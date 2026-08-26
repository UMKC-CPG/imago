!! SPDX-License-Identifier: ECL-2.0
!! Copyright (c) 2026 Paul Rulis

module O_SecularEquation

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Define module data.
   real (kind=double), allocatable, dimension (:,:,:) :: energyEigenValues
         ! States, KPoints, Spin

#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:) :: valeVale
   complex (kind=double), allocatable, dimension (:,:)   :: valeValeOL
   complex (kind=double), allocatable, dimension (:,:,:) :: valeValeMM
   complex (kind=double), allocatable, dimension (:,:,:) :: valeValeKO
#else
   real (kind=double), allocatable, dimension (:,:,:) :: valeValeGamma
   real (kind=double), allocatable, dimension (:,:)   :: valeValeOLGamma
   real (kind=double), allocatable, dimension (:,:,:) :: valeValeMMGamma
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains


subroutine secularEqnSCF(spinDirection, numStates)

   ! Import necessary modules.
   use HDF5
   use O_Kinds
   use O_TimeStamps
   use O_KPoints, only: numKPoints
   use O_Potential, only: rel, spin, potDim, potCoeffs, numPlusUJAtoms, &
         & currIteration
   use O_AtomicSites, only: valeDim
   use O_SCFEigValHDF5, only: eigenValues_did, states
   use O_SCFEigVecHDF5, only: eigenVectors_did, eigenVectors_aid, valeStates
   use O_SCFIntegralsHDF5, only: packedVVDims, atomOverlap_did, &
         & atomKEOverlap_did, atomMVOverlap_did, atomNPOverlap_did, &
         & atomPotOverlap_did
   use O_MPI, only: mpiSize, sendCtrlMPI, sendPackedMPI, &
         & recvDblVecMPI, solveTask, mpiTagHam, mpiTagOvlp, &
         & mpiTagVals, mpiTagVecs
   use O_ELPASolve, only: elpaAvailable, rootCollectiveSolve
#ifndef GAMMA
   use O_MPI, only: recvCmplxBlockMPI
   use O_LAPACKZHEGV
   use O_MatrixSubs, only: readPackedMatrix, readPackedMatrixAccum,&
         & unpackMatrix
#else
   use O_MPI, only: recvPackedMPI
   use O_LAPACKDSYGV
   use O_MatrixSubs, only: readPackedMatrix, readPackedMatrixAccum, &
         & unpackMatrixGamma
#endif

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define the passed parameters.
   integer :: spinDirection
   integer :: numStates

   ! Define the local variables used in this subroutine.
   integer :: i,j ! Loop index variables
!integer :: k, l, m, n, o, p
   integer :: hdferr
   integer :: hdf5Status
   integer(hsize_t), dimension (1) :: attribIntDims ! Attribute dataspace dim
   integer :: dim1
   real    (kind=double), allocatable, dimension (:,:)   :: packedValeVale
   real    (kind=double), allocatable, dimension (:,:)   :: tempPackedValeVale
!complex(kind=double), allocatable, dimension(:,:) :: identity

   ! The k-point deal (PSEUDOCODE 26.2/26.3): which k-points an
   !   earlier pass found complete, which rank owns each undone one,
   !   and the buffers the dispatch and collection use. A dealt
   !   k-point's packed overlap travels in its own buffer because
   !   the Hamiltonian occupies packedValeVale at the same time.
   real    (kind=double), allocatable, dimension (:,:)   :: packedOverlap
   integer, allocatable, dimension (:) :: kPointDone
   integer, allocatable, dimension (:) :: kPointOwner
   integer, allocatable, dimension (:) :: kPointRound
   integer :: numShipped ! K-points dealt to worker ranks this call.
   integer :: dealWidth  ! Ranks in the deal: mpiSize, or 1 when the
         ! deal is off (one rank, one k-point, or Hubbard-U).
   integer :: numRounds  ! Deal rounds; each hands every rank at most
         ! one k-point (the deadlock-free discipline of 26.3).
   integer :: roundIdx
   logical :: useCollective ! The stage-B arm (PSEUDOCODE 27): one
         ! k-point, several ranks, ELPA built -- the single solve is
         ! distributed instead of dealt.
   real (kind=double), allocatable, dimension (:) :: valsPlusTime
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:) :: eigVecBlock
#else
   real (kind=double), allocatable, dimension (:,:) :: eigVecBlockGamma
#endif
   real (kind=double) :: minSolveSeconds, maxSolveSeconds

   ! Record the date and time that we start.
   call timeStampStart (15)

   ! Initialize the dimension of the packed matrices to include two components
   !   (real,imaginary) or just a real component.
#ifndef GAMMA
   dim1 = 2
#else
   dim1 = 1
#endif

   ! Only allocate for first spin to prevent double allocation.  These arrays
   !   will be deallocated in the makeValenceRho subroutine to accomodate the
   !   common case of one kpoint SCF calculations where it is not necessary
   !   (or efficient) to write the energy eigen values and wave function to
   !   disk with each iteration.  These matrices can simply be kept in memory
   !   and used directly when needed in makeValenceRho.  Note that only 1
   !   kpoint is needed at a time, that is the meaning of the "1" in the
   !   valeVale and valeValeOL allocation statements.
   if (spinDirection == 1) then
      allocate(energyEigenValues (numStates,numKPoints,spin))
#ifndef GAMMA
      allocate(valeVale          (valeDim,valeDim,spin)) ! Complex
#else
      allocate(valeValeGamma     (valeDim,valeDim,spin)) ! Real
#endif
   endif

   ! These matrices are deallocated at the end of this subroutine.
#ifndef GAMMA
   allocate(valeValeOL         (valeDim,valeDim)) ! Complex
#else
   allocate(valeValeOLGamma    (valeDim,valeDim)) ! Real
#endif
   allocate(packedValeVale     (dim1,valeDim*(valeDim+1)/2))
   allocate(tempPackedValeVale (dim1,valeDim*(valeDim+1)/2))
   allocate(packedOverlap      (dim1,valeDim*(valeDim+1)/2))
   allocate(kPointDone  (numKPoints))
   allocate(kPointOwner (numKPoints))
   allocate(kPointRound (numKPoints))

   ! Read every k-point's completion attribute up front: the deal
   !   below must know the undone set before anything is shipped
   !   (PSEUDOCODE 26.2). The attribute handle is NOT closed for a
   !   completed k-point -- the reset loop at the end of this
   !   routine writes through the same handle, so closing it here
   !   would make that write fail on a mid-iteration restart.
   attribIntDims(1) = 1
   do i = 1, numKPoints
      hdf5Status = 0
      call h5aread_f(eigenVectors_aid(i,spinDirection),H5T_NATIVE_INTEGER,&
            & hdf5Status,attribIntDims,hdferr)
      if (hdferr /= 0) stop 'Failed to read eigen vector SCF status.'
      kPointDone(i) = hdf5Status
      if (hdf5Status == 1) then
         write(20,*) "Wave function for kpoint ",i," already computed."
      endif
   enddo

   ! The deal (PSEUDOCODE 26.2/26.3): round-robin over the undone
   !   k-points, root taking a share like any rank, in ROUNDS of at
   !   most one k-point per rank -- the discipline that keeps the
   !   blocking transport deadlock-free (a worker never holds a
   !   second pending task while its reply is unsent; the first
   !   acceptance run proved the alternative deadlocks once the
   !   matrices pass MPI's eager threshold). K-points carrying the
   !   Hubbard-U correction are never dealt (the UJ path couples to
   !   the previous iteration's density matrix held here on root),
   !   and a width of one -- one rank, or one k-point -- degenerates
   !   to today's all-on-root path.
   kPointOwner(:) = 0
   kPointRound(:) = 0
   numShipped = 0
   dealWidth = 1
   if ((mpiSize > 1) .and. (numKPoints > 1) .and. &
         & (numPlusUJAtoms == 0)) then
      dealWidth = mpiSize
   endif

   ! The arm choice (PSEUDOCODE 27.1): at one k-point the deal cannot
   !   help, so the single solve is distributed across every rank
   !   with ELPA -- when it is built, there are ranks to use, and no
   !   Hubbard-U coupling holds the solve on root. Everything outside
   !   the solver call itself is identical between the arms.
   useCollective = elpaAvailable .and. (mpiSize > 1) .and. &
         & (numKPoints == 1) .and. (numPlusUJAtoms == 0) .and. &
         & (valeDim >= mpiSize)
   j = 0
   do i = 1, numKPoints
      if (kPointDone(i) == 1) cycle
      kPointOwner(i) = mod(j, dealWidth)
      kPointRound(i) = j / dealWidth
      if (kPointOwner(i) /= 0) numShipped = numShipped + 1
      j = j + 1
   enddo
   numRounds = (j + dealWidth - 1) / dealWidth

   if (numShipped > 0) then
      allocate (valsPlusTime (numStates + 1))
#ifndef GAMMA
      allocate (eigVecBlock (valeDim, numStates))
#else
      allocate (eigVecBlockGamma (valeDim, numStates))
#endif
      minSolveSeconds = huge (minSolveSeconds)
      maxSolveSeconds = 0.0_double
   endif

   ! The round loop (PSEUDOCODE 26.3): ship this round's worker
   !   k-points (each worker is idle in recvCtrl, so the sends
   !   complete), solve root's own k-point while they work, then
   !   collect the round's replies before dealing the next.
   do roundIdx = 0, numRounds - 1

   ! Dispatch this round's worker-owned k-points.
   do i = 1, numKPoints
      if ((kPointDone(i) == 1) .or. (kPointOwner(i) == 0) .or. &
            & (kPointRound(i) /= roundIdx)) cycle
      call assembleSecularSCF (i, spinDirection, packedValeVale,&
            & packedOverlap, tempPackedValeVale)
      call sendCtrlMPI (solveTask, i, kPointOwner(i))
      call sendPackedMPI (packedValeVale, kPointOwner(i), mpiTagHam)
      call sendPackedMPI (packedOverlap, kPointOwner(i), mpiTagOvlp)
   enddo

   ! Root's own k-point this round: today's serial body, unchanged
   !   in effect. (The loop finds at most one match per round; it
   !   stays a loop so the width-one path reads as today's code.)
   do i = 1,numKPoints

      ! Skip k-points already computed (the mid-iteration restart),
      !   k-points dealt to a worker (collected below), and other
      !   rounds' k-points.
      if ((kPointDone(i) == 1) .or. (kPointOwner(i) /= 0) .or. &
            & (kPointRound(i) /= roundIdx)) cycle

      ! Prepare the matrices.
      packedValeVale(:,:) = 0.0_double
      tempPackedValeVale(:,:) = 0.0_double

#ifndef GAMMA
      valeVale(:,:,spinDirection) = cmplx(0.0_double,0.0_double,double)
      valeValeOL(:,:) = cmplx(0.0_double,0.0_double,double)
#else
      valeValeGamma(:,:,spinDirection) = 0.0_double
      valeValeOLGamma(:,:) = 0.0_double
#endif

      ! Assemble this k-point's packed Hamiltonian and overlap --
      !   the SAME assembly the dispatch loop ships to workers
      !   (nuclear + kinetic [+ mass velocity] + the potDim potential
      !   terms weighted by this iteration's potCoeffs; see
      !   assembleSecularSCF).
      call assembleSecularSCF (i, spinDirection, packedValeVale,&
            & packedOverlap, tempPackedValeVale)

      ! Unpack the hamiltonian and overlap matrices. The LAPACK arm
      !   reads only the upper triangle its solvers consume; the
      !   collective arm unpacks BOTH triangles, because ELPA takes
      !   the full Hermitian matrix (an upper-only scatter fed it
      !   zeros for half the matrix -- measured 2026-08-22 as
      !   garbage eigenvectors and a positive total energy).
      if (useCollective) then
#ifndef GAMMA
         call unpackMatrix(valeVale(:,:,spinDirection),&
               & packedValeVale,valeDim,1)
         call unpackMatrix(valeValeOL(:,:),packedOverlap,valeDim,1)
#else
         call unpackMatrixGamma(valeValeGamma(:,:,spinDirection),&
               & packedValeVale,valeDim,1)
         call unpackMatrixGamma(valeValeOLGamma(:,:),packedOverlap,&
               & valeDim,1)
#endif
      else
#ifndef GAMMA
      call unpackMatrix(valeVale(:,:,spinDirection),packedValeVale,valeDim,0)
      call unpackMatrix(valeValeOL(:,:),packedOverlap,valeDim,0)
#else
      call unpackMatrixGamma(valeValeGamma(:,:,spinDirection),&
            & packedValeVale,valeDim,0)
      call unpackMatrixGamma(valeValeOLGamma(:,:),packedOverlap,valeDim,0)
#endif
      endif

      ! For each atom with a Hubbard U and Hund J term, we need to apply its
      !   effect on relevant matrix elements of the Hamiltonian. This is only
      !   needed in the event that we actually have atoms with these terms.
      ! Prior to application of the UJ effect on the Hamiltonian, we need to
      !   complete the update process that was started in the makeValenceRho
      !   subroutine of the valeCharge.F90 file. At that point we had the charge
      !   density matrix in a convenient structure (i.e. an unpacked matrix as
      !   opposed to a packed matrix), but we didn't have the overlap matrix in
      !   a convenient form (i.e. it would have been packed). Therefore, to
      !   multiply the charge density matrix terms by the appropriate overlap
      !   matrix terms would have required us to do a bit of annoying triangle
      !   math which at the moment I don't have time to develop. (Although it is
      !   probably fairly easy-ish.) At any rate, we have access to the overlap
      !   matrix elements in an unpacked matrix now. So, we will multiply them
      !   against the stored charge density matrix elements from the previous
      !   iteration before we apply the result to the Hamiltonian.
      ! On the first iteration of an SCF calculation there are no stored results
      !   from the previous iteration. Therefore, we will just have zeros for
      !   the plusUJ elements and there will be no need to modify the
      !   Hamiltonian.
      ! As a reminder, the reason that we need two phases for the update process
      !   is that although the charge density matrix is obtained through the
      !   Psi* * Psi operation 
      if ((numPlusUJAtoms > 0) .and. (currIteration > 1)) then
         call update2AndApplyUJ(i,spinDirection)
      endif

      ! Solve the eigen problem: collectively with ELPA when the arm
      !   policy chose it (the workers were woken by
      !   rootCollectiveSolve; the eigenvectors land in the leading
      !   columns exactly as the LAPACK contract puts them), or with
      !   the serial LAPACK routine as always.
      if (useCollective) then
#ifndef GAMMA
         call rootCollectiveSolve (valeVale(:,:,spinDirection),&
               & valeValeOL(:,:),&
               & energyEigenValues(:,i,spinDirection))
#else
         call rootCollectiveSolve (valeValeGamma(:,:,spinDirection),&
               & valeValeOLGamma(:,:),&
               & energyEigenValues(:,i,spinDirection))
#endif
      else
#ifndef GAMMA
!write(20,*) "valeValeSCF"
!do j = 1,valeDim
!   do k = 1,valeDim
!      write(20,fmt="(2e12.4)",advance="NO") valeVale(k,j,spinDirection)
!   enddo
!   write(20,*)
!enddo
      call solveZHEGV(valeDim,numStates,valeVale(:,:,spinDirection),&
            & valeValeOL(:,:),energyEigenValues(:,i,spinDirection))
#else
      call solveDSYGV(valeDim,numStates,valeValeGamma(:,:,spinDirection),&
            & valeValeOLGamma(:,:),energyEigenValues(:,i,spinDirection))
#endif
      endif

!! Test normalization.
!allocate (identity(numStates,numStates))
!! Read the atomic overlap matrix. 
!call readPackedMatrix(atomOverlap_did(i),packedValeVale,&
!      & packedVVDims,dim1,valeDim)
!
!      ! Unpack the overlap matrix.
!#ifndef GAMMA
!call unpackMatrix(valeValeOL(:,:,1,1),packedValeVale,valeDim,0)
!#else
!call unpackMatrixGamma(valeValeOLGamma(:,:,1),packedValeVale,valeDim,0)
!#endif
!
!identity(:numStates,:numStates) = &
!      & matmul(&
!      & matmul(transpose(conjg(valeVale(:valeDim,:numStates,1,spinDirection))),&
!             & valeValeOL(:,:,1,1)),valeVale(:,:numStates,1,spinDirection))
!write(24,*) "Normalization i = ",i
!do k = 1, num_states
!write(24,*) "k,I(k,k) = ",k,identity(k,k)
!enddo
!
!deallocate (identity)

      ! Write the energy eigenValues onto disk in HDF5 format in a.u.
      call h5dwrite_f (eigenValues_did(i,spinDirection),H5T_NATIVE_DOUBLE,&
            & energyEigenValues(:,i,spinDirection),states,hdferr)
      if (hdferr /= 0) stop 'Cannot write energy eigen values SCF.'

      ! In the event that we have some atoms with plusUJ terms, we need to
      !   update the terms. The update depends on the charge in each of the
      !   highest d or f orbitals of the affected atoms. Therefore, we will
      !   need to read in the overlap matrix again because it was just
      !   destroyed in the ZHEGV solution. An alternative approach might be
      !   to make a copy of the overlap matrix and hold it in reserve until we
      !   need it for this calculation.
      if (numPlusUJAtoms > 0) then

         ! Read the atomic overlap matrix. 
         call readPackedMatrix(atomOverlap_did(i),packedValeVale,&
               & packedVVDims,dim1,valeDim)

         ! Unpack the overlap matrix.
#ifndef GAMMA
         call unpackMatrix(valeValeOL(:,:),packedValeVale,valeDim,0)
#else
         call unpackMatrixGamma(valeValeOLGamma(:,:),packedValeVale,valeDim,0)
#endif
      endif

#ifndef GAMMA
      ! Write the eigenVectors onto disk in HDF5 format for this
      !   kpoint and spin direction.
      call h5dwrite_f(eigenVectors_did(1,i,spinDirection),&
            & H5T_NATIVE_DOUBLE,real(valeVale(:,:numStates,spinDirection),&
            & double),valeStates,hdferr)
      if (hdferr /= 0) stop 'Cannot write real energy eigen vectors SCF.'
      call h5dwrite_f(eigenVectors_did(2,i,spinDirection),&
            & H5T_NATIVE_DOUBLE,aimag(valeVale(:,:numStates,spinDirection)),&
            & valeStates,hdferr)
      if (hdferr /= 0) stop 'Cannot write imag energy eigen vectors SCF.'
#else
      call h5dwrite_f(eigenVectors_did(1,i,spinDirection),&
            & H5T_NATIVE_DOUBLE,valeValeGamma(:,:numStates,spinDirection),&
            & valeStates,hdferr)
      if (hdferr /= 0) stop 'Cannot write real energy eigen vectors SCF.'
#endif

      ! Record that this calculation is complete.
      call h5awrite_f(eigenVectors_aid(i,spinDirection),H5T_NATIVE_INTEGER,&
            & 1,attribIntDims,hdferr)
      if (hdferr /= 0) stop 'Failed to record eigenvector success SCF.'
   enddo ! Loop i over kpoints.

   ! Collect this round's dealt k-points (PSEUDOCODE 26.3): each
   !   owner returns the eigenvalues (with its solve seconds
   !   appended) and the lowest numStates eigenvector columns, and
   !   root writes them exactly where its own solves write --
   !   populate, the density, and every later consumer see the
   !   serial state.
   if (numShipped > 0) then
      do i = 1, numKPoints
         if ((kPointDone(i) == 1) .or. (kPointOwner(i) == 0) .or. &
               & (kPointRound(i) /= roundIdx)) cycle

         call recvDblVecMPI (valsPlusTime, kPointOwner(i), mpiTagVals)
         energyEigenValues(:,i,spinDirection) = &
               & valsPlusTime(1:numStates)
         minSolveSeconds = min (minSolveSeconds,&
               & valsPlusTime(numStates+1))
         maxSolveSeconds = max (maxSolveSeconds,&
               & valsPlusTime(numStates+1))

         ! Write the energy eigenValues onto disk in a.u.
         call h5dwrite_f (eigenValues_did(i,spinDirection),&
               & H5T_NATIVE_DOUBLE,energyEigenValues(:,i,spinDirection),&
               & states,hdferr)
         if (hdferr /= 0) stop 'Cannot write energy eigen values SCF.'

#ifndef GAMMA
         call recvCmplxBlockMPI (eigVecBlock, kPointOwner(i),&
               & mpiTagVecs)
         call h5dwrite_f(eigenVectors_did(1,i,spinDirection),&
               & H5T_NATIVE_DOUBLE,real(eigVecBlock(:,:),double),&
               & valeStates,hdferr)
         if (hdferr /= 0) stop &
               & 'Cannot write real energy eigen vectors SCF.'
         call h5dwrite_f(eigenVectors_did(2,i,spinDirection),&
               & H5T_NATIVE_DOUBLE,aimag(eigVecBlock(:,:)),&
               & valeStates,hdferr)
         if (hdferr /= 0) stop &
               & 'Cannot write imag energy eigen vectors SCF.'
#else
         call recvPackedMPI (eigVecBlockGamma, kPointOwner(i),&
               & mpiTagVecs)
         call h5dwrite_f(eigenVectors_did(1,i,spinDirection),&
               & H5T_NATIVE_DOUBLE,eigVecBlockGamma(:,:),&
               & valeStates,hdferr)
         if (hdferr /= 0) stop &
               & 'Cannot write real energy eigen vectors SCF.'
#endif

         ! Record that this calculation is complete.
         call h5awrite_f(eigenVectors_aid(i,spinDirection),&
               & H5T_NATIVE_INTEGER,1,attribIntDims,hdferr)
         if (hdferr /= 0) stop 'Failed to record eigenvector success.'
      enddo
   endif

   enddo ! roundIdx over the deal rounds.

   if (numShipped > 0) then
      ! One summary line: the k-point deal's record for this call.
      write (20,*) 'Secular k-point deal: ', numShipped,&
            & ' shipped over ', mpiSize - 1, ' workers; solve s ',&
            & minSolveSeconds, maxSolveSeconds
      deallocate (valsPlusTime)
#ifndef GAMMA
      deallocate (eigVecBlock)
#else
      deallocate (eigVecBlockGamma)
#endif
   endif

   ! Once all kpoints are done, we need to reset the attributes so that the
   !   next iteration doesn't think that all the kpoints are done already.
   do i = 1, numKPoints
      call h5awrite_f(eigenVectors_aid(i,spinDirection),H5T_NATIVE_INTEGER,&
            & 0,attribIntDims,hdferr)
      if (hdferr /= 0) stop 'Failed to reset eigenvector success to 0 SCF.'
   enddo

   ! Deallocate unnecessary arrays and matrices.
   deallocate (packedValeVale)
   deallocate (tempPackedValeVale)
   deallocate (packedOverlap)
   deallocate (kPointDone)
   deallocate (kPointOwner)
   deallocate (kPointRound)
#ifndef GAMMA
   deallocate(valeValeOL)
#else
   deallocate(valeValeOLGamma)
#endif

   ! Record the date and time that we finish.
   call timeStampEnd (15)

end subroutine secularEqnSCF


subroutine assembleSecularSCF (kPointIndex, spinDirection, packedHam,&
      & packedOvlp, scratch)

   ! Assemble one k-point's packed Hamiltonian and packed overlap
   !   from root's open HDF5 handles (PSEUDOCODE 26.3): the nuclear
   !   potential, the kinetic energy, the mass velocity when the
   !   calculation is scalar relativistic, and the potDim
   !   three-centre potential terms weighted by this iteration's
   !   potential coefficients; then the overlap matrix. This is THE
   !   assembly of the secular problem, factored out so that root's
   !   own solves and the k-point deal's shipments cannot drift
   !   apart. It runs on root only -- workers have no file handles.

   ! Import necessary modules.
   use HDF5
   use O_Kinds
   use O_Potential, only: rel, potDim, potCoeffs
   use O_AtomicSites, only: valeDim
   use O_SCFIntegralsHDF5, only: packedVVDims, atomOverlap_did, &
         & atomKEOverlap_did, atomMVOverlap_did, atomNPOverlap_did, &
         & atomPotOverlap_did
   use O_MatrixSubs, only: readPackedMatrix, readPackedMatrixAccum

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define the passed parameters.
   integer, intent (in) :: kPointIndex
   integer, intent (in) :: spinDirection
   real (kind=double), dimension (:,:), intent (out) :: packedHam
   real (kind=double), dimension (:,:), intent (out) :: packedOvlp
   real (kind=double), dimension (:,:), intent (inout) :: scratch

   ! Define the local variables used in this subroutine.
   integer :: j
   integer :: dim1

   dim1 = size (packedHam, 1)
   packedHam(:,:) = 0.0_double
   scratch(:,:) = 0.0_double

   ! Read the nuclear potential term into the packed hamiltonian.
   call readPackedMatrix(atomNPOverlap_did(kPointIndex),packedHam,&
         & packedVVDims,dim1,valeDim)

   ! Read the kinetic energy term into the still packed hamiltonian.
   call readPackedMatrixAccum(atomKEOverlap_did(kPointIndex),&
         & packedHam,scratch,packedVVDims,0.0_double,dim1,valeDim)

   ! Read the mass velocity term into the still packed hamiltonian
   !   if we are doing a scalar relativistic calculation.
   if (rel == 1) then
      ! Note that the -1.0 introduce a negative sign to the term. In
      !   the future, the sign should be incorporated into the matrix
      !   calculation itself to avoid the extra work here.
      call readPackedMatrixAccum(atomMVOverlap_did(kPointIndex),&
            & packedHam,scratch,packedVVDims,0.0_double,dim1,valeDim)
   endif

   ! Read the atomic potential terms into the still packed
   !   hamiltonian, each weighted by this iteration's coefficient.
   do j = 1, potDim
      call readPackedMatrixAccum(atomPotOverlap_did(kPointIndex,j),&
            & packedHam,scratch,packedVVDims,&
            & potCoeffs(j,spinDirection),dim1,valeDim)
   enddo

   ! Read the atomic overlap matrix into its own packed buffer.
   call readPackedMatrix(atomOverlap_did(kPointIndex),packedOvlp,&
         & packedVVDims,dim1,valeDim)

end subroutine assembleSecularSCF


subroutine solveServerLoop

   ! The worker side of the k-point deal (PSEUDOCODE 26.4). Between
   !   the term stage and the certificate barrier, every worker rank
   !   sits here: receive a control message from root; on a TASK,
   !   receive the k-point's packed Hamiltonian and overlap, unpack,
   !   solve with the SAME serial backend root uses, and return the
   !   eigenvalues (solve seconds appended) and the lowest numStates
   !   eigenvector columns; on SHUTDOWN, return to the caller, which
   !   falls to the certificate barrier. Everything needed beyond
   !   the messages -- valeDim, numStates, the LAPACK block size --
   !   is replicated setup state this rank already holds.
   !
   ! The control-code dispatch is deliberately extensible: stage B
   !   (ELPA) adds a collective-solve code here rather than a new
   !   worker structure.

   ! Import necessary modules.
   use O_Kinds
   use O_Input, only: numStates
   use O_AtomicSites, only: valeDim
   use O_LAPACKParameters, only: setBlockSize
   use O_MPI, only: recvCtrlMPI, recvPackedMPI, sendDblVecMPI, &
         & solveShutdown, solveTask, solveCollective, elecStatTask, &
         & stopMPI, mpiTagHam, mpiTagOvlp, mpiTagVals, mpiTagVecs
   use O_ELPASolve, only: workerCollectiveSolve
   use O_ElectroStatics, only: neutralAndNuclearQPot, residualQ
#ifndef GAMMA
   use O_MPI, only: sendCmplxBlockMPI
   use O_LAPACKZHEGV
   use O_MatrixSubs, only: unpackMatrix
#else
   use O_MPI, only: sendPackedMPI
   use O_LAPACKDSYGV
   use O_MatrixSubs, only: unpackMatrixGamma
#endif

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define the local variables used in this subroutine.
   integer :: code
   integer :: kPointIndex
   integer :: dim1
   integer :: clockBefore, clockAfter, clockRate
   real (kind=double), allocatable, dimension (:,:) :: packedHam
   real (kind=double), allocatable, dimension (:,:) :: packedOvlp
   real (kind=double), allocatable, dimension (:)   :: valsPlusTime
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:) :: hamiltonian
   complex (kind=double), allocatable, dimension (:,:) :: overlap
#else
   real (kind=double), allocatable, dimension (:,:) :: hamiltonian
   real (kind=double), allocatable, dimension (:,:) :: overlap
#endif

#ifndef GAMMA
   dim1 = 2
#else
   dim1 = 1
#endif

   ! Set up LAPACK machine parameters, mirroring mainSCF.
   call setBlockSize (valeDim)

   allocate (packedHam  (dim1, valeDim*(valeDim+1)/2))
   allocate (packedOvlp (dim1, valeDim*(valeDim+1)/2))
   allocate (hamiltonian (valeDim, valeDim))
   allocate (overlap     (valeDim, valeDim))
   allocate (valsPlusTime (numStates + 1))

   do
      call recvCtrlMPI (code, kPointIndex, 0)
      if (code == solveShutdown) exit
      if (code == solveCollective) then
         ! Join the collective ELPA solve (PSEUDOCODE 27.4) and
         !   return to the loop -- the extensible-control design of
         !   26.4 doing exactly what it was built for.
         call workerCollectiveSolve (valeDim, numStates)
         write (20,*) 'Joined collective solve for k-point ',&
               & kPointIndex
         cycle
      endif
      if (code == elecStatTask) then
         ! Join a dealt electrostatic-setup sub-stage (PSEUDOCODE
         !   28): the same subroutine root is inside, whose
         !   root-only entry and exit blocks skip themselves on a
         !   worker; the second control integer names the
         !   sub-stage.
         if (kPointIndex == 1) then
            call neutralAndNuclearQPot
         else
            call residualQ
         endif
         write (20,*) 'Joined electrostatics sub-stage ', kPointIndex
         cycle
      endif
      if (code /= solveTask) then
         call stopMPI ('Solve server received an unknown code.')
      endif

      call recvPackedMPI (packedHam, 0, mpiTagHam)
      call recvPackedMPI (packedOvlp, 0, mpiTagOvlp)

#ifndef GAMMA
      hamiltonian(:,:) = cmplx (0.0_double, 0.0_double, double)
      overlap(:,:)     = cmplx (0.0_double, 0.0_double, double)
      call unpackMatrix (hamiltonian, packedHam, valeDim, 0)
      call unpackMatrix (overlap, packedOvlp, valeDim, 0)
#else
      hamiltonian(:,:) = 0.0_double
      overlap(:,:)     = 0.0_double
      call unpackMatrixGamma (hamiltonian, packedHam, valeDim, 0)
      call unpackMatrixGamma (overlap, packedOvlp, valeDim, 0)
#endif

      ! Solve with the same backend root uses; the eigenvectors land
      !   in the hamiltonian's columns, exactly as on root.
      call system_clock (clockBefore, clockRate)
#ifndef GAMMA
      call solveZHEGV (valeDim, numStates, hamiltonian, overlap,&
            & valsPlusTime(1:numStates))
#else
      call solveDSYGV (valeDim, numStates, hamiltonian, overlap,&
            & valsPlusTime(1:numStates))
#endif
      call system_clock (clockAfter)
      valsPlusTime(numStates+1) = real (clockAfter - clockBefore,&
            & double) / real (clockRate, double)

      ! Log the service in this rank's own log (visible under
      !   IMAGO_RANK_LOGS; discarded by default).
      write (20,*) 'Served k-point ', kPointIndex, ' in ',&
            & valsPlusTime(numStates+1), ' s'

      call sendDblVecMPI (valsPlusTime, 0, mpiTagVals)
#ifndef GAMMA
      call sendCmplxBlockMPI (hamiltonian(:,1:numStates), 0,&
            & mpiTagVecs)
#else
      call sendPackedMPI (hamiltonian(:,1:numStates), 0, mpiTagVecs)
#endif
   enddo

   deallocate (packedHam)
   deallocate (packedOvlp)
   deallocate (hamiltonian)
   deallocate (overlap)
   deallocate (valsPlusTime)

end subroutine solveServerLoop


subroutine secularEqnPSCF(spinDirection,numStates,numComponents,ol_did,&
      & ham_did,eVal_did,eVec_did,eVec_aid)

   ! Import necessary modules.
   use HDF5
   use O_Kinds
   use O_TimeStamps
   use O_KPoints, only: numKPoints
   use O_Potential, only: spin, numPlusUJAtoms, currIteration
   use O_AtomicSites, only: valeDim
   use O_PSCFEigValHDF5, only: states
   use O_PSCFEigVecHDF5, only: valeStatesPSCF
   use O_PSCFIntegralsHDF5, only: packedVVDimsPSCF
#ifndef GAMMA
   use O_LAPACKZHEGV
   use O_MatrixSubs, only: readPackedMatrix, readPackedMatrixAccum,&
         & unpackMatrix
#else
   use O_LAPACKDSYGV
   use O_MatrixSubs, only: readPackedMatrix, readPackedMatrixAccum, &
         & unpackMatrixGamma
#endif

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define the passed parameters.
   integer, intent(in) :: spinDirection
   integer, intent(in) :: numStates
   integer, intent(in) :: numComponents
   integer(hid_t), dimension(numKPoints), intent(in) :: ol_did
   integer(hid_t), dimension(numKPoints,spin), intent(in) :: ham_did
   integer(hid_t), dimension(numKPoints,spin), intent(in) :: eVal_did
   integer(hid_t), dimension(numComponents,numKPoints,spin), &
         & intent(in) :: eVec_did
   integer(hid_t), dimension(numKPoints,spin), intent(in) :: eVec_aid

   ! Define the local variables used in this subroutine.
   integer :: i ! Loop index variables
!integer :: k, l, m, n, o, p
   integer :: hdferr
   integer :: hdf5Status
   integer(hsize_t), dimension (1) :: attribIntDims ! Attribute dataspace dim
   integer :: dim1
   real    (kind=double), allocatable, dimension (:,:)   :: packedValeVale
   real    (kind=double), allocatable, dimension (:,:)   :: tempPackedValeVale
!complex(kind=double), allocatable, dimension(:,:) :: identity

   ! Record the date and time that we start.
   call timeStampStart (15)

   ! Initialize the dimension of the packed matrices to include two components
   !   (real,imaginary) or just a real component.
#ifndef GAMMA
   dim1 = 2
#else
   dim1 = 1
#endif

   ! Only allocate for first spin to prevent double allocation because we do
   !   spin sequentially. All these arrays are deallocated !  . These arrays
   !   will be deallocated in the makeValenceRho subroutine to accomodate the
   !   common case of one kpoint SCF calculations where it is not necessary
   !   (or efficient) to write the energy eigen values and wave function to
   !   disk with each iteration.  These matrices can simply be kept in memory
   !   and used directly when needed in makeValenceRho.  Note that only 1
   !   kpoint is needed at a time, that is the meaning of the "1" in the
   !   valeVale and valeValeOL allocation statements.
   if (spinDirection == 1) then
      allocate(energyEigenValues (numStates,numKPoints,spin))
#ifndef GAMMA
      allocate(valeVale          (valeDim,valeDim,spin)) ! Complex
#else
      allocate(valeValeGamma     (valeDim,valeDim,spin)) ! Real
#endif
   endif

   ! These matrices are deallocated at the end of this subroutine.
#ifndef GAMMA
   allocate(valeValeOL         (valeDim,valeDim)) ! Complex
#else
   allocate(valeValeOLGamma    (valeDim,valeDim)) ! Real
#endif
   allocate(packedValeVale     (dim1,valeDim*(valeDim+1)/2))
   allocate(tempPackedValeVale (dim1,valeDim*(valeDim+1)/2))


   ! Begin loop over all kpoints.
   do i = 1,numKPoints

      ! If the eigenvectors for this kpoint and spin direction are already
      !   computed, then skip.
      hdf5Status = 0
      attribIntDims(1) = 1
      call h5aread_f(eVec_aid(i,spinDirection),H5T_NATIVE_INTEGER,&
            & hdf5Status,attribIntDims,hdferr)
      if (hdferr /= 0) stop 'Failed to read eigen vector PSCF status.'
      if (hdf5Status == 1) then
         write(20,*) "Wave function for kpoint ",i," already computed."
         call h5aclose_f(eVec_aid(i,spinDirection),hdferr)
         if (hdferr /= 0) stop 'Failed to close eigen vector PSCF status'
         call h5dread_f(eVal_did(i,spinDirection),H5T_NATIVE_DOUBLE,&
               & energyEigenValues(:numStates,i,spinDirection),states,hdferr)
         if (hdferr /= 0) stop 'Failed to read EVals PSCF'
         cycle
      endif

      ! Prepare the matrices.
      packedValeVale(:,:) = 0.0_double
      tempPackedValeVale(:,:) = 0.0_double

#ifndef GAMMA
      valeVale(:,:,spinDirection) = cmplx(0.0_double,0.0_double,double)
      valeValeOL(:,:) = cmplx(0.0_double,0.0_double,double)
#else
      valeValeGamma(:,:,spinDirection) = 0.0_double
      valeValeOLGamma(:,:) = 0.0_double
#endif

      ! Read the combined Hamiltonian term into the packed valeVale.
      call readPackedMatrix(ham_did(i,spinDirection),&
            & packedValeVale,packedVVDimsPSCF,dim1,valeDim)

      ! Unpack the hamiltonian matrix.
#ifndef GAMMA
      call unpackMatrix(valeVale(:,:,spinDirection),packedValeVale,valeDim,0)
#else
      call unpackMatrixGamma(valeValeGamma(:,:,spinDirection),&
            & packedValeVale,valeDim,0)
#endif

      ! Read the atomic overlap matrix. 
      call readPackedMatrix(ol_did(i),packedValeVale,packedVVDimsPSCF,dim1,&
            & valeDim)

      ! Unpack the overlap matrix.
#ifndef GAMMA
      call unpackMatrix(valeValeOL(:,:),packedValeVale,valeDim,0)
#else
      call unpackMatrixGamma(valeValeOLGamma(:,:),packedValeVale,valeDim,0)
#endif

      ! For each atom with a Hubbard U and Hund J term, we need to apply its
      !   effect on relevant matrix elements of the Hamiltonian. This is only
      !   needed in the event that we actually have atoms with these terms.
      ! Prior to application of the UJ effect on the Hamiltonian, we need to
      !   complete the update process that was started in the makeValenceRho
      !   subroutine of the valeCharge.F90 file. At that point we had the charge
      !   density matrix in a convenient structure (i.e. an unpacked matrix as
      !   opposed to a packed matrix), but we didn't have the overlap matrix in
      !   a convenient form (i.e. it would have been packed). Therefore, to
      !   multiply the charge density matrix terms by the appropriate overlap
      !   matrix terms would have required us to do a bit of annoying triangle
      !   math which at the moment I don't have time to develop. (Although it is
      !   probably fairly easy-ish.) At any rate, we have access to the overlap
      !   matrix elements in an unpacked matrix now. So, we will multiply them
      !   against the stored charge density matrix elements from the previous
      !   iteration before we apply the result to the Hamiltonian.
      ! On the first iteration of an SCF calculation there are no stored results
      !   from the previous iteration. Therefore, we will just have zeros for
      !   the plusUJ elements and there will be no need to modify the
      !   Hamiltonian.
      ! As a reminder, the reason that we need two phases for the update process
      !   is that although the charge density matrix is obtained through the
      !   Psi* * Psi operation 
      if ((numPlusUJAtoms > 0) .and. (currIteration > 1)) then
         call update2AndApplyUJ(i,spinDirection)
      endif

      ! Solve the eigen problem with a LAPACK routine.
#ifndef GAMMA
!write(20,*) "valeValePSCF"
!do j = 1,valeDim
!   do k = 1,valeDim
!      write(20,fmt="(2e12.4)",advance="NO") valeVale(k,j,spinDirection)
!   enddo
!   write(20,*)
!enddo
      call solveZHEGV(valeDim,numStates,valeVale(:,:,spinDirection),&
            & valeValeOL(:,:),energyEigenValues(:,i,spinDirection))
#else
      call solveDSYGV(valeDim,numStates,valeValeGamma(:,:,spinDirection),&
            & valeValeOLGamma(:,:),energyEigenValues(:,i,spinDirection))
#endif

!! Test normalization.
!allocate (identity(numStates,numStates))
!! Read the atomic overlap matrix. 
!call readPackedMatrix(atomOverlap_did(i),packedValeVale,&
!      & packedVVDims,dim1,valeDim)
!
!      ! Unpack the overlap matrix.
!#ifndef GAMMA
!call unpackMatrix(valeValeOL(:,:,1,1),packedValeVale,valeDim,0)
!#else
!call unpackMatrixGamma(valeValeOLGamma(:,:,1),packedValeVale,valeDim,0)
!#endif
!
!identity(:numStates,:numStates) = &
!      & matmul(&
!      & matmul(transpose(conjg(valeVale(:valeDim,:numStates,1,spinDirection))),&
!             & valeValeOL(:,:,1,1)),valeVale(:,:numStates,1,spinDirection))
!write(24,*) "Normalization i = ",i
!do k = 1, num_states
!write(24,*) "k,I(k,k) = ",k,identity(k,k)
!enddo
!
!deallocate (identity)

      ! Write the energy eigenValues onto disk in HDF5 format in a.u.
      call h5dwrite_f (eVal_did(i,spinDirection),H5T_NATIVE_DOUBLE,&
            & energyEigenValues(:,i,spinDirection),states,hdferr)
      if (hdferr /= 0) stop 'Cannot write energy eigen values PSCF.'

      ! If we have some atoms with plusUJ terms, we need to update the just
      !   obtained eVectors. The update depends on the charge in each of the
      !   highest d or f orbitals of the affected atoms. Therefore, we will
      !   need to reread in the overlap matrix again because it was just
      !   destroyed in the ZHEGV solution. An alternative approach might be
      !   to make a copy of the overlap matrix and hold it in reserve until we
      !   need it for this calculation.
      if (numPlusUJAtoms > 0) then

         ! Read the atomic overlap matrix. 
         call readPackedMatrix(ol_did(i),packedValeVale,packedVVDimsPSCF,dim1,&
               & valeDim)

         ! Unpack the overlap matrix.
#ifndef GAMMA
         call unpackMatrix(valeValeOL(:,:),packedValeVale,valeDim,0)
#else
         call unpackMatrixGamma(valeValeOLGamma(:,:),packedValeVale,valeDim,0)
#endif
      endif

#ifndef GAMMA
      ! Write the eigenVectors onto disk in HDF5 format for this
      !   kpoint and spin direction.
      call h5dwrite_f(eVec_did(1,i,spinDirection),&
            & H5T_NATIVE_DOUBLE,real(valeVale(:,:numStates,spinDirection),&
            & double), valeStatesPSCF,hdferr)
      if (hdferr /= 0) stop 'Cannot write real energy eigen vectors PSCF.'
      call h5dwrite_f(eVec_did(2,i,spinDirection),&
            & H5T_NATIVE_DOUBLE,aimag(valeVale(:,:numStates,spinDirection)),&
            & valeStatesPSCF,hdferr)
      if (hdferr /= 0) stop 'Cannot write imag energy eigen vectors PSCF.'
#else
      call h5dwrite_f(eVec_did(1,i,spinDirection),&
            & H5T_NATIVE_DOUBLE,real(valeValeGamma(:,:numStates,&
            & spinDirection),double),valeStatesPSCF,hdferr)
      if (hdferr /= 0) stop 'Cannot write real energy eigen vectors PSCF.'
#endif

      ! Record that this kpoint has been finished.
      if (mod(i,10) .eq. 0) then
         write (20,ADVANCE="NO",FMT="(a1)") "|"
      else
         write (20,ADVANCE="NO",FMT="(a1)") "."
      endif
      if (mod(i,50) .eq. 0) then
         write (20,*) " ",i
      endif
      call flush (20)

      ! Record that this calculation is complete.
      call h5awrite_f(eVec_aid(i,spinDirection),H5T_NATIVE_INTEGER,&
            & 1,attribIntDims,hdferr)
      if (hdferr /= 0) stop 'Failed to record eigenvector success PSCF.'
      call h5aclose_f(eVec_aid(i,spinDirection),hdferr)
      if (hdferr /= 0) stop 'Failed to close eigenvector attribute PSCF.'
   enddo ! Loop i over kpoints.

   ! Deallocate unnecessary arrays and matrices.
   deallocate (packedValeVale)
   deallocate (tempPackedValeVale)
#ifndef GAMMA
   deallocate(valeVale)
   deallocate(valeValeOL)
#else
   deallocate(valeValeGamma)
   deallocate(valeValeOLGamma)
#endif

   ! Record the date and time that we finish.
   call timeStampEnd (15)

end subroutine secularEqnPSCF


! This subroutine will update the plusUJ term values on the basis of the charge
!   density matrix (obtained from the product of wave function coefficients
!   (psi* * psi) but before the overlap matrix elements have been multiplied
!   against the charge density matrix elements. That multiplication will be done
!   later (in update2AndApplyUJ) when the overlap matrix is available in an
!   unpacked form.
! The action of this algorithm will be to multiply each appropriate 5x5 or 7x7
!   block of the charge density matrix by the U coefficient
subroutine update1UJ (currKPoint, valeValeRho)

   ! Use necessary modules.
   use O_Kinds
   use O_Potential, only: numPlusUJAtoms, plusUJAtomSize, plusUJAtomValeIndex, &
         & plusUJ, spin

   ! Make sure that no unnecessary variables are declared.
   implicit none

   ! Declare the passed parameters.
   integer, intent(in) :: currKPoint
#ifndef GAMMA
   complex (kind=double), intent(inout), dimension (:,:,:) :: valeValeRho
#else
   real (kind=double),intent(inout), dimension (:,:,:) :: valeValeRho
#endif

   ! Declare local variables.
   integer :: i,j,k ! Loop index variables.
   integer :: currUJIndex ! The valence dimension number that is one before the
         ! first d or f orbital of whichever atom is being treated at the time.
   integer :: currUJSize ! The number of orbitals (5 for d, 7 for f) for
         ! whichever atom is being treated at the time.


   ! The essential physical goal that we need to perform for each atom with a
   !   plusUJ contribution follows equation (9) from Anisimov which is: V_ms =
   !   U*SUM_m'(n_m'(-s) - n^0) + (U-J)*SUM_m'/=m(n_m's - n^0) + V^LDA. Note
   !   that m = d orbital magnetic quantum number (1-5) or f orbital magnetic
   !   quantum number (1-7); m' is a dummy index number that runs over the
   !   magnetic quantum numbers; n_m's and n_m'(-s) are the occupation numbers
   !   of the indexed magnetic quantum number and spin states; n^0 is a
   !   reference point defined to assume even distribution of all d (f)
   !   electrons across d (f) orbitals; V^LDA is the already applied LDA
   !   potential; and s = one spin direction and (-s) is the other.
   ! We will compute a spin-up V_ms and a spin-down V_m(-s) that is constructed
   !   on the basis of the actual occupation of the spin-up and spin-down
   !   orbitals.
   ! See: Anisimov VI, Zaanen J, Andersen OK. Band theory and Mott insulators:
   !   Hubbard U instead of Stoner I. Physical Review B, 1991;44(3):943.
   !   Available from: http://dx.doi.org/10.1103/PhysRevB.44.943

   ! Start collecting the update to plusUJ from each atom with a UJ term.
   do i = 1, spin
      do j = 1, numPlusUJAtoms

         currUJIndex = plusUJAtomValeIndex(j)
         currUJSize = plusUJAtomSize(j)

         do k = 1, currUJSize ! Either 5 or 7
            ! A DIAGONAL element of a Hermitian density matrix is
            !   real by construction, so discarding the imaginary
            !   part is exact rather than approximate. Said with an
            !   explicit real() so a reader can see that is the
            !   claim being made, instead of inferring it from an
            !   implicit conversion.
            !
            ! The two arms of the old #ifndef GAMMA here were
            !   character-for-character identical, so the branch
            !   said nothing and is gone. In the real build
            !   valeValeRho is already real and real() is the
            !   identity.
            plusUJ(k,j,i,currKPoint) = &
                  & real(valeValeRho(currUJIndex+k,currUJIndex+k,i),&
                  & double)
         enddo
      enddo
   enddo

end subroutine update1UJ

! The big picture for this subroutine is that we need to modify specifc matrix
!   elements of the Hamiltonian matrix.
! Before we can do that though, we need to complete the construction of the
!   plusUJ terms. If this is the first iteration of an SCF cycle, then this
!   subroutine should not be called because the plusUJ terms are all zero. Only
!   on the second iteration will we have any knowledge about the plusUJ terms.
!   (Basically this is because I am assuming that the scfV.dat input file will
!   not have any initial values for the plusUJ terms and that the plusUJ terms
!   can only be created *after* we have the charge density matrix. The charge
!   density matrix can only be created from the single particle wave functions.)
! However, if this is part of a non-SCF calculation, then there are no
!   "iterations" and we need to complete the plusUJ right away.
subroutine update2AndApplyUJ(currKPoint,spinDirection)

   ! Use necessary modules.
   use O_Kinds
!   use O_KPoints, only: numKPoints, kPointWeight
   use O_KPoints, only: kPointWeight
   use O_Potential, only: plusUJ, plusUJAtomSize, numPlusUJAtoms, &
         & plusUJAtomValeIndex, plusUJAtomValue, plusUJAtomGSElectrons, spin

   ! Make sure that no funny variables are accidentally defined.
   implicit none

   ! Define passed parameters.
   integer, intent(in) :: currKPoint
   integer, intent(in) :: spinDirection

   ! Define local variables.
   integer :: i,j,k ! Loop index variables.
   integer :: currUJIndex ! The valence dimension number that is one before the
         ! first d or f orbital of whichever atom is being treated at the time.
   integer :: currUJSize ! The number of orbitals (5 for d, 7 for f) for
         ! whichever atom is being treated at the time.
   integer :: oppositeSpin ! If spinDirection==1 then this is 2 and vice-versa.
   real (kind=double) :: sum1, sum2 ! Two intermediate summations used in the
         ! computation of the plusUJ term for localized d- or f-orbitals in
         ! select atoms.

   ! An important thing to point out is that the plusUJ terms depend on *both*
   !   spin-up and spin-down charge densities. At this point we have access to
   !   both spin-up and spin-down charge densities in the plusUJ matrix from the
   !   update1 subroutine called in the previous SCF iteration (in the
   !   makeValenceRho subroutine). So, we *could* finalize the plusUJ term for
   !   both spin directions, but we will only do one (the current spinDirection)
   !   at a time. The reason is that at this point in the program, the valeVale
   !   matrix holds the wave function coefficients for only one spin. So, we
   !   can only apply the contributions to the up or down Hamiltonian.

   ! The essential physical goal that we need to perform for each atom with a
   !   plusUJ contribution follows equation (9) from Anisimov which is: V_ms =
   !   U*SUM_m'(n_m'(-s) - n^0) + (U-J)*SUM_m'/=m(n_m's - n^0) + V^LDA. Note
   !   that m = d orbital magnetic quantum number (1-5) or f orbital magnetic
   !   quantum number (1-7); m' is a dummy index number that runs over the
   !   magnetic quantum numbers; n_m's and n_m'(-s) are the occupation numbers
   !   of the indexed magnetic quantum number and spin states; n^0 is a
   !   reference point defined to assume even distribution of all d (f)
   !   electrons across d (f) orbitals; V^LDA is the already applied LDA
   !   potential; and s = one spin direction and (-s) is the other;
   ! We will compute a spin-up V_ms and a spin-down V_m(-s) that is constructed
   !   on the basis of the actual occupation of the spin-up and spin-down
   !   orbitals.
   ! See: Anisimov VI, Zaanen J, Andersen OK. Band theory and Mott insulators:
   !   Hubbard U instead of Stoner I. Physical Review B, 1991;44(3):943.
   !   Available from: http://dx.doi.org/10.1103/PhysRevB.44.943

   ! As another note, the division of plusUJAtomGSElectrons by (2.0*currUJSize)
   !   is designed to be a division by either 10 or 14 depending on whether
   !   this is a d or f orbital so that the value of n^0 will be (Actual number
   !   of electrons) / 10 or (Actual number of electrons) / 14. Of course, this
   !   quantity is also divided by the current kPoint weight to scale the number
   !   of electrons to the current kPoint.

   ! Establisht he value of the opposite spin direction.
   oppositeSpin = mod(spinDirection,2) + 1

   ! Update the plusUJ term and apply it to the Hamiltonian for each atom with
   !   some plusUJ contribution.
   do i = 1, numPlusUJAtoms

      ! Get the valence dimension index and size of the plusUJ term for the
      !   current atom.
      currUJIndex = plusUJAtomValeIndex(i)
      currUJSize = plusUJAtomSize(i)

      ! Update the plusUJ term by multiplying the charge density matrix by the
      !   overlap matrix elements to get the actual charge density. (Previously,
      !   the terms only held the products of wave function coefficients. This
      !   is close to the actual charge, but it lacks the overlap.)
      do j = 1, currUJSize

         ! Compute the first temporary sum for this orbital of this atom.
         sum1 = 0.0_double
         do k = 1, currUJSize
#ifndef GAMMA
            ! Diagonal of a Hermitian overlap: real by construction,
            !   so real() here is exact and states that claim.
            sum1 = sum1 + plusUJ(k,i,oppositeSpin,currKPoint) * &
                  & real(valeValeOL(currUJIndex+k,currUJIndex+k),&
                  & double) - &
                  & plusUJAtomGSElectrons(i) / (2.0_double * &
                  & real(currUJSize,double)) * kPointWeight(currKPoint) / &
                  & real(spin,double)
#else
            sum1 = sum1 + plusUJ(k,i,oppositeSpin,currKPoint) * &
                  & valeValeOLGamma(currUJIndex+k,currUJIndex+k) - &
                  & plusUJAtomGSElectrons(i) / (2.0_double * &
                  & real(currUJSize,double)) * kPointWeight(currKPoint) / &
                  & real(spin,double)
#endif
         enddo
         sum1 = sum1 * plusUJAtomValue(1,i)

         ! Compute the second temporary sum for this orbital of this atom.
         sum2 = 0.0_double
         do k = 1, currUJSize

            ! Don't include m=m' terms in the summation.
            if (k == j) then
               cycle
            endif

            ! Accumulate summation terms.
#ifndef GAMMA
            ! Diagonal of a Hermitian overlap, as above.
            sum2 = sum2 + plusUJ(k,i,spinDirection,currKPoint) * &
                  & real(valeValeOL(currUJIndex+k,currUJIndex+k),&
                  & double) - &
                  & plusUJAtomGSElectrons(i) / (2.0_double * &
                  & real(currUJSize,double)) * kPointWeight(currKPoint) / &
                  & real(spin,double)
#else
            sum2 = sum2 + plusUJ(k,i,spinDirection,currKPoint) * &
                  & valeValeOLGamma(currUJIndex+k,currUJIndex+k) - &
                  & plusUJAtomGSElectrons(i) / (2.0_double * &
                  & real(currUJSize,double)) * kPointWeight(currKPoint) / &
                  & real(spin,double)
#endif
         enddo
         sum2 = sum2 * (plusUJAtomValue(1,i) - plusUJAtomValue(2,i))

      ! Now that the plusUJ term has been fully updated (for the current
      !   spinDirection) we will apply it to the Hamiltonian with the current
      !   spinDirection.
#ifndef GAMMA
         valeVale(currUJIndex+j,currUJIndex+j,spinDirection) = &
               & valeVale(currUJIndex+j,currUJIndex+j,spinDirection) + &
               & sum1 + sum2
#else
         valeValeGamma(currUJIndex+j,currUJIndex+j,spinDirection) = &
               & valeValeGamma(currUJIndex+j,currUJIndex+j,spinDirection) + &
               & sum1 + sum2
#endif
      enddo
   enddo

end subroutine update2AndApplyUJ


subroutine shiftEnergyEigenValues(energyShift)

!use O_KPoints
   ! Make sure that no funny variables are defined.
   implicit none

   ! Define dummy variables passed to this subroutine.
   real (kind=double) :: energyShift
!integer :: i,j,k

   ! Shift the energyEigenValues down by the requested about.
   energyEigenValues(:,:,:) = energyEigenValues(:,:,:) - energyShift

!write(20,*) "energyShift:  ", energyShift
!do i = 1, 1
!do j = 1, numKPoints
!do k = 1, numStates
!write(20,*) i,j,k, energyEigenValues(k,j,i)
!enddo
!enddo
!enddo

end subroutine shiftEnergyEigenValues


!subroutine readEnergyEigenValuesBand(numStates)
!
!   ! Use necessary modules
!   use HDF5
!   use O_Kinds
!   use O_KPoints, only: numKPoints
!   use O_Potential, only: spin
!   use O_PSCFBandHDF5, only: statesBand, eigenValuesBand_did
!
!   ! Make sure that no funny variables are defined.
!   implicit none
!
!   ! define variables passed to this subroutine.
!   integer :: numStates
!
!   ! Define local variables used in this subroutine.
!   integer :: i,j
!   integer :: hdferr
!   real (kind=double), allocatable, dimension (:)   :: energyValuesTemp
!
!   ! Allocate space for a reading buffer.
!   allocate (energyValuesTemp (numStates))
!
!   ! Read the ground state energy eigen values.
!
!   ! Loop over each kpoint and spin direction to read the energy values.
!   do i = 1, numKPoints
!      do j = 1, spin
!         call h5dread_f (eigenValuesBand_did(i,j),H5T_NATIVE_DOUBLE,&
!               & energyValuesTemp(:),statesBand,hdferr)
!         if (hdferr /= 0) stop 'Failed to read energy eigen values'
!
!         ! Copy the necessary values into the final array.
!         energyEigenValues(:,i,j) = energyValuesTemp(:)
!      enddo
!   enddo
!
!   ! Deallocate the reading buffer
!   deallocate (energyValuesTemp)
!
!end subroutine readEnergyEigenValuesBand
!
!
!subroutine appendExcitedEValsBand (firstStateIndex,numStates)
!
!   ! Use necessary modules
!   use HDF5
!   use O_Kinds
!   use O_KPoints, only: numKPoints
!   use O_Potential, only: spin
!   use O_PSCFBandHDF5, only: statesBand, eigenValuesBand2_did
!
!   ! Make sure that no funny variables are defined.
!   implicit none
!
!   ! Define variables passed to this subroutine.
!   integer :: firstStateIndex
!   integer :: numStates
!
!   ! Define local variables used in this subroutine.
!   integer :: i,j
!   integer :: hdferr
!   real (kind=double), allocatable, dimension (:)   :: energyValuesTemp
!
!   ! Allocate space for a reading buffer.
!   allocate (energyValuesTemp (numStates))
!
!   ! It is assumed that we know the highest occupied state for each kpoint.
!   !   We now read the excited state energy eigen values and merge them into
!   !   the ground state ones with the merge points being the highest occupied
!   !   states.  We must be careful when dealing with degenerate states.
!   do i = 1, numKPoints
!      do j = 1, spin
!
!         ! Get the excited state energy values.
!         call h5dread_f (eigenValuesBand2_did(i,j),H5T_NATIVE_DOUBLE,&
!               & energyValuesTemp(:),statesBand,hdferr)
!         if (hdferr /= 0) stop 'Failed to read excited energy eigen values'
!
!         energyEigenValues(firstStateIndex:numStates,i,j) = &
!               & energyValuesTemp(firstStateIndex:numStates)
!
!      enddo
!   enddo
!
!   ! Deallocate the reading buffer
!   deallocate (energyValuesTemp)
!
!end subroutine appendExcitedEValsBand


!subroutine preserveValeValeOL
!
!   ! Use necessary modules.
!   use O_Potential, only: spin, numPlusUJAtoms
!   use O_AtomicSites, only: valeDim
!
!   ! Make sure that no funny variables are defined.
!   implicit none
!
!   ! In some cases (e.g. during a spin polarized calculation) it is necessary
!   !   to preserve the valeValeOL during the diagonalization.  This is
!   !   because the LAPACK routine will destroy the overlap matrix because it
!   !   is the same for spin up or spin down (only one copy is made).  We will
!   !   copy the valeValeOL (spin index 1) into the spin index 2 part of the
!   !   same matrix (valeValeOL).
!   if ((spin == 2) .or. (numPlusUJAtoms > 0)) then
!#ifndef GAMMA
!      valeValeOL(:valeDim,:valeDim,1,2) = valeValeOL(:valeDim,:valeDim)
!#else
!      valeValeOLGamma(:valeDim,:valeDim,2) = &
!            & valeValeOLGamma(:valeDim,:valeDim,1)
!#endif
!   endif
!
!   ! IMPORTANT NOTE:  This whole thing could be "improved" in terms of memory
!   !   usage if necessary.  Instead of copying to the spin index 2 of the
!   !   valeValeOL we could copy to spin index 2 of the valeVale because it
!   !   has not yet been used.  Along with this, the valeValeOL would only need
!   !   to be allocated 1 spin index worth of memory.  This could be important,
!   !   but it will not be pursued now because it is an ugly and potentially
!   !   more confusing thing to do.  (The next subroutine would also need to be
!   !   changed accordingly.)
!
!end subroutine preserveValeValeOL
!
!
!subroutine restoreValeValeOL
!
!   ! Use necessary modules.
!   use O_Potential, only: spin
!   use O_AtomicSites, only: valeDim
!
!   ! Make sure that no funny variables are defined.
!   implicit none
!
!   ! As with the above subroutine, the valeValeOL matrix must sometimes be
!   !   saved from destruction during the diagonalization.  There is no spin
!   !   aspect to the overlap matrix and so the matrix was copied from the
!   !   spin index 1 to the spin index 2 for safe keeping.  This subroutine
!   !   will restore the matrix from spin index 2 to spin index 1.
!   if (spin == 2) then
!#ifndef GAMMA
!      valeValeOL(:valeDim,:valeDim) = valeValeOL(:valeDim,:valeDim,1,2)
!#else
!      valeValeOLGamma(:valeDim,:valeDim) = &
!            & valeValeOLGamma(:valeDim,:valeDim,2)
!#endif
!   endif
!end subroutine restoreValeValeOL



subroutine readDataSCF(h,i,numStates,matrixCode,slab)

   ! Deliver one k-point's eigenvectors, and the integral matrices the
   !   caller names by matrixCode, from the SCF data structures.  The
   !   contract for the eigenvectors (DESIGN 2.8, PSEUDOCODE 33): after
   !   this call the eigenvectors of spin h at k-point i occupy slab
   !   `slab` of valeVale (Gamma build: valeValeGamma), whatever the
   !   build and whichever routine allocated that array.  The slab is
   !   the CALLER's to name because only the caller knows the shape it
   !   allocated: a consumer that holds one slab per spin (the valence
   !   density, dimo, field, mtop) omits the argument and receives slab
   !   h; a consumer that processes the spins one at a time and holds a
   !   single slab (the DOS, bond order and optical programs) passes
   !   slab = 1.  Before this argument existed those consumers read
   !   slab 1 while the eigenvectors of spin two were delivered to slab
   !   2, so every spin-polarized result for spin two was wrong or an
   !   out-of-bounds write (DEBUG.md BUG-028).

   ! Import necessary data modules.
   use O_AtomicSites, only: valeDim
   use O_SCFIntegralsHDF5, only: packedVVDims,atomOverlap_did,&
         & atomMMOverlap_did
#ifndef GAMMA
   ! The k-point overlap datasets (matrix codes 3 through 8) are full,
   !   unpacked matrices read only on the multi-k path; the gamma build
   !   reads nothing but the packed datasets above, so these names and
   !   the full-matrix dimensions exist only in the complex build.
   use O_SCFIntegralsHDF5, only: fullVVDims, atomKOverlap_did, &
         & atomKOverlapPlusG_did
   use O_SCFEigVecHDF5, only: valeStates,eigenVectors_did
   use O_MatrixSubs, only: readMatrix,readPackedMatrix,unpackMatrix
#else
   use O_MatrixSubs, only: readPackedMatrix,unpackMatrixGamma
#endif

   ! Define passed parameters.  In the multi-k build h and numStates
   !   drive the eigenvector re-read from the file; in the gamma build
   !   the single k-point's eigenvectors stay resident in the solver's
   !   array, and the two arguments drive only the copy that delivers
   !   them to the requested slab.
   integer, intent(in) :: h ! Spin variable.
   integer, intent(in) :: i ! KPoint variable
   integer, intent(in) :: numStates
   integer, intent(in) :: matrixCode
   integer, intent(in), optional :: slab ! Destination slab; default h.

   ! Define local variables.
   integer :: destSlab ! Where the eigenvectors of spin h are delivered.
   integer :: dim1
   integer :: j ! Loop index (usually xyz).
   real (kind=double), allocatable, dimension (:,:) :: packedValeVale
#ifndef GAMMA
   real (kind=double), allocatable, dimension (:,:) :: tempRealValeVale
   real (kind=double), allocatable, dimension (:,:) :: tempImagValeVale
#endif

   ! Initialize the dimension of the packed matrices to include two components
   !   (real,imaginary) or just a real component.
#ifndef GAMMA
   dim1 = 2
#else
   dim1 = 1
#endif

   if (matrixCode > 0) then

      ! Allocate space to read a packed matrix.
      allocate (packedValeVale(dim1,valeDim*(valeDim+1)/2))

      if (matrixCode == 1) then
         ! Read the overlap matrix.  The tempPackedMatrix is not used.
         call readPackedMatrix(atomOverlap_did(i),packedValeVale,&
               & packedVVDims,dim1,valeDim)

         ! Unpack the matrix.
#ifndef GAMMA
         call unpackMatrix(valeValeOL(:,:),packedValeVale,valeDim,1)
#else
         call unpackMatrixGamma(valeValeOLGamma(:,:),packedValeVale,valeDim,1)
#endif
      elseif (matrixCode == 2) then
         do j = 1, 3
            ! Read the xyz momentum matrix elements.
            call readPackedMatrix(atomMMOverlap_did(i,j),packedValeVale,&
                  & packedVVDims,dim1,valeDim)

            ! Unpack the matrix.
#ifndef GAMMA
           call unpackMatrix(valeValeMM(:,:,j),packedValeVale,valeDim,1)
#else
           call unpackMatrixGamma(valeValeMMGamma(:,:,j),packedValeVale,valeDim,1)
#endif
         enddo
      elseif ((matrixCode >= 3) .and. (matrixCode <= 5)) then

#ifndef GAMMA
         ! Allocate space to read the complex KOverlap matrix
         allocate (tempRealValeVale (valeDim,valeDim))
         allocate (tempImagValeVale (valeDim,valeDim))

         ! Read the complex KOverlap matrix from the datasets.
         call readMatrix(atomKOverlap_did(i,matrixCode-2,:),&
            & valeValeKO(:,:,matrixCode-2),&
            & tempRealValeVale(:,:),tempImagValeVale(:,:),&
            & fullVVDims,valeDim,valeDim)

         ! Deallocate the space to read the complex wave function.
         deallocate (tempRealValeVale)
         deallocate (tempImagValeVale)
#endif
      elseif ((matrixCode >= 6) .and. (matrixCode <= 8)) then

#ifndef GAMMA
         ! Allocate space to read the complex KOverlapPlusG matrix
         allocate (tempRealValeVale (valeDim,valeDim))
         allocate (tempImagValeVale (valeDim,valeDim))

         ! Read the complex KOverlap matrix from the datasets.
         call readMatrix(atomKOverlapPlusG_did(i,matrixCode-5,:),&
            & valeValeKO(:,:,matrixCode-5),&
            & tempRealValeVale(:,:),tempImagValeVale(:,:),&
            & fullVVDims,valeDim,valeDim)

         ! Deallocate the space to read the complex wave function.
         deallocate (tempRealValeVale)
         deallocate (tempImagValeVale)
#endif
      endif

      ! Deallocate the packed matrix used to read and unpack the data.
      deallocate (packedValeVale)

   endif

   ! The destination slab for the eigenvectors: the caller's choice, or
   !   the spin index when the caller holds one slab per spin.
   destSlab = h
   if (present(slab)) destSlab = slab

#ifndef GAMMA
   ! Read the wave functions for this kpoint from the datasets into
   !   the requested slab of the valeVale matrix.
!   if (numKPoints > 1) then

      ! Allocate space to read the wave functions.
      allocate (tempRealValeVale (valeDim,numStates))
      allocate (tempImagValeVale (valeDim,numStates))

      call readMatrix(eigenVectors_did(:,i,h),&
            & valeVale(:,:numStates,destSlab),&
            & tempRealValeVale(:,:),tempImagValeVale(:,:),&
            & valeStates,valeDim,numStates)

      ! Deallocate the space to read the wave functions.
      deallocate (tempRealValeVale)
      deallocate (tempImagValeVale)
!   endif

#else
   ! The gamma build never writes its single k-point's eigenvectors to
   !   the file during the SCF: they stay in slab h of the solver's
   !   valeValeGamma(valeDim,valeDim,spin) from the last solve.  To
   !   deliver them to another slab is therefore a copy, not a read.
   !   It is a valeDim by numStates copy, once per spin per consumer,
   !   and it happens only for a one-slab consumer on spin two.
   if (destSlab /= h) then
      valeValeGamma(:,:numStates,destSlab) = valeValeGamma(:,:numStates,h)
   endif
#endif

end subroutine readDataSCF



subroutine readDataPSCF(h,i,numStates,matrixCode,slab)

   ! The post-SCF twin of readDataSCF, with the same eigenvector
   !   contract (DESIGN 2.8, PSEUDOCODE 33): the eigenvectors of spin h
   !   at k-point i are delivered to slab `slab` of valeVale (Gamma
   !   build: valeValeGamma), default slab h.  See readDataSCF for why
   !   the caller names the slab.  Here both builds read the vectors
   !   from the post-SCF file, so the slab is simply the read target.

   ! Use necessary modules.
!   use O_KPoints, only: numKPoints
   use O_AtomicSites, only: valeDim
   use O_PSCFIntegralsHDF5, only: packedVVDimsPSCF,&
         & atomOverlapPSCF_did,atomMMOverlapPSCF_did
   use O_PSCFEigVecHDF5, only: valeStatesPSCF,eigenVectorsPSCF_did
#ifndef GAMMA
   ! The k-point overlap datasets (matrix codes 3 through 8) are full,
   !   unpacked matrices read only on the multi-k path; the gamma build
   !   reads nothing but the packed datasets above, so these names and
   !   the full-matrix dimensions exist only in the complex build.
   use O_PSCFIntegralsHDF5, only: fullVVDimsPSCF, atomKOverlapPSCF_did, &
         & atomKOverlapPlusGPSCF_did
   use O_MatrixSubs, only: readMatrix, readPackedMatrix, unpackMatrix
#else
   use O_MatrixSubs, only: readMatrixGamma, readPackedMatrix, &
         & unpackMatrixGamma
#endif

   ! Define passed parameters.
   integer, intent(in) :: h ! Spin variable.
   integer, intent(in) :: i ! KPoint variable
   integer, intent(in) :: numStates
   integer, intent(in) :: matrixCode
   integer, intent(in), optional :: slab ! Destination slab; default h.

   ! Define local variables.
   integer :: destSlab ! Where the eigenvectors of spin h are delivered.
   integer :: dim1
   integer :: j ! Loop index (usually xyz).
   real (kind=double), allocatable, dimension (:,:) :: packedValeVale
#ifndef GAMMA
   real (kind=double), allocatable, dimension (:,:) :: tempRealValeVale
   real (kind=double), allocatable, dimension (:,:) :: tempImagValeVale
#endif

   ! The destination slab for the eigenvectors: the caller's choice, or
   !   the spin index when the caller holds one slab per spin.
   destSlab = h
   if (present(slab)) destSlab = slab

   ! Initialize the dimension of the packed matrices to include two components
   !   (real,imaginary) or just a real component.
#ifndef GAMMA
   dim1 = 2
#else
   dim1 = 1
#endif

   if (matrixCode > 0) then

      ! Allocate space to read a packed matrix.
      allocate (packedValeVale(dim1,valeDim*(valeDim+1)/2))

      if (matrixCode == 1) then
         ! Read the overlap matrix.  The tempPackedMatrix is not used.
         call readPackedMatrix(atomOverlapPSCF_did(i),packedValeVale,&
               & packedVVDimsPSCF,dim1,valeDim)

         ! Unpack the matrix.
#ifndef GAMMA
         call unpackMatrix(valeValeOL(:,:),packedValeVale,valeDim,1)
#else
         call unpackMatrixGamma(valeValeOLGamma(:,:),packedValeVale,valeDim,1)
#endif

      elseif (matrixCode == 2) then
         do j = 1, 3
            ! Read the xyz momentum matrix elements.
            call readPackedMatrix(atomMMOverlapPSCF_did(i,j),packedValeVale,&
                  & packedVVDimsPSCF,dim1,valeDim)

            ! Unpack the matrix.
#ifndef GAMMA
            call unpackMatrix(valeValeMM(:,:,j),packedValeVale,valeDim,1)
#else
            call unpackMatrixGamma(valeValeMMGamma(:,:,j),packedValeVale,valeDim,1)
#endif
         enddo
      elseif ((matrixCode >= 3) .and. (matrixCode <= 5)) then

#ifndef GAMMA
         ! Allocate space to read the complex KOverlap matrix
         allocate (tempRealValeVale (valeDim,valeDim))
         allocate (tempImagValeVale (valeDim,valeDim))

         ! Read the complex KOverlap matrix from the datasets.
         call readMatrix(atomKOverlapPSCF_did(i,matrixCode-2,:),&
            & valeValeKO(:,:,matrixCode-2),&
            & tempRealValeVale(:,:),tempImagValeVale(:,:),&
            & fullVVDimsPSCF,valeDim,valeDim)

         ! Deallocate the space to read the complex wave function.
         deallocate (tempRealValeVale)
         deallocate (tempImagValeVale)
#endif

      elseif ((matrixCode >= 6) .and. (matrixCode <= 8)) then

#ifndef GAMMA
         ! Allocate space to read the complex KOverlapPlusG matrix
         allocate (tempRealValeVale (valeDim,valeDim))
         allocate (tempImagValeVale (valeDim,valeDim))

         ! Read the complex KOverlap matrix from the datasets.
         call readMatrix(atomKOverlapPlusGPSCF_did(i,matrixCode-5,:),&
            & valeValeKO(:,:,matrixCode-5),&
            & tempRealValeVale(:,:),tempImagValeVale(:,:),&
            & fullVVDimsPSCF,valeDim,valeDim)

         ! Deallocate the space to read the complex wave function.
         deallocate (tempRealValeVale)
         deallocate (tempImagValeVale)
#endif

      endif

      ! Deallocate the packed matrix used to read and unpack the data.
      deallocate (packedValeVale)
   endif

#ifndef GAMMA
   ! Read the wave functions for this kpoint from the datasets into
   !   the valeVale matrix.  If numKPoints==1, the wave functions should
   !   already be in the valeVale(1:valeDim,1:numStates) matrix, unless it
   !   was already computed once in which case we just have not yet read it
   !   in (same for the eigen values).
!   if (numKPoints > 1) then

      ! Allocate space to read the complex wave function.
      allocate (tempRealValeVale (valeDim,numStates))
      allocate (tempImagValeVale (valeDim,numStates))

      ! Read the complex wave function from the datasets into the
      !   requested slab.
      call readMatrix(eigenVectorsPSCF_did(:,i,h),&
         & valeVale(:,:numStates,destSlab),&
         & tempRealValeVale(:,:),tempImagValeVale(:,:),&
         & valeStatesPSCF,valeDim,numStates)

      ! Deallocate the space to read the complex wave function.
      deallocate (tempRealValeVale)
      deallocate (tempImagValeVale)
!   endif
#else
   ! Read the real wave function from the datasets into the requested
   !   slab.
   call readMatrixGamma(eigenVectorsPSCF_did(1,i,h),&
         & valeValeGamma(:,:numStates,destSlab),valeStatesPSCF,valeDim,&
         & numStates)
#endif

end subroutine readDataPSCF


subroutine cleanUpSecularEqn

   ! Make sure that no funny variables are defined.
   implicit none

   ! Guard each deallocation because this subroutine
   !   is called from multiple program paths (SCF
   !   cleanup, PSCF cleanup) where the allocation
   !   state depends on which code path ran. For
   !   example, secularEqnPSCF deallocates valeVale
   !   after writing eigenvectors to HDF5, but the
   !   LAT TDOS path never re-allocates it (unlike
   !   the Gaussian PDOS path which does).
   if (allocated(energyEigenValues)) then
      deallocate (energyEigenValues)
   endif
#ifndef GAMMA
   if (allocated(valeVale)) then
      deallocate (valeVale)
   endif
#else
   if (allocated(valeValeGamma)) then
      deallocate (valeValeGamma)
   endif
#endif

end subroutine cleanUpSecularEqn

end module O_SecularEquation

!! SPDX-License-Identifier: ECL-2.0
!! Copyright (c) 2026 Paul Rulis

module O_ValeCharge

   ! Import necessary modules.
   use O_Kinds
   use O_Constants

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Charge fitting structures.  The second index for potRho and the only
   !   index for the others is reserved for spin.
   real (kind=double), allocatable, dimension (:,:) :: potRho
   real (kind=double), allocatable, dimension (:)   :: chargeDensityTrace
   real (kind=double), allocatable, dimension (:)   :: nucPotTrace
   real (kind=double), allocatable, dimension (:)   :: kineticEnergyTrace
   real (kind=double), allocatable, dimension (:)   :: massVelocityTrace
   real (kind=double), allocatable, dimension (:,:) :: dipoleMomentTrace

   real (kind=double), allocatable, dimension (:,:)   :: packedValeVale
   real (kind=double), allocatable, dimension (:,:,:) :: packedValeValeRho


   ! Split measurement of the valence charge density stage (PSEUDOCODE
   !   30).  Building the density matrix from the eigenvectors and
   !   contracting it against the stored integral matrices are two kinds
   !   of work with two different limits: the first is arithmetic over a
   !   resident matrix, the second is dominated by reading matrices back
   !   from the HDF5 file.  How the stage divides between them, and
   !   whether the building half is limited by arithmetic or by memory
   !   bandwidth, decides which parallel decomposition is worth building
   !   (DESIGN 9.6).  Nothing in the program reads these values; they
   !   exist to be logged once per call and read by a person.

   ! An integer kind wide enough to hold a clock tick count for the
   !   duration of a run.  Eighteen decimal digits is far more than a
   !   nanosecond-resolution counter needs, so no tick rollover handling
   !   is required.  Declared here rather than in the shared O_Kinds,
   !   which every subprogram in the project compiles against and which
   !   should not change to serve one measurement.
   integer, parameter :: clockIntKind = selected_int_kind(18)

   ! Ticks per second, as reported by system_clock at the top of each
   !   call.  A value of zero means the system has no usable clock, in
   !   which case the seconds below stay zero and the log line says so
   !   rather than dividing by it.
   integer (kind=clockIntKind) :: valeRhoClockRate

   ! Wall-clock seconds accumulated in each of the four regions, summed
   !   over the k-points of one call.  Wall clock, not processor time:
   !   the measurement these serve varies the BLAS thread count and
   !   watches the duration respond, and a processor-time total would
   !   grow with the thread count by construction and report perfect
   !   insensitivity no matter what the code did.
   real (kind=double) :: valeRhoReadVectorSeconds   ! Eigenvector reads
   real (kind=double) :: valeRhoAccumulateSeconds   ! Build the density
   real (kind=double) :: valeRhoReadIntegralSeconds ! Integral reads
   real (kind=double) :: valeRhoContractSeconds     ! Contract arithmetic

   ! Rank-1 update calls actually issued, and k-points not skipped.  The
   !   first turns a load-bearing but assumed quantity into a measured
   !   one: the traffic estimate that orders the candidate
   !   decompositions multiplies the size of the density matrix by the
   !   number of states whose occupation clears the negligibility
   !   threshold, and that state count has never been counted.  The
   !   second records how many k-points the skip test let through, which
   !   silently sets the denominator of every per-k-point figure derived
   !   from the log.
   integer :: valeRhoRankUpdateCount
   integer :: valeRhoKPointsProcessed

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains


! Read the wall clock into a caller-held marker, opening a timed region.
!   The marker is a local of the calling routine rather than module data
!   so that the regions cannot accidentally share one start time.
subroutine beginTimedRegion (regionStartTicks)

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define passed parameters.
   integer (kind=clockIntKind), intent(out) :: regionStartTicks

   call system_clock (count=regionStartTicks)

end subroutine beginTimedRegion


! Close a timed region opened by beginTimedRegion and add its duration to
!   a running total.  Regions do not nest, but one region may be entered
!   and left many times per call (the integral reads and the contractions
!   alternate, once per stored matrix), which is why the total
!   accumulates here instead of being assigned.
subroutine endTimedRegion (regionStartTicks,accumulatedSeconds)

   ! Import the precision variables.
   use O_Kinds

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define passed parameters.
   integer (kind=clockIntKind), intent(in) :: regionStartTicks
   real (kind=double), intent(inout) :: accumulatedSeconds

   ! Define the local variables.
   integer (kind=clockIntKind) :: regionEndTicks

   ! Without a usable clock there is nothing to add.  The counters in
   !   this module remain meaningful and are still reported.
   if (valeRhoClockRate <= 0) then
      return
   endif

   call system_clock (count=regionEndTicks)

   accumulatedSeconds = accumulatedSeconds + &
         & real(regionEndTicks - regionStartTicks,double) / &
         & real(valeRhoClockRate,double)

end subroutine endTimedRegion

subroutine makeValenceRho(inSCF)

   ! Import the necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_CommandLine, only: doForce_SCF, doForce_PSCF
   use O_AtomicSites, only: valeDim
   use O_Input, only: numStates
   use O_KPoints, only: numKPoints
   use O_Constants, only: smallThresh
   use O_Potential, only: rel,spin,potDim,potCoeffs,&
         & numPlusUJAtoms, converged
   use O_Populate, only: electronPopulation,cleanUpPopulation, &
         & electronPopulation_LAT
   use O_KPoints, only: kPointIntgCode
   use O_SCFIntegralsHDF5, only: atomOverlap_did,atomKEOverlap_did, &
         & atomMVOverlap_did,atomNPOverlap_did,atomPotOverlap_did,packedVVDims
   use O_PSCFIntegralsHDF5, only: atomOverlapPSCF_did,packedVVDimsPSCF
#ifndef GAMMA
   use O_BLASZHER
   use O_SecularEquation, only: valeVale,cleanUpSecularEqn,energyEigenValues,&
         & update1UJ, readDataSCF, readDataPSCF
   use O_MatrixSubs, only: readPackedMatrix,matrixElementMult,packMatrix
   use O_Force, only: computeForce
#else
   use O_BLASDSYR
   use O_SecularEquation, only: valeValeGamma, cleanUpSecularEqn, &
         & energyEigenValues, update1UJ, readDataSCF, readDataPSCF
   use O_MatrixSubs, only: readPackedMatrix, &
         & matrixElementMultGamma,packMatrixGamma
   use O_Force, only: computeForceGamma
#endif

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define passed parameters.
   integer, intent(in) :: inSCF

   ! Define the local variables used in this subroutine.
   integer :: i,j,k ! Loop index variables
! The l index belongs to the commented-out force-matrix symmetrization
!   code in the force block near the end of this subroutine; restore
!   this declaration together with that code.
!   integer :: l
#ifndef GAMMA
   ! These two exist only for the multi-k eigenvector re-read below:
   !   valeVale holds one k-point at a time, so each k-point's wave
   !   functions must be read back from disk (h indexes the spin of
   !   that read), and k-points with negligible occupation are skipped
   !   (skipKP). The gamma build keeps its single k-point resident in
   !   memory, so it never re-reads. Spin polarization itself is NOT
   !   variant-dependent: the do-k=1,spin loops below are shared.
   integer :: h
   integer :: skipKP
#endif
   integer :: dim1
   integer :: energyLevelCounter
   real (kind=double) :: sumElecEnergy
   real (kind=double), allocatable, dimension (:)     :: tempDensity
   real (kind=double), allocatable, dimension (:)     :: electronEnergy
   real (kind=double), allocatable, dimension (:)     :: currentPopulation
   real (kind=double), allocatable, dimension (:,:,:) :: &
         & structuredElectronPopulation
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:) :: valeValeRho
#else
   real    (kind=double), allocatable, dimension (:,:,:) :: valeValeRhoGamma
#endif

   ! Marker holding the tick count at which the currently open timed
   !   region began (PSEUDOCODE 30).  The four regions run one after
   !   another and never nest, so a single marker serves all of them.
   integer (kind=clockIntKind) :: clockAtRegionStart

   ! Log the date and time that we start.
   call timeStampStart (17)

   ! Begin the split measurement of this call (PSEUDOCODE 30).  The
   !   totals are cleared here rather than at module load so that every
   !   call reports its own iteration instead of a running sum over the
   !   whole self-consistency loop.  The clock rate is asked for once,
   !   because it does not change during a run.
   call system_clock (count_rate=valeRhoClockRate)
   valeRhoReadVectorSeconds   = 0.0_double
   valeRhoAccumulateSeconds   = 0.0_double
   valeRhoReadIntegralSeconds = 0.0_double
   valeRhoContractSeconds     = 0.0_double
   valeRhoRankUpdateCount     = 0
   valeRhoKPointsProcessed    = 0

   ! Define whether the packed arrays have two rows (complex) or one (real).
   !   This also defines the size of some other small arrays.
#ifndef GAMMA
   dim1 = 2
#else
   dim1 = 1
#endif

   ! Allocate the main valence valence density matrix.
#ifndef GAMMA
   allocate (valeValeRho(valeDim,valeDim,spin)) ! Complex
#else
   allocate (valeValeRhoGamma(valeDim,valeDim,spin)) ! Real
#endif

   ! Allocate a temporary packed valeVale matrix with a spin component.
   allocate (packedValeValeRho(dim1,valeDim*(valeDim+1)/2,spin))

   ! Allocate a temporary holder for the charge density for the case that the
   !   calculation is spin polarized.  This is needed when rewriting the charge
   !   in a up+down, up-down form.
   allocate (tempDensity(dim1))

   ! Allocate space for the valence charge density (as represented by a
   !   summation of atom centered Gaussian functions in the same way as the
   !   potential function except with different coefficients).
   allocate (potRho (potDim,spin)) ! This will be deallocated in the makeSCFPot
         ! subroutine since after its values are copied to a local array there
         ! it will not be needed again until here.

   ! Allocate space to hold the trace of the charge density, nuclear potential,
   !   kinetic energy, mass velocity, and dipole moment.
   if (inSCF == 1) then
      allocate (nucPotTrace(spin))
      nucPotTrace(:) = 0.0_double
      allocate (chargeDensityTrace(spin))
      chargeDensityTrace(:) = 0.0_double
      allocate (kineticEnergyTrace(spin))
      kineticEnergyTrace(:) = 0.0_double
      if (rel == 1) then
         allocate (massVelocityTrace(spin))
         massVelocityTrace(:) = 0.0_double
      endif
   else
#ifndef GAMMA
      allocate (valeVale(valeDim,numStates,spin))
      valeVale(:,:,:) = cmplx(0.0_double,0.0_double,double)
#else
      allocate (valeValeGamma(valeDim,numStates,spin))
      valeValeGamma(:,:,:) = 0.0_double
#endif
   endif

   ! Allocate space to hold the currentPopulation based on spin
   allocate (currentPopulation (spin))
   allocate (electronEnergy (spin))
   allocate (structuredElectronPopulation (numStates,numKPoints,spin))

   ! Initialize local variables
   electronEnergy(:) = 0.0_double
   potRho(:,:) = 0.0_double


   ! Fill a matrix of electron populations from whichever integration
   !   method produced them (DESIGN 1.6c).
   !
   ! Under LAT the occupations already have the (state, kpoint, spin)
   !   shape this matrix wants, so the LAT branch is a straight copy and
   !   the unpack below is skipped entirely rather than being fed a
   !   differently-ordered array.  The one conversion needed is of
   !   CONVENTION: electronPopulation carries the kPointWeight
   !   convention, whose weights sum to 2.0 so that a non-polarized
   !   calculation holds two electrons per state, while
   !   electronPopulation_LAT holds pure Brillouin-zone volume fractions
   !   summing to 1.0 per occupied band per spin.  The factor 2/spin
   !   converts between them, the same conversion computeBond applies at
   !   its own point of use (DESIGN 1.6d).
   if (kPointIntgCode == 1) then
      structuredElectronPopulation(:,:,:) = &
            & electronPopulation_LAT(:,:,:) &
            & * 2.0_double / real(spin, double)
   else

      ! Note that electronPopulation is a one dimensional array that has
      !   some order, but is not sorted in the way that the energy eigen
      !   values were sorted.  Please read the comments in the
      !   populateLevels subroutine to understand the order.
      !   (You can also probably get it from the loop order here ;)
      energyLevelCounter=0
      do i = 1, numKPoints
         do j = 1, spin
            do k = 1, numStates
               energyLevelCounter = energyLevelCounter + 1
               structuredElectronPopulation (k,i,j) = &
                     & electronPopulation(energyLevelCounter)
            enddo
         enddo
      enddo
   endif

   do i = 1, numKPoints

      ! Initialize space to read the wave functions.  Note that this matrix is
      !   a double complex matrix because the wave function data contains both
      !   real and imaginary parts and has no symmetry that might permit
      !   packing or use of triangular form (if it were a Hermitian matrix
      !   instead).
      ! Note that it is only necessary to go through the initialization action
      !   if there are more than 1 kpoint.  For the 1 kpoint case, the valeVale
      !   matrix was not changed so it can still be used here.  Also note
      !   that the gammaKPoint option will never enter the "if" block.
#ifndef GAMMA
!      if (numKPoints > 1) then

         ! Skip any kpoints with a negligable contribution for each state.
         skipKP = 1 ! Assume that we will skip this kpoint.
         do j = 1, numStates
            if (sum(abs(structuredElectronPopulation(j,i,:)))>smallThresh) then
               skipKP = 0 ! Enough contribution to not skip.
               exit
            endif
         enddo
         if (skipKP == 1) then
            cycle
         endif

         ! This k-point survived the skip test above and its work will be
         !   done, so count it (PSEUDOCODE 30).
         valeRhoKPointsProcessed = valeRhoKPointsProcessed + 1

         ! Determine if we are doing the valeCharge in a post-SCF calculation
         !   or within an SCF calculation.
         ! This read is timed on its own (PSEUDOCODE 30 region R) rather
         !   than being folded into the accumulation that follows it.
         !   Only the complex build reads wave functions back: it holds
         !   one k-point at a time.  Folding a file read into the
         !   arithmetic total would report reading as computing on
         !   precisely the multi-k-point decks where the two must be
         !   told apart.
         call beginTimedRegion (clockAtRegionStart)
         do h = 1, spin
            if (inSCF == 1) then
               call readDataSCF(h,i,numStates,0) ! Read wave functions only.
            else
               call readDataPSCF(h,i,numStates,0) ! Read wave functions only.
            endif
         enddo
         call endTimedRegion (clockAtRegionStart,valeRhoReadVectorSeconds)
!      endif

      ! Open the accumulation region (PSEUDOCODE 30 region A, first
      !   span).  It opens before the zeroing rather than after it
      !   because zeroing writes the whole valence-by-valence matrix,
      !   which is the same memory traffic the accumulation itself is
      !   made of and belongs with the half it serves.
      call beginTimedRegion (clockAtRegionStart)

      ! Initialize matrix to receive the valeVale density matrix (square of the
      !   wave function).
      valeValeRho(:,:,:) = cmplx(0.0_double,0.0_double,double)
#else
      ! The gamma build keeps its single k-point resident, so there is no
      !   read to time here and region R stays exactly zero for the whole
      !   run.  That zero is worth seeing in the log rather than being
      !   inferred from the build name.
      valeRhoKPointsProcessed = valeRhoKPointsProcessed + 1

      ! Open the accumulation region (PSEUDOCODE 30 region A, first
      !   span).  See the note in the complex arm above.
      call beginTimedRegion (clockAtRegionStart)

      ! All of the information we need is already available in system memory.
      !   The only thing we need to do is initialize this matrix to zero.
      valeValeRhoGamma(:,:,:) = 0.0_double
#endif


      ! Accumulate the valeValeRho matrix upper triangle and electron energy.
#ifndef GAMMA
      do j = 1, numStates
         currentPopulation(:) = structuredElectronPopulation(j,i,:)

         if (sum(abs(currentPopulation(:))) < smallThresh) cycle

         electronEnergy(:) = electronEnergy(:) + currentPopulation(:) * &
               & energyEigenValues(j,i,:)

         do k = 1, spin
            ! Count the update before issuing it (PSEUDOCODE 30).  Each
            !   one streams the whole upper triangle of the density
            !   matrix through memory, so this count times that
            !   triangle's size is the memory traffic of the
            !   accumulation -- the quantity DESIGN 9.6 estimates and
            !   this measurement exists to check.
            valeRhoRankUpdateCount = valeRhoRankUpdateCount + 1
            call zher('U',valeDim,currentPopulation(k),valeVale(:,j,k),1,&
                  & valeValeRho(:,:,k),valeDim)
         enddo
      enddo

      ! In the event that the calculation includes plusUJ terms, then we need
      !   to compute the plusUJ potential from each atom with such a
      !   contribution. Of course, we also want to do it as efficiently as
      !   possible. The computation follows equation #9 from Anisimov VI,
      !   Zaanen J, Andersen OK. Band theory and Mott insulators: Hubbard U
      !   instead of Stoner I. Physical Review B, 1991;44(3):943. Available
      !   from: http://dx.doi.org/10.1103/PhysRevB.44.943.
      ! We will perform the update in two phases so that we are only doing work
      !   when the appropriate data structures are most available (to avoid
      !   having to do any extra data-reading or data-transfer just for the
      !   purpose of updating the plusUJ terms). We will do phase 1 here
      !   because at this point in the program we have access to the charge
      !   density matrix with all kpoint and electron energy level population
      !   effects accounted for already. The charge density matrix is not
      !   currently packed so it is fairly easy to reference the matrix
      !   elements. Also, at this point the charge denstiy matrix still has a
      !   up and down spin representation instead of a total and up minus down
      !   representation.

      ! Start the computation of the plusUJ terms if there are any.
      if (numPlusUJAtoms > 0) then
         call update1UJ(i,valeValeRho)
      endif

      ! Pack the matrix for easy comparison with the hamiltonian terms
      !   to be read in next.  Store the result in the appropriate packed
      !   valeVale spin array.
      do j = 1, spin
         call packMatrix(valeValeRho(:,:,j),packedValeValeRho(:,:,j),&
               & valeDim)
      enddo
#else
      do j = 1, numStates
         currentPopulation(:) = structuredElectronPopulation(j,i,:)
         if (sum(abs(currentPopulation(:))) < smallThresh) cycle

         electronEnergy(:) = electronEnergy(:) + currentPopulation(:) * &
               & energyEigenValues(j,i,:)

         do k = 1, spin
            ! Count the update before issuing it (PSEUDOCODE 30).  See
            !   the note in the complex arm above for what the count is
            !   used to compute.
            valeRhoRankUpdateCount = valeRhoRankUpdateCount + 1
            call dsyr('U',valeDim,currentPopulation(k),&
                  & valeValeGamma(:,j,k),1,valeValeRhoGamma(:,:,k),valeDim)
         enddo
      enddo

      ! Note: the documentation written above for the plusUJ applies here too.
      if (numPlusUJAtoms > 0) then
         call update1UJ(i,valeValeRhoGamma)
      endif

      ! Pack the matrix for easy comparison with the hamiltonian terms
      !   to be read in next.  Store the result in the appropriate packed
      !   valeVale spin array.
      do j = 1, spin
         call packMatrixGamma(valeValeRhoGamma(:,:,j),&
               & packedValeValeRho(:,:,j),valeDim)
      enddo
#endif

      ! Close the accumulation region (PSEUDOCODE 30 region A, first
      !   span).  Placed after the build-variant block so that one
      !   statement closes the span for both the complex and the real
      !   arm.
      call endTimedRegion (clockAtRegionStart,valeRhoAccumulateSeconds)


      ! Allocate space to hold the overlap matrix. Later, this will also be
      !   used to read in the Hamiltonian matrix terms (KE, nuclear, electronic
      !   potential, and (if needed) the mass velocity).
      ! The allocation belongs to neither region and is deliberately left
      !   outside both, so the four regions sum to slightly less than the
      !   stage as a whole rather than appearing to account for all of it.
      allocate (packedValeVale(dim1,valeDim*(valeDim+1)/2))

      ! Read the overlap matrix into the packedValeVale representation.
      ! Reading a stored matrix and contracting against it alternate from
      !   here to the end of the k-point body, and they are timed apart
      !   (PSEUDOCODE 30 regions I and M).  The ratio between them is
      !   what says whether holding these matrices in memory would be
      !   worth more than dividing the work of reading them, which a
      !   single combined total could not distinguish.
      call beginTimedRegion (clockAtRegionStart)
      if (inSCF == 1) then
         call readPackedMatrix (atomOverlap_did(i),packedValeVale,&
               & packedVVDims,dim1,valeDim)
      else
         call readPackedMatrix (atomOverlapPSCF_did(i),packedValeVale,&
               & packedVVDimsPSCF,dim1,valeDim)
      endif
      call endTimedRegion (clockAtRegionStart,valeRhoReadIntegralSeconds)

      ! In the case that the calculation is spin polarized (spin=2) then we
      !   need to convert the values in the packedValeValeRho density matrix
      !   from being spin up and spin down to being spin up + spin down and
      !   spin up - spin down.  This is accomplished with a temporary variable.
      ! This rewrite of the density matrix is part of building it, not of
      !   contracting it, so it is timed into the accumulation total as
      !   its second span (PSEUDOCODE 30 region A).  It sits here, rather
      !   than beside the packing above, only because the overlap read is
      !   what the recombined density is first contracted against.
      !   Leaving it untimed would cost nothing on every deck in the
      !   benchmark set, all of which are unpolarized and skip it
      !   entirely, and would then quietly misreport the first
      !   spin-polarized deck anyone measured.
      if (spin == 2) then
         call beginTimedRegion (clockAtRegionStart)
         do j = 1, valeDim*(valeDim+1)/2
            tempDensity(:) = packedValeValeRho(:,j,1)
            packedValeValeRho(:,j,1)=tempDensity(:)+packedValeValeRho(:,j,2)
            packedValeValeRho(:,j,2)=tempDensity(:)-packedValeValeRho(:,j,2)
         enddo
         call endTimedRegion (clockAtRegionStart,valeRhoAccumulateSeconds)
      endif

      ! In the next section, we will read in the hamiltonian terms (and
      !   overlap) and perform an element by element multiplication with the
      !   density matrix for each term.

      ! Compute the integration of the charge density to show that <Psi|Psi> is
      !   actually equal to the expected number of eletrons. Note that in the
      !   spin polarized case the first index will refer to the total number of
      !   electrons and the second index will refer to the spin difference.

      if (inSCF == 1) then

         call beginTimedRegion (clockAtRegionStart)
         do j = 1, spin ! j=1 -> Total; j=2 -> Difference
#ifndef GAMMA
            call matrixElementMult (chargeDensityTrace(j),packedValeVale,&
                  & packedValeValeRho(:,:,j),dim1,valeDim)
#else
            call matrixElementMultGamma (chargeDensityTrace(j),&
                  & packedValeVale,packedValeValeRho(:,:,j),dim1,valeDim)
#endif
         enddo
         call endTimedRegion (clockAtRegionStart,valeRhoContractSeconds)

         ! Compute the nuclear contribution to the fitted potential first.
         call beginTimedRegion (clockAtRegionStart)
         call readPackedMatrix (atomNPOverlap_did(i),packedValeVale,&
               & packedVVDims,dim1,valeDim)
         call endTimedRegion (clockAtRegionStart,valeRhoReadIntegralSeconds)
         call beginTimedRegion (clockAtRegionStart)
         do j = 1, spin ! j=1 -> Total; j=2 -> Difference
#ifndef GAMMA
            call matrixElementMult (nucPotTrace(j),packedValeVale,&
                  & packedValeValeRho(:,:,j),dim1,valeDim)
#else
            call matrixElementMultGamma (nucPotTrace(j),packedValeVale,&
                  & packedValeValeRho(:,:,j),dim1,valeDim)
#endif
         enddo
         call endTimedRegion (clockAtRegionStart,valeRhoContractSeconds)

         ! Now compute the kinetic energy.
         call beginTimedRegion (clockAtRegionStart)
         call readPackedMatrix (atomKEOverlap_did(i),packedValeVale,&
               & packedVVDims,dim1,valeDim)
         call endTimedRegion (clockAtRegionStart,valeRhoReadIntegralSeconds)
         call beginTimedRegion (clockAtRegionStart)
         do j = 1, spin ! j=1 -> Total; j=2 -> Difference
#ifndef GAMMA
            call matrixElementMult (kineticEnergyTrace(j),packedValeVale,&
                  & packedValeValeRho(:,:,j),dim1,valeDim)
#else
            call matrixElementMultGamma (kineticEnergyTrace(j),packedValeVale,&
                  & packedValeValeRho(:,:,j),dim1,valeDim)
#endif
         enddo
         call endTimedRegion (clockAtRegionStart,valeRhoContractSeconds)

         ! If needed, compute the mass velocity.
         if (rel == 1) then
            call beginTimedRegion (clockAtRegionStart)
            call readPackedMatrix (atomMVOverlap_did(i),packedValeVale,&
                  & packedVVDims,dim1,valeDim)
            call endTimedRegion (clockAtRegionStart,&
                  & valeRhoReadIntegralSeconds)
            call beginTimedRegion (clockAtRegionStart)
            do j = 1, spin ! j=1 -> Total; j=2 -> Difference
#ifndef GAMMA
               call matrixElementMult (massVelocityTrace(j),packedValeVale,&
                     & packedValeValeRho(:,:,j),dim1,valeDim)
#else
               call matrixElementMultGamma (massVelocityTrace(j),packedValeVale,&
                     & packedValeValeRho(:,:,j),dim1,valeDim)
#endif
            enddo
            call endTimedRegion (clockAtRegionStart,valeRhoContractSeconds)
         endif

         ! Loop over atomic potential terms next.
         ! These are the bulk of both regions: one stored matrix per
         !   potential term, read and then contracted, so the two timed
         !   regions below are opened and closed potDim times each.
         do j = 1, potDim
            call beginTimedRegion (clockAtRegionStart)
            call readPackedMatrix (atomPotOverlap_did(i,j),packedValeVale,&
                  & packedVVDims,dim1,valeDim)
            call endTimedRegion (clockAtRegionStart,&
                  & valeRhoReadIntegralSeconds)
            call beginTimedRegion (clockAtRegionStart)
            do k = 1, spin ! j=1 -> Total; j=2 -> Difference
#ifndef GAMMA
               call matrixElementMult (potRho(j,k),packedValeVale,&
                     & packedValeValeRho(:,:,k),dim1,valeDim)
#else
               call matrixElementMultGamma (potRho(j,k),packedValeVale,&
                     & packedValeValeRho(:,:,k),dim1,valeDim)
#endif
            enddo
            call endTimedRegion (clockAtRegionStart,valeRhoContractSeconds)
         enddo
      endif ! inSCF == 1

      ! If needed, compute the forces
      if (((doForce_SCF == 1) .and. (converged == 1)) .or. &
            & ((doForce_PSCF == 1) .and. (inSCF == 0))) then
         do j = 1, 3 ! xyz directions
            do k = 1, spin
#ifndef GAMMA
!               call packMatrix(valeValeF(:,:,i,k,j),packedValeVale(:,:),&
!                     & valeDim)
!               valeValeRho(:,:,k) = valeValeRho(:,:,k) + transpose(valeValeRho(:,:,k))
!               do l = 1, valeDim
!                  valeValeRho(l,l,k) = valeValeRho(l,l,k) / 2.0_double
!               enddo
!               valeValeF(:,:,i,k,j) = valeValeF(:,:,i,k,j) * valeValeRho(:,:,k)
!   packedValeVale(1,:) = packedValeVale(1,:)*packedValeValeRho(1,:,spin) + &
!         & packedValeVale(2,:)*packedValeValeRho(2,:,spin)
               call computeForce(valeValeRho,i,k,j)
#else
!               valeValeRhoGamma(:,:,k) = valeValeRhoGamma(:,:,k) + &
!                     & transpose(valeValeRhoGamma(:,:,k))
!               do l = 1, valeDim
!                  valeValeRhoGamma(l,l,k) = valeValeRhoGamma(l,l,k) / 2.0_double
!               enddo
!               valeValeFGamma(:,:,i,k,j) = valeValeFGamma(:,:,i,k,j) * valeValeRhoGamma(:,:,k)
!               call packMatrixGamma(valeValeFGamma(:,:,k,j),&
!                     & packedValeVale(:,:),valeDim)
               call computeForceGamma(valeValeRhoGamma,k,j)
#endif
            enddo
         enddo
      endif

      ! Deallocate the space no longer needed.
      deallocate (packedValeVale)

   enddo ! i   numKPoints

   ! Now that all k-points have been accumulated, record the total and spin
   !   difference number of electrons in the system.
   if (inSCF == 1) then
      write (20,*) "Total number of electrons from <psi|psi> = ",&
            & chargeDensityTrace(1)
      if (spin == 2) then
         write (20,*) "Electron spin difference from <psi|psi> = ",&
               & chargeDensityTrace(2)
      endif

      if (spin == 1) then

         ! Obtain a summation of the electrons.
         if (rel == 0) then
            sumElecEnergy = sum(potCoeffs(:,1) * potRho(:,1)) + &
                  & kineticEnergyTrace(1) + nucPotTrace(1)
         else
            sumElecEnergy = sum(potCoeffs(:,1) * potRho(:,1)) + &
                  & kineticEnergyTrace(1) + massVelocityTrace(1) + &
                  & nucPotTrace(1)
         endif


         ! Write the two electron numbers (both are energy values, but they
         !   are computed through different means.)
         write (20,fmt='(a17,f18.8)') 'Electron Energy = ',electronEnergy(1)
         write (20,fmt='(a17,f18.8)') 'Electron Sum    = ',sumElecEnergy
      else

         ! Obtain a summation of the electrons for each spin.

         ! Do total plus the spin difference to obtain the spin up.
         if (rel == 0) then
            sumElecEnergy = 0.5_double * (sum(potCoeffs(:,1) * &
                  & (potRho(:,1) + potRho(:,2))) + kineticEnergyTrace(1) + &
                  & kineticEnergyTrace(2) + nucPotTrace(1) + nucPotTrace(2))
         else
            sumElecEnergy = 0.5_double * (sum(potCoeffs(:,1) * &
                  & (potRho(:,1) + potRho(:,2))) + kineticEnergyTrace(1) + &
                  & kineticEnergyTrace(2) + massVelocityTrace(1) + &
                  & massVelocityTrace(2) + nucPotTrace(1) + nucPotTrace(2))
         endif
         ! Write the two electron numbers. (Both are energy values, but they
         !   are computed thorugh different means.
         write (20,fmt='(a24,f18.8)') '(UP)   Electron Energy = ', &
               & electronEnergy(1)
         write (20,fmt='(a24,f18.8)') '(UP)   Electron Sum    = ', &
               & sumElecEnergy


         ! Do total minus the spin difference to obtain the spin down.
         if (rel == 0) then
            sumElecEnergy = 0.5_double * (sum(potCoeffs(:,2) * &
                  & (potRho(:,1) - potRho(:,2))) + kineticEnergyTrace(1) - &
                  & kineticEnergyTrace(2) + nucPotTrace(1) - nucPotTrace(2))
         else
            sumElecEnergy = 0.5_double * (sum(potCoeffs(:,2) * &
                  & (potRho(:,1) - potRho(:,2))) + kineticEnergyTrace(1) - &
                  & kineticEnergyTrace(2) + massVelocityTrace(1) - &
                  & massVelocityTrace(2) + nucPotTrace(1) - nucPotTrace(2))
         endif

         ! Write the spin down electronEnergy???
         write (20,fmt='(a24,f18.8)') '(DOWN) Electron Energy = ', &
               & electronEnergy(2)
         write (20,fmt='(a24,f18.8)') '(DOWN) Electron Sum    = ', &
               & sumElecEnergy
      endif
   endif

   ! Report how this call divided between building the density matrix
   !   and contracting it (PSEUDOCODE 30).  One line per call, written
   !   inside the stage's own start/end stamps so a reader finds it where
   !   the stage is, and carrying a fixed leading label so that a whole
   !   run's worth can be pulled out of the log with a single search.
   !
   ! Raw values only.  Every ratio, rate and projection a reader wants
   !   from these -- the share spent accumulating, the memory traffic
   !   implied by the update count, the read rate per stored matrix -- is
   !   computed outside the program, because a derived number written
   !   here is one that nobody re-checks when the derivation behind it
   !   turns out to be wrong.
   !
   ! inSCF is reported because the post-SCF path calls this same routine
   !   with a different set of stored matrices, and a reading that
   !   silently mixed the two would be worse than no reading.
   if (valeRhoClockRate > 0) then
      write (20,fmt='(a,6(1x,i0),4(1x,es14.6),2(1x,i0))') &
            & 'VALEDENSITY SPLIT', inSCF, numKPoints, valeDim, &
            & numStates, potDim, spin, &
            & valeRhoReadVectorSeconds, valeRhoAccumulateSeconds, &
            & valeRhoReadIntegralSeconds, valeRhoContractSeconds, &
            & valeRhoRankUpdateCount, valeRhoKPointsProcessed
   else
      write (20,fmt='(a,6(1x,i0),a,2(1x,i0))') &
            & 'VALEDENSITY SPLIT', inSCF, numKPoints, valeDim, &
            & numStates, potDim, spin, &
            & '  seconds unavailable (no system clock)', &
            & valeRhoRankUpdateCount, valeRhoKPointsProcessed
   endif
   call flush (20)

   ! Deallocate arrays associated with the electron population.
   call cleanUpPopulation

   ! Deallocate arrays associated with the secular equation and its solution.
   call cleanUpSecularEqn

   ! Deallocate other arrays and matrices defined in this subroutine.
   deallocate (structuredElectronPopulation)
   deallocate (packedValeValeRho)
   deallocate (currentPopulation)
   deallocate (electronEnergy)
   deallocate (tempDensity)
#ifndef GAMMA
   deallocate (valeValeRho)
#else
   deallocate (valeValeRhoGamma)
#endif

   ! Log the date and time we end.
   call timeStampEnd (17)

end subroutine makeValenceRho


end module O_ValeCharge

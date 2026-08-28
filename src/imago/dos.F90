!! SPDX-License-Identifier: ECL-2.0
!! Copyright (c) 2026 Paul Rulis

module O_DOS

   ! Import the necessary modules.
   use O_Kinds

   ! Make sure that no variables are declared accidentally.
   implicit none

   ! Define module data specifically for TDOS calculation for each SCF cycle.
   real (kind=double), allocatable, dimension (:,:,:) :: tdos
   real (kind=double), allocatable, dimension (:)     :: currentEnergyValues

   ! Define module data for the general PDOS/TDOS calculation.
   integer :: numEnergyPoints
   real (kind=double), allocatable, dimension (:)     :: energyScale

   contains

subroutine computeIterationTDOS

   ! Import the necessary modules.
   use O_Kinds
   use O_Constants,       only: pi, hartree
   use O_Populate,        only: occupiedEnergy
   use O_SecularEquation, only: energyEigenValues
   use O_KPoints,         only: numKPoints, kPointWeight
   use O_Potential,       only: spin, lastIteration, currIteration
   use O_Input, only: sigmaDOS, eminDOS, emaxDOS, deltaDOS, numStates

   ! Make sure that no variables are declared accidentally.
   implicit none

   ! Define local variables
   integer :: i,j,k,l
   real (kind=double) :: sigmaSqrtPi
   real (kind=double) :: expTerm
   real (kind=double) :: expFactor

   ! Normalization for Gaussian broadening.
   sigmaSqrtPi = sqrt(pi) * sigmaDOS

   if (.not. allocated(tdos)) then

      ! Determine the number of energy buckets to be computed for.
      numEnergyPoints = int((emaxDOS - eminDOS ) / deltaDOS)

      ! Allocate memory to hold the data for those points. The funny-stuff
      !   with last iteration: If lastIteration is zero, that is a special
      !   case. We still just compute one iteration, but it is also a signal
      !   to not re-compute the wave function later.
      if (lastIteration > 0) then
         allocate (tdos(numEnergyPoints,spin,lastIteration))
      else
         allocate (tdos(numEnergyPoints,spin,1))
      endif
      allocate (energyScale(numEnergyPoints))
      allocate (currentEnergyValues(numStates))

      ! Initialize the data for accumulation.
      tdos(:,:,:) = 0.0_double

      ! Compute the values for the energy scale.
      do i = 1, numEnergyPoints
         energyScale(i) = eminDOS + (i-1) * deltaDOS
      enddo
   endif


   do i = 1, spin
      do j = 1, numKPoints

         ! Obtain a copy of the current energy values for this spin and kpoint,
         !   and shift them in accordance  with the highest occupied energy
         !   level.
         currentEnergyValues(:) = (energyEigenValues(:,j,i) - occupiedEnergy)

         do k = 1, numStates

            ! Apply broadening to each state.
            do l = 1, numEnergyPoints

               ! Compute the exponential term for the broadening of this point.
               !   Note that from the usual Gaussian term of exp(-alpha*x^2)
               !   we are have (eV-eS)^2 = x^2 and (1/sigma)^2 = alpha.
               expTerm = ((currentEnergyValues(k)-energyScale(l))/sigmaDOS)**2

               ! If the exponential term is less than 50 we apply the
               !   broadening.
               if (expTerm < 50.0_double) then

                  ! Compute the exponential factor.  It is at this point that
                  !   the kpoint weighting factor is applied to make the energy
                  !   bucket we are considering be filled with the right number
                  !   of electrons.  Note that we must account for the value
                  !   of "spin".  This will put either 1 or 2 electrons in each
                  !   state.
                  expFactor = exp(-expTerm) / sigmaSqrtPi / hartree * &
                        & kPointWeight(j) / real(spin,double)

                  ! Store the broadened TDOS.
                  tdos(l,i,currIteration) = tdos(l,i,currIteration) + expFactor
               endif
            enddo
         enddo
      enddo
   enddo

end subroutine computeIterationTDOS

subroutine printIterationTDOS

   ! Import the necessary modules.
   use O_Kinds
   use O_Constants, only: hartree
   use O_Potential, only: spin, currIteration

   ! Make sure that no variables are declared accidentally.
   implicit none

   ! Define local variables
   integer :: i,j
!   integer :: pointCount
   integer :: fileID
   character*9 :: fileName

   ! Note that the current iteration was already incremented by one when this
   !   subroutine is called so we must refer back to the previous iteration
   !   when accessing data.

   do i = 1, spin

      ! Define the file ID to open.
      fileID = currIteration+i*1000

      ! Define the file name.
      write (fileName,fmt="(a5,i4)") "fort.",fileID

      ! Open the file for creating the opendx output.
      open (unit=fileID,file=fileName,status='unknown',form='formatted')


do j = 1, numEnergyPoints
! Make sure to convert to eV upon output.
write (fileID,fmt="(1x,3f12.8)") energyScale(j)*hartree,tdos(j,1,1),&
      & tdos(j,1,currIteration-1)

enddo

   enddo


   ! Deallocate unused memory.
   deallocate (tdos)
   deallocate (energyScale)
   deallocate (currentEnergyValues)

end subroutine printIterationTDOS


subroutine computeDOS(inSCF)

   ! Import the necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Potential,   only: spin
   use O_Populate,    only: electronPopulation
   use O_KPoints, only: numKPoints, kPointWeight, &
         & kPointIntgCode, numPointOps, &
         & symmetrizeLATPartials
   use O_Constants, only: pi, hartree, lAngMomCount
   use O_AtomicSites, only: valeDim, &
         & numAtomSites, atomSites, invAtomPerm
   use O_AtomicTypes, only: numAtomTypes, &
         & atomTypes, maxNumValeStates
   use O_Input, only: numStates, sigmaDOS, &
         & eminDOS, emaxDOS, deltaDOS, &
         & detailCodePDOS
#ifndef GAMMA
   use O_SecularEquation, only: valeValeOL, &
         & valeVale, energyEigenValues, &
         & readDataSCF, readDataPSCF
#else
   use O_SecularEquation, only: &
         & valeValeOLGamma, valeValeGamma, &
         & energyEigenValues, &
         & readDataSCF, readDataPSCF
#endif
   ! The shared state-projection routine forms T = S C, every state
   !   carried through the overlap, as one matrix product per k-point
   !   (PSEUDOCODE 34).
   use O_StateProjection, only: projectStatesOntoBasis


   ! Make sure that no funny variables are used.
   implicit none

   ! Define passed parameters.
   integer, intent(in) :: inSCF

   ! Define local variables.
   integer :: h, i, j, k, l  ! Loop index variables.
   character*17 :: formatString
   character*1, dimension (lAngMomCount) :: &
         & QN_lLetter
   character*14, dimension (4,7) :: QN_mLetter
   integer, allocatable, dimension (:) :: cumulNumDOS
   integer, allocatable, dimension (:) :: pdosIndex
         !   Given a valence state, this stores the place within pdosAccum that
         !   this state should be stored.
   integer, allocatable, dimension (:) :: &
         & numAtomStates
   integer :: cumulDOSTotal
   integer :: numSQN_l
   integer :: numPQN_l
   integer :: numDQN_l
   integer :: numFQN_l
   integer :: initSIndex
   integer :: initPIndex
   integer :: initDIndex
   integer :: initFIndex
   integer :: currentType
   integer :: valeDimIndex
   integer :: initIndex
   integer :: finIndex
   integer :: energyLevelCounter
   integer :: numCols
   integer :: stateSpinKPointIndex
   real (kind=double) :: occupancyNumber
   real (kind=double) :: oneValeRealAccum
   real (kind=double) :: expTerm
   real (kind=double) :: expFactor
   real (kind=double) :: sigmaSqrtPi
   real (kind=double) :: integratedArea
   real (kind=double) :: numStatesInRange
   real (kind=double) :: totalElectronsComputed
   real (kind=double) :: electronFactor
   real (kind=double) :: currentPopulation
   real (kind=double), allocatable, &
         & dimension (:)   :: pdosAccum
   real (kind=double), allocatable, &
         & dimension (:)   :: localizationIndex
   real (kind=double), allocatable, &
         & dimension (:)   :: totalSystemDos
   real (kind=double), allocatable, &
         & dimension (:)   :: energyValuesAvg
   real (kind=double), allocatable, &
         & dimension (:,:) :: pdosComplete
   real (kind=double), allocatable, &
         & dimension (:,:) :: electronNumber
   ! The state projection T = S C (PSEUDOCODE 34): every state carried
   !   through the overlap, formed once per k-point as one matrix
   !   product.  It is complex in the multi-k build and real in the
   !   gamma build.
#ifndef GAMMA
   complex (kind=double), allocatable, &
         & dimension (:,:) :: statesProjected
#else
   real (kind=double), allocatable, &
         & dimension (:,:) :: statesProjected
#endif
   ! The Mulliken projection P(mu,j) = Re(conjg(C(mu,j)) T(mu,j)),
   !   real in both builds, formed once per k-point from the
   !   eigenvectors and T.  The atom-by-orbital binning loop reads one
   !   entry of it per basis function and state (PSEUDOCODE 34.5).
   real (kind=double), allocatable, &
         & dimension (:,:) :: mullikenProj

   ! LAT PDOS variables (used only when kPointIntgCode == 1).
   real (kind=double), allocatable, &
         & dimension(:,:,:) :: projArray

   ! The channel permutation table is needed by BOTH pathways, for
   !   different reasons (PSEUDOCODE 23.1), so it is not a LAT variable.
   integer, allocatable, &
         & dimension(:,:) :: channelPermTable
   logical :: buildPermTable ! Whether this run needs that table.

   ! Log the date and time we start.
   call timeStampStart (19)

   ! Mode 3 restriction: per-atom per-lm PDOS with LAT requires D^l rotation
   !   matrices for the channel permutation, which are not available. Stop with
   !   a clear error (DESIGN 1.4).
   if (kPointIntgCode == 1 .and. &
         & detailCodePDOS == 3) then
      write (20, *) 'ERROR: LAT PDOS (kPointIntg Code=1) with per-atom'
      write (20, *) 'per-lm detail (detailCodePDOS= 3) requires D^l rotation'
      write (20, *) 'matrices that are not available. Use mode 0, 1, or 2'
      write (20, *) 'with LAT integration.'
      stop 'computeDOS: mode 3 + LAT unsupported'
   endif

   ! Normalization for Gaussian broadening.
   sigmaSqrtPi = sqrt(pi) * sigmaDOS

   ! Allocate arrays and matrices for this computation.
   allocate (pdosIndex     (valeDim))
   allocate (numAtomStates (numAtomSites))
   ! Store DOS for each type's orbital sum.
   if     (detailCodePDOS == 0) then
      allocate (cumulNumDOS (numAtomTypes + 1))
   ! Store DOS for each atom's TDOS.
   elseif (detailCodePDOS == 1) then
      allocate (cumulNumDOS (1))
   ! Store DOS for each QN_nl resolved atom.
   elseif (detailCodePDOS == 2) then
      allocate (cumulNumDOS (numAtomSites + 1))
   ! Store DOS for each QN_nlm resolved atom.
   elseif (detailCodePDOS == 3) then
      allocate (cumulNumDOS (numAtomSites + 1))
   endif

   ! The pdosAccum array, the projection matrices, and the overlap are
   !   only needed for the Gaussian path. The LAT path handles its own
   !   allocations inside computeProjections_LAT.
   ! The eigenvector array is ONE slab on the post-SCF path: the spins
   !   are processed one at a time and the reader delivers each spin's
   !   vectors to slab 1 on request (slab = 1 at the read).  On the SCF
   !   path the solver's own spin-slab array is still allocated and the
   !   reader delivers into its slab 1 the same way (DESIGN 2.8).
   if (kPointIntgCode /= 1) then
      if     (detailCodePDOS == 0) then
         allocate (pdosAccum (valeDim))
      elseif (detailCodePDOS == 1) then
         allocate (pdosAccum (numAtomSites))
      elseif (detailCodePDOS == 2) then
         allocate (pdosAccum (valeDim))
      elseif (detailCodePDOS == 3) then
         allocate (pdosAccum (valeDim))
      endif
      allocate (statesProjected(valeDim, numStates))
      allocate (mullikenProj   (valeDim, numStates))
#ifndef GAMMA
      allocate (valeValeOL(valeDim, valeDim))
      if (inSCF == 0) then
         allocate (valeVale(valeDim, numStates, 1))
      endif
#else
      allocate (valeValeOLGamma(valeDim, valeDim))
      if (inSCF == 0) then
         allocate (valeValeGamma( &
               & valeDim, numStates, 1))
      endif
#endif
   endif

   ! Define the QN_l letters.
   QN_lLetter(1) = 's'
   QN_lLetter(2) = 'p'
   QN_lLetter(3) = 'd'
   QN_lLetter(4) = 'f'

   ! Define the QN_m resolved letters.
   QN_mLetter(1,1) = 'r'
   QN_mLetter(2,1) = 'x'
   QN_mLetter(2,2) = 'y'
   QN_mLetter(2,3) = 'z'
   QN_mLetter(3,1) = 'xy'
   QN_mLetter(3,2) = 'xz'
   QN_mLetter(3,3) = 'yz'
   QN_mLetter(3,4) = 'xx~yy'
   QN_mLetter(3,5) = '2zz~xx~yy'
   QN_mLetter(4,1) = 'xyz'
   QN_mLetter(4,2) = 'xxz~yyz'
   QN_mLetter(4,3) = 'xxx~3yyx'
   QN_mLetter(4,4) = '3xxy~yyy'
   QN_mLetter(4,5) = '2zzz~3xxz~3yyz'
   QN_mLetter(4,6) = '4zzx~xxx~yyx'
   QN_mLetter(4,7) = '4zzy~xxy~yyy'

   ! Initialize other variables.
   cumulNumDOS(:) = 0
   cumulDOSTotal  = 0

   if (detailCodePDOS == 0) then

      ! Initialize counter to index the cumulative sum of QN_l orbitals for all
      !   atomic types.  (This PDOS will give a DOS for each QN_nl pair of each
      !   atomic type.)
      cumulNumDOS(1) = 0

      ! Loop to record the number of orbitals that each type contributes.  (An
      !   orbital is just a QN_nl pair.)
      do i = 1, numAtomTypes
         cumulNumDOS(i+1) = cumulNumDOS(i) + &
               & sum(atomTypes(i)%numQN_lValeRadialFns(:))
      enddo

      ! Record the total number of orbitals summed over all types.
      cumulDOSTotal = cumulNumDOS(numAtomTypes+1)

   elseif (detailCodePDOS == 1) then

      ! There will be one DOS curve for each atom and that is it.
      cumulDOSTotal = numAtomSites

   elseif (detailCodePDOS == 2) then

      ! Initialize counter to index the cumulative sum of QN_l orbitals for all
      !   atoms.  (This PDOS will give a DOS for each QN_nl pair of each atom.)
      cumulNumDOS(1) = 0

      ! Loop to record the number of orbitals that for each atom contributes.
      !   (An orbital is just a QN_nl pair.)
      do i = 1, numAtomSites
         cumulNumDOS(i+1) = cumulNumDOS(i) + sum( &
               & atomTypes(atomSites(i)%atomTypeAssn)%numQN_lValeRadialFns(:))
      enddo

      ! Record the total number of orbitals summed over all atoms.
      cumulDOSTotal = cumulNumDOS(numAtomSites+1)

   elseif (detailCodePDOS == 3) then

      ! Initialize counter to index the cumulative sum of QN_l orbitals for
      !   all atoms.  (This PDOS will give a DOS for each QN_nlm set for each
      !   atom.)
      cumulNumDOS(1) = 0

      ! Loop to record the index number for each atom's orbitals.
      do i = 1, numAtomSites
         cumulNumDOS(i+1) = cumulNumDOS(i) + &
            & atomTypes(atomSites(i)%atomTypeAssn)%numQN_lValeRadialFns(1)*1 + &
            & atomTypes(atomSites(i)%atomTypeAssn)%numQN_lValeRadialFns(2)*3 + &
            & atomTypes(atomSites(i)%atomTypeAssn)%numQN_lValeRadialFns(3)*5 + &
            & atomTypes(atomSites(i)%atomTypeAssn)%numQN_lValeRadialFns(4)*7 
      enddo

      ! Record the total number of QN_m resolved orbitals summed over all atoms.
      cumulDOSTotal = cumulNumDOS(numAtomSites+1)

   endif

   ! Initialize valeDimIndex to record the index number in valeDim that each
   !   pdos state is at.
   valeDimIndex = 0

   ! Loop over every atom in the system to index where the pdos values for
   !   each atom should be stored.
   do i = 1, numAtomSites

      ! Obtain the type of the current atom.
      currentType = atomSites(i)%atomTypeAssn

      ! Identify and store the number of valence states for this atom.
      numAtomStates(i) = atomTypes(currentType)%numValeStates

      ! In the case where the PDOS should be collected by types we loop
      !   through each QN_nl pair for this atom and record the index where its
      !   PDOS should be recorded according to its type.  In the case where
      !   the PDOS is collected by atoms we link each QN_nl pair of this atom
      !   with the index of this atom.  In the case where the PDOS is collected
      !   for each QN_nl pair of each atom the index is as in the first case
      !   that now we index for each atom as opposed to each type.
      if (detailCodePDOS == 0) then

         numSQN_l = atomTypes(currentType)%numQN_lValeRadialFns(1)
         numPQN_l = atomTypes(currentType)%numQN_lValeRadialFns(2)
         numDQN_l = atomTypes(currentType)%numQN_lValeRadialFns(3)
         numFQN_l = atomTypes(currentType)%numQN_lValeRadialFns(4)

         initSIndex = cumulNumDOS(currentType)
         initPIndex = cumulNumDOS(currentType) + numSQN_l
         initDIndex = cumulNumDOS(currentType) + numSQN_l + numPQN_l
         initFIndex = cumulNumDOS(currentType) + numSQN_l + numPQN_l + numDQN_l

         do j = 1, numSQN_l
            valeDimIndex = valeDimIndex + 1
            pdosIndex(valeDimIndex) = initSIndex + j
         enddo
         do j = 1, numPQN_l
            do k = 1,3
               valeDimIndex = valeDimIndex + 1
               pdosIndex(valeDimIndex) = initPIndex + j
            enddo
         enddo
         do j = 1, numDQN_l
            do k = 1,5
               valeDimIndex = valeDimIndex + 1
               pdosIndex(valeDimIndex) = initDIndex + j
         enddo
         enddo
         do j = 1, numFQN_l
            do k = 1,7
               valeDimIndex = valeDimIndex + 1
               pdosIndex(valeDimIndex) = initFIndex + j
            enddo
         enddo
      elseif (detailCodePDOS == 1) then
         do j = 1, numAtomStates(i)
            valeDimIndex = valeDimIndex + 1
            pdosIndex(valeDimIndex) = i  ! NOTE THAT THIS IS 'i', NOT 'j'.
         enddo
      elseif (detailCodePDOS == 2) then ! Consider spdf for each atom too.

         numSQN_l = atomTypes(currentType)%numQN_lValeRadialFns(1)
         numPQN_l = atomTypes(currentType)%numQN_lValeRadialFns(2)
         numDQN_l = atomTypes(currentType)%numQN_lValeRadialFns(3)
         numFQN_l = atomTypes(currentType)%numQN_lValeRadialFns(4)

         initSIndex = cumulNumDOS(i)
         initPIndex = cumulNumDOS(i) + numSQN_l
         initDIndex = cumulNumDOS(i) + numSQN_l + numPQN_l
         initFIndex = cumulNumDOS(i) + numSQN_l + numPQN_l + numDQN_l

         do j = 1, numSQN_l
            valeDimIndex = valeDimIndex + 1
            pdosIndex(valeDimIndex) = initSIndex + j
         enddo
         do j = 1, numPQN_l
            do k = 1,3
               valeDimIndex = valeDimIndex + 1
               pdosIndex(valeDimIndex) = initPIndex + j
            enddo
         enddo
         do j = 1, numDQN_l
            do k = 1,5
               valeDimIndex = valeDimIndex + 1
               pdosIndex(valeDimIndex) = initDIndex + j
         enddo
         enddo
         do j = 1, numFQN_l
            do k = 1,7
               valeDimIndex = valeDimIndex + 1
               pdosIndex(valeDimIndex) = initFIndex + j
            enddo
         enddo
      elseif (detailCodePDOS == 3) then
         do j = 1, numAtomStates(i)
            valeDimIndex = valeDimIndex + 1
            pdosIndex(valeDimIndex) = valeDimIndex ! Each QN_nlm is saved.
         enddo
      endif
   enddo

   ! Determine the number of energy buckets to be computed for.
   numEnergyPoints = int((emaxDOS - eminDOS ) / deltaDOS)

   ! Allocate space to hold the pdos and localization index results
   allocate (localizationIndex (numStates))
   allocate (energyScale       (numEnergyPoints))
   allocate (totalSystemDos    (numEnergyPoints))
   allocate (electronNumber    (maxNumValeStates,numAtomSites))
   allocate (energyValuesAvg   (numStates))
   allocate (pdosComplete      (cumulDOSTotal,numEnergyPoints))

   ! Assign values to the energy scale.
   do i = 1, numEnergyPoints
      energyScale(i) = eminDOS + (i-1) * deltaDOS
   enddo

   ! Build the channel permutation table before the spin loop (it is
   !   spin-independent). This maps PDOS channel indices through the inverse
   !   atom permutation.
   !
   ! BOTH pathways need it, for different reasons. The tetrahedron path
   !   permutes once per tetrahedron corner as it assembles them (DESIGN 1.4,
   !   PSEUDOCODE 8.1). The Gaussian path never visits a star member at all,
   !   so it instead group-averages the finished accumulation, which
   !   PSEUDOCODE 23.1 shows is exactly the same thing.
   !
   ! Only the atom-resolved modes need it. Mode 0 is type-resolved and a sum
   !   over the atoms of a type is already invariant, because an operation
   !   carries an atom onto one of the same type (DESIGN 2.3). Mode 3 is not
   !   covered by the table at all -- see the note at the symmetrization call.
   buildPermTable = (detailCodePDOS == 1) .or. (detailCodePDOS == 2)
   if (kPointIntgCode == 1) then
      buildPermTable = .true.
   endif

   ! invAtomPerm is built for every run except SYBD, which needs no unfolding
   !   and never initializes the point-operation machinery (DESIGN 2.6). The
   !   guard is on the array rather than on the run type so that this cannot
   !   disagree with whoever decides that elsewhere.
   if (.not. allocated(invAtomPerm)) then
      buildPermTable = .false.
   endif

   if (buildPermTable) then
      call buildChannelPermTable(detailCodePDOS, &
            & numPointOps, cumulDOSTotal, &
            & cumulNumDOS, numAtomSites, &
            & invAtomPerm, channelPermTable)
   endif

   do h = 1, spin

      ! Track the stateSpinKPoint index number.
      stateSpinKPointIndex = (h-1) * numStates

      ! Record which calculation is being done.
      if (spin == 2) then
         if (h == 1) then
            write (20,*) "Computing spin up DOS."
         else
            write (20,*) "Computing spin down DOS."
         endif
      endif

      ! Initialize various arrays and matrices.
      electronNumber  (:,:) = 0.0_double
      pdosComplete    (:,:) = 0.0_double
      localizationIndex (:) = 0.0_double

      ! Branch on integration method. The LAT path uses a two-pass design
      !   (project then integrate via tetrahedra); the Gaussian path uses a
      !   single-pass design (project and broaden simultaneously). Both fill
      !   pdosComplete, electronNumber, and localizationIndex. The output phase
      !   that follows is shared.
      if (kPointIntgCode == 1) then

         !! - LAT two-pass PDOS (DESIGN 1.4). -

         ! Pass 1: stream IBZ k-points, compute Mulliken projections, store in
         !   projArray. Also accumulates electronNumber and localizationIndex
         !   for the diagnostic.
         write (20, *) "LAT PDOS Pass 1: projections"
         call computeProjections_LAT(inSCF, h, &
               & numKPoints, numStates, &
               & numAtomSites, numAtomStates, &
               & pdosIndex, valeDim, &
               & cumulDOSTotal, spin, projArray, &
               & electronNumber, localizationIndex)

         ! Pass 2: tetrahedron integration with Bloechl corner weights and
         !   channel permutation for IBZ unfolding.
         write (20, *) "LAT PDOS Pass 2: integration"
         call integratePDOS_LAT(projArray, &
               & channelPermTable, pdosComplete, h, &
               & numStates, cumulDOSTotal, &
               & numEnergyPoints, energyScale)

         ! Free the projection array (large memory).
         deallocate (projArray)

         ! Pass 3: average the atom-resolved result over the point group
         !   (DESIGN 1.7). The tetrahedron decomposition is not carried onto
         !   itself by every operation of every point group, so k-points
         !   related by symmetry can receive unrelated weights and atoms that
         !   must be equivalent come out unequal. Averaging over the group
         !   removes precisely that part of the error, and leaves the total
         !   DOS untouched because it is the same operation as averaging each
         !   k-point's weight over its star.
         !
         ! Detail code 0 needs nothing: a type-level sum is already invariant,
         !   since every operation carries an atom onto an atom of the same
         !   type. Detail code 3 never reaches here, being refused on this
         !   pathway earlier.
         if ((symmetrizeLATPartials == 1) .and. &
               & (detailCodePDOS >= 1)) then
            call symmetrizePDOS_LAT(pdosComplete, &
                  & channelPermTable, cumulDOSTotal, &
                  & numEnergyPoints)
         endif

      else

      ! ----------------------------------------- Gaussian single-pass PDOS
      ! (original code). -----------------------------------------

      ! Begin accumulating the DOS values
      do i = 1, numKPoints

         ! Track the stateSpinKPoint index number.
         if ((spin == 2) .and. (i /= 1)) then
            stateSpinKPointIndex = stateSpinKPointIndex + numStates
         endif

         ! Determine if we are doing the DOS in a post-SCF calculation, or
         !   within an SCF calculation.  Either way the eigenvectors of
         !   this spin are delivered to slab 1 (slab = 1): this routine
         !   processes the spins one at a time and holds a single slab,
         !   and it reads slab 1 below (DESIGN 2.8, PSEUDOCODE 33).
         if (inSCF == 1) then
            ! Read necessary data from SCF (setup,main) data structures.
            call readDataSCF(h,i,numStates,1,slab=1) ! 1 = Overlap matrixCode
         else
            ! Read necessary data from post SCF (intg,band) data structures.
            call readDataPSCF(h,i,numStates,1,slab=1) ! 1 = Overlap matrixCode
         endif

         ! Project every state through the overlap in one matrix
         !   product, T = S C, then form the Mulliken projection
         !   P(mu,j) = Re(conjg(C(mu,j)) T(mu,j)) for all states at
         !   once (PSEUDOCODE 34).  This replaces the per-state,
         !   per-basis-function inner loop that walked a whole column
         !   of the overlap for every (state, basis function) pair;
         !   the atom-by-orbital loops below now only read P and place
         !   each value in its channel.  In the gamma build the
         !   eigenvectors are real, so the projection is C(mu,j) *
         !   T(mu,j) with no conjugate.  The eigenvector columns are
         !   sliced explicitly as 1:numStates, because on the SCF path
         !   the solver's slab is (valeDim, valeDim, spin) and only its
         !   first numStates columns are the states of interest; the
         !   post-SCF slab already has exactly numStates columns.
#ifndef GAMMA
         call projectStatesOntoBasis(valeValeOL, &
               & valeVale(:,1:numStates,1), valeDim, numStates, &
               & statesProjected)
         mullikenProj(:,:) = real(conjg(valeVale(:,1:numStates,1)) &
               & * statesProjected(:,:), double)
#else
         call projectStatesOntoBasis(valeValeOLGamma, &
               & valeValeGamma(:,1:numStates,1), valeDim, &
               & numStates, statesProjected)
         mullikenProj(:,:) = valeValeGamma(:,1:numStates,1) * &
               & statesProjected(:,:)
#endif

         do j = 1, numStates

            ! Track the stateSpinKPoint index number.
            stateSpinKPointIndex = stateSpinKPointIndex + 1

!            occupancyNumber = electronPopulation(stateSpinKPointIndex)

            ! Determine the occupancy number for this state.  The default
            !   assumption is that each state contains 2 electrons.  This is
            !   because the sum of all kpoint weighting factors is equal to 2.
            !   This will make spin polarized states half occupied and spin
            !   non-polarized states fully occupied.
            if (energyEigenValues(j,i,h) <= 0.0_double) then
               occupancyNumber = 1.0_double/real(spin,double)
            else
               occupancyNumber = 0.0_double
            endif

            ! Initialize the accumlator for the wave function--overlap
            !   interaction.
            pdosAccum(:) = 0.0_double

            ! Initialize a counter that tracks the index of valeDim (total
            !   number of states used in the system).
            valeDimIndex = 0

            ! Begin a loop over all atoms
            do k = 1, numAtomSites

               ! Loop over all the valence states for this atom
               do l = 1, numAtomStates(k)

                  ! Increment the valeDimIndex
                  valeDimIndex = valeDimIndex + 1

                  ! Read the Mulliken projection for this basis
                  !   function and state, formed above by the single
                  !   matrix product (PSEUDOCODE 34.5).  This value
                  !   equals the direct per-column overlap sum to the
                  !   reassociation floor (PSEUDOCODE 34.2, 34.7), so
                  !   the three accumulations that follow are unchanged.
                  oneValeRealAccum = mullikenProj(valeDimIndex, j)

                  ! Store the current electron number assignment.
                  electronNumber(l,k) = electronNumber(l,k) + &
                        & oneValeRealAccum * kPointWeight(i) * occupancyNumber
!                  electronNumber(l,k) = electronNumber(l,k) + &
!                        & oneValeRealAccum * occupancyNumber


                  pdosAccum(pdosIndex(valeDimIndex)) = &
                        & pdosAccum(pdosIndex(valeDimIndex)) + &
                        & oneValeRealAccum / real(spin,double)


                  ! Store the square of the accumulation as the localization
                  !   index for this state.
                  localizationIndex(j) = localizationIndex(j) + &
                        & oneValeRealAccum * oneValeRealAccum * &
                        & kPointWeight(i) / real(spin,double)
!                  localizationIndex(j) = localizationIndex(j) + &
!                        & oneValeRealAccum * oneValeRealAccum * & &
!                        occupancyNumber
               enddo
            enddo


            ! Apply broadening to the pdos values according to the given
            !   pdos sigma factor.
            do k = 1, numEnergyPoints

               ! Compute the exponential term for the broadening of this point.
               !   Note that from the usual Gaussian term of exp(-alpha*x^2)
               !   we are have (eV-eS)^2 = x^2 and (1/sigma)^2 = alpha.
               expTerm = ((energyEigenValues(j,i,h) - &
                     & energyScale(k))/sigmaDOS)**2

               ! If the exponential term is less than 50 we apply the
               !   broadening.
               if (expTerm < 50.0_double) then

                  ! Compute the exponential factor.  It is at this point that
                  !   the kpoint weighting factor is applied to make the energy
                  !   bucket we are considering be filled with the right number
                  !   of electrons.  The pdosAccum below already has the
                  !   1/spin factor included.  Note that the hartree conversion
                  !   must be included here because the sigmaDOS (that was used
                  !   to make sigmaSqrtPi) is in units of hartree and the
                  !   exponential and kPointWeights are unitless. Thus, if we
                  !   want to get States / [eV Cell] then we need to convert
                  !   the sigmaSqrtPi back to eV. Because it is in the
                  !   denominator we need to divide the whole expression by
                  !   eV/hartree to get a final energy unit of eV. Recall that
                  !   the hartree variable equals 27.211... eV/hartree.
                  expFactor = exp(-expTerm) / sigmaSqrtPi / hartree * &
                        & kpointWeight(i)

                  ! Broaden and store the complete pdos
                  pdosComplete(1:cumulDOSTotal,k) = &
                        & pdosComplete(1:cumulDOSTotal,k) &
                        & + pdosAccum(1:cumulDOSTotal) * expFactor
               endif
            enddo
         enddo ! (j numStates)

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

      enddo ! (i numKPoints)

      ! Reset the line pointer to the beginning of the next line if necessary.
      if (mod(numKPoints,10) .ne. 0) then
         write (20,*)
         call flush (20)
      endif

      ! Unfold the Gaussian per-atom PDOS over the star of each irreducible
      !   k-point (PSEUDOCODE 23, TODO C148).
      !
      ! This pathway never visits a star member: it accumulates over
      !   irreducible points with the multiplicity carried inside
      !   kPointWeight, so every member of a star was credited with the
      !   REPRESENTATIVE's atom-by-atom breakdown. On cubic KNbO3 that left
      !   three symmetry-equivalent oxygens differing by seventy percent.
      !
      ! One group average of the finished accumulation is EXACTLY the
      !   per-k-point star average, and PSEUDOCODE 23.1 derives why: the
      !   projection at an irreducible point is invariant under that point's
      !   little group, so a star sum equals the whole-group average times
      !   the multiplicity, and the group average commutes with the sum over
      !   k-points because it is linear. So no loop moves and nothing inside
      !   the accumulation changes.
      !
      ! Unlike symmetrizePDOS_LAT below, this is NOT optional cleanup. There
      !   the per-corner permutation already gives the right answer and the
      !   averaging only makes equal things exactly equal. Here the average
      !   IS the correctness mechanism, so it is deliberately not wired to
      !   symmetrizeLATPartials -- that switch would turn off a fix rather
      !   than a polish.
      !
      ! Mode 3 is NOT covered, because the m components of an orbital mix
      !   under a rotation and a channel permutation cannot express that.
      !   It needs the D^l(R) matrices, which is why DESIGN 1.4 refuses mode
      !   3 on the tetrahedron path. The Gaussian path does not refuse it, so
      !   mode 3 on a reduced mesh stays wrong here (PSEUDOCODE 23.4).
      ! No test on kPointIntgCode: this sits inside the Gaussian arm of that
      !   branch already. The table exists only for the atom-resolved modes,
      !   so asking whether it was built also selects the modes that need it.
      if (allocated(channelPermTable)) then
         call symmetrizePDOS_Gaussian(pdosComplete, channelPermTable, &
               & cumulDOSTotal, numEnergyPoints)
      endif

      endif ! (kPointIntgCode branch)


      ! Compute the total DOS for the whole system.
      totalSystemDos(:) = 0.0_double
      do i = 1, numEnergyPoints
         totalSystemDos(i) = totalSystemDos(i) + sum(pdosComplete(:,i))
      enddo


      ! Compute the number of electron states present in the range of eminDOS
      !   to emaxDOS for this spin orientation or combination.
      numStatesInRange = 0.0_double
      do i = 1, numKPoints
         do j = 1, numStates
            if ((energyEigenValues(j,i,h) > eminDOS) .and. &
                  & (energyEigenValues(j,i,h) < emaxDOS)) then
               numStatesInRange = numStatesInRange + kPointWeight(i) / &
                     & real(spin,double)
            endif
         enddo
      enddo


      ! Compute the total number of electrons determined from the pdos
      !   calculation.
      totalElectronsComputed = sum(electronNumber(:,:))


      ! Compute the total number of electrons in each spin direction (up,down).
      !   This number is then be compared to the value for the
      !   totalElectronsComputed to normalize the spectra.
      currentPopulation = 0.0_double
      energyLevelCounter = 0
      do i = 1, numKPoints
         do j = 1, spin
            do k = 1, numStates
               energyLevelCounter = energyLevelCounter + 1
               ! Accumulate only electrons for the current spin orientation.
               !   In the case of spin non-polarized calculations we will
               !   accumulate all populations.
               if (j == h) then
                  currentPopulation = currentPopulation + &
                        & electronPopulation(energyLevelCounter)
               endif
            enddo
         enddo
      enddo


      ! Compute the ratio of the exact electron count to the computed electron
      !   count. For the Gaussian path this corrects broadening-tail truncation.
      !   For LAT, tetrahedron corner weights provide exact BZ integration so
      !   electronFactor should be ~1.0 (DESIGN 1.4).
      electronFactor = currentPopulation &
            & / totalElectronsComputed

      if (kPointIntgCode == 1) then
         ! LAT path: log the ratio as a diagnostic but do NOT apply it. A ratio
         !   far from 1.0 signals an integration bug.
         write (20, *) 'LAT electronFactor: ', &
               & electronFactor
      else
         ! Gaussian path: apply the correction to electronNumber, TDOS, and
         !   PDOS.
         electronNumber(:,:) = &
               & electronNumber(:,:) * electronFactor
         totalSystemDos(:) = &
               & totalSystemDos(:) * electronFactor
         pdosComplete(:,:) = &
               & pdosComplete(:,:) * electronFactor
      endif


      ! Find the total area under the computed TDOS. This should equal the
      !   number of electron spin states available in the system (occupied +
      !   unoccupied) OVER THE REQUESTED ENERGY RANGE ONLY.
      integratedArea = sum( &
            & totalSystemDos(1:numEnergyPoints - 1) &
            & + totalSystemDos(2:numEnergyPoints)) &
            & * deltaDOS * 0.5_double * hartree

      ! Record the exact and calculated values for electrons and states.
      write (20, *) 'Electrons Calculated:    ', &
            & totalElectronsComputed
      write (20, *) 'Electrons Expected:      ', &
            & currentPopulation
      write (20, *) 'Spin States Calculated:  ', &
            & integratedArea
      write (20, *) 'Spin States Expected:    ', &
            & numStatesInRange


      ! Begin recording the results to disk.

      ! Record the total system DOS, converting the scale to eV. For the LAT
      !   path, TDOS has already been written by computeTDOS_LAT to fort.60/61,
      !   so we skip the TDOS write here.
      if (kPointIntgCode /= 1) then
         write (59+h, fmt="(a14,a14)") &
               & "    ENERGY(eV)", "          TDOS"
         do i = 1, numEnergyPoints
            write (59+h, fmt="(f14.4,f14.6)") &
                  & energyScale(i) * hartree, &
                  & totalSystemDos(i)
         enddo
      endif

      ! Loop over types for the types-based DOS
      if (detailCodePDOS == 0) then

         ! Print the key bits of information for the PDOS output.
         write (69+h,fmt="(a7)") 'STYLE 1'
         write (69+h,fmt="(a10,i6)") 'NUM_UNITS ',numAtomTypes
         write (69+h,fmt="(a11,i9)") 'NUM_POINTS ',numEnergyPoints

         ! Print the energy scale used by all types, converting to eV.
         do i = 1, numEnergyPoints
            write (69+h,fmt="(f16.8)") energyScale(i) * hartree
         enddo

         do i = 1, numAtomTypes

            ! Determine the indices of the pdosComplete to use for this type,
            !   beginning with the first orbital and ending with the last for
            !   this type.
            initIndex = cumulNumDOS(i) + 1
            finIndex  = cumulNumDOS(i+1)


            ! Write the label and orbital information for this type.
            write (69+h,fmt="(a13,i5)") 'SEQUENCE_NUM ',i
            write (69+h,fmt="(a13,a3)") 'ELEMENT_NAME ',atomTypes(i)%elementName
            write (69+h,fmt="(a11,i5)") 'SPECIES_ID ',atomTypes(i)%speciesID
            write (69+h,fmt="(a8,i5)") 'TYPE_ID ',atomTypes(i)%typeID
            write (69+h,fmt="(a10,1x,i2)") 'COL_LABELS',finIndex-initIndex+2
            write (69+h,ADVANCE='NO',fmt="(a6)") 'TOTAL '
            numCols = 1
            do j = 1, lAngMomCount  ! 1=s; 2=p; 3=d; 4=f
               do k = 1, atomTypes(i)%numQN_lValeRadialFns(j)

                  numCols = numCols + 1

                  ! Write the orbital definition in the standard 1s,2s,3s,2p,3p
                  !   notation.  The number is the n quantum number and we
                  !   exclude the orbitals that are in the core and begin
                  !   counting with the orbitals in the valence.  Note that the
                  !   s orbitals ascend like 1s, 2s, 3s while the p orbitals
                  !   ascend like 2p, 3p.  The (+j-1) accounts for this need.
                  write (69+h,ADVANCE='NO',fmt="(i1,a1,1x)") &
                        & atomTypes(i)%numQN_lCoreRadialFns(j)+k+j-1, &
                        & QN_lLetter(j)

                  if (mod(numCols,6) == 0) then
                     write (69+h,*)
                  endif
               enddo
            enddo

            ! Advance the file pointer if the data fit was not exact.
            if (mod(numCols,6) /= 0) then
               write (69+h,*)
            endif

            ! Write the raw data for this unit.  The first value is the sum of
            !   all orbitals for this type.  (TDOS for this type.)  The next
            !   values are the QN_nl resolved PDOS for this type.
            do j = 1, numEnergyPoints
               write (69+h,ADVANCE="NO",fmt="(e13.6,1x)") &
                  & sum(pdosComplete(initIndex:finIndex,j))
               numCols = 1
               do k = initIndex,finIndex
                  numCols = numCols + 1
                  write (69+h,ADVANCE="NO",fmt="(e13.6,1x)") pdosComplete(k,j)
                  if (mod(numCols,6) == 0) then
                    write (69+h,*)
                  endif
               enddo

               ! Advance the file pointer if the data fit was not exact.
               if (mod(numCols,6) /= 0) then
                  write (69+h,*)
               endif
            enddo
         enddo
      ! Loop over atoms for the atom-based DOS
      elseif (detailCodePDOS == 1) then

         ! Print the key bits of information for the PDOS output.
         write (69+h,fmt="(a7)") 'STYLE 1'
         write (69+h,fmt="(a10,i6)") 'NUM_UNITS ', numAtomSites
         write (69+h,fmt="(a11,i9)") 'NUM_POINTS ', numEnergyPoints

         ! Print the energy scale used by all atoms, converting to eV.
         do i = 1, numEnergyPoints
            write (69+h,fmt="(f16.8)") energyScale(i) * hartree
         enddo

         do i = 1, numAtomSites

            ! Obtain the type of the current atom.
            currentType = atomSites(i)%atomTypeAssn

            ! Write the label information for this atom.
            write (69+h,fmt="(a13,i5,1x)") 'SEQUENCE_NUM ',i
            write (69+h,fmt="(a13,a3)") 'ELEMENT_NAME ',&
                  & atomTypes(currentType)%elementName
            write (69+h,fmt="(a11,i5)") 'SPECIES_ID ',&
                  & atomTypes(currentType)%speciesID
            write (69+h,fmt="(a8,i5)") 'TYPE_ID ',&
                  & atomTypes(currentType)%typeID
            write (69+h,fmt="(a10,1x,i2)") 'COL_LABELS',1
            write (69+h,fmt="(a6)") 'TOTAL '

            ! Write the raw data for this unit.
            do j = 1, numEnergyPoints
               write (69+h,fmt="(f16.8)") pdosComplete(i,j)
            enddo
         enddo
      ! Loop over atoms & orbitals (valeDim).
      elseif (detailCodePDOS == 2) then

         ! Print the key bits of information for the PDOS output.
         write (69+h,fmt="(a7)") 'STYLE 1'
         write (69+h,fmt="(a10,i6)") 'NUM_UNITS ',numAtomSites
         write (69+h,fmt="(a11,i9)") 'NUM_POINTS ', numEnergyPoints

         ! Print the energy scale used by all atoms, converting to eV.
         do i = 1, numEnergyPoints
            write (69+h,fmt="(f16.8)") energyScale(i) * hartree
         enddo

         ! Initialize a record of which valence index the loop is on.
         valeDimIndex = 0

         do i = 1, numAtomSites

            ! Obtain the type of the current atom.
            currentType = atomSites(i)%atomTypeAssn

            ! Determine the indices of the pdosComplete to use for this atom,
            !   beginning with the first orbital and ending with the last for
            !   this atom.
            initIndex = cumulNumDOS(i)+1
            finIndex  = cumulNumDOS(i+1)


            ! Write the label information for this atom.
            write (69+h,fmt="(a13,i5,1x)") 'SEQUENCE_NUM ',i
            write (69+h,fmt="(a13,a3)") 'ELEMENT_NAME ',&
                  & atomTypes(currentType)%elementName
            write (69+h,fmt="(a11,i5)") 'SPECIES_ID ',&
                  & atomTypes(currentType)%speciesID
            write (69+h,fmt="(a8,i5)") 'TYPE_ID ',&
                  & atomTypes(currentType)%typeID
            write (69+h,fmt="(a10,1x,i2)") 'COL_LABELS',finIndex-initIndex+2
            write (69+h,ADVANCE='NO',fmt="(a6)") 'TOTAL '
            numCols = 1
            do j = 1, lAngMomCount
               do k = 1, atomTypes(currentType)%numQN_lValeRadialFns(j)

                  numCols = numCols + 1

                  ! Write the orbital definition in the standard 1s,2s,3s,2p,3p
                  !   notation.  The number is the n quantum number and we
                  !   exclude the orbitals that are in the core and begin
                  !   counting with the orbitals in the valence.  Note that the
                  !   s orbitals ascend like 1s, 2s, 3s while the p orbitals
                  !   ascend like 2p, 3p.  The (+j-1) accounts for this need.
                  write (69+h,ADVANCE='NO',fmt="(i1,a1,1x)") &
                        & atomTypes(currentType)%numQN_lCoreRadialFns(j)+k+j-1,&
                        & QN_lLetter(j)

                  if (mod(numCols,6) == 0) then
                     write (69+h,*)
                  endif
               enddo
            enddo

            ! Advance the file pointer if the data fit was not exact.
            if (mod(numCols,6) /= 0) then
               write (69+h,*)
            endif

            ! Write the raw data for this unit.  The first value is the TDOS
            !   for this atom.  The other values are the QN_nl resolved PDOS
            !   for this atom.
            do j = 1, numEnergyPoints
               write (69+h,ADVANCE="NO",fmt="(e13.6,1x)") &
                     & sum(pdosComplete(initIndex:finIndex,j))
               numCols = 1
               do k = initIndex,finIndex
                  numCols = numCols + 1
                  write (69+h,ADVANCE="NO",fmt="(e13.6,1x)") pdosComplete(k,j)
                  if (mod(numCols,6) == 0) then
                    write (69+h,*)
                  endif
               enddo

               ! Advance the file pointer if the data fit was not exact.
               if (mod(numCols,6) /= 0) then
                  write (69+h,*)
               endif
            enddo
         enddo
      ! Loop over atoms and QN_nlm orbitals.
      elseif (detailCodePDOS == 3) then

         ! Print the key bits of information for the PDOS output.
         write (69+h,fmt="(a7)") 'STYLE 1'
         write (69+h,fmt="(a10,i6)") 'NUM_UNITS ',numAtomSites
         write (69+h,fmt="(a11,i9)") 'NUM_POINTS ', numEnergyPoints

         ! Print the energy scale used by all atoms, converting to eV.
         do i = 1, numEnergyPoints
            write (69+h,fmt="(f16.8)") energyScale(i) * hartree
         enddo

         ! Initialize a record of which valence index the loop is on.
         valeDimIndex = 0

         do i = 1, numAtomSites

            ! Obtain the type of the current atom.
            currentType = atomSites(i)%atomTypeAssn

            ! Determine the indices of the pdosComplete to use for this atom,
            !   beginning with the first orbital and ending with the last for
            !   this atom.
            initIndex = cumulNumDOS(i)+1
            finIndex  = cumulNumDOS(i+1)


            ! Write the label information for this atom.
            write (69+h,fmt="(a13,i5,1x)") 'SEQUENCE_NUM ',i
            write (69+h,fmt="(a13,a3)") 'ELEMENT_NAME ',&
                  & atomTypes(currentType)%elementName
            write (69+h,fmt="(a11,i5)") 'SPECIES_ID ',&
                  & atomTypes(currentType)%speciesID
            write (69+h,fmt="(a8,i5)") 'TYPE_ID ',&
                  & atomTypes(currentType)%typeID
            write (69+h,fmt="(a10,1x,i2)") 'COL_LABELS',finIndex-initIndex+2
            write (69+h,ADVANCE='NO',fmt="(a6)") 'TOTAL '
            numCols = 1
            do j = 1, lAngMomCount
               do k = 1, atomTypes(currentType)%numQN_lValeRadialFns(j)
                  do l = 1, (j-1)*2+1

                     numCols = numCols + 1

                     ! Write the QN_nlm orbital definition in the following
                     !   notation:  1s,2s,3s,2px,2py,2pz,3px,3py,3pz,...  The
                     !   first number is the n QN (excluding all orbitals from
                     !   the core).  The first letter is the l QN, and the
                     !   string after that is the m QN in xyz notation.  Note
                     !   that the s orbitals ascend like 1s, 2s, 3s, while the
                     !   p orbitals ascend like 2p, 3p.  The (+j-1) accounts
                     !   for this need.
                     write (formatString,fmt="(a11,i2.2,a4)") "(i1,a1,a1,a", &
                           & len_trim(QN_mLetter(j,l)),",1x)"

                     write (69+h,ADVANCE='NO',fmt=formatString) &
                           & atomTypes(currentType)%&
                           & numQN_lCoreRadialFns(j)+k+j-1, &
                           & QN_lLetter(j),"_",QN_mLetter(j,l)

                     if (mod(numCols,6) == 0) then
                        write (69+h,*)
                     endif
                  enddo
               enddo
            enddo

            ! Advance the file pointer if the data fit was not exact.
            if (mod(numCols,6) /= 0) then
               write (69+h,*)
            endif

            ! Write the raw data for this unit.  The first value is the TDOS
            !   for this atom.  The other values are the QN_nlm resolved PDOS
            !   for this atom.
            do j = 1, numEnergyPoints
               write (69+h,ADVANCE="NO",fmt="(e13.6,1x)") &
                     & sum(pdosComplete(initIndex:finIndex,j))
               numCols = 1
               do k = initIndex,finIndex
                  numCols = numCols + 1
                  write (69+h,ADVANCE="NO",fmt="(e13.6,1x)") pdosComplete(k,j)
                  if (mod(numCols,6) == 0) then
                    write (69+h,*)
                  endif
               enddo

               ! Advance the file pointer if the data fit was not exact.
               if (mod(numCols,6) /= 0) then
                  write (69+h,*)
               endif
            enddo
         enddo
      endif


      ! Compute the weighted average of the energy of a given band across all
      !   kpoints. This energy value is used as the energy for the localization
      !   index.

      ! Initialize the averageEnergy (stored in energyValuesAvg).
      energyValuesAvg(:) = 0.0_double

      do i = 1, numKPoints

         ! Collect weighted sum of the energy for this kpoint over all states.
         energyValuesAvg(:) = energyValuesAvg(:) + &
               & energyEigenValues(:numStates,i,h) * kPointWeight(i)
      enddo

      ! Divide the result by two since the sum of kpoint weights is two and not
      !   one.
      energyValuesAvg(:) = energyValuesAvg(:) / 2.0_double
      localizationIndex(:) = localizationIndex(:) / 2.0_double

      ! Record the results to disk as the localization index for each state,
      !   converting the energy scale to eV.
      do i = 1, numStates
         write (79+h,fmt="(1x,f24.8,2x,f24.8)") energyValuesAvg(i) * hartree,&
               & localizationIndex(i)
      enddo
   enddo ! (h spin)

   ! Deallocate all the unnecessary matrices and arrays. The Gaussian-only
   !   arrays (pdosAccum, the projection matrices, overlap) are guarded because
   !   the LAT path handles its own allocations.
   deallocate (cumulNumDOS)
   deallocate (pdosIndex)
   deallocate (numAtomStates)
   deallocate (energyScale)
   deallocate (localizationIndex)
   deallocate (totalSystemDos)
   deallocate (pdosComplete)
   deallocate (electronNumber)

   ! Gaussian-path-only deallocations.
   if (kPointIntgCode /= 1) then
      deallocate (pdosAccum)
      deallocate (statesProjected)
      deallocate (mullikenProj)
#ifndef GAMMA
      if (inSCF == 0) then
         deallocate (valeValeOL)
      endif
#else
      if (inSCF == 0) then
         deallocate (valeValeOLGamma)
      endif
#endif
   endif ! (kPointIntgCode /= 1)

   ! The channel permutation table was built before the spin loop, on whichever
   !   pathway needed it. Released on the SAME condition it was built under,
   !   expressed by asking the array rather than by repeating the condition:
   !   widening one guard and not the other either leaks the table or frees
   !   something that was never allocated, and a permutation table whose
   !   lifetime disagreed with its consumer is exactly what made O9's fix
   !   silently do nothing (PSEUDOCODE 19.2.1a, 23.2).
   if (allocated(channelPermTable)) then
      deallocate (channelPermTable)
   endif

   ! Log the date and time we end.
   call timeStampEnd (19)

end subroutine computeDOS


! Compute the total density of states using the Linear Analytic Tetrahedron
!   (LAT) method (Bloechl, Jepsen, & Andersen, PRB 49, 16223, 1994). This
!   approach decomposes the Brillouin zone into tetrahedra and integrates the
!   DOS analytically within each one, eliminating the need for a broadening
!   parameter.
!
!   The algorithm loops over bands and tetrahedra. For each tetrahedron, the
!   four corner eigenvalues are looked up via the fullKPToIBZKPMap (which maps
!   full-mesh kpoint indices to IBZ kpoint indices), sorted, and the Bloechl
!   analytic formulas are applied to each energy grid point.
!
!   Output is written to the same TDOS file (unit 60/61) in the same format as
!   the Gaussian broadening path. PDOS is not computed in this subroutine
!   (future work).
subroutine computeTDOS_LAT

   ! Import the necessary modules.
   use O_Kinds
   use O_MathSubs, only: bloechlCornerDOSWt
   use O_TimeStamps
   use O_Potential, only: spin
   use O_KPoints, only: numKPoints, numTetrahedra, &
         & tetraVol, tetrahedra, &
         & fullKPToIBZKPMap, kPointWeight
   use O_Constants, only: hartree
   use O_Input, only: numStates, eminDOS, emaxDOS, &
         & deltaDOS
   use O_SecularEquation, only: energyEigenValues

   ! Make sure that no funny variables are used.
   implicit none

   ! Define local variables.
   integer :: h         ! Spin loop index.
   integer :: n         ! Band (state) loop index.
   integer :: t         ! Tetrahedron loop index.
   integer :: iE        ! Energy grid loop index.
   integer :: corner    ! Corner loop index.
   integer :: ibzK      ! IBZ kpoint index for corner.
   integer :: minIdx    ! Index of minimum for sort.
   real (kind=double) :: energy ! Current grid point.
   real (kind=double) :: integratedArea
   real (kind=double) :: tempVal
   real (kind=double) :: kpWtSum ! sum(kPointWeight)
   real (kind=double), allocatable, dimension(:) :: &
         & totalSystemDos
   real (kind=double), dimension(4) :: eps
   real (kind=double), dimension(4) :: cornerDOSWt

   ! Log the date and time we start.
   call timeStampStart (19)

   ! Determine the number of energy grid points.
   numEnergyPoints = int((emaxDOS - eminDOS) / deltaDOS)

   ! Allocate and fill the energy scale.
   if (allocated(energyScale)) deallocate(energyScale)
   allocate (energyScale(numEnergyPoints))
   do iE = 1, numEnergyPoints
      energyScale(iE) = eminDOS + (iE - 1) * deltaDOS
   enddo

   ! Allocate the TDOS accumulation array.
   allocate (totalSystemDos(numEnergyPoints))

   ! Normalization factor: the tetrahedron BZ integration sums tetraVol to 1.0
   !   over the full BZ, but the Gaussian path uses kPointWeight (which sums to
   !   2.0) as its BZ integration weight. The factor of 2 in kPointWeight
   !   accounts for the two electron spin states per band in a spin-unpolarized
   !   calculation (spin=1). Both paths divide by spin separately, so to match
   !   conventions the LAT path must include this same weight sum. For spin=1:
   !   2/1 = 2 spin states per band. For spin=2: 2/2 = 1 spin state per band.
   kpWtSum = sum(kPointWeight(1:numKPoints))

   ! Loop over spin orientations.
   do h = 1, spin

      ! Initialize the TDOS array.
      totalSystemDos(:) = 0.0_double

      ! Record which calculation is being done.
      if (spin == 2) then
         if (h == 1) then
            write (20,*) "Computing spin up LAT TDOS."
         else
            write (20,*) "Computing spin dn LAT TDOS."
         endif
      endif

      ! Loop over bands (states).
      do n = 1, numStates

         ! Loop over tetrahedra.
         do t = 1, numTetrahedra

            ! Look up the four corner kpoint indices from the tetrahedra array
            !   (full mesh indices) and map them to IBZ kpoint indices.
            do corner = 1, 4
               ibzK = fullKPToIBZKPMap(tetrahedra(corner, t))
               eps(corner) = &
                     & energyEigenValues(n, ibzK, h)
            enddo

            ! Sort the four eigenvalues in ascending order using a simple
            !   selection sort (4 elements).
            do corner = 1, 3
               minIdx = corner
               do ibzK = corner + 1, 4
                  if (eps(ibzK) < eps(minIdx)) then
                     minIdx = ibzK
                  endif
               enddo
               if (minIdx /= corner) then
                  tempVal = eps(corner)
                  eps(corner) = eps(minIdx)
                  eps(minIdx) = tempVal
               endif
            enddo

            ! Loop over the energy grid and accumulate the DOS contribution from
            !   this tet.
            do iE = 1, numEnergyPoints
               energy = energyScale(iE)

               ! Skip if outside eigenvalue range.
               if (energy < eps(1) .or. &
                     & energy >= eps(4)) cycle

               ! Compute per-corner DOS weights. The TDOS uses only their sum.
               call bloechlCornerDOSWt( &
                     & energy, eps, cornerDOSWt)

               ! Accumulate into the TDOS. tetraVol is the BZ fraction for one
               !   tet. kpWtSum matches the kPointWeight convention (sums to 2).
               !   1/spin accounts for spin degeneracy. 1/hartree converts
               !   states/Hartree to states/eV.
               totalSystemDos(iE) = &
                     & totalSystemDos(iE) &
                     & + sum(cornerDOSWt) &
                     & * tetraVol * kpWtSum &
                     & / real(spin, double) &
                     & / hartree

            enddo ! iE (energy grid)
         enddo ! t (tetrahedra)

         ! Progress indicator.
         if (mod(n, 50) == 0) then
            write (20, ADVANCE="NO", FMT="(a1)") "."
            call flush(20)
         endif

      enddo ! n (states)

      write (20, *) ""

      ! Compute the integrated area under the TDOS curve using the trapezoidal
      !   rule.
      integratedArea = sum( &
            & totalSystemDos(1:numEnergyPoints - 1) &
            & + totalSystemDos(2:numEnergyPoints)) &
            & * deltaDOS * 0.5_double * hartree
      write (20, *) "LAT TDOS integrated area: ", &
            & integratedArea

      ! Write the TDOS to file (same format as Gaussian).
      write (59+h, fmt="(a14,a14)") &
            & "    ENERGY(eV)", "          TDOS"
      do iE = 1, numEnergyPoints
         write (59+h, fmt="(f14.4,f14.6)") &
               & energyScale(iE) * hartree, &
               & totalSystemDos(iE)
      enddo

   enddo ! h (spin)

   ! Clean up. Deallocate both the local TDOS array and the module-level
   !   energyScale so that computeDOS can re-allocate it independently.
   deallocate (totalSystemDos)
   deallocate (energyScale)

   ! Log the date and time we end.
   call timeStampEnd (19)

end subroutine computeTDOS_LAT



! Build the channel permutation lookup table for LAT PDOS IBZ unfolding
!   (PSEUDOCODE 8.1). For each point group operation R and PDOS channel alpha,
!   channelPermTable(R, alpha) gives the channel index at the IBZ k-point whose
!   projection should be used for channel alpha at the full-mesh k-point related
!   by R.
!
!   The permutation depends on the PDOS detail mode: Mode 0 (per-type, per-l):
!     identity. Type-level
!       sums are invariant under R because R permutes atoms within each type.
!     Mode 1 (per-atom total): channel = atom index, so permute via invAtomPerm.
!     Mode 2 (per-atom, per-l): remap the atom index via invAtomPerm while
!       preserving the l-shell offset within the atom (same species => same
!       orbital structure).
!     Mode 3 (per-atom, per-lm): not supported with LAT (requires D^l rotation
!       matrices).
!
!   Built once before the spin loop since channel permutation is
!   spin-independent.
subroutine buildChannelPermTable(detailCode, &
      & numOps, totalChannels, cumulDOS, &
      & numSites, invPerm, channelPermTable)

   use O_Kinds

   implicit none

   ! Passed parameters.
   integer, intent(in) :: detailCode
   integer, intent(in) :: numOps
   integer, intent(in) :: totalChannels
   integer, dimension(:), intent(in) :: cumulDOS
   integer, intent(in) :: numSites
   integer, dimension(:,:), intent(in) :: invPerm
   integer, allocatable, dimension(:,:), &
         & intent(out) :: channelPermTable

   ! Local variables.
   integer :: opR       ! Point group operation index.
   integer :: alpha     ! Channel loop index.
   integer :: atomA     ! Atom index decoded from alpha.
   integer :: permAtom  ! Permuted atom index.
   integer :: baseOld   ! Cumulative offset for atomA.
   integer :: baseNew   ! Cumulative offset for permAtom.
   integer :: nOrbitals ! Number of l-shells for atomA.
   integer :: orbOff    ! Orbital offset loop index.

   allocate (channelPermTable(numOps, totalChannels))

   ! Mode 0: per-type, per-l. Type-level sums are invariant under R (R permutes
   !   atoms within each type, so the type sum is unchanged). Channel
   !   permutation is the identity.
   if (detailCode == 0) then
      do opR = 1, numOps
         do alpha = 1, totalChannels
            channelPermTable(opR, alpha) = alpha
         enddo
      enddo
      return
   endif

   ! Mode 1: per-atom total. Channel index = atom index. Permute directly via
   !   invAtomPerm.
   if (detailCode == 1) then
      do opR = 1, numOps
         do atomA = 1, numSites
            channelPermTable(opR, atomA) = &
                  & invPerm(opR, atomA)
         enddo
      enddo
      return
   endif

   ! Mode 2: per-atom, per-l. Decode alpha into (atom, l-offset), permute the
   !   atom via invAtomPerm, re-encode using the permuted atom's cumulative
   !   offset. The l-shell offset is unchanged because same species implies
   !   identical orbital structure.
   if (detailCode == 2) then
      do opR = 1, numOps
         do atomA = 1, numSites
            permAtom = invPerm(opR, atomA)
            baseOld = cumulDOS(atomA)
            baseNew = cumulDOS(permAtom)
            nOrbitals = cumulDOS(atomA + 1) &
                  & - cumulDOS(atomA)
            do orbOff = 1, nOrbitals
               channelPermTable(opR, &
                     & baseOld + orbOff) = &
                     & baseNew + orbOff
            enddo
         enddo
      enddo
      return
   endif

end subroutine buildChannelPermTable


! Pass 1 of the LAT PDOS computation (PSEUDOCODE 8.2). Stream through IBZ
!   k-points, read eigenvectors and overlap from HDF5, compute Mulliken
!   projections for all bands and channels, and store into projArray.
!
!   The Mulliken decomposition is identical to the Gaussian path: the states
!   are carried through the overlap in one matrix product, T = S C, and the
!   projection P(mu,j) = Re(conjg(C(mu,j)) T(mu,j)) gives the weight of basis
!   function mu in band j (PSEUDOCODE 34). These are summed by channel
!   according to pdosIndex.
!
!   In addition, electronNumber and localizationIndex are accumulated using the
!   same formulas as the Gaussian path (kPointWeight * step-function occupancy)
!   to provide the normalization diagnostic.
!
!   projArray is allocated here and must be deallocated by the caller after
!   integratePDOS_LAT completes.
subroutine computeProjections_LAT(inSCF, spinIdx, &
      & numKP, numSt, numSites, numAtmSt, pIndex, &
      & vDim, cumDOSTotal, spinCount, projArr, &
      & elecNum, locIdx)

   use O_Kinds
   use O_KPoints, only: kPointWeight
   use O_SecularEquation, only: energyEigenValues, &
         & readDataSCF, readDataPSCF
#ifndef GAMMA
   use O_SecularEquation, only: valeValeOL, valeVale
#else
   use O_SecularEquation, only: valeValeOLGamma, &
         & valeValeGamma
#endif
   ! The shared state-projection routine (PSEUDOCODE 34), the same
   !   T = S C used on the Gaussian path in computeDOS.
   use O_StateProjection, only: projectStatesOntoBasis

   implicit none

   ! Passed parameters.
   integer, intent(in) :: inSCF
   integer, intent(in) :: spinIdx
   integer, intent(in) :: numKP
   integer, intent(in) :: numSt
   integer, intent(in) :: numSites
   integer, dimension(:), intent(in) :: numAtmSt
   integer, dimension(:), intent(in) :: pIndex
   integer, intent(in) :: vDim
   integer, intent(in) :: cumDOSTotal
   integer, intent(in) :: spinCount
   real (kind=double), allocatable, &
         & dimension(:,:,:), intent(out) :: projArr
   real (kind=double), dimension(:,:), &
         & intent(inout) :: elecNum
   real (kind=double), dimension(:), &
         & intent(inout) :: locIdx

   ! Local variables.
   integer :: i, j, k, l       ! Loop indices.
   integer :: valeDimIdx        ! Tracks position in the full valence dimension.
   real (kind=double) :: occupNum  ! Step-function occupancy for normalization
         !   diagnostic.
   real (kind=double) :: oneValeRA ! Mulliken projection for one basis function.
   ! The state projection T = S C and the Mulliken projection
   !   P(mu,j) = Re(conjg(C(mu,j)) T(mu,j)), formed once per k-point by
   !   one matrix product (PSEUDOCODE 34), exactly as on the Gaussian
   !   path in computeDOS.  T is complex in the multi-k build and real
   !   in the gamma build; P is real in both.
#ifndef GAMMA
   complex (kind=double), allocatable, &
         & dimension (:,:) :: statesProjected
#else
   real (kind=double), allocatable, &
         & dimension (:,:) :: statesProjected
#endif
   real (kind=double), allocatable, &
         & dimension (:,:) :: mullikenProj

   ! Allocate the projection array: (channel, band, IBZ kpoint). This is the
   !   main memory cost of the LAT PDOS computation.
   allocate (projArr(cumDOSTotal, numSt, numKP))
   projArr(:,:,:) = 0.0_double

   ! Allocate work arrays for the Mulliken product.  The eigenvector
   !   array is one slab: this pass runs for one spin and the reader
   !   delivers that spin's vectors to slab 1 on request (DESIGN 2.8).
   allocate (statesProjected(vDim, numSt))
   allocate (mullikenProj   (vDim, numSt))
#ifndef GAMMA
   allocate (valeValeOL(vDim, vDim))
   if (inSCF == 0) then
      allocate (valeVale(vDim, numSt, 1))
   endif
#else
   allocate (valeValeOLGamma(vDim, vDim))
   if (inSCF == 0) then
      allocate (valeValeGamma(vDim, numSt, 1))
   endif
#endif

   ! Stream through IBZ k-points, reading eigenvectors and overlap one k-point
   !   at a time.
   do i = 1, numKP

      ! Read eigenvectors + overlap for this IBZ kpoint and spin orientation.
      !   The eigenvectors are delivered to slab 1 (slab = 1): this pass
      !   holds a single slab and reads slab 1 below (DESIGN 2.8,
      !   PSEUDOCODE 33).
      if (inSCF == 1) then
         call readDataSCF(spinIdx, i, numSt, 1, slab=1)
      else
         call readDataPSCF(spinIdx, i, numSt, 1, slab=1)
      endif

      ! Project every state through the overlap in one matrix product,
      !   T = S C, then form the Mulliken projection P(mu,j) for all
      !   states at once (PSEUDOCODE 34), the same recast used on the
      !   Gaussian path.  In the gamma build the eigenvectors are real,
      !   so the projection is C(mu,j) * T(mu,j) with no conjugate.
      !   The eigenvector columns are sliced explicitly as 1:numSt (see
      !   the note on the Gaussian path in computeDOS): the SCF slab is
      !   (vDim, vDim, spin) and only its first numSt columns are the
      !   states of interest.
#ifndef GAMMA
      call projectStatesOntoBasis(valeValeOL, &
            & valeVale(:,1:numSt,1), vDim, numSt, statesProjected)
      mullikenProj(:,:) = real(conjg(valeVale(:,1:numSt,1)) * &
            & statesProjected(:,:), double)
#else
      call projectStatesOntoBasis(valeValeOLGamma, &
            & valeValeGamma(:,1:numSt,1), vDim, numSt, &
            & statesProjected)
      mullikenProj(:,:) = valeValeGamma(:,1:numSt,1) * &
            & statesProjected(:,:)
#endif

      do j = 1, numSt

         ! Step-function occupancy for the electron number diagnostic. States
         !   below the Fermi level (shifted to 0) get 1/spin; above get 0.
         if (energyEigenValues(j,i,spinIdx) &
               & <= 0.0_double) then
            occupNum = 1.0_double &
                  & / real(spinCount, double)
         else
            occupNum = 0.0_double
         endif

         ! Reset the valence dimension tracker.
         valeDimIdx = 0

         ! Loop over all atoms and their valence states to compute Mulliken
         !   projections.
         do k = 1, numSites
            do l = 1, numAtmSt(k)

               valeDimIdx = valeDimIdx + 1

               ! Read the Mulliken projection for this basis function
               !   and state, formed above by the single matrix product
               !   (PSEUDOCODE 34.5).  This value equals the direct
               !   per-column overlap sum to the reassociation floor
               !   (PSEUDOCODE 34.2, 34.7), so the three accumulations
               !   that follow are unchanged.
               oneValeRA = mullikenProj(valeDimIdx, j)

               ! Accumulate into the channel determined by pdosIndex. The 1/spin
               !   factor ensures that the projection sums correctly for
               !   spin-polarized calculations.
               projArr(pIndex(valeDimIdx), j, i) &
                     & = projArr(pIndex( &
                     & valeDimIdx), j, i) &
                     & + oneValeRA &
                     & / real(spinCount, double)

               ! Accumulate electron number for the normalization diagnostic
               !   (same formula as the Gaussian path).
               elecNum(l, k) = elecNum(l, k) &
                     & + oneValeRA &
                     & * kPointWeight(i) * occupNum

               ! Accumulate localization index (same formula as Gaussian path).
               locIdx(j) = locIdx(j) &
                     & + oneValeRA * oneValeRA &
                     & * kPointWeight(i) &
                     & / real(spinCount, double)
            enddo
         enddo
      enddo ! j (states)

      ! Progress indicator.
      if (mod(i, 10) == 0) then
         write (20, ADVANCE="NO", &
               & FMT="(a1)") "|"
      else
         write (20, ADVANCE="NO", &
               & FMT="(a1)") "."
      endif
      if (mod(i, 50) == 0) then
         write (20, *) " ", i
      endif
      call flush(20)

   enddo ! i (kpoints)

   ! Newline after progress if needed.
   if (mod(numKP, 50) /= 0) then
      write (20, *)
      call flush(20)
   endif

   ! Clean up work arrays. projArray is kept alive for the caller to use in Pass
   !   2.
   deallocate (statesProjected)
   deallocate (mullikenProj)
#ifndef GAMMA
   if (inSCF == 0) then
      deallocate (valeValeOL)
   endif
#else
   if (inSCF == 0) then
      deallocate (valeValeOLGamma)
   endif
#endif

end subroutine computeProjections_LAT


! Pass 2 of the LAT PDOS computation (PSEUDOCODE 8.3). Loop over bands and
!   tetrahedra. For each tetrahedron, look up the four corner eigenvalues via
!   the IBZ map, sort them with a tracked permutation, compute Bloechl corner
!   weights at each energy grid point, and accumulate weighted projections into
!   pdosComplete using the channel permutation table for IBZ unfolding.
!
!   The corner weight at each energy E distributes spectral density among the
!   four tetrahedron corners. The channel permutation maps the projection stored
!   at the IBZ k-point to the correct channel at the full-mesh k-point.
!
!   pdosComplete must be initialized to zero by the caller before this
!   subroutine is called.
subroutine integratePDOS_LAT(projArr, &
      & channelPermTable, pdosComp, spinIdx, &
      & numSt, cumDOSTotal, numEPts, eScale)

   use O_Kinds
   use O_MathSubs, only: bloechlCornerDOSWt
   use O_Constants, only: hartree
   use O_KPoints, only: numKPoints, numTetrahedra, &
         & tetraVol, tetrahedra, &
         & fullKPToIBZKPMap, fullKPToIBZOpMap, &
         & kPointWeight
   use O_SecularEquation, only: energyEigenValues

   implicit none

   ! Passed parameters.
   real (kind=double), dimension(:,:,:), &
         & intent(in) :: projArr
   integer, dimension(:,:), intent(in) :: &
         & channelPermTable
   real (kind=double), dimension(:,:), &
         & intent(inout) :: pdosComp
   integer, intent(in) :: spinIdx
   integer, intent(in) :: numSt
   integer, intent(in) :: cumDOSTotal
   integer, intent(in) :: numEPts
   real (kind=double), dimension(:), &
         & intent(in) :: eScale

   ! Local variables.
   integer :: n         ! Band (state) loop index.
   integer :: t         ! Tetrahedron loop index.
   integer :: c         ! Corner loop index (1-4).
   integer :: iE        ! Energy grid loop index.
   integer :: alpha     ! Channel loop index.
   integer :: orig      ! Original corner before sort.
   integer :: minIdx    ! Sort: index of minimum.
   integer :: opR       ! Point group op for corner.
   integer :: kIBZc     ! IBZ kpoint for corner.
   integer :: permAlpha ! Permuted channel index.
   real (kind=double) :: energy   ! Current grid pt.
   real (kind=double) :: tempVal  ! Sort swap temp.
   real (kind=double) :: kpWtSum  ! sum(kPointWeight)
   integer :: tempInt             ! Sort swap temp.

   ! Per-tetrahedron arrays.
   real (kind=double), dimension(4) :: eps
   real (kind=double), dimension(4) :: cornerDOSWt
   integer, dimension(4) :: kFull  ! Full-mesh corners.
   integer, dimension(4) :: kIBZ   ! IBZ corners.
   integer, dimension(4) :: opIdx  ! Op index per corner.
   integer, dimension(4) :: sigma  ! Sort permutation.

   ! Normalization: match the kPointWeight convention used by the Gaussian path.
   !   See the comment in computeTDOS_LAT for the full explanation.
   kpWtSum = sum(kPointWeight(1:numKPoints))

   ! Loop over bands (states).
   do n = 1, numSt

      ! Loop over tetrahedra.
      do t = 1, numTetrahedra

         ! Look up the four corner k-point indices from the tetrahedra array
         !   (full mesh) and map them to IBZ k-point indices and operation
         !   indices.
         do c = 1, 4
            kFull(c) = tetrahedra(c, t)
            kIBZ(c) = &
                  & fullKPToIBZKPMap(kFull(c))
            opIdx(c) = &
                  & fullKPToIBZOpMap(kFull(c))
            eps(c) = energyEigenValues( &
                  & n, kIBZ(c), spinIdx)
            sigma(c) = c
         enddo

         ! Sort the four eigenvalues in ascending order, tracking the
         !   permutation sigma so we can map sorted corners back to their
         !   original IBZ kpoint and operation.
         do c = 1, 3
            minIdx = c
            do iE = c + 1, 4
               if (eps(iE) < eps(minIdx)) then
                  minIdx = iE
               endif
            enddo
            if (minIdx /= c) then
               tempVal = eps(c)
               eps(c) = eps(minIdx)
               eps(minIdx) = tempVal
               tempInt = sigma(c)
               sigma(c) = sigma(minIdx)
               sigma(minIdx) = tempInt
            endif
         enddo

         ! Loop over the energy grid and accumulate weighted projections into
         !   pdosComplete.
         do iE = 1, numEPts
            energy = eScale(iE)

            ! Skip if outside eigenvalue range.
            if (energy < eps(1) .or. &
                  & energy >= eps(4)) cycle

            ! Compute per-corner DOS density weights (not the cumulative corner
            !   weights from bloechlCornerWeights, which are for integrated
            !   properties only).
            call bloechlCornerDOSWt( &
                  & energy, eps, cornerDOSWt)

            ! Accumulate weighted projections. Each sorted corner c maps back to
            !   its original corner sigma(c), whose IBZ kpoint and operation
            !   index determine the projection lookup. The channel permutation
            !   table handles the IBZ unfolding of the channel index.
            do c = 1, 4
               if (abs(cornerDOSWt(c)) &
                     & < 1.0d-30) cycle
               orig = sigma(c)
               opR = opIdx(orig)
               kIBZc = kIBZ(orig)

               do alpha = 1, cumDOSTotal
                  permAlpha = channelPermTable( &
                        & opR, alpha)
                  pdosComp(alpha, iE) = &
                        & pdosComp(alpha, iE) &
                        & + cornerDOSWt(c) &
                        & * tetraVol * kpWtSum &
                        & / hartree &
                        & * projArr(permAlpha, &
                        & n, kIBZc)
               enddo
            enddo
         enddo ! iE (energy grid)
      enddo ! t (tetrahedra)

      ! Progress indicator for long computations.
      if (mod(n, 50) == 0) then
         write (20, ADVANCE="NO", &
               & FMT="(a1)") "."
         call flush(20)
      endif

   enddo ! n (states)

   write (20, *) ""

end subroutine integratePDOS_LAT


! Average the atom-resolved partial DOS over the crystal point group and say
!   in the log what that did (DESIGN 1.7; PSEUDOCODE 20).
!
!   The averaging itself is generic and lives in O_MathSubs. What belongs here
!   is everything specific to this consumer: whether the permutation table
!   exists at all, and reporting the result in terms a reader of fort.20 can
!   act on.
!
!   The equality this imposes must be visible. A run whose equivalent atoms
!   agree because they were averaged and a run whose equivalent atoms agree
!   because the integration was sound produce identical output files, and only
!   the log can tell them apart. The spread reported here is also a free
!   measurement of the residual asymmetry, which is why turning the averaging
!   off is a rarely-needed diagnostic rather than the only way to see it.
subroutine symmetrizePDOS_LAT(pdosComp, channelPermTable, &
      & cumDOSTotal, numEPts)

   use O_Kinds
   use O_MathSubs, only: symmetrizeChannels
   use O_KPoints, only: numPointOps

   implicit none

   ! Passed parameters.
   real (kind=double), dimension(:,:), &
         & intent(inout) :: pdosComp
   integer, allocatable, dimension(:,:), &
         & intent(in) :: channelPermTable
   integer, intent(in) :: cumDOSTotal
   integer, intent(in) :: numEPts

   ! Local variables.
   real (kind=double) :: largestSpread ! Biggest gap made equal.
   real (kind=double) :: largestValue  ! Peak of the spectrum.
   real (kind=double) :: relativeSpread

   ! Without the symmetry maps there is nothing to average over. That happens
   !   for kpoint style code 0, an explicit list of kpoints, where Imago never
   !   builds the full mesh and so cannot construct atomPerm or anything
   !   derived from it. Say so rather than returning quietly: a silent skip
   !   looks exactly like a completed job in every output file.
   if (.not. allocated(channelPermTable)) then
      write (20,*) "PDOS symmetrization SKIPPED: no point group maps are"
      write (20,*) "available. This happens with an explicit kpoint list"
      write (20,*) "(style code 0). Symmetry-equivalent atoms may not"
      write (20,*) "agree in the partial DOS below."
      call flush (20)
      return
   endif

   call symmetrizeChannels(pdosComp, channelPermTable, &
         & numPointOps, cumDOSTotal, numEPts, &
         & largestSpread, largestValue)

   ! Report the disagreement that was averaged away, as a fraction of the
   !   spectrum it sits in. An absolute number alone says nothing: the same
   !   gap is negligible under a tall peak and damning under a small one.
   if (largestValue > 0.0_double) then
      relativeSpread = largestSpread / largestValue
   else
      relativeSpread = 0.0_double
   endif

   write (20,*) "PDOS averaged over ",numPointOps," point group operations."
   write (20,fmt="(a,e12.5,a,e12.5)") &
         & " Largest disagreement made equal: ",largestSpread, &
         & " of peak ",largestValue
   write (20,fmt="(a,e12.5)") &
         & " That is a relative spread of: ",relativeSpread
   write (20,*) "This equality is IMPOSED by averaging, not earned by the"
   write (20,*) "integration. A large value means the tetrahedron"
   write (20,*) "decomposition is far from point-group invariant on this"
   write (20,*) "lattice (DESIGN 1.2 and 1.7)."
   call flush (20)

end subroutine symmetrizePDOS_LAT


! Unfold the Gaussian per-atom PDOS over the star of each irreducible
!   k-point (PSEUDOCODE 23, TODO C148).
!
! The counterpart of symmetrizePDOS_LAT above, sharing its machinery and
!   differing in what the operation MEANS. There, the per-corner permutation
!   has already produced the right answer and this averaging only makes
!   equal things exactly equal, so a user may switch it off. Here there is
!   no other unfolding anywhere on the pathway: the accumulation credits
!   every star member with the representative's atom-by-atom breakdown, and
!   this call is what repairs it. It is therefore not optional and takes no
!   control setting.
!
! Why one average at the end is exact, rather than a permutation inside the
!   k-point loop, is derived in PSEUDOCODE 23.1. In short: the projection at
!   an irreducible k-point is invariant under that point's little group, so
!   summing over a star equals averaging over the whole point group times
!   the star multiplicity; and the group average commutes with the sum over
!   k-points because averaging is linear.
subroutine symmetrizePDOS_Gaussian(pdosComp, channelPermTable, &
      & cumDOSTotal, numEPts)

   use O_Kinds
   use O_MathSubs, only: symmetrizeChannels
   use O_KPoints, only: numPointOps

   implicit none

   ! Passed parameters.
   real (kind=double), dimension(:,:), &
         & intent(inout) :: pdosComp
   integer, allocatable, dimension(:,:), &
         & intent(in) :: channelPermTable
   integer, intent(in) :: cumDOSTotal
   integer, intent(in) :: numEPts

   ! Local variables.
   real (kind=double) :: largestSpread ! Biggest gap made equal.
   real (kind=double) :: largestValue  ! Peak of the spectrum.
   real (kind=double) :: relativeSpread

   ! The caller already tests this, but a routine that would read an
   !   unallocated array if called wrongly should say so itself rather than
   !   rely on every future caller repeating the guard.
   if (.not. allocated(channelPermTable)) then
      write (20,*) "PDOS unfolding SKIPPED: no point group maps are"
      write (20,*) "available. This happens with an explicit kpoint list"
      write (20,*) "(style code 0), where no mesh was folded and so none"
      write (20,*) "needs unfolding. The partial DOS below is unaffected."
      call flush (20)
      return
   endif

   call symmetrizeChannels(pdosComp, channelPermTable, &
         & numPointOps, cumDOSTotal, numEPts, &
         & largestSpread, largestValue)

   if (largestValue > 0.0_double) then
      relativeSpread = largestSpread / largestValue
   else
      relativeSpread = 0.0_double
   endif

   write (20,*) "PDOS unfolded over ",numPointOps," point group operations."
   write (20,fmt="(a,e12.5,a,e12.5)") &
         & " Disagreement between symmetry-equivalent channels: ", &
         & largestSpread," of peak ",largestValue
   write (20,fmt="(a,e12.5)") &
         & " That is a relative spread of: ",relativeSpread
   write (20,*) "Unlike the tetrahedron case, this equality is EARNED and"
   write (20,*) "not imposed: the number above is the error the reduction"
   write (20,*) "introduced by crediting a whole star with one k-point's"
   write (20,*) "orientation, and averaging removes it. On an unreduced"
   write (20,*) "mesh it is zero (PSEUDOCODE 23.5)."
   call flush (20)

end subroutine symmetrizePDOS_Gaussian


end module O_DOS

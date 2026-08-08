!! SPDX-License-Identifier: ECL-2.0
!! Copyright (c) 2026 Paul Rulis

module O_OptcTransitions

   ! Import necessary modules
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define module data.
   real (kind=double) :: orbitalDiff ! Energy difference between highest
         !   occupied state and lowest unoccupied state. (Normal energy onset.)
   real (kind=double) :: energyCutoff ! The energy cut off value determined
         !   from the input data and the type of calculation.  This represents
         !   the highest energy level above 0 included in the calc.
   real (kind=double) :: energyMin ! The minimum energy in the energy window
         !   to be computed for.
   real (kind=double) :: maxTransEnergy ! The maximum transition energy
         !   determined from the input data and the type of calculation.  This
         !   represents the largest transition energy included in the calc.
   real (kind=double), allocatable, dimension (:) :: energyScale ! The fine
         !   scale used for the energy axis.

   integer :: maxPairs ! The largest value in the transCounter array defined
         !   below.
   integer, allocatable, dimension (:,:) :: transCounter ! A count of the
         !   number of transitions present for each kpoint of both spin dirs.

   integer, allocatable, dimension (:,:) :: firstOccupiedState ! Array with the
         !   index number of the first occupied state for each kpoint and spin.
   integer, allocatable, dimension (:,:) :: lastOccupiedState ! Array with the
         !   index number of the last occupied state for each kpoint and spin.
   integer, allocatable, dimension (:,:) :: firstUnoccupiedState ! Array with
         !   the index number of the first unoccupied state for each
         !   kpoint and spin.
   integer, allocatable, dimension (:,:) :: lastUnoccupiedState ! Array with the
         !   index number of the last unoccupied state for each kpoint and spin.

   real (kind=double), allocatable, dimension (:) :: indirectGap ! This is the
         !   smallest energy difference between the highest occupied state and
         !   the lowest unoccupied state for all kpoints.
   real (kind=double), allocatable, dimension (:,:) :: directGap ! This is the
         !   smallest energy difference between the highest occupied state and
         !   the lowest unoccupied state for any single kpoint.
   real (kind=double), allocatable, dimension (:,:) :: indirectGapEnergies !
         !   This is the list of upper and lower energy values for the
         !   directGap of each kpoint and spin.

   real (kind=double), allocatable, dimension (:,:,:,:) :: transitionProb
   real (kind=double), allocatable, &
         & dimension (:,:,:,:,:,:) :: transitionProbPOPTC
   real (kind=double), allocatable, dimension (:,:,:)   :: energyDiff

   ! How wide the stores and the printed records are (DESIGN 13.7,
   !   PSEUDOCODE 21.3 and 21.7). These are derived once from the
   !   direction codes read out of fort.5, and every allocation and
   !   every loop bound below is taken from them rather than from the
   !   literal 3 that used to be written everywhere.
   !
   !   direction code   stored   accumulated   printed
   !        0              1          1           1   TOTAL
   !        1              3          3           4   TOTAL x y z
   !        2              3          6           7   TOTAL xx yy zz xy xz yz
   !
   ! Three counts rather than one because the widths genuinely differ at
   !   each stage. "Stored" is what the producer keeps per transition,
   !   "accumulated" is how many spectra the deposit builds, and
   !   "printed" adds the leading isotropic column, which is formed at
   !   the moment of writing from the three diagonal entries rather than
   !   accumulated in its own slot.
   !
   ! Codes 1 and 2 store the same three numbers and differ only in what
   !   is formed from them at the deposit, because the off diagonal
   !   tensor entries are products of the SAME three components. Code 0
   !   stores one number because the direction average is invariant
   !   under every rotation, so it can be collapsed at the producer and
   !   never needs the components again.
   integer :: numStoredCompOPTC   ! Components kept per transition.
   integer :: numStoredCompPOPTC  ! The same, for the decomposed store.
   integer :: numAccumCompOPTC    ! Spectra built by the deposit.
   integer :: numAccumCompPOPTC   ! The same, for the decomposed spectra.
   integer :: numPrintColOPTC     ! Columns written per spectrum.
   integer :: numPrintColPOPTC    ! The same, for a decomposed unit.

   ! The momentum matrix element itself, kept COMPLEX rather than
   !   squared on the spot (DESIGN 13.3). The momentum operator is a
   !   vector, so unfolding a matrix element from an irreducible k-point
   !   onto a member of its star mixes the three Cartesian components:
   !
   !     P_i(Rk) = sum over j of R_ij P_j(k)
   !
   !   A squared modulus cannot be rotated -- the information needed to
   !   form the rotated square, namely the cross products between
   !   different components, is destroyed by the squaring. So the square
   !   has to happen AFTER the rotation, and that means the element must
   !   survive the producer intact. These arrays exist only when the
   !   matching direction code is 1 or 2; at code 0 the collapse is
   !   still done in the producer and the real stores above carry it.
   !
   ! Indices mirror transitionProb and transProbBanded exactly, so that
   !   the accumulation reads them with the same subscripts it already
   !   uses. See the index-order argument on transProbBanded below; it
   !   applies here unchanged.
#ifndef GAMMA
   complex (kind=double), allocatable, &
         & dimension (:,:,:,:) :: transitionMoment
   complex (kind=double), allocatable, &
         & dimension (:,:,:,:,:) :: transMomentBanded
#else
   ! At the Gamma point the matrix element is real, and there is a
   !   second reason this build needs nothing more: a Gamma-only
   !   calculation has a single k-point that every symmetry operation
   !   holds fixed, so its star has one member and no rotation is ever
   !   applied. The arrays exist so the shared code below compiles and
   !   indexes identically, not because the physics differs.
   real (kind=double), allocatable, &
         & dimension (:,:,:,:) :: transitionMoment
   real (kind=double), allocatable, &
         & dimension (:,:,:,:,:) :: transMomentBanded
#endif

   ! The occupancy weight that used to be folded into the stored square,
   !   now carried alongside it (PSEUDOCODE 21.3).
   !
   ! Why it cannot stay folded in. The stored quantity is now the
   !   element rather than its square, so folding in the occupancy would
   !   mean multiplying by the SQUARE ROOT of the initial occupancy
   !   times the final vacancy. Near a Fermi surface that vacancy is a
   !   difference of two nearly equal numbers and can land a hair below
   !   zero through round-off, and the square root of a negative number
   !   is a crash rather than a small error. Kept separate it is applied
   !   at the deposit as the plain real scalar it has always been.
   real (kind=double), allocatable, dimension (:,:,:) :: pairOccupancy
         !   (pair, kPoint, spin), matching energyDiff.
   real (kind=double), allocatable, &
         & dimension (:,:,:,:) :: bandedOccupancy
         !   (kPoint, initial band, final band, spin).

   ! Tetrahedron (LAT) integration pathway. DESIGN 12, PSEUDOCODE 19.
   !
   ! The Gaussian pathway stores each transition at a position in a list
   !   that has been SORTED BY TRANSITION ENERGY, which discards the band
   !   pair that produced it. A tetrahedron needs the same band pair at
   !   all four of its corners, and sorted position p is a different pair
   !   at each corner, so that storage cannot be reused. These arrays are
   !   the same strengths held under a band-pair index instead.
   real (kind=double), allocatable, &
         & dimension (:,:,:,:,:) :: transProbBanded
         !   (component, kPoint, initial band, final band, spin).
         !
         ! INDEX ORDER IS A PERFORMANCE DECISION, not an accident, and it
         !   is chosen against the loop that reads it rather than the one
         !   that fills it. Fortran stores the leftmost index fastest.
         !   The accumulation holds a band pair fixed and walks every
         !   tetrahedron, fetching four different k-points per
         !   tetrahedron, so the k-point index sits immediately after the
         !   component. The block (:,:,i,j,h) is then contiguous -- three
         !   components by the k-point count, tens of kilobytes -- and
         !   stays cache resident across every tetrahedron for that pair.
         !   The more obvious (component, i, j, kPoint) would stride by
         !   the component count times the two band counts on each corner
         !   fetch, and miss on every one. Do not reorder this without
         !   reordering the accumulation with it.

   ! The decomposed counterpart, carrying the pair matrix over partials
   !   alongside the band pair. Indexed (initial partial, final partial,
   !   component, kPoint, initial band, final band, spin).
   real (kind=double), allocatable, &
         & dimension (:,:,:,:,:,:,:) :: transProbPOPTCBanded
         !
         ! The same ordering argument applies but reaches a weaker
         !   conclusion, because this array is far too large to hold one
         !   band pair's slice in cache. The aim here is only that the
         !   pair matrix for one component at one corner be contiguous,
         !   and that the partial loops run the leftmost index innermost.
         !   The destination of the accumulation stays scattered whatever
         !   is done, because permuting the two partial indices is the
         !   entire purpose of the operation (DESIGN 12.4).

   ! The band ranges spanned by transProbBanded. These are the UNION over
   !   k-points and spins of the per-k-point occupied and unoccupied
   !   ranges, not any single k-point's range. A tetrahedron corner may
   !   need a band pair that another corner does not enumerate, so a
   !   store built from one k-point's range would leave holes at exactly
   !   the k-points a Fermi surface passes through.
   integer :: bandedInitLo, bandedInitHi ! Occupied band range.
   integer :: bandedFinLo,  bandedFinHi  ! Unoccupied band range.

   ! Which band pairs are worth storing at all. The Gaussian producer
   !   drops pairs failing the transition-energy cutoff as it goes, which
   !   is safe when each k-point stands alone and unsafe here, since a
   !   pair may fail at one tetrahedron corner and pass at three. The
   !   mask below is therefore taken over ALL k-points at once: a pair is
   !   kept if it is in range anywhere. Sized (initial, final, spin).
   logical, allocatable, dimension (:,:,:) :: pairIsWanted

   ! POPTC specific variables. The decomposition these describe is the
   !   two by two grid of DESIGN 11: a grouping (by type or by atomic
   !   site) crossed with a resolution (the whole group in one partial,
   !   or one partial per QN_nl radial function). A "segment" is one
   !   member of the grouping, so it is a type for detail codes 1 and 2
   !   and a site for codes 3 and 4, and each segment owns a contiguous
   !   block of partials.
   integer, allocatable, dimension (:) :: pOptcIndex ! For each valence
         !   basis function, the partial it contributes to. This array is
         !   the whole of the decomposition; everything downstream simply
         !   accumulates through it.
   integer, allocatable, dimension (:) :: segmentBase ! Where each
         !   segment's block of partials begins, as a zero based offset,
         !   with one extra final entry holding sumNumPartials. Sized
         !   (numSegments + 1).
   integer, allocatable, dimension (:) :: slotsPerSegment ! How many
         !   partials each segment owns, which is the difference of
         !   consecutive segmentBase entries. Stored rather than
         !   recomputed because the star unfolding walks it directly.
   integer, allocatable, dimension (:,:) :: partialPerm ! The image of
         !   each partial under each point operation, indexed
         !   (numPointOps, sumNumPartials). Built from atomPerm and the
         !   partial layout, and used only by the atom grouped codes to
         !   carry the pair matrix across the star of an IBZ k-point.
   integer :: sumNumPartials ! Total number of POPTC partials. The stored
         !   pair matrix is this squared, so it is also the cost driver
         !   described in DESIGN 11.4.
   integer :: initVDBI ! Index for initial state valeDim basis fns
   integer :: finVDBI  ! Index for final state valeDim basis fns



!#ifndef GAMMA
!   complex (kind=double), allocatable, dimension (:,:,:,:) :: valeValeMom
!   complex (kind=double), allocatable, dimension (:,:,:)   :: coreValeOL
!#else
!   real (kind=double), allocatable, dimension (:,:,:) :: valeValeMomGamma
!   real (kind=double), allocatable, dimension (:,:)   :: coreValeOLGamma
!#endif

   real (kind=double), allocatable, dimension (:,:) :: sigmaEAccumulator

contains

! Turn the two direction codes read from fort.5 into the two sizes they
!   govern: how many Cartesian components each store keeps per
!   transition, and how many columns each spectrum prints
!   (PSEUDOCODE 21.3 and 21.7).
!
! This must run before ANY optical store is sized. The call sits at the
!   top of getEnergyStatistics because that is the first routine to
!   allocate one, and because the codes themselves are read much
!   earlier, in parseInput. A later reader moving an allocation ahead of
!   this call would size it from an undefined width, so the ordering is
!   stated here rather than left to be inferred.
subroutine setOptcStoreSize

   ! Import necessary modules.
   use O_Input, only: cartesianCodeOPTC, cartesianCodePOPTC

   implicit none

   ! Codes 1 and 2 both keep the three Cartesian components; they part
   !   company only at the deposit, where code 1 forms the three squared
   !   magnitudes and code 2 additionally forms the three off-diagonal
   !   products. Code 0 keeps the direction average alone, which is
   !   invariant and so needs no components at all.
   if (cartesianCodeOPTC == 0) then
      numStoredCompOPTC = 1
      numAccumCompOPTC  = 1
      numPrintColOPTC   = 1
   elseif (cartesianCodeOPTC == 1) then
      numStoredCompOPTC = 3
      numAccumCompOPTC  = 3
      numPrintColOPTC   = 4
   else
      numStoredCompOPTC = 3
      numAccumCompOPTC  = 6
      numPrintColOPTC   = 7
   endif

   if (cartesianCodePOPTC == 0) then
      numStoredCompPOPTC = 1
      numAccumCompPOPTC  = 1
      numPrintColPOPTC   = 1
   elseif (cartesianCodePOPTC == 1) then
      numStoredCompPOPTC = 3
      numAccumCompPOPTC  = 3
      numPrintColPOPTC   = 4
   else
      numStoredCompPOPTC = 3
      numAccumCompPOPTC  = 6
      numPrintColPOPTC   = 7
   endif

   write (20,fmt="(a,i2,a,i2,a,i2,a)") &
         & " Optical directions, totals:   store ",numStoredCompOPTC, &
         & ", accumulate ",numAccumCompOPTC,", print ", &
         & numPrintColOPTC
   write (20,fmt="(a,i2,a,i2,a,i2,a)") &
         & " Optical directions, partials: store ",numStoredCompPOPTC,&
         & ", accumulate ",numAccumCompPOPTC,", print ", &
         & numPrintColPOPTC
   call flush (20)

end subroutine setOptcStoreSize


subroutine getEnergyStatistics(doOPTC)

   ! Import necessary modules.
   use HDF5
   use O_Kinds
   use O_Potential,       only: spin
   use O_SecularEquation, only: energyEigenValues
   use O_Populate,        only: electronPopulation
   use O_KPoints,         only: numKPoints, kPointWeight
   use O_Constants,       only: dim3, smallThresh, bigThresh, hartree
   use O_Input, only: numStates, cutoffEnOPTC, maxTransEnOPTC, &
         & totalEnergyDiffPACS, energyWindowPACS, firstInitStatePACS, &
         & lastInitStatePACS, onsetEnergySlackPACS, cutoffEnSIGE, &
         & maxTransEnSIGE, cutoffEnNLOP, maxTransEnNLOP

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define passed parameters.
   integer, intent(in) :: doOPTC

   ! Define local variables.
   integer :: h,i,j,k ! Loop index variables.
   real (kind=double) :: currentEnergyDiff
   real (kind=double) :: currentGap
   integer :: firstInit
   integer :: lastInit
   integer :: firstFin
   integer :: lastFin
   integer :: orderedIndex

   ! Settle how wide the optical stores are before any of them is sized.
   !   The allocation at the end of this routine is the first one in the
   !   program, so this is the last moment the widths can be fixed.
   call setOptcStoreSize

   ! Pull variables out of imported modules.
   if (doOPTC == 1) then   ! Standard optical properties calculation.
      ! The energy onset for standard optical properties calculations is the
      !   band gap width in eV.  The energy scale for evaluating the
      !   accumulated (and broadened) transitions will begin as close to 0 eV
      !   as conveniently possible.
      energyCutoff   = cutoffEnOPTC
      maxTransEnergy = maxTransEnOPTC
   elseif (doOPTC == 2) then  ! PACS type calculation
      ! PACS calculations have the interesting feature that the energy onset is
      !   at some very high energy that is dependent on the particular target
      !   atom being excited (elemental dependency).  This is because the
      !   transitions we are considering are from a deep core state to the
      !   conduction band.  We want that energy onset to be defined by the
      !   difference in total energy between the ground state and the excited
      !   state instead of the orbital energy difference.  (WHY?)

      !   To make that happen we must first subtract away the original energy
      !   onset which is defined by the orbital energy difference between the
      !   initial state and the lowest energy state that the electron is
      !   excited into (bottom of the conduction band).
      energyCutoff   = bigThresh

      ! Compute the max transition energy as the total energy difference
      !   between ground and excited states.  Since we want the first value to
      !   be a nice round number (multiple of 5) we subtract out so-called
      !   "slack" which is the remainder of an integer division by 5.  This is
      !   the lowest energy in the output data set (no transitions at this
      !   energy though).  Then, we add the energy window we want to compute
      !   for to get the maximum transition energy.
      energyMin = totalEnergyDiffPACS - &
            & mod(totalEnergyDiffPACS,onsetEnergySlackPACS)
      maxTransEnergy = energyMin + energyWindowPACS
   elseif (doOPTC == 3) then ! Sigma(E) type calculation
      energyCutoff   = cutoffEnSIGE
      maxTransEnergy = maxTransEnSIGE
   elseif (doOPTC == 4) then ! Nonlinear optical properties calculation
      ! Just as for the linear optical properties, the energy onset is the
      !   band gap width in eV. The energy scale for evaluating the
      !   accumulated (and broadened) transitions will begin as close to 0
      !   eV as conveniently possible.
      energyCutoff   = cutoffEnNLOP
      maxTransEnergy = maxTransEnNLOP
   endif


   ! Allocate arrays
   allocate (transCounter         (numKPoints,spin))
   allocate (firstOccupiedState   (numKPoints,spin))
   allocate (lastOccupiedState    (numKPoints,spin))
   allocate (firstUnoccupiedState (numKPoints,spin))
   allocate (lastUnoccupiedState  (numKPoints,spin))
   allocate (directGap            (numKPoints,spin))
   allocate (indirectGapEnergies  (spin,2)) ! Last index: Upper=1, Lower=2
   allocate (indirectGap          (spin))

   ! Initialize the arrays.
   transCounter(:,:)         = 0
   firstOccupiedState(:,:)   = 0
   lastOccupiedState(:,:)    = 0
   firstUnoccupiedState(:,:) = 0
   lastUnoccupiedState(:,:)  = 0
   directGap(:,:)            = bigThresh
   indirectGapEnergies(:,1)  = bigThresh
   indirectGapEnergies(:,2)  = -bigThresh
   inDirectGap(:)            = 0.0_double

   ! Initialize the running maximum for the transition pair count. This is
   !   a module variable that is only ever built up with max() below, so it
   !   needs a starting value smaller than any count it could be compared
   !   against, and zero is the natural floor for a count. Setting it here
   !   rather than leaving it to whatever the module variable happens to
   !   hold keeps the result independent of how a given compiler chooses
   !   to lay out and pre-fill module storage.
   maxPairs = 0

   ! The purpose of this subroutine is to gather important statistics and
   !   indices for use later on.  The important values that will be determined
   !   are:  1) The first and last occupied state and the first and last 
   !   unoccupied state for each kpoint and each spin that each different type
   !   of calculation cares about (e.g. PACS transitions are from one core
   !   state to numerous CB states; OPTC transitions are from many VB states to
   !   numerous CB states; and SIGE considers a finite range near the fermi
   !   level); 2) The direct band gap for each spin; 3) The minimum indirect
   !   band gap for each spin; 4) The number of transitions for each kpoint and
   !   spin so that we can allocate memory easily; 5) whatever ...
   ! We must take some extra consideration of the situation in which thermal
   !   smearing is present. It should be understood that in this case the
   !   population statistics are a bit different than the 0K case where the
   !   Fermi edge is a flat step function. For thermal smearing we have the
   !   condition that states near the Fermi level that are fully occupied at
   !   0K will be partially occupied at a finite temperature. Further, states
   !   that are totally unoccupied at 0K will also be partially occupied at a
   !   finite temperature. Thus we must consider the case of transitions
   !   between occupied and unoccupied states a bit carefully.
   ! Consider first the PACS case. Here we have the initial state(s) being core
   !   states that will always be fully occupied in the ground state regardless
   !   of temperature (becaues they are so deep). The electron from the core
   !   state may transition into any state with partial or zero occupation.
   !   Thus, when making a list of the states to transition into, we must find
   !   the first one that has a non-negligable portion that is unoccupied (i.e.
   !   a mostly occupied VB state that is a bit far-ish from the traditional
   !   top of the VB.) The core electron may transition into this state, but
   !   the probability of doing so must be affected by (of course) the
   !   momentum matrix element that represents conservation of angular momentum
   !   selection rules (i.e. s->p and p->s,d etc) *and* the fact that the
   !   state at that energy is mostly occupied already (e.g. 99%). Therefore,
   !   the transition calculation proceeds as normal except that the intensity
   !   must be scaled by a factor of 0.01. Without thermal smearing, the
   !   probability of transitioning into this state would be 0%. As the higher
   !   energy states are considered, the so-called occupation scaling factor
   !   will increase to a maximum of 1.0 in the case that the core electron is
   !   transitioning into a state that has zero initial occupation.
   ! Consider second the case of a traditional VB optical properties
   !   calculation. Now the situation is even more complicated because the
   !   initial state may be partially occupied as well as the final state.
   !   However, similar physical principles will apply and the probability of
   !   making a transition will be scaled by two occupation scaling factors,
   !   one for the originating state and one for the final state.


   do h = 1, spin

      ! Pacs calculations need to have the resultant spectra shifted according
      !   to the difference in orbital energies from the ground and excited
      !   states.  To find the amount of the shift we initialize the search
      !   number.
      if (doOPTC == 2) then  ! Doing PACS calculation.
         orbitalDiff = bigThresh
      endif

      do i = 1, numKPoints
         do j = 2, numStates

            ! Find the last occupied state for this KPoint and spin. Also note
            !   that this will deal with degenerate highest occupied states.
            ! Note: In the case that thermal smearing is turned on, then the
            !   zero of energy is the Fermi level (which may be in the gap for
            !   insulators). There may be partially occupied orbitals above
            !   the Fermi level and partially occupid orbitals below the Fermi
            !   level because of the smearing.

            ! Transitions are from the core states into any of the partially
            !   occupied states scaled by the degree of population for PACS
            !   calculations.
            ! Transitions are from the valence states (including any partially
            !   occupied states above or below the fermi level) into any of the
            !   partially occupied states above or below the fermi level (that
            !   are also *above* the initial state) for traditional optical
            !   properties calculations.
            ! Thermal excitations in conducting materials are taken to occur
            !   between any occupied (or partially occupied) state and any
            !   unoccupied (or partially occupied state of higher energy) that
            !   is within a designated range (usually on the order of a few eV
            !   or tenths of an eV). This is for sigma(E) calculations.

            ! Determine the array index value of the current spin-kpoint-state
            !   as defined by the tempEnergyEigenValues loop near the beginning
            !   of the population subroutine.
            orderedIndex = j+numStates*(h-1)+numStates*spin*(i-1)

            ! Note that we should anticipate that the following "if" statement
            !   will never be true when j=1. This means that the electron
            !   population of the first (lowest) state will always be the
            !   maximum possible (the kPointWeight). This condition might be
            !   broken if the thermal smearing was set so insanely high that
            !   it smeared all the way to the lowest state (some -20 eV). This
            !   would correspond to rediculous temperatures. Just to be sure
            !   though, the j-loop starts at 2. (The first state will not be
            !   the first (partially) unoccupied state.) This is important
            !   because some j-1 and j-k with k/=0 type calculations are done
            !   inside this if-statement.
            ! This "if" statement will be true any time that we find a state
            !   with less than full occupation.
            if (abs(electronPopulation(orderedIndex)-kPointWeight(i) / &
                  & real(spin,double)) > smallThresh) then

               ! In the condition that we are doing a PACS calculation then
               !   the initial state core state(s) *will* be partially
               !   occupied because we have pulled an electron out of them.
               !   These states will never be final states and so we have to
               !   skip them.
               if (doOPTC == 2) then
                  if (j <= lastInitStatePACS) cycle
               endif

               ! We have found a state that is either partially or fully
               !   unoccupied. Therefore, this state may be the lowest energy
               !   state that an electron *could* transfer into (indirect gap)
               !   considering that each spin and kpoint has its own set of
               !   states. If this state is lower in energy than any other
               !   previously found state with non-occupation then we record it
               !   as the final state for the indirect gap for this spin. (Note
               !   that for metals with thermal smearing turned on, this state
               !   will be below the fermi level while there will be partially
               !   occupied states above the fermi level. Thus, there will be
               !   no gap. Metals should have no gap.)
               if (indirectGapEnergies(h,1) > energyEigenValues(j,i,h)) then
                  indirectGapEnergies(h,1) = energyEigenValues(j,i,h)
               endif

               ! Store the index values for the first unoccupied state for this
               !   kpoint and spin.
               if (firstUnoccupiedState(i,h) == 0) then
                  firstUnoccupiedState(i,h) = j
               endif

               ! In the event that the previous state is partially or fully
               !   occupied then we want to check the energy difference
               !   between it and the current state (which is partially
               !   occupied or fully unoccupied) to see if it is a "gap". For
               !   the 0K case with insulators, this is an obvious calculation,
               !   but for the finite temperature case with thermal smearing
               !   the situation is more complex.  Essentailly, for the current
               !   kpoint, we are looking for the largest separation between a
               !   state with some electrons in it and a state that isn't fully
               !   occupied. At the moment we are at a state that has less than
               !   full occupation and we are checking the previous state to
               !   see if it has any electrons.
               if (electronPopulation(orderedIndex-1) > 0.0_double) then

                  ! Obviously, this is also an occupied state and the last time
                  !   that we get inside this "if" statement (inside the upper
                  !   one too) we will have found the last occupied state.
                  !   Thus, we assume that every time is the last time and
                  !   eventually it will be correct.  (This is only useful for
                  !   traditional VB optical properties calculations and
                  !   sigma(E) calculations. It is not useful for PACS
                  !   calculatoins because they will use the user specified
                  !   core states by overriding this determination later on.)
                  lastOccupiedState(i,h) = j-1

                  ! Similarly, we will find the highest occupied state energy.
                  if (indirectGapEnergies(h,2)<energyEigenValues(j-1,i,h)) then
                     indirectGapEnergies(h,2) = energyEigenValues(j-1,i,h)
                  endif

                  ! Compute the size of the energy difference between the
                  !   current partially occupied or fully unoccupied state and
                  !   the previous (j-1) state which is either partially
                  !   occupied or fully occupied. Note that we had to check
                  !   that the previous state was at least partially occupied.
                  currentGap = abs(energyEigenValues(j,i,h) - &
                        & energyEigenValues(j-1,i,h))

                  ! This so-called currentGap will become the direct gap *for
                  !   this kpoint and spin* if it is minimal compared to other
                  !   previous currentGap values for this kpoint (and spin).
                  !   Note that the comparison between kpoints will be done
                  !   later. Also, just remember that we are always only
                  !   comparing the energy difference between adjacent states
                  !   for the current kpoint.
                  if (currentGap < directGap(i,h)) then
                     ! The largest currentGap for this kpoint might be a
                     !   directGap when compared with other directGap values
                     !   from other kpoints (to be determined later).
                     directGap(i,h) = currentGap
                  endif
               else ! The previous state has no electrons in it.
                  ! No need to search through any higher states in the outer j
                  !   loop.  All relevant information has been obtained.
                  exit
               endif ! The previous state has at least some electrons in it.

            endif ! (Found a partially occupied or fully unoccupied state)
         enddo ! (j = numStates)

         ! If we are doing a sige calculation then we are only concerned
         !   with states that are a few eV from the Fermi energy.  This
         !   means that in addition to having to seek out the first
         !   unoccupied state and the last occupied state we also need to
         !   identify the last unoccupied state (which will be *just*
         !   above the last occupied state according to the
         !   maxTransEnergy given in the input file) and the first
         !   occupied state (which will be *just* below the first
         !   unoccupied state according to the maxTransEnergy given in
         !   the input file (imago.dat)).
         ! Note that this search only needs to be done once per kpoint so this
         !   is why we are doing it after the numStates loop.
         if (doOPTC == 3) then

            ! Loop higher than the last occupied state to find the
            !   lowest unoccupied state that is *greater* than the maximum
            !   transition energy from the Fermi level.
            do k = lastOccupiedState(i,h), numStates
               if (energyEigenValues(k,i,h) > &
                     & (energyEigenValues(lastOccupiedState(i,h),i,h) + &
                     & maxTransEnergy)) then
                  lastUnoccupiedState(i,h) = k-1
                  exit
               endif
               if (k == numStates) then
                  lastUnoccupiedState(i,h) = numStates
               endif
            enddo

            ! Loop lower than the first unoccupied state to find the
            !   highest occupied state that is *less* than the maximum
            !   transition energy from the Fermi level.
            do k = 1, firstUnoccupiedState(i,h)-1
               if (energyEigenValues(firstUnoccupiedState(i,h)-k,i,h) < &
                        & (energyEigenValues(firstUnoccupiedState(i,h),i,h) - &
                        & maxTransEnergy)) then
                  firstOccupiedState(i,h) = firstUnoccupiedState(i,h)-k+1
                  exit
               endif
               if (k == firstUnoccupiedState(i,h)-1) then
                  firstOccupiedState(i,h) = 1
               endif
            enddo
         endif ! doOPTC 3

         ! Initialize the counter for the number of transitions for this kpoint.
         transCounter(i,h) = 0

         ! For normal optical properties calculations the last unoccupied state
         !   we care about is the last (highest) state in the calculation, and
         !   the first occupied state is always the first (lowest) state in the
         !   calculation.  For PACS, these values depend on which core state
         !   has been excited.  For Sigma(E) these values depend on the range
         !   around the Fermi level to consider (defined by the user) and the
         !   determination above.

         if (doOPTC == 1) then ! Normal optical properties calculation.
            firstOccupiedState(i,h) = 1
            ! lastOccupiedState determined above.
            ! firstUnoccupiedState determined above.
            lastUnoccupiedState(i,h) = numStates
         elseif (doOPTC == 2) then ! Doing a PACS calculation.
            firstOccupiedState(i,h) = firstInitStatePACS ! From O_Input
            lastOccupiedState(i,h)  = lastInitStatePACS  ! From O_Input
            ! firstUnoccupiedState determined above.
            lastUnoccupiedState(i,h) = numStates
         elseif (doOPTC == 3) then ! Doing a Sigma(E) calculation.
            ! firstOccupiedState determined above.
            ! lastOccupiedState determined above.
            ! firstUnoccupiedState determined above.
            ! lastUnoccupiedState determined above.
         elseif (doOPTC == 4) then ! Non-linear optical properties.
            ! The non-linear properties have an input block of their own
            !   (NLOP_INPUT_DATA, read by readNlopControl) and a job ID
            !   that reaches this far, but no routine anywhere computes
            !   them: there is no counterpart to computePairs or
            !   computeSigmaE for the second order response. Stopping here
            !   is deliberate. The alternative is to fall through and
            !   silently emit whichever spectrum the surrounding code
            !   happens to produce, labelled as a non-linear result.
            stop "Non-linear optical properties are not implemented."
         else
            ! Error, no other options.
            stop "Check optical properties command line parameter: doOPTC"
         endif

         ! Store the state variables for temporary use as loop indices.
         firstInit = firstOccupiedState(i,h)
         lastInit  = lastOccupiedState(i,h)
         firstFin  = firstUnoccupiedState(i,h)
         lastFin   = lastUnoccupiedState(i,h)

         if (doOPTC == 2) then ! Doing PACS calculation.
            ! Determine the orbital energy difference for this kpoint and
            !   compare it to the smallest difference yet obtained.
            orbitalDiff = min(orbitalDiff,&
                  & energyEigenValues(firstFin,i,h) - &
                  & energyEigenValues(lastInit,i,h))
         endif

         ! Loop over all the possible transitions for this kpoint to determine
         !   the number of accepted transitions for this kpoint.  We will also
         !   refine the value for the lastUnoccupied state to consider.
         do j = firstInit, lastInit
            do k = firstFin, lastFin

               ! An important note is that with thermal smearing turned on,
               !   there is the potential for a state to be *both* an initial
               !   and a finel state because it is partially occupied. If we
               !   encounter any such states we will not count the case where
               !   the final state is lower in energy that the initial state.
               if (j >= k) cycle

               ! If the energy of the final state is higher than the requested
               !   cut-off then we adjust the record for the last unoccupied
               !   state and go to the next initial state because all the
               !   remaining final states will be greater. This should never be
               !   entered for the sigma(E) (doOPTC==3) case because the
               !   energyCutoff should always be set much larger than the
               !   maxTransEnergy. Thus, the first/last Occ./Unocc. states
               !   index numbers will all be clustered around the Fermi energy
               !   index number where the quality of the state functions is
               !   expected to be high (and so we should have no reason to
               !   discard them).
               if (energyEigenValues(k,i,h) > energyCutoff) then
                  lastUnoccupiedState(i,h) = k-1
                  exit
               endif

               ! Compute the energy of transition between the current states.
               currentEnergyDiff = energyEigenValues(k,i,h) - &
                     & energyEigenValues(j,i,h)

               ! Check if the energy difference is less than the maximum
               !   transition energy that the input file requested computation
               !   for. If it fails, then we go to the next initial state
               !   because all the remaining final states for this energy will
               !   be greater.
               if (currentEnergyDiff > maxTransEnergy) then
                  exit
               endif

               ! At this point the transition is valid and one we would want to
               !   compute.  Unfortunately there isn't a good way to save the
               !   energyDiff computation that we did above for later use, and
               !   it will have to be done again.  (BOO)  (Unless we do some
               !   static memory allocation. (BOO))
               transCounter(i,h) = transCounter(i,h) + 1
            enddo ! (k First to Last Fin)
         enddo ! (j First to Last Init)
      enddo ! (i kpoints)

      ! Determine the indirect band gap.
      indirectGap(h) = indirectGapEnergies(h,1) - indirectGapEnergies(h,2)
      if (indirectGap(h) < 0.0_double) then
         indirectGap(h) = 0.0_double
      endif

      ! Write the indirect band gap and determine+write the direct band gap.
      if (spin == 1) then
         write (20,*) "Indirect Band Gap(eV) = ",indirectGap(h)*hartree
         write (20,*) "Direct Band Gap(eV)   = ",minval(directGap(:,h))*hartree
         call flush (20)
      elseif (h == 1) then
         write (20,*) "(Up) Indirect Band Gap(eV) = ",indirectGap(h)*hartree
         write (20,*) "(Up) Direct Band Gap(eV)   = ",minval(directGap(:,h))* &
               & hartree
         call flush (20)
      else ! spin == 2 and h == 2
         write (20,*) "(Dn) Indirect Band Gap(eV) = ",indirectGap(h)*hartree
         write (20,*) "(Dn) Direct Band Gap(eV)   = ",minval(directGap(:,h))* &
               & hartree
         call flush (20)
      endif

      ! Determine the maximum number of transition pairs of all the kpoints and
      !   for both spins.
      maxPairs = max(maxPairs,maxval(transCounter(:,h)))
   enddo ! (h spin)

   ! Now that the number of transitions for each kpoint are known we can 
   !   allocate space to hold information based on the number transitions.
   !   Note that we allocate transitionProb here regardless of whether we
   !   are doing POPTC or not, because we will *always* do a total optc.
   if (doOPTC /= 3) then  ! Not doing a Sigma(E)
      allocate (energyDiff (maxPairs,numKPoints,spin))
      allocate (transitionProb (numStoredCompOPTC,maxPairs,numKPoints,&
            & spin))

      ! Initialize these arrays.
      energyDiff(:maxPairs,:numKPoints,:) = 0.0_double
      transitionProb(:,:,:,:) = 0.0_double

      ! At direction codes 1 and 2 the producer keeps the complex matrix
      !   element instead of its square, so that the star unfolding can
      !   rotate the components before they are squared, and the
      !   occupancy weight travels beside it rather than inside it
      !   (PSEUDOCODE 21.3). At code 0 neither array is needed: the
      !   direction average is rotation invariant, so the producer
      !   collapses it on the spot into transitionProb above.
      if (numStoredCompOPTC > 1) then
         allocate (transitionMoment (numStoredCompOPTC,maxPairs, &
               & numKPoints,spin))
         allocate (pairOccupancy (maxPairs,numKPoints,spin))
         transitionMoment(:,:,:,:) = 0.0_double
         pairOccupancy(:,:,:) = 0.0_double
      endif
   endif

end subroutine getEnergyStatistics



subroutine computeTransitions(inSCF,doOPTC)

   ! Import the necessary modules.
   use HDF5
   use O_Kinds
   use O_TimeStamps
   use O_Potential,     only: spin
   use O_Constants,     only: dim3
   use O_KPoints,       only: numKPoints, kPointIntgCode
   use O_AtomicSites,   only: coreDim, valeDim
   use O_CommandLine,   only: serialXYZ
   use O_Input,         only: numStates, totalEnergyDiffPACS, detailCodePOPTC
#ifndef GAMMA
   use O_SecularEquation, only: valeVale, valeValeMM, readDataSCF, &
         & readDataPSCF
#else
   use O_SecularEquation, only: valeValeGamma, valeValeMMGamma, readDataSCF, &
         & readDataPSCF
#endif

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define passed parameters.
   integer, intent(in) :: inSCF
   integer, intent(in) :: doOPTC

   ! Define local variables.
   integer :: h,i,j ! Loop index variables
integer :: k,l
   real (kind=double) :: energyShift
!   real (kind=double), allocatable, dimension (:,:) :: tempRealValeVale
!#ifndef GAMMA
!   real (kind=double), allocatable, dimension (:,:) :: tempImagValeVale
!#endif


   ! Log the date and time we start.
   call timeStampStart (23)

   ! Allocate the matrix to hold the wave function and momentum matrix
   !   elements and initialize them.
#ifndef GAMMA
   if (inSCF == 0) then
      allocate (valeVale (valeDim,numStates,1))
      valeVale(:,:,1) = cmplx(0.0_double,0.0_double,double)
   endif
   allocate (valeValeMM (valeDim,valeDim,3))
   valeValeMM(:,:,:) = cmplx(0.0_double,0.0_double,double)
#else
   if (inSCF == 0) then
      allocate (valeValeGamma(valeDim,numStates,1))
      valeValeGamma(:,:,1) = 0.0_double
   endif
   allocate (valeValeMMGamma (valeDim,valeDim,3))
   valeValeMMGamma(:,:,:) = 0.0_double
#endif

   ! Build the decomposition index once, before any producer runs. It is
   !   fixed by the structure and the detail code, so it does not vary
   !   over the loop below, and the tetrahedron store cannot even be
   !   sized without the partial count it yields (PSEUDOCODE 18).
   if (detailCodePOPTC /= 0) then
      call buildPOPTCIndex
   endif

   ! For the tetrahedron pathway, decide the band-pair store's shape and
   !   allocate it before the loop below begins. It cannot be a per-call
   !   array the way the Gaussian temporaries are: the accumulation reads
   !   four k-points at a time, so the whole store has to be resident at
   !   once (DESIGN 12.4, PSEUDOCODE 19.2).
   if (kPointIntgCode == 1) then
      call selectBandedPairs

      allocate (transProbBanded (numStoredCompOPTC,numKPoints, &
            & bandedInitLo:bandedInitHi,bandedFinLo:bandedFinHi,spin))

      ! Zeroed because the pair mask leaves gaps: a pair no tetrahedron
      !   wants is never written, and the accumulation must read a zero
      !   there rather than whatever the allocation happened to contain.
      transProbBanded(:,:,:,:,:) = 0.0_double

      ! At direction codes 1 and 2 the element and its occupancy weight
      !   are held separately so the star unfolding can rotate before it
      !   squares (PSEUDOCODE 21.3). The occupancy carries no component
      !   index: it depends on the two bands and the k-point alone, so
      !   one number serves all three components.
      if (numStoredCompOPTC > 1) then
         allocate (transMomentBanded (numStoredCompOPTC,numKPoints, &
               & bandedInitLo:bandedInitHi,bandedFinLo:bandedFinHi, &
               & spin))
         allocate (bandedOccupancy (numKPoints, &
               & bandedInitLo:bandedInitHi,bandedFinLo:bandedFinHi, &
               & spin))
         transMomentBanded(:,:,:,:,:) = 0.0_double
         bandedOccupancy(:,:,:,:) = 0.0_double
      endif

      ! The decomposed store, when a decomposition was requested. This
      !   is the array whose size DESIGN 11.4 warns about: it carries
      !   the partial count SQUARED on top of everything the total store
      !   holds, so the atom-nl cell reaches tens of gigabytes on a cell
      !   of a few tens of atoms.
      if (detailCodePOPTC /= 0) then
         allocate (transProbPOPTCBanded (sumNumPartials,sumNumPartials, &
               & dim3,numKPoints,bandedInitLo:bandedInitHi, &
               & bandedFinLo:bandedFinHi,spin))
         transProbPOPTCBanded(:,:,:,:,:,:,:) = 0.0_double
      endif
   endif


   do h = 1, spin

      ! Record the fact that we are starting the k-point loop so that when
      !   someone looks at the output file as the job is running they will
      !   know that all the little dots represent a count of the number of
      !   k-points. Then they can figure out the progress and progress rate.
      write (20,*) "Beginning k-point loop."
      if (numKPoints > 1) write (20,*) "Expecting ",numKPoints," iterations."
      call flush (20)

      ! Begin a loop over the number of kpoints
      do i = 1, numKPoints

         ! Determine if we are doing the OPTC in a post-SCF calculation, or
         !   within an SCF calculation. Generally, we will not need to read in
         !   the energy eigenvalues here because we already did it.
         if (inSCF == 1) then
            ! Read necessary data from SCF (setup,main) data structures.
            if (doOPTC /= 2) then ! Not doing a PACS calculation
               call readDataSCF(h,i,numStates,2) ! 2 = regular MME matrixCode
            else
               call readDataSCF(h,i,numStates,3) ! 3 = PACS MME matrixCode
            endif
         else
            ! Read necessary data from post SCF data structures.
            if (doOPTC /= 2) then ! Not doing a PACS calculation
               call readDataPSCF(h,i,numStates,2) ! 2 = regular MME matrixCode
            else
               call readDataPSCF(h,i,numStates,3) ! 3 = PACS MME matrixCode
            endif
         endif


!         ! Allocate temporary reading matrices.
!#ifndef GAMMA
!         allocate (tempRealValeVale(valeDim,numStates))
!         allocate (tempImagValeVale(valeDim,numStates))
!#else
!         allocate (tempRealValeVale(valeDim,numStates))
!#endif
!
!         if (doOPTC /= 2) then  ! Not doing a PACS calculation.
!
!            ! Read the datasets for this kpoint.
!#ifndef GAMMA
!            call readMatrix(eigenVectorsBand_did(:,i,h),valeVale(:,:,1,1),&
!                  & tempRealValeVale(:,:),tempImagValeVale(:,:),&
!                  & valeStatesBand,valeDim,numStates)
!#else
!            call readMatrixGamma(eigenVectorsBand_did(1,i,h),&
!                  & valeValeGamma(:,:,1),valeStatesBand,valeDim,numStates)
!#endif
!         else
!
!#ifndef GAMMA
!            ! Read the data for the ground state for this kpoint.
!            call readPartialWaveFns(eigenVectorsBand_did(:,i,h),&
!                  & valeVale(:,:,1,1),tempRealValeVale(:,:),&
!                  & tempImagValeVale(:,:),valeStatesBand,&
!                  & firstOccupiedState(i,h),lastOccupiedState(i,h),&
!                  & valeDim,numStates)
!
!            ! Read the data for the excited state for this kpoint.
!            call readPartialWaveFns(eigenVectorsBand2_did(:,i,h),&
!                  & valeVale(:,:,1,1),tempRealValeVale(:,:),&
!                  & tempImagValeVale(:,:),valeStatesBand,&
!                  & lastOccupiedState(i,h)+1,lastUnoccupiedState(i,h),&
!                  & valeDim,numStates)
!#else
!            ! Read the data for the ground state for this kpoint.
!            call readPartialWaveFnsGamma(eigenVectorsBand_did(1,i,h),&
!                  & valeValeGamma(:,:,1),tempRealValeVale(:,:),valeStatesBand,&
!                  & firstOccupiedState(i,h),lastOccupiedState(i,h),&
!                  & valeDim,numStates)
!
!            ! Read the data for the excited state for this kpoint.
!            call readPartialWaveFnsGamma(eigenVectorsBand2_did(1,i,h),&
!                  & valeValeGamma(:,:,1),tempRealValeVale(:,:),valeStatesBand,&
!                  & lastOccupiedState(i,h)+1,lastUnoccupiedState(i,h),&
!                  & valeDim,numStates)
!#endif
!         endif
!
!
!         ! Read the orthogonalization coefficients after allocating space to
!         !   hold them (and the temp reading matrix).
!#ifndef GAMMA
!         deallocate (tempRealValeVale)
!         deallocate (tempImagValeVale)
!         allocate   (tempRealValeVale (coreDim,valeDim))
!         allocate   (tempImagValeVale (coreDim,valeDim))
!         allocate   (coreValeOL (coreDim,valeDim,1))
!         if (coreDim /= 0) then
!            call readMatrix(coreValeBand_did(:,i),coreValeOL(:,:,1),&
!                  & tempRealValeVale(:,:),tempImagValeVale(:,:),&
!                  & coreValeBand,coreDim,valeDim)
!         endif
!         deallocate (tempRealValeVale)
!         deallocate (tempImagValeVale)
!#else
!         deallocate (tempRealValeVale)
!         allocate   (coreValeOLGamma  (coreDim,valeDim))
!         if (coreDim /= 0) then
!            call readMatrixGamma(coreValeBand_did(1,i),coreValeOLGamma(:,:),&
!                  & coreValeBand,coreDim,valeDim)
!         endif
!#endif


         ! Perform the computations in serial or all together.
         if (serialXYZ == 0) then
!#ifndef GAMMA
!do j = 1, 3
!do k = 1, valeDim
!do l = 1, valeDim
!write(20,*) k,l,valeValeMM(l,k,j)
!enddo
!enddo
!enddo
!#endif
!#ifndef GAMMA
!
!            allocate   (valeValeMom (valeDim,valeDim,1,3))
!            ! Get the integral results for the x, y, z momentum matrices.
!            ! Runcode:  3 = XMom, 4 = YMom, 5 = ZMom; (j+2)
!            do j = 1, 3
!               call getIntgResults (valeValeMom(:,:,:,j),coreValeOL,&
!                     & i,j+2,valeValeBand_did(i),valeValeBand,1,1)
!            enddo
!            ! Deallocate matrices that are no longer needed in this iteration to
!            !   make space for those that are needed.
!            deallocate (coreValeOL)
!#else
!            allocate   (valeValeMomGamma (valeDim,valeDim,3))
!            ! Get the integral results for the x, y, z momentum matrices.
!            ! Runcode:  3 = XMom, 4 = YMom, 5 = ZMom; (j+2)
!            do j = 1, 3
!               call getIntgResults (valeValeMomGamma(:,:,j),coreValeOLGamma,&
!                     & j+2,valeValeBand_did(i),valeValeBand,1,1)
!            enddo
!            ! Deallocate matrices that are no longer needed in this iteration to
!            !   make space for those that are needed.
!            deallocate (coreValeOLGamma)
!#endif

            if (doOPTC /= 3) then  ! Not doing a Sigma(E) calculation.

               ! The two Brillouin-zone integration methods need the
               !   transition strengths held under different indices, so
               !   the producer is chosen here rather than inside one
               !   (DESIGN 12.2). Both are called at the same point in
               !   the loop, immediately after the read that loaded this
               !   k-point's wave functions and momentum matrix.
               if (kPointIntgCode == 1) then

                  ! The total spectra are needed whether or not a
                  !   decomposition was asked for, exactly as on the
                  !   Gaussian side, so the undecomposed store is filled
                  !   in both cases and the decomposed one in addition.
                  call computeTransProbBanded (i,0,h)
                  if (detailCodePOPTC /= 0) then
                     call computeTransProbPOPTCBanded (i,0,h)
                  endif
               elseif (detailCodePOPTC == 0) then ! Standard OPTC calc.
                  call computePairs (i,0,h,doOPTC)
               else ! Doing pOptc calculation
                  call computePOPTCPairs (i,0,h,doOPTC)
               endif
            else
               call computeSigmaE (i,0,h)
            endif

         else
!            do j = 1, 3
!#ifndef GAMMA
!               allocate   (valeValeMom (valeDim,valeDim,1,1))
!               ! Runcode:  3 = XMom, 4 = YMom, 5 = ZMom; (j+2)
!               call getIntgResults (valeValeMom(:,:,:,1),coreValeOL,&
!                     & i,j+2,valeValeBand_did(i),valeValeBand,1,1)
!#else
!               allocate   (valeValeMomGamma (valeDim,valeDim,1))
!               ! Runcode:  3 = XMom, 4 = YMom, 5 = ZMom; (j+2)
!               call getIntgResults (valeValeMomGamma(:,:,1),coreValeOLGamma,&
!                     & j+2,valeValeBand_did(i),valeValeBand,1,1)
!#endif
!
!               if (doOPTC /= 3) then  ! Not doing a Sigma(E) calc.
!                  if (detailCodePOPTC == 0) then ! Doing standard OPTC calc.
!                     call computePairs (i,j,h,doOPTC)
!                  else ! Doing pOptc calculation
!                     call computePOPTCPairs (i,j,h,doOPTC)
!                  endif
!               else
!                  call computeSigmaE (i,j,h)
!               endif
!            enddo
!
!#ifndef GAMMA
!            deallocate (coreValeOL)
!#else
!            deallocate (coreValeOLGamma)
!#endif
         endif


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

      enddo ! i numKPoints

      ! Add a final return if one was not already made.
      if (mod(numKPoints,10) /= 0) then
         write (20,*)
      endif


      ! In the case of XANES calculations we need to shift the calculated
      !   energy values.  The amount of the shift is equal to the difference
      !   between (the total energy difference between the ground and excited
      !   states) and (the calculated LUMO,core difference).
      if (doOPTC == 2) then  ! Doing PACS calculation.
         energyShift = totalEnergyDiffPACS - orbitalDiff
         energyDiff(:,:,h) = energyDiff(:,:,h) + energyShift
      endif
   enddo ! (h spin)

   ! Deallocate unnecessary matrices and arrays
#ifndef GAMMA
!   deallocate (valeVale)
   deallocate (valeValeMM)
#else
!   deallocate (valeValeGamma)
   deallocate (valeValeMMGamma)
#endif
   deallocate (firstOccupiedState)
   deallocate (lastOccupiedState)
   deallocate (firstUnoccupiedState)
   deallocate (lastUnoccupiedState)
   deallocate (directGap)
   deallocate (indirectGapEnergies)
   deallocate (indirectGap)

   ! NOTE that neither the decomposition index nor pairIsWanted is
   !   released here, even though this routine allocated both. The
   !   accumulation runs after this routine returns and reads them: the
   !   mask to skip band pairs that were never filled, and partialPerm
   !   to permute the two partial indices at each tetrahedron corner. A
   !   lifetime spans every consumer, so both belong to subroutine optc,
   !   which owns the whole optical phase (PSEUDOCODE 19.2.1).
   !
   ! Releasing the index here instead was wrong, and wrong silently.
   !   The accumulation guards its use of partialPerm on
   !   allocated(partialPerm), because the type grouped detail codes
   !   never build the table at all. That guard cannot tell a table
   !   that was legitimately never built from one freed too early, so
   !   the permutation was simply skipped and a plausible spectrum came
   !   out with the wrong per-atom values.

   ! Log the date and time we end.
   call timeStampEnd (23)

end subroutine computeTransitions


! Decide the shape of the band-pair store before any of it is filled.
!
! Two jobs, both of which need only eigenvalues and so are essentially
!   free beside the matrix element work that follows. The first is to
!   find the band ranges the store must span. The second is to decide
!   which band pairs inside those ranges are worth storing at all.
!
! Why the pruning has to be redone rather than inherited. The Gaussian
!   producer drops a pair the moment its transition energy exceeds the
!   requested maximum, and does so k-point by k-point. That is correct
!   when each k-point is integrated on its own. It is wrong here: the
!   tetrahedron that owns a k-point also owns three others, and a pair
!   that is out of range at one corner may be in range at the other
!   three. Dropping it at the first corner would leave the tetrahedron
!   unable to form its energy difference. So the test below asks whether
!   a pair is in range at ANY k-point, and keeps it everywhere if so.
!
! Why prune at all, given that the expensive work is unpruned anyway.
!   The dominant cost of the producer is building conjWaveMomSum, which
!   runs over the whole final-state range before any cutoff applies and
!   is therefore identical in both pathways. What pruning saves is not
!   time but STORAGE, and that is worth having: every retained pair
!   costs three doubles per k-point here, and three times the square of
!   the partial count in the decomposed case.
subroutine selectBandedPairs

   ! Import the necessary modules.
   use O_Kinds
   use O_Potential,       only: spin
   use O_KPoints,         only: numKPoints
   use O_SecularEquation, only: energyEigenValues

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define local variables.
   integer :: h,j,k ! Loop indices: spin, initial band, final band.
   real (kind=double) :: smallestEnergyDiff

   ! The union of the per-k-point ranges over every k-point and spin.
   !   getEnergyStatistics has already filled these arrays, one entry per
   !   k-point and spin, so this is a scan rather than a computation.
   bandedInitLo = minval(firstOccupiedState  (:numKPoints,:spin))
   bandedInitHi = maxval(lastOccupiedState   (:numKPoints,:spin))
   bandedFinLo  = minval(firstUnoccupiedState(:numKPoints,:spin))
   bandedFinHi  = maxval(lastUnoccupiedState (:numKPoints,:spin))

   allocate (pairIsWanted(bandedInitLo:bandedInitHi, &
         & bandedFinLo:bandedFinHi, spin))
   pairIsWanted(:,:,:) = .false.

   do h = 1, spin
      do j = bandedInitLo, bandedInitHi
         do k = bandedFinLo, bandedFinHi

            ! A final state below an initial state is not a transition.
            !   The Gaussian producer applies the same rule.
            if (j >= k) cycle

            ! The smallest transition energy this pair reaches anywhere
            !   in the zone. If even that exceeds the requested maximum
            !   then the pair is out of range at every k-point and no
            !   tetrahedron can want it.
            smallestEnergyDiff = minval( &
                  & energyEigenValues(k,:numKPoints,h) &
                  & - energyEigenValues(j,:numKPoints,h))

            if (smallestEnergyDiff <= maxTransEnergy) then
               pairIsWanted(j,k,h) = .true.
            endif
         enddo
      enddo
   enddo

end subroutine selectBandedPairs


! Build the sum, over all basis functions, of the conjugated wave function
!   of each final state against the momentum matrix. This is the expensive
!   part of forming a transition probability, and it is factored out here
!   because it does not depend on the initial state: every initial state
!   that pairs with a given final state reuses the same sum, so it is
!   computed once per final state rather than once per pair.
!
! It lives on its own for a second reason that matters more as the code
!   grows. The Brillouin-zone integration methods differ in how they
!   INDEX and filter the transitions they store, not in the physics of
!   the matrix element itself (DESIGN 12.2). Keeping that physics in one
!   routine means the several producers that will consume it cannot drift
!   apart in the one place where a discrepancy would be hardest to see.
!
! NOTE for the Gamma-point build: this routine MUTATES valeValeMMGamma,
!   negating its upper triangle to restore Hermiticity. That is a side
!   effect on a module array, and the invariant it needs is ONE CALL PER
!   READ of that array -- not one call per k-point, which for a Gamma
!   build would be vacuous since there is only the zone-centre point.
!
! It holds today because computeTransitions re-reads the momentum matrix
!   through readDataSCF or readDataPSCF at the top of every spin and
!   k-point iteration, and calls this routine exactly once afterwards. A
!   spin-polarized Gamma run therefore calls it twice against two
!   separate reads, which is correct. What would break it is a caller
!   that reads once and then calls more than once -- for instance a
!   restored serial-XYZ path looping the three Cartesian components
!   around a single read. The triangle would be negated back and every
!   transition probability after it would change, with no symptom.
subroutine buildConjWaveMomSum (firstFin,lastFin,initComponent, &
      & finComponent,conjWaveMomSum)

   ! Import the necessary modules.
   use O_Kinds
   use O_AtomicSites, only: valeDim
#ifndef GAMMA
   use O_SecularEquation, only: valeVale, valeValeMM
#else
   use O_SecularEquation, only: valeValeGamma, valeValeMMGamma
#endif

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   integer, intent(in) :: firstFin
   integer, intent(in) :: lastFin
   integer, intent(in) :: initComponent
   integer, intent(in) :: finComponent
#ifndef GAMMA
   complex (kind=double), dimension (:,:,:), intent(out) :: conjWaveMomSum
#else
   real    (kind=double), dimension (:,:,:), intent(out) :: conjWaveMomSum
#endif

   ! Define local variables.
   integer :: i,j,k ! Loop index variables.
   integer :: finalStateIndex

#ifndef GAMMA

   ! Compute the sum over the final states. The 1 for the valeVale is for
   !   the 1 kpoint. The finComponent is 3 for all three components at
   !   once, and 1 for when X, Y, Z are done separately.
   do i = initComponent, finComponent
      finalStateIndex = 0
      do j = firstFin, lastFin
         ! Define the final index for conjWaveMomSum
         finalStateIndex = finalStateIndex + 1
         do k = 1, valeDim
            conjWaveMomSum(k,finalStateIndex,i) = &
                  & sum(conjg(valeVale(:,j,1)) * valeValeMM(:,k,i))
         enddo
      enddo
   enddo

#else

   ! Documentation similar to the above non-gamma case.
   do i = initComponent, finComponent

      ! Make the upper triangle correct for Hermiticity.  Recall that for
      !   the Gamma K Point all the matrices are real (except the momentum
      !   matrix which was multiplied by a -i and is hence imaginary).
      !   Since it must be Hermitian we need to apply that now. See the
      !   note above this routine: this write is why it may be called only
      !   once per k-point per component set.
      do j = 1, valeDim
         valeValeMMGamma(1:j,j,i) = -valeValeMMGamma(1:j,j,i)
      enddo

      finalStateIndex = 0
      do j = firstFin, lastFin

         ! Increment the finalStateIndex for conjWaveMomSum
         finalStateIndex = finalStateIndex + 1

         do k = 1, valeDim
            conjWaveMomSum(k,finalStateIndex,i) = &
                  & sum(valeValeGamma(:,j,1) * valeValeMMGamma(:,k,i))
         enddo
      enddo
   enddo

#endif

end subroutine buildConjWaveMomSum



subroutine computePairs (currentKPoint,xyzComponents,spinDirection,doOPTC)

   ! Import the necessary modules.
   use O_Kinds
   use O_Constants,   only: dim3
   use O_Potential,   only: spin
   use O_AtomicSites, only: valeDim
   use O_SortSubs,    only: mergeSort
   use O_Input,       only: numStates
   use O_KPoints,     only: kPointWeight
   use O_Populate,    only: electronPopulation
#ifndef GAMMA
   ! The momentum matrix itself is read by buildConjWaveMomSum rather than
   !   here, so only the wave functions are needed at this level.
   use O_SecularEquation, only: energyEigenValues, valeVale
#else
   use O_SecularEquation, only: energyEigenValues, valeValeGamma
#endif

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   integer, intent(in) :: currentKPoint
   integer, intent(in) :: xyzComponents ! 0=all, 1=x, 2=y, 3=z
   integer, intent(in) :: spinDirection
   integer, intent(in) :: doOPTC

   ! Define local variables.
   integer :: i,j,k ! Loop index variables
   integer :: initComponent
   integer :: finComponent
   integer :: transPairCount
   integer :: firstInit
   integer :: lastInit
   integer :: firstFin
   integer :: lastFin
   integer :: finalStateIndex
   integer :: orderedIndex
   integer, allocatable, dimension (:) :: sortOrder
   integer, allocatable, dimension (:) :: segmentBorders
   real    (kind=double) :: initStateFactor
   real    (kind=double) :: finStatefactor
   real    (kind=double) :: currentEnergyDiff
   real    (kind=double), allocatable, dimension (:)     :: energyDiffTemp
   real    (kind=double), allocatable, dimension (:,:)   :: transitionProbTemp
   ! The unsorted companions of the two direction-resolved stores. Every
   !   quantity this routine produces has to be held under the order the
   !   pairs were MET and then copied into place under the order they
   !   sort into, so each store needs its own temporary here.
   real    (kind=double), allocatable, dimension (:) :: occupancyTemp
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:) :: conjWaveMomSum
   complex (kind=double),              dimension (dim3)  :: valeValeXMom
   complex (kind=double), allocatable, dimension (:,:) :: transMomentTemp
#else
   real    (kind=double), allocatable, dimension (:,:,:) :: conjWaveMomSumGamma
   real    (kind=double),              dimension (dim3)  :: valeValeXMomGamma
   real    (kind=double), allocatable, dimension (:,:) :: transMomentTemp
#endif

   ! Initialize a counter for the current number of transition pairs
   transPairCount = 0

   ! Make shorthand for the state indices.
   firstInit = firstOccupiedState(currentKPoint,spinDirection)
   lastInit  = lastOccupiedState(currentKPoint,spinDirection)
   firstFin  = firstUnoccupiedState(currentKPoint,spinDirection)
   lastFin   = lastUnoccupiedState(currentKPoint,spinDirection)

   ! Determine the range of components (xyz) that should be considered.
   if (xyzComponents == 0) then
      initComponent = 1
      finComponent = 3
   else
      initComponent = 1
      finComponent = 1
   endif


   ! Allocate space to hold the sum(conjg(valeVale(:,j)) * valeVale_Mom(:,k,1))
   !   for each of the possible final states.  This is done since the values
   !   are independent of the initial states.  The finComponent is 3 for all
   !   three at once, and 1 for when X, Y, Z are done separately.
   !
   ! The sum itself is built by buildConjWaveMomSum, which both Brillouin
   !   -zone integration methods share (DESIGN 12.2).
#ifndef GAMMA
   allocate (conjWaveMomSum (valeDim,lastFin-firstFin+1,finComponent))

   call buildConjWaveMomSum (firstFin,lastFin,initComponent,finComponent, &
         & conjWaveMomSum)
#else
   allocate (conjWaveMomSumGamma (valeDim,lastFin-firstFin+1,finComponent))

   call buildConjWaveMomSum (firstFin,lastFin,initComponent,finComponent, &
         & conjWaveMomSumGamma)
#endif


   ! Allocate space for the energy difference.
   allocate (energyDiffTemp (maxPairs))

   ! Only one of the two temporaries below is ever filled, chosen by the
   !   direction code: at code 0 the square is formed here and there is
   !   nothing to rotate later, while at codes 1 and 2 the element is
   !   carried out intact for the star unfolding to rotate. Both are
   !   allocated regardless so that the copy-into-place block at the end
   !   has a valid array to reference on either branch, and because the
   !   unused one costs a single k-point's worth of pairs.
   allocate (transitionProbTemp (finComponent,maxPairs))
   allocate (transMomentTemp (finComponent,maxPairs))
   allocate (occupancyTemp (maxPairs))

   ! Initialize the temporary energy transition array.
   energyDiffTemp(:) = 0.0_double
   occupancyTemp(:) = 0.0_double

   ! Allocate space to hold the indices for each segment of the energyDiff
   !   array.
   allocate (segmentBorders (lastInit-firstInit+2))

   ! Initialize the first index since it will always be 0.
   segmentBorders(1) = 0

   ! Begin the double loop to determine the transition energies.
   do i = firstInit, lastInit
      do j = firstFin, lastFin

         ! Index into conjWaveMomSum, which was filled above for *every*
         !   final state in the firstFin to lastFin range. The index is
         !   derived from j rather than accumulated by a counter because
         !   the loop below skips some j: with thermal smearing a state
         !   can be both initial and final, and those skipped j still
         !   occupy their slot in conjWaveMomSum. A counter incremented
         !   only on accepted pairs would fall behind at the first skip
         !   and read a different final state's momentum sum from then
         !   on, silently, for the rest of this initial state.
         finalStateIndex = j - firstFin + 1

         ! Recall that thermal smearing may allow some states to be both
         !   initial and final. We do not consider transitions where the final
         !   state has an energy less than the initial.
         if (i >= j) cycle

         ! If the energy of the final state is higher than the requested
         !   cut-off we go to the next initial state.
         if (energyEigenValues(j,currentKPoint,spinDirection) > &
               & energyCutoff) exit

         ! Compute the energy of the transition from the current states.
         currentEnergyDiff = energyEigenValues(j,currentKPoint,spinDirection)-&
               & energyEigenValues(i,currentKPoint,spinDirection)

         ! Check if the energy difference is less than the maximum
         !   transition energy that the input file requested computation
         !   for.  If it fails, then we go to the next initial state because
         !   all the remaining final states for this energy will be greater.
         if (currentEnergyDiff > maxTransEnergy) exit

         ! Increment the number of transition pairs counted so far.
         transPairCount = transPairCount + 1

         ! Store the transition energy for the current pair.
         energyDiffTemp(transPairCount) = currentEnergyDiff

         ! In the event that thermal smearing is turned on. The state that the
         !   e- comes from and goes into may be fully, partially, or not
         !   occupied. We will scale the probability of a transition linearly
         !   according to the percent occupation of both the initial and final
         !   states.

         ! Determine the array index value of the current initial (index i)
         !   spin-kpoint-state as defined by the tempEnergyEigenValues loop
         !   near the beginning of the population subroutine.
         orderedIndex = i + numStates*(spinDirection-1) + &
               & numStates*spin*(currentKPoint-1)

         ! Use the normal state factor for non-PACS calculations. For PACS
         !   calculations the initStateFactor is always 1 even though the
         !   initial core state(s) will have an electron missing.
         if (doOPTC /= 2) then
            initStateFactor = electronPopulation(orderedIndex) / &
                  & (kPointWeight(currentKPoint)/real(spin,double))
         else
            initStateFactor = 1.0_double
         endif

         ! Determine the array index value of the current final (index j)
         !   spin-kpoint-state as defined by the tempEnergyEigenValues loop
         !   near the beginning of the population subroutine.
         orderedIndex = j + numStates*(spinDirection-1) + &
               & numStates*spin*(currentKPoint-1)

         finStateFactor = 1.0_double - electronPopulation(orderedIndex) / &
               & (kPointWeight(currentKPoint)/real(spin,double))

         ! The occupancy weight for this pair. At direction code 0 it is
         !   folded into the square immediately below, exactly as it
         !   always has been; at codes 1 and 2 it is carried out of here
         !   in its own array so that it never has to be square rooted
         !   (PSEUDOCODE 21.3).
         occupancyTemp(transPairCount) = initStateFactor * finStateFactor

#ifndef GAMMA

         ! Loop to obtain the wave function times the momentum integral.
         do k = initComponent,finComponent
             valeValeXMom(k) = sum(valeVale(:,i,1) * &
                   & conjWaveMomSum(:,finalStateIndex,k))
         enddo

         if (numStoredCompOPTC == 1) then

            ! Direction code 0. Collapse to the direction average here
            !   and keep nothing else: the sum of the three squared
            !   magnitudes is unchanged by any rotation, so no later
            !   stage ever needs the components back. Note that the
            !   squared magnitude is what the physics asks for even
            !   though the element looks real, because the momentum
            !   matrix carries a factor of -i applied in getIntgResults
            !   -- what is stored as a real number is the y in (x+iy).
            transitionProbTemp(1,transPairCount) = &
                  & sum(real(valeValeXMom(initComponent:finComponent), &
                  & double)**2 &
                  & + aimag(valeValeXMom(initComponent:finComponent)) &
                  & **2) * initStateFactor * finStateFactor
         else

            ! Direction codes 1 and 2. Carry the element out intact. It
            !   is squared only after the star unfolding has rotated it,
            !   because squaring first destroys the cross products
            !   between components that a rotated square is built from.
            transMomentTemp(initComponent:finComponent,transPairCount) &
                  & = valeValeXMom(initComponent:finComponent)
         endif
#else

         ! Loop to get the wave function times the momentum matrix element.
         do k = initComponent,finComponent
            valeValeXMomGamma(k) = sum(valeValeGamma(:,i,1) * &
                  & conjWaveMomSumGamma(:,finalStateIndex,k))
         enddo

         if (numStoredCompOPTC == 1) then

            ! Direction code 0, as above but with a real element.
            transitionProbTemp(1,transPairCount) = &
                  & sum(valeValeXMomGamma(initComponent:finComponent) &
                  & **2) * initStateFactor * finStateFactor
         else
            transMomentTemp(initComponent:finComponent,transPairCount) &
                  & = valeValeXMomGamma(initComponent:finComponent)
         endif
#endif
      enddo ! Fin loop j

      ! Save the index for the end border of this segment.
      segmentBorders(i - firstInit + 2) = transPairCount
   enddo ! Init loop i

   ! Deallocate unnecessary matrix
#ifndef GAMMA
   deallocate (conjWaveMomSum)
#else
   deallocate (conjWaveMomSumGamma)
#endif

   ! Determine if there was only one segment.  In this case we don't have to
   !   sort anything.

   ! Sort energyDiffTemp into energyDiff, and obtain the indices for the
   !   correct sorted order of energyDiff so that we can copy the energy
   !   momentum directly.

   allocate (sortOrder (transPairCount))

   call mergeSort (energyDiffTemp,energyDiff(:,currentKPoint,spinDirection),&
         & sortOrder,segmentBorders,transPairCount)

   ! Copy the temporaries into the real stores using the sorting order
   !   determined in the mergeSort subroutine. Which store is filled
   !   follows the direction code, and it is the same choice the deposit
   !   loop above made.
   if (numStoredCompOPTC == 1) then

      ! Direction code 0. One collapsed number per transition.
      if (xyzComponents == 0) then
         do i = 1, transPairCount
            transitionProb(1,i,currentKPoint,spinDirection) = &
                  & transitionProbTemp(1,sortOrder(i))
         enddo
      else

         ! One component per call, so the direction average has to be
         !   built up across the three calls rather than written. The
         !   sort order is a function of the transition energies alone,
         !   which do not depend on the component, so position i means
         !   the same pair on every call and the sum is well defined.
         do i = 1, transPairCount
            transitionProb(1,i,currentKPoint,spinDirection) = &
                  & transitionProb(1,i,currentKPoint,spinDirection) &
                  & + transitionProbTemp(1,sortOrder(i))
         enddo
      endif
   else

      ! Direction codes 1 and 2. The complex element and its occupancy
      !   weight travel separately, and the weight is written rather
      !   than accumulated because it is the same on all three calls.
      if (xyzComponents == 0) then
         do i = 1, transPairCount
            transitionMoment(:,i,currentKPoint,spinDirection) = &
                  & transMomentTemp(:,sortOrder(i))
            pairOccupancy(i,currentKPoint,spinDirection) = &
                  & occupancyTemp(sortOrder(i))
         enddo
      else
         do i = 1, transPairCount
            transitionMoment(xyzComponents,i,currentKPoint, &
                  & spinDirection) = transMomentTemp(1,sortOrder(i))
            pairOccupancy(i,currentKPoint,spinDirection) = &
                  & occupancyTemp(sortOrder(i))
         enddo
      endif
   endif

   ! Deallocate unnecessary arrays and matrices
   deallocate (energyDiffTemp)
   deallocate (transitionProbTemp)
   deallocate (transMomentTemp)
   deallocate (occupancyTemp)
   deallocate (segmentBorders)
   deallocate (sortOrder)

end subroutine computePairs


! Fill one k-point's slice of the band-pair transition store, for the
!   tetrahedron integration pathway. PSEUDOCODE 19.2.
!
! This is the counterpart of computePairs, and the physics inside the
!   two is the same: both build conjWaveMomSum through the shared
!   routine and both form the squared momentum matrix element weighted
!   by the initial state's occupancy and the final state's vacancy. What
!   differs is entirely bookkeeping, and in five ways. This routine
!   indexes by band pair rather than by position in a list; it does not
!   sort; it spans the union of the band ranges rather than this
!   k-point's own; it prunes from a mask taken over all k-points at once
!   rather than dropping pairs as it meets them; and it takes its
!   occupancies from the tetrahedron scheme rather than the Gaussian
!   one.
!
! Why it is called per k-point rather than owning a k-point loop. The
!   wave functions and the momentum matrix are module arrays holding ONE
!   k-point: readDataSCF and readDataPSCF overwrite them on every pass,
!   and valeValeMM has no k-point dimension at all. A routine that
!   looped over k-points itself would see only whichever k-point was
!   read last, so it would have to perform the reads too -- duplicating
!   the logic that chooses between the SCF and post-SCF sources and
!   between the regular and PACS matrix codes. Sitting inside the
!   existing loop avoids that second copy, and satisfies for free the
!   rule that buildConjWaveMomSum is called exactly once per read.
subroutine computeTransProbBanded (currentKPoint,xyzComponents, &
      & spinDirection)

   ! Import the necessary modules.
   use O_Kinds
   use O_AtomicSites, only: valeDim
   use O_KPoints,     only: kPointWeight
   use O_Populate,    only: electronPopulation_LAT
#ifndef GAMMA
   use O_SecularEquation, only: valeVale
#else
   use O_SecularEquation, only: valeValeGamma
#endif

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   integer, intent(in) :: currentKPoint
   integer, intent(in) :: xyzComponents ! 0=all, 1=x, 2=y, 3=z
   integer, intent(in) :: spinDirection

   ! Define local variables.
   integer :: i,j,k ! Loop indices: initial band, final band, component.
   integer :: initComponent
   integer :: finComponent
   integer :: storeComponent  ! Where in dim3 this component belongs.
   integer :: finalStateIndex ! Position of band j within conjWaveMomSum.
   real (kind=double) :: initStateFactor
   real (kind=double) :: finStateFactor
   real (kind=double) :: fullOccupancy
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:) :: conjWaveMomSum
   complex (kind=double) :: valeValeXMom
#else
   real (kind=double), allocatable, dimension (:,:,:) :: conjWaveMomSum
   real (kind=double) :: valeValeXMom
#endif

   ! Determine the range of components (xyz) that should be considered.
   if (xyzComponents == 0) then
      initComponent = 1
      finComponent = 3
   else
      initComponent = 1
      finComponent = 1
   endif

   ! The occupancy denominator. electronPopulation_LAT holds a
   !   Brillouin-zone VOLUME FRACTION rather than an occupancy: summed
   !   over every k-point it reaches one for a fully occupied band. The
   !   share belonging to a single k-point is that point's own volume
   !   fraction, which is half its weight because the weights sum to two
   !   by Imago's convention. Dividing by it leaves the pure zero-to-one
   !   occupancy the transition probability needs.
   !
   ! The zone measure must be removed here rather than carried along,
   !   because on this pathway it re-enters through tetraVol during the
   !   accumulation. Leaving it in would apply the measure twice. The
   !   Gaussian producer performs the matching division by
   !   kPointWeight/spin for the same reason, against its own array's
   !   convention (DESIGN 1.6d, and the conversion in computeBond).
   fullOccupancy = kPointWeight(currentKPoint) * 0.5_double

   ! Build the shared momentum sum over the FULL final-state range. It is
   !   deliberately not restricted to the wanted pairs: the cost is set
   !   by the range rather than by the pairs, and every wanted pair needs
   !   a final state somewhere in it.
   allocate (conjWaveMomSum (valeDim,bandedFinHi-bandedFinLo+1, &
         & finComponent))

#ifndef GAMMA
   call buildConjWaveMomSum (bandedFinLo,bandedFinHi,initComponent, &
         & finComponent,conjWaveMomSum)
#else
   call buildConjWaveMomSum (bandedFinLo,bandedFinHi,initComponent, &
         & finComponent,conjWaveMomSum)
#endif

   do i = bandedInitLo, bandedInitHi

      ! The occupancy of the initial state at this k-point.
      initStateFactor = electronPopulation_LAT(i,currentKPoint, &
            & spinDirection) / fullOccupancy

      do j = bandedFinLo, bandedFinHi

         ! Skip pairs that no tetrahedron anywhere in the zone can want.
         !   The mask already excludes the case of a final state below an
         !   initial one.
         if (.not. pairIsWanted(i,j,spinDirection)) cycle

         ! Where band j sits within conjWaveMomSum. Derived from the band
         !   index rather than counted, so that skipping a pair cannot
         !   put this out of step with the array it addresses.
         finalStateIndex = j - bandedFinLo + 1

         ! The vacancy of the final state at this k-point.
         finStateFactor = 1.0_double - electronPopulation_LAT(j, &
               & currentKPoint,spinDirection) / fullOccupancy

         ! The occupancy weight for this band pair at this k-point. At
         !   direction code 0 it is folded into the square below, as it
         !   always has been; at codes 1 and 2 it is stored beside the
         !   element so that it never has to be square rooted
         !   (PSEUDOCODE 21.3). It carries no component index because
         !   occupancy does not depend on direction.
         if (numStoredCompOPTC > 1) then
            bandedOccupancy(currentKPoint,i,j,spinDirection) = &
                  & initStateFactor * finStateFactor
         endif

         do k = initComponent, finComponent

            ! When the components are done one at a time the computed
            !   component always lands in slot 1 of conjWaveMomSum but
            !   belongs in slot xyzComponents of the store. At direction
            !   code 0 there is only one slot and every component adds
            !   into it, so the position is 1 whichever way the
            !   components were computed.
            if (numStoredCompOPTC == 1) then
               storeComponent = 1
            elseif (xyzComponents == 0) then
               storeComponent = k
            else
               storeComponent = xyzComponents
            endif

#ifndef GAMMA
            valeValeXMom = sum(valeVale(:,i,1) &
                  & * conjWaveMomSum(:,finalStateIndex,k))

            if (numStoredCompOPTC == 1) then

               ! Direction code 0. Accumulate the three squared
               !   magnitudes into the single slot to form the direction
               !   average, which no rotation can change. The square is
               !   the right quantity even though the element looks
               !   real, because the momentum matrix carries a factor of
               !   -i applied in getIntgResults: what is stored as a
               !   real number is the y in (x+iy).
               transProbBanded(1,currentKPoint,i,j,spinDirection) = &
                     & transProbBanded(1,currentKPoint,i,j, &
                     & spinDirection) + (real(valeValeXMom,double)**2 &
                     & + aimag(valeValeXMom)**2) &
                     & * initStateFactor * finStateFactor
            else

               ! Direction codes 1 and 2. Keep the element itself; it is
               !   squared only after the star unfolding has rotated it.
               transMomentBanded(storeComponent,currentKPoint,i,j, &
                     & spinDirection) = valeValeXMom
            endif
#else
            valeValeXMom = sum(valeValeGamma(:,i,1) &
                  & * conjWaveMomSum(:,finalStateIndex,k))

            if (numStoredCompOPTC == 1) then
               transProbBanded(1,currentKPoint,i,j,spinDirection) = &
                     & transProbBanded(1,currentKPoint,i,j, &
                     & spinDirection) + valeValeXMom**2 &
                     & * initStateFactor * finStateFactor
            else
               transMomentBanded(storeComponent,currentKPoint,i,j, &
                     & spinDirection) = valeValeXMom
            endif
#endif
         enddo
      enddo
   enddo

   deallocate (conjWaveMomSum)

end subroutine computeTransProbBanded


! Fill one k-point's slice of the DECOMPOSED band-pair store, for the
!   tetrahedron pathway. The counterpart of computeTransProbBanded, and
!   the same relationship to computePOPTCPairs that that routine has to
!   computePairs: identical physics, different bookkeeping.
!
! What the extra work here is. Instead of collapsing the momentum matrix
!   element to a single number per transition, the element is resolved
!   into a matrix over partial pairs, and the transition probability is
!   distributed over that matrix so that summing it reproduces the
!   total. The construction is the "sum squared to sum of squares"
!   arrangement: each entry is multiplied by the SUM over all entries,
!   separately for the real and imaginary parts, which is what makes the
!   partials add up to the undecomposed answer.
subroutine computeTransProbPOPTCBanded (currentKPoint,xyzComponents, &
      & spinDirection)

   ! Import the necessary modules.
   use O_Kinds
   use O_Constants,   only: dim3
   use O_AtomicSites, only: valeDim
   use O_KPoints,     only: kPointWeight
   use O_Populate,    only: electronPopulation_LAT
#ifndef GAMMA
   use O_SecularEquation, only: valeVale
#else
   use O_SecularEquation, only: valeValeGamma
#endif

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   integer, intent(in) :: currentKPoint
   integer, intent(in) :: xyzComponents ! 0=all, 1=x, 2=y, 3=z
   integer, intent(in) :: spinDirection

   ! Define local variables.
   integer :: i,j,k ! Loop indices: initial band, final band, component.
   integer :: l,n,o ! Loop indices: basis function, final and initial
         !   partial. The pair matrix is addressed (o,n) so that the
         !   initial partial is leftmost and therefore fastest.
   integer :: initComponent
   integer :: finComponent
   integer :: storeComponent  ! Where in dim3 this component belongs.
   integer :: finalStateIndex ! Position of band j within conjWaveMomSum.
   real (kind=double) :: initStateFactor
   real (kind=double) :: finStateFactor
   real (kind=double) :: fullOccupancy
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:,:) :: conjWaveMomSum
   complex (kind=double), allocatable, dimension (:,:,:) :: valeValeXMom
   real (kind=double) :: valeValeXMomSumReal
   real (kind=double) :: valeValeXMomSumImag
#else
   real (kind=double), allocatable, dimension (:,:,:,:) :: conjWaveMomSum
   real (kind=double), allocatable, dimension (:,:,:) :: valeValeXMom
   real (kind=double) :: valeValeXMomSum
#endif

   ! Determine the range of components (xyz) that should be considered.
   if (xyzComponents == 0) then
      initComponent = 1
      finComponent = 3
   else
      initComponent = 1
      finComponent = 1
   endif

   ! See computeTransProbBanded for why the zone measure is divided out
   !   here: it re-enters through tetraVol during the accumulation, and
   !   applying it twice would be silent.
   fullOccupancy = kPointWeight(currentKPoint) * 0.5_double

   allocate (conjWaveMomSum (valeDim,sumNumPartials, &
         & bandedFinHi-bandedFinLo+1,finComponent))
   allocate (valeValeXMom (sumNumPartials,sumNumPartials,dim3))

   call buildConjWaveMomSumPOPTC (bandedFinLo,bandedFinHi,initComponent, &
         & finComponent,conjWaveMomSum)

   do i = bandedInitLo, bandedInitHi

      initStateFactor = electronPopulation_LAT(i,currentKPoint, &
            & spinDirection) / fullOccupancy

      do j = bandedFinLo, bandedFinHi

         ! Skip pairs that no tetrahedron anywhere in the zone can want.
         if (.not. pairIsWanted(i,j,spinDirection)) cycle

         ! Derived rather than counted, so that a skipped pair cannot put
         !   this out of step with the array it addresses.
         finalStateIndex = j - bandedFinLo + 1

         finStateFactor = 1.0_double - electronPopulation_LAT(j, &
               & currentKPoint,spinDirection) / fullOccupancy

         do k = initComponent, finComponent

            if (xyzComponents == 0) then
               storeComponent = k
            else
               storeComponent = xyzComponents
            endif

#ifndef GAMMA
            valeValeXMom(:,:,k) = cmplx(0.0_double,0.0_double,double)

            ! Resolve the matrix element by the partial its initial-state
            !   basis function belongs to. The final-state partial is
            !   already carried by conjWaveMomSum's second index.
            do l = 1, valeDim
               do n = 1, sumNumPartials
                  valeValeXMom(pOptcIndex(l),n,k) = &
                        & valeValeXMom(pOptcIndex(l),n,k) &
                        & + valeVale(l,i,1) &
                        & * conjWaveMomSum(l,n,finalStateIndex,k)
               enddo
            enddo

            ! The totals that turn a sum of squares back into a squared
            !   sum, so that the partials reproduce the undecomposed
            !   transition probability when added together.
            valeValeXMomSumReal = sum(real(valeValeXMom(:,:,k),double))
            valeValeXMomSumImag = sum(aimag(valeValeXMom(:,:,k)))

            ! n outer, o inner: o is the leftmost index of both the
            !   scratch matrix and the store, so it must run innermost.
            do n = 1, sumNumPartials
               do o = 1, sumNumPartials
                  transProbPOPTCBanded(o,n,storeComponent,currentKPoint, &
                        & i,j,spinDirection) = &
                        & ((real(valeValeXMom(o,n,k),double) &
                        & * valeValeXMomSumReal) &
                        & + (aimag(valeValeXMom(o,n,k)) &
                        & * valeValeXMomSumImag)) &
                        & * initStateFactor * finStateFactor
               enddo
            enddo
#else
            valeValeXMom(:,:,k) = 0.0_double

            do l = 1, valeDim
               do n = 1, sumNumPartials
                  valeValeXMom(pOptcIndex(l),n,k) = &
                        & valeValeXMom(pOptcIndex(l),n,k) &
                        & + valeValeGamma(l,i,1) &
                        & * conjWaveMomSum(l,n,finalStateIndex,k)
               enddo
            enddo

            valeValeXMomSum = sum(valeValeXMom(:,:,k))

            do n = 1, sumNumPartials
               do o = 1, sumNumPartials
                  transProbPOPTCBanded(o,n,storeComponent,currentKPoint, &
                        & i,j,spinDirection) = &
                        & valeValeXMom(o,n,k) * valeValeXMomSum &
                        & * initStateFactor * finStateFactor
               enddo
            enddo
#endif
         enddo
      enddo
   enddo

   deallocate (conjWaveMomSum)
   deallocate (valeValeXMom)

end subroutine computeTransProbPOPTCBanded


! The decomposed counterpart of buildConjWaveMomSum. It answers the same
!   question -- for each final state, the conjugated wave function summed
!   against the momentum matrix -- but keeps the answer resolved by
!   partial rather than collapsed to a single number, so that the pair
!   matrix of DESIGN 11 can be formed from it afterwards.
!
! The extra index is supplied by pOptcIndex, which section 18 of the
!   pseudocode fills: it sends each basis function to the partial it
!   belongs to. Every basis function therefore ACCUMULATES into its
!   partial's slot rather than assigning to its own, which is why this
!   array is zeroed first and the undecomposed one is not.
!
! The same Gamma-point caution applies as for buildConjWaveMomSum: the
!   Gamma build negates the upper triangle of valeValeMMGamma in place,
!   so this routine may be called only once per READ of that matrix. See
!   the note above that routine for why the read, rather than the
!   k-point, is the unit that matters.
subroutine buildConjWaveMomSumPOPTC (firstFin,lastFin,initComponent, &
      & finComponent,conjWaveMomSum)

   ! Import the necessary modules.
   use O_Kinds
   use O_AtomicSites, only: valeDim
#ifndef GAMMA
   use O_SecularEquation, only: valeVale, valeValeMM
#else
   use O_SecularEquation, only: valeValeGamma, valeValeMMGamma
#endif

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   integer, intent(in) :: firstFin
   integer, intent(in) :: lastFin
   integer, intent(in) :: initComponent
   integer, intent(in) :: finComponent
#ifndef GAMMA
   complex (kind=double), dimension (:,:,:,:), intent(out) :: conjWaveMomSum
#else
   real    (kind=double), dimension (:,:,:,:), intent(out) :: conjWaveMomSum
#endif

   ! Define local variables.
   integer :: i,j,k ! Loop index variables.
   integer :: basisFn ! The basis function currently being deposited.
   integer :: finalStateIndex

   ! Every basis function adds into the slot its partial owns, so the
   !   destination must start empty.
#ifndef GAMMA
   conjWaveMomSum = cmplx(0.0_double,0.0_double,double)
#else
   conjWaveMomSum = 0.0_double
#endif

#ifndef GAMMA

   ! Compute the sum over the final states.
   do i = initComponent, finComponent
      finalStateIndex = 0

      do j = firstFin, lastFin
         ! Define the final index for conjWaveMomSum
         finalStateIndex = finalStateIndex + 1

         do k = 1, valeDim

            ! Walk every basis function and send its contribution to the
            !   partial that pOptcIndex assigns it. The basis functions
            !   are laid out site by site and, within a site, state by
            !   state, so walking valeDim directly visits them in exactly
            !   that order.
            do basisFn = 1, valeDim
               conjWaveMomSum(k,pOptcIndex(basisFn),finalStateIndex,i) = &
                     & conjWaveMomSum(k,pOptcIndex(basisFn), &
                     & finalStateIndex,i) &
                     & + (conjg(valeVale(basisFn,j,1)) &
                     & * valeValeMM(basisFn,k,i))
            enddo
         enddo
      enddo
   enddo

#else

   ! Documentation similar to the above non-gamma case.
   do i = initComponent, finComponent

      ! Make the upper triangle correct for Hermiticity.  Recall that for
      !   the Gamma K Point all the matrices are real (except the momentum
      !   matrix which was multiplied by a -i and is hence imaginary).
      !   Since it must be Hermitian we need to apply that now.
      do j = 1, valeDim
         valeValeMMGamma(1:j,j,i) = -valeValeMMGamma(1:j,j,i)
      enddo

      finalStateIndex = 0
      do j = firstFin, lastFin

         ! Increment the finalStateIndex for conjWaveMomSum
         finalStateIndex = finalStateIndex + 1

         do k = 1, valeDim
            do basisFn = 1, valeDim
               conjWaveMomSum(k,pOptcIndex(basisFn),finalStateIndex,i) = &
                     & conjWaveMomSum(k,pOptcIndex(basisFn), &
                     & finalStateIndex,i) &
                     & + (valeValeGamma(basisFn,j,1) &
                     & * valeValeMMGamma(basisFn,k,i))
            enddo
         enddo
      enddo
   enddo

#endif

end subroutine buildConjWaveMomSumPOPTC



! Build the decomposition index that assigns every basis function to a
!   partial, and the tables that describe the resulting layout.
!   PSEUDOCODE 18, DESIGN 11.
!
! The answer does not depend on the k-point, the spin, or anything else
!   that varies during a run: it is fixed by the structure and the
!   requested detail code. So this is called ONCE, before the k-point
!   loop, and the arrays it fills are read by whichever transition
!   producer the integration method selected. It used to be rebuilt on
!   every call of computePOPTCPairs, which was wasted work even then and
!   is not available at all to the tetrahedron pathway, whose store must
!   be sized from sumNumPartials before the loop begins.
subroutine buildPOPTCIndex

   ! Import the necessary modules.
   use O_Kinds
   use O_Constants,   only: lAngMomCount
   use O_KPoints,     only: numPointOps
   use O_AtomicSites, only: valeDim, numAtomSites, atomSites, atomPerm
   use O_AtomicTypes, only: numAtomTypes, atomTypes
   use O_Input,       only: detailCodePOPTC

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define local variables.
   integer :: i,j,k,l ! Loop index variables.
   integer :: currentType
   integer :: valeDimIndex
   integer :: opIdx  ! Point operation being tabulated.
   integer :: siteRot ! Image of a site under that operation.

   ! Variables that resolve the decomposition request into the two
   !   independent parameters of DESIGN 11.2 and then lay the partials
   !   out. Writing the assignment as one parameterized walk rather than
   !   one branch per detail code keeps the grid of DESIGN 11 visible in
   !   the code, so that a cell added later is a new parameter value
   !   rather than a new branch that has to remember to do everything
   !   the other branches do.
   logical :: groupByType  ! Codes 1,2 group by type; codes 3,4 by site
   logical :: resolveTotal ! Codes 1,3 give a segment one shared slot
   integer :: numSegments  ! Type count when grouping by type, else sites
   integer :: typeOfSegment    ! The type whose basis a segment carries
   integer :: slotsThisSegment ! Partials the segment being laid out owns
   integer :: currentSegment   ! The segment that the current site feeds
   integer :: currentSlot      ! Slot within that segment
   integer :: currentPartial   ! The partial those two resolve to

   ! How many basis functions feed each partial. Local because only the
   !   Kramers-Kronig factor written below consumes it.
   real (kind=double), allocatable, dimension (:) :: partialsIndex

   ! Resolve the requested decomposition into the two independent
   !   parameters that DESIGN 11.2 defines a cell by. (Detail code 0
   !   means no decomposition at all and never reaches this routine.)
   if (detailCodePOPTC <= 2) then
      groupByType = .true.  ! Codes 1 and 2: a segment is an atomic type.
   else
      groupByType = .false. ! Codes 3 and 4: a segment is an atomic site.
   endif
   if ((detailCodePOPTC == 1) .or. (detailCodePOPTC == 3)) then
      resolveTotal = .true.  ! Codes 1 and 3: one partial per segment.
   else
      resolveTotal = .false. ! Codes 2 and 4: one per radial function.
   endif

   ! The number of segments follows from the grouping alone.
   if (groupByType) then
      numSegments = numAtomTypes
   else
      numSegments = numAtomSites
   endif

   allocate (segmentBase     (numSegments + 1))
   allocate (slotsPerSegment (numSegments))

   ! Lay the partials out. Each segment owns a contiguous block of them:
   !   segmentBase records where that block starts as a zero based
   !   offset, and slotsPerSegment records how long it is. Both outlive
   !   this walk because the star unfolding needs them to carry a partial
   !   from one site onto the block of the site's image under a symmetry
   !   operation.
   segmentBase(1) = 0
   do i = 1, numSegments

      ! When grouping by type the segment index already is a type index.
      !   When grouping by site it is a site index, and the basis that
      !   site carries is determined by its type.
      if (groupByType) then
         typeOfSegment = i
      else
         typeOfSegment = atomSites(i)%atomTypeAssn
      endif

      ! A total resolved segment holds a single partial that everything
      !   in the segment shares. An nl resolved one holds a partial per
      !   radial function, summed over the s, p, d and f shells. Note
      !   that this counts QN_nl pairs and not basis functions: all m
      !   components of a shell feed one partial, which is exactly the
      !   property that lets the result be unfolded from an irreducible
      !   wedge with the atom permutation alone (DESIGN 2.5).
      if (resolveTotal) then
         slotsThisSegment = 1
      else
         slotsThisSegment = &
               & sum(atomTypes(typeOfSegment)%numQN_lValeRadialFns(:))
      endif

      slotsPerSegment(i) = slotsThisSegment
      segmentBase(i+1)   = segmentBase(i) + slotsThisSegment
   enddo

   ! The final base sits one past the last block, so it is the count.
   sumNumPartials = segmentBase(numSegments + 1)

   ! For each sub-group, record the number of basis functions that
   !   contribute to it so that we can properly normalize the KKC.
   allocate (partialsIndex(sumNumPartials))
   partialsIndex(:) = 0

   ! Allocate storage for pOptcIndex so that each basis function can be
   !   "sent" (or indexed) to the correct accumulation group.
   allocate (pOptcIndex (valeDim))

   ! Track which basis function is currently under consideration for the
   !   mapping into pOptcIndex. The valeDim is "valence dimension", one
   !   for each basis function (orbital).
   valeDimIndex = 0

   ! Loop over every atom in the system to index where the pOptc values
   !   for each atom should be stored.
   do i = 1, numAtomSites

      ! Obtain the type of the current atom.
      currentType = atomSites(i)%atomTypeAssn

      ! Identify the segment that this site feeds. Grouping by type
      !   sends every atom of a type to one segment; grouping by site
      !   gives each atom its own. The walk is over sites either way,
      !   because that is the order the basis functions are laid out in,
      !   so only the destination changes between the two.
      if (groupByType) then
         currentSegment = currentType
      else
         currentSegment = i
      endif

      ! One loop over the angular momentum shells covers s, p, d and f
      !   together: the number of m components of shell l is (l-1)*2+1,
      !   which gives 1, 3, 5 and 7 for l = 1 through 4. A total
      !   resolved segment holds every radial function on its one shared
      !   slot, while an nl resolved one advances the slot per radial
      !   function, so its slots run in s, p, d, f order. That ordering
      !   matters beyond tidiness: printSpectrumPOPTC walks the partials
      !   in this same layout order, and the sequence numbers written
      !   into the output file are the only thing tying a spectrum to
      !   its label.
      currentSlot = 0
      do j = 1, lAngMomCount ! 1=s; 2=p; 3=d; 4=f
         do k = 1, atomTypes(currentType)%numQN_lValeRadialFns(j)

            if (resolveTotal) then
               currentSlot = 1
            else
               currentSlot = currentSlot + 1
            endif

            currentPartial = segmentBase(currentSegment) + currentSlot

            do l = 1, (j-1)*2 + 1
               valeDimIndex = valeDimIndex + 1
               pOptcIndex(valeDimIndex) = currentPartial

               ! Record how many basis functions feed this partial, which
               !   is what imagoKKc normalizes the additive constant of
               !   eps1 with. Counting by an increment per assigned basis
               !   function is correct for every cell, including the type
               !   grouped ones where many sites feed one partial.
               partialsIndex(currentPartial) = &
                     & partialsIndex(currentPartial) + 1.0_double
            enddo
         enddo
      enddo
   enddo ! i = 1, numAtomSites

   ! Build the table that carries each partial onto its image under each
   !   point operation (PSEUDOCODE 7a). Only the atom grouped codes need
   !   it. The type grouped codes are already correct on a reduced mesh,
   !   because every operation carries an atom onto an atom of the same
   !   type and a type level sum therefore maps onto itself. The guard on
   !   atomPerm covers style code 0, where an explicit k-point list
   !   leaves Imago no symmetry from which to build the maps at all.
   !
   ! A partial is a slot within a segment, and for these codes a segment
   !   is an atomic site, so carrying a partial through an operation
   !   means re-basing its slot onto the block that belongs to the
   !   site's image. The slot number itself survives unchanged, and that
   !   is what makes this work: buildAtomPerm only ever maps an atom onto
   !   an atom of the same type, and atoms of one type share a basis, so
   !   the image site's slots stand in one to one correspondence with the
   !   original's and carry the same QN_nl meaning. For detail code 3
   !   there is exactly one slot per segment and the table reduces to
   !   atomPerm itself, which is why the star average is written once for
   !   both codes rather than special cased.
   if ((detailCodePOPTC >= 3) .and. (allocated(atomPerm))) then
      allocate (partialPerm (numPointOps, sumNumPartials))

      do opIdx = 1, numPointOps
         do i = 1, numAtomSites
            siteRot = atomPerm(opIdx,i)
            do j = 1, slotsPerSegment(i)
               partialPerm(opIdx, segmentBase(i) + j) = &
                     & segmentBase(siteRot) + j
            enddo
         enddo
      enddo
   endif

   ! Write the KKC factor to file for use in imagoKKc. Written once,
   !   which is now simply a consequence of this routine being called
   !   once rather than something a k-point test has to arrange.
   write (209,fmt="(a17,a5)") 'POPTC_KKC_FACTOR ', '    1'
   do j = 1, sumNumPartials
      do  k = 1, sumNumPartials
      write (209,fmt="(a17,1e15.7)") 'POPTC_KKC_FACTOR ', &
            & (partialsIndex(j) * partialsIndex(k) / (valeDimIndex**2))
      enddo
   enddo

   deallocate (partialsIndex)

end subroutine buildPOPTCIndex


! Release what buildPOPTCIndex allocated. Separate from the build so
!   that the caller owning the k-point loop owns the lifetime, rather
!   than a producer inside the loop having to guess whether it is the
!   last one to run.
subroutine cleanUpPOPTCIndex

   implicit none

   if (allocated(segmentBase))     deallocate (segmentBase)
   if (allocated(slotsPerSegment)) deallocate (slotsPerSegment)
   if (allocated(pOptcIndex))      deallocate (pOptcIndex)

   ! Built only for the atom grouped codes, and only when the symmetry
   !   maps exist, so its release is guarded the same way.
   if (allocated(partialPerm))     deallocate (partialPerm)

end subroutine cleanUpPOPTCIndex


subroutine computePOPTCPairs(currentKPoint,xyzComponents,spinDirection,doOPTC)

   ! Import the necessary modules.
   use O_Kinds
   use O_TimeStamps
   use O_Potential,   only: spin
   use O_SortSubs,    only: mergeSort
   use O_Populate,    only: electronPopulation
   use O_KPoints,     only: kPointWeight, numKPoints, numFullMeshKP, &
         & fullKPToIBZKPMap, fullKPToIBZOpMap, kPointIntgCode
   use O_Constants,   only: pi, hartree, dim3
   use O_AtomicSites, only: valeDim, atomPerm
   use O_Input,       only: numStates, detailCodePOPTC

   ! The momentum matrix itself is read by buildConjWaveMomSumPOPTC rather
   !   than here, so only the wave functions are needed at this level.
#ifndef GAMMA
   use O_SecularEquation, only: valeVale, energyEigenValues
#else
   use O_SecularEquation, only: valeValeGamma, energyEigenValues
#endif


   ! Make sure that no funny variables.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   integer :: currentKPoint
   integer :: xyzComponents ! 0=all, 1=x, 2=y, 3=z
   integer :: spinDirection
   integer, intent(in) :: doOPTC

   ! Define local variables specific to POPTC. The decomposition index
   !   and its layout tables are built by buildPOPTCIndex before the
   !   k-point loop, so none of the variables that construct them appear
   !   here any more.
   integer :: i,j,k,l,n,o ! Loop index variables

   ! Variables for the IBZ star unfolding of the atom grouped pair
   !   matrix (PSEUDOCODE 7a). These are used only for detail codes 3
   !   and 4, and the indices they carry are PARTIALS rather than atoms:
   !   for code 3 a partial is exactly an atomic site, but for code 4 a
   !   partial is one QN_nl slot within a site, so the permutation that
   !   acts on the pair matrix is partialPerm and not atomPerm.
   integer :: starSize   ! Full-mesh k-points folding onto this IBZ point
   integer :: fullIdx    ! Index of a full-mesh k-point in that star
   integer :: opIdx      ! Point operation carrying the IBZ point to it
   integer :: pairIndex  ! Transition pair within this k-point
   integer :: component  ! Cartesian component of the momentum operator
   integer :: partialA   ! Initial-state partial index of the pair
   integer :: partialB   ! Final-state partial index of the pair
   integer :: permPartialA ! Image of partialA under the operation
   integer :: permPartialB ! Image of partialB under the operation
   real (kind=double), allocatable, dimension (:,:) :: pairSlabSym
   real (kind=double), allocatable, &
         & dimension (:,:,:,:) :: transitionProbTemp

   ! Define local variables that are the same as computePairs.
   integer :: initComponent
   integer :: finComponent
   integer :: transPairCount
   integer :: firstInit
   integer :: lastInit
   integer :: firstFin
   integer :: lastFin
   integer :: finalStateIndex
   integer :: orderedIndex
   integer, allocatable, dimension (:) :: sortOrder
   integer, allocatable, dimension (:) :: segmentBorders
   real    (kind=double) :: initStateFactor
   real    (kind=double) :: finStatefactor
   real    (kind=double) :: currentEnergyDiff
   real    (kind=double), allocatable, dimension (:)       :: energyDiffTemp
#ifndef GAMMA
   real    (kind=double) :: valeValeXMomSumReal
   real    (kind=double) :: valeValeXMomSumImag
   complex (kind=double), allocatable, dimension (:,:,:,:) :: conjWaveMomSum
   complex (kind=double), allocatable, dimension (:,:,:)   :: valeValeXMom
#else
   real    (kind=double) :: valeValeXMomGammaSum
   real  (kind=double), allocatable, dimension (:,:,:,:) :: conjWaveMomSumGamma
   real    (kind=double), allocatable, dimension (:,:,:)   :: valeValeXMomGamma
#endif


   ! The decomposition index is built once before the k-point loop, by
   !   buildPOPTCIndex, so this routine reads pOptcIndex, segmentBase,
   !   slotsPerSegment, sumNumPartials and partialPerm rather than
   !   constructing them. It does not own them and must not free them.

   ! Storage for the transition probabilities of all partial pairs. This
   !   is the GAUSSIAN pathway's store; the tetrahedron pathway fills
   !   transProbPOPTCBanded instead, under a band-pair index.
   if (.not. allocated(transitionProbPOPTC)) then
      allocate(transitionProbPOPTC(sumNumPartials,sumNumPartials,dim3,&
            & maxPairs,numKPoints,spin))
      transitionProbPOPTC (:,:,:,:,:,:) = 0.0_double
   endif

   ! Make shorthand for the state indices.
   firstInit = firstOccupiedState(currentKPoint,spinDirection)
   lastInit  = lastOccupiedState(currentKPoint,spinDirection)
   firstFin  = firstUnoccupiedState(currentKPoint,spinDirection)
   lastFin   = lastUnoccupiedState(currentKPoint,spinDirection)

   ! Initialize a counter for the current number of transition pairs
   transPairCount = 0

   ! Determine the range of components (xyz) that should be considered.
   if (xyzComponents == 0) then
      initComponent = 1
      finComponent = 3
   else
      initComponent = 1
      finComponent = 1
   endif

   ! Allocate space to hold the sum(conjg(valeVale(:,j)) * valeVale_Mom(:,k,1))
   !   for each of the possible final states, resolved by partial.  This is
   !   done since the values are independent of the initial states.  The
   !   finComponent is 3 for all three at once, and 1 for when X, Y, Z are
   !   done separately.
   !
   ! The sum itself is built by buildConjWaveMomSumPOPTC, which both
   !   Brillouin-zone integration methods share (DESIGN 12.2).
#ifndef GAMMA
   allocate (conjWaveMomSum (valeDim,sumNumPartials,lastFin-firstFin+1,&
                             & finComponent))

   call buildConjWaveMomSumPOPTC (firstFin,lastFin,initComponent, &
         & finComponent,conjWaveMomSum)
#else
   allocate (conjWaveMomSumGamma (valeDim,sumNumPartials,lastFin-firstFin+1,&
                                  & finComponent))

   call buildConjWaveMomSumPOPTC (firstFin,lastFin,initComponent, &
         & finComponent,conjWaveMomSumGamma)
#endif


   ! Allocate space for the energy difference.
   allocate (energyDiffTemp (maxPairs))
   allocate (transitionProbTemp (sumNumPartials,sumNumPartials,finComponent,&
                                 & maxPairs))

   ! Initialize the temporary energy transition array.
   energyDiffTemp(:) = 0.0_double

   ! Allocate space to hold the indices for each segment of the energyDiff
   !   array.
   allocate (segmentBorders (lastInit-firstInit+2))

   ! Initialize the first index since it will always be 0.
   segmentBorders(1) = 0

#ifndef GAMMA
     allocate (valeValeXMom (sumNumPartials, sumNumPartials, finComponent))
#else
     allocate (valeValeXMomGamma (sumNumPartials,sumNumPartials,finComponent))
#endif


   ! Begin the double loop to determine the transition energies.
   do i = firstInit, lastInit
      do j = firstFin, lastFin

         ! Index into conjWaveMomSum, which was filled above for *every*
         !   final state in the firstFin to lastFin range. The index is
         !   derived from j rather than accumulated by a counter because
         !   the loop below skips some j: with thermal smearing a state
         !   can be both initial and final, and those skipped j still
         !   occupy their slot in conjWaveMomSum. A counter incremented
         !   only on accepted pairs would fall behind at the first skip
         !   and read a different final state's momentum sum from then
         !   on, silently, for the rest of this initial state.
         finalStateIndex = j - firstFin + 1

         ! Recall that thermal smearing may allow some states to be both
         !   initial and final. We do not consider transitions where the final
         !   state has an energy less than the initial.
         if (i >= j) cycle

         ! If the energy of the final state is higher than the requested
         !   cut-off we go to the next initial state.
         if (energyEigenValues(j,currentKPoint,spinDirection) > &
               & energyCutoff) exit

         ! Compute the energy of the transition from the current states.
         currentEnergyDiff = energyEigenValues(j,currentKPoint,spinDirection)-&
               & energyEigenValues(i,currentKPoint,spinDirection)

         ! Check if the energy difference is less than the maximum
         !   transition energy that the input file requested computation
         !   for.  If it fails, then we go to the next initial state because
         !   all the remaining final states for this energy will be greater.
         if (currentEnergyDiff > maxTransEnergy) exit

         ! Increment the number of transition pairs counted so far.
         transPairCount = transPairCount + 1

         ! Store the transition energy for the current pair.
         energyDiffTemp(transPairCount) = currentEnergyDiff

         ! In the event that thermal smearing is turned on. The state that the
         !   e- comes from and goes into may be fully, partially, or not
         !   occupied. We will scale the probability of a transition linearly
         !   according to the percent occupation of both the initial and final
         !   states.

         ! Determine the array index value of the current initial (index i)
         !   spin-kpoint-state as defined by the tempEnergyEigenValues loop
         !   near the beginning of the population subroutine.
         orderedIndex = i + numStates*(spinDirection-1) + &
               & numStates*spin*(currentKPoint-1)

         ! Use the normal state factor for non-PACS calculations. For PACS
         !   calculations the initStateFactor is always 1 even though the
         !   initial core state(s) will have an electron missing.
         if (doOPTC /= 2) then
            initStateFactor = electronPopulation(orderedIndex) / &
                  & (kPointWeight(currentKPoint)/real(spin,double))
         else
            initStateFactor = 1.0_double
         endif

         ! Determine the array index value of the current final (index j)
         !   spin-kpoint-state as defined by the tempEnergyEigenValues loop
         !   near the beginning of the population subroutine.
         orderedIndex = j + numStates*(spinDirection-1) + &
               & numStates*spin*(currentKPoint-1)

         finStateFactor = 1.0_double - electronPopulation(orderedIndex) / &
               & (kPointWeight(currentKPoint)/real(spin,double))


#ifndef GAMMA
        
         ! Loop to obtain the wave function times the momentum integral.
         do k = initComponent,finComponent

            valeValeXMom (:,:,:) = cmplx(0.0_double,0.0_double,double)
                   
            initVDBI = 0

            do l = 1, valeDim
               initVDBI = initVDBI + 1

               do n = 1,sumNumPartials
                  valeValeXMom(pOptcIndex(initVDBI),n,k) = &
                        & valeValeXMom(pOptcIndex(initVDBI),n,k) &
                        & + valeVale(initVDBI,i,1) &
                        & * conjWaveMomSum(initVDBI,n,finalStateIndex,k)
               enddo
            enddo

            ! Initialize 
            valeValeXMomSumReal = 0
            valeValeXMomSumImag = 0

            ! The sum of the real and imaginary parts of valeValeXMom. These
            !   are for the transition probability to account for the change 
            !   from sum squared to sum of squares
            valeValeXMomSumReal = sum(real(valeValeXMom(:,:,k),double))
            valeValeXMomSumImag = sum(aimag(valeValeXMom(:,:,k)))
            do n = 1, sumNumPartials
               do o = 1,sumNumPartials

                  transitionProbTemp(o,n,k,transPairCount) = &
                        & ((real(valeValeXMom(o,n,k))*valeValeXMomSumReal) &
                        & + (aimag(valeValeXMom(o,n,k))*valeValeXMomSumImag)) &
                        & * initStateFactor*finStateFactor
               enddo
            enddo
         enddo

#else
         ! Loop to get the wave function times the momentum matrix element.

         do k = initComponent,finComponent

            valeValeXMomGamma (:,:,:) = 0.0_double
            initVDBI = 0

            do l = 1, valeDim
               initVDBI = initVDBI + 1

               do n = 1,sumNumPartials
                  valeValeXMomGamma(pOptcIndex(initVDBI),n,k) &
                        & = valeValeXMomGamma(pOptcIndex(initVDBI),n,k) &
                        & + valeValeGamma(initVDBI,i,1) &
                        & * conjWaveMomSumGamma(initVDBI,n,finalStateIndex,k)

               enddo
            enddo

            ! Initialize
            valeValeXMomGammaSum = 0
           
            ! The sum of all valeValeXMomGamma. Used in the transition 
            !   probability to account from the change from sum squared 
            !   to the sum of squares.
            valeValeXMomGammaSum = sum(valeValeXMomGamma(:,:,k))
            do n = 1, sumNumPartials
               do o = 1,sumNumPartials

                  transitionProbTemp(o,n,k,transPairCount) = &
                        & valeValeXMomGamma(o,n,k)*valeValeXMomGammaSum &
                        & * initStateFactor*finStateFactor
               enddo
            enddo
         enddo
#endif
      enddo ! Fin loop j

      ! Save the index for the end border of this segment.
      segmentBorders(i - firstInit + 2) = transPairCount
   enddo ! Init loop i

   ! Deallocate the matrices that the transition loop needed. Both are
   !   finished with: the loop above closed at the "Init loop i" line and
   !   nothing below this point reads either of them. The release is
   !   written out explicitly even though a local allocatable is handed
   !   back automatically when the subroutine returns, because this
   !   routine is called once per k-point per spin and pairing each
   !   deallocate with its allocate is what lets a reader confirm the
   !   per-call cost is actually returned without having to rely on
   !   knowing that language rule.
#ifndef GAMMA
   deallocate (conjWaveMomSum)
   deallocate (valeValeXMom)
#else
   deallocate (conjWaveMomSumGamma)
   deallocate (valeValeXMomGamma)
#endif

   ! ----------------------------------------------------------------------
   ! Unfold the atom grouped pair matrix over the star of this IBZ
   !   k-point. (PSEUDOCODE 7a, DESIGN 2.5.)
   !
   ! Detail codes 3 and 4 need this and codes 1 and 2 do not. The
   !   numbering of DESIGN 11.3 puts grouping ahead of resolution, so the
   !   test is a threshold rather than a memorized set of cases: the low
   !   codes group by TYPE, and every operation that Imago reduces the
   !   mesh by carries each atom onto an atom of the same type, so a type
   !   level sum already maps onto itself and is correct as computed.
   !   (That closure is enforced at startup by buildAtomPerm rather than
   !   assumed from a type being a symmetry orbit, which in an amorphous
   !   cell or a defect supercell it is not. See DESIGN 2.3.)
   !
   ! The resolution axis does not enter. Both offered resolutions sum
   !   over complete shells -- a whole group, or a QN_nl radial function
   !   summed over its m components -- and that is exactly the condition
   !   the invariance argument of DESIGN 2.3 needs. A QN_nlm resolution
   !   would break it and require the deferred D^l(R) representation
   !   matrices, which is one of the two reasons DESIGN 11.2 does not
   !   offer it. So nothing here is deferred: every offered cell is
   !   either correct as computed or made correct in this block.
   !
   ! What the star average is exact for. The momentum operator is a
   !   vector, so an operation both relabels the atoms and mixes the
   !   Cartesian components, P_c(Rk) = sum_d R_cd P_d(k). The permutation
   !   below handles the relabeling and not the mixing. The mixing cancels
   !   when the three components are summed, because R is orthogonal:
   !   sum_c R_cd R_ce = delta_de kills every cross term, and the
   !   component-summed pair matrix therefore transforms by pure index
   !   permutation. So the isotropic column that printSpectrumPOPTC writes
   !   becomes correct per atom pair on a reduced mesh, while the separate
   !   x, y and z columns stay exactly as unverified as they already are
   !   for the total spectra. Repairing those is a separate question that
   !   cannot be answered on this quantity at all -- it has to be answered
   !   on the complex matrix element before the probability is formed.
   !
   ! Note that the total spectra do not move by so much as a bit. The
   !   total is the sum of this matrix over both indices, and permuting or
   !   averaging permuted copies does not change a sum. That is also why
   !   the identity that the partials sum to the total cannot be used to
   !   check any of this: it holds equally whether the unfolding is
   !   present, absent, or wrong. Check instead by running a structure
   !   whose symmetry-equivalent atoms are inequivalently oriented on a
   !   full mesh and on a reduced mesh, and requiring the two to agree per
   !   atom pair.
   !
   ! This belongs to the GAUSSIAN pathway alone, which is why the guard
   !   below tests the integration code. That pathway visits only
   !   irreducible k-points and must spread each one's contribution over
   !   the members of its star. The tetrahedron pathway visits full-mesh
   !   corners directly and permutes once per corner as it fetches, so
   !   this block does not run there at all -- doing both would count the
   !   symmetry twice. Both pathways are correct under reduction; they
   !   reach it by different routes (PSEUDOCODE 19.5, DESIGN 12.6).
   ! ----------------------------------------------------------------------
   if ((kPointIntgCode == 0) .and. (detailCodePOPTC >= 3) .and. &
         & (allocated(atomPerm))) then

      ! Count the star: the number of full-mesh k-points that fold onto
      !   this IBZ representative. Counted the same way as in computeBond.
      starSize = 0
      do fullIdx = 1, numFullMeshKP
         if (fullKPToIBZKPMap(fullIdx) == currentKPoint) then
            starSize = starSize + 1
         endif
      enddo

      ! Scratch space for one symmetrized slab. The averaging cannot be
      !   done in place because each star member reads the whole original
      !   slab while writing scattered elements of the result. The slab
      !   is indexed by PARTIAL rather than by site: the two coincide for
      !   detail code 3, where each site owns one partial, but for code 4
      !   a site owns one partial per radial function and a slab bounded
      !   by the site count would be far too small.
      allocate (pairSlabSym(sumNumPartials,sumNumPartials))

      ! The star sum is a fixed linear map on the pair matrix: it does not
      !   depend on energy, and the Gaussian broadening applied later in
      !   getOptcCondPOPTC is linear. So it is applied once per transition
      !   per k-point, here, and the weighted accumulation over IBZ points
      !   is left alone. Doing it inside the energy loop instead would
      !   multiply the innermost work by the reduction factor of 4 to 48
      !   for the very same answer.
      do pairIndex = 1, transPairCount
         do component = initComponent, finComponent

            pairSlabSym(:,:) = 0.0_double

            ! Walk the star. Each member contributes the IBZ slab with
            !   both of its partial indices carried through that member's
            !   operation, and the division by starSize turns the sum into
            !   an average, so the totals are preserved exactly.
            do fullIdx = 1, numFullMeshKP
               if (fullKPToIBZKPMap(fullIdx) /= currentKPoint) cycle
               opIdx = fullKPToIBZOpMap(fullIdx)

               do partialA = 1, sumNumPartials
                  permPartialA = partialPerm(opIdx,partialA)

                  do partialB = 1, sumNumPartials
                     permPartialB = partialPerm(opIdx,partialB)

                     pairSlabSym(permPartialA,permPartialB) = &
                           & pairSlabSym(permPartialA,permPartialB) &
                           & + transitionProbTemp(partialA,partialB, &
                           & component,pairIndex) / real(starSize,double)
                  enddo ! partialB
               enddo ! partialA
            enddo ! fullIdx (star members)

            transitionProbTemp(:,:,component,pairIndex) = pairSlabSym(:,:)
         enddo ! component
      enddo ! pairIndex

      deallocate (pairSlabSym)
   endif

   ! Determine if there was only one segment.  In this case we don't have to
   !   sort anything.

   ! Sort energyDiffTemp into energyDiff, and obtain the indices for the
   !   correct sorted order of energyDiff so that we can copy the energy
   !   momentum directly.

   allocate (sortOrder (transPairCount))

   call mergeSort (energyDiffTemp,energyDiff(:,currentKPoint,spinDirection),&
         & sortOrder,segmentBorders,transPairCount)

   ! Copy transitionProbTemp to the real transitionProb data structure using
   !   the sorting order determined in the mergeSort subroutine.
   if (xyzComponents == 0) then
      do i = 1, transPairCount
         transitionProbPOPTC(:,:,:,i,currentKPoint,spinDirection) = &
               & transitionProbTemp(:,:,:,sortOrder(i))
         do j = 1, dim3
            transitionProb(j,i,currentKPoint,spinDirection) = &
                  & sum(transitionProbPOPTC(:,:,j,i,currentKPoint,&
                  & spinDirection))
         enddo
      enddo
   else
      do i = 1, transPairCount
         transitionProbPOPTC(:,:,xyzComponents,i, &
               & currentKPoint,spinDirection) = &
               & transitionProbTemp(:,:,1,sortOrder(i))
         transitionProb(xyzComponents,i,currentKPoint,spinDirection) = &
               & sum(transitionProbPOPTC(:,:,xyzComponents,i,currentKPoint,&
               & spinDirection))
      enddo
   endif

   ! Deallocate unnecessary arrays and matrices. The decomposition index
   !   is NOT among them: it is built once before the k-point loop and
   !   released after it, because both this routine and the tetrahedron
   !   producer read it and neither owns it.
   deallocate (energyDiffTemp)
   deallocate (transitionProbTemp)
   deallocate (segmentBorders)
   deallocate (sortOrder)

end subroutine computePOPTCPairs


subroutine computeSigmaE (currentKPoint,xyzComponents,spinDirection)

   ! Import the necessary data modules.
   use O_Kinds
   use O_Constants, only: dim3, pi, auTime, lightFactor, hartree
   use O_KPoints, only: numKPoints, kPointWeight
   use O_SortSubs
   use O_Input, only: maxTransEnSIGE, deltaSIGE, sigmaSIGE, numStates
   use O_Populate, only: occupiedEnergy, electronPopulation
   use O_Lattice, only: realCellVolume
   use O_Potential, only: spin
   use O_AtomicSites, only: valeDim
#ifndef GAMMA
   use O_SecularEquation, only: valeVale, valeValeMM, energyEigenValues
#else
   use O_SecularEquation, only: valeValeGamma, valeValeMMGamma, &
         & energyEigenValues
#endif

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   integer :: currentKPoint
   integer :: xyzComponents
   integer :: spinDirection

   ! Define local variables.
   integer :: i,j,k ! Loop index variables
   integer :: initComponent
   integer :: finComponent
   integer :: firstInit
   integer :: lastInit
   integer :: firstFin
   integer :: lastFin
   integer :: numInitEnergyStates
   integer :: numFinEnergyStates
   integer :: numTotalEnergyStates
   integer :: finalStateIndex
   integer :: initialStateIndex
   integer :: numEnergyPoints
   integer :: orderedIndex
   real (kind=double) :: alphaFactor
   real (kind=double) :: initStateFactor
   real (kind=double) :: finStateFactor
   real (kind=double) :: kPointFactor
   real (kind=double) :: sigmaSqrt2Pi
   real (kind=double) :: conversionFactor
   real (kind=double) :: totalSigma
   real    (kind=double), allocatable, dimension (:,:,:) :: transitionProb
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:,:) :: conjWaveMomSum
   complex (kind=double)                                 :: valeValeXMom
#else
   real    (kind=double), allocatable, dimension (:,:,:) :: conjWaveMomSumGamma
   real    (kind=double)                                 :: valeValeXMomGamma
#endif

   ! Compute the number of energy points to evaluate.
   numEnergyPoints = int((maxTransEnSIGE * 2.0_double) / deltaSIGE + 1)

   ! Define constants for normalizing the broadening gaussian.
   sigmaSqrt2Pi = sigmaSIGE * sqrt(2.0_double * pi)

   ! Obtain the state bounds for this kpoint.
   firstInit =   firstOccupiedState(currentKPoint,spinDirection)
   lastInit  =    lastOccupiedState(currentKPoint,spinDirection)
   firstFin  = firstUnoccupiedState(currentKPoint,spinDirection)
   lastFin   =  lastUnoccupiedState(currentKPoint,spinDirection)

   ! Get the number of energy states in the system.
   numFinEnergyStates   = lastFin  - firstFin  + 1
   numInitEnergyStates  = lastInit - firstInit + 1
   numTotalEnergyStates = lastFin  - firstInit + 1

   ! Allocate space for the resulting sigmaE and initialize it to zero during
   !   the first kpoint iteration.  Also create the energy scale and initialize
   !   it.  NOTE:  The sigmaEAccumulator
   if (.not. allocated(sigmaEAccumulator)) then
      allocate (sigmaEAccumulator (numEnergyPoints,dim3))
      sigmaEAccumulator (:,:) = 0.0_double

      allocate (energyScale (numEnergyPoints))
      do i = 1, numEnergyPoints
         energyScale(i) = -maxTransEnSIGE + deltaSIGE * (i-1) + occupiedEnergy
      enddo
   endif

   ! Determine the range of components (xyz) that should be considered.
   if (xyzComponents == 0) then
      initComponent = 1
      finComponent = 3
   else
      initComponent = 1
      finComponent = 1
   endif

#ifndef GAMMA
   ! Allocate space to hold the sum(conjg(valeVale(:,j)) * valeVale_Mom(:,k,1))
   !   for each of the possible energy states.
   allocate (conjWaveMomSum (valeDim,numFinEnergyStates,finComponent))

   ! Compute the sum over all the energy states.
   do i = initComponent,finComponent

      finalStateIndex = 0
      do j = firstFin, lastFin

         ! Increment the final state index for conjWaveMomSum. Note that this
         !   will compute the conjWaveMomSumGamma for every possible final
         !   state wihtin the firstFin - lastFin range. Later on, we may find
         !   that for certain initial states we don't actually need every final
         !   state. In those cases we will cycle past the final states.
         finalStateIndex = finalStateIndex + 1

         do k = 1, valeDim
            conjWaveMomSum(k,finalStateIndex,i) = &
                  & sum(conjg(valeVale(:,j,1) * valeValeMM(:,k,i)))
         enddo
      enddo
   enddo

#else

   allocate (conjWaveMomSumGamma (valeDim,numFinEnergyStates,finComponent))

   ! Compute the sum over all the energy states.
   do i = initComponent,finComponent

      finalStateIndex = 0
      do j = firstFin, lastFin

         ! Increment the final state index for conjWaveMomSumGamma. Note that
         !   this will compute the conjWaveMomSumGamma for every possible final
         !   state wihtin the firstFin - lastFin range. Later on, we may find
         !   that for certain initial states we don't actually need every final
         !   state. In those cases we will cycle past the final states.
         finalStateIndex = finalStateIndex + 1

         do k = 1, valeDim
            conjWaveMomSumGamma(k,finalStateIndex,i) = &
                  & sum(valeValeGamma(:,j,1) * valeValeMMGamma(:,k,i))
         enddo
      enddo
   enddo

#endif


   ! Allocate space for the fully computed set of transition probabilties
   !   between all of the involved states.  FinComponent is either 1 or 3.
   allocate (transitionProb (finComponent,numTotalEnergyStates,&
         & numTotalEnergyStates))

   ! Begin the double loop to determine the momentum matrix elements between
   !   each possible state.
   ! Initialize the counter for the initial state index number. (The point is
   !   that the firstInit and lastInit values could be something like 345 and
   !   400.  We want to index an array from 1 to 56 so we use counters like
   !   this.)
   initialStateIndex = 0
   do i = firstInit, lastInit

      ! Increment the index number for the initial states.
      initialStateIndex = initialStateIndex + 1

      ! Initialize the counter for the final state index number.
      finalStateIndex = 0
      do j = firstFin, lastFin

         ! Define the indices for the conjWaveMomSum and the momentum matrix
         !   elements. Note that we need to increment the index for every
         !   j-loop iteration because we computed the conjWaveMomSum and
         !   conjWaveMomSumGamma for every final state. We have to increment
         !   the counter for every case, even if we don't use it.  (This note
         !   is here because this was a point of confusion in the past when the
         !   code was being developed.)
         finalStateIndex = finalStateIndex + 1

         ! We don't want double counting or self interactions so we skip the
         !   i>=j cases.
         if (i >= j) cycle

         ! In the event that thermal smearing is turned on. The state that the
         !   e- comes from and goes into may be fully, partially, or not
         !   occupied. We will scale the probability of a transition linearly
         !   according to the percent occupation of both the initial and final
         !   states.

         ! Determine the array index value of the current initial (index i)
         !   spin-kpoint-state as defined by the tempEnergyEigenValues loop
         !   near the beginning of the population subroutine.
         orderedIndex = i + numStates*(spinDirection-1) + &
               & numStates*spin*(currentKPoint-1)

         ! The initial state factor is the population of this state divided by
         !   the kpoint weight times the spin parameter. (I.e., the k and spin
         !   weighted electron population.)
         initStateFactor = electronPopulation(orderedIndex) / &
               & (kPointWeight(currentKPoint)/real(spin,double))

         ! Determine the array index value of the current final (index j)
         !   spin-kpoint-state as defined by the tempEnergyEigenValues loop
         !   near the beginning of the population subroutine.
         orderedIndex = j + numStates*(spinDirection-1) + &
               & numStates*spin*(currentKPoint-1)

         ! The final state factor is just like the initial state factor
         !   except that here we care about the size of the "hole". Hence,
         !   we compute 1 - the k and spin weighted electron population.
         finStateFactor = 1.0_double - electronPopulation(orderedIndex) / &
               & (kPointWeight(currentKPoint)/real(spin,double))

#ifndef GAMMA
         ! Loop to obtain the wave function times the momentum integral.
         do k = initComponent,finComponent
             valeValeXMom = sum(valeVale(:,i,1) * &
                   & conjWaveMomSum(:,finalStateIndex,k))

            ! Compute the imaginary component of the square of the
            !   valeValeXMom to obtain the transition probabilty element.
            ! Note that the reason it is imaginary is because of the negative
            !   sign included in the getIntgResults subroutine for the
            !   momentum matrix. (See notes in that code.)
            transitionProb(k,initialStateIndex,finalStateIndex) = &
                  & (real(valeValeXMom,double)**2 + aimag(valeValeXMom)**2) * &
                  & (initStateFactor * finStateFactor)
         enddo
#else

         ! Loop to obtain the wave function times the momentum integral.
         do k = initComponent, finComponent
            valeValeXMomGamma = sum(valeValeGamma(:,i,1) * &
                  & conjWaveMomSumGamma(:,finalStateIndex,k))

            ! Compute the real component of the square of the valeValeXMom to
            !   obtain the transition probability.
            transitionProb(k,initialStateIndex,finalStateIndex) = &
                  & valeValeXMomGamma**2 * initStateFactor*finStateFactor
         enddo
#endif
      enddo
   enddo

   ! Determine the weighting effect of this kpoint and include the normalizaion
   !   factor for the gaussian. We will actually be multiplying two Gaussians
   !   together so we need this squared.
   kPointFactor = kPointWeight(currentKPoint)/real(spin,double) / &
         & (sigmaSqrt2Pi)**2

   ! Now compute the exponential alpha factor which is -1/(2 * sigma^2).
   alphaFactor = -1.0_double / (2.0_double * sigmaSIGE**2)

   ! Now we fill up the sigmaEAccumulator.  There are two scenarios in which
   !   this may happen, the xyz all at once (xyzComponents==0) case, and the
   !   each axis (x, y, z) one at a time case (xyzComponents/=0). The only
   !   real difference is the k-loop in the first case and the specific index
   !   access in the second.
   ! The basic approach is to loop over all pairs of states and for each pair
   !   create a Gaussian (of unit area) for each, multiply them together to
   !   get a new Gaussian between them, and then numerically evaluate and
   !   accumulate those on a mesh.
   if (xyzComponents == 0) then

      ! Initialize the counter for the index number of the initial states.
      initialStateIndex = 0

      ! Loop over the set of initial states from first to last. Recall that
      !   this list will include fully occupied states a little below the first
      !   unoccupied state up to the last state with any occupation (as long
      !   as it is within range of the fermi level).
      do i = firstInit, lastInit

         ! Increment the index number for the initial states.
         initialStateIndex = initialStateIndex + 1

         ! Initialize the counter for the index number of the final states.
         finalStateIndex = 0
         do j = firstFin, lastFin

            ! Increment the index number for the final states. Note that we
            !   need to increment the index for every j-loop iteration because
            !   we computed the conjWaveMomSum and conjWaveMomSumGamma for
            !   every final state. We have to increment the counter for every
            !   case, even if we don't use it.  (This note is here because this
            !   was a point of confusion in the past when the code was being
            !   developed.)
            finalStateIndex = finalStateIndex + 1

            ! We don't want double counting or self interactions so we skip the
            !   i>=j cases.
            if (i >= j) cycle

            ! Here we compute the product of two gaussians evaluated on a mesh.
            !   The kPointFactor is the square of (the kpoint weight divided by
            !   the spin (2 or 1 depending on spin polarized or not), divided
            !   by sqrt(2Pi)*sigmaSIGE broadening factor). The exponential has
            !   the alphaFactor of -1/(2*sigmaSIGE^2) times the r^2 value where
            !   the r^2 value is essentially the distance between the current
            !   mesh point and the center of the new Gaussian that is formed
            !   by the product of the Gaussians for each transition state. The
            !   r^2 can be expressed as r1^2 +r2^2 where r1=e-e1 and r2=e-e2
            !   with e=the current mesh point energy and e1=the current initial
            !   state energy value and e2=the current final state energy value.
            do k = 1, 3
               sigmaEAccumulator(:,k) = sigmaEAccumulator(:,k) + &
                     & kPointFactor * exp(alphaFactor * ( &
                     & (energyScale(:) - energyEigenValues(i,currentKPoint,&
                     & spinDirection))**2 + &
                     & (energyScale(:) - energyEigenValues(j,currentKPoint,&
                     & spinDirection))**2)) * transitionProb(k,&
                     & initialStateIndex,finalStateIndex)
            enddo
         enddo
      enddo
   else

      ! Now we create a Gaussin of unit area for each pair of states and
      !   multiply them together to create a new Gaussian (somewhere between
      !   them). These product Gaussians will be evaluated and accumulated on a
      !   numerical mesh.
      ! Initialize the counter for the index number of the initial states.
      initialStateIndex = 0
      do i = firstInit, lastInit

         ! Increment the index number for the initial states.
         initialStateIndex = initialStateIndex + 1

         ! Initialize the counter for the index number of the final states.
         finalStateIndex = 0
         do j = firstFin, lastFin

            ! Increment the index number for the final states. Note that we
            !   need to increment the index for every j-loop iteration because
            !   we computed the conjWaveMomSum and conjWaveMomSumGamma for
            !   every final state. We have to increment the counter for every
            !   case, even if we don't use it.  (This note is here because this
            !   was a point of confusion in the past when the code was being
            !   developed.)
            finalStateIndex = finalStateIndex + 1

            ! We don't want double counting or self interactions so we skip the
            !   i>=j cases.
            if (i >= j) cycle

            ! See notes for the above case.
            sigmaEAccumulator(:,xyzComponents) = &
                  & sigmaEAccumulator(:,xyzComponents) + &
                  & kPointFactor * exp(alphaFactor * ( &
                  & (energyScale(:) - energyEigenValues(i,currentKPoint,&
                  & spinDirection))**2 + &
                  & (energyScale(:) - energyEigenValues(j,currentKPoint,&
                  & spinDirection))**2)) * transitionProb(1,&
                  & initialStateIndex,finalStateIndex)
         enddo
      enddo
   endif

   ! Print the results during the last KPoint iteration for the last component
   !   or set of components.
   if ((currentKPoint == numKPoints) .and. ((xyzComponents == 0) .or. &
         & (xyzComponents == 3))) then

      ! This is a bit tricky and so it will help a lot to know something
      !   about cgs esu units before you start.  The Kubo-Greenwood formula
      !   (KGF) given in this reference:  Greenwood, Proc. Phys. Soc. London,
      !   71,585, (1958) is in cgs esu units.  We did the calculation in a.u.
      !   and the result must be presented in SI.  Arrrg.  Note some deviations
      !   from the formula as presented and as we calculated.  (1)  The
      !   formula used the velocity operator while this calculation used the
      !   momentum operator, so our equations have the electron mass squared in
      !   the denominator to compensate.  (2)  The formula should be considered
      !   as a vector outer product that produces a rank-2 tensor.  Our
      !   computed result is an x, y, z vector that is then averaged when it is
      !   printed out.  That is where the 1/3 factor comes from in our equation.

      ! The equation used is as follows:
      !    2 * Pi * hbar * e^2 / (3 * m^2 * Omega) *
      !    Sum ( |Pnm|^2 * delta(E-En) * delta(E-Em))

      ! In a.u. hbar = 1, e = 1, and m = 1 (They still have units though)

      ! In a.u. we have E=energy, T=time, M=mass, L=length.

      ! The first step is to apply the terms in the coefficient that are not
      !   equal to one.  (Note that we do not divide by three here since we
      !   will perform the averaging later.  We also do not change the units of
      !   the cell volume since that will be accounted for next.)
      conversionFactor = 2.0_double * pi / realCellVolume


      ! The units of the KGF in a.u. are seen as:
      ! SigmaE = (ET EL) / (M^2 L^3) * (M^2 L^2)/(T^2) * 1/(E^2)
      ! SigmaE = 1/T
      ! hBar is ET, charge^2 (Q^2) is ML^3/T^2 = EL, we divide by the volume
      !   explicitly unlike in the KGF.  The tricky part is the charge because
      !   it is not a separate fundamental unit in cgs.  In SI it is in a
      !   separate fundamental unit called Coulombs, but in cgs it is not.  In
      !   cgs the units of Q are sqrt(g cm^3 / s^2).

      ! The second step is then to convert the 1/T in a.u. to 1/sec in cgs.
      !   This can be done by the conversion factor of 2.418884326505x10-17 s
      !   equals one atomic unit of time taken from NIST (auTime).  Note that
      !   the factor of 1e-17 is not included in the constant auTime and it
      !   must be accounted for later.

      ! This will convert the result from (a.u. 1/T) to (cgs emu 1/s).  Note
      !   that it still needs to be multiplied by 1e+17.
      conversionFactor = conversionFactor / auTime

      ! The third step is then to convert the 1/sec in cgs to 1/(ohm m) in
      !   SI.  This can be done with the knowledge that ohm = ET/Q^2 = Js/C^2.
      !   Writing that in cgs esu we have:  (g cm^2 / s^2) s s^2/(g cm^3) = 
      !   s / cm.  Then 1/(ohm cm) = cm / (s cm) = 1/s in cgs units.  This
      !   shows that the two units are equivalent.  Now we just need the
      !   conversion factor between them.
      !   We have from Jackson:  1 / (ohm m) = (1e-16 * c^2) * 1e9 1/s
      !   Therefor: 1/s (cgs emu) = 1/(1e-16 * c^2) * 1e-9 1/(ohm m)
      !   The value of 1e-16 * c^2 is DEFINED as:  8.9875517873681764.

      ! This will convert the result from (cgs emu 1/s) to (SI 1/(ohm m)).
      !   Note that it still needs to be multiplied by an additional 1e-9.
      !   This creates a total exponent multiplication factor of 1e+8 so far.
      conversionFactor = conversionFactor / lightFactor

      ! Finally, we want to express the result in (micro ohm cm)^-1 which is
      !   equal to (1e-8 ohm m)^-1 = 1e8(ohm m)^-1.  This creates a total
      !   exponent multiplication factor of 1 so we don't have to make any more
      !   modifications.

      ! Adjust the sigmaE to have the correct units.
      sigmaEAccumulator(:,:) = sigmaEAccumulator(:,:) * conversionFactor

      ! Now we can compute the final DC conductivity.
      do i = 1, numEnergyPoints
!         dcCond = dcCond + 
      enddo

      do i = 1, numEnergyPoints

         ! Write the averaged total, and the x, y, z components, making sure to
         !   convert the energy scale into eV.
         write (49+spinDirection,fmt="(5e20.8e4)") energyScale(i)*hartree,&
            & sum(sigmaEAccumulator(i,:))/3.0_double,sigmaEAccumulator(i,:)
      enddo

      ! Accumulate the total electronic contribution to thermal conductivity.
      totalSigma = sum(sigmaEAccumulator(:,:)) / 3.0_double
      write (20,*) "The total sigma is: ",totalSigma

      deallocate (sigmaEAccumulator)
      deallocate (energyScale)
   endif

   ! Deallocate unnecessary matrices.
#ifndef GAMMA
   deallocate (conjWaveMomSum)
#else
   deallocate (conjWaveMomSumGamma)
#endif
   deallocate (transitionProb)

end subroutine computeSigmaE


end module O_OptcTransitions

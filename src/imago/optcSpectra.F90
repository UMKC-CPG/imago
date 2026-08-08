!! SPDX-License-Identifier: ECL-2.0
!! Copyright (c) 2026 Paul Rulis

module O_OptcSpectra

   use O_Kinds

   implicit none

   ! The accumulated spectra, held at module scope because accumulating
   !   them and writing them are now separate routines, called in turn
   !   from subroutine optc exactly as subroutine dos calls its own pair
   !   (DESIGN 12.2). This is the pattern O_OptcTransitions already
   !   follows with transitionProb; these were local only by habit.
   real (kind=double), allocatable, dimension (:,:,:) :: optcCond
   real (kind=double), allocatable, &
         & dimension (:,:,:,:,:) :: optcCondPOPTC

   ! Set while the spectra are computed and read again while they are
   !   written, so they travel with the arrays rather than being derived
   !   twice from the input.
   integer :: numEnergyPoints
   real (kind=double) :: conversionFactor
   real (kind=double) :: conversionFactorEps2

contains

! Compute the broadened spectra and leave them in the module arrays
!   above. The Brillouin-zone integration method is selected here, at
!   the top of the work rather than inside it, because the two methods
!   do not share a loop structure: the Gaussian accumulation walks
!   k-points and the tetrahedron accumulation walks tetrahedra, so an
!   internal branch would be a branch around the entire body.
subroutine computeOptcSpectra(doOPTC)

   ! Import necessary data modules.
   use O_Kinds
   use O_Potential,       only: spin
   use O_Lattice,         only: realCellVolume
   use O_KPoints,         only: numKPoints, kPointWeight
   use O_OptcTransitions, only: maxTransEnergy, energyMin, energyScale,&
                                & sumNumPartials, numAccumCompOPTC, &
                                & numAccumCompPOPTC
   use O_Input,           only: sigmaOPTC, deltaOPTC, sigmaPACS, deltaPACS,&
                                & detailCodePOPTC
   use O_Constants,       only: dim3, pi, auTime, eCharge, hPlanck, hartree
   use O_KPoints,         only: kPointIntgCode, symmetrizeLATPartials

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define passed parameters.
   integer, intent(in) :: doOPTC

   ! Define local variables
   integer :: i ! Loop index variables
   real (kind=double) :: sigma
   real (kind=double) :: energyDelta
   real (kind=double) :: sigmaSqrtPi
   real (kind=double), allocatable, dimension (:)     :: kPointFactor
   real (kind=double) :: latFactor ! Tetrahedron pathway scale factor.


   ! Initialize variables.
   if (doOPTC == 1) then ! Standard optical properties calculation.
      sigma           = sigmaOPTC
      energyDelta     = deltaOPTC
      energyMin       = deltaOPTC ! Start as close to 0 as possible.
      numEnergyPoints = int(maxTransEnergy / energyDelta) + 1
   elseif (doOPTC == 2) then ! PACS calculation
      sigma           = sigmaPACS
      energyDelta     = deltaPACS
      ! The energyMin was already determined for PACS calculations.
      numEnergyPoints = int((maxTransEnergy - energyMin) / energyDelta) + 1
   endif ! Sigma(E) calculations never call this subroutine.

   ! Initialize local conversion parameters

   ! Written to agree with Cohen and Chelikowsky, "Electronic Structure and
   !   Optical Properties of Semiconductors" section 4.1 equation 4.10.  Note
   !   that we first compute sigma = e2(w) * w / (4*pi)  Note that equation
   !   4.10 comes from Ehrenreich and Cohen, Phys. Rev., 115, 786, (1959).
   ! w=omega

   ! The equation that we use is:
   !   Sigma(omega) = 2 * pi * e^2 / (3 * m^2 * omega) * (1/(2*Pi)^3) *
   !      sum(ij) [Int(BZ) [delta(Ej - Ei - hbar*omega) * |Mij(k)|^2 d3k]]
   !   omega = angular frequency
   !   Mij = Momentum matrix element
   !   Note that the integration over BZ will produce a factor of (2*Pi)^3 / 
   !      Omega.  Where Omega equals the cell volume.
   ! 

   ! In a.u. hbar=1, e=1, and m=1 (Though they still have units)

   ! In a.u. we have E=energy, T=time, M=mass, L=length.

   ! The first step is to apply the terms in the coefficient that are not
   !   equal to one.  (Note that we do not divide by three here since we
   !   will perform the averaging later.  We also do not change the units of
   !   the cell volume since that will be accounted for later.)
   conversionFactor = (2.0_double*pi) / realCellVolume

   ! The next step is to note that our calculation uses energy instead of
   !   angular frequency (omega), and that our energy is in eV.  We must
   !   convert it back from eV to a.u., then we use the fact that hbar in a.u.
   !   equals 1 to relate energy to frequency.  Since we will divide by eV we
   !   convert it to a.u. by multiplying by the hartree factor.
!   conversionFactor = conversionFactor * hartree

   ! Along a similar vein, the delta term (delta(Ej - Ei - hbar*omega)) has
   !   units of inverse energy, but this too is in eV and must be converted
   !   back to a.u. by another multiplication by the hartree factor.
!   conversionFactor = conversionFactor * hartree

   ! The next step is to do dimensional analysis and understand the unit
   !   conversion.  The equation has the following units in a.u.:
   
   ! Sigma = EL / (M^2 * T^-1) * 1/L^3 * 1/E * (ML/T)^2
   ! Sigma = 1/T
   ! Charge^2 (Q^2) is ML^3/T^2 = EL.  (In cgs charge is sqrt(g cm^3 / s^2).)
   ! Convert 1/T in a.u. to 1/s in cgs by applying the conversion factor of
   !   2.418884326505x10-17 s = 1 a.u. of time (taken from NIST).  Note that
   !   the factor of 1e-17 is not included in the constant auTime and it must
   !   be accounted for somewhere down the road.
   conversionFactor = conversionFactor/(auTime)

   ! The last step for the conductivity is to adjust the result to have the
   !   appropriate order of magnitude.  The result should be in units of
   !   1e15 * 1/sec.  The 1/1e-17 from auTime gives us 1e17 1/s.  To put the
   !   result in the right units we must multiply the answer by 100.
   conversionFactor = conversionFactor * 100.0_double

   ! We also want to produce a result in terms of the unitless epsilon 2.  To
   !   do this we will multiply the conductivity by 4pi and divide by the
   !   frequency.  This means that we have to multiply by hbar in eV to
   !   have the proper units since we are giving the value in eV and it
   !   needs to be a frequency in 1/s.  The 2pi from the hbar will cancel
   !   with the 4pi to leave just 2.  We also divide by eCharge to put the
   !   value of hPlanck in eV s instead of Js.
   conversionFactorEps2 = 2.0_double * hPlanck / eCharge

   ! Again we must adjust for the units so we will multiply the result by
   !   1x10^-15 because of the difference between 1d-34 and 1d-19 for hPlanck
   !   and eCharge.  I'm not yet sure why we don't apply this factor.  CHECK!
!   conversionFactorEps2 = conversionFactorEps2 * 1.0d-15
   conversionFactorEps2 = conversionFactorEps2

   ! Used for normalization of the convoluted Gaussian
   sigmaSqrtPi = sqrt(pi) * sigma

   ! Allocate space for local arrays and matrices.
   allocate (kPointFactor (numKPoints))

   ! Begin setting up and initializing the optical conductivity parameters.

   ! Allocate space to hold the energy scale.
   allocate (energyScale (numEnergyPoints))

   ! Allocate space to hold the appropriate optical conductivity and then
   !   initialize.
   ! How many spectra each array carries is set by the direction code
   !   rather than fixed at three (PSEUDOCODE 21.3): one for the
   !   direction average alone, three for x, y and z, six for the
   !   symmetric tensor. The leading isotropic column is NOT among them;
   !   it is formed from the diagonal entries when the file is written.
   if (detailCodePOPTC == 0) then
      allocate (optcCond(numAccumCompOPTC,numEnergyPoints,spin))
      optcCond(:,:,:) = 0.0_double
   else
      allocate (optcCond(numAccumCompOPTC,numEnergyPoints,spin))
      allocate (optcCondPOPTC(sumNumPartials,sumNumPartials, &
            & numAccumCompPOPTC,numEnergyPoints,spin))
      optcCond(:,:,:) = 0.0_double
      optcCondPOPTC(:,:,:,:,:) = 0.0_double
   endif

   ! Assign values to the energy range
   do i = 1, numEnergyPoints
      energyScale(i) = energyMin + energyDelta * (i-1)
   enddo

   ! Fill in factor for broadening based on kpoint weight.  The 0.5 must be
   !   included because we have already accounted for the fact that each state
   !   (in the spin non-polarized case) has two electrons and when we consider
   !   the kpoint weighting we don't want to re-multiply by 2.0.
   !   Recall that the sum(kPointWeight(:)) == 2.  In the spin polarized case
   !   the division by "spin" will divide by two because now we say that there
   !   is only one electron per state.  Also note that we must divide by
   !   hartree to have the right units for sigmaSqrtPi.
   kPointFactor(:) = kPointWeight(:) * 0.5_double / sigmaSqrtPi / hartree / &
         & real(spin,double)

   ! The tetrahedron pathway's scale factor, which replaces
   !   kPointFactor rather than adjusting it. Derived term by term in
   !   DESIGN 12.4 against the expression just above:
   !
   !   - kPointWeight summed over the irreducible points becomes tetraVol
   !     summed over tetrahedra. The first sums to 2 and the second to 1,
   !     so sum(kPointWeight) restores the scale.
   !   - 1/sigmaSqrtPi is DROPPED. It normalizes a Gaussian standing in
   !     for a delta function, and the corner density weight is that
   !     delta function evaluated exactly, already carrying units of
   !     inverse energy.
   !   - the 0.5, the 1/hartree and the 1/spin belong to the quantity
   !     rather than to the integration, so all three survive. The 0.5 in
   !     particular is not geometric: the transition probabilities
   !     already account for two electrons per state in the
   !     spin-unpolarized case, and it stops the weights from counting
   !     them a second time. Dropping it would halve every spectrum.
   latFactor = sum(kPointWeight(:)) * 0.5_double / hartree &
         & / real(spin,double)

   ! Accumulate the spectra by whichever method was requested. Both
   !   leave the same module arrays filled, so nothing downstream learns
   !   which one ran.
   if (kPointIntgCode == 1) then
      call accumulateOptcCond_LAT (latFactor)
      if (detailCodePOPTC /= 0) then
         call accumulateOptcCondPOPTC_LAT (latFactor)
      endif
   else
      call accumulateOptcCond (kPointFactor,sigma)
      if (detailCodePOPTC /= 0) then
         call accumulateOptcCondPOPTC (kPointFactor,sigma)
      endif
   endif

   ! Average the atom-resolved spectra over the point group (DESIGN 1.7).
   !   Only the tetrahedron pathway needs it: the Gaussian one spreads each
   !   irreducible point's contribution evenly across its star and so is
   !   already symmetric. Only the atom grouped codes need it either, since
   !   every operation carries an atom onto an atom of the same type and a
   !   type level sum therefore maps onto itself.
   !
   ! This must run BEFORE printOptcSpectra, and it does so while partialPerm
   !   is still alive -- subroutine optc releases the decomposition index
   !   after the spectra are written, precisely so that this window exists.
   if ((kPointIntgCode == 1) .and. (symmetrizeLATPartials == 1) &
         & .and. (detailCodePOPTC >= 3)) then
      call symmetrizeOptcPOPTC_LAT
   endif

   deallocate (kPointFactor)

end subroutine computeOptcSpectra


! Write the spectra that computeOptcSpectra left in the module arrays,
!   and release them. Kept apart from the computation so that each
!   routine can be named for the one job it does, and so that the
!   integration method is invisible here: by this point the two
!   pathways have produced the same arrays and differ in nothing the
!   printer can see.
subroutine printOptcSpectra(doOPTC)

   ! Import necessary data modules.
   use O_Kinds
   use O_OptcTransitions, only: energyScale
   use O_Input,           only: detailCodePOPTC

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define passed parameters.
   integer, intent(in) :: doOPTC

   if (detailCodePOPTC == 0) then ! Regular total optical properties.

      if (doOPTC == 2) then ! Doing PACS calculation.
         call printSpectrum(0,numEnergyPoints,optcCond,conversionFactor)
      else ! Do optical conductivity followed by epsilon 2.
         call printSpectrum(1,numEnergyPoints,optcCond,conversionFactor)
         call printSpectrum(2,numEnergyPoints,optcCond,conversionFactorEps2)
      endif

   else ! Partial optical properties

      if (doOPTC == 2) then ! Doing PACS calculation.
         call printSpectrum(0,numEnergyPoints,optcCond,conversionFactor)
         call printSpectrumPOPTC(0,numEnergyPoints,optcCondPOPTC,&
               & conversionFactor)
      else ! Do optical conductivity followed by epsilon 2.
         call printSpectrum(1,numEnergyPoints,optcCond,conversionFactor)
         call printSpectrumPOPTC(1,numEnergyPoints,optcCondPOPTC,&
               & conversionFactor)
         call printSpectrum(2,numEnergyPoints,optcCond,conversionFactorEps2)
         call printSpectrumPOPTC(2,numEnergyPoints,optcCondPOPTC,&
               & conversionFactorEps2)
      endif
   endif

   ! Deallocate arrays.
   deallocate (energyScale)
   deallocate (optcCond)
   if (detailCodePOPTC /= 0) then
      deallocate (optcCondPOPTC)
   endif

end subroutine printOptcSpectra


! Accumulate the broadened total spectrum by Gaussian smearing over the
!   irreducible k-points. Named for what it does: it adds into an array
!   rather than returning anything, and the unsuffixed name marks it as
!   the Gaussian member of the pair, following the convention DESIGN 1.5
!   sets with electronPopulation and electronPopulation_LAT. That the
!   plain name means Gaussian is a convention rather than something the
!   name states.
subroutine accumulateOptcCond (kPointFactor, sigma)

   ! Import the necessary data modules.
   use O_Kinds
   use O_Potential,       only: spin
   use O_KPoints,         only: numKPoints
   use O_OptcTransitions, only: energyScale, energyDiff, transCounter, &
         & transitionProb, transitionMoment, pairOccupancy, &
         & numStoredCompOPTC, numAccumCompOPTC

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), dimension (:)     :: kPointFactor
   real (kind=double) :: sigma

   ! Define local variables
   real (kind=double) :: broadenEnergyDiff
   real (kind=double) :: expAlpha
   real (kind=double) :: broadenWeight ! Everything but the strength.
   integer :: h,i,j,k
   integer :: c,d,e ! Cartesian component indices.
   real (kind=double), allocatable, dimension (:,:,:,:,:) :: starRotAvg
   real (kind=double), dimension (numAccumCompOPTC) :: strength
         !   This transition's contribution before broadening.
   real (kind=double), dimension (3,3) :: momentProduct
         !   Re[ M_d conjg(M_e) ], the only combination of the stored
         !   element any direction code needs.

   ! At direction codes 1 and 2 the star-averaged rotation products have
   !   to be in hand before the loop begins. See the routine's own
   !   commentary for why an average and not a sum.
   if (numStoredCompOPTC > 1) then
      allocate (starRotAvg (3,3,3,3,numKPoints))
      call buildStarRotationAverage (starRotAvg)
   endif

   do h = 1, spin
      do i = 1, numKPoints
         do j = 1, transCounter(i,h)

            ! Form this transition's strength once per pair rather than
            !   once per energy point. It does not depend on the energy
            !   point at all, and the energy loop below runs to
            !   thousands of points, so building it inside would repeat
            !   the whole rotation contraction for every one of them.
            if (numStoredCompOPTC == 1) then

               ! Direction code 0. The producer already collapsed to the
               !   rotation-invariant direction sum, so there is nothing
               !   to rotate and nothing to combine.
               strength(1) = transitionProb(1,j,i,h)
            else

               ! Codes 1 and 2. Build the outer product of the matrix
               !   element with its own conjugate. Only the real part
               !   survives into any spectrum: the imaginary part is
               !   antisymmetric in d and e while the rotation product
               !   contracted against it is symmetric, so it cancels.
               do d = 1, 3
                  do e = 1, 3
#ifndef GAMMA
                     momentProduct(d,e) = &
                           & real(transitionMoment(d,j,i,h) &
                           & * conjg(transitionMoment(e,j,i,h)),double)
#else
                     momentProduct(d,e) = transitionMoment(d,j,i,h) &
                           & * transitionMoment(e,j,i,h)
#endif
                  enddo
               enddo

               ! Contract with the star average. For the diagonal
               !   entries this is the rotated squared magnitude
               !   averaged over the star; for the off-diagonal entries
               !   of code 2 it is the corresponding averaged product.
               call contractStarRotation (starRotAvg(:,:,:,:,i), &
                     & momentProduct,numAccumCompOPTC,strength)

               ! The occupancy weight, held out of the element so that
               !   it never had to be square rooted, rejoins here.
               strength(:) = strength(:) * pairOccupancy(j,i,h)
            endif

            do k = 1, numEnergyPoints

               ! Determine the energy difference between the current transition
               !   energy and the current energy scale point.
               broadenEnergyDiff = energyDiff(j,i,h) - energyScale(k)

               ! Determine the exponential alpha for the broadening factor.
               expAlpha = broadenEnergyDiff * broadenEnergyDiff / &
                     & (sigma * sigma)

               ! If the exponential alpha is too large, then we don't have to
               !   complete the broadening for this set because it will not
               !   have a significant effect.
               if (expAlpha < 50.0_double) then
                  broadenWeight = exp(-expAlpha) * kPointFactor(i)
                  optcCond(:,k,h) = optcCond(:,k,h) &
                        & + strength(:) * broadenWeight
               endif
            enddo
         enddo
      enddo
   enddo

   if (allocated(starRotAvg)) then
      deallocate (starRotAvg)
   endif

end subroutine accumulateOptcCond


! Average the product of two rotation entries over the star of each
!   irreducible k-point (DESIGN 13.5, PSEUDOCODE 21.5).
!
! Why this exists at all. The Gaussian pathway has NO star loop: it
!   accumulates over irreducible points and carries the star
!   multiplicity inside kPointWeight. So there is nowhere for a
!   per-star-member rotation to live, and adding such a loop would
!   multiply the cost of the accumulation by the reduction factor --
!   four to forty-eight. What the accumulation actually needs from the
!   star is only the AVERAGE of R(c,d)*R(c2,e) over its members, which
!   depends on the star alone: not on the band, not on the energy, not
!   on the matrix element. So it is computed once here and reused.
!
! It must be the AVERAGE and never the sum. kPointWeight already
!   carries the star multiplicity, so summing here would apply the
!   symmetry twice -- the same double-counting trap PSEUDOCODE 19.5
!   records for the tetrahedron pathway. On an unreduced mesh every
!   star has one member, the average is the identity's own product, and
!   the accumulation reduces exactly to what it was before this change.
!
! The four index form generalizes PSEUDOCODE 21.5, which writes three
!   indices because it describes direction code 1. Code 1 needs only
!   the entries with c2 equal to c, and reads them from this array
!   unchanged; code 2 needs the off diagonal c2 as well.
subroutine buildStarRotationAverage (starRotAvg)

   ! Import the necessary data modules.
   use O_Kinds
   use O_KPoints, only: numKPoints, numFullMeshKP, fullKPToIBZKPMap, &
         & fullKPToIBZOpMap, xyzRealPointOps

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), dimension (3,3,3,3,numKPoints) :: starRotAvg

   ! Define local variables.
   integer :: kFull, kIBZ, c, c2, d, e
   integer, allocatable, dimension (:) :: starSize
   real (kind=double), dimension (3,3) :: rotation

   allocate (starSize (numKPoints))
   starSize(:) = 0
   starRotAvg(:,:,:,:,:) = 0.0_double

   ! Without a folded mesh there is no star to average over. This
   !   happens for an explicitly listed k-point set, where each point
   !   stands for itself; the identity below then leaves the
   !   accumulation exactly as it was.
   if ((.not. allocated(fullKPToIBZKPMap)) .or. &
         & (.not. allocated(xyzRealPointOps))) then
      do kIBZ = 1, numKPoints
         do c = 1, 3
            starRotAvg(c,c,c,c,kIBZ) = 1.0_double
         enddo
      enddo
      deallocate (starSize)
      return
   endif

   do kFull = 1, numFullMeshKP
      kIBZ = fullKPToIBZKPMap(kFull)
      rotation(:,:) = xyzRealPointOps(:,:,fullKPToIBZOpMap(kFull))
      starSize(kIBZ) = starSize(kIBZ) + 1

      do c = 1, 3
         do c2 = 1, 3
            do d = 1, 3
               do e = 1, 3
                  starRotAvg(c,c2,d,e,kIBZ) = &
                        & starRotAvg(c,c2,d,e,kIBZ) &
                        & + rotation(c,d) * rotation(c2,e)
               enddo
            enddo
         enddo
      enddo
   enddo

   do kIBZ = 1, numKPoints
      if (starSize(kIBZ) > 0) then
         starRotAvg(:,:,:,:,kIBZ) = starRotAvg(:,:,:,:,kIBZ) &
               & / real(starSize(kIBZ),double)
      endif
   enddo

   deallocate (starSize)

end subroutine buildStarRotationAverage


! Contract one k-point's star-averaged rotation products against one
!   transition's moment product, producing the spectra entries this
!   direction code asks for.
!
! The tensor entry order is the one PSEUDOCODE 21.7 declares and the
!   printer relies on: xx, yy, zz, then xy, xz, yz. Changing it here
!   without changing the header there would mislabel every column, in a
!   way no sum rule could detect.
subroutine contractStarRotation (starRotAvg,momentProduct,numWanted, &
      & strength)

   ! Import the necessary data modules.
   use O_Kinds

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), dimension (3,3,3,3) :: starRotAvg
   real (kind=double), dimension (3,3) :: momentProduct
   integer :: numWanted
   real (kind=double), dimension (numWanted) :: strength

   ! Define local variables.
   integer :: slot, d, e
   integer, dimension (6), parameter :: firstOfPair = &
         & (/ 1, 2, 3, 1, 1, 2 /)
   integer, dimension (6), parameter :: secondOfPair = &
         & (/ 1, 2, 3, 2, 3, 3 /)

   strength(:) = 0.0_double

   do slot = 1, numWanted
      do d = 1, 3
         do e = 1, 3
            strength(slot) = strength(slot) &
                  & + starRotAvg(firstOfPair(slot), &
                  & secondOfPair(slot),d,e) * momentProduct(d,e)
         enddo
      enddo
   enddo

end subroutine contractStarRotation


! The decomposed counterpart, identical in structure and looping over
!   the pair matrix instead of the single total.
subroutine accumulateOptcCondPOPTC (kPointFactor, sigma)

   ! Import the necessary data modules.
   use O_Kinds
   use O_Potential,       only: spin
   use O_KPoints,         only: numKPoints
   use O_OptcTransitions, only: energyScale, energyDiff, transCounter, &
         & transitionProbPOPTC

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   real (kind=double), dimension (:)     :: kPointFactor
   real (kind=double) :: sigma

   ! Define local variables
   real (kind=double) :: broadenEnergyDiff
   real (kind=double) :: expAlpha
   integer :: h,i,j,k


   do h = 1, spin
      do i = 1, numKPoints
         do j = 1, transCounter(i,h)
            do k = 1, numEnergyPoints

               ! Determine the energy difference between the current transition
               !   energy and the current energy scale point.
               broadenEnergyDiff = energyDiff(j,i,h) - energyScale(k)

               ! Determine the exponential alpha for the broadening factor.
               expAlpha = broadenEnergyDiff * broadenEnergyDiff / &
                     & (sigma * sigma)

               ! If the exponential alpha is too large, then we don't have to
               !   complete the broadening for this set because it will not
               !   have a significant effect.
               if (expAlpha < 50.0_double) then
                  optcCondPOPTC(:,:,:,k,h) = optcCondPOPTC(:,:,:,k,h) &
                        & + transitionProbPOPTC(:,:,:,j,i,h) &
                        & * exp(-expAlpha) * kPointFactor(i)
               endif
            enddo
         enddo
      enddo
   enddo

end subroutine accumulateOptcCondPOPTC


! Accumulate the broadened total spectrum by tetrahedron integration.
!   PSEUDOCODE 19.3.
!
! The loop structure inverts relative to the Gaussian routine above:
!   outer over band pairs and tetrahedra rather than over k-points and
!   the transitions found at each. What is integrated is a JOINT density
!   of states -- the surface where the difference of two bands equals
!   the output energy -- so the corner values fed to the Bloechl weights
!   are differences of eigenvalues rather than eigenvalues, and the
!   weight routine is reused untouched because it is a function of four
!   corner values and does not care what produced them.
subroutine accumulateOptcCond_LAT (latFactor)

   ! Import the necessary data modules.
   use O_Kinds
   use O_Potential,       only: spin
   use O_MathSubs,        only: bloechlCornerDOSWt
   use O_SecularEquation, only: energyEigenValues
   use O_KPoints,         only: numTetrahedra, tetraVol, tetrahedra, &
         & fullKPToIBZKPMap, fullKPToIBZOpMap, xyzRealPointOps
   use O_OptcTransitions, only: energyScale, transProbBanded, &
         & bandedInitLo, bandedInitHi, bandedFinLo, bandedFinHi, &
         & pairIsWanted, transMomentBanded, bandedOccupancy, &
         & numStoredCompOPTC, numAccumCompOPTC

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   real (kind=double) :: latFactor

   ! Define local variables.
   integer :: h,i,j    ! Spin, initial band, final band.
   integer :: t        ! Tetrahedron.
   integer :: c,iE     ! Corner, energy point.
   integer :: minIdx, tempInt, orig
   integer :: firstPoint, lastPoint ! Energy range this tetrahedron spans.
   integer, dimension (4) :: sortOrder ! Sorted corner to original corner.
   integer, dimension (4) :: cornerKP  ! IBZ k-point of each corner.
   integer, dimension (4) :: cornerOp  ! Operation carrying IBZ to corner.
   real (kind=double), dimension (4) :: epsDiff
   real (kind=double), dimension (4) :: cornerDOSWt
   real (kind=double) :: tempVal
   real (kind=double) :: energyStep
   real (kind=double) :: energy
   integer :: slot, d, e ! Spectrum slot and Cartesian component indices.
   real (kind=double), dimension (numAccumCompOPTC) :: cornerStrength
         !   One corner's contribution before the tetrahedron weight.
   ! The tensor entry order declared by PSEUDOCODE 21.7 and relied on by
   !   the printed header: xx, yy, zz, xy, xz, yz.
   integer, dimension (6), parameter :: firstOfPair = &
         & (/ 1, 2, 3, 1, 1, 2 /)
   integer, dimension (6), parameter :: secondOfPair = &
         & (/ 1, 2, 3, 2, 3, 3 /)
#ifndef GAMMA
   complex (kind=double), dimension (3) :: rotatedMoment
#else
   real (kind=double), dimension (3) :: rotatedMoment
#endif

   ! The grid is uniform, so its step can be recovered from it rather
   !   than passed in. Taking it from the scale itself means the bounds
   !   computed below cannot disagree with the points they index.
   if (size(energyScale) > 1) then
      energyStep = energyScale(2) - energyScale(1)
   else
      energyStep = 1.0_double
   endif

   do h = 1, spin
      do i = bandedInitLo, bandedInitHi
         do j = bandedFinLo, bandedFinHi

            ! Pairs no tetrahedron can want were dropped when the store
            !   was built and hold nothing but zeros.
            if (.not. pairIsWanted(i,j,h)) cycle

            do t = 1, numTetrahedra

               ! Gather the four corner values of the band-energy
               !   DIFFERENCE. Both eigenvalues are read at the same IBZ
               !   representative, and e(Rk) = e(k), so no permutation
               !   enters here.
               do c = 1, 4
                  cornerKP(c) = fullKPToIBZKPMap(tetrahedra(c,t))
                  cornerOp(c) = fullKPToIBZOpMap(tetrahedra(c,t))
                  epsDiff(c) = energyEigenValues(j,cornerKP(c),h) &
                        & - energyEigenValues(i,cornerKP(c),h)
                  sortOrder(c) = c
               enddo

               ! Sort ascending, tracking the permutation so that a
               !   sorted corner can be mapped back to the k-point it
               !   came from. Same selection sort as integratePDOS_LAT.
               do c = 1, 3
                  minIdx = c
                  do iE = c + 1, 4
                     if (epsDiff(iE) < epsDiff(minIdx)) then
                        minIdx = iE
                     endif
                  enddo
                  if (minIdx /= c) then
                     tempVal = epsDiff(c)
                     epsDiff(c) = epsDiff(minIdx)
                     epsDiff(minIdx) = tempVal
                     tempInt = sortOrder(c)
                     sortOrder(c) = sortOrder(minIdx)
                     sortOrder(minIdx) = tempInt
                  endif
               enddo

               ! Bound the energy loop instead of sweeping the whole
               !   grid and testing inside it. This is not an
               !   optimization to leave for later: a tetrahedron spans
               !   only the few energy points between its lowest and
               !   highest corner difference, there are six tetrahedra
               !   per full-mesh k-point against one pass per
               !   irreducible k-point on the Gaussian side, and the
               !   grid is thousands of points long. Sweeping it all
               !   would multiply work that is already the largest term
               !   by the ratio of the two (PSEUDOCODE 19.2.2).
               !
               ! The range tests are KEPT inside the loop even so. They
               ! cost almost nothing and they mean the bound arithmetic
               ! can be off by one at either end without changing the
               ! answer, only the speed.
               firstPoint = int((epsDiff(1) - energyScale(1)) &
                     & / energyStep) + 1
               lastPoint  = int((epsDiff(4) - energyScale(1)) &
                     & / energyStep) + 2
               if (firstPoint < 1) firstPoint = 1
               if (lastPoint > size(energyScale)) then
                  lastPoint = size(energyScale)
               endif

               do iE = firstPoint, lastPoint
                  energy = energyScale(iE)
                  if (energy < epsDiff(1)) cycle
                  if (energy >= epsDiff(4)) cycle

                  call bloechlCornerDOSWt (energy,epsDiff,cornerDOSWt)

                  ! Each sorted corner carries its weight back to the
                  !   corner it belongs to, whose k-point supplies the
                  !   matrix element. Pairing a sorted weight with an
                  !   unsorted corner is the mistake that yields a
                  !   plausible spectrum rather than a broken one.
                  do c = 1, 4
                     if (abs(cornerDOSWt(c)) < 1.0d-30) cycle
                     orig = sortOrder(c)

                     if (numStoredCompOPTC == 1) then

                        ! Direction code 0. The stored number is the
                        !   rotation-invariant direction sum, so this
                        !   corner needs no operation applied to it.
                        cornerStrength(1) = &
                              & transProbBanded(1,cornerKP(orig),i,j,h)
                     else

                        ! Codes 1 and 2. This corner is a full-mesh
                        !   k-point reached from its irreducible
                        !   representative by cornerOp, and the momentum
                        !   operator is a vector, so the element that
                        !   belongs at this corner is the stored one
                        !   ROTATED by that operation. Rotate first and
                        !   square afterwards; doing it the other way
                        !   round is precisely the defect this change
                        !   exists to remove, because squaring destroys
                        !   the cross terms a rotated square is built
                        !   from.
                        rotatedMoment(:) = matmul( &
                              & xyzRealPointOps(:,:,cornerOp(orig)), &
                              & transMomentBanded(:,cornerKP(orig), &
                              & i,j,h))

                        do slot = 1, numAccumCompOPTC
                           d = firstOfPair(slot)
                           e = secondOfPair(slot)
#ifndef GAMMA
                           cornerStrength(slot) = &
                                 & real(rotatedMoment(d) &
                                 & * conjg(rotatedMoment(e)),double)
#else
                           cornerStrength(slot) = rotatedMoment(d) &
                                 & * rotatedMoment(e)
#endif
                        enddo

                        ! The occupancy weight, kept out of the element
                        !   so it never had to be square rooted, rejoins
                        !   here. It belongs to the irreducible
                        !   representative because occupancy is a
                        !   property of the band and the k-point, both
                        !   of which the operation leaves unchanged.
                        cornerStrength(:) = cornerStrength(:) &
                              & * bandedOccupancy(cornerKP(orig),i,j,h)
                     endif

                     optcCond(:,iE,h) = optcCond(:,iE,h) &
                           & + cornerDOSWt(c) * tetraVol * latFactor &
                           & * cornerStrength(:)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo

end subroutine accumulateOptcCond_LAT


! The decomposed counterpart. Identical in structure, with the pair
!   matrix carried through and BOTH of its partial indices permuted at
!   the corner. This is where the IBZ unfolding happens on this pathway:
!   the star average of PSEUDOCODE 7a does not run here, because corner
!   assembly visits full-mesh points directly and applying both would
!   count the symmetry twice (PSEUDOCODE 19.5).
subroutine accumulateOptcCondPOPTC_LAT (latFactor)

   ! Import the necessary data modules.
   use O_Kinds
   use O_Potential,       only: spin
   use O_MathSubs,        only: bloechlCornerDOSWt
   use O_SecularEquation, only: energyEigenValues
   use O_Input,           only: detailCodePOPTC
   use O_KPoints,         only: numTetrahedra, tetraVol, tetrahedra, &
         & fullKPToIBZKPMap, fullKPToIBZOpMap, xyzRealPointOps
   use O_OptcTransitions, only: energyScale, transProbPOPTCBanded, &
         & bandedInitLo, bandedInitHi, bandedFinLo, bandedFinHi, &
         & pairIsWanted, sumNumPartials, partialPerm, &
         & transMomentBanded, numStoredCompPOPTC, numAccumCompPOPTC

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   real (kind=double) :: latFactor

   ! Define local variables.
   integer :: h,i,j
   integer :: t
   integer :: c,iE
   integer :: minIdx, tempInt, orig, opR
   integer :: firstPoint, lastPoint
   integer :: a,b          ! Initial and final partial of the pair.
   integer :: aRot, bRot   ! Their images under the corner's operation.
   integer :: slot, d      ! Spectrum slot and Cartesian component.
   integer, dimension (4) :: sortOrder
   integer, dimension (4) :: cornerKP
   integer, dimension (4) :: cornerOp ! Operation reaching each corner.
   real (kind=double), dimension (4) :: epsDiff
   real (kind=double), dimension (4) :: cornerDOSWt
   real (kind=double) :: tempVal
   real (kind=double) :: energyStep
   real (kind=double) :: energy
   real (kind=double) :: weight ! The full per-corner scale factor.
   ! The tensor entry order of PSEUDOCODE 21.7, named apart from the
   !   total spectrum's copy because the two widths are set by different
   !   direction codes and must not be made to share a variable.
   integer, dimension (6), parameter :: firstOfPairP = &
         & (/ 1, 2, 3, 1, 1, 2 /)
   integer, dimension (6), parameter :: secondOfPairP = &
         & (/ 1, 2, 3, 2, 3, 3 /)
#ifndef GAMMA
   complex (kind=double), dimension (3) :: rotatedTotal
   complex (kind=double), dimension (6,3) :: crossFactor
#else
   real (kind=double), dimension (3) :: rotatedTotal
   real (kind=double), dimension (6,3) :: crossFactor
#endif

   if (size(energyScale) > 1) then
      energyStep = energyScale(2) - energyScale(1)
   else
      energyStep = 1.0_double
   endif

   do h = 1, spin
      do i = bandedInitLo, bandedInitHi
         do j = bandedFinLo, bandedFinHi

            if (.not. pairIsWanted(i,j,h)) cycle

            do t = 1, numTetrahedra

               do c = 1, 4
                  cornerKP(c) = fullKPToIBZKPMap(tetrahedra(c,t))
                  cornerOp(c) = fullKPToIBZOpMap(tetrahedra(c,t))
                  epsDiff(c) = energyEigenValues(j,cornerKP(c),h) &
                        & - energyEigenValues(i,cornerKP(c),h)
                  sortOrder(c) = c
               enddo

               do c = 1, 3
                  minIdx = c
                  do iE = c + 1, 4
                     if (epsDiff(iE) < epsDiff(minIdx)) then
                        minIdx = iE
                     endif
                  enddo
                  if (minIdx /= c) then
                     tempVal = epsDiff(c)
                     epsDiff(c) = epsDiff(minIdx)
                     epsDiff(minIdx) = tempVal
                     tempInt = sortOrder(c)
                     sortOrder(c) = sortOrder(minIdx)
                     sortOrder(minIdx) = tempInt
                  endif
               enddo

               firstPoint = int((epsDiff(1) - energyScale(1)) &
                     & / energyStep) + 1
               lastPoint  = int((epsDiff(4) - energyScale(1)) &
                     & / energyStep) + 2
               if (firstPoint < 1) firstPoint = 1
               if (lastPoint > size(energyScale)) then
                  lastPoint = size(energyScale)
               endif

               do iE = firstPoint, lastPoint
                  energy = energyScale(iE)
                  if (energy < epsDiff(1)) cycle
                  if (energy >= epsDiff(4)) cycle

                  call bloechlCornerDOSWt (energy,epsDiff,cornerDOSWt)

                  do c = 1, 4
                     if (abs(cornerDOSWt(c)) < 1.0d-30) cycle
                     orig = sortOrder(c)
                     opR  = cornerOp(orig)
                     weight = cornerDOSWt(c) * tetraVol * latFactor

                     ! At partial direction code 1 or 2 the cross term
                     !   has to be built HERE, because both of its
                     !   factors are vectors that this corner's
                     !   operation rotates. Everything that does not
                     !   depend on the partial pair is formed once per
                     !   corner, outside the pair loops below:
                     !
                     !     rotatedTotal(c) = sum_e R(c,e) S_e
                     !     crossFactor(slot,d) = R(c1,d)
                     !                           * conjg(rotatedTotal(c2))
                     !
                     !   after which one partial pair's entry for a slot
                     !   is just sum_d Re[ M_d(a,b) * crossFactor ].
                     !   Hoisting it matters: the pair loops run over
                     !   the partial count SQUARED, which DESIGN 11.4
                     !   names as the cost driver of the whole method.
                     if (numStoredCompPOPTC > 1) then
                        do d = 1, 3
                           rotatedTotal(d) = &
                                 & sum(xyzRealPointOps(d,:,opR) &
                                 & * transMomentBanded(:, &
                                 & cornerKP(orig),i,j,h))
                        enddo
                        do slot = 1, numAccumCompPOPTC
                           do d = 1, 3
                              crossFactor(slot,d) = &
                                    & xyzRealPointOps( &
                                    & firstOfPairP(slot),d,opR) &
#ifndef GAMMA
                                    & * conjg(rotatedTotal( &
                                    & secondOfPairP(slot)))
#else
                                    & * rotatedTotal( &
                                    & secondOfPairP(slot))
#endif
                           enddo
                        enddo
                     endif

                     ! The type grouped codes need no permutation: every
                     !   operation carries an atom onto an atom of the
                     !   same type, so a type level sum maps onto itself
                     !   (DESIGN 2.5). partialPerm is not even built for
                     !   them, so it must not be indexed here.
                     if ((detailCodePOPTC >= 3) .and. &
                           & (allocated(partialPerm))) then

                        ! b outer, a inner: a is the leftmost index of
                        !   the store and must run innermost.
                        do b = 1, sumNumPartials
                           bRot = partialPerm(opR,b)
                           do a = 1, sumNumPartials
                              aRot = partialPerm(opR,a)
                              if (numStoredCompPOPTC == 1) then
                                 optcCondPOPTC(aRot,bRot,:,iE,h) = &
                                     & optcCondPOPTC(aRot,bRot,:,iE,h) &
                                     & + weight &
                                     & * transProbPOPTCBanded(a,b,:, &
                                     & cornerKP(orig),i,j,h)
                              else
                                 call depositRotatedPair (a,b,aRot, &
                                     & bRot,iE,h,i,j,cornerKP(orig), &
                                     & weight,crossFactor)
                              endif
                           enddo
                        enddo
                     else
                        do b = 1, sumNumPartials
                           do a = 1, sumNumPartials
                              if (numStoredCompPOPTC == 1) then
                                 optcCondPOPTC(a,b,:,iE,h) = &
                                     & optcCondPOPTC(a,b,:,iE,h) &
                                     & + weight &
                                     & * transProbPOPTCBanded(a,b,:, &
                                     & cornerKP(orig),i,j,h)
                              else
                                 call depositRotatedPair (a,b,a,b,iE, &
                                     & h,i,j,cornerKP(orig),weight, &
                                     & crossFactor)
                              endif
                           enddo
                        enddo
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo

end subroutine accumulateOptcCondPOPTC_LAT


! Deposit one partial pair's rotated cross term at one tetrahedron
!   corner (PSEUDOCODE 21.6).
!
! Everything that does not depend on the partial pair has already been
!   folded into crossFactor by the caller, so what is left here is the
!   contraction over the Cartesian index of the pair's own element:
!
!     entry(slot) = sum over d of Re[ M_d(a,b) * crossFactor(slot,d) ]
!
!   where crossFactor carries R(c1,d) times the conjugate of the ROTATED
!   total. Both factors of the cross term are therefore rotated, which
!   is the whole point: rotating only one of them would leave a quantity
!   that is not the projection of anything.
!
! Split out of the accumulation rather than written inline because it
!   appears twice there -- once for the permuted atom grouped codes and
!   once for the type grouped codes that need no permutation -- and two
!   hand-copied contractions would be two places to get the slot order
!   wrong.
subroutine depositRotatedPair (a,b,aRot,bRot,iE,h,initBand,finBand, &
      & kIBZ,weight,crossFactor)

   use O_Kinds
   use O_OptcTransitions, only: transMomentPOPTCBanded, &
         & bandedOccupancy, numAccumCompPOPTC

   implicit none

   ! Define the dummy variables passed to this subroutine.
   integer :: a,b          ! The pair's partials, as stored.
   integer :: aRot, bRot   ! Where they land after the permutation.
   integer :: iE, h        ! Energy point and spin.
   integer :: initBand, finBand
   integer :: kIBZ         ! The corner's irreducible k-point.
   real (kind=double) :: weight
#ifndef GAMMA
   complex (kind=double), dimension (6,3) :: crossFactor
#else
   real (kind=double), dimension (6,3) :: crossFactor
#endif

   ! Define local variables.
   integer :: slot, d
   real (kind=double) :: slotValue ! Not named "entry": that is a
         !   Fortran statement keyword, and shadowing it reads as a
         !   declaration to anyone skimming.

   do slot = 1, numAccumCompPOPTC
      slotValue = 0.0_double
      do d = 1, 3
#ifndef GAMMA
         slotValue = slotValue + real(transMomentPOPTCBanded(a,b,d, &
               & kIBZ,initBand,finBand,h) * crossFactor(slot,d),double)
#else
         slotValue = slotValue + transMomentPOPTCBanded(a,b,d,kIBZ, &
               & initBand,finBand,h) * crossFactor(slot,d)
#endif
      enddo

      ! The occupancy weight, held out of the element so that it never
      !   had to be square rooted, rejoins here. It belongs to the
      !   irreducible representative, since occupancy depends on the
      !   band and the k-point and the operation changes neither.
      optcCondPOPTC(aRot,bRot,slot,iE,h) = &
            & optcCondPOPTC(aRot,bRot,slot,iE,h) + weight * slotValue &
            & * bandedOccupancy(kIBZ,initBand,finBand,h)
   enddo

end subroutine depositRotatedPair


! Write the two column-description lines that head every unit of a
!   decomposed spectrum file (PSEUDOCODE 21.7).
!
! The declared width follows the partial direction code rather than
!   being the literal 4 it used to be, and the consumers -- imagoKKc and
!   processPOPTC.py -- are expected to read COL_LABELS rather than
!   assume. The isotropic column stays FIRST in every case, so a reader
!   that only wants it needs to know nothing about the code at all.
subroutine writePOPTCColumnLabels (outputUnit)

   use O_OptcTransitions, only: numPrintColPOPTC

   implicit none

   integer :: outputUnit

   write (outputUnit,fmt="(a11,i2)") 'COL_LABELS ',numPrintColPOPTC
   if (numPrintColPOPTC == 1) then
      write (outputUnit,fmt="(a)") 'TOTAL'
   elseif (numPrintColPOPTC == 4) then
      write (outputUnit,fmt="(a)") 'TOTAL x y z'
   else
      write (outputUnit,fmt="(a)") 'TOTAL xx yy zz xy xz yz'
   endif

end subroutine writePOPTCColumnLabels


! Write one energy point of a decomposed spectrum: the isotropic column
!   followed by whichever direction resolved columns this run carries.
!
! The isotropic value is formed here from the first three slots rather
!   than accumulated in a slot of its own. At direction code 0 the one
!   stored slot is already the SUM over the three components, so the
!   same division by three applies; at code 2 the off diagonal slots
!   must be left out of it.
subroutine writePOPTCRecord (outputUnit,slotValues)

   use O_Kinds
   use O_OptcTransitions, only: numAccumCompPOPTC, numPrintColPOPTC

   implicit none

   integer :: outputUnit
   real (kind=double), dimension (:) :: slotValues

   character*20 :: recordFormat
   real (kind=double) :: isotropic

   isotropic = sum(slotValues(1:min(3,numAccumCompPOPTC))) &
         & / 3.0_double

   write (recordFormat,fmt="(a,i2,a)") "(",numPrintColPOPTC,"e15.7)"

   if (numPrintColPOPTC == 1) then
      write (outputUnit,fmt=recordFormat) isotropic
   else
      write (outputUnit,fmt=recordFormat) isotropic, &
            & slotValues(1:numAccumCompPOPTC)
   endif

end subroutine writePOPTCRecord


! Average the atom-resolved optical spectra over the crystal point group and
!   record in the log what that did (DESIGN 1.7; PSEUDOCODE 20).
!
!   The averaging is generic and lives in O_MathSubs; what belongs here is
!   whether the permutation table exists, and reporting the outcome so that an
!   imposed equality can be told apart from an earned one.
!
!   Both indices of each pair are permuted by the same operation. That is what
!   preserves the meaning of a pair: the partial the transition starts on and
!   the one it ends on are carried to their images together, exactly as the
!   atoms they sit on are.
subroutine symmetrizeOptcPOPTC_LAT

   use O_Kinds
   use O_MathSubs,        only: symmetrizePairs
   use O_Potential,       only: spin
   use O_KPoints,         only: numPointOps
   use O_Constants,       only: dim3
   use O_OptcTransitions, only: sumNumPartials, partialPerm

   implicit none

   ! Local variables.
   integer :: h                        ! Spin loop index.
   real (kind=double) :: largestSpread ! Biggest gap made equal.
   real (kind=double) :: largestValue  ! Peak of the spectra.
   real (kind=double) :: worstSpread   ! Largest over both spins.
   real (kind=double) :: worstValue    ! Peak over both spins.
   real (kind=double) :: relativeSpread

   ! Without the symmetry maps there is nothing to average over. partialPerm
   !   is built only for the atom grouped codes and only when atomPerm exists,
   !   which it does not for an explicit kpoint list (style code 0). Say so
   !   rather than returning quietly, because a skipped symmetrization and a
   !   completed one are indistinguishable in every output file.
   if (.not. allocated(partialPerm)) then
      write (20,*) "Optical symmetrization SKIPPED: no point group maps"
      write (20,*) "are available. This happens with an explicit kpoint"
      write (20,*) "list (style code 0). Symmetry-equivalent atoms may"
      write (20,*) "not agree in the partial spectra."
      call flush (20)
      return
   endif

   worstSpread = 0.0_double
   worstValue  = 0.0_double

   ! Each spin channel is averaged on its own. The point group acts on
   !   positions and carries a spin channel onto itself, so mixing the two
   !   would average unrelated quantities together.
   do h = 1, spin
      call symmetrizePairs(optcCondPOPTC(:,:,:,:,h), &
            & partialPerm, numPointOps, sumNumPartials, &
            & dim3, numEnergyPoints, largestSpread, &
            & largestValue)
      worstSpread = max(worstSpread, largestSpread)
      worstValue  = max(worstValue,  largestValue)
   enddo

   ! Report the disagreement that was averaged away relative to the spectra it
   !   sits in. The absolute number alone would say nothing, the same gap
   !   being negligible under a tall peak and damning under a small one.
   if (worstValue > 0.0_double) then
      relativeSpread = worstSpread / worstValue
   else
      relativeSpread = 0.0_double
   endif

   write (20,*) "Partial optical spectra averaged over ",numPointOps,&
         & " point group operations."
   write (20,fmt="(a,e12.5,a,e12.5)") &
         & " Largest disagreement made equal: ",worstSpread, &
         & " of peak ",worstValue
   write (20,fmt="(a,e12.5)") &
         & " That is a relative spread of: ",relativeSpread
   write (20,*) "This equality is IMPOSED by averaging, not earned by the"
   write (20,*) "integration. A large value means the tetrahedron"
   write (20,*) "decomposition is far from point-group invariant on this"
   write (20,*) "lattice (DESIGN 1.2 and 1.7)."
   call flush (20)

end subroutine symmetrizeOptcPOPTC_LAT


subroutine printSpectrum (specType,numEnergyPoints,spectrum,conversionFactor)

   ! Include the modules we need.
   use O_Kinds
   use O_Potential,       only: spin
   use O_Constants,       only: hartree
   use O_OptcTransitions, only: energyScale, numAccumCompOPTC, &
         & numPrintColOPTC

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   integer :: specType ! 0 = XAS; 1 = conductivity; 2 = epsilon2
   integer :: numEnergyPoints
   real (kind=double), dimension (:,:,:) :: spectrum
   real (kind=double) :: conversionFactor

   ! Define the local variables
   ! Wide enough for the largest record the direction codes allow: the
   !   energy, the isotropic column and the six tensor entries, each in
   !   a fifteen character field.
   character*120 :: header
   character*20  :: headerFormat, recordFormat
   character*10  :: specTag ! "Cond", "Eps2" or "XANES".
   integer :: h,i
   integer :: unitBase
   real (kind=double) :: isotropic
   ! The suffixes naming each accumulated column, in the order
   !   PSEUDOCODE 21.7 declares and the deposit fills. Only the first
   !   numAccumCompOPTC of them are ever used.
   !
   ! Two sets, because the diagonal is named differently depending on
   !   whether anything else is present. At direction code 1 the three
   !   columns are directions and are called x, y and z. At code 2 they
   !   are the diagonal ENTRIES of a tensor whose off-diagonal entries
   !   are also present, so they are called xx, yy and zz -- and calling
   !   them x, y and z there would put this header at odds with the one
   !   imagoKKc writes for the same quantity.
   character*3, dimension (6), parameter :: axisTag = &
         & (/ "x  ", "y  ", "z  ", "   ", "   ", "   " /)
   character*3, dimension (6), parameter :: tensorTag = &
         & (/ "xx ", "yy ", "zz ", "xy ", "xz ", "yz " /)
   character*3, dimension (6) :: columnTag



   ! Customize the output for the current spectrum type. The column
   !   count follows the direction code rather than being fixed at five,
   !   so both the header and the record are built rather than written
   !   as literals (PSEUDOCODE 21.7).
   if (specType == 0) then ! XANES/ELNES
      unitBase = 49
      specTag = "XANES"
   elseif (specType == 1) then ! Optical Conductivity
      unitBase = 39
      specTag = "Cond"
   elseif (specType == 2) then ! Epsilon2
      unitBase = 49
      specTag = "Eps2"
   else
      stop "Need a spectrum type of 0, 1, or 2"
   endif

   ! One field for the energy, then the isotropic column, then whichever
   !   direction resolved columns this run carries. At direction code 0
   !   the loop below adds nothing and the record is energy plus the
   !   isotropic value alone.
   write (headerFormat,fmt="(a,i2,a)") "(",numPrintColOPTC+1,"a15)"
   write (recordFormat,fmt="(a,i2,a)") "(",numPrintColOPTC+1,"e15.7)"

   if (numAccumCompOPTC == 6) then
      columnTag(:) = tensorTag(:)
   else
      columnTag(:) = axisTag(:)
   endif

   header = ""
   if (numPrintColOPTC == 1) then
      write (header,fmt=headerFormat) "Energy","total"//trim(specTag)
   else
      write (header,fmt=headerFormat) "Energy","total"//trim(specTag),&
            & (trim(columnTag(i))//trim(specTag), &
            & i = 1, numAccumCompOPTC)
   endif

   ! Print the total (if spin == 1) or spin up and spin down (if spin == 2).
   do h = 1, spin

      ! Insert the appropriate header. Trimmed rather than written at
      !   the full buffer width so that a narrow direction code does not
      !   leave a trail of blanks behind the last column name.
      write(unitBase+h,fmt="(a)") trim(header)

      do i = 1, numEnergyPoints

         ! Adjust the spectrum for the correct units. NOTE: When this is
         !   called for printing either XANES/ELNES or the optical
         !   conductivity then the optcCond data structure is modified. For
         !   printing epsilon2, that modification side-effect is an expected
         !   prerequisite operation.
         spectrum(:,i,h) = spectrum(:,i,h) * conversionFactor / energyScale(i)

         ! The isotropic column is the average over the three Cartesian
         !   directions, and it is formed here rather than accumulated.
         !   At direction code 0 the single stored number is already the
         !   SUM over the three, so the same division by three applies;
         !   at codes 1 and 2 the first three accumulated entries are
         !   the diagonal, and the off diagonal entries of code 2 must
         !   be left out of it.
         isotropic = sum(spectrum(1:min(3,numAccumCompOPTC),i,h)) &
               & / 3.0_double

         ! Record the spectra to disk, making sure to convert the scale to eV.
         if (numPrintColOPTC == 1) then
            write (unitBase+h,fmt=recordFormat) &
                  & energyScale(i)*hartree,isotropic
         else
            write (unitBase+h,fmt=recordFormat) &
                  & energyScale(i)*hartree,isotropic, &
                  & spectrum(1:numAccumCompOPTC,i,h)
         endif
      enddo
   enddo
end subroutine printSpectrum


subroutine printSpectrumPOPTC (specType,numEnergyPoints,spectrumPOPTC,&
      & conversionFactor)

   ! Include the modules we need.
   use O_Kinds
   use O_Potential,       only: spin
   use O_Constants,       only: hartree, lAngMomCount
   use O_AtomicSites,     only: numAtomSites, atomSites
   use O_AtomicTypes,     only: numAtomTypes, atomTypes
   use O_OptcTransitions, only: energyScale, sumNumPartials, &
         & numAccumCompPOPTC
   use O_Input,           only: detailCodePOPTC

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the dummy variables passed to this subroutine.
   integer :: specType ! 0 = XAS; 1 = conductivity; 2 = epsilon2
   integer :: numEnergyPoints
   real (kind=double), dimension (:,:,:,:,:) :: spectrumPOPTC
   real (kind=double) :: conversionFactor

   ! Define the local variables
   integer :: h,i,j,k,l,m,n,o
   integer :: unitBase
   integer :: sequenceNum
   integer :: currentTypeI ! Initial
   integer :: currentTypeF ! Final
   integer :: poptcI ! pOptc initial
   integer :: poptcF ! pOptc final
   integer :: slot   ! Spectrum slot within one partial pair.
   real (kind=double), dimension (6) :: slotTotals
         !   The spectrum summed over every partial pair, one entry per
         !   slot. Sized for the widest direction code so that it does
         !   not have to be allocatable for the sake of one array.
   character*1, dimension (lAngMomCount) :: QN_lLetter

   ! Customize the output for the current spectrum type.
   if (specType == 0) then ! XANES/ELNES
      unitBase = 249
   elseif (specType == 1) then ! Optical Conductivity
      unitBase = 239
   elseif (specType == 2) then ! Epsilon2
      unitBase = 249
   else
      stop "Need a spectrum type of 0, 1, or 2"
   endif

   ! Define the QN_l letters.
   QN_lLetter(1) = 's'
   QN_lLetter(2) = 'p'
   QN_lLetter(3) = 'd'
   QN_lLetter(4) = 'f'


   ! Print the total (if spin == 1) or spin up and spin down (if spin == 2).
   do h = 1, spin

      ! Print the key bits of information for the pOptc output.
      write (unitBase+h,fmt="(a7)") 'STYLE 2'
      write (unitBase+h,fmt="(a10,i6)") 'NUM_UNITS ', (sumNumPartials**2)+1
      write (unitBase+h,fmt="(a11,i9)") 'NUM_POINTS ', numEnergyPoints
 
      ! Print the energy scale used by all atoms, converting to eV.
      do i = 1, numEnergyPoints
         write (unitBase+h,fmt="(f16.8)") energyScale(i) * hartree
      enddo

      ! Regardless of the decomposition style, we will always print the
      !   total spectrum first.
      sequenceNum = 1
      write (unitBase+h,fmt="(a13,i5)") 'SEQUENCE_NUM ',sequenceNum
      write (unitBase+h,fmt="(a20)") 'ELEMENT_1_NAME total'
      write (unitBase+h,fmt="(a20)") 'ELEMENT_2_NAME total'
      call writePOPTCColumnLabels (unitBase+h)

      do i = 1, numEnergyPoints

         ! Adjust the spectrum for the correct units. NOTE: When this is
         !   called for printing either XANES/ELNES or the optical
         !   conductivity then the optcCond data structure is *modified*. For
         !   printing epsilon2, that modification side-effect is an expected
         !   prerequisite operation. 
         spectrumPOPTC(:,:,:,i,h) = spectrumPOPTC(:,:,:,i,h) &
               & * conversionFactor / energyScale(i)

         ! Record the total spectrum to disk.
         ! The total over every partial pair, one value per slot.
         do slot = 1, numAccumCompPOPTC
            slotTotals(slot) = sum(spectrumPOPTC(:,:,slot,i,h))
         enddo
         call writePOPTCRecord (unitBase+h,slotTotals)
      enddo


      if (detailCodePOPTC == 1) then ! Decomposed by type.

         ! Print the Partial contributions to the optical conductivity.
         do i = 1, numAtomTypes
            do k  = 1, numAtomTypes

               sequenceNum = sequenceNum + 1

               ! Print the total for each pair
               write (unitBase+h,fmt="(a13,i5)") 'SEQUENCE_NUM ',sequenceNum
               write (unitBase+h,fmt="(a15,a3)") 'ELEMENT_1_NAME ',&
                     & atomTypes(i)%elementName
               write (unitBase+h,fmt="(a15,a3)") 'ELEMENT_2_NAME ',&
                     & atomTypes(k)%elementName
               write (unitBase+h,fmt="(a7,i6)") 'TYPE_1 ',i
               write (unitBase+h,fmt="(a7,i6)") 'TYPE_2 ',k
               call writePOPTCColumnLabels (unitBase+h)


               do l = 1, numEnergyPoints
                  
                  ! Record the optical conductivity to disk, making sure to
                  !   use au instead of eV.
                  call writePOPTCRecord (unitBase+h, &
                        & spectrumPOPTC(i,k,:,l,h))
               enddo
            enddo
         enddo

      ! Decompose by type and QN_nl. The pair of loop nests below walks
      !   types, then l shells, then the radial functions of that shell,
      !   which reproduces exactly the ordering that the pOptcIndex
      !   construction in O_Optc laid down, so poptcI and poptcF address
      !   the partial that was accumulated for the same (type, QN_nl)
      !   pair. Every atom carrying a type has already been summed into
      !   it, which is why no site loop appears here.
      elseif (detailCodePOPTC == 2) then

         ! Print the partial contributions to the spectrum.
         poptcI = 0
         do i = 1, numAtomTypes
            do j = 1, lAngMomCount  ! 1=s; 2=p; 3=d; 4=f
               do k = 1, atomTypes(i)%numQN_lValeRadialFns(j)
                  poptcI = poptcI + 1
                  poptcF = 0
                  do l = 1, numAtomTypes
                     do m = 1, lAngMomCount  ! 1=s; 2=p; 3=d; 4=f
                        do n = 1, atomTypes(l)%numQN_lValeRadialFns(m)
                           poptcF = poptcF + 1

                           sequenceNum = sequenceNum + 1
                   

                           ! Print the total for each pair
                           write (unitBase+h,fmt="(a13,i5)") 'SEQUENCE_NUM ',&
                                 & sequenceNum
                           write (unitBase+h,fmt="(a15,a,a1,i1,a1)") &
                                 & 'ELEMENT_1_NAME ', &
                                 & trim(atomTypes(i)%elementName), &
                                 & "_",atomTypes(i)%numQN_lCoreRadialFns(j) &
                                 & + k + j - 1, QN_lLetter(j)
                           write (unitBase+h,fmt="(a15,a,a1,i1,a1)") &
                                 & 'ELEMENT_2_NAME ', &
                                 & trim(atomTypes(l)%elementName), &
                                 & "_",atomTypes(l)%numQN_lCoreRadialFns(m) &
                                 & + n + m - 1, QN_lLetter(m)
                           write (unitBase+h,fmt="(a7,i6)") 'TYPE_1 ',i
                           write (unitBase+h,fmt="(a7,i6)") 'TYPE_2 ',l
                           call writePOPTCColumnLabels (unitBase+h)

                           do o = 1, numEnergyPoints

                              ! Record the spectrum to disk.
                              call writePOPTCRecord (unitBase+h, &
                                    & spectrumPOPTC(poptcI,poptcF,:,o,h))
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo

      ! Decompose by atom. Each atomic site owns one partial, so the
      !   loops walk sites directly and the partial index is the site
      !   index. The label carries the site number as well as the
      !   element, because two atoms of one element are separate
      !   partials here and would otherwise be indistinguishable in the
      !   output.
      elseif (detailCodePOPTC == 3) then

         ! Print the partial contributions.
         do i = 1, numAtomSites

            ! Obtain the type of the current initial state atom.
            currentTypeI = atomSites(i)%atomTypeAssn

            do k  = 1, numAtomSites

               ! Obtain the type of the current final state atom.
               currentTypeF = atomSites(k)%atomTypeAssn

               sequenceNum = sequenceNum + 1

               ! Print the total for each pair.
               write (unitBase+h,fmt="(a13,i5)") 'SEQUENCE_NUM ',sequenceNum
               write (unitBase+h,fmt="(a15,a,i0)") 'ELEMENT_1_NAME ',&
                     & trim(atomTypes(currentTypeI)%elementName), i
               write (unitBase+h,fmt="(a15,a,i0)") 'ELEMENT_2_NAME ',&
                     & trim(atomTypes(currentTypeF)%elementName), k
               write (unitBase+h,fmt="(a7,i6)") 'TYPE_1 ',currentTypeI
               write (unitBase+h,fmt="(a7,i6)") 'TYPE_2 ',currentTypeF
               call writePOPTCColumnLabels (unitBase+h)


               do l = 1, numEnergyPoints

                  ! Record the spectrum to disk.
                  call writePOPTCRecord (unitBase+h, &
                        & spectrumPOPTC(i,k,:,l,h))
               enddo
            enddo
         enddo

      ! Decompose by atom and QN_nl. This is the finest decomposition on
      !   offer: a partial is one radial function of one site, summed
      !   over the m components of its shell. The loop nests walk sites,
      !   then l shells, then the radial functions of that shell, which
      !   is exactly the order the index construction in O_Optc laid the
      !   partials down in, so poptcI and poptcF address the partial
      !   that was accumulated for the same (site, QN_nl) pair. Note
      !   that the pair count here is the square of the total radial
      !   function count over all sites, so this cell grows quickly with
      !   the cell size (DESIGN 11.4).
      elseif (detailCodePOPTC == 4) then

         ! Print the partial contributions to the spectrum.
         poptcI = 0
         do i = 1, numAtomSites

            ! Obtain the type of the current initial state atom.
            currentTypeI = atomSites(i)%atomTypeAssn

            do j = 1, lAngMomCount  ! 1=s; 2=p; 3=d; 4=f
               do k = 1, atomTypes(currentTypeI)%numQN_lValeRadialFns(j)
                  poptcI = poptcI + 1
                  poptcF = 0
                  do l = 1, numAtomSites

                     ! Obtain the type of the current final state atom.
                     currentTypeF = atomSites(l)%atomTypeAssn

                     do m = 1, lAngMomCount  ! 1=s; 2=p; 3=d; 4=f
                        do n = 1, atomTypes(currentTypeF)% &
                              & numQN_lValeRadialFns(m)
                           poptcF = poptcF + 1

                           sequenceNum = sequenceNum + 1

                           ! Print the total for each pair. The label
                           !   names the element and the site, then the
                           !   principal quantum number and QN_l letter
                           !   of the radial function, giving names such
                           !   as O3_2p.
                           write (unitBase+h,fmt="(a13,i5)") 'SEQUENCE_NUM ',&
                                 & sequenceNum
                           write (unitBase+h,fmt="(a15,a,i0,a1,i1,a1)") &
                                 & 'ELEMENT_1_NAME ', &
                                 & trim(atomTypes(currentTypeI)%elementName), &
                                 & i, "_", atomTypes(currentTypeI)% &
                                 & numQN_lCoreRadialFns(j) + k + j - 1, &
                                 & QN_lLetter(j)
                           write (unitBase+h,fmt="(a15,a,i0,a1,i1,a1)") &
                                 & 'ELEMENT_2_NAME ', &
                                 & trim(atomTypes(currentTypeF)%elementName), &
                                 & l, "_", atomTypes(currentTypeF)% &
                                 & numQN_lCoreRadialFns(m) + n + m - 1, &
                                 & QN_lLetter(m)
                           write (unitBase+h,fmt="(a7,i6)") 'TYPE_1 ',&
                                 & currentTypeI
                           write (unitBase+h,fmt="(a7,i6)") 'TYPE_2 ',&
                                 & currentTypeF
                           call writePOPTCColumnLabels (unitBase+h)

                           do o = 1, numEnergyPoints

                              ! Record the spectrum to disk.
                              call writePOPTCRecord (unitBase+h, &
                                    & spectrumPOPTC(poptcI,poptcF,:,o,h))
                           enddo ! o
                        enddo ! n
                     enddo ! m
                  enddo ! l
               enddo ! k
            enddo ! j
         enddo ! i
      endif
   enddo ! h

end subroutine printSpectrumPOPTC

end module O_OptcSpectra

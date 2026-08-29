!! SPDX-License-Identifier: ECL-2.0
!! Copyright (c) 2026 Paul Rulis

!-----------------------------------------------------------------------
! The distributed valence density build (PSEUDOCODE 36).
!
! Under the ranks-only ruling (ARCHITECTURE 6.5) the accumulate half of
!   the valence density -- the rank-k update that builds the density
!   matrix from the occupied eigenvectors (PSEUDOCODE 31) -- is
!   distributed over ranks rather than accelerated by threads.  This
!   module owns the one routine that does it, dealtValenceRankUpdate.
!
! It lives in its OWN module, not in O_ValeCharge, for the reason
!   O_ELPASolve does: O_ValeCharge already uses O_SecularEquation, and a
!   worker calls this routine from solveServerLoop, which is IN
!   O_SecularEquation -- so a routine placed in O_ValeCharge would close
!   a circular module dependency.  This module depends only on O_Kinds,
!   O_MPI and the zherk/dsyrk interface; everything else it needs -- the
!   scaled columns, the counts, the accumulator -- arrives as an
!   argument (PSEUDOCODE 36.2).
!
! Root and the workers call the SAME routine.  Root, in makeValenceRho,
!   passes the scaled columns and the counts (it alone holds them); the
!   workers, in solveServerLoop, omit those optional arguments and
!   receive the counts by broadcast and their column block by scatter.
!   Every rank does the rank-k update over ITS block into a private
!   accumulator, and the accumulators sum in place onto root
!   (MPI_IN_PLACE), where the accumulator IS the density matrix.
!
! The deal is NOT bit-exact: the reduce regroups the column sum the
!   serial zherk performs in one pass, a measured floor (DESIGN 9.6).
!-----------------------------------------------------------------------
module O_ValenceDeal

   ! Import the precision kinds shared across the code base.
   use O_Kinds

   ! Make sure that no variables are declared accidentally.
   implicit none

   contains

subroutine dealtValenceRankUpdate (rho, scaledCols, numPos, numNeg)

   ! Distribute one (k-point, spin) rank-k density update over the
   !   ranks (PSEUDOCODE 36.3): scatter each rank its column block, let
   !   each build its own partial density with a rank-k update, and sum
   !   the partials in place onto root.
   use O_Kinds
   use O_MPI, only: loadBalMPI, scattervColumnsMPI, reduceSumMPI, &
         & bcastIntVecMPI
#ifndef GAMMA
   use zherkInterface
#else
   use dsyrkInterface
#endif

   ! Make sure that no variables are declared accidentally.
   implicit none

   ! Define the passed parameters.
   !   rho        this rank's valeDim x valeDim accumulator.  On ROOT
   !              it IS valeValeRho(:,:,spin) and ends holding the full
   !              summed density; on a WORKER a scratch buffer whose
   !              only purpose is the in-place reduce.
   !   scaledCols the scaled occupied columns, valeDim by at least
   !              totalColumns; OPTIONAL -- present on ROOT, which
   !              scaled them, absent on the WORKERS.
   !   numPos     the count of positive-weight columns; OPTIONAL, as
   !              scaledCols is.
   !   numNeg     the count of negative-weight columns (possible only
   !              under the linear tetrahedron method); OPTIONAL.
#ifndef GAMMA
   complex (kind=double), dimension (:,:), intent (inout) :: rho
   complex (kind=double), dimension (:,:), intent (in), &
         & optional :: scaledCols
#else
   real (kind=double), dimension (:,:), intent (inout) :: rho
   real (kind=double), dimension (:,:), intent (in), &
         & optional :: scaledCols
#endif
   integer, intent (in), optional :: numPos
   integer, intent (in), optional :: numNeg

   ! Define local variables.
   integer :: valeDim           ! Rows of the density and the columns.
   integer :: counts(2)         ! [numPos, numNeg], broadcast to all.
   integer :: totalColumns      ! numPos + numNeg.
   integer :: firstCol, lastCol ! This rank's GLOBAL column range.
   integer :: myWidth           ! Columns in this rank's block.
   integer :: localPosCount     ! Positive columns in this rank's block.
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:) :: myColumns
#else
   real (kind=double), allocatable, dimension (:,:) :: myColumns
#endif

   ! The density's leading dimension IS the valence dimension.
   valeDim = size (rho, 1)

   ! Root has the counts; broadcast them so every rank learns the
   !   positive/negative split and the total column count.
   if (present (numPos)) then
      counts(1) = numPos
      counts(2) = numNeg
   endif
   call bcastIntVecMPI (counts)
   totalColumns = counts(1) + counts(2)

   ! Every rank learns its own column range by the SAME rule the
   !   scatter prices the blocks with (PSEUDOCODE 36.2).
   call loadBalMPI (totalColumns, firstCol, lastCol)
   myWidth = lastCol - firstCol + 1

   ! Receive this rank's column block.  Root's send buffer is
   !   scaledCols; a worker, holding no columns, passes none.  The
   !   allocation floors at one column so an empty block still has a
   !   valid array to name (its zero-column section carries nothing).
   allocate (myColumns (valeDim, max (myWidth, 1)))
   if (present (scaledCols)) then
      call scattervColumnsMPI (myColumns(:, 1:myWidth), &
            & totalColumns, scaledCols)
   else
      call scattervColumnsMPI (myColumns(:, 1:myWidth), totalColumns)
   endif

   ! Build this rank's contribution to the density.  A column whose
   !   GLOBAL index is <= numPositiveColumns is in the positive group
   !   and is ADDED; one above it is in the negative group (tetrahedron
   !   method only) and is SUBTRACTED (PSEUDOCODE 36.4).  Because the
   !   positive columns are stored first, the positive columns in a
   !   contiguous block are its first `localPosCount` columns.
   rho(:,:) = 0.0_double
   localPosCount = max (0, min (lastCol, counts(1)) - firstCol + 1)
   if (localPosCount > 0) then
#ifndef GAMMA
      call zherk ('U', 'N', valeDim, localPosCount, 1.0_double, &
            & myColumns(1,1), valeDim, 1.0_double, rho, valeDim)
#else
      call dsyrk ('U', 'N', valeDim, localPosCount, 1.0_double, &
            & myColumns(1,1), valeDim, 1.0_double, rho, valeDim)
#endif
   endif
   if (myWidth > localPosCount) then
#ifndef GAMMA
      call zherk ('U', 'N', valeDim, myWidth - localPosCount, &
            & -1.0_double, myColumns(1, localPosCount+1), valeDim, &
            & 1.0_double, rho, valeDim)
#else
      call dsyrk ('U', 'N', valeDim, myWidth - localPosCount, &
            & -1.0_double, myColumns(1, localPosCount+1), valeDim, &
            & 1.0_double, rho, valeDim)
#endif
   endif

   deallocate (myColumns)

   ! Sum every rank's partial in place onto root (MPI_IN_PLACE).  On
   !   root, rho is valeValeRho(:,:,spin), so it ends holding the full
   !   density; on a worker the summed copy is simply discarded.
   call reduceSumMPI (rho)

end subroutine dealtValenceRankUpdate

end module O_ValenceDeal

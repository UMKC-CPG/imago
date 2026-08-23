!! SPDX-License-Identifier: ECL-2.0
!! Copyright (c) 2026 Paul Rulis

!****************************************************************************
!
! The collective ELPA solve (PSEUDOCODE 27; DESIGN 9.6 stage B, 9.3).
!
! When the k-point deal cannot fill the ranks -- the one-k-point case,
!   the program's dominant use -- the single secular solve is
!   distributed: every rank holds a block-cyclic piece of H and S,
!   ELPA solves the generalized problem H c = e S c cooperatively
!   (internal Cholesky; its two-stage tridiagonalization is why it
!   beats a threaded LAPACK driver, which measurement showed to be
!   nearly thread-proof), and the eigenvector columns return to root.
!   Root's callers cannot tell which arm ran: the same full matrices
!   go in, the same eigenvalues and eigenvector columns come out.
!
! Everything ELPA-specific lives HERE, behind IMAGO_ELPA: the
!   most-square process grid with its ONE coordinate convention, the
!   block-cyclic local-extent arithmetic of DESIGN 9.3 (one spelling
!   of the layout drives the scatter, the gather, and the handle
!   setup, so the three cannot drift), and the once-per-run handle
!   lifecycle. Without IMAGO_ELPA the module compiles to
!   elpaAvailable = .false. and stubs, and the arm policy in
!   secularEqnSCF never selects the collective path.
!
!****************************************************************************
module O_ELPASolve

   ! Import the precision kinds and the MPI layer (identity, the
   !   solve-server protocol, and the block transport the scatter and
   !   gather ride on -- PSEUDOCODE 27.1: no new transport).
   use O_Kinds
   use O_MPI

#ifdef IMAGO_ELPA
   ! The ELPA Fortran API (the cpgp toolchain provides the module and
   !   library; elpa_init(20241105) is the handshake this cluster has
   !   verified).
   use elpa
#endif

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

#ifdef IMAGO_ELPA
   logical, parameter :: elpaAvailable = .true.
#else
   logical, parameter :: elpaAvailable = .false.
#endif

   ! The block-cyclic block edge (DESIGN 9.3): 64 by default, shrunk
   !   at setup when the matrix is small enough that a fixed 64 would
   !   leave a grid coordinate owning NOTHING -- the generalized
   !   solve returns garbage on such a degenerate layout (measured
   !   2026-08-22 on a 60-row problem over a 1x2 grid). Production
   !   sizes (valeDim 1800-5184 measured) keep the full 64.
   integer, parameter :: elpaBlockDefault = 64

#ifdef IMAGO_ELPA
   ! The process grid and this rank's place in it, set once per run.
   !   The coordinate convention -- column-major, processRow =
   !   mod(rank, numProcRows) -- is defined ONLY here; every layout
   !   computation below goes through it.
   logical :: elpaReady = .false.
   integer :: elpaBlockUsed ! The adaptive block edge of this run.
   integer :: numProcRows
   integer :: numProcCols
   integer :: myProcRow
   integer :: myProcCol
   integer :: myLocalRows ! Local extent of valeDim over my row coord.
   integer :: myLocalCols ! Local extent of valeDim over my col coord.
   integer :: blacsContext ! The BLACS grid over the same layout: the
         ! GENERALIZED solve path runs ScaLAPACK operations for its
         ! Cholesky transformation and refuses to run without one
         ! ("BLACS context has not been set beforehand", measured
         ! 2026-08-22 -- the standard path needs only the plain
         ! communicator, the generalized path needs both).

   ! The handle, configured once (na = valeDim, nev = numStates, the
   !   local extents, the block size, the communicator and grid
   !   coordinates) and reused every iteration.
   class (elpa_t), pointer :: elpaHandle => null ()

   ! This rank's persistent local pieces, sized at setup and reused
   !   every solve: H (overwritten), S (re-scattered and re-decomposed
   !   every iteration -- reusing the constant overlap's Cholesky via
   !   is_already_decomposed is a recorded later optimization), and
   !   the eigenvector result Q.
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:) :: myHLocal
   complex (kind=double), allocatable, dimension (:,:) :: mySLocal
   complex (kind=double), allocatable, dimension (:,:) :: myQLocal
#else
   real (kind=double), allocatable, dimension (:,:) :: myHLocal
   real (kind=double), allocatable, dimension (:,:) :: mySLocal
   real (kind=double), allocatable, dimension (:,:) :: myQLocal
#endif
   real (kind=double), allocatable, dimension (:) :: allEigenValues
#endif

contains


#ifdef IMAGO_ELPA

function localExtent (globalCount, myCoord, numCoords)

   ! How many of `globalCount` global indices land on the process at
   !   coordinate `myCoord` (0-based) of `numCoords`, under the
   !   block-cyclic deal with elpaBlockUsed (DESIGN 9.3): whole
   !   rounds of numCoords blocks give every process elpaBlockUsed
   !   indices, and the remainder falls block by block on the lowest
   !   coordinates.

   implicit none

   ! Define passed parameters and the return value.
   integer, intent (in) :: globalCount
   integer, intent (in) :: myCoord
   integer, intent (in) :: numCoords
   integer :: localExtent

   ! Define local variables.
   integer :: fullRounds
   integer :: remainder

   fullRounds = globalCount / (elpaBlockUsed * numCoords)
   remainder  = globalCount - fullRounds * elpaBlockUsed * numCoords
   localExtent = fullRounds * elpaBlockUsed &
         & + max (0, min (elpaBlockUsed,&
         &                remainder - myCoord * elpaBlockUsed))

end function localExtent


function globalIndex (localIdx, myCoord, numCoords)

   ! The global index (1-based) of local index `localIdx` on the
   !   process at coordinate `myCoord` of `numCoords` -- the inverse
   !   walk of the same block-cyclic deal localExtent counts.

   implicit none

   ! Define passed parameters and the return value.
   integer, intent (in) :: localIdx
   integer, intent (in) :: myCoord
   integer, intent (in) :: numCoords
   integer :: globalIndex

   ! Define local variables.
   integer :: localBlock
   integer :: offsetInBlock

   localBlock    = (localIdx - 1) / elpaBlockUsed
   offsetInBlock = localIdx - 1 - localBlock * elpaBlockUsed
   globalIndex = 1 + (localBlock * numCoords + myCoord)&
         & * elpaBlockUsed + offsetInBlock

end function globalIndex


subroutine ensureELPA (valeDim, numStates)

   ! First collective solve of the run: choose the most-square
   !   process grid, size this rank's locals, and configure the ELPA
   !   handle. Every later call returns immediately -- the dimensions
   !   never change within a run, so the handle is reused.

   implicit none

   ! Define passed parameters.
   integer, intent (in) :: valeDim
   integer, intent (in) :: numStates

   ! Define local variables.
   integer :: elpaErr

   if (elpaReady) return

   ! The most-square factorization of the rank count, rows from
   !   below: square grids minimize the communication volume of the
   !   factorization sweeps (DESIGN 9.3).
   numProcRows = int (sqrt (real (mpiSize, double)))
   do
      if (mod (mpiSize, numProcRows) == 0) exit
      numProcRows = numProcRows - 1
   enddo
   numProcCols = mpiSize / numProcRows
   myProcRow = mod (mpiRank, numProcRows)
   myProcCol = mpiRank / numProcRows

   ! The adaptive block edge: never so large that a grid coordinate
   !   owns nothing (see elpaBlockDefault's header; the arm policy
   !   guarantees valeDim >= mpiSize, so this is always >= 1 with
   !   every coordinate covered).
   elpaBlockUsed = max (1, min (elpaBlockDefault,&
         & valeDim / numProcRows, valeDim / numProcCols))

   myLocalRows = localExtent (valeDim, myProcRow, numProcRows)
   myLocalCols = localExtent (valeDim, myProcCol, numProcCols)

   allocate (myHLocal (max (1, myLocalRows), max (1, myLocalCols)))
   allocate (mySLocal (max (1, myLocalRows), max (1, myLocalCols)))
   allocate (myQLocal (max (1, myLocalRows), max (1, myLocalCols)))
   allocate (allEigenValues (valeDim))

   ! Initialize the library and configure the one handle of the run.
   if (elpa_init (20241105) /= ELPA_OK) then
      call stopMPI ('ELPA API version handshake failed.')
   endif
   elpaHandle => elpa_allocate (elpaErr)
   if (elpaErr /= ELPA_OK) then
      call stopMPI ('Failed to allocate the ELPA handle.')
   endif
   call elpaHandle%set ("na", valeDim, elpaErr)
   if (elpaErr /= ELPA_OK) call stopMPI ('ELPA set na failed.')
   call elpaHandle%set ("nev", numStates, elpaErr)
   if (elpaErr /= ELPA_OK) call stopMPI ('ELPA set nev failed.')
   call elpaHandle%set ("local_nrows", myLocalRows, elpaErr)
   if (elpaErr /= ELPA_OK) call stopMPI ('ELPA set nrows failed.')
   call elpaHandle%set ("local_ncols", myLocalCols, elpaErr)
   if (elpaErr /= ELPA_OK) call stopMPI ('ELPA set ncols failed.')
   call elpaHandle%set ("nblk", elpaBlockUsed, elpaErr)
   if (elpaErr /= ELPA_OK) call stopMPI ('ELPA set nblk failed.')
   call elpaHandle%set ("mpi_comm_parent", MPI_COMM_WORLD%MPI_VAL,&
         & elpaErr)
   if (elpaErr /= ELPA_OK) call stopMPI ('ELPA set comm failed.')
   call elpaHandle%set ("process_row", myProcRow, elpaErr)
   if (elpaErr /= ELPA_OK) call stopMPI ('ELPA set row failed.')
   call elpaHandle%set ("process_col", myProcCol, elpaErr)
   if (elpaErr /= ELPA_OK) call stopMPI ('ELPA set col failed.')

   ! The BLACS grid over the SAME layout, for the generalized path's
   !   internal ScaLAPACK operations. 'C' (column-major) maps rank r
   !   to (mod(r, rows), r / rows) -- exactly this module's one
   !   convention, so the BLACS view and the scattered data agree.
   !   blacs_gridinit is collective; every rank reaches this point
   !   together (the wake-before-ensure ordering guarantees it).
   call blacs_get (0, 0, blacsContext)
   call blacs_gridinit (blacsContext, 'C', numProcRows, numProcCols)
   call elpaHandle%set ("blacs_context", blacsContext, elpaErr)
   if (elpaErr /= ELPA_OK) call stopMPI ('ELPA set blacs failed.')

   elpaErr = elpaHandle%setup ()
   if (elpaErr /= ELPA_OK) call stopMPI ('ELPA setup failed.')

   ! The two-stage solver is the algorithm that makes ELPA worth the
   !   trip for LARGE dense problems (ARCHITECTURE 6.6); on small
   !   matrices its banded intermediate stage has nothing to work
   !   with and the solve errors outright (measured 2026-08-22 on a
   !   ~130-row test problem). Small problems take the one-stage
   !   solver; the threshold is generous -- every production target
   !   of this arm (valeDim 1800-5184 measured) sits far above it.
   if (valeDim >= 512) then
      call elpaHandle%set ("solver", ELPA_SOLVER_2STAGE, elpaErr)
   else
      call elpaHandle%set ("solver", ELPA_SOLVER_1STAGE, elpaErr)
   endif
   if (elpaErr /= ELPA_OK) call stopMPI ('ELPA set solver failed.')

   elpaReady = .true.

end subroutine ensureELPA


subroutine scatterLocals (fullMatrix, destRank, buffer)

   ! Fill `buffer` with `destRank`'s block-cyclic locals of
   !   `fullMatrix` -- the one walk of the layout, used for both H
   !   and S, for every rank including root itself.

   implicit none

   ! Define passed parameters.
#ifndef GAMMA
   complex (kind=double), dimension (:,:), intent (in) :: fullMatrix
   complex (kind=double), dimension (:,:), intent (out) :: buffer
#else
   real (kind=double), dimension (:,:), intent (in) :: fullMatrix
   real (kind=double), dimension (:,:), intent (out) :: buffer
#endif
   integer, intent (in) :: destRank

   ! Define local variables.
   integer :: destRow, destCol
   integer :: localRows, localCols
   integer :: il, jl

   destRow = mod (destRank, numProcRows)
   destCol = destRank / numProcRows
   localRows = localExtent (size (fullMatrix, 1), destRow,&
         & numProcRows)
   localCols = localExtent (size (fullMatrix, 2), destCol,&
         & numProcCols)
   do jl = 1, localCols
      do il = 1, localRows
         buffer(il, jl) = fullMatrix(&
               & globalIndex (il, destRow, numProcRows),&
               & globalIndex (jl, destCol, numProcCols))
      enddo
   enddo

end subroutine scatterLocals


subroutine rootCollectiveSolve (fullHQ, fullS, eigenValues)

   ! Root's side of the collective solve (PSEUDOCODE 27.3). On
   !   entry, fullHQ holds the unpacked Hamiltonian and fullS the
   !   unpacked overlap, exactly as the local LAPACK arm receives
   !   them; on return, eigenValues holds the lowest numStates
   !   eigenvalues and fullHQ's first numStates COLUMNS hold the
   !   eigenvectors -- the same contract as solveZHEGV/solveDSYGV,
   !   so the caller's write path is untouched. The gather
   !   overwrites fullHQ only after every scatter read of it is
   !   done.

   implicit none

   ! Define passed parameters.
#ifndef GAMMA
   complex (kind=double), dimension (:,:), intent (inout) :: fullHQ
   complex (kind=double), dimension (:,:), intent (inout) :: fullS
#else
   real (kind=double), dimension (:,:), intent (inout) :: fullHQ
   real (kind=double), dimension (:,:), intent (inout) :: fullS
#endif
   real (kind=double), dimension (:), intent (out) :: eigenValues

   ! Define local variables.
   integer :: valeDim, numStates
   integer :: r, il, jl, jg
   integer :: destRow, destCol
   integer :: localRows, localCols, colsUsed
   integer :: elpaErr
#ifndef GAMMA
   complex (kind=double), allocatable, dimension (:,:) :: buffer
#else
   real (kind=double), allocatable, dimension (:,:) :: buffer
#endif

   valeDim = size (fullHQ, 1)
   numStates = size (eigenValues)

   ! Wake every worker into the collective FIRST (PSEUDOCODE 27.3):
   !   the handle setup inside ensureELPA is COLLECTIVE over the
   !   world communicator (it splits the row and column
   !   communicators), so every worker must be inside its own ensure
   !   call before root enters the setup -- the first acceptance run
   !   hung on the reversed order. The k-point index is 1: the
   !   collective arm exists only at one k-point.
   do r = 1, mpiSize - 1
      call sendCtrlMPI (solveCollective, 1, r)
   enddo
   call ensureELPA (valeDim, numStates)

   ! Scatter each rank's locals of H and S (root fills its own
   !   directly). Each worker is idle in its receive when its
   !   messages arrive -- the round discipline's guarantee, width
   !   one deep here.
   do r = 0, mpiSize - 1
      destRow = mod (r, numProcRows)
      destCol = r / numProcRows
      localRows = localExtent (valeDim, destRow, numProcRows)
      localCols = localExtent (valeDim, destCol, numProcCols)
      if (r == 0) then
         call scatterLocals (fullHQ, r, myHLocal)
         call scatterLocals (fullS, r, mySLocal)
      else
         allocate (buffer (max (1, localRows), max (1, localCols)))
         call scatterLocals (fullHQ, r, buffer)
#ifndef GAMMA
         call sendCmplxBlockMPI (buffer, r, mpiTagHam)
#else
         call sendPackedMPI (buffer, r, mpiTagHam)
#endif
         call scatterLocals (fullS, r, buffer)
#ifndef GAMMA
         call sendCmplxBlockMPI (buffer, r, mpiTagOvlp)
#else
         call sendPackedMPI (buffer, r, mpiTagOvlp)
#endif
         deallocate (buffer)
      endif
   enddo

   ! The cooperative solve. ELPA returns the FULL eigenvalue vector
   !   on every rank and this rank's locals of the eigenvectors.
   call elpaHandle%generalized_eigenvectors (myHLocal, mySLocal,&
         & allEigenValues, myQLocal, .false., elpaErr)
   if (elpaErr /= ELPA_OK) then
      call stopMPI ('ELPA generalized eigensolve failed on root.')
   endif
   eigenValues(:) = allEigenValues(1:numStates)

   ! Gather the strips of the eigenvector columns with global index
   !   <= numStates back into fullHQ's leading columns. The local
   !   column -> global column map is monotone, so each rank's strip
   !   is its FIRST colsUsed local columns.
   do r = 0, mpiSize - 1
      destRow = mod (r, numProcRows)
      destCol = r / numProcRows
      localRows = localExtent (valeDim, destRow, numProcRows)
      localCols = localExtent (valeDim, destCol, numProcCols)
      colsUsed = 0
      do jl = 1, localCols
         if (globalIndex (jl, destCol, numProcCols) > numStates) exit
         colsUsed = jl
      enddo
      if (colsUsed == 0) cycle
      if (r == 0) then
         do jl = 1, colsUsed
            jg = globalIndex (jl, destCol, numProcCols)
            do il = 1, localRows
               fullHQ(globalIndex (il, destRow, numProcRows), jg) = &
                     & myQLocal(il, jl)
            enddo
         enddo
      else
         allocate (buffer (max (1, localRows), colsUsed))
#ifndef GAMMA
         call recvCmplxBlockMPI (buffer, r, mpiTagVecs)
#else
         call recvPackedMPI (buffer, r, mpiTagVecs)
#endif
         do jl = 1, colsUsed
            jg = globalIndex (jl, destCol, numProcCols)
            do il = 1, localRows
               fullHQ(globalIndex (il, destRow, numProcRows), jg) = &
                     & buffer(il, jl)
            enddo
         enddo
         deallocate (buffer)
      endif
   enddo

end subroutine rootCollectiveSolve


subroutine workerCollectiveSolve (valeDim, numStates)

   ! The worker's side (PSEUDOCODE 27.4), entered from the solve
   !   server on the solveCollective control code: receive my
   !   locals, join the cooperative solve, return my strips of the
   !   leading eigenvector columns, and go back to the server loop.

   implicit none

   ! Define passed parameters.
   integer, intent (in) :: valeDim
   integer, intent (in) :: numStates

   ! Define local variables.
   integer :: jl, colsUsed
   integer :: elpaErr

   call ensureELPA (valeDim, numStates)

#ifndef GAMMA
   call recvCmplxBlockMPI (myHLocal, 0, mpiTagHam)
   call recvCmplxBlockMPI (mySLocal, 0, mpiTagOvlp)
#else
   call recvPackedMPI (myHLocal, 0, mpiTagHam)
   call recvPackedMPI (mySLocal, 0, mpiTagOvlp)
#endif

   call elpaHandle%generalized_eigenvectors (myHLocal, mySLocal,&
         & allEigenValues, myQLocal, .false., elpaErr)
   if (elpaErr /= ELPA_OK) then
      call stopMPI ('ELPA generalized eigensolve failed on worker.')
   endif

   colsUsed = 0
   do jl = 1, myLocalCols
      if (globalIndex (jl, myProcCol, numProcCols) > numStates) exit
      colsUsed = jl
   enddo
   if (colsUsed > 0) then
#ifndef GAMMA
      call sendCmplxBlockMPI (myQLocal(:, 1:colsUsed), 0, mpiTagVecs)
#else
      call sendPackedMPI (myQLocal(:, 1:colsUsed), 0, mpiTagVecs)
#endif
   endif

end subroutine workerCollectiveSolve


subroutine teardownELPA

   ! Release the handle and the library at the end of the run. Safe
   !   to call when no collective solve ever happened.

   implicit none

   ! Define local variables.
   integer :: elpaErr

   if (.not. elpaReady) return
   call elpa_deallocate (elpaHandle, elpaErr)
   call elpa_uninit (elpaErr)
   call blacs_gridexit (blacsContext)
   deallocate (myHLocal)
   deallocate (mySLocal)
   deallocate (myQLocal)
   deallocate (allEigenValues)
   elpaReady = .false.

end subroutine teardownELPA


#else
! ---------------------------------------------------------------------
! The IMAGO_ELPA-less stubs: elpaAvailable is .false., so the arm
!   policy in secularEqnSCF never selects the collective path and
!   none of these can be reached; they exist so every caller compiles
!   identically in both builds.
! ---------------------------------------------------------------------

subroutine rootCollectiveSolve (fullHQ, fullS, eigenValues)

   implicit none

   ! Define passed parameters.
#ifndef GAMMA
   complex (kind=double), dimension (:,:), intent (inout) :: fullHQ
   complex (kind=double), dimension (:,:), intent (inout) :: fullS
#else
   real (kind=double), dimension (:,:), intent (inout) :: fullHQ
   real (kind=double), dimension (:,:), intent (inout) :: fullS
#endif
   real (kind=double), dimension (:), intent (out) :: eigenValues

   eigenValues(:) = 0.0_double
   associate (unusedA => fullHQ, unusedB => fullS)
   end associate
   call stopMPI ('Collective solve requested without ELPA built.')

end subroutine rootCollectiveSolve


subroutine workerCollectiveSolve (valeDim, numStates)

   implicit none

   ! Define passed parameters.
   integer, intent (in) :: valeDim
   integer, intent (in) :: numStates

   associate (unusedA => valeDim, unusedB => numStates)
   end associate
   call stopMPI ('Collective solve requested without ELPA built.')

end subroutine workerCollectiveSolve


subroutine teardownELPA

   implicit none

end subroutine teardownELPA

#endif

end module O_ELPASolve

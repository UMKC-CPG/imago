!! SPDX-License-Identifier: ECL-2.0
!! Copyright (c) 2026 Paul Rulis

!****************************************************************************
!
! The MPI lifecycle module (PSEUDOCODE 24; DESIGN 9.1/9.2; ARCH 6.5).
!
! This one module owns everything about the MPI runtime: starting it,
!   ending it, knowing this process's rank and the world size, dealing
!   work out across ranks, and aborting the whole job when one rank must
!   die. It compiles in BOTH builds from the same source. Under the
!   parallel build (-DIMAGO_MPI, compiler wrapper h5pfc) it is a real
!   MPI client; under the serial build it is a set of constants and
!   no-ops. That single-source property is the point: every call site
!   elsewhere in the program writes plain calls with no preprocessor
!   guard of its own, and the serial build keeps meaning exactly what
!   it meant before this module existed.
!
! The serial truths are the defaults: mpiRank = 0 and mpiSize = 1. The
!   serial build never changes them, so any consumer may read them
!   unconditionally -- "am I root" is (mpiRank == 0) in both builds,
!   and a loop balanced over mpiSize ranks degenerates to the whole
!   range at size 1. This is what lets the one-rank parallel binary and
!   the serial binary take literally identical paths, which is the
!   acceptance test for every parallel increment.
!
! Abort discipline: a plain Fortran `stop` on one rank of a parallel
!   run leaves the other ranks alive inside their next collective call,
!   and the job then hangs until the scheduler kills it. Parallel code
!   therefore aborts through stopMPI, which logs the message and takes
!   down the WHOLE job via MPI_Abort. The hundreds of existing serial
!   `stop` statements are deliberately left alone: until a path is
!   distributed, every rank executes it identically and stops together.
!   A `stop` is converted to stopMPI only when its path is parallelized
!   (PA3, PA4), never wholesale.
!
!****************************************************************************
module O_MPI

#ifdef IMAGO_MPI
   ! The Fortran 2008 MPI interface. It is imported at module scope so
   !   that every guarded procedure body below sees it; the serial
   !   build never reads this line.
   use mpi_f08
#endif

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   ! The identity of this process within the parallel run. The values
   !   below are the SERIAL truths and double as the defaults: the
   !   serial build never touches them, and the parallel build
   !   overwrites them in initMPI before any other code can look.
   integer :: mpiRank = 0 ! Rank of this process, 0-based; 0 is root.
   integer :: mpiSize = 1 ! Number of ranks in the world.

contains


subroutine initMPI

   ! Start the MPI runtime and learn who we are. This must be the FIRST
   !   call of the run (made from the program unit imagoWrapper, before
   !   Imago itself), because MPI implementations reserve the right to
   !   do anything from remapping file descriptors to re-executing the
   !   process image during initialization. In the serial build this is
   !   a no-op and the defaults above stand.

   implicit none

#ifdef IMAGO_MPI
   call MPI_Init ()
   call MPI_Comm_rank (MPI_COMM_WORLD, mpiRank)
   call MPI_Comm_size (MPI_COMM_WORLD, mpiSize)
#endif

end subroutine initMPI


subroutine closeMPI

   ! End the MPI runtime cleanly. This must be the LAST call of the run
   !   (again from imagoWrapper, after Imago returns): MPI_Finalize is
   !   collective, so every rank must reach it, and nothing MPI-related
   !   may happen after it. In the serial build this is a no-op.

   implicit none

#ifdef IMAGO_MPI
   call MPI_Finalize ()
#endif

end subroutine closeMPI


subroutine barrierMPI

   ! Wait here until every rank arrives. Serial build: a no-op, because
   !   one process is always at a barrier with itself. The one caller
   !   that must exist is the end of subroutine Imago, where the barrier
   !   gives the fort.2 success certificate its collective meaning: root
   !   writes it only after EVERY rank has finished its work.

   implicit none

#ifdef IMAGO_MPI
   call MPI_Barrier (MPI_COMM_WORLD)
#endif

end subroutine barrierMPI


subroutine stopMPI (message)

   ! The abort path for parallel code, replacing `stop` on distributed
   !   paths (see the module header for the conversion doctrine). The
   !   message is recorded first, then the whole job is taken down. We
   !   use MPI_Abort rather than a Finalize-and-stop sequence because
   !   Finalize is collective: calling it on one rank while the others
   !   are working does not end the job, it deadlocks it.

   use iso_fortran_env, only: error_unit

   implicit none

   ! Define passed parameters.
   character (len=*), intent (in) :: message

   ! Define local variables.
   logical :: logIsOpen

   ! Say it on standard error first. The MPI launcher aggregates every
   !   rank's stderr into the single job error file, so this line is
   !   the shared, rank-stamped error record (PSEUDOCODE 24.4): one
   !   destination for every rank with no file locking, and it works
   !   even when the abort comes before any log is open.
   write (error_unit, fmt = '(a,i0,2a)') "ABORT (rank ", mpiRank,&
         & "): ", message

   ! Write the same line to the log if the log is open -- on root that
   !   is the real fort.20; on a worker it is /dev/null or the rank
   !   file, per parseCommandLine. The inquire guard matters: a write
   !   to a CLOSED unit does not fail, it silently reconnects unit 20
   !   to a fresh fort.20 and truncates the real log -- exactly the
   !   evidence a person would need to see why the job died.
   inquire (unit = 20, opened = logIsOpen)
   if (logIsOpen) then
      write (20, fmt = '(a,i0,2a)') "ABORT (rank ", mpiRank, "): ",&
            & message
   endif

#ifdef IMAGO_MPI
   ! Take down every rank of the job. The error code is what the
   !   launcher (mpirun/srun) reports to the scheduler.
   call MPI_Abort (MPI_COMM_WORLD, 1)
#else
   stop 1
#endif

end subroutine stopMPI


subroutine loadBalMPI (toBalance, initialIdx, finalIdx)

   ! Deal `toBalance` items out over mpiSize ranks in contiguous
   !   ranges (DESIGN 9.2). Every rank receives jobsPer items, and the
   !   highest `remainder` ranks each take one additional item so that
   !   no work is dropped when the division is uneven. The caller then
   !   loops over its [initialIdx, finalIdx] range.
   !
   ! At mpiSize 1 this returns (1, toBalance): the serial build and the
   !   one-rank parallel build both walk the full range, which is what
   !   the acceptance checks lean on.

   implicit none

   ! Define passed parameters.
   integer, intent (in)  :: toBalance  ! Number of items to divide.
   integer, intent (out) :: initialIdx ! First item owned by this rank.
   integer, intent (out) :: finalIdx   ! Last item owned by this rank.

   ! Define local variables.
   integer :: jobsPer   ! Items every rank gets before the remainder.
   integer :: remainder ! Items left over after the even division.
   integer :: shift     ! How many extra items sit below this rank.

   jobsPer   = toBalance / mpiSize
   remainder = mod (toBalance, mpiSize)

   ! The even deal: rank r takes items [jobsPer*r + 1, jobsPer*(r+1)].
   initialIdx = jobsPer * mpiRank + 1
   finalIdx   = jobsPer * (mpiRank + 1)

   ! The highest `remainder` ranks take one extra item each. A rank in
   !   that group is shifted up by the number of extras dealt to the
   !   ranks below it, and its own range grows by one.
   if (mpiRank >= mpiSize - remainder) then
      shift = remainder - (mpiSize - mpiRank)
      initialIdx = initialIdx + shift
      finalIdx   = finalIdx + shift + 1
   endif

end subroutine loadBalMPI


end module O_MPI

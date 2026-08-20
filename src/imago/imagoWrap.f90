!! SPDX-License-Identifier: ECL-2.0
!! Copyright (c) 2026 Paul Rulis

program imagoWrapper

   use O_Imago
   use O_MPI, only: initMPI, closeMPI

   ! The MPI lifecycle brackets the entire run and lives HERE, in the
   !   only program unit, so that nothing inside O_Imago ever starts or
   !   ends the runtime (PSEUDOCODE 24.1). In the serial build both
   !   calls are no-ops and this program means what it always meant.
   call initMPI

   call Imago

   call closeMPI

end program imagoWrapper

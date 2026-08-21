!! SPDX-License-Identifier: ECL-2.0
!! Copyright (c) 2026 Paul Rulis

!****************************************************************************
!
! Generic HDF5 file-handling utilities (ARCHITECTURE 2; PSEUDOCODE 25).
!
! This is the small shared layer BENEATH the per-file owner modules
!   (O_SCFHDF5 today; the PSCF owner when post-SCF parallelization
!   arrives). It holds only operations that are about HDF5 files in
!   general, never about the layout of any particular file:
!
!   - openFileWithRetry: acquire a file whose write lock another
!     process may hold. The serial HDF5 library permits one writer
!     per file and enforces it with its own lock, and the parallel
!     term-distribution design (DESIGN 9.5) uses that lock AS the
!     write mutex -- so a failed open here is an expected, ordinary
!     event, retried until the lock is granted.
!
!   - closeAllFileObjects: close every object still open on a file
!     id. HDF5 defers the real file close until the LAST handle on
!     it dies, so a single leaked attribute would silently keep the
!     file open -- and its lock held -- after h5fclose returns
!     success. This sweep turns "the file is truly closed" from a
!     hope into a guarantee.
!
!****************************************************************************
module O_HDF5Utils

   ! Import the HDF5 interface and the precision kinds.
   use HDF5
   use O_Kinds

   ! Import the abort discipline: these routines run on parallel
   !   paths, where a plain stop would strand the other ranks.
   use O_MPI, only: stopMPI

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

contains


subroutine openFileWithRetry (fileName, fid, waitedSeconds)

   ! Open `fileName` read-write, retrying until HDF5's file lock is
   !   granted. Each miss spins about a tenth of a second on our own
   !   core -- the rank has nothing else to do while it waits for its
   !   turn at the file. The error stack is silenced around the retry
   !   because a lock-refused open is EXPECTED here, and restored
   !   before returning so that real errors elsewhere stay loud. The
   !   cap turns a pathological stall (a dead process holding the
   !   lock, a filesystem that stopped honoring locks) into a visible
   !   abort instead of a silent hang.
   !
   ! The time spent waiting is reported to the caller: it is one of
   !   the acceptance numbers of the term-distribution design
   !   (PSEUDOCODE 25.6), measured rather than assumed.

   implicit none

   ! Define passed dummy variables.
   character (len=*), intent (in) :: fileName
   integer(hid_t), intent (out) :: fid
   real (kind=double), intent (out) :: waitedSeconds

   ! Declare local variables.
   integer :: hdferr
   integer :: clockStart, clockNow, clockRate
   integer :: spinStart

   call system_clock (clockStart, clockRate)

   call h5eset_auto_f (0, hdferr)
   do
      call h5fopen_f (fileName, H5F_ACC_RDWR_F, fid, hdferr)
      if (hdferr == 0) exit

      ! Not our turn yet. Check the cap, then spin ~0.1 s and retry.
      call system_clock (clockNow)
      waitedSeconds = real (clockNow - clockStart, double)&
            & / real (clockRate, double)
      if (waitedSeconds > 3600.0_double) then
         call stopMPI ('Could not acquire the HDF5 file lock within&
               & one hour: ' // fileName)
      endif
      call system_clock (spinStart)
      do
         call system_clock (clockNow)
         if (real (clockNow - spinStart, double)&
               & / real (clockRate, double) >= 0.1_double) exit
      enddo
   enddo
   call h5eset_auto_f (1, hdferr)

   call system_clock (clockNow)
   waitedSeconds = real (clockNow - clockStart, double)&
         & / real (clockRate, double)

end subroutine openFileWithRetry


subroutine closeAllFileObjects (fid)

   ! Close every attribute, dataset, group and named datatype still
   !   open on the file `fid`. Called by a file owner just before its
   !   final h5fclose when the file must GENUINELY close -- e.g. the
   !   suspend that precedes the distributed term stage, where a
   !   leaked handle would keep the file lock held and starve every
   !   other rank's retry loop. Attributes go first (they hang off
   !   the other objects), then datasets, then groups, then named
   !   datatypes; the file handle itself is not touched.

   implicit none

   ! Define passed dummy variables.
   integer(hid_t), intent (in) :: fid

   ! Declare local variables.
   integer :: i
   integer :: hdferr
   integer :: idType
   integer(size_t) :: objCount, numReturned
   integer(hid_t), allocatable, dimension (:) :: objIds

   ! Enumerate EVERYTHING still open on the file in one query.
   !   H5F_OBJ_ALL_F covers attributes as well as datasets, groups
   !   and named datatypes (the Fortran interface offers no separate
   !   per-attribute query), so each returned id is classified with
   !   h5iget_type_f and closed by the call its class requires. The
   !   file's own handle comes back in the same list and is skipped:
   !   the owner closes it, after this sweep, with h5fclose.
   call h5fget_obj_count_f (fid, H5F_OBJ_ALL_F, objCount, hdferr)
   if (hdferr /= 0) then
      call stopMPI ('Failed to count open HDF5 objects.')
   endif

   if (objCount > 0) then
      allocate (objIds (objCount))
      call h5fget_obj_ids_f (fid, H5F_OBJ_ALL_F, objCount, objIds,&
            & hdferr, numReturned)
      if (hdferr /= 0) then
         call stopMPI ('Failed to list open HDF5 objects.')
      endif

      do i = 1, int (numReturned)
         call h5iget_type_f (objIds(i), idType, hdferr)
         if (hdferr /= 0) then
            call stopMPI ('Failed to classify an open HDF5 object.')
         endif
         if (idType == H5I_FILE_F) cycle ! The owner closes the file.
         if (idType == H5I_ATTR_F) then
            call h5aclose_f (objIds(i), hdferr)
         elseif (idType == H5I_DATASET_F) then
            call h5dclose_f (objIds(i), hdferr)
         elseif (idType == H5I_GROUP_F) then
            call h5gclose_f (objIds(i), hdferr)
         elseif (idType == H5I_DATATYPE_F) then
            call h5tclose_f (objIds(i), hdferr)
         else
            call stopMPI ('An open HDF5 object has an unknown class.')
         endif
         if (hdferr /= 0) then
            call stopMPI ('Failed to close a swept HDF5 object.')
         endif
      enddo
      deallocate (objIds)
   endif

   ! Verify the sweep left nothing behind but the file handle itself.
   !   More than that means an object slipped through, which must be
   !   surfaced rather than allowed to hold the file open.
   call h5fget_obj_count_f (fid, H5F_OBJ_ALL_F, objCount, hdferr)
   if (hdferr /= 0) then
      call stopMPI ('Failed to recount open HDF5 objects.')
   endif
   if (objCount > 1) then
      call stopMPI ('HDF5 objects remain open after the sweep.')
   endif

end subroutine closeAllFileObjects


end module O_HDF5Utils

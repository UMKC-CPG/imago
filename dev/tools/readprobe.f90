! readprobe -- time the read of ONE named HDF5 dataset, repeatedly.
!
! This is the instrument for TODO PF8 (dev/PERFORMANCE.md, "Integral
!   storage layout and read cost"; DESIGN 9.7).  The self-consistency
!   loop reads every stored integral matrix twice per iteration, and
!   the open question is what that read actually costs: the work of
!   INFLATING the DEFLATE-compressed chunk, or the time to fetch the
!   bytes from STORAGE.  The two look alike in a whole-program stamp
!   and divide across processes very differently, so the design of the
!   dealt integral cache waits on telling them apart.
!
! The method turns the page cache into the instrument instead of
!   fighting it.  The FIRST read of a dataset on a node that has never
!   touched the file is cold and includes the storage; every later
!   repeat is served from memory and leaves only inflation plus the
!   copy.  Run on two datasets of identical logical size and very
!   different on-disk size, plus uncompressed copies of each, and the
!   four warm times separate the mechanisms: if inflation dominates,
!   the compressed reads stay far above the uncompressed floor and the
!   dense dataset costs more than the sparse one; if storage dominates,
!   the warm reads all collapse toward the floor and the spread shows
!   up only in the cold read, tracking on-disk size.
!
! It deliberately uses the same library call imago's own
!   readPackedMatrix makes (h5dread_f of the whole dataset into a
!   contiguous double array), linked against the same HDF5 the serial
!   build links, so the number it reports is the number the program
!   pays and not a command-line tool's overhead.
!
! Usage:  readprobe FILE DATASET_PATH [REPEATS]
!   e.g.  readprobe gs_scf-fb.hdf5 \
!           /00001_00001_00001/atomIntgGroup/atomPotOverlap/0000001/0000001 5
!
! Output, one line per repeat, fixed leading label so a whole job's
!   worth can be pulled out with a single search:
!   READPROBE <path> repeat <k> seconds <s> logicalMB <m> MB/s <r> sum <c>
! The checksum is the plain sum of every element read; it proves the
!   read happened and that repeats returned the same data.
program readprobe

   use HDF5
   use, intrinsic :: iso_fortran_env, only: int64, real64

   implicit none

   ! Command-line arguments.
   character (len=1024) :: fileName
   character (len=1024) :: datasetPath
   character (len=32)   :: repeatArgument
   integer :: numRepeats

   ! HDF5 handles and the dataset's shape.  Fortran-order extents, so
   !   a packed real matrix arrives as (1, n(n+1)/2) and a packed
   !   complex one as (2, n(n+1)/2), exactly as imago allocates them.
   integer (hid_t) :: fileID
   integer (hid_t) :: datasetID
   integer (hid_t) :: dataspaceID
   integer (hsize_t), dimension (2) :: datasetDims
   integer (hsize_t), dimension (2) :: datasetMaxDims
   integer :: datasetRank
   integer :: hdferr

   ! The read buffer and the clock.
   real (kind=real64), allocatable, dimension (:,:) :: buffer
   integer (kind=int64) :: ticksAtStart
   integer (kind=int64) :: ticksAtEnd
   integer (kind=int64) :: ticksPerSecond
   real (kind=real64) :: elapsedSeconds
   real (kind=real64) :: logicalMegabytes
   real (kind=real64) :: checksum
   integer :: repeat

   ! Read the arguments; the repeat count defaults to one.
   call get_command_argument (1, fileName)
   call get_command_argument (2, datasetPath)
   call get_command_argument (3, repeatArgument)
   if (len_trim (fileName) == 0 .or. len_trim (datasetPath) == 0) then
      write (*,'(a)') 'usage: readprobe FILE DATASET_PATH [REPEATS]'
      stop 1
   endif
   numRepeats = 1
   if (len_trim (repeatArgument) > 0) read (repeatArgument, *) numRepeats

   ! Open the file read-only and the dataset by its full path, then ask
   !   the dataspace for the extents so the buffer matches the dataset
   !   whatever its shape.
   call h5open_f (hdferr)
   call h5fopen_f (trim (fileName), H5F_ACC_RDONLY_F, fileID, hdferr)
   if (hdferr /= 0) stop 'readprobe: cannot open file'
   call h5dopen_f (fileID, trim (datasetPath), datasetID, hdferr)
   if (hdferr /= 0) stop 'readprobe: cannot open dataset'
   call h5dget_space_f (datasetID, dataspaceID, hdferr)
   call h5sget_simple_extent_ndims_f (dataspaceID, datasetRank, hdferr)
   if (datasetRank /= 2) stop 'readprobe: expected a rank-2 dataset'
   call h5sget_simple_extent_dims_f (dataspaceID, datasetDims, &
         & datasetMaxDims, hdferr)
   allocate (buffer (datasetDims(1), datasetDims(2)))
   logicalMegabytes = real (datasetDims(1) * datasetDims(2), real64) &
         & * 8.0_real64 / 1.0e6_real64

   ! Time each whole-dataset read with the wall clock.  The first
   !   repeat is the cold read if this node has never touched the file;
   !   the rest are warm.
   call system_clock (count_rate=ticksPerSecond)
   do repeat = 1, numRepeats
      call system_clock (ticksAtStart)
      call h5dread_f (datasetID, H5T_NATIVE_DOUBLE, buffer, datasetDims, &
            & hdferr)
      call system_clock (ticksAtEnd)
      if (hdferr /= 0) stop 'readprobe: read failed'
      elapsedSeconds = real (ticksAtEnd - ticksAtStart, real64) &
            & / real (ticksPerSecond, real64)
      checksum = sum (buffer)
      write (*,'(a,1x,a,1x,a,1x,i0,1x,a,1x,f10.4,1x,a,1x,f10.2,1x,a,1x,' // &
            & 'f10.1,1x,a,1x,es22.14)') 'READPROBE', trim (datasetPath), &
            & 'repeat', repeat, 'seconds', elapsedSeconds, 'logicalMB', &
            & logicalMegabytes, 'MB/s', logicalMegabytes / elapsedSeconds, &
            & 'sum', checksum
   enddo

   call h5sclose_f (dataspaceID, hdferr)
   call h5dclose_f (datasetID, hdferr)
   call h5fclose_f (fileID, hdferr)
   call h5close_f (hdferr)

end program readprobe

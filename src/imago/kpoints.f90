module O_KPoints

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

   integer :: kPointStyleCode ! In all cases the kpoint information is read from
         !   a separate file than the imago.dat file. This separates the atomic
         !   description (basis and potential) in the imago.dat file and the
         !   kpoint mesh in a separate file. The code value in that file will
         !   cause the kpoints to be set up in different ways. If 0 then the
         !   kpoints will be read as an explicit list. If 1 then the kpoints
         !   will be read as a mesh definition along
         !      with a specific shift.
         !   If 2 then the kpoints will be read as a minimum density that
         !      will be applied along each axis. Additionally a specific
         !      shift is also read in.
   integer :: kPointIntgCode ! For each kpoint mesh it may be possible to use a
         !   variety of different integration methods from simple histograms to
         !   linear analytic tetrahedron (LAT) integration, etc. There are some
         !   constraints to be aware of. If the kpoints are given in the form of
         !   a regular grid (not necessarily square) then any LAT method can be
         !   used. A regular grid will certainly be present for kPointStyleCode
         !   values of 1 and 2, but for kPointStyleCode 0 it is possible that
         !   the mesh will not be regular. It is the responsibilty of the user
         !   to know whether or not a LAT method can be used if they provide the
         !   kPoints as an explicit list with kPointStyleCode 0. If
         !   kPointIntgCode=0 then the histogram method will be used. If
         !   kPointIntgCode=1 then the LAT method will be used.
   integer :: numKPoints ! The number of kpoints in the system.
   integer :: numPaths ! The number of discontinuous high-sym. KP paths.
   integer :: numPathKP ! Number of kpoints that will be used to create the
         !   path between the high symmetry kpoints.
   integer :: isCartesian ! 1 = yes, 0 = no (Are the high symmetry kpoints
         !   given in cartesian coordinates?)
   integer :: numTotalHighSymKP ! Total number of high symmetry kpoints over
         !   all paths.
   integer, dimension (dim3) :: numAxialKPoints ! Number of kpoints along each
         !   of the a, b, c axes defining the kpoint mesh. This is only used
         !   if the kPointStyleCode equals 1 or 2.
   integer, dimension (dim3) :: numAxialMTOP_KP ! Number of kpoints along
         !   each of the a,b,c axes defining the kpoint mesh. This is only
         !   used if either doMTOP_SCF or doMTOP_PSCF flag is set to 1.
   integer, allocatable, dimension (:) :: numHighSymKP ! Number of high
         !   symmetry kpoints that define the vertices of each path to be taken
         !   for the band diagram.
   integer, allocatable, dimension (:,:) :: mtopKPMap ! Map between the
         !   desired sequence of kpoints for MTOP calculations and the one-
         !   dimensional list of kpoints.
   real (kind=double) :: minKPointDensity ! The minimum number of kpoints per
         !   unit reciprocal-space volume (Bohr^-3). The total kpoint count is
         !   at least minKPointDensity * recipCellVolume, distributed as
         !   uniformly as possible across the three axes. Only used when
         !   kPointStyleCode equals 2.
   real (kind=double), dimension (dim3) :: kPointShift ! Fractional amount to
         !   shift the kpoint mesh by along each of the a,b,c axes. This is
         !   only used if the kPointStyleCode equals 1 or 2.
   real (kind=double), dimension (dim3) :: kPointShiftMTOP ! Fractional amount
         !   to shift the kpoint mesh by along each of the a,b,c axes. This is
         !   only used if either doMTOP_SCF or doMTOP_PSCF flag is set to 1.
   real (kind=double), allocatable, dimension (:) :: kPointWeight ! The
         !   weight assigned to each kpoint.
   real (kind=double), allocatable, dimension (:,:) :: kPoints ! The acutal
         !   kpoints.  The first dimension holds the three cartesian xyz
         !   coordinates.  The second dimension is the index over the number
         !   of kpoints.
   real (kind=double), allocatable, dimension (:,:,:) :: highSymKP ! These
         !   are the high symmetry kpoints of the requested path.  The first
         !   dimension holds the three space coordinates of the kpoints
         !   to be indexed by the second dimension (numHighSymKP). The last
         !   dimension indexes which discrete path these points define.
   real (kind=double), allocatable, dimension (:) :: pathKPointMag ! Distance
         !   magnitudes of kpoints on the path for band diagrams.
   complex (kind=double), allocatable, dimension (:,:) :: phaseFactor ! The
         !   cosine and sine of the dot product between each kpoint and each
         !   real space superlattice vector.
   character*1, allocatable, dimension (:,:) :: highSymKPChar ! The characters
         !   that identify each point on the high symmetry kpoint paths.
   integer :: numPointOps ! Number of point group operations for IBZ
         !   symmetry reduction. Read from the kpoint input file for style
         !   codes 1 and 2. These are the rotational parts of the space
         !   group operations (translations stripped).
   real (kind=double), allocatable, dimension (:,:,:) :: convAbcPointOps
         !   The 3x3 rotation matrices for each point group operation,
         !   stored in their on-disk conventional-cell-abc fractional
         !   form -- byte-identical to the entries in
         !   share/spaceDB/<sg>.  The kp file carries them verbatim;
         !   makeinput.py applies no producer-side similarity transform.
         !   computeRealPointOps and computeRecipPointOps rebase them at
         !   runtime via the change-of-basis matrix
         !   C = invRealVectors^T * convLattice (full conjugation in
         !   prim mode, identity-shortcut copy in full mode).
         !   Dimensions: (3, 3, numPointOps).  See DESIGN 2.7.
   real (kind=double), allocatable, dimension (:,:) :: convAbcFracTrans
         !   The translation vector for each point group operation, in
         !   the same on-disk conventional-cell-abc fractional form as
         !   convAbcPointOps above.  For symmorphic space groups all
         !   entries are zero; for non-symmorphic groups (e.g. Fd-3m)
         !   some operations carry non-zero translations (screw axes,
         !   glide planes).  Conjugated into the loaded real lattice's
         !   abc basis by computeRealPointOps and consumed by
         !   buildAtomPerm to form the real-space atom mapping
         !   r_B = R * r_A + t.  Dimensions: (3, numPointOps).
   real (kind=double), dimension (3, 3) :: convLattice
         !   The conventional-cell matrix in Bohr, read from the kp
         !   file's CONV_LATTICE block.  Rows are conventional lattice
         !   vectors, columns are xyz, matching the row layout of
         !   O_Lattice's realVectors.  Required by computeRealPointOps
         !   and computeRecipPointOps to form the change-of-basis
         !   matrix C = invRealVectors^T * convLattice that carries
         !   conv-abc operations into the basis of the loaded cell.
         !   For style code 0 (explicit kpoint list) this defaults to
         !   realVectors so the identity shortcut applies trivially.
         !   See DESIGN 2.7.
   character (len=4) :: cellMode
         !   Cell-mode flag read from the kp file's CELL_MODE block --
         !   either 'full' (the loaded cell IS the conventional cell:
         !   M_loaded == M_conv, so C = I and the transform routines
         !   take the identity shortcut) or 'prim' (the loaded cell is
         !   a primitive reduction of the conventional cell: C != I
         !   and the full similarity transform runs).  This mirrors the
         !   skeleton's full/prim lattice-mode flag.  Defaults to
         !   'full' for style code 0 so the trivial-identity path goes
         !   through the same shortcut.  See DESIGN 2.7.
   real (kind=double), allocatable, dimension (:,:,:) :: abcRecipPointOps
         !   The point group operations converted to reciprocal-space abc
         !   coordinates using the real and reciprocal lattice vectors.
         !   Computed by computeRecipPointOps. Dims: (3,3,numPointOps).
   real (kind=double), allocatable, dimension (:,:,:) :: abcRealPointOps
         !   The point group operations transformed into the basis of the
         !   real-space lattice currently loaded in O_Lattice. That
         !   lattice may be the full conventional cell or a primitive
         !   reduction of it -- whichever the skeleton's `full`/`prim`
         !   flag and applySpaceGroup produced. Computed by
         !   computeRealPointOps. This is the real-space sibling of
         !   abcRecipPointOps: rotations stored here act on atom
         !   positions expressed in fractional coordinates of that real
         !   lattice; rotations stored in abcRecipPointOps act on kpoints
         !   expressed in fractional coordinates of the matching
         !   reciprocal lattice. buildAtomPerm consumes these.
         !   Dimensions: (3, 3, numPointOps).
   real (kind=double), allocatable, dimension (:,:) :: abcRealFracTrans
         !   Per-operation fractional translation vectors transformed
         !   into the same real-space abc basis as abcRealPointOps -- the
         !   basis of whichever cell O_Lattice currently holds, not
         !   necessarily a primitive one. buildAtomPerm forms
         !   r_B = R * r_A + t entirely in that basis, so both the
         !   rotation matrix and the translation vector must already live
         !   in it before the routine runs. Dimensions: (3, numPointOps).
   integer :: numFullMeshKP ! The total number of kpoints in the full uniform
         !   mesh BEFORE IBZ reduction. This is the product of
         !   numAxialKPoints(1) * numAxialKPoints(2) * numAxialKPoints(3). Only
         !   set when the mesh is built internally (style codes 1, 2).
   integer, allocatable, dimension (:) :: fullKPToIBZKPMap
         !   Maps each full-mesh kpoint index (1..numFullMeshKP) to its IBZ
         !   representative kpoint index (1..numKPoints). Used by LAT
         !   integration to unfold eigenvalues from the IBZ to the full mesh:
         !   eigenValues(band, fullKPToIBZKPMap(k), spin). Only allocated when
         !   the mesh is built internally (style codes 1, 2).
   integer, allocatable, dimension (:) :: &
         & fullKPToIBZOpMap
         !   For each full-mesh kpoint index (1..numFullMeshKP), the index of
         !   the point group operation R such that R(k_IBZ) = k_full -- i.e.,
         !   the forward-direction operation that maps the IBZ representative to
         !   this full-mesh kpoint. For an IBZ representative itself, the stored
         !   value is 1 (the identity operation, which is always first in the
         !   space group database). Used together with atomPerm to distribute
         !   eigenvector-dependent properties (Q*, bond order) across the star
         !   of each IBZ kpoint. Only allocated for style codes 1 and 2.
   integer :: numTetrahedra ! The total number of tetrahedra tiling the
         !   reciprocal cell. Equal to 6 * nA * nB * nC.
   real (kind=double) :: tetraVol ! The BZ fraction for each tetrahedron: 1 /
         !   numTetrahedra.
   integer, allocatable, dimension (:,:) :: tetrahedra
         !   The four corner kpoint indices for each tetrahedron. Dimensions:
         !   (4, numTetrahedra). Indices reference the full uniform mesh
         !   (1..numFullMeshKP), NOT the IBZ-reduced list.

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine readKPoints(readUnit, writeUnit)

   ! Include the modules we need
   use O_Kinds
   use O_Constants, only: dim3

   ! Import necessary subroutine modules.
   use O_ReadDataSubs

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from
                                        ! which we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Define the local variables used in this subroutine.
   integer :: i ! Loop index variable.
   integer :: counter ! Dummy variable to read the kpoint index number.

   ! Read the kpoint definition style code. I.e., how the mesh should be
   !   constructed.
   call readData(readUnit,writeUnit,kPointStyleCode,len('KPOINT_STYLE_CODE'),&
         & 'KPOINT_STYLE_CODE')

   ! Read the code that defines the type of integration method to be used.
   !   I.e., Weighted sums, tetrahedral, etc.
   call readData(readUnit,writeUnit,kPointIntgCode,len('KPOINT_INTG_CODE'),&
         & 'KPOINT_INTG_CODE')

   ! Depending on kPointStyleCode we will read the associated relevant data.
   if (kPointStyleCode == 0) then ! Read an explicit list of kpoints.

      ! Read the number of kpoints.
      call readData(readUnit,writeUnit,numKPoints,len('NUM_BLOCH_VECTORS'),&
            & 'NUM_BLOCH_VECTORS')

      ! Read the number of kpoints requested on each axis.
      call readData(readUnit,writeUnit,3,numAxialKPoints,&
            & len('NUM_AXIAL_KPOINTS'),'NUM_AXIAL_KPOINTS')

      ! Allocate space to hold the kpoints and their weighting factors.
      allocate (kPoints(dim3,numKPoints))
      allocate (kPointWeight(numKPoints))

      ! Read past the 'NUM_WEIGHT_KA_KB_KC' header.
      call readAndCheckLabel(readUnit,writeUnit,len('NUM_WEIGHT_KA_KB_KC'),&
            & 'NUM_WEIGHT_KA_KB_KC')

      ! Note that the values read in here are in terms of a,b,c fractional
      !   coordinates. Later, after the reciprocal lattice has been initialized
      !   these values will be changed into x,y,z cartesian coordinates.
      do i = 1, numKPoints
         read (readUnit,*) counter, kPointWeight(i), kPoints(1:dim3,i)
         write (writeUnit,fmt="(i5,1x,4f15.8)") counter, kPointWeight(i), &
               & kPoints(1:dim3,i)
      enddo
      call flush (writeUnit)

   elseif (KPointStyleCode == 1) then
      ! Read axial kpoint counts, a shift, and the point group operations needed
      !   for IBZ symmetry reduction. This format is produced by makeinput.py
      !   mesh mode (-kp/-scfkp/-pscfkp).

      ! Read the number of kpoints along each a,b,c axis.
      call readData(readUnit,writeUnit,3,&
            & numAxialKPoints,len('NUM_KP_A_B_C'),&
            & 'NUM_KP_A_B_C')

      ! Read the fractional distance that the uniform mesh should be shifted
      !   away from the origin along each a,b,c axis.
      call readData(readUnit,writeUnit,3,kPointShift,&
            & len('KP_SHIFT_A_B_C'),&
            & 'KP_SHIFT_A_B_C')

      ! Read the number of point group operations for IBZ reduction.
      !   These are the rotational parts of the space group symmetry
      !   operations, expressed in their on-disk conventional-cell-abc
      !   fractional form -- byte-identical to share/spaceDB/<sg>.
      !   computeRealPointOps and computeRecipPointOps rebase them
      !   into the loaded-cell abc basis at runtime.
      call readData(readUnit,writeUnit,numPointOps,&
            & len('NUM_POINT_OPS'),'NUM_POINT_OPS')

      ! Read the point group operation matrices and their per-operation
      !   translation vectors. Each operation is a 3x3 conv-abc rotation
      !   matrix followed by a 3-component conv-abc fractional
      !   translation. The file may contain blank lines between
      !   operations for readability.
      allocate (convAbcPointOps(3,3,numPointOps))
      allocate (convAbcFracTrans(3,numPointOps))
      call readAndCheckLabel(readUnit,writeUnit,&
            & len('POINT_OPS'),'POINT_OPS')
      do i = 1, numPointOps
         ! Skip any blank line before this operation. The first
         !   operation has no blank line before it, but subsequent
         !   ones do.
         if (i > 1) read (readUnit,*)
         ! Read each row of the matrix from one file line. Each
         !   line fills one column of the Fortran array (column-major
         !   storage); see computeRealPointOps in this file for the
         !   downstream consumer that depends on this convention.
         read (readUnit,*) convAbcPointOps(:,1,i)
         read (readUnit,*) convAbcPointOps(:,2,i)
         read (readUnit,*) convAbcPointOps(:,3,i)
         ! Read the per-operation conv-abc fractional translation.
         !   Symmorphic groups carry all-zero translations; non-
         !   symmorphic groups (e.g. Fd-3m) carry non-zero values
         !   here (screw axes, glide planes).
         read (readUnit,*) convAbcFracTrans(:,i)
         do counter = 1, 3
            write (writeUnit,fmt="(3f14.8)") &
                  & convAbcPointOps(:,counter,i)
         enddo
         write (writeUnit,fmt="(3f14.8)") &
               & convAbcFracTrans(:,i)
         write (writeUnit,*)
      enddo
      call flush (writeUnit)

      ! Read the conventional-cell matrix (Bohr).  Rows are
      !   conventional lattice vectors, columns are xyz, matching
      !   the row layout of O_Lattice's realVectors.  Required by
      !   the per-routine change-of-basis matrix C =
      !   invRealVectors^T * convLattice (see DESIGN 2.7 and the
      !   computeRealPointOps docstring).
      call readAndCheckLabel(readUnit,writeUnit,&
            & len('CONV_LATTICE'),'CONV_LATTICE')
      do counter = 1, 3
         read (readUnit,*) convLattice(counter,:)
         write (writeUnit,fmt="(3f18.12)") &
               & convLattice(counter,:)
      enddo
      call flush (writeUnit)

      ! Read the cell-mode flag.  'full' means the loaded cell IS
      !   the conventional cell (so C = I and the transform routines
      !   take the identity shortcut); 'prim' means the loaded cell
      !   is a primitive reduction (so C != I and the full
      !   conjugation runs).
      call readAndCheckLabel(readUnit,writeUnit,&
            & len('CELL_MODE'),'CELL_MODE')
      read (readUnit,*) cellMode
      write (writeUnit,fmt="(a)") trim(cellMode)
      call flush (writeUnit)

      ! The kpoint mesh will be constructed later, once implicit information
      !   such as the reciprocal lattice has been obtained.

   elseif (KPointStyleCode == 2) then
      ! Read density-based kpoint specification: a minimum volume density,
      !   a shift, and the point group operations needed for IBZ symmetry
      !   reduction.

      ! Read the minimum volume density of kpoints (kpoints per unit reciprocal-
      !   space volume). Note: the file label says LINE_DENSITY for historical
      !   reasons but the value is actually a volume density.
      call readData(readUnit,writeUnit,minKPointDensity,&
            & len('MIN_KP_LINE_DENSITY'),&
            & 'MIN_KP_LINE_DENSITY')

      ! Read the fractional shift for the mesh along each a,b,c axis.
      call readData(readUnit,writeUnit,3,kPointShift,&
            & len('KP_SHIFT_A_B_C'),&
            & 'KP_SHIFT_A_B_C')

      ! Read the number of point group operations for IBZ reduction.
      !   These are the rotational parts of the space group symmetry
      !   operations along with their per-operation translations,
      !   expressed in their on-disk conventional-cell-abc fractional
      !   form -- byte-identical to share/spaceDB/<sg>.
      !   computeRealPointOps and computeRecipPointOps rebase them
      !   into the loaded-cell abc basis at runtime.
      call readData(readUnit,writeUnit,numPointOps,&
            & len('NUM_POINT_OPS'),'NUM_POINT_OPS')

      ! Read the point group operation matrices and their per-operation
      !   translation vectors. Each operation is a 3x3 conv-abc rotation
      !   matrix followed by a 3-component conv-abc fractional
      !   translation.
      allocate (convAbcPointOps(3,3,numPointOps))
      allocate (convAbcFracTrans(3,numPointOps))
      call readAndCheckLabel(readUnit,writeUnit,&
            & len('POINT_OPS'),'POINT_OPS')
      do i = 1, numPointOps
         ! Skip any blank line before this operation. The first operation has
         !   no blank line before it, but subsequent ones do.
         if (i > 1) read (readUnit,*)
         ! Read each row of the matrix from one file line. Each line fills
         !   one column of the Fortran array (column-major storage).
         read (readUnit,*) convAbcPointOps(:,1,i)
         read (readUnit,*) convAbcPointOps(:,2,i)
         read (readUnit,*) convAbcPointOps(:,3,i)
         ! Read the conv-abc fractional translation for this operation.
         read (readUnit,*) convAbcFracTrans(:,i)
         do counter = 1, 3
            write (writeUnit,fmt="(3f14.8)") &
                  & convAbcPointOps(:,counter,i)
         enddo
         write (writeUnit,fmt="(3f14.8)") &
               & convAbcFracTrans(:,i)
         write (writeUnit,*)
      enddo
      call flush (writeUnit)

      ! Read the conventional-cell matrix (Bohr).  See the matching
      !   block in the style-code-1 branch above for the layout
      !   convention and the role of this block under DESIGN 2.7.
      call readAndCheckLabel(readUnit,writeUnit,&
            & len('CONV_LATTICE'),'CONV_LATTICE')
      do counter = 1, 3
         read (readUnit,*) convLattice(counter,:)
         write (writeUnit,fmt="(3f18.12)") &
               & convLattice(counter,:)
      enddo
      call flush (writeUnit)

      ! Read the cell-mode flag.  See the matching block in the
      !   style-code-1 branch above for the role of this flag.
      call readAndCheckLabel(readUnit,writeUnit,&
            & len('CELL_MODE'),'CELL_MODE')
      read (readUnit,*) cellMode
      write (writeUnit,fmt="(a)") trim(cellMode)
      call flush (writeUnit)

      ! The kpoint mesh will be constructed later, once implicit information
      !   such as the reciprocal lattice has been obtained.
   endif


   ! Note: tetrahedra generation for LAT (kPointIntgCode==1) is deferred to
   !   initializeKPoints, after the mesh is built and numAxialKPoints is
   !   fully determined.

end subroutine readKPoints


subroutine readSYBDKPoints(readUnit, writeUnit)

   ! Include the modules we need
   use O_Kinds
   use O_Constants, only: dim3
   use O_ReadDataSubs
   use O_WriteDataSubs

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Passed parameters
   integer, intent(in)    :: readUnit   ! The unit number of the file from which
                                        ! we are reading.
   integer, intent(in)    :: writeUnit  ! The unit number of the file to which
                                        ! we are writing.

   ! Define the local variables used in this subroutine.
   integer :: i, j ! Loop index variable.

   ! A band structure is often plotted as a contiguous path of points in k-space
   !   connecting specific high-symmetry k-points. However, it is not a
   !   requirement. If one so chooses, the band structure can be plotted using a
   !   series of paths that do not have to be contiguous with each other. See
   !   Setyawan and Curtarolo, "High-throughput electronic band structure
   !   calculations: Challenges and tools", Comp. Mat. Sci., vol. 49, pages
   !   299-312, (2010).

   ! Hence, we first read in the number of distcontinuous paths to plot, the
   !   number of kpoints to use across all paths, and whether or not the set
   !   of all kpoints are given in Cartesian coordinates or not.
   call readData(readUnit,writeUnit,numPaths,numPathKP,isCartesian,&
                    len('SYBD_INPUT_DATA'),'SYBD_INPUT_DATA')
   write (writeUnit,*) 'Number of discrete paths         = ',numPaths
   write (writeUnit,*) 'Number of path K-Points          = ',numPathKP
   write (writeUnit,*) '1 = cart, 0 = fract:             = ',isCartesian
   call flush (writeUnit)

   ! Allocate space to hold the number of high symmetry kpoints that will
   !   exist for each path. E.g., there may be 6 kpoints that define the first
   !   path and only 2 kpoints that define the second path.
   allocate (numHighSymKP(numPaths))

   ! Read the number of high symmetry kpoints that define the vertices for
   !   each path.
   call readData(readUnit,writeUnit,numPaths,numHighSymKP(:),0,'')

   ! Compute the total number of high symmetry kpoints over all paths.
   numTotalHighSymKP = sum(numHighSymKP)

   write (writeUnit,*) 'Number of high symmetry K-Points = ',numTotalHighSymKP
   call flush (writeUnit)
   
   ! Allocate space to hold the high symmetry kpoints and their identifying
   !   characters.
   allocate (highSymKP(dim3,maxval(numHighSymKP),numPaths))
   allocate (highSymKPChar(maxval(numHighSymKP),numPaths))

   ! Read the high symmetry kpoints for each path.
   do i = 1, numPaths
   
      ! Read the coordinates of the high symmetry kpoints.
      do j = 1, numHighSymKP(i)
         call readData(readUnit,writeUnit,3,highSymKP(:,j,i),0,'')
      enddo
   enddo

   ! Read the set of characters identifying the high symmetry kpoints.
   do i = 1, numPaths
      do j = 1, numHighSymKP(i)
         call readData(readUnit,writeUnit,highSymKPChar(j,i),0)
      enddo
   enddo
   call writeData(writeUnit)

end subroutine readSYBDKPoints


subroutine readMTOPKPoints(readUnit, writeUnit)

   ! Include the necessary modules
   use O_Kinds
   use O_Constants, only: dim3
   use O_ReadDataSubs

   implicit none

   ! Passed parameters
   integer, intent(in) :: readUnit   ! The unit number of the file to read
   integer, intent(in) :: writeUnit  ! The unit number of the file to write

   ! Read input data
   call readAndCheckLabel(readUnit, writeUnit, len('MTOP_INPUT_DATA'), &
         & 'MTOP_INPUT_DATA')
   call readData(readUnit, writeUnit, 3, numAxialMTOP_KP(:), 0, '')
   call readData(readUnit, writeUnit, 3, kPointShiftMTOP(:), 0, '')
   call flush(writeUnit)

end subroutine readMTOPKPoints


! This converts kpoints from abc fractional (of the recip space cell) to
!   cartesian xyz (of the recip space cell).
subroutine convertKPointsToXYZ

   ! Include the modules we need
   use O_Kinds
   use O_Lattice, only: recipVectors

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the local variables used in this subroutine.
   integer :: i ! Loop index variables.
   real (kind=double) :: kPointX, kPointY, kPointZ

   ! Do the conversion. We must use a temp variable since each kpoint component
   !   is used to calculate the other components of that kpoint. I.e., we are
   !   *overwriting* the abc fractional coordinates with xyz Cartesian. The
   !   units are inverse bohr radii.
   write (20,*) 'Kpoints in x,y,z cartesian form in recip space cell are:'
   do i = 1, numKPoints
      kPointX = dot_product(kPoints(:,i),recipVectors(1,:))
      kPointY = dot_product(kPoints(:,i),recipVectors(2,:))
      kPointZ = dot_product(kPoints(:,i),recipVectors(3,:))
      kPoints(1,i) = kPointX
      kPoints(2,i) = kPointY
      kPoints(3,i) = kPointZ

      write (20,100) i,kPoints(:,i)
   enddo

   100 format (i5,3f15.8)
end subroutine convertKPointsToXYZ


! This subroutine will compute the k dot r phase factors for each kpoint and
!   real space cell in the superlattice.
subroutine computePhaseFactors

   ! Include the modules we need
   use O_Kinds
   use O_Lattice, only: numCellsReal, cellDimsReal

   ! Make sure that there are not accidental variable declarations.
   implicit none

   ! Define the local variables used in this subroutine.
   integer :: i,j
   real (kind=double) :: dotProduct

   ! First we allocate space to hold the matrix
   allocate (phaseFactor(numKPoints,numCellsReal))

   ! Assign the values as the dot product of kpoints and celldims.
!write(22,*) "Real Imaginary"
!write(24,*) "Real Imaginary"
   do i = 1, numCellsReal
      do j = 1, numKPoints
         dotProduct = dot_product(kPoints(:,j),cellDimsReal(:,i))
         phaseFactor(j,i) = cmplx(cos(dotProduct),sin(dotProduct),double)
!if ((j==2) .or. (j==26)) then
!write(20,*) "i,j",i,j,real(phaseFactor(j,i),double),aimag(phaseFactor(j,i))
!write(20,*) "cDR",cellDimsReal(:,i)
!endif
      enddo
   enddo

end subroutine computePhaseFactors


subroutine makePathKPoints

   ! Import necessary modules
   use O_Kinds
   use O_Constants, only: dim3
   use O_Lattice, only: recipVectors
   use O_TimeStamps

   ! Make sure that there are not accidental variable declarations.
   implicit none


   ! Define local parameters used in this subroutine
   integer :: i,j,k ! Loop index variables
   integer :: kPointCounter
   integer :: numSegmentKPoints
   real (kind=double) :: averageDelta
   real (kind=double), dimension (3) :: tempKPoint
   real (kind=double), dimension (3) :: segmentDelta
   real (kind=double), allocatable, dimension (:,:,:) :: symKPDistVect
   real (kind=double), allocatable, dimension (:,:)   :: symKPDistMag

   ! Begin the kpoint path creation.
   call timeStampStart(21)


   ! The main purpose of this subroutine is going to be to replace the kpoints
   !   specified in fort.15 with the kpoints defined by the path in the input
   !   file.  So, we begin by redefining the number of kpoints in the system.
   numKPoints = numPathKP

   ! Decallocate the kpoints and weights that may have been read in.
   if (allocated(kPointWeight)) then
      deallocate (kPointWeight)
      deallocate (kPoints)
   endif

   ! Reallocate the kpoint storage with the new parameters.  Note that the
   !   weights are going to store the magnitudes and not the usual kpoint
   !   weights.
   allocate (kPointWeight (numPathKP))
   allocate (kPoints      (dim3,numPathKP))

   ! Allocate space to hold the magnitudes along the path for each point.
   allocate (pathKPointMag (numPathKP))

   ! Allocate space for details of highly symmetric kpoints.
   allocate (symKPDistVect (3,maxval(numHighSymKP),numPaths))
   allocate (symKPDistMag  (maxval(numHighSymKP),numPaths))


   ! Convert the kpoints from fractional coordinates to cartesian if needed.
   !   Normally, the kpoints are given in fractional coordinates of the
   !   reciprocal space cell. This will convert the coordinates to Cartesian
   !   coordinates.
   if (isCartesian /= 1) then

      do i = 1, numPaths
         do j = 1, numHighSymKP(i)
            tempKPoint(:) = highSymKP(1,j,i) * recipVectors(:,1) + &
                          & highSymKP(2,j,i) * recipVectors(:,2) + &
                          & highSymKP(3,j,i) * recipVectors(:,3)
            highSymKP(:,j,i) = tempKPoint(:)
!write (20,*) "j,i = ",j,i
!write (20,fmt="(a,3f15.8)") "highSymKP(:,j,i)",highSymKP(:,j,i)
         enddo
      enddo
   endif

   ! Identify distances between each consecutive symmetric kpoint on the path.
   !   Note that when transitioning from one path to the next, the distance
   !   between those two kpoints should be zero.
   do i = 1, numPaths
      do j = 1, numHighSymKP(i)
         if ((i == 1) .and. (j == 1)) then ! First point on first path
            symKPDistVect(:,1,i) = 0.0_double
            symKPDistMag(1,i)    = 0.0_double
         elseif (j == 1) then ! First point on any path after the first path
               ! needs to be equal in position and magnitude to the last point
               ! on the previous path. This will signify to later programs and
               ! subroutines that we have started a new path.
            symKPDistVect(:,1,i) = symKPDistVect(:,numHighSymKP(i-1),i-1)
            symKPDistMag(1,i)    = symKPDistMag(numHighSymKP(i-1),i-1)
         else
            symKPDistVect(:,j,i) = highSymKP(:,j,i) - highSymKP(:,j-1,i)
            symKPDistMag(j,i) = symKPDistMag(j-1,i) &
                  + sqrt(sum(symKPDistVect(:,j,i)**2))
         endif
!write (20,*) "j,i = ",j,i
!write (20,fmt='(a,3f15.8)') "symKPDistVect(:,j,i)",symKPDistVect(:,j,i)
!write (20,*) "symKPDistMag(j,i)",symKPDistMag(j,i)
      enddo
   enddo

   ! Consider the special case where the number of requested path kpoints is
   !   equal to the number of highly symmetric kpoints provided.
   if (numPathKP == sum(numHighSymKP(:))) then

      ! Copy the highly symmetric kpoint coordinates to the path kpoints.
      kPointCounter = 0
      do i = 1, numPaths
         do j = 1, numHighSymKP(i)
            kPointCounter = kPointCounter + 1
            kPoints(:,kPointCounter)     = highSymKP(:,j,i)
            pathKPointMag(kPointCounter) = symKPDistMag(j,i)

            ! Record the index number of the other high symmetry k points.
            write (20,*) 'HIGH SYMMETRY K POINT INDEX NUMBER :',i
         enddo
      enddo

   else

      ! Determine the scaling factor for the distance between kpoints in a
      !   given segment between highly symmetric kpoints. This is the most
      !   most distant kpoint magnitude divided by the number of kpoints.
      averageDelta = symKPDistMag(numHighSymKP(numPaths),numPaths) / &
            & (numPathKP - 1)
!write (20,*) "numPaths = ",numPaths
!write (20,*) "numHighSymKP = ",numHighSymKP(:)
!write (20,*) "averageDelta = ",averageDelta

      ! Initialize a kpoint counter
      kPointCounter = 0

      do i = 1, numPaths
         do j = 1, numHighSymKP(i)
            if (j /= 1) then

               ! Determine the number of k points to use for the current
               !   segment between highly symmetric k points (j) and (j-1) of
               !   this path.
               numSegmentKPoints = int((symKPDistMag(j,i) - &
                     & symKPDistMag(j-1,i)) / averageDelta)

               ! Check that this number of segment kpoints will not take us
               !   over the limit. If it does, then reduce the number of
               !   points in this segment.
               if (kPointCounter + numSegmentKPoints > numPathKP) then
                  numSegmentKPoints = numPathKP - kPointCounter
               endif
!write (20,*) "j,i = ",j,i
!write (20,*) "numSegmentKPoints = ",numSegmentKPoints
!write (20,*) "symKPDistMag(j,i) = ",symKPDistMag(j,i)
!write (20,*) "symKPDistMag(j-1,i) = ",symKPDistMag(j-1,i)
!call flush(20)
   
               ! If there are not going to be any segment kpoints for this
               !   segment then cycle to the next symmetric kpoint. This is
               !   most commonly the case when we start a new path.
               if (numSegmentKPoints >= 1) then
   
                  ! Get the size of the x,y,z deltas for this segment
                  segmentDelta(:) = symKPDistVect(:,j,i) / numSegmentKPoints
      
                  ! Loop to assign the k point positions for this segment
                  do k = 1, numSegmentKPoints - 1
      
                     ! Increment the K point counter
                     kPointCounter = kPointCounter + 1
      
                     ! Store the position of the next kpoint on this segment by
                     !   adding the above determined delta to the last known
                     !   path kpoint.
                     kPoints(:,kPointCounter) = &
                           & kPoints(:,kPointCounter - 1) + segmentDelta(:)
      
                     ! Calculate the magnitude of the vector for the above path
                     !   KPoint.
                     pathKPointMag(kPointCounter) = &
                           & pathKPointMag(kPointCounter - 1) + &
                           & sqrt(sum(segmentDelta(:)**2))
                  enddo
               endif
            endif
   
            ! Now, record the j-th high symmetry k point. For the first kpoint
            !   (where i==1 and j==1) this will happen first, before the code
            !   segment above.
   
            ! Increment the K point counter.
            kPointCounter = kPointCounter + 1
   
            ! Record the index number of the high symmetry k point.
            write (20,*) 'HIGH SYMMETRY K POINT INDEX NUMBER :',kPointCounter
   
            ! Store the position of the high symmetry k point in the path array.
            kPoints(:,kPointCounter) = highSymKP(:,j,i)
   
            ! Calculate the magnitude of the vector for this high symm KPoint.
            if ((i==1) .and. (j==1)) then
               pathKPointMag(kPointCounter) = 0.0_double
            elseif (j==1) then
               pathKPointMag(kPointCounter) = pathKPointMag(kPointCounter - 1)
            else
               pathKPointMag(kPointCounter) = &
                     & pathKPointMag(kPointCounter - 1) + &
                     & sqrt(sum(segmentDelta(:)**2))
            endif
         enddo
      enddo
   endif

   ! Assign weighting factors to each kpoint on the path so that they are all
   !   equal.
   kPointWeight(:) = 2.0_double / kPointCounter

   ! Use the newly calculated number of kpoints for all future program parts.
   numPathKP   = kPointCounter
   numKPoints  = kPointCounter

   ! Record the index, magnitude and position of each kpoint on the path.
   do i = 1, kPointCounter
      write (20,fmt="(i5,4f8.3,f18.12)") i,pathKPointMag(i),kPoints(:,i),&
            & kPointWeight(i)
   enddo

   deallocate (symKPDistVect)
   deallocate (symKPDistMag)

   ! Finish the kpoint path creation.
   call timeStampEnd(21)

end subroutine makePathKPoints


subroutine initializeKPointMesh(applySymmetry)

   ! Define used modules.
   use O_Kinds
   use O_Constants

   ! Make sure no funny variables are created.
   implicit none

   ! Define passed parameters.
   integer, intent (in) :: applySymmetry

   ! Define local variables.
   integer :: i, j, k, m, n
   integer :: numMeshKPoints
   integer :: numFoldedKPoints
   integer :: isMatch
   integer, dimension (3) :: loopIndex
   integer, allocatable, dimension(:) :: kPointTracker
   real (kind=double) :: kpThresh
   real (kind=double) :: weightSum
   real (kind=double) :: initWeight
   real (kind=double), dimension(3) :: abcDelta
   real (kind=double), dimension(3) :: foldedKPoint
   real (kind=double), allocatable, dimension(:,:) :: &
         & abcMeshKPoints
   real (kind=double), allocatable, dimension(:) :: &
         & foldedWeight
   real (kind=double), allocatable, dimension(:,:) :: &
         & abcFoldedKPoints

   ! Define the threshold for when two kpoints are considered symmetrically
   !   equivalent under point group operations.
   kpThresh = 0.00001_double

   ! Compute the number of uniform mesh points.
   numMeshKPoints = product(numAxialKPoints(:))

   ! Initialize the total weight to 2. (Two electrons per state by default.
   !   Spin-polarized calculations account for the single electron per state
   !   separately via the 'spin' variable.)
   weightSum = 2.0_double

   ! Compute the mesh step size in fractional units.
   abcDelta(:) = 1.0_double / numAxialKPoints(:)

   ! Allocate space for the uniform mesh.
   allocate (abcMeshKPoints(3,numMeshKPoints))

   ! Create the initial uniform mesh in fractional units of the abc reciprocal
   !   space cell.
   numMeshKPoints = 0
   do i = 1, numAxialKPoints(1)
      loopIndex(1) = i
      do j = 1, numAxialKPoints(2)
         loopIndex(2) = j
         do k = 1, numAxialKPoints(3)
            numMeshKPoints = numMeshKPoints + 1
            loopIndex(3) = k
            abcMeshKPoints(:,numMeshKPoints) = &
                  & -0.5_double + (loopIndex(:) &
                  & - 1.0_double + kPointShift(:)) &
                  & * abcDelta(:)
         enddo
      enddo
   enddo

   ! Deallocate any previously allocated kpoint arrays.
   if (allocated(kPoints)) then
      deallocate(kPoints)
      deallocate(kPointWeight)
   endif
   if (allocated(fullKPToIBZKPMap)) then
      deallocate(fullKPToIBZKPMap)
   endif
   if (allocated(fullKPToIBZOpMap)) then
      deallocate(fullKPToIBZOpMap)
   endif

   ! Record the full mesh size for tetrahedra generation and eigenvalue
   !   unfolding.
   numFullMeshKP = numMeshKPoints

   ! If symmetry reduction is not requested, copy the full mesh directly. Build
   !   identity mappings (every point is its own IBZ representative under the
   !   identity operation).
   if (applySymmetry == 0) then
      allocate(kPoints(3,numMeshKPoints))
      allocate(kPointWeight(numMeshKPoints))
      kPoints(:,:) = abcMeshKPoints(:,:)
      kPointWeight(:) = weightSum &
            & / real(numMeshKPoints,double)
      numKPoints = numMeshKPoints
      allocate(fullKPToIBZKPMap(numMeshKPoints))
      allocate(fullKPToIBZOpMap(numMeshKPoints))
      do i = 1, numMeshKPoints
         fullKPToIBZKPMap(i) = i
         fullKPToIBZOpMap(i) = 1
      enddo
      deallocate(abcMeshKPoints)
      return
   endif

   ! ----------------------------------------------- IBZ symmetry reduction
   ! (fold the mesh). ----------------------------------------------- The
   ! algorithm applies each point group operation to each mesh kpoint.
   !   If the rotated kpoint matches another mesh kpoint (within kpThresh),
   !   those two points are symmetry-equivalent. The matching point is
   !   removed from further consideration and its weight is added to the
   !   irreducible representative.

   ! Allocate working arrays for the folding process. We allocate the full
   !   mesh size because we do not yet know how many irreducible points
   !   there will be.
   allocate (kPointTracker(numMeshKPoints))
   allocate (abcFoldedKPoints(3,numMeshKPoints))
   allocate (foldedWeight(numMeshKPoints))
   allocate (fullKPToIBZOpMap(numFullMeshKP))

   ! Compute the initial weight of each uniform mesh kpoint.
   initWeight = weightSum / real(numMeshKPoints,double)

   ! Initialize the kpoint tracker. Each mesh point starts as its own
   !   representative (tracker(i) == i). When a point is found to be equivalent
   !   to an earlier one, its tracker is set to the negative of the IBZ index of
   !   that representative.
   do i = 1, numMeshKPoints
      kPointTracker(i) = i
   enddo

   ! Initialize the count of irreducible kpoints.
   numFoldedKPoints = 0

   ! Consider each mesh kpoint in turn.
   do i = 1, numMeshKPoints

      ! Skip points already matched to an earlier IBZ representative.
      if (kPointTracker(i) /= i) cycle

      ! This is a new irreducible kpoint. The identity operation (always index 1
      !   in the space group database) maps this IBZ representative to itself.
      numFoldedKPoints = numFoldedKPoints + 1
      kPointTracker(i) = -numFoldedKPoints
      fullKPToIBZOpMap(i) = 1
      foldedWeight(numFoldedKPoints) = initWeight
      abcFoldedKPoints(:,numFoldedKPoints) = &
            & abcMeshKPoints(:,i)

      ! If this is the last mesh point, no further comparisons are needed.
      if (i == numMeshKPoints) cycle

      ! Apply each point group operation to this kpoint and look for matches
      !   among the remaining unmatched mesh points.
      do m = 1, numPointOps

         ! Compute the folded (rotated) kpoint by applying point group operation
         !   m.
         do n = 1, 3
            foldedKPoint(n) = &
                  & sum(abcRecipPointOps(n,:,m) &
                  & * abcMeshKPoints(:,i))
         enddo

         ! Compare against remaining unmatched points.
         do j = i + 1, numMeshKPoints
            if (kPointTracker(j) /= j) cycle

            ! Check if the difference is negligible along all three axes.
            isMatch = 1
            do k = 1, 3
               if (abs(foldedKPoint(k) &
                     & - abcMeshKPoints(k,j)) &
                     & > kpThresh) then
                  isMatch = 0
                  exit
               endif
            enddo

            ! If matched, absorb this point's weight and mark it as folded.
            !   Record that operation m maps the IBZ representative to this
            !   full-mesh kpoint (forward direction: m(k_IBZ) = k_full(j)).
            if (isMatch == 1) then
               foldedWeight(numFoldedKPoints) = &
                     & foldedWeight(numFoldedKPoints) &
                     & + initWeight
               kPointTracker(j) = -numFoldedKPoints
               fullKPToIBZOpMap(j) = m
            endif

         enddo ! j (remaining mesh points)
      enddo ! m (point group operations)
   enddo ! i (mesh points)

   ! Copy the irreducible kpoints into the module-level kPoints and kPointWeight
   !   arrays.
   allocate(kPoints(3,numFoldedKPoints))
   allocate(kPointWeight(numFoldedKPoints))
   kPoints(:,1:numFoldedKPoints) = &
         & abcFoldedKPoints(:,1:numFoldedKPoints)
   kPointWeight(1:numFoldedKPoints) = &
         & foldedWeight(1:numFoldedKPoints)
   numKPoints = numFoldedKPoints

   ! Build the full-to-IBZ mapping from kPointTracker. kPointTracker(i) ==
   !   -ibzIndex for every mesh point. Convert to positive IBZ indices.
   allocate(fullKPToIBZKPMap(numFullMeshKP))
   do i = 1, numFullMeshKP
      fullKPToIBZKPMap(i) = -kPointTracker(i)
   enddo

   ! Clean up working arrays.
   deallocate(abcMeshKPoints)
   deallocate(kPointTracker)
   deallocate(abcFoldedKPoints)
   deallocate(foldedWeight)

end subroutine initializeKPointMesh


subroutine initializeKPoints (inSCF)

   ! Define used modules.
   use O_Kinds
   use O_Constants, only: dim3
   use O_CommandLine, only: doSYBD_SCF, doSYBD_PSCF, doMTOP_SCF, doMTOP_PSCF
   use O_Lattice, only: realVectors
         ! Used by the style-code-0 trivial-identity setup to seed
         !   convLattice with the currently-loaded cell, so that
         !   computeRealPointOps' identity shortcut applies (see
         !   DESIGN 2.7).

   ! Prevent accidental variable definition.
   implicit none

   ! Declare passed parameters.
   integer, intent(in) :: inSCF

   ! Define local variables.
   integer :: i ! Loop index for identity-map setup (style code 0).

   ! The considerations for this operation are:

   ! At this point an SCF or PSCF kpoint input file has been read in. The
   !   imago.dat file has been read in. Implicit information has been computed.
   !   The lattices have been initialized.

   ! If either doSYBD_SCF or doSYBD_PSCF variable was set, then we need to
   !   shift to use the defined SYBD kpoints regardless of any other setting.

   ! If either doMTOP_SCF or doMTOP_PSCF variable was set, then we need to
   !   shift to use the defined MTOP kpoints regardless of any other setting.
   !   The MTOP kpoints are a uniform mesh given as a number of kpoints in
   !   each of the a,b,c directions along with the amount of shift from the
   !   origin.

   ! If none of doSYBD_SCF, doSYBD_PSCF, doMTOP_SCF, and doMTOP_PSCF are set,
   !   then we turn to the kpoints as defined in the kpoint input file.
   !   Otherwise, if one of the above situations was true, then there is
   !   nothing left to do and we can leave this subroutine.

   ! If the kpoint style was set to zero, then the kpoints were explicitly
   !   provided in the kpoint input file, the number of kpoints was set,
   !   the kpoint weights were set, and the kpoint values themselves were
   !   set. Therefore, there is nothing left to do and we can leave this
   !   subroutine.

   ! If the kpoint style was set to one, then the number of kpoints in each
   !   a,b,c direction were provided along with the amount of shift that
   !   should be applied to the mesh. Therefore, it is necessary to create
   !   a mesh from this information and reduce it by symmetry.

   ! If the kpoint style was set to two, then the density of kpoints was
   !   provided along with the amount of shift that should be applied to the
   !   mesh. Therefore, it is necessary to create a mesh from this information
   !   and reduce it by symmetry.

   ! Proceed:

   ! Check for the above sequence of situations.
   if ((doSYBD_SCF == 1) .and. (inSCF == 1)) then
      call makePathKPoints
   elseif ((doSYBD_PSCF == 1) .and. (inSCF == 0)) then
      call makePathKPoints
   elseif ((doMTOP_SCF == 1) .and. (inSCF == 1)) then
      numAxialKPoints(:) = numAxialMTOP_KP(:)
      kPointShift(:) = kPointShiftMTOP(:)
      call initializeKPointMesh(0) ! Do not apply symmetry.
      call convertKPointsToXYZ
      call makeMTOPIndexMap
   elseif ((doMTOP_PSCF == 1) .and. (inSCF == 0)) then
      numAxialKPoints(:) = numAxialMTOP_KP(:)
      kPointShift(:) = kPointShiftMTOP(:)
      call initializeKPointMesh(0) ! Do not apply symmetry.
      call convertKPointsToXYZ
      call makeMTOPIndexMap
   elseif (kPointStyleCode == 0) then

      ! Warn that decomposition properties may not be correct with an explicit
      !   k-point list. The code sets up trivial identity IBZ maps (each k-point
      !   is its own IBZ representative with star size 1 and identity
      !   permutation), so the machinery works, but the results are only
      !   physically correct if the user-supplied list is a complete symmetric
      !   mesh.
      write (20,*)
      write (20,*) '================================'
      write (20,*) '=          WARNING             ='
      write (20,*) '================================'
      write (20,*) 'Style-code-0 k-points (explicit'
      write (20,*) 'list) are in use. Decomposition'
      write (20,*) 'properties (effective charge,'
      write (20,*) 'bond order, partial DOS) may'
      write (20,*) 'NOT be correct unless the list'
      write (20,*) 'is a complete symmetric mesh.'
      write (20,*) 'Use -kp or -kpd in makeinput.py'
      write (20,*) 'for guaranteed correct results.'
      write (20,*) '================================'
      write (20,*)

      call convertKPointsToXYZ

      ! Set up trivial identity IBZ infrastructure so that the star-distribution
      !   code in computeBond works uniformly for all style codes. Each k-point
      !   is its own IBZ representative (star size 1) with the identity point
      !   group operation.
      numPointOps = 1
      allocate (convAbcPointOps(3, 3, 1))
      convAbcPointOps(:,:,1) = 0.0_double
      convAbcPointOps(1,1,1) = 1.0_double
      convAbcPointOps(2,2,1) = 1.0_double
      convAbcPointOps(3,3,1) = 1.0_double
      allocate (convAbcFracTrans(3, 1))
      convAbcFracTrans(:,1) = 0.0_double

      ! Style code 0 never sees a CONV_LATTICE / CELL_MODE block on
      !   disk because there are no real symmetry operations to
      !   rebase. Default to 'full' with convLattice = realVectors so
      !   computeRealPointOps' identity shortcut applies trivially --
      !   M_loaded == M_conv gives C = I and the single identity op
      !   passes through unchanged. See DESIGN 2.7 for the shortcut
      !   convention.
      cellMode = 'full'
      convLattice(:,:) = realVectors(:,:)

      ! Build the real-abc twin of the (trivial) identity operation in
      !   the basis of whichever cell is currently loaded. For a single
      !   identity op the transform is a no-op, but invoking
      !   computeRealPointOps unconditionally keeps every code path that
      !   consumes abcRealPointOps free of style-code special cases.
      call computeRealPointOps

      numFullMeshKP = numKPoints
      allocate (fullKPToIBZKPMap(numKPoints))
      allocate (fullKPToIBZOpMap(numKPoints))
      do i = 1, numKPoints
         fullKPToIBZKPMap(i) = i
         fullKPToIBZOpMap(i) = 1
      enddo
   elseif (kPointStyleCode == 1) then
      call computeRecipPointOps
      call computeRealPointOps
      call initializeKPointMesh(1) ! Apply symmetry.
      call convertKPointsToXYZ
   elseif (kPointStyleCode == 2) then
      call computeAxialKPoints
      call computeRecipPointOps
      call computeRealPointOps
      call initializeKPointMesh(1) ! Apply symmetry.
      call convertKPointsToXYZ
   endif

   ! If the LAT integration method was requested, generate the tetrahedra
   !   and compute the BZ integration weight per tetrahedron. This must
   !   happen after the mesh is built (so that numAxialKPoints is set).
   !   Tetrahedra reference the full uniform mesh, not the IBZ-reduced
   !   kpoints.
   if (kPointIntgCode == 1) then
      call generateTetrahedra
      call computeTetraVol
   endif

   ! Compute the phase factors.
   call computePhaseFactors

end subroutine initializeKPoints


subroutine computeAxialKPoints

   ! Define used modules.
   use O_Kinds
   use O_Constants, only:dim3
   use O_Lattice, only: recipCellVolume, recipMag

   ! Define local variables.
   integer :: i
   integer :: nearIndex
   real (kind=double) :: stepSize
   real (kind=double) :: numKPointsEstimate
   real (kind=double), dimension (dim3) :: fractionalPoints

   ! The user specifies a volume density (kpoints per unit reciprocal-space
   !   volume in Bohr^-3). Multiply by the reciprocal cell volume to get the
   !   target total number of kpoints in the full mesh.
   numKPointsEstimate = minKPointDensity * recipCellVolume

   ! Distribute the total count across the three axes so that the spacing is
   !   as uniform as possible (equal step size along each axis).

   ! To force the step sizes to be equal for each axis, the ratio of the
   !   axial magnitude and the number of kpoints along that axis should be
   !   equal for all axes. I.e., recipMag(1) / numAxialKPoints(1) = stepSize,
   !   and so forth for 2, and 3. With a fixed stepSize this would yield
   !   real numbers for the number of kpoints. We will assume this for now
   !   and then make integer later.

   ! Now, establish a point density constraint: recipCellVolume /
   !   product(numAxialKPoints(:)) = 1 / minKPointsDensity

   ! Therefore, the step size is:
   stepSize = (product(recipMag(:)) / &
         & (recipCellVolume * minKPointDensity))**(1.0_double/3.0_double)

   ! Now, solve for the number of axial points:
   numAxialKPoints(:) = int(recipMag(:) / stepSize)

   ! Compute the residual fractional part.
   fractionalPoints(:) = (recipMag(:) / stepSize) &
         & - real(numAxialKPoints(:),double)

   ! Finally, adjust for the shift to integer values.

   ! First, increment the number of kpoints along axes that are closest to
   !   the next integer number of points until the density exceeds the target.
   do while (minKPointDensity * recipCellVolume > product(numAxialKPoints(:)))
      nearIndex = maxloc(fractionalPoints(:),1)
      numAxialKPoints(nearIndex) = numAxialKPoints(nearIndex) + 1
      fractionalPoints(nearIndex) = 0.0_double
   enddo

   ! Then, decrement the number of kpoints along the axes that are closest to
   !   the previous integer (floor) until the density is less than the target.
   !   First, fix any fractionalPoints that are zero to be one so that they
   !   are never reduced.
   do i = 1, 3
      if (fractionalPoints(i) == 0.0_double) then
         fractionalPoints(i) = 1.0_double
      endif
   enddo
   do while (minKPointDensity * recipCellVolume < product(numAxialKPoints(:)))
      nearIndex = minloc(fractionalPoints(:),1)
      numAxialKPoints(nearIndex) = numAxialKPoints(nearIndex) - 1
      if (numAxialKPoints(nearIndex) == 0) then
         numAxialKPoints(nearIndex) = 1
      endif
      fractionalPoints(nearIndex) = 1.0_double
   enddo

   ! At this point, the numAxialKPoints should be optimized to be as close as
   !   possible to the target density while keeping inter-point spacing as
   !   constant as possible (comparing axes).

end subroutine computeAxialKPoints


subroutine makeMTOPIndexMap

   ! Use necessary modules

   ! Make sure no variables are accidentally declared.
   implicit none

   ! Declare local variables.
   integer :: kPointCount
   integer :: i,j,k

   ! Allocate space to hold the MTOP index map.
   allocate (mtopKPMap(3,numKPoints))

   ! Build a map between the linear list of kpoints as they were created in
   !   the initializeKPointMesh subroutine and the order in which the c-axis
   !   strings of kpoints are needed for MTOP calculations.
   kPointCount = 0
   do i = 1, numAxialKPoints(1)
      do j = 1, numAxialKPoints(2)
         do k = 1, numAxialKPoints(3)
            kPointCount = kPointCount + 1
            mtopKPMap(3,kPointCount) = kPointCount
         enddo
      enddo
   enddo

   ! Build a map between the linear list of kpoints as they were created in
   !   the initializeKPointMesh subroutine and the order in which the b-axis
   !   strings of kpoints are needed for MTOP calculations.
   kPointCount = 0
   do i = 1, numAxialKPoints(1)
      do j = 1, numAxialKPoints(3)
         do k = 1, numAxialKPoints(2)
            kPointCount = kPointCount + 1
            mtopKPMap(2,kPointCount) = getIndexFromIndices(i,k,j)
         enddo
      enddo
   enddo

   ! Build a map between the linear list of kpoints as they were created in
   !   the initializeKPointMesh subroutine and the order in which the a-axis
   !   strings of kpoints are needed for MTOP calculations.
   kPointCount = 0
   do i = 1, numAxialKPoints(2)
      do j = 1, numAxialKPoints(3)
         do k = 1, numAxialKPoints(1)
            kPointCount = kPointCount + 1
            mtopKPMap(1,kPointCount) = getIndexFromIndices(k,i,j)
         enddo
      enddo
   enddo

!write(20,*) "AXIS1"
!do i = 1,kPointCount
!   write(20,*) "i map kp",i,mtopKPMap(1,i), kPoints(:,i)
!enddo
!write(20,*) "AXIS2"
!do i = 1,kPointCount
!   write(20,*) "i map kp",i,mtopKPMap(2,i), kPoints(:,i)
!enddo
!write(20,*) "AXIS3"
!do i = 1,kPointCount
!   write(20,*) "i map kp",i,mtopKPMap(3,i), kPoints(:,i)
!enddo

end subroutine makeMTOPIndexMap


! Given the a,b,c indices of a kpoint in the full kpoint mesh, this function
!   will return then index of the same kpoint in the linear list of kpoints.
function getIndexFromIndices(a,b,c)

   implicit none

   ! Passed parameters.
   integer, intent(in) :: a, b, c
   integer :: getIndexFromIndices

   ! The linear list of kpoints is created by nested loops over kpoints in
   !   the a, b, c axes in that order. Hence, the c-axis kpoints iterate
   !   fastest through the linear list, followed by the b-axis kpoints. So
   !   the linear list index is the current c index + nAKP(2) times the
   !   number of full lines on the c-axis we passed to get to the current
   !   b-axis point, followed by nAKP(2)*nAKP(3) times the number of full
   !   c,b planes we passed to get to the current a-axis point.
   getIndexFromIndices = c + (b-1)*numAxialKPoints(3) + &
         & (a-1)*numAxialKPoints(2)*numAxialKPoints(3)

end function getIndexFromIndices


! Generate the tetrahedra that tile the reciprocal cell. The uniform
!   Monkhorst-Pack mesh defines a grid of nA x nB x nC parallelepipeds. Each
!   parallelepiped has 8 corners and is decomposed into exactly 6 tetrahedra
!   that share the main diagonal M1-M8 (Bloechl 1994). The mesh is periodic, so
!   indices wrap with modular arithmetic.
!
!   Parallelepiped corners at grid position (a, b, c): M1 = (a, b, c ) M5 =
!     (a+1, b+1, c ) M2 = (a+1, b, c ) M6 = (a+1, b, c+1) M3 = (a, b+1, c ) M7 =
!     (a, b+1, c+1) M4 = (a, b, c+1) M8 = (a+1, b+1, c+1)
!
!   Six tetrahedra sharing diagonal M1-M8: T1: M1, M2, M5, M8 T4: M1, M4, M7, M8
!     T2: M1, M3, M5, M8 T5: M1, M4, M6, M8 T3: M1, M3, M7, M8 T6: M1, M2, M6,
!     M8
!
!   The tetrahedra indices reference the full uniform mesh (1..numFullMeshKP),
!   not the IBZ-reduced list.
subroutine generateTetrahedra

   implicit none

   ! Define local variables.
   integer :: a, b, c, t
   integer :: nA, nB, nC
   integer :: M1, M2, M3, M4, M5, M6, M7, M8

   ! Local shorthand for the axial kpoint counts.
   nA = numAxialKPoints(1)
   nB = numAxialKPoints(2)
   nC = numAxialKPoints(3)

   ! Compute the total number of tetrahedra: 6 per parallelepiped, one
   !   parallelepiped per grid cell.
   numTetrahedra = 6 * nA * nB * nC

   ! Allocate the tetrahedra array.
   if (allocated(tetrahedra)) deallocate(tetrahedra)
   allocate (tetrahedra(4, numTetrahedra))

   ! Build the tetrahedra by iterating over every grid cell and decomposing it
   !   into 6 tetrahedra.
   t = 0
   do a = 1, nA
      do b = 1, nB
         do c = 1, nC

            ! Compute the 8 corner indices of this parallelepiped with periodic
            !   wrapping.
            M1 = getIndexFromIndices( &
                  & a, b, c)
            M2 = getIndexFromIndices( &
                  & mod(a, nA) + 1, b, c)
            M3 = getIndexFromIndices( &
                  & a, mod(b, nB) + 1, c)
            M4 = getIndexFromIndices( &
                  & a, b, mod(c, nC) + 1)
            M5 = getIndexFromIndices( &
                  & mod(a, nA) + 1, &
                  & mod(b, nB) + 1, c)
            M6 = getIndexFromIndices( &
                  & mod(a, nA) + 1, &
                  & b, mod(c, nC) + 1)
            M7 = getIndexFromIndices( &
                  & a, mod(b, nB) + 1, &
                  & mod(c, nC) + 1)
            M8 = getIndexFromIndices( &
                  & mod(a, nA) + 1, &
                  & mod(b, nB) + 1, &
                  & mod(c, nC) + 1)

            ! Six tetrahedra sharing diagonal M1-M8.
            tetrahedra(:, t+1) = (/M1, M2, M5, M8/)
            tetrahedra(:, t+2) = (/M1, M3, M5, M8/)
            tetrahedra(:, t+3) = (/M1, M3, M7, M8/)
            tetrahedra(:, t+4) = (/M1, M4, M7, M8/)
            tetrahedra(:, t+5) = (/M1, M4, M6, M8/)
            tetrahedra(:, t+6) = (/M1, M2, M6, M8/)
            t = t + 6

         enddo ! c
      enddo ! b
   enddo ! a

end subroutine generateTetrahedra


! Compute the BZ integration weight for each tetrahedron. All tetrahedra tile
!   the BZ uniformly, so each one represents an equal fraction 1/numTetrahedra
!   of the full zone.
subroutine computeTetraVol

   use O_Kinds

   implicit none

   tetraVol = 1.0_double &
         & / real(numTetrahedra, double)

end subroutine computeTetraVol


! Convert the on-disk point group operations (read into
!   convAbcPointOps from the kp file in their native conventional-
!   cell-abc fractional form) into the basis of the reciprocal
!   lattice currently held in O_Lattice. The resulting
!   abcRecipPointOps act directly on kpoints expressed in
!   fractional reciprocal coordinates of that lattice and are used
!   by initializeKPointMesh for runtime IBZ symmetry reduction.
!
! Input convention. The operations arrive here byte-identical to
!   share/spaceDB/<sg>: each rotation matrix is in the natural
!   axes of the conventional crystallographic setting and acts on
!   conv-abc fractional vectors directly. No producer-side
!   similarity transform has been applied; the entire basis change
!   happens here on the consumer side. See DESIGN 2.7 for the
!   motivation behind locating the conjugation here.
!
! Output convention. Whichever cell ended up in O_Lattice after
!   the input pipeline (the full conventional cell when the
!   skeleton's flag is `full`, a primitive reduction when it is
!   `prim`), the conjugation below uses realVectors and
!   invRealVectors as they currently stand together with the
!   conventional cell read from the kp file's CONV_LATTICE block.
!   The result is correct for any cell imago can describe.
!
! Mathematics. The change-of-basis matrix that carries a vector
!   from the conv-abc basis to the loaded-cell abc basis is
!
!     C = M_loaded^{-1} * M_conv
!       = invRealVectors^T * convLattice,
!
!   where the equality uses the orthogonality identity
!   realVectors^T * recipVectors = 2*pi*I, which makes
!   realVectors^{-1} equal to invRealVectors^T (no extra
!   inversion of realVectors at runtime).  The point-group
!   rotation transforms by the similarity
!
!     R_recip_abc = C * R_conv_abc * C^{-1}.
!
!   Point-group rotations transform identically under the dual
!   real- and reciprocal-space abc bases, so this matches the
!   real-space transform in computeRealPointOps below.  There is
!   no translation field to transform on the reciprocal side.
!
! Identity shortcut. In `full` mode the loaded cell IS the
!   conventional cell, so M_loaded == M_conv and C = I; the
!   conjugation collapses to a copy.  The cellMode flag (read
!   from the kp file's CELL_MODE block) selects between the
!   shortcut and the full conjugation path.
subroutine computeRecipPointOps

   ! Import necessary modules.
   use O_Kinds
   use O_Lattice, only: invRealVectors

   ! Make sure no funny variables are defined.
   implicit none

   ! Define local variables.
   integer :: i ! Loop index over point group operations.
   real (kind=double), dimension (3,3) :: changeBasis
         !   Composed change-of-basis matrix C = invRealVectors^T
         !   * convLattice.  Computed once before the per-op loop;
         !   acts as the conjugating factor on the left of each
         !   similarity transform.
   real (kind=double), dimension (3,3) :: changeBasisInv
         !   Inverse of changeBasis (C^{-1}), used on the right
         !   of each similarity transform.  Computed once before
         !   the per-op loop.

   ! Allocate storage for the reciprocal-space operations.
   allocate (abcRecipPointOps(3,3,numPointOps))

   if (cellMode == 'full') then
      ! Identity shortcut: M_loaded == M_conv, so C = I and the
      !   conjugation collapses to a copy from the on-disk array
      !   into the runtime array.  Skips the matrix inverse and
      !   the per-op kernel for the most common cell mode.
      do i = 1, numPointOps
         abcRecipPointOps(:,:,i) = convAbcPointOps(:,:,i)
      enddo
   else
      ! Full conjugation path.  Form C once, then apply the
      !   similarity R_recip_abc = C * R_conv_abc * C^{-1} per
      !   operation.  invRealVectors^T equals realVectors^{-1}
      !   by the orthogonality identity (see header), so the
      !   factor on the left of C uses pre-computed O_Lattice
      !   state without a runtime inversion of realVectors.
      changeBasis = matmul(transpose(invRealVectors), convLattice)
      call invert3x3(changeBasis, changeBasisInv)

      do i = 1, numPointOps
         abcRecipPointOps(:,:,i) = matmul( &
               & changeBasis, &
               & matmul(convAbcPointOps(:,:,i), &
               &        changeBasisInv))
      enddo
   endif

end subroutine computeRecipPointOps


! Convert the on-disk point group operations and their fractional
!   translation vectors (read into convAbcPointOps and
!   convAbcFracTrans from the kp file in conv-abc form) into the
!   basis of the real-space lattice currently held in O_Lattice.
!   The resulting abcRealPointOps and abcRealFracTrans act directly
!   on atom positions expressed in fractional real-space
!   coordinates of that same lattice, which is how buildAtomPerm
!   represents them internally (see atomicSites.f90 step 1).
!
! Input convention. The operations arrive here byte-identical to
!   share/spaceDB/<sg>: each rotation matrix and per-operation
!   translation vector is in the natural axes of the conventional
!   crystallographic setting and acts on conv-abc fractional
!   vectors directly. No producer-side similarity transform has
!   been applied; the entire basis change happens here on the
!   consumer side. See DESIGN 2.7 for the motivation behind
!   locating the conjugation here.
!
! Output convention. Whichever cell ended up in O_Lattice after
!   the input pipeline (the full conventional cell when the
!   skeleton's flag is `full`, a primitive reduction when it is
!   `prim`), the conjugation below uses realVectors and
!   invRealVectors as they currently stand together with the
!   conventional cell read from the kp file's CONV_LATTICE block.
!   The result is correct for any cell imago can describe.
!
! Mathematics. The change-of-basis matrix that carries a vector
!   from the conv-abc basis to the loaded-cell abc basis is
!
!     C = M_loaded^{-1} * M_conv
!       = invRealVectors^T * convLattice,
!
!   where the equality uses the orthogonality identity
!   realVectors^T * recipVectors = 2*pi*I, which makes
!   realVectors^{-1} equal to invRealVectors^T (no extra
!   inversion of realVectors at runtime).  The rotation and the
!   translation then transform as
!
!     R_real_abc = C * R_conv_abc * C^{-1},
!     t_real_abc = C * t_conv_abc.
!
!   Both identities hold for any non-singular realVectors and
!   convLattice; no cell-shape, centering, or symmetry assumption
!   is built into the derivation, and the routine is therefore
!   correct for every cell type imago can describe.
!
! Identity shortcut. In `full` mode the loaded cell IS the
!   conventional cell, so M_loaded == M_conv and C = I; the
!   conjugation collapses to a copy of both the rotation matrices
!   and the translation vectors.  The cellMode flag (read from
!   the kp file's CELL_MODE block) selects between the shortcut
!   and the full conjugation path.
!
! Relationship to computeRecipPointOps above. That routine
!   performs the matching transform for the reciprocal-space side
!   using the same change-of-basis matrix C, so kpoint folding
!   via abcRecipPointOps and atom permutation via abcRealPointOps
!   descend from a single source of operations and a single
!   change-of-basis step.
subroutine computeRealPointOps

   ! Import necessary modules.
   use O_Kinds
   use O_Lattice, only: invRealVectors

   ! Make sure no funny variables are defined.
   implicit none

   ! Define local variables.
   integer :: i ! Loop index over point group operations.
   real (kind=double), dimension (3,3) :: changeBasis
         !   Composed change-of-basis matrix C = invRealVectors^T
         !   * convLattice.  Computed once before the per-op loop;
         !   acts as the conjugating factor on the left of each
         !   similarity transform and as the multiplier on each
         !   conv-abc translation.
   real (kind=double), dimension (3,3) :: changeBasisInv
         !   Inverse of changeBasis (C^{-1}), used on the right
         !   of each rotation's similarity transform.

   ! Allocate storage for both transformed quantities. The arrays
   !   are freed in cleanUpKPoints alongside the other point-group
   !   tables.
   allocate (abcRealPointOps(3, 3, numPointOps))
   allocate (abcRealFracTrans(3, numPointOps))

   if (cellMode == 'full') then
      ! Identity shortcut: M_loaded == M_conv, so C = I and the
      !   conjugation collapses to a copy from the on-disk arrays
      !   into the runtime arrays.  Skips the matrix inverse and
      !   the per-op kernel for the most common cell mode.
      do i = 1, numPointOps
         abcRealPointOps(:,:,i)  = convAbcPointOps(:,:,i)
         abcRealFracTrans(:,i)   = convAbcFracTrans(:,i)
      enddo
   else
      ! Full conjugation path.  Form C once, then apply the
      !   similarity per operation:
      !     R_real_abc = C * R_conv_abc * C^{-1}
      !     t_real_abc = C * t_conv_abc
      !   invRealVectors^T equals realVectors^{-1} by the
      !   orthogonality identity (see header), so the factor on
      !   the left of C uses pre-computed O_Lattice state without
      !   a runtime inversion of realVectors.
      changeBasis = matmul(transpose(invRealVectors), convLattice)
      call invert3x3(changeBasis, changeBasisInv)

      do i = 1, numPointOps
         abcRealPointOps(:,:,i) = matmul( &
               & changeBasis, &
               & matmul(convAbcPointOps(:,:,i), &
               &        changeBasisInv))
         abcRealFracTrans(:,i) = matmul( &
               & changeBasis, convAbcFracTrans(:,i))
      enddo
   endif

end subroutine computeRealPointOps


! Invert a 3x3 real matrix using cofactor expansion.  Used by
!   computeRealPointOps and computeRecipPointOps to form
!   C^{-1} from C = invRealVectors^T * convLattice when the
!   `prim` cell-mode path is taken (the `full` path skips this
!   entirely via the identity shortcut).
!
! A non-singular convLattice and realVectors are preconditions of
!   the calling code, so a singular C here means the kp file's
!   CONV_LATTICE block is corrupt or the loaded real lattice has
!   collapsed -- either points to an input problem rather than a
!   bug in this routine.  The check uses |det| < 1.0e-14 as the
!   singularity threshold, matching the producer-side check that
!   used to live in makeinput.py's _inv_3x3 helper.
subroutine invert3x3(matrixIn, matrixOut)

   ! Import necessary modules.
   use O_Kinds

   implicit none

   ! Passed parameters.
   real (kind=double), dimension(3,3), intent(in)  :: matrixIn
   real (kind=double), dimension(3,3), intent(out) :: matrixOut

   ! Local variables.
   real (kind=double) :: det
         !   Determinant of matrixIn.  Used to check for
         !   singularity and to divide out of the cofactor matrix.
   real (kind=double) :: invDet
         !   Reciprocal of det, cached so each output entry is a
         !   single multiply.

   ! Compute the determinant via expansion along the first row.
   det = matrixIn(1,1) * (matrixIn(2,2)*matrixIn(3,3) &
       &                - matrixIn(2,3)*matrixIn(3,2)) &
       & - matrixIn(1,2) * (matrixIn(2,1)*matrixIn(3,3) &
       &                  - matrixIn(2,3)*matrixIn(3,1)) &
       & + matrixIn(1,3) * (matrixIn(2,1)*matrixIn(3,2) &
       &                  - matrixIn(2,2)*matrixIn(3,1))

   if (abs(det) < 1.0e-14_double) then
      write (20,*) 'invert3x3: singular matrix (|det| =', &
            & abs(det), '); check kp file CONV_LATTICE and the'
      write (20,*) 'loaded real lattice.  This is an input/state'
      write (20,*) 'problem rather than a bug in invert3x3.'
      stop
   endif

   invDet = 1.0_double / det

   ! Cofactor / adjugate construction: matrixOut(i,j) = (-1)^{i+j}
   !   * minor(j,i) / det.  Note the transpose: the (i,j) entry
   !   of the inverse comes from the (j,i) cofactor.
   matrixOut(1,1) =  (matrixIn(2,2)*matrixIn(3,3) &
                  & - matrixIn(2,3)*matrixIn(3,2)) * invDet
   matrixOut(1,2) = -(matrixIn(1,2)*matrixIn(3,3) &
                  & - matrixIn(1,3)*matrixIn(3,2)) * invDet
   matrixOut(1,3) =  (matrixIn(1,2)*matrixIn(2,3) &
                  & - matrixIn(1,3)*matrixIn(2,2)) * invDet
   matrixOut(2,1) = -(matrixIn(2,1)*matrixIn(3,3) &
                  & - matrixIn(2,3)*matrixIn(3,1)) * invDet
   matrixOut(2,2) =  (matrixIn(1,1)*matrixIn(3,3) &
                  & - matrixIn(1,3)*matrixIn(3,1)) * invDet
   matrixOut(2,3) = -(matrixIn(1,1)*matrixIn(2,3) &
                  & - matrixIn(1,3)*matrixIn(2,1)) * invDet
   matrixOut(3,1) =  (matrixIn(2,1)*matrixIn(3,2) &
                  & - matrixIn(2,2)*matrixIn(3,1)) * invDet
   matrixOut(3,2) = -(matrixIn(1,1)*matrixIn(3,2) &
                  & - matrixIn(1,2)*matrixIn(3,1)) * invDet
   matrixOut(3,3) =  (matrixIn(1,1)*matrixIn(2,2) &
                  & - matrixIn(1,2)*matrixIn(2,1)) * invDet

end subroutine invert3x3


subroutine cleanUpKPoints

   implicit none

   deallocate (kPointWeight)
   deallocate (kPoints)
   deallocate (numHighSymKP)
   deallocate (highSymKP)
   deallocate (highSymKPChar)

   if (allocated(pathKPointMag)) then
      deallocate(pathKPointMag)
   endif

   ! Only allocated in setup.
   if (allocated(phaseFactor)) then
      deallocate (phaseFactor)
   endif

   ! Only allocated for MTOP.
   if (allocated(mtopKPMap)) then
      deallocate (mtopKPMap)
   endif

   ! On-disk conv-abc operations and their per-op translations.
   !   Allocated for every style code (style 0 sets up a single
   !   identity operation in memory; styles 1 and 2 read them from
   !   the kp file).
   if (allocated(convAbcPointOps)) then
      deallocate (convAbcPointOps)
   endif
   if (allocated(convAbcFracTrans)) then
      deallocate (convAbcFracTrans)
   endif
   if (allocated(abcRecipPointOps)) then
      deallocate (abcRecipPointOps)
   endif
   if (allocated(abcRealPointOps)) then
      deallocate (abcRealPointOps)
   endif
   if (allocated(abcRealFracTrans)) then
      deallocate (abcRealFracTrans)
   endif

   ! Only allocated when the mesh is built internally (style codes 1, 2) with or
   !   without IBZ reduction.
   if (allocated(fullKPToIBZKPMap)) then
      deallocate (fullKPToIBZKPMap)
   endif
   if (allocated(fullKPToIBZOpMap)) then
      deallocate (fullKPToIBZOpMap)
   endif

   ! Only allocated when kPointIntgCode == 1 (LAT).
   if (allocated(tetrahedra)) then
      deallocate (tetrahedra)
   endif

end subroutine cleanUpKPoints


end module O_KPoints

!! SPDX-License-Identifier: ECL-2.0
!! Copyright (c) 2026 Paul Rulis

program imagoKkc
!Written by Paul Rulis.  Adapted from EPS1 by Dr. Y.-N. Xu.
!Last modified on April 1, 2022 by Alysse Weigand and Paul Rulis

!The purpose of the program is to calculate epsilon1 and the energy loss
! function from the given epsilon2 using the Kramers-Kronig conversion
! relation.

!The usual procedure for calculating the energy loss function (ELF) is as
! follows.  First form a new finer grained set of energies by linear
! interpolation.  Form a finer grained version of the epsilon2 data.  Perform
! Kramers Kronig conversion to determine epsilon1.  Use epsilon1 and epsilon2
! to find the ELF.

   use O_Kinds
   use O_Constants

   implicit none

   real (kind=double), allocatable, dimension(:) :: energy, fine_energy
   real (kind=double), allocatable, dimension(:) :: totalEps1,totalEps1i
   real (kind=double), allocatable, dimension(:) :: totalEps2,totalElf
   real (kind=double), allocatable, dimension(:) :: totalRefractIndex
   real (kind=double), allocatable, dimension(:) :: totalExtnctCoeff
   real (kind=double), allocatable, dimension(:) :: totalReflectivity
   real (kind=double), allocatable, dimension(:) :: totalAbsorpCoeff
   real (kind=double), allocatable, dimension(:,:) :: eps1,eps1i,elf
   real (kind=double), allocatable, dimension(:,:) :: refractIndex,extnctCoeff
   real (kind=double), allocatable, dimension(:,:) :: reflectivity,absorpCoeff
   real (kind=double), allocatable, dimension(:,:) :: eps2,fine_eps2
   
   real (kind=double) :: coarseEnergyDiff,fineEnergyDiff,pOptcFactor
   integer :: usePOPTCFiles ! Flag for using POPTC file numbers.
   integer :: numValues=0  ! This value will be set later to length*grain.
   integer :: grain=10     ! Grainularity factor for refining the integration.

   ! The variables that will be defined by a command line parameters
   integer :: length  ! Length of the optical conductivity file
   integer :: spin    ! 1 for spin up or total, 2 for spin down.

   ! How many DIRECTION resolved columns the input carries, and how many
   !   columns this program therefore works in. Imago's optical writers
   !   choose their record width from the direction code of DESIGN 13.7,
   !   so it is no longer fixed at three and has to be discovered from
   !   the file rather than assumed.
   !
   !   code 0   energy, total                        numDirCols = 0
   !   code 1   energy, total, x, y, z               numDirCols = 3
   !   code 2   energy, total, xx yy zz xy xz yz     numDirCols = 6
   !
   ! TWO widths are needed, and confusing them is the trap this program
   !   has to avoid (DESIGN 13.7, PSEUDOCODE 21.7).
   !
   ! numWorkCols sizes epsilon-2 and epsilon-1, which carry every
   !   component the input holds. Kramers-Kronig acts on each component
   !   independently, because causality holds element by element, so
   !   epsilon-1 generalizes to all six without argument.
   !
   ! numScalarCols sizes the five quantities that are functions of a
   !   SCALAR dielectric function -- the energy loss function, the
   !   refractive index, the extinction coefficient, the reflectivity
   !   and the absorption. None of those has a meaning applied to an
   !   off-diagonal entry, so they are computed from the DIAGONAL
   !   entries alone and numScalarCols never exceeds three.
   !
   ! The isotropic average runs over numScalarCols, never over
   !   numWorkCols. At code 2 the two differ, and averaging over all six
   !   would fold the off-diagonal entries into a quantity that is meant
   !   to be one third of the TRACE -- finite, plausible, and wrong.
   integer :: numDirCols
   integer :: numWorkCols
   integer :: numScalarCols
   integer :: numFields    ! Numeric fields on one input data line.
   character*400 :: probeLine

   character*40 :: buffer ! Character buffer for command line arguments

   ! Define the subroutine interfaces for proper passing of allocatable arrays.
   interface
      subroutine readOPTCData (length, energy, totalEps2, eps2)
         use O_Kinds
         use O_Constants
         integer :: length
         real (kind=double), dimension (:) :: energy
         real (kind=double), dimension (:) :: totalEps2
         real (kind=double), dimension (:,:) :: eps2
      end subroutine readOPTCData
      subroutine getCoarseEnergyDiff(length, energy, coarseEnergyDiff)
         use O_Kinds
         use O_Constants
         integer :: length
         real (kind=double), dimension (:) :: energy
         real (kind=double) :: coarseEnergyDiff
      end subroutine getCoarseEnergyDiff
      subroutine makeFineGrainEnergy(length,grain,energy,fine_energy)
         use O_Kinds
         use O_Constants
         integer :: length, grain
         real (kind=double), dimension (:) :: energy, fine_energy
      end subroutine makeFineGrainEnergy
      subroutine getFineEnergyDiff(numValues,grain,fine_energy,fineEnergyDiff)
         use O_Kinds
         use O_Constants
         integer :: numValues, grain
         real (kind=double), dimension (:) :: fine_energy
         real (kind=double) :: fineEnergyDiff
      end subroutine getFineEnergyDiff
      subroutine makeFineEps2(length,grain,energy,fine_energy,eps2,&
                             &fine_eps2)
         use O_Kinds
         use O_Constants
         integer :: length,grain
         real (kind=double), dimension (:) :: energy,fine_energy
         real (kind=double), dimension (:,:) :: eps2,fine_eps2
      end subroutine makeFineEps2
      subroutine kramersKronig(length,grain,numValues,fineEnergyDiff,&
                              &energy,fine_energy,eps1,eps1i,fine_eps2,&
                              &pOptcFactor)
         use O_Kinds
         use O_Constants
         integer :: length,grain,numValues
         real (kind=double) :: fineEnergyDiff,pOptcFactor
         real (kind=double), dimension (:) :: energy,fine_energy
         real (kind=double), dimension (:,:) :: eps1,eps1i,fine_eps2
      end subroutine kramersKronig
      subroutine get_n_k_R_a(length,energy,eps1,eps2,refractIndex,&
            & extnctCoeff,reflectivity,absorpCoeff)
         use O_Kinds
         use O_Constants
         integer :: length
         real (kind=double), dimension (:) :: energy
         real (kind=double), dimension (:,:) :: eps1,eps2
         real (kind=double), dimension (:,:) :: refractIndex,extnctCoeff
         real (kind=double), dimension (:,:) :: reflectivity,absorpCoeff
      end subroutine get_n_k_R_a
      subroutine getELF(length,eps1,eps2,elf)
         use O_Kinds
         use O_Constants
         integer :: length
         real (kind=double), dimension (:,:) :: eps1,eps2,elf
      end subroutine getELF
      subroutine averageFunctions(length,eps1,eps1i,refractIndex,extnctCoeff,&
                                & reflectivity,absorpCoeff,totalEps1,&
                                & totalEps1i,totalEps2,totalElf,&
                                & totalRefractIndex,totalExtnctCoeff,&
                                & totalReflectivity,totalAbsorpCoeff)
         use O_Kinds
         use O_Constants
         integer :: length
         real (kind=double), dimension (:,:) :: eps1,eps1i
         real (kind=double), dimension (:,:) :: refractIndex,extnctCoeff
         real (kind=double), dimension (:,:) :: reflectivity,absorpCoeff
         real (kind=double), dimension (:)   :: totalEps1,totalEps1i
         real (kind=double), dimension (:)   :: totalEps2,totalElf
         real (kind=double), dimension (:)   :: totalRefractIndex
         real (kind=double), dimension (:)   :: totalExtnctCoeff
         real (kind=double), dimension (:)   :: totalReflectivity
         real (kind=double), dimension (:)   :: totalAbsorpCoeff
      end subroutine averageFunctions
      subroutine printData(length,energy,totalEps1,totalEps1i,totalEps2,&
                         & totalElf,totalRefractIndex,totalExtnctCoeff,&
                         & totalReflectivity,totalAbsorpCoeff,eps1,eps1i,elf,&
                         & refractIndex,extnctCoeff,reflectivity,absorpCoeff)
         use O_Kinds
         use O_Constants
         integer :: length
         real (kind=double), dimension (:)   :: energy
         real (kind=double), dimension (:)   :: totalEps1,totalEps1i
         real (kind=double), dimension (:)   :: totalEps2,totalElf
         real (kind=double), dimension (:)   :: totalRefractIndex
         real (kind=double), dimension (:)   :: totalExtnctCoeff
         real (kind=double), dimension (:)   :: totalReflectivity
         real (kind=double), dimension (:)   :: totalAbsorpCoeff
         real (kind=double), dimension (:,:) :: eps1,eps1i,elf
         real (kind=double), dimension (:,:) :: refractIndex,extnctCoeff
         real (kind=double), dimension (:,:) :: reflectivity,absorpCoeff
      end subroutine printData
      integer function countNumericFields (textLine)
         character*(*) :: textLine
      end function countNumericFields
   end interface

   ! Get command line parameter for the number of lines. Note that the
   !   caller passes a count of LINES in the eps2 file, header included,
   !   not a count of data points. Both callers produce it that way:
   !   imago.py runs "wc -l" on fort.50 for the total spectrum, and
   !   processPOPTC.py counts the lines of the fort.450 that makePDOS.py
   !   wrote for one partial pair. Each of those files carries exactly
   !   one header line, so one line is subtracted just below and every
   !   later use of "length" means the number of data points.
   call getarg(1,buffer)
   read (buffer,*) length

   ! Adjust the length to subtract the header
   length = length - 1

   ! Determine the number of values from the number of lines and the
   !   granularity we desire. This is an upper bound rather than the
   !   working size: makeFineGrainEnergy interpolates between adjacent
   !   pairs of data points, and there are length-1 such gaps, so only
   !   the first grain*(length-1) entries are ever filled or read. That
   !   count is written throughout as numValues-grain, which is where
   !   that otherwise puzzling expression comes from. The final grain
   !   entries of the fine arrays stay untouched.
   numValues=length*grain

   ! Determine if this is a spin polarized calculation or not.
   call getarg(2,buffer)
   read (buffer,*) spin

   ! Read in a flag indicating that POPTC input file numbers should be used.
   call getarg(3,buffer)
   read (buffer,*) usePOPTCFiles

   ! Read in the POPTC factor for this KKC.
   call getarg(4,buffer)
   read (buffer,*) pOptcFactor

   if (usePOPTCFiles == 0) then
      if (spin == 1) then
         open (unit=50, file='fort.50', status="old",form="formatted") ! Eps2
         open (unit=100,file='fort.100',status="new",form="formatted") ! Total
         open (unit=110,file='fort.110',status="new",form="formatted") ! Eps1
         open (unit=120,file='fort.120',status="new",form="formatted") ! ELF 
         open (unit=130,file='fort.130',status="new",form="formatted") ! Refrct
         open (unit=140,file='fort.140',status="new",form="formatted") ! Extnct
         open (unit=150,file='fort.150',status="new",form="formatted") ! Eps1i
         open (unit=160,file='fort.160',status="new",form="formatted") ! Reflct
         open (unit=170,file='fort.170',status="new",form="formatted") ! Absorp
      else
         open (unit=50, file='fort.51', status="old",form="formatted") ! Eps2
         open (unit=100,file='fort.101',status="new",form="formatted") ! Total
         open (unit=110,file='fort.111',status="new",form="formatted") ! Eps1
         open (unit=120,file='fort.121',status="new",form="formatted") ! ELF 
         open (unit=130,file='fort.131',status="new",form="formatted") ! Refrct
         open (unit=140,file='fort.141',status="new",form="formatted") ! Extnct
         open (unit=150,file='fort.151',status="new",form="formatted") ! Eps1i
         open (unit=160,file='fort.161',status="new",form="formatted") ! Reflct
         open (unit=170,file='fort.171',status="new",form="formatted") ! Absorp
      endif
   else
      ! The processPOPTC script automatically calls imagoKkc multiple times
      !   for all necessary spectra and it uses the same input in all cases.
      !   (I.e., processPOPTC does not produce different input files for
      !   imagoKkc for spin up vs. spin down.)
      open (unit=50, file='fort.450',status="old",form="formatted") ! Eps2
      open (unit=100,file='fort.500',status="new",form="formatted") ! Total
      open (unit=110,file='fort.510',status="new",form="formatted") ! Eps1
      open (unit=120,file='fort.520',status="new",form="formatted") ! ELF 
      open (unit=130,file='fort.530',status="new",form="formatted") ! Refrct
      open (unit=140,file='fort.540',status="new",form="formatted") ! Extnct
      open (unit=150,file='fort.550',status="new",form="formatted") ! Eps1i
      open (unit=160,file='fort.560',status="new",form="formatted") ! Reflct
      open (unit=170,file='fort.570',status="new",form="formatted") ! Absorp
   endif

   ! Discover the record width before anything is sized by it. Counting
   !   the numeric fields on the first DATA line works for both kinds of
   !   input this program is given: the total spectrum files, whose
   !   header names its columns but declares no count, and the
   !   per-partial files that makePDOS.py extracts, whose header is
   !   reduced to a single line by the time it arrives here. Reading the
   !   count from the data is therefore the one method that serves both,
   !   and it cannot disagree with the data it describes.
   read (50,fmt="(a)") probeLine
   read (50,fmt="(a)") probeLine
   rewind (50)

   numFields = countNumericFields(probeLine)

   ! Drop the energy column and the leading isotropic column; whatever
   !   remains is direction resolved.
   numDirCols = numFields - 2

   if ((numDirCols == 0) .or. (numDirCols == 3) .or. &
         & (numDirCols == 6)) then
      numWorkCols = max(numDirCols,1)
      numScalarCols = min(3,numWorkCols)
   else
      write (6,*) 'imagoKKc: this input carries ',numFields, &
            & ' fields per line, giving ',numDirCols, &
            & ' direction resolved columns.'
      write (6,*) 'Only 0 (the direction average alone), 3 (x, y and z)'
      write (6,*) 'and 6 (the symmetric tensor) are meaningful. A count'
      write (6,*) 'outside those means the file was not written by an'
      write (6,*) 'optical writer, or was truncated.'
      stop 'imagoKKc: unsupported number of direction columns'
   endif

   allocate (energy(length))
   allocate (fine_energy(numValues))
   allocate (totalEps1(length))
   allocate (totalEps1i(length))
   allocate (totalEps2(length))
   allocate (totalElf(length))
   allocate (totalRefractIndex(length))
   allocate (totalExtnctCoeff(length))
   allocate (totalReflectivity(length))
   allocate (totalAbsorpCoeff(length))
   ! Epsilon-2 and epsilon-1 carry every component the input holds.
   allocate (eps1(numWorkCols,length))
   allocate (eps1i(numWorkCols,length))
   allocate (eps2(numWorkCols,length))
   allocate (fine_eps2(numWorkCols,numValues))

   ! The scalar-derived five carry the diagonal entries only.
   allocate (elf(numScalarCols,length))
   allocate (refractIndex(numScalarCols,length))
   allocate (extnctCoeff(numScalarCols,length))
   allocate (reflectivity(numScalarCols,length))
   allocate (absorpCoeff(numScalarCols,length))

   ! Read the energy and epsilon2 data (output of the optc program).
   call readOPTCData(length, energy, totalEps2, eps2)

   ! Determine the average energy difference for the coarse grain input.
   call getCoarseEnergyDiff(length,energy,coarseEnergyDiff)

   ! Generate a fine grained version of the energy data.
   call makeFineGrainEnergy(length,grain,energy,fine_energy)

   ! Determine the average energy difference for the fine grain energies.
   call getFineEnergyDiff(numValues,grain,fine_energy,fineEnergyDiff)

   ! Generate a fine grained version of the eps2 data. This will be used to
   !  improve the accuracy of the integration as well as to avoid the poles in
   !  the integral.
   call makeFineEps2(length,grain,energy,fine_energy,eps2,fine_eps2)

   ! Do Kramers Kronig Conversion to get Epsilon1 from Epsilon2, and Epsilon1i
   !   from Epsiton2.
   call kramersKronig(length,grain,numValues,fineEnergyDiff,&
                      &energy,fine_energy,eps1,eps1i,fine_eps2,&
                      &pOptcFactor)

   ! Calculate the refractive index, extinction coefficient, reflectivity,
   !   and the absorption coefficient from epsilon1.
   call get_n_k_R_a(length,energy,eps1,eps2,refractIndex,extnctCoeff,&
         & reflectivity,absorpCoeff)

   ! Calculate the Energy Loss Function from Epsilon1 and Epsilon2
   call getELF(length,eps1,eps2,elf)

   ! Average Epsilon1, Epsilon1i, Epsilon2, and the ELF over 3 dimensions
   call averageFunctions(length,eps1,eps1i,refractIndex,extnctCoeff,&
                       & reflectivity,absorpCoeff,totalEps1,totalEps1i,&
                       & totalEps2,totalElf,totalRefractIndex,&
                       & totalExtnctCoeff,totalReflectivity,totalAbsorpCoeff)

   ! Print a data file of Epsilon1, Epsilon1i, Epsilon2 and the ELF
   call printData(length,energy,totalEps1,totalEps1i,totalEps2,totalElf,&
                & totalRefractIndex,totalExtnctCoeff,totalReflectivity,&
                & totalAbsorpCoeff,eps1,eps1i,elf,refractIndex,extnctCoeff,&
                & reflectivity,absorpCoeff)

   deallocate (energy)
   deallocate (fine_energy)
   deallocate (totalEps1)
   deallocate (totalEps1i)
   deallocate (totalEps2)
   deallocate (totalElf)
   deallocate (totalRefractIndex)
   deallocate (totalExtnctCoeff)
   deallocate (totalReflectivity)
   deallocate (totalAbsorpCoeff)
   deallocate (eps1)
   deallocate (eps1i)
   deallocate (eps2)
   deallocate (fine_eps2)
   deallocate (elf)
   deallocate (refractIndex)
   deallocate (extnctCoeff)
   deallocate (reflectivity)
   deallocate (absorpCoeff)

   write(0,*) 'Kramers-Kronig conversion relation applied.'

end program imagoKkc

subroutine readOPTCData (length, energy, totalEps2, eps2)

   use O_Kinds
   use O_Constants

   implicit none

   integer :: length
   real (kind=double), dimension (:) :: energy
   real (kind=double), dimension (:) :: totalEps2
   real (kind=double), dimension (:,:) :: eps2

   integer :: i

   ! Read past the header
   read (50,*)

   ! Read the energy and epsilon2 data (output of the optc program).
   !
   ! At direction code 0 the record holds the energy and one value, and
   !   that value IS the working column, so it is copied into the single
   !   slot rather than read separately. Reading eps2(:,i) as well would
   !   run the list-directed read on into the NEXT record, quietly
   !   shifting every subsequent line by one -- a misalignment that
   !   produces a complete, plausible spectrum rather than an error.
   do i=1,length
      if (size(eps2,1) == 1) then
         read (50,*) energy(i),totalEps2(i)
         eps2(1,i) = totalEps2(i)
      else
         read (50,*) energy(i),totalEps2(i),eps2(:,i)
      endif
   enddo

end subroutine readOPTCData

subroutine getCoarseEnergyDiff(length, energy, coarseEnergyDiff)

   use O_Kinds
   use O_Constants

   implicit none

   integer :: length
   real (kind=double), dimension (:) :: energy
   real (kind=double) :: coarseEnergyDiff

   integer :: i

   !Determine the average energy difference for the coarse grain input.
   coarseEnergyDiff=0.0_double
   do i=1,length-1
      coarseEnergyDiff=coarseEnergyDiff+(energy(i+1)-energy(i))
   enddo
   coarseEnergyDiff=coarseEnergyDiff/(length-1)

end subroutine getCoarseEnergyDiff

subroutine makeFineGrainEnergy(length,grain,energy,fine_energy)

   use O_Kinds
   use O_Constants

   implicit none

   integer :: length, grain
   real (kind=double), dimension (:) :: energy, fine_energy

   integer :: i,j,k

   ! Generate a fine grained version of the energy data.
   do i=1,length-1
      do j=1,grain
         k=int((i-1)*grain+j)
         fine_energy(k) = energy(i)+(energy(i+1)-energy(i))*(j-1)/grain
      enddo
   enddo

end subroutine makeFineGrainEnergy

subroutine getFineEnergyDiff(numValues,grain,fine_energy,fineEnergyDiff)

   use O_Kinds
   use O_Constants

   implicit none

   integer :: numValues, grain
   real (kind=double), dimension (:) :: fine_energy
   real (kind=double) :: fineEnergyDiff

   integer :: i

   fineEnergyDiff=0.0_double
   do i=1,numValues-grain-1
      fineEnergyDiff=fineEnergyDiff+fine_energy(i+1)-fine_energy(i)
   enddo
   fineEnergyDiff=fineEnergyDiff/(numValues-grain-1)

end subroutine getFineEnergyDiff

subroutine makeFineEps2(length,grain,energy,fine_energy,eps2,fine_eps2)

   use O_Kinds
   use O_Constants

   implicit none

   integer :: length,grain
   real (kind=double), dimension (:) :: energy,fine_energy
   real (kind=double), dimension (:,:) :: eps2,fine_eps2

   integer :: i,j,k
   real (kind=double), allocatable, dimension(:,:) :: diffRatio

   ! Sized from the array it differentiates rather than from a constant,
   !   so the width follows the caller's direction code automatically.
   allocate(diffRatio(size(eps2,1),length))

   do i=1,length-1
      diffRatio(:,i)=(eps2(:,i+1)-eps2(:,i))/(energy(i+1)-energy(i))
      do j=1,grain
         k=int((i-1)*grain+j)
         fine_eps2(:,k) = eps2(:,i)+diffRatio(:,i)*(fine_energy(k)-energy(i))
      enddo
   enddo

   deallocate(diffRatio)
end subroutine makeFineEps2

subroutine kramersKronig(length,grain,numValues,fineEnergyDiff,&
                        &energy,fine_energy,eps1,eps1i,fine_eps2,&
                        &pOptcFactor)

   use O_Constants
   use O_Kinds

   implicit none

   integer :: length,grain,numValues
   real (kind=double) :: fineEnergyDiff,pOptcFactor
   real (kind=double), dimension (:) :: energy,fine_energy
   real (kind=double), dimension (:,:) :: eps1,eps1i,fine_eps2

   integer :: i,j
   real (kind=double) :: fine_e,original_e
   real (kind=double), allocatable, dimension (:) :: evenSum,oddSum
   real (kind=double), allocatable, dimension (:) :: evenSumi,oddSumi
   real (kind=double), allocatable, dimension (:) :: vacuumTerm
   real (kind=double), allocatable, dimension (:,:) :: integrand,integrandi

   real (kind=double) :: multFactor

   ! This constant carries TWO unrelated factors that happen to be
   !   applied at the same moment, and reading it as one number is
   !   misleading enough to be worth spelling out. The real part of the
   !   dielectric function is recovered from the imaginary part by the
   !   Kramers-Kronig relation, evaluated here as a composite Simpson
   !   integration over the fine-grained energy scale:
   !
   !! eps1(w) = 1 + (2/pi) * P Int[ eps2(w') w' / (w'^2 - w^2) dw' ]
   !!
   !!           |<-- KK prefactor        |<-- Simpson's rule
   !!           |    2/pi                |    (h/3)*[f_0 + 4*odd
   !!           |    one per Cartesian   |     + 2*even + f_n]
   !!           |    component           |
   !
   !   The 2/pi is the Kramers-Kronig prefactor and belongs to each
   !   Cartesian component separately. The 1/3 is Simpson's third: the
   !   quadrature expression below multiplies by fineEnergyDiff, which
   !   is the plain fine-grid spacing h and not h/3, so the missing
   !   third is supplied from here instead.
   !
   !   The 1/3 is therefore NOT a directional average. The only average
   !   over x, y and z is the division by three in averageFunctions,
   !   which forms totalEps1 from the three component eps1 values. A
   !   reader who takes this 1/3 for that average will conclude the
   !   code divides by three twice and that every eps1, ELF, n, k, R
   !   and alpha Imago prints is too small by a factor of three. It
   !   does not, and they are not.
   multFactor = 2.0_double/3.0_double/pi

   allocate (evenSum    (size(eps1,1)))
   allocate (oddSum     (size(eps1,1)))
   allocate (integrand  (size(eps1,1),numValues))
   allocate (evenSumi   (size(eps1,1)))
   allocate (oddSumi    (size(eps1,1)))
   allocate (integrandi (size(eps1,1),numValues))

   ! The additive term of the dispersion relation,
   !
   !   eps1_ij(w) = delta_ij + (2/pi) P Int[ ... ]
   !
   !   is a KRONECKER DELTA: it is the vacuum contribution, present on
   !   the diagonal and absent off it. So it goes on the first three
   !   entries and nowhere else. Adding it to every column would leave
   !   each off-diagonal epsilon-1 sitting at exactly one where it
   !   should sit at zero -- a constant offset, uniform across the
   !   spectrum, that looks like a plausible baseline rather than an
   !   error.
   !
   ! At direction codes 0 and 1 every column IS diagonal, so this is
   !   the unconditional addition it has always been.
   allocate (vacuumTerm (size(eps1,1)))
   vacuumTerm(:) = 0.0_double
   vacuumTerm(1:min(3,size(eps1,1))) = pOptcFactor

   do i=1,length
      oddSum(:)    = 0.0_double
      evenSum(:)   = 0.0_double
      oddSumi(:)   = 0.0_double
      evenSumi(:)  = 0.0_double

      ! The first and last nodes each carry weight one in the sum below,
      !   and each is skipped when it lands on the pole at energy(i).
      !   Clear them before the guards so that a skipped node adds
      !   nothing, rather than adding whatever the array happened to
      !   hold. On the first pass through i it holds nothing at all,
      !   because an allocate does not initialize.
      integrand(:,1)                = 0.0_double
      integrandi(:,1)               = 0.0_double
      integrand(:,numValues-grain)  = 0.0_double
      integrandi(:,numValues-grain) = 0.0_double

      ! Do the first point explicitly to save computation in the loop.
      !   The pole sits at energy(i), the energy this eps1 value is
      !   being computed for, and not at either end of the range.
      fine_e=fine_energy(1)
      original_e=energy(i)
      if (abs(fine_e - original_e) .gt. 0.00001_double) then
         integrand(:,1)=fine_eps2(:,1)*fine_e/ &
               & (original_e + fine_e)/(fine_e - original_e)
         integrandi(:,1)=fine_eps2(:,1)*fine_e/ &
               & (fine_e*fine_e + original_e*original_e)
      endif

      ! Loop through the middle part
      do j=2,numValues-grain-1
         !Exclude all points from energy(i) in the summation to avoid
         !  poles in the integration.
         fine_e=fine_energy(j)
         original_e=energy(i)
         if (abs(fine_e - original_e) .gt. 0.00001_double) then
            integrand(:,j)=fine_eps2(:,j)*fine_e / &
                  & (original_e + fine_e)/(fine_e - original_e)
            integrandi(:,j)=fine_eps2(:,j)*fine_e / &
                  & (fine_e*fine_e + original_e*original_e)
            if ((j/2)*2 .eq. j) then
               evenSum(:)=evenSum(:)+integrand(:,j)
               evenSumi(:)=evenSumi(:)+integrandi(:,j)
            else
               oddSum(:)=oddSum(:)+integrand(:,j)
               oddSumi(:)=oddSumi(:)+integrandi(:,j)
            endif
         endif
      enddo

      ! Do the last point explicitly to save computation in the loop.
      !   The pole is energy(i) here too. Using the top of the range
      !   instead put this node a single fine step from its own pole for
      !   every i, so it evaluated to roughly -eps2(top)/(2h). The h
      !   then cancelled against the h the sum is multiplied by, leaving
      !   a spurious constant -eps2(top)/(3*pi) subtracted from every
      !   eps1 value -- an error that no amount of grid refinement could
      !   reduce, because it did not depend on the grid.
      fine_e=fine_energy(numValues-grain)
      original_e=energy(i)
      if (abs(fine_e - original_e) .gt. 0.00001_double) then
         integrand(:,numValues-grain)=fine_eps2(:,numValues-grain)*fine_e/ &
               & (original_e + fine_e)/(fine_e - original_e)
         integrandi(:,numValues-grain)=fine_eps2(:,numValues-grain)*fine_e/ &
               & (fine_e*fine_e + original_e*original_e)
      endif

      ! Assemble the Simpson sum. The leading fineEnergyDiff is the bare
      !   fine-grid spacing h, so the third that completes Simpson's h/3
      !   arrives with multFactor instead; see the note at its
      !   assignment. The trailing pOptcFactor is the additive 1 of
      !   eps1 = 1 + (2/pi)*Int[...], divided among the partials by the
      !   caller so that a full set of them still sums to that single 1
      !   rather than to one per partial.
      ! Simpson's rule weights the even-numbered nodes by four and the
      !   odd-numbered interior ones by two, counting the first node as
      !   number one. The two weights were previously the other way
      !   round, which left the total weight short and doubled the
      !   error of the whole integration.
      eps1(:,i)=fineEnergyDiff*(integrand(:,int(numValues-grain)) + &
            & integrand(:,1) + 4.0_double * evenSum(:) + 2.0_double * &
            & oddSum(:)) * multFactor + vacuumTerm(:)
      eps1i(:,i)=fineEnergyDiff*(integrandi(:,int(numValues-grain)) + &
            & integrandi(:,1) + 4.0_double * evenSumi(:) + 2.0_double * &
            & oddSumi(:)) * multFactor + vacuumTerm(:)
   enddo

   deallocate (evenSum)
   deallocate (oddSum)
   deallocate (integrand)
   deallocate (evenSumi)
   deallocate (oddSumi)
   deallocate (integrandi)
   deallocate (vacuumTerm)

end subroutine kramersKronig

subroutine get_n_k_R_a(length,energy,eps1,eps2,refractIndex,extnctCoeff,&
      & reflectivity,absorpCoeff)

   use O_Kinds
   use O_Constants

   implicit none

   ! Define passed parameters
   integer :: length
   real (kind=double), dimension (:) :: energy
   real (kind=double), dimension (:,:) :: eps1,eps2
   real (kind=double), dimension (:,:) :: refractIndex,extnctCoeff
   real (kind=double), dimension (:,:) :: reflectivity,absorpCoeff

   ! Define local variables.
   integer :: i,j
   complex (kind=double) :: reflecTemp

   ! Bounded by the DESTINATION arrays rather than by three, and not by
   !   eps1 either. At direction code 0 there is a single column, and
   !   running to three here writes past the end of every one of these
   !   -- which does not fail at the write but corrupts the heap and
   !   aborts later, somewhere with no connection to the cause. At code
   !   2 eps1 is six wide while these are three, so bounding by eps1
   !   would make the same mistake in the other direction.
   do i = 1,length
      do j = 1,size(refractIndex,1)

         ! ~eps = eps1 + i eps2
         ! eps1 = n^2 - k^2  (n=refractive index, k=extinction coeff)
         ! eps2 = 2nk
         ! n = sqrt[(sqrt(eps1^2 + eps2^2)+eps1)/2]
         ! k = sqrt[(sqrt(eps1^2 + eps2^2)-eps1)/2] ; kappa
!         ! R = ((n-1)^2 + k) / ((n+1)^2 + k) ! Reflectivity
         ! R = ((n-1)^2) / ((n+1)^2 + k) ! Reflectivity
         ! alpha = 4 * pi * kappa / lambda = 2 * kappa * omega / c
         ! omega = (energy in eV) / hPlanck -> angular frequency
         refractIndex(j,i) = sqrt((sqrt(eps1(j,i)*eps1(j,i) + &
               & eps2(j,i)*eps2(j,i))+eps1(j,i))/2.0_double)

         extnctCoeff(j,i) = sqrt((sqrt(eps1(j,i)*eps1(j,i) + &
               & eps2(j,i)*eps2(j,i))-eps1(j,i))/2.0_double)

!         reflectivity(j,i) = ((refractIndex(j,i)-1)**2 + extnctCoeff(j,i)) / &
!                         & + ((refractIndex(j,i)+1)**2 + extnctCoeff(j,i))
         reflecTemp = (refractIndex(j,i) + &
               & cmplx(0.0_double,extnctCoeff(j,i),double) - 1.0_double)
         reflectivity(j,i) = real(reflecTemp * conjg(reflecTemp),double)
         reflecTemp = (refractIndex(j,i) + &
               & cmplx(0.0_double,extnctCoeff(j,i),double) + 1.0_double)
         reflectivity(j,i) = reflectivity(j,i) / &
               & real(reflecTemp * conjg(reflecTemp),double)

         absorpCoeff(j,i) = 2.0_double * extnctCoeff(j,i) * energy(i) / &
                         & hPlanck / lightSpeed
      enddo
   enddo
end subroutine get_n_k_R_a

subroutine getELF(length,eps1,eps2,elf)

   use O_Kinds
   use O_Constants

   implicit none

   integer :: length
   real (kind=double), dimension (:,:) :: eps1,eps2,elf

   integer :: i

   ! ELF = Im(-1/eps) = eps2 / (eps1^2 + eps2^2)
   !
   ! Restricted to the DIAGONAL entries. This expression is the loss
   !   function of a scalar dielectric function; the tensor form is
   !   Im(-eps^-1), a matrix inverse, which is not this applied
   !   element-wise. So elf is narrower than eps1 and eps2 whenever the
   !   input carries off-diagonal components (DESIGN 13.7).
   do i = 1,length
      elf(:,i) = eps2(1:size(elf,1),i) &
            & / (eps1(1:size(elf,1),i)**2 + eps2(1:size(elf,1),i)**2)
   enddo


end subroutine getELF

subroutine averageFunctions(length,eps1,eps1i,refractIndex,extnctCoeff,&
                          & reflectivity,absorpCoeff,totalEps1,totalEps1i,&
                          & totalEps2,totalElf,totalRefractIndex,&
                          & totalExtnctCoeff,totalReflectivity,&
                          & totalAbsorpCoeff)

   use O_Kinds
   use O_Constants

   implicit none

   integer :: length
   real (kind=double), dimension (:,:) :: eps1,eps1i
   real (kind=double), dimension (:,:) :: refractIndex,extnctCoeff
   real (kind=double), dimension (:,:) :: reflectivity,absorpCoeff
   real (kind=double), dimension (:) :: totalEps1,totalEps1i,totalEps2,&
         & totalElf,totalRefractIndex,totalExtnctCoeff,totalReflectivity,&
         & totalAbsorpCoeff

   integer :: i
   integer :: numDiagonal

   ! The divisor is the number of DIAGONAL columns carried, not the
   !   literal three. At direction code 1 that IS three and every result
   !   is unchanged. At direction code 0 there is a single column which
   !   is already the direction average, so dividing by one returns it
   !   untouched -- dividing by three there would shrink every derived
   !   spectrum to a third of its value, silently.
   !
   ! Note that totalElf is NOT an average of the per-column values. It
   !   is built from the total eps1 and eps2, because the energy loss
   !   function is a non-linear function of the dielectric function and
   !   the average of Im(-1/eps) is not Im(-1/eps_avg).
   ! The isotropic value is one third of the TRACE, so the average runs
   !   over the DIAGONAL entries only. At direction code 2 epsilon-1 is
   !   six columns wide while this average must still take three;
   !   averaging over all six would fold the off-diagonal entries into
   !   the trace and produce a number that is finite, plausible and
   !   wrong. `refractIndex` is the diagonal width by construction, so
   !   it is what the bound is taken from.
   numDiagonal = size(refractIndex,1)

   do i = 1,length
      totalEps1(i)         = sum(eps1(1:numDiagonal,i)) &
            & /real(numDiagonal,double)
      totalEps1i(i)        = sum(eps1i(1:numDiagonal,i)) &
            & /real(numDiagonal,double)
      totalElf(i)          = totalEps2(i)/(totalEps1(i)**2+totalEps2(i)**2)
      totalRefractIndex(i) = sum(refractIndex(:,i)) &
            & /real(numDiagonal,double)
      totalExtnctCoeff(i)  = sum(extnctCoeff(:,i)) &
            & /real(numDiagonal,double)
      totalReflectivity(i) = sum(reflectivity(:,i)) &
            & /real(numDiagonal,double)
      totalAbsorpCoeff(i)  = sum(absorpCoeff(:,i)) &
            & /real(numDiagonal,double)
   enddo

end subroutine averageFunctions

subroutine printData(length,energy,totalEps1,totalEps1i,totalEps2,totalElf,&
      & totalRefractIndex,totalExtnctCoeff,totalReflectivity,&
      & totalAbsorpCoeff,eps1,eps1i,elf,refractIndex,extnctCoeff,&
      & reflectivity,absorpCoeff)

   use O_Kinds
   use O_Constants

   implicit none

   integer :: length
   real (kind=double), dimension (:)   :: energy,totalEps1,totalEps1i
   real (kind=double), dimension (:)   :: totalEps2,totalElf
   real (kind=double), dimension (:)   :: totalRefractIndex
   real (kind=double), dimension (:)   :: totalExtnctCoeff
   real (kind=double), dimension (:)   :: totalReflectivity
   real (kind=double), dimension (:)   :: totalAbsorpCoeff
   real (kind=double), dimension (:,:) :: eps1,eps1i,elf
   real (kind=double), dimension (:,:) :: refractIndex,extnctCoeff
   real (kind=double), dimension (:,:) :: reflectivity,absorpCoeff

   integer :: i

   write (100,*) "Energy   Epsilon1   Epsilon2   ELF"
   do i = 1,length
      write (100,701) energy(i),totalEps1(i),totalEps2(i),totalELF(i)
   enddo

   ! Each of these writes a header naming its columns and then one
   !   record per energy point, at whatever width the input carried.
   !   The per-direction columns disappear entirely at direction code 0,
   !   where the single value is the direction average and repeating it
   !   under an "x" heading would invite a reader to treat one number as
   !   three independent ones.
   call writeSpectrum (110,"Eps1",length,size(eps1,1),energy, &
         & totalEps1,eps1)
   call writeSpectrum (120,"ELF",length,size(elf,1),energy, &
         & totalELF,elf)
   call writeSpectrum (130,"n",length,size(refractIndex,1),energy, &
         & totalRefractIndex,refractIndex)
   call writeSpectrum (140,"k",length,size(extnctCoeff,1),energy, &
         & totalExtnctCoeff,extnctCoeff)
   call writeSpectrum (150,"Eps1i",length,size(eps1i,1),energy, &
         & totalEps1i,eps1i)
   call writeSpectrum (160,"Reflectivity",length, &
         & size(reflectivity,1),energy,totalReflectivity,reflectivity)

   call writeSpectrum (170,"Alpha",length,size(absorpCoeff,1),energy, &
         & totalAbsorpCoeff,absorpCoeff)

   701 format (4e15.7)

end subroutine printData


! Write one derived spectrum: a header naming its columns, then one
!   record per energy point.
!
! The width follows the array rather than a literal, because Imago's
!   optical writers choose their record width from the direction code
!   of DESIGN 13.7. At code 1 this produces exactly the five field
!   record and the "total x y z" heading that this program has always
!   written. At code 0 there is one working column which IS the
!   direction average, so only the energy and that average are written.
! The arrays are declared with EXPLICIT shape rather than assumed
!   shape, so that no explicit interface is needed to call this. That
!   matters here because the caller is itself an external subroutine
!   and cannot see the interface block the main program carries.
subroutine writeSpectrum (outputUnit,quantityName,length,numCols, &
      & energy,totalValues,columnValues)

   use O_Kinds

   implicit none

   integer :: outputUnit
   character*(*) :: quantityName
   integer :: length
   integer :: numCols
   real (kind=double), dimension (length) :: energy, totalValues
   real (kind=double), dimension (numCols,length) :: columnValues

   integer :: i, col
   character*20 :: recordFormat
   character*200 :: header
   ! Three columns are named for the Cartesian axes; six are named for
   !   the entries of the symmetric tensor, in the order PSEUDOCODE
   !   21.7 declares and the engine's writers fill.
   character*3, dimension (6), parameter :: axisTag = &
         & (/ "x  ", "y  ", "z  ", "   ", "   ", "   " /)
   character*3, dimension (6), parameter :: tensorTag = &
         & (/ "xx ", "yy ", "zz ", "xy ", "xz ", "yz " /)

   header = " Energy   total"//trim(quantityName)
   if (numCols > 1) then
      do col = 1, numCols
         if (numCols == 6) then
            header = trim(header)//"   "//trim(tensorTag(col))// &
                  & trim(quantityName)
         else
            header = trim(header)//"   "//trim(axisTag(col))// &
                  & trim(quantityName)
         endif
      enddo
   endif
   write (outputUnit,fmt="(a)") trim(header)

   if (numCols == 1) then
      write (recordFormat,fmt="(a)") "(2e15.7)"
      do i = 1,length
         write (outputUnit,fmt=recordFormat) energy(i),totalValues(i)
      enddo
   else
      write (recordFormat,fmt="(a,i2,a)") "(",numCols+2,"e15.7)"
      do i = 1,length
         write (outputUnit,fmt=recordFormat) energy(i),totalValues(i),&
               & columnValues(:,i)
      enddo
   endif

end subroutine writeSpectrum


! Count the numeric fields on one line of an input record.
!
! Used to discover how many direction resolved columns the input
!   carries, since the total spectrum files name their columns in a
!   header but declare no count. Counting transitions from blank to
!   non-blank is enough: every field on these lines is a number written
!   in a fixed width numeric format, so there are no quoted strings or
!   embedded spaces to confuse it.
integer function countNumericFields (textLine)

   implicit none

   character*(*) :: textLine

   integer :: i
   logical :: insideField

   countNumericFields = 0
   insideField = .false.

   do i = 1, len_trim(textLine)
      if (textLine(i:i) == ' ') then
         insideField = .false.
      else
         if (.not. insideField) then
            countNumericFields = countNumericFields + 1
         endif
         insideField = .true.
      endif
   enddo

end function countNumericFields

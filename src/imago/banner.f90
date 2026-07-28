!! SPDX-License-Identifier: ECL-2.0
!! Copyright (c) 2026 Paul Rulis

!****************************************************************************
!
! The identity block: the logo, the wordmark, and Imago's own citation,
!   printed at the head of the log.
!
! A license header and a CITATION.cff file are addressed at whoever
!   redistributes or packages the code. Neither reaches the person writing
!   the paper, who has results in front of them and is composing a methods
!   section, and who will cite what they can see. Printing the citation
!   into the output they are already reading is the one attribution
!   mechanism this project has that arrives at the moment the decision is
!   actually made, which is why LAMMPS, VASP, and Quantum ESPRESSO all do
!   the same thing.
!
! The artwork and the citation text are not compiled in. They live in
!   banner.txt, installed alongside elements.dat and found through the
!   IMAGO_DATA environment variable, so updating the DOI at release time is
!   an edit to one text file rather than to a Fortran literal.
!
! This module deliberately depends on nothing but O_Verboseness. The
!   citations for the methods a run exercised are a separate module,
!   O_MethodCitations, because those predicates must read engine state and
!   this block must print before any of that state exists. See PSEUDOCODE
!   17.1: a single module holding both halves closes a module cycle through
!   O_CommandLine and O_KPoints that no build order resolves.
!
!****************************************************************************
module O_Banner

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains


! Print the identity block, if this run wants it.
!
! The gate is the whole of what this routine adds, and it is the reason
!   initVerboseness must already have run. Reversed, the mask is still zero,
!   isVerbose returns false, and this returns without printing or
!   complaining: the run succeeds and the banner is simply absent, with
!   nothing in the compiler or the output pointing at why.
subroutine printIdentityBlock

   ! Import necessary modules.
   use O_Verboseness, only: isVerbose, VERB_BANNER

   ! Make sure that no funny variables are defined.
   implicit none

   if (.not. isVerbose(VERB_BANNER)) return

   call echoBannerFile

   ! Separate the block from the first timestamp rule the same way every
   !   later pair of blocks is separated.
   write (20,*)
   call flush (20)

end subroutine printIdentityBlock


! Copy banner.txt to the log, line for line.
!
! The formats are the point of this routine. Reading under fmt='(a)'
!   preserves leading blanks, which carry the kerning of the wordmark and
!   the centring of the butterfly, and writing under fmt='(a)' starts at
!   column one. A list directed write (20,*) would insert a leading blank
!   on every line and shift the whole figure one column right of the
!   character(len=51) rules that O_TimeStamps prints directly beneath it.
!   The result would look almost right, which is worse than looking wrong.
!
! For the same reason trim is used and adjustl is not. Trailing blanks are
!   insignificant and trim removes only those; adjustl would left justify
!   every line and destroy the artwork.
subroutine echoBannerFile

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define the local variables used in this subroutine.
   character(len=100) :: dataDirectory  ! Where the databases were installed.
   character(len=100) :: bannerFileName ! Full path to banner.txt.

   ! The buffer is 132 characters and not 51. The artwork obeys the 51
   !   column limit that matches the opLabels of O_TimeStamps, its widest
   !   line being 50, but the citation lines beneath it neither do nor need
   !   to: the longest today is the 64 character line naming the repository,
   !   and the DOI line will be longer still. A 51 character buffer would
   !   silently truncate the very text this block exists to deliver.
   character(len=132) :: lineBuffer

   ! Chosen to avoid unit 313 (elementData), unit 9 (potential), and unit
   !   20 (the log). Any free number would do; the point is that it was
   !   checked.
   integer, parameter :: bannerUnit = 314

   integer :: envStatus  ! Whether IMAGO_DATA could be read.
   integer :: openStatus ! Whether banner.txt could be opened.
   integer :: readStatus ! Non-zero at end of file.

   ! A missing IMAGO_DATA cannot happen here in practice: Imago calls
   !   initElementData before parseCommandLine, and that routine stops the
   !   run outright when the variable is unreadable, so a run that reaches
   !   the banner has already proved it good. The test exists so this
   !   routine is honest on its own terms instead of relying on code two
   !   files away, and it is silent because anything it could say has
   !   already been said by the code that actually failed.
   call get_environment_variable (NAME="IMAGO_DATA", &
         & VALUE=dataDirectory, STATUS=envStatus)
   if (envStatus /= 0) return

   bannerFileName = trim(dataDirectory)//"/banner.txt"

   open (unit=bannerUnit, file=bannerFileName, form='formatted', &
         & status='old', IOSTAT=openStatus)
   if (openStatus /= 0) then
      write (20,*) "Could not open ",trim(bannerFileName)
      write (20,*) "Continuing without the identity block."
      return
   endif

   do
      read (bannerUnit, fmt='(a)', IOSTAT=readStatus) lineBuffer
      if (readStatus /= 0) exit
      write (20, fmt='(a)') trim(lineBuffer)
   enddo

   close (bannerUnit)

end subroutine echoBannerFile


end module O_Banner

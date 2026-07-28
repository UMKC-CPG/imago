!! SPDX-License-Identifier: ECL-2.0
!! Copyright (c) 2026 Paul Rulis

!****************************************************************************
!
! The methods block: citations for the methods a run actually exercised,
!   printed at the end of the run.
!
! This cannot go at the head of the log beside the identity block, because
!   at startup the program does not yet know whether the tetrahedron
!   integration will be used. Printing early would mean either listing
!   everything Imago could conceivably do, which trains the reader to skip
!   it, or guessing from the job code, which is wrong whenever a branch is
!   not taken.
!
! The structure is a registry pairing each reference with a predicate
!   answering "did this run use it," and only the references whose
!   predicate is true are printed. Keeping the two beside each other is
!   what makes adding a method to Imago and adding its citation one task
!   rather than two, and makes a reference that no predicate can ever
!   select visibly dead.
!
! The predicates read state the engine already holds, so nothing new has to
!   be tracked to support this. That is also why this module is separate
!   from O_Banner: reading engine state is exactly what the identity block
!   may not do, since it prints before that state exists.
!
! Note that the registry is narrower than the reference list in the design
!   documents, and must be. The UFF parameters, and the force field
!   references beside them, belong to the Python path that builds LAMMPS
!   inputs; they appear nowhere in this engine, so an entry here could
!   never be selected. If that path should announce its own citations, it
!   announces them itself.
!
!****************************************************************************
module O_MethodCitations

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! The number of references in the registry, and the greatest number of
   !   lines any one of them occupies. Adding a reference means raising the
   !   first, declaring one parameter, filling in its text, and adding one
   !   case to methodWasUsed.
   integer, parameter :: numMethodRefs = 2
   integer, parameter :: maxRefLines   = 4

   ! The registry index of each reference. Callers and the predicate use
   !   these names rather than bare numbers, for the same reason the
   !   verboseness categories do.
   integer, parameter :: METHOD_MONKHORST_PACK = 1
   integer, parameter :: METHOD_BLOECHL_LAT    = 2

   ! The citation text itself, and how many lines of it each reference uses.
   character(len=76), dimension(maxRefLines,numMethodRefs) :: refText
   integer, dimension(numMethodRefs) :: refNumLines

   ! The rule printed above the block, matching the width of the operation
   !   labels in O_TimeStamps so the two line up.
   character(len=51), parameter :: sectionRule = &
         & '***************************************************'

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains


! Fill in the citation text for every reference in the registry. These
!   restate what the design documents already carry, which is duplication,
!   accepted because the engine cannot read a design document. The two must
!   be edited together.
subroutine initMethodCitations

   ! Make sure that no funny variables are defined.
   implicit none

   ! The uniform reciprocal space mesh, its shift, and the symmetry folding
   !   to the irreducible zone.
   refNumLines(METHOD_MONKHORST_PACK) = 3
   refText(1,METHOD_MONKHORST_PACK) = &
         & "  H. J. Monkhorst, J. D. Pack, ""Special points for"
   refText(2,METHOD_MONKHORST_PACK) = &
         & "  Brillouin-zone integrations,"" Phys. Rev. B 13, 5188"
   refText(3,METHOD_MONKHORST_PACK) = &
         & "  (1976).  DOI: 10.1103/PhysRevB.13.5188"

   ! The linear analytic tetrahedron method used for the density of states
   !   and for integrated properties.
   refNumLines(METHOD_BLOECHL_LAT) = 4
   refText(1,METHOD_BLOECHL_LAT) = &
         & "  P. E. Bloechl, O. Jepsen, O. K. Andersen, ""Improved"
   refText(2,METHOD_BLOECHL_LAT) = &
         & "  tetrahedron method for Brillouin-zone integrations,"""
   refText(3,METHOD_BLOECHL_LAT) = &
         & "  Phys. Rev. B 49, 16223 (1994)."
   refText(4,METHOD_BLOECHL_LAT) = &
         & "  DOI: 10.1103/PhysRevB.49.16223"

end subroutine initMethodCitations


! Answer whether this run exercised the given method.
function methodWasUsed (methodRef)

   ! Import necessary modules.
   use O_KPoints, only: kPointStyleCode, kPointIntgCode

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define the function result.
   logical :: methodWasUsed

   ! Define the passed dummy variables.
   integer, intent(in) :: methodRef ! A registry index, e.g. METHOD_*

   select case (methodRef)

   case (METHOD_MONKHORST_PACK)
      ! Style code 1 is an explicit mesh plus a shift and style code 2 is a
      !   minimum density; both build a Monkhorst-Pack mesh. Style code 0
      !   is a bare list of k-points the user supplied, which is not
      !   necessarily one, so it earns no citation.
      methodWasUsed = ((kPointStyleCode == 1) .or. &
                     & (kPointStyleCode == 2))

   case (METHOD_BLOECHL_LAT)
      ! Integration code 1 selects the linear analytic tetrahedron method;
      !   code 0 is the histogram method, which cites nothing.
      methodWasUsed = (kPointIntgCode == 1)

   case default
      methodWasUsed = .false.

   end select

end function methodWasUsed


! Print the citations for whichever methods this run exercised.
!
! This block is not governed by IMAGO_VERBOSENESS. It is a handful of lines,
!   it is the part a reader is meant to copy into a paper, and a flight
!   that has turned the artwork off has no reason to discard the citations
!   its runs earned.
subroutine printMethodsBlock

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define the local variables used in this subroutine.
   integer :: methodRef     ! Loop index over the registry.
   integer :: line          ! Loop index over one reference's text.
   logical :: anyMethodUsed ! Whether anything at all will be printed.
   logical :: logIsOpen     ! Whether unit 20 is still connected.

   call initMethodCitations

   ! The log must still be open. Writing to a closed unit does not fail:
   !   Fortran reconnects it to fort.20 and truncates, which would take the
   !   run's entire log with it and leave behind only the citation that
   !   destroyed it. Complain where it can still be seen, and write
   !   nothing. This is deliberately not a repair. Closing the log once at
   !   the end of Imago is what keeps the unit open, and if some later
   !   change closes it early again then that change should cost the
   !   citations and announce itself, rather than silently costing the run
   !   everything the way it did the first time.
   inquire (unit=20, opened=logIsOpen)
   if (.not. logIsOpen) then
      write (6,*) "The log was closed before the citations could be"
      write (6,*) "  written; they are omitted."
      return
   endif

   ! Say nothing at all rather than print an empty invitation to cite.
   anyMethodUsed = .false.
   do methodRef = 1, numMethodRefs
      if (methodWasUsed(methodRef)) then
         anyMethodUsed = .true.
      endif
   enddo
   if (.not. anyMethodUsed) return

   write (20,*)
   write (20,fmt='(a51)') sectionRule
   write (20,*) "Methods exercised by this run.  Please cite:"
   write (20,*)

   do methodRef = 1, numMethodRefs
      if (.not. methodWasUsed(methodRef)) cycle
      do line = 1, refNumLines(methodRef)
         write (20,fmt='(a)') trim(refText(line,methodRef))
      enddo
      write (20,*)
   enddo

   call flush (20)

end subroutine printMethodsBlock


end module O_MethodCitations

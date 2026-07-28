!! SPDX-License-Identifier: ECL-2.0
!! Copyright (c) 2026 Paul Rulis

!****************************************************************************
!
! Runtime output control.
!
! This module decides how much the engine is willing to say. A run states
!   what it wants through the IMAGO_VERBOSENESS environment variable, whose
!   value is a comma separated list of category names:
!
!     IMAGO_VERBOSENESS=normal          (the default when unset)
!     IMAGO_VERBOSENESS=none
!     IMAGO_VERBOSENESS=banner,timing   (order is irrelevant)
!
! Names are the public contract and bit positions are an implementation
!   detail. A number would have been easier to parse and much worse to
!   live with: this project records every invocation into a "command" file,
!   so a bit assignment written into a job script persists in the permanent
!   record, and renumbering afterward would silently change what those old
!   scripts do. Names can be reordered freely, and an unrecognized name can
!   be reported where a stray bit cannot.
!
! Callers ask "may I print this?" by passing one of the named parameters
!   below to isVerbose. They never pass a bare number. If any call site
!   hardcodes a position then renumbering stops being free, which is the
!   one property that makes this design safe to extend.
!
!****************************************************************************
module O_Verboseness

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! The number of output categories that exist. Adding a category means
   !   raising this number, declaring one parameter below, and filling in
   !   one row of initCategoryTable. Nothing else moves.
   integer, parameter :: numCategories = 1

   ! The public contract. The value of each parameter is a row of the
   !   category table, NOT a bit position. Nothing outside this module ever
   !   sees a bit position at all.
   integer, parameter :: VERB_BANNER = 1

   ! The category table. This is the single source of truth: a category's
   !   bit IS its row index minus one, so there is no second column of bit
   !   numbers that could drift out of step with the names, and the
   !   translation from name to bit exists in exactly one expression.
   character(len=20), dimension(numCategories) :: categoryName

   ! Whether each category is included in the "normal" default. This is
   !   what lets the decision "should a fresh installation print this?" be
   !   made on the same line that names the category.
   logical, dimension(numCategories) :: categoryInNormal

   ! The parsed selection, one bit per category. Private, because the whole
   !   point of the design is that the bit layout is nobody else's business.
   !   A default integer gives 31 usable positions, which is far beyond any
   !   plausible category count but is a real ceiling all the same.
   integer, private :: verbosenessMask = 0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains


! Fill in the category table. One row per category, naming it and saying
!   whether the "normal" default includes it.
subroutine initCategoryTable

   ! Make sure that no funny variables are defined.
   implicit none

   categoryName(VERB_BANNER)     = "banner"
   categoryInNormal(VERB_BANNER) = .true.

end subroutine initCategoryTable


! Read IMAGO_VERBOSENESS and set the mask from it. This must be called
!   before any isVerbose query, and in particular before the citation
!   banner, which is gated on the mask this sets. It is called from
!   parseCommandLine immediately after the log file is opened, because the
!   log unit does not exist any earlier and every complaint below is
!   written to it.
!
! Nothing here ever stops the run. A mistyped environment variable must not
!   kill a cluster job hours after it was queued. That is the opposite of
!   how elementData treats a missing IMAGO_DATA, and the difference is
!   deliberate: missing element data makes the run impossible, whereas a
!   bad verboseness request only makes the log wrong.
subroutine initVerboseness

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define the local variables used in this subroutine.
   character(len=1024) :: requestString ! The raw environment value.
   integer :: readStatus    ! 0 = read; 1 = unset; -1 = truncated.
   integer :: valueLength   ! True length, even when truncated.
   integer :: requestLength ! Length actually available to parse.
   integer :: tokenStart    ! First character of the current token.
   integer :: tokenEnd      ! Last character of the current token.
   integer :: commaOffset   ! Where the next comma sits, or zero.

   call initCategoryTable

   verbosenessMask = 0

   call get_environment_variable (NAME="IMAGO_VERBOSENESS", &
         & VALUE=requestString, LENGTH=valueLength, STATUS=readStatus)

   ! Unset is not the same as silent. A fresh installation with nothing
   !   configured must still produce a useful log, so an absent variable
   !   (status 1) and a variable holding only blanks both mean "normal".
   if ((readStatus == 1) .or. (len_trim(requestString) == 0)) then
      call applyNormalDefault
      return
   endif

   ! A status of -1 means the value was longer than the buffer and has been
   !   truncated, which would quietly drop whatever categories came last.
   !   Say so, then parse the part that did arrive.
   if (readStatus == -1) then
      write (20,*) "IMAGO_VERBOSENESS is ",valueLength," characters long"
      write (20,*) "  and was truncated to ",len(requestString),"."
      write (20,*) "  Trailing categories were ignored."
   elseif (readStatus /= 0) then
      write (20,*) "Could not read IMAGO_VERBOSENESS (status ", &
            & readStatus,")."
      write (20,*) "  Using the normal default."
      call applyNormalDefault
      return
   endif

   ! Walk the value, splitting it at commas and applying one token at a
   !   time. A zero length slice is legal in Fortran, so an empty token
   !   from ",," or a trailing comma passes through harmlessly.
   requestLength = len_trim(requestString)
   tokenStart = 1
   do while (tokenStart <= requestLength)

      commaOffset = index(requestString(tokenStart:requestLength),",")
      if (commaOffset == 0) then
         tokenEnd = requestLength
      else
         tokenEnd = tokenStart + commaOffset - 2
      endif

      call applyCategoryToken (requestString(tokenStart:tokenEnd))

      tokenStart = tokenEnd + 2
   enddo

end subroutine initVerboseness


! Turn on every category that the "normal" default includes.
subroutine applyNormalDefault

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define the local variables used in this subroutine.
   integer :: row ! Loop index over the category table.

   do row = 1, numCategories
      if (categoryInNormal(row)) then
         verbosenessMask = ibset(verbosenessMask, row - 1)
      endif
   enddo

end subroutine applyNormalDefault


! Apply one comma delimited token from IMAGO_VERBOSENESS.
!
! Tokens combine by union and nothing subtracts, so "none,banner" is
!   exactly "banner" rather than a mute switch that a later token overrides.
!   A caller wanting silence writes "none" on its own.
subroutine applyCategoryToken (rawToken)

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define the passed dummy variables.
   character(len=*), intent(in) :: rawToken

   ! Define the local variables used in this subroutine.
   character(len=64) :: candidate ! The token, trimmed and lowercased.
   integer :: row     ! Loop index over the category table.
   logical :: matched ! Whether the token named a real category.

   ! An empty token is a trailing or doubled comma. That is a typo with no
   !   ambiguity about intent, and warning about it would only train the
   !   reader to ignore the warnings that do matter.
   if (len_trim(rawToken) == 0) return

   ! Matching is case insensitive. A difference in case is not the kind of
   !   mistake the "unrecognized name" warning exists for, and reporting it
   !   as one would name a category the user can plainly see they set.
   candidate = toLowerCase(adjustl(rawToken))

   ! "normal" and "none" are set valued aliases rather than categories, and
   !   are reserved: no category may ever be given either name.
   if (trim(candidate) == "none") return
   if (trim(candidate) == "normal") then
      call applyNormalDefault
      return
   endif

   matched = .false.
   do row = 1, numCategories
      if (trim(candidate) == trim(categoryName(row))) then
         verbosenessMask = ibset(verbosenessMask, row - 1)
         matched = .true.
         exit
      endif
   enddo

   ! A name nobody recognizes must not silently do nothing, or a typo looks
   !   exactly like a feature that failed to work.
   if (.not. matched) then
      write (20,*) "Unrecognized IMAGO_VERBOSENESS category '", &
            & trim(candidate),"' ignored."
   endif

end subroutine applyCategoryToken


! Answer whether one category is active. Callers pass a named parameter,
!   such as VERB_BANNER, and never a literal number.
!
! The bounds test below catches an out of range value and nothing more. A
!   hardcoded 1 is indistinguishable from VERB_BANNER at run time, so the
!   rule that callers use names is a review discipline rather than
!   something this function can enforce. Worth saying plainly, because a
!   guard that catches part of a problem invites the belief that it catches
!   all of it.
function isVerbose (category)

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define the function result.
   logical :: isVerbose

   ! Define the passed dummy variables.
   integer, intent(in) :: category ! A category table row, e.g. VERB_BANNER

   if ((category < 1) .or. (category > numCategories)) then
      write (20,*) "isVerbose called with out of range category ", &
            & category," -- treating it as off."
      isVerbose = .false.
      return
   endif

   isVerbose = btest(verbosenessMask, category - 1)

end function isVerbose


! Return a copy of the given text with every upper case letter folded to
!   lower case. Written against the ASCII collating sequence through iachar
!   and achar rather than against the processor's native character set, so
!   the result does not depend on which one that is.
function toLowerCase (inputText)

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define the passed dummy variables.
   character(len=*), intent(in) :: inputText

   ! Define the function result.
   character(len=len(inputText)) :: toLowerCase

   ! Define the local variables used in this function.
   integer :: i        ! Loop index over the characters of the input.
   integer :: charCode ! ASCII code of the character being examined.

   do i = 1, len(inputText)
      charCode = iachar(inputText(i:i))
      if ((charCode >= iachar("A")) .and. (charCode <= iachar("Z"))) then
         toLowerCase(i:i) = achar(charCode + 32)
      else
         toLowerCase(i:i) = inputText(i:i)
      endif
   enddo

end function toLowerCase


end module O_Verboseness

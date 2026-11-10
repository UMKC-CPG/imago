!! SPDX-License-Identifier: ECL-2.0
!! Copyright (c) 2026 Paul Rulis

!****************************************************************************
!
! Stage stamps and the stage timing report.
!
! A *stage* is one named operation of a run -- the overlap integrals, the
!   electrostatic matrices, the secular equation, and so on. Each stage is
!   bracketed in the log by a start stamp and an end stamp carrying the
!   date and time (timeStampStart and timeStampEnd below).
!
! The stamps say WHEN a stage ran. The stage timing report says where its
!   time went, which is what a user needs in order to decide how many
!   ranks the next run should have (DESIGN 9.11, PSEUDOCODE 40). Inside a
!   stage three kinds of wall-clock time are told apart:
!
!!   dealt phase        Work every rank does a share of at once, such as
!!                      a loop over atom pairs. Each rank times its own
!!                      share, so the phase has one number per rank.
!!   rank-0-only phase  Work only rank 0 does while the others wait, such
!!                      as orthogonalizing and writing a matrix. One
!!                      number.
!!   stage wall time    Rank 0's clock from the start stamp to the end
!!                      stamp.
!
! A stage that wants a report records its phases between its two stamps,
!   with recordDealtPhase and recordRootPhase. The end stamp then writes
!   the report just above itself. A stage that records nothing prints its
!   stamps exactly as it always has.
!
! Every stage's calls and wall time are also added up over the whole run,
!   together with what its reports found, and printRunTiming writes those
!   totals as one table at the end of the log.
!
! What is printed is governed by two O_Verboseness categories: "timing"
!   (the short summary of each stage and the end-of-run table; on by
!   default) and "ranktiming" (one row per rank; off by default).
!
!****************************************************************************
module O_TimeStamps

   ! Import necessary modules.
   use O_Kinds
   use, intrinsic :: iso_fortran_env, only: int64

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Define local variables.
   character(len=8)  :: tempDate
   character(len=10) :: date
   character(len=10) :: tempTime
   character(len=12) :: time
   character(len=51), dimension(34) :: opLabels
   integer :: numOpCodes

   ! The most phases one stage may record between its stamps. Six dealt
   !   phases is what one row of the per-rank table can carry inside 80
   !   columns; the other two leave room for rank-0-only phases.
   integer, parameter :: maxPhases = 8

   ! The three kinds of phase. A dealt COUNT is a per-rank number that is
   !   not a time (how many potential terms each rank was dealt, say); it
   !   is summarized like a time but takes no part in the time arithmetic.
   integer, parameter :: PHASE_DEALT_TIME  = 1
   integer, parameter :: PHASE_DEALT_COUNT = 2
   integer, parameter :: PHASE_ROOT_ONLY   = 3

   ! The clocks, one entry per operation code. These accumulate over the
   !   whole run, because a stage inside the SCF cycle runs once per
   !   iteration; the totals are what an end-of-run table is built from.
   integer (kind=int64) :: clockTicksPerSecond ! system_clock's rate.
   integer (kind=int64) :: runStartTicks ! Tick when the run began.
   integer (kind=int64), dimension(34) :: opStartTicks ! Latest start.
   integer, dimension(34) :: opCalls      ! Completed start/end pairs.
   integer, dimension(34) :: opFirstOrder ! 1, 2, ... in the order the
         !   stages first started; 0 for a stage that has never run.
   integer :: numOpsStarted ! How many stages have a first order.
   real (kind=double), dimension(34) :: opWallSeconds   ! Summed wall.
   real (kind=double), dimension(34) :: opSerialSeconds ! Summed
         !   rank-0-only seconds.
   real (kind=double), dimension(34) :: opWorstBalance  ! The largest
         !   max/mean of any dealt time phase in any call; 0 if none.
   logical, dimension(34) :: opHasPhases ! True once a phase is recorded.
   logical, dimension(34) :: opOverflowReported ! True once the "too many
         !   phases" notice has been written for the stage now open.

   ! The phases recorded for stages that are open now. Only rank 0 keeps
   !   these, since only rank 0 writes the log. A slot is filled by one of
   !   the record routines and emptied by the end stamp of its stage.
   integer :: numPending
   integer, dimension(maxPhases) :: pendingOp   ! Owning operation code.
   integer, dimension(maxPhases) :: pendingKind ! One of the PHASE_ kinds.
   character(len=24), dimension(maxPhases) :: pendingLabel
   real (kind=double), dimension(maxPhases) :: pendingSeconds ! The one
         !   value of a rank-0-only phase.
   real (kind=double), allocatable, dimension(:,:) :: pendingRanks ! The
         !   per-rank values of a dealt phase; (slot, rank + 1).

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

! Name the stages and start the clocks. Called once, first thing in the
!   run and after MPI has started, on every rank.
subroutine initOperationLabels

   ! Import necessary modules.
   use O_MPI, only: mpiRank, mpiSize

   ! Make sure that no funny variables are defined.
   implicit none

   ! Set the number of operation codes.
   numOpCodes = 34

   ! Zero every accumulator and start the run clock. The tick rate is
   !   asked for once and kept, so that every later reading is turned into
   !   seconds by the same divisor.
   call system_clock (runStartTicks, clockTicksPerSecond)
   opStartTicks(:)       = 0
   opCalls(:)            = 0
   opFirstOrder(:)       = 0
   numOpsStarted         = 0
   opWallSeconds(:)      = 0.0_double
   opSerialSeconds(:)    = 0.0_double
   opWorstBalance(:)     = 0.0_double
   opHasPhases(:)        = .false.
   opOverflowReported(:) = .false.
   numPending            = 0

   ! Only rank 0 keeps recorded phases, so only rank 0 needs the store
   !   for their per-rank values. It lives for the whole run.
   if ((mpiRank == 0) .and. (.not. allocated(pendingRanks))) then
      allocate (pendingRanks(maxPhases,mpiSize))
   endif

   opLabels(1)  = '***************  Parse Input Files  ***************'
   opLabels(2)  = '***********  Get Implicit Information  ************'
   opLabels(3)  = '************  Lattice Initialization  *************'
   opLabels(4)  = '****  Init. Lattice Vector Finding Structures  ****'
   opLabels(5)  = '******  Get Exchange Correlation Mesh Counts  *****'
   opLabels(6) = '***********  Initialize SCF HDF5 File  ************'
   opLabels(7)  = '*********  Make Exchange Correlation Mesh  ********'
   opLabels(8)  = '***************  Overlap Integrals  ***************'
   opLabels(9)  = '************  Kinetic Energy Integrals  ***********'
   opLabels(10) = '**********  Nuclear Potential Integrals  **********'
   opLabels(11) = '*********  Electronic Potential Integrals  ********'
   opLabels(12) = '**************  Momentum Integrals  ***************'
   opLabels(13) = '**************  Core Charge Density  **************'
   opLabels(14) = '*************  Electrostatic Matrices  ************'
   opLabels(15) = '****************  Secular Equation  ***************'
   opLabels(16) = '************  Populate Electron Levels  ***********'
   opLabels(17) = '*************  Valence Charge Density  ************'
   opLabels(18) = '***************  Make SCF Potential  **************'
   opLabels(19) = '************  PDOS/TDOS/LI Calculation  ***********'
   opLabels(20) = '*************  Hamiltonian Integrals  *************'
   opLabels(21) = '*****  Symmetric Band Structure Path KPoints  *****'
   opLabels(22) = '********  Bond Order and Effective Charge  ********'
   opLabels(23) = '*********  Optical Transition Calculation  ********'
   opLabels(24) = '*************  Parse The Command Line  ************'
   opLabels(25) = '*************  Evaluate Data on Mesh  *************'
   opLabels(26) = '************  Three Center Bond Order  ************'
   opLabels(27) = '***********  Initialize Field HDF5 File  **********'
   opLabels(28) = '**********  Compute Bispectrum Component  *********'
   opLabels(29) = '************  Dipole Moment Integrals  ************'
   opLabels(30) = '************  Mass Velocity Integrals  ************'
   opLabels(31) = '***********  Initialize PSCF HDF5 File  ***********'
   opLabels(32) = '**************  KOverlap Integrals  ***************'
   opLabels(33) = '**** Modern Theory of Polarization Calculation ****'
   opLabels(34) = '*****************  Dipole Moment  *****************'

end subroutine initOperationLabels

subroutine timeStampStart (opCode)

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare passed dummy variables.
   integer :: opCode

   ! Make sure we are given a valid opCode.
   if ((opCode > numOpCodes) .or. (opCode < 1)) then
      write (20,*) "opCode range is (1,",numOpCodes,")"
      write (20,*) "Given opCode is:  ",opCode
      stop
   endif

   ! Start this stage's clock, and note its place in the order in which
   !   the stages of this run first began (PSEUDOCODE 40.5).
   call system_clock (opStartTicks(opCode))
   opOverflowReported(opCode) = .false.
   if (opFirstOrder(opCode) == 0) then
      numOpsStarted = numOpsStarted + 1
      opFirstOrder(opCode) = numOpsStarted
   endif

   ! Print the leading banner followed by the appropriate label.
   write (20,*)
   write (20,fmt='(a51)') '***************************************************'
   write (20,fmt='(a51)') opLabels(opCode)

   ! Log the date and time
   call date_and_time(tempDate,tempTime)
   date = tempDate(1:4)//'/'//tempDate(5:6)//'/'//tempDate(7:)
   time = tempTime(1:2)//':'//tempTime(3:4)//':'//tempTime(5:)
   write (20,100) date,time
   write (20,*)
   call flush (20)

   100 format ('Date is: ',a10,' Time is: ',a12) ! Date and time output format

end subroutine timeStampStart

subroutine timeStampEnd (opCode)

   ! Import necessary modules.
   use O_Kinds
   use O_MPI, only: mpiRank

   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare passed dummy variables.
   integer :: opCode

   ! Define local variables.
   integer (kind=int64) :: ticksNow ! The clock as this stage ends.
   real (kind=double) :: wallSeconds ! This call's start-to-end time.
   integer :: slot      ! Loop index over the recorded phases.
   integer :: keptSlots ! How many phases remain after this stage's go.

   ! Make sure we are given a valid opCode.
   if ((opCode > numOpCodes) .or. (opCode < 1)) then
      write (20,*) "opCode range is between 0 and ",numOpCodes
      write (20,*) "Given opCode is:  ",opCode
      stop
   endif

   ! Stop this stage's clock and add the call to the stage's run totals
   !   (PSEUDOCODE 40.5). A stage ended without ever being started has no
   !   start tick to measure from and is given no time.
   wallSeconds = 0.0_double
   if (opStartTicks(opCode) /= 0) then
      call system_clock (ticksNow)
      wallSeconds = real (ticksNow - opStartTicks(opCode), double) &
            & / real (clockTicksPerSecond, double)
   endif
   opCalls(opCode)       = opCalls(opCode) + 1
   opWallSeconds(opCode) = opWallSeconds(opCode) + wallSeconds

   ! If this stage recorded any phases, write its timing report now, so
   !   that it sits directly above the end stamp it explains, and then
   !   release the slots. Phases belonging to some other stage that is
   !   still open are kept, moved down to close the gaps.
   if ((mpiRank == 0) .and. (any(pendingOp(1:numPending) == opCode))) then
      call writeStageReport (opCode, wallSeconds)
      keptSlots = 0
      do slot = 1, numPending
         if (pendingOp(slot) == opCode) cycle
         keptSlots = keptSlots + 1
         pendingOp(keptSlots)      = pendingOp(slot)
         pendingKind(keptSlots)    = pendingKind(slot)
         pendingLabel(keptSlots)   = pendingLabel(slot)
         pendingSeconds(keptSlots) = pendingSeconds(slot)
         pendingRanks(keptSlots,:) = pendingRanks(slot,:)
      enddo
      numPending = keptSlots
   endif

   ! Log the date and time.
   write (20,*)
   call date_and_time(tempDate,tempTime)
   date = tempDate(1:4)//'/'//tempDate(5:6)//'/'//tempDate(7:)
   time = tempTime(1:2)//':'//tempTime(3:4)//':'//tempTime(5:)
   write (20,100) date,time

   ! Print the appropriate label followed by the trailing banner.
   write (20,fmt='(a51)') opLabels(opCode)
   write (20,fmt='(a51)') '***************************************************'
   write (20,*)
   call flush (20)

   100 format ('Date is: ',a10,' Time is: ',a12) ! Date and time output format

end subroutine timeStampEnd


! Claim the next free slot for a phase of the stage opCode, and return its
!   index, or zero when the stage has already recorded maxPhases phases.
!   Running out of slots must never stop a run -- this is only a report --
!   so the phase is dropped, and the log says so once for the stage.
function claimPhaseSlot (opCode, label)

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define the function result.
   integer :: claimPhaseSlot

   ! Declare passed dummy variables.
   integer, intent(in) :: opCode
   character(len=*), intent(in) :: label

   claimPhaseSlot = 0
   if (numPending == maxPhases) then
      if (.not. opOverflowReported(opCode)) then
         write (20,fmt='(a,i0,a)') ' Stage timing: more than ', &
               & maxPhases, ' phases recorded; "'//trim(label)// &
               & '" and any later ones are dropped.'
         opOverflowReported(opCode) = .true.
      endif
      return
   endif

   numPending = numPending + 1
   pendingOp(numPending)      = opCode
   pendingLabel(numPending)   = label
   pendingSeconds(numPending) = 0.0_double
   pendingRanks(numPending,:) = 0.0_double
   claimPhaseSlot = numPending

end function claimPhaseSlot


! Record one dealt phase of the open stage opCode: work that every rank
!   did a share of. Each rank passes the time its own share took (or, with
!   isCount true, a count of the things it was dealt).
!
! THIS CALL IS COLLECTIVE. Every rank that took part in the phase must
!   make it, and where a stage records several phases every rank must
!   record them in the same order -- the values travel to rank 0 through
!   an MPI gather, which only completes when all ranks have joined it. So
!   it belongs at a point every rank reaches, before any rank leaves the
!   stage.
!
! The label names the phase in the summary and heads its column in the
!   per-rank table. Keep it to 11 characters so that it fits the column.
subroutine recordDealtPhase (opCode, label, myValue, isCount)

   ! Import necessary modules.
   use O_MPI, only: mpiRank, mpiSize, gatherTimesMPI

   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare passed dummy variables.
   integer, intent(in) :: opCode ! The stage this phase belongs to.
   character(len=*), intent(in) :: label
   real (kind=double), intent(in) :: myValue ! This rank's own value.
   logical, intent(in), optional :: isCount  ! True for a count; a time
         !   in seconds when false or absent.

   ! Define local variables.
   real (kind=double), dimension(1) :: myValues ! myValue, as the vector
         !   that the gather routine expects.
   real (kind=double), dimension(1,mpiSize) :: allValues ! Every rank's
         !   value, side by side; filled on rank 0 only.
   integer :: slot ! The slot this phase is stored in.

   ! Bring every rank's value to rank 0. In the serial build, and with
   !   one rank, this copies rank 0's own value across.
   myValues(1) = myValue
   call gatherTimesMPI (myValues, allValues)

   ! Only rank 0 keeps what was gathered.
   if (mpiRank /= 0) return

   slot = claimPhaseSlot (opCode, label)
   if (slot == 0) return

   pendingKind(slot) = PHASE_DEALT_TIME
   if (present(isCount)) then
      if (isCount) pendingKind(slot) = PHASE_DEALT_COUNT
   endif
   pendingRanks(slot,:) = allValues(1,:)

end subroutine recordDealtPhase


! Record one rank-0-only phase of the open stage opCode: work that rank 0
!   alone did while the other ranks waited, and which adding ranks
!   therefore cannot shorten.
!
! This call is NOT collective, and it may be made without a rank test: on
!   any rank other than 0 it does nothing.
!
! Recording the same label twice within one open stage ADDS the times
!   into a single line. A stage made of several parts that each end with,
!   say, a write to the file can so report one "write" figure for the
!   whole stage.
subroutine recordRootPhase (opCode, label, seconds)

   ! Import necessary modules.
   use O_MPI, only: mpiRank

   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare passed dummy variables.
   integer, intent(in) :: opCode ! The stage this phase belongs to.
   character(len=*), intent(in) :: label
   real (kind=double), intent(in) :: seconds ! Rank 0's time in it.

   ! Define local variables.
   integer :: slot ! Loop index, then the slot the phase is stored in.

   if (mpiRank /= 0) return

   ! Add to a phase of this stage that already carries this label.
   do slot = 1, numPending
      if ((pendingOp(slot) == opCode) .and. &
            & (pendingKind(slot) == PHASE_ROOT_ONLY) .and. &
            & (trim(pendingLabel(slot)) == trim(label))) then
         pendingSeconds(slot) = pendingSeconds(slot) + seconds
         return
      endif
   enddo

   ! Otherwise this is a new phase.
   slot = claimPhaseSlot (opCode, label)
   if (slot == 0) return
   pendingKind(slot)    = PHASE_ROOT_ONLY
   pendingSeconds(slot) = seconds

end subroutine recordRootPhase


! Write the timing report of the stage opCode from the phases it has
!   recorded, and fold them into the stage's run totals (PSEUDOCODE 40.6).
!   Called by timeStampEnd, on rank 0, with the wall time of the call
!   that is ending.
!
! The numbers, and how to read them:
!
!!   min, mean, max   Of a dealt phase, over the ranks.
!!   max/mean         The balance of the deal. 1.0 is a perfectly even
!!                    deal; 2.0 means the slowest rank took twice the
!!                    average, so the phase ran at half the speed an
!!                    even deal would have given.
!!   rank 0 only      Seconds of work that only rank 0 did.
!!   not timed        The stage wall time, less rank 0's own share of
!!                    every dealt phase, less every rank-0-only phase.
!!                    Small means the phases account for the stage;
!!                    large means the stage holds work that nobody put
!!                    a clock on, or rank 0 waiting on slower ranks.
!!   serial share     The rank-0-only seconds as a percentage of the
!!                    stage wall time: the part that more ranks cannot
!!                    shorten. "Not timed" is NOT counted in it, since
!!                    what that time holds is unknown.
subroutine writeStageReport (opCode, wallSeconds)

   ! Import necessary modules.
   use O_MPI, only: mpiSize
   use O_Verboseness, only: isVerbose, VERB_TIMING, VERB_RANKTIMING

   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare passed dummy variables.
   integer, intent(in) :: opCode
   real (kind=double), intent(in) :: wallSeconds

   ! Define local variables.
   integer :: slot ! Loop index over the recorded phases.
   integer :: rank ! Loop index over the ranks, counted from one.
   integer :: numDealt ! How many dealt phases this stage recorded.
   real (kind=double), dimension(maxPhases) :: lowest, mean, highest
   real (kind=double), dimension(maxPhases) :: balance ! highest / mean.
   real (kind=double) :: rankZeroDealtSeconds ! Rank 0's own time, summed
         !   over the dealt time phases.
   real (kind=double) :: serialSeconds ! Summed rank-0-only phases.
   real (kind=double) :: worstBalance  ! Largest balance of a time phase.
   real (kind=double) :: notTimedSeconds
   integer :: serialSharePercent
   character(len=16) :: wallText ! The wall time, for the heading.
   character(len=5)  :: rankWord ! "rank" or "ranks".
   character(len=24) :: rowLabel ! A phase label, as the summary shows it.
   character(len=37) :: wideLabel ! The label of a line with one value.
         !   (A character variable is used for these labels, not a
         !   literal, because an "a" edit descriptor wider than its text
         !   pads on the LEFT; a variable is already padded on the right,
         !   which keeps the labels flush with the left margin.)
   character(len=12) :: columnHead ! A phase label, as a column heading.

   ! Reduce each phase to the numbers the report needs.
   rankZeroDealtSeconds = 0.0_double
   serialSeconds        = 0.0_double
   worstBalance         = 0.0_double
   numDealt             = 0
   do slot = 1, numPending
      if (pendingOp(slot) /= opCode) cycle
      if (pendingKind(slot) == PHASE_ROOT_ONLY) then
         serialSeconds = serialSeconds + pendingSeconds(slot)
         cycle
      endif

      ! A dealt phase, whether a time or a count.
      numDealt      = numDealt + 1
      lowest(slot)  = minval(pendingRanks(slot,:))
      highest(slot) = maxval(pendingRanks(slot,:))
      mean(slot)    = sum(pendingRanks(slot,:)) / real(mpiSize,double)
      if (mean(slot) > 0.0_double) then
         balance(slot) = highest(slot) / mean(slot)
      else
         balance(slot) = 1.0_double
      endif

      ! Only times enter the time arithmetic and the balance record.
      if (pendingKind(slot) == PHASE_DEALT_TIME) then
         rankZeroDealtSeconds = rankZeroDealtSeconds + pendingRanks(slot,1)
         worstBalance = max(worstBalance, balance(slot))
      endif
   enddo

   ! "Not timed" is a difference of clock readings, so when the phases
   !   account for the whole stage it lands a rounding error either side
   !   of zero. Anything smaller than half the last printed digit is set
   !   to exactly zero so that the line never reads "-0.000". A larger
   !   negative value is left alone: it says the phases claim more time
   !   than the stage took, and that fault in the clocks must be seen.
   notTimedSeconds = wallSeconds - rankZeroDealtSeconds - serialSeconds
   if (abs(notTimedSeconds) < 0.0005_double) notTimedSeconds = 0.0_double
   if (wallSeconds > 0.0_double) then
      serialSharePercent = nint(100.0_double * serialSeconds / wallSeconds)
   else
      serialSharePercent = 0
   endif

   ! Fold this call into the stage's run totals. This is done whatever is
   !   printed below: what a run chooses to print must not change what
   !   its totals say.
   opSerialSeconds(opCode) = opSerialSeconds(opCode) + serialSeconds
   opWorstBalance(opCode)  = max(opWorstBalance(opCode), worstBalance)
   opHasPhases(opCode)     = .true.

   ! The summary: one line per phase.
   if (isVerbose(VERB_TIMING)) then
      write (wallText,fmt='(f16.1)') wallSeconds
      rankWord = 'ranks'
      if (mpiSize == 1) rankWord = 'rank'
      write (20,*)
      write (20,fmt='(a,i0,a)') ' Stage timing: ', mpiSize, ' '// &
            & trim(rankWord)//', '//trim(adjustl(wallText))//' s wall'
      if (numDealt > 0) then
         rowLabel = 'phase'
         write (20,fmt='(3x,a24,a9,2a11,a10)') rowLabel, 'min', 'mean', &
               & 'max', 'max/mean'
      endif
      do slot = 1, numPending
         if (pendingOp(slot) /= opCode) cycle
         if (pendingKind(slot) == PHASE_DEALT_TIME) then
            write (20,fmt='(3x,a24,f9.3,2f11.3,f10.1)') &
                  & pendingLabel(slot), lowest(slot), mean(slot), &
                  & highest(slot), balance(slot)
         elseif (pendingKind(slot) == PHASE_DEALT_COUNT) then
            rowLabel = trim(pendingLabel(slot))//' (count)'
            write (20,fmt='(3x,a24,i9,f11.1,i11,f10.1)') rowLabel, &
                  & nint(lowest(slot)), mean(slot), nint(highest(slot)), &
                  & balance(slot)
         endif
      enddo
      do slot = 1, numPending
         if (pendingOp(slot) /= opCode) cycle
         if (pendingKind(slot) /= PHASE_ROOT_ONLY) cycle
         write (20,fmt='(3x,a13,a24,f18.3)') 'rank 0 only: ', &
               & pendingLabel(slot), pendingSeconds(slot)
      enddo
      wideLabel = 'not timed'
      write (20,fmt='(3x,a37,f18.3)') wideLabel, notTimedSeconds
      wideLabel = 'serial share of this stage'
      write (20,fmt='(3x,a37,i16,a2)') wideLabel, serialSharePercent, ' %'
   endif

   ! The per-rank table: one row per rank, one 12-character column per
   !   dealt phase, in the order the phases were recorded.
   if (isVerbose(VERB_RANKTIMING) .and. (numDealt > 0)) then
      write (20,*)
      write (20,fmt='(a)') ' Stage timing by rank:'
      write (20,fmt='(3x,a4)',advance='no') 'rank'
      do slot = 1, numPending
         if (pendingOp(slot) /= opCode) cycle
         if (pendingKind(slot) == PHASE_ROOT_ONLY) cycle
         columnHead = pendingLabel(slot)(1:12)
         write (20,fmt='(a12)',advance='no') adjustr(columnHead)
      enddo
      write (20,*)
      do rank = 1, mpiSize
         write (20,fmt='(i7)',advance='no') rank - 1
         do slot = 1, numPending
            if (pendingOp(slot) /= opCode) cycle
            if (pendingKind(slot) == PHASE_DEALT_TIME) then
               write (20,fmt='(f12.3)',advance='no') &
                     & pendingRanks(slot,rank)
            elseif (pendingKind(slot) == PHASE_DEALT_COUNT) then
               write (20,fmt='(i12)',advance='no') &
                     & nint(pendingRanks(slot,rank))
            endif
         enddo
         write (20,*)
      enddo
   endif

end subroutine writeStageReport


! Return the plain name of the stage opCode: its banner label with the
!   rows of asterisks and the blanks around the words taken off, so that
!   '***************  Overlap Integrals  ***************' comes back as
!   'Overlap Integrals'. The result is 41 characters long because that is
!   the length of the longest name among the labels.
function stageName (opCode)

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define the function result.
   character(len=41) :: stageName

   ! Declare passed dummy variables.
   integer, intent(in) :: opCode

   ! Define local variables.
   integer :: firstChar ! Where the first word of the name begins.
   integer :: lastChar  ! Where the last word of the name ends.

   ! Step inward from each end of the label, past asterisks and blanks.
   !   verify returns the position of the first character (or, with
   !   back=.true., the last) that is NOT one of those two.
   firstChar = verify(opLabels(opCode), '* ')
   lastChar  = verify(opLabels(opCode), '* ', back=.true.)

   ! A label holding nothing but asterisks and blanks has no name.
   stageName = ' '
   if (firstChar > 0) stageName = opLabels(opCode)(firstChar:lastChar)

end function stageName


! Write the end-of-run timing table: one row for every stage that ran,
!   gathering in one place what the stage summaries said along the way
!   (PSEUDOCODE 40.7). This is the table to read when choosing how many
!   ranks the next run of a similar system should have. Called once, by
!   Imago, after all the work is done and before the log is closed.
!
! The columns, and how to read them:
!
!!   calls      How many times the stage ran. A stage inside the SCF cycle
!!              runs once per iteration.
!!   wall s     The stage's wall-clock seconds, summed over its calls.
!!   % run      That time as a percentage of the whole run. The stages
!!              with large numbers here are the ones worth speeding up.
!!   serial %   The seconds of work that only rank 0 did, as a percentage
!!              of the stage's wall seconds: the part of the stage that
!!              more ranks cannot shorten.
!!   balance    The largest max/mean that any dealt phase of the stage
!!              showed in any call. 1.0 is a perfectly even deal.
!
! A dash in the last two columns means the stage recorded no phases, so
!   nothing is known about how its time divides. That is NOT the same as
!   a serial share of zero: most stages are simply not instrumented yet.
!
! Rows come in the order the stages first started, which is the order a
!   reader meets them in the log. The last two rows are the sum of the
!   wall column ("all stages") and the run's own clock ("run"). They are
!   printed side by side so that time spent outside every stage shows as
!   a difference between them instead of hiding.
subroutine printRunTiming

   ! Import necessary modules.
   use O_MPI, only: mpiRank, mpiSize
   use O_Verboseness, only: isVerbose, VERB_TIMING

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define local variables.
   integer (kind=int64) :: ticksNow ! The clock as the run ends.
   real (kind=double) :: runSeconds ! The whole run, start to now.
   real (kind=double) :: allStagesSeconds ! The sum of the wall column.
   real (kind=double) :: runSharePercent  ! A wall time as a % of the run.
   integer :: order  ! Loop index over the order of first starts.
   integer :: opCode ! Loop index over the operation codes.
   logical :: logIsOpen ! Whether unit 20 is still connected.
   character(len=16) :: wallText ! The run's wall time, for the heading.
   character(len=5)  :: rankWord ! "rank" or "ranks".
   character(len=41) :: rowLabel ! The first column of a row. (A variable
         !   and not a literal, so that it is padded on the right; see
         !   the same note in writeStageReport.)
   character(len=9)  :: serialText  ! The serial share, or a dash.
   character(len=8)  :: balanceText ! The worst balance, or a dash.

   ! Only rank 0 writes the log, and the table belongs to the "timing"
   !   category along with the stage summaries it gathers.
   if (mpiRank /= 0) return
   if (.not. isVerbose(VERB_TIMING)) return

   ! The log must still be open. Writing to a closed unit does not fail:
   !   Fortran reconnects it to fort.20 and truncates, which would replace
   !   the run's entire log with this table. Complain where it can still
   !   be seen, and write nothing.
   inquire (unit=20, opened=logIsOpen)
   if (.not. logIsOpen) then
      write (6,*) "The log was closed before the run timing table could"
      write (6,*) "  be written; it is omitted."
      return
   endif

   ! Read the run clock that initOperationLabels started.
   call system_clock (ticksNow)
   runSeconds = real (ticksNow - runStartTicks, double) &
         & / real (clockTicksPerSecond, double)

   ! The heading and the column titles.
   write (wallText,fmt='(f16.1)') runSeconds
   rankWord = 'ranks'
   if (mpiSize == 1) rankWord = 'rank'
   write (20,*)
   write (20,fmt='(a,i0,a)') ' Run timing: ', mpiSize, ' '// &
         & trim(rankWord)//', '//trim(adjustl(wallText))//' s wall'
   rowLabel = 'stage'
   write (20,fmt='(1x,a41,a6,a9,a6,a9,a8)') rowLabel, 'calls', 'wall s', &
         & '% run', 'serial %', 'balance'

   ! One row per stage, in the order the stages first started. Each pass
   !   looks up the one operation code holding that place in the order.
   allStagesSeconds = 0.0_double
   do order = 1, numOpsStarted
      do opCode = 1, numOpCodes
         if (opFirstOrder(opCode) /= order) cycle

         ! A stage that was started but never ended has no completed
         !   call to report.
         if (opCalls(opCode) == 0) cycle

         allStagesSeconds = allStagesSeconds + opWallSeconds(opCode)
         runSharePercent  = percentOf (opWallSeconds(opCode), runSeconds)

         ! The last two columns are known only for a stage that recorded
         !   phases. The balance needs a dealt time phase besides: a stage
         !   with rank-0-only phases alone has a serial share and no deal
         !   to be balanced, which its zero worst balance marks.
         serialText  = adjustr('-        ')
         balanceText = adjustr('-       ')
         if (opHasPhases(opCode)) then
            write (serialText,fmt='(i9)') nint (percentOf ( &
                  & opSerialSeconds(opCode), opWallSeconds(opCode)))
            if (opWorstBalance(opCode) > 0.0_double) then
               write (balanceText,fmt='(f8.1)') opWorstBalance(opCode)
            endif
         endif

         write (20,fmt='(1x,a41,i6,f9.1,f6.1,a9,a8)') stageName(opCode), &
               & opCalls(opCode), opWallSeconds(opCode), runSharePercent, &
               & serialText, balanceText
      enddo
   enddo

   ! The two closing rows: what the stages add up to, and the run itself.
   rowLabel = 'all stages'
   write (20,fmt='(1x,a41,6x,f9.1,f6.1)') rowLabel, allStagesSeconds, &
         & percentOf (allStagesSeconds, runSeconds)
   rowLabel = 'run'
   write (20,fmt='(1x,a41,6x,f9.1,f6.1)') rowLabel, runSeconds, &
         & 100.0_double
   call flush (20)

end subroutine printRunTiming


! Return part as a percentage of whole, or zero when whole is zero (a
!   stage, or a whole run, too short for the clock to have ticked).
function percentOf (part, whole)

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define the function result.
   real (kind=double) :: percentOf

   ! Declare passed dummy variables.
   real (kind=double), intent(in) :: part
   real (kind=double), intent(in) :: whole

   percentOf = 0.0_double
   if (whole > 0.0_double) percentOf = 100.0_double * part / whole

end function percentOf

end module O_TimeStamps

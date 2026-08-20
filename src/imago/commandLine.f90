!! SPDX-License-Identifier: ECL-2.0
!! Copyright (c) 2026 Paul Rulis

module O_CommandLine

   ! Import necessary modules.
   use O_Kinds

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Integer to track which number command line argument is to be read in next.
   integer :: nextArg

   ! Flags that will be set to indicate that SCF and/or PSCF calculations
   !   should be done.
   integer :: doSCF
   integer :: doPSCF

   ! Choice of basis for SCF and PSCF stages.
   integer :: basisCode_SCF ! 0=NO; 1=MB; 2=FB; 3=EB
   integer :: basisCode_PSCF ! 0=NO; 1=MB; 2=FB; 3=EB

   ! Core electron excitation control paramters.
   integer :: excitedQN_n ! 1, 2, 3, 4, ...
   integer :: excitedQN_l ! 0=s, 1=p, 2=d, 3=f, ...

   ! Overall job ID.
   integer :: jobID

   ! Properties to compute in the SCF stage.
   integer :: doDOS_SCF  ! Include DOS/PDOS calculation.
   integer :: doBond_SCF ! Include bond order calculation.
   integer :: doDIMO_SCF ! Include dipole moment matrix elements.
   integer :: doOPTC_SCF ! Include optical properties. -1=none; 1=valence band
         ! optical properties; 2=core level XANES/ELNES; 3=sigma(E) electronic
         ! contribution to thermal conductivity; 4=non-linear valence band
         ! optical properties.
   integer :: doSYBD_SCF ! Shift to using the defined path of high-symmetry
         ! k-points (1) AND produce a band structure diagram (2), AND include
         ! partial band structure data (3).
   integer :: doForce_SCF ! Include computation of the force between atoms.
   integer :: doField_SCF ! Compute a charge, potential, or wave fn field.
   integer :: doMTOP_SCF ! Compute the polarization using the modern theory.

   ! Properties to compute in the PSCF stage.
   integer :: doDOS_PSCF  ! Include DOS/PDOS calculation.
   integer :: doBond_PSCF ! Include bond order calculation.
   integer :: doDIMO_PSCF ! Include dipole moment matrix elements.
   integer :: doOPTC_PSCF ! Include optical properties. -1=none; 1=valence band
         ! optical properties; 2=core level XANES/ELNES; 3=sigma(E) electronic
         ! contribution to thermal conductivity; 4=non-linear valence band
         ! optical properties.
   integer :: doSYBD_PSCF ! Shift to using the defined path of high-symmetry
         ! k-points (1) AND produce a band structure diagram (2), AND include
         ! partial band structure data (3).
   integer :: doForce_PSCF ! Include computation of the force between atoms.
   integer :: doField_PSCF ! Compute a charge, potential, or wave fn field.
   integer :: doMTOP_PSCF ! Compute the polarization using the modern theory.

   ! Properties to compute based only on geometry and not the wave function.
   integer :: doLoEn ! Compute a metric that quantifies the local environment.

   ! Compute any XYZ components in serial.
   integer :: serialXYZ

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains


subroutine parseCommandLine

   ! Use necessary modules.
   use O_TimeStamps
   use O_Verboseness, only: initVerboseness
   use O_Banner, only: printIdentityBlock
   use O_MPI, only: mpiRank

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define the local variables that will be used to parse the command line.
   character*25 :: commandBuffer
   character*16 :: rankLogName
   integer :: rankLogEnvLen  ! Length of $IMAGO_RANK_LOGS; 0 if unset.

   ! Open the file that will be written to as output for this program.
   !   Under MPI every rank executes this subroutine, and two hazards
   !   shape the branch below (PSEUDOCODE 24.4). First, N ranks
   !   appending to one file would interleave and corrupt the log.
   !   Second, a rank may NOT simply skip the open: a Fortran write to
   !   a unit that was never opened does not fail -- it silently
   !   auto-connects unit 20 to a file named fort.20 -- so a worker
   !   with no connection of its own would tear root's log with its
   !   first stray write. Every rank therefore connects unit 20 here.
   !
   !   Root keeps fort.20: byte-for-byte the serial log, and the file
   !   imago.py reads. A worker's routine output is a replica of
   !   root's, and a run on thousands of ranks must not shed
   !   thousands of log files, so by default the workers connect to
   !   /dev/null and their chatter is discarded. Setting
   !   IMAGO_RANK_LOGS (to any non-empty value) switches the workers
   !   to per-rank fort.20.rNNNN files for debugging sessions where
   !   one rank's view of events is the question. Errors never depend
   !   on the worker log: stopMPI writes rank-stamped messages to
   !   standard error, which the MPI launcher aggregates into the one
   !   job error file. In the serial build mpiRank is 0 and this
   !   collapses to the original open.
   if (mpiRank == 0) then
      open(20,file='fort.20',status='unknown',form='formatted')
   else
      call get_environment_variable('IMAGO_RANK_LOGS',&
            & length=rankLogEnvLen)
      if (rankLogEnvLen > 0) then
         write(rankLogName,fmt='(a9,i4.4)') 'fort.20.r', mpiRank
         open(20,file=trim(rankLogName),status='unknown',&
               & form='formatted')
      else
         open(20,file='/dev/null',status='old',form='formatted')
      endif
   endif

   ! Establish what this run is willing to print, and then print the
   !   identity block. This slot is forced from both sides: nothing may
   !   print earlier because unit 20 does not exist earlier, and nothing
   !   may print later because the timestamp below writes the first rule of
   !   the log, which the banner must sit above.
   !
   ! The order of these two calls is forced as well, and a mistake here is
   !   invisible. Reversed, printIdentityBlock tests a mask that is still
   !   zero and returns without printing or complaining, so the run
   !   succeeds and the banner is simply missing.
   call initVerboseness
   call printIdentityBlock

   ! Record the date and time that we start.
   call timeStampStart (24)

   ! Initialize all the command line parameters.
   call initCLP

   ! Begin Parsing the command line.

   call readBasisCodes

   call readExcitedQN

   call readJobID

   ! Read a flag to request that any XYZ based calculation be done in serial.
   call getarg(nextArg,commandBuffer)
   nextArg = nextArg + 1
   read (commandBuffer,*) serialXYZ
   write (20,*) "serialXYZ = ",serialXYZ

   ! Record the date and time that we end.
   call timeStampEnd (24)

end subroutine parseCommandLine


subroutine readBasisCodes

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define the local variables that will be used to parse the command line.
   character*25 :: commandBuffer

   ! Get the command line argument that defines the basis set to use for the
   !   SCF portion of the calculation.  0=NO; 1=MB; 2=FB; 3=EB
   call getarg(nextArg,commandBuffer)
   nextArg = nextArg + 1
   read (commandBuffer,*) basisCode_SCF

   ! If the SCF basis code is a 0, that indicates that we will not perform an
   !   SCF calculation. If it is non-zero, then we do an SCF calculation.
   if (basisCode_SCF /= 0) then
      doSCF = 1
   endif

   ! Record the basis code in the output.
   write (20,*) "basisCode_SCF = ",basisCode_SCF
   write (20,*) "0=NO; 1=MB; 2=FB; 3=EB"
   write (20,*)

   ! Get the command line argument that defines the basis set to use for the
   !   PSCF portion of the calculation.  1=MB; 2=FB; 3=EB
   call getarg(nextArg,commandBuffer)
   nextArg = nextArg + 1
   read (commandBuffer,*) basisCode_PSCF

   ! If the PSCF basis code is a 0, that indicates that we will not perform a
   !   PSCF calculation. If it is non-zero, then we do some PSCF calculation.
   if (basisCode_PSCF /= 0) then
      doPSCF = 1
   endif

   ! Record the basis code in the output.
   write (20,*) "basisCode_PSCF = ",basisCode_PSCF
   write (20,*) "0=NO; 1=MB; 2=FB; 3=EB"
   write (20,*)

   ! If both basis codes are 0, then set the PSCF code to 1 so that input
   !   files can be easily read for non-SCF and non-PSCF calculations.
   if ((basisCode_SCF == 0) .and. (basisCode_PSCF == 0)) then
      basisCode_PSCF = 1
      write (20,*) "Both basis code requested to be zero."
      write (20,*) "Setting basisCode_PSCF to 1."
      write (20,*)
   endif


end subroutine readBasisCodes



subroutine readExcitedQN

   ! Make sure that there are no accidental variable declarations.
   implicit none

   ! Define the local variables that will be used to parse the command line.
   character*25 :: commandBuffer

   ! Store the command line argument that defines the QN_n of which electron
   !   will be excited.  (0=ground state; 1=K; 2=L; 3=M; 4=N; ...)
   call getarg(nextArg,commandBuffer)
   nextArg = nextArg + 1
   read (commandBuffer,*) excitedQN_n

   ! Store the command line argument that defines the QN_l of which electron
   !   will be excited.  (0=s; 1=p; 2=d; 3=f; ...)
   call getarg(nextArg,commandBuffer)
   nextArg = nextArg + 1
   read (commandBuffer,*) excitedQN_l

   ! Check to make sure that QN_l < QN_n (or that QN_l == QN_n == 0 (gs)).
   if ((excitedQN_l < excitedQN_n) .or. &
         & ((excitedQN_n == 0) .and. (excitedQN_l==0))) then
      write (20,*) "excitedQN_n = ",excitedQN_n
      write (20,*) "excitedQN_l = ",excitedQN_l
      write (20,*)
   else
      write (20,*) "CLP 'excitedQN_n' = ",excitedQN_n
      write (20,*) "CLP 'excitedQN_l' = ",excitedQN_l
      write (20,*) "QN_l should be < QN_n"
      stop
   endif

end subroutine readExcitedQN


subroutine readJobID

   ! Make sure that nothing funny is declared.
   implicit none

   ! Define the local variables that will be used to parse the command line.
   character*25 :: commandBuffer

   ! Read a flag indicating that a dipole moment calculation should be tacked
   !   on to the end of the SCF iterations.
   call getarg(nextArg,commandBuffer)
   nextArg = nextArg + 1
   read (commandBuffer,*) jobID
   write (20,*) "jobID = ",jobID

   if (jobID == 0) then
      write (20,*) "Doing SCF Total Energy Only"
   elseif (jobID == 101) then
      doDOS_SCF = 1
      write (20,*) "Doing SCF Density of States"
   elseif (jobID == 102) then
      doBond_SCF = 1
      write (20,*) "Doing SCF Bond Order and Q*"
   elseif (jobID == 103) then
      doDIMO_SCF = 1
      write (20,*) "Doing SCF Dipole Moment"
   elseif (jobID == 104) then
      doOPTC_SCF = 1
      write (20,*) "Doing SCF Valence Band Optical Properties"
   elseif (jobID == 105) then
      doOPTC_SCF = 2
      write (20,*) "Doing SCF Photo-Absorption Cross Section"  ! XANES/ELNES
   elseif (jobID == 106) then
      ! The doOPTC codes are the ones declared above and acted on inside
      !   the optc module: 3 selects sigma(E) and 4 selects the non-linear
      !   properties. They are not the same numbering as the job ID, where
      !   106 is non-linear and 107 is sigma(E), so the two cross over here.
      doOPTC_SCF = 4
      write (20,*) "Doing SCF Non-Linear Optical Properties"
   elseif (jobID == 107) then
      doOPTC_SCF = 3
      write (20,*) "Doing SCF Sigma(E)"
   elseif (jobID == 108) then
      doSYBD_SCF = 1
      write (20,*) "Doing SCF Symmetric Band Structure"
   elseif (jobID == 109) then
      doForce_SCF = 1
      write (20,*) "Doing SCF Force"
   elseif (jobID == 110) then
      doField_SCF = 1
      write (20,*) "Doing SCF Field"
   elseif (jobID == 111) then
      doMTOP_SCF = 1
      write (20,*) "Doing SCF Modern Polarization"
   elseif (jobID == 201) then
      doDOS_PSCF = 1
      write (20,*) "Doing PSCF Density of States"
   elseif (jobID == 202) then
      doBond_PSCF = 1
      write (20,*) "Doing PSCF Bond Order and Q*"
   elseif (jobID == 203) then
      doDIMO_PSCF = 1
      write (20,*) "Doing PSCF Dipole Moment"
   elseif (jobID == 204) then
      doOPTC_PSCF = 1
      write (20,*) "Doing PSCF Valence Band Optical Properties"
   elseif (jobID == 205) then
      doOPTC_PSCF = 2
      write (20,*) "Doing PSCF Photo-Absorption Cross Section"  ! XANES/ELNES
   elseif (jobID == 206) then
      ! See the note on the matching SCF case above: the doOPTC code and
      !   the job ID number the last two optical properties in opposite
      !   order, so 206 (non-linear) maps to code 4 and 207 to code 3.
      doOPTC_PSCF = 4
      write (20,*) "Doing PSCF Non-Linear Optical Properties"
   elseif (jobID == 207) then
      doOPTC_PSCF = 3
      write (20,*) "Doing PSCF Sigma(E)"
   elseif (jobID == 208) then
      doSYBD_PSCF = 1
      write (20,*) "Doing PSCF Symmetric Band Structure"
   elseif (jobID == 209) then
      doForce_PSCF = 1
      write (20,*) "Doing PSCF Force"
   elseif (jobID == 210) then
      doField_PSCF = 1
      write (20,*) "Doing PSCF Field"
   elseif (jobID == 211) then
      doMTOP_PSCF = 1
      write (20,*) "Doing PSCF Modern Polarization"
   elseif (jobID == 311) then
      doLoEn = 1
      write (20,*) "Doing Local Environment"

      ! loen is a standalone, geometry-only job and never shares a
      !   run with an SCF or post-SCF pass (DESIGN 5.10.2).  The test
      !   is on the doSCF/doPSCF decisions, NOT on the basis codes:
      !   when both codes arrive as zero, readBasisCodes substitutes
      !   a valid post-SCF code so the input readers have a slot to
      !   index, and that substitution must not be mistaken for a
      !   requested pass.  With a pass present, loen's parseInput
      !   would either index the per-type arrays with a zero basis
      !   code (SCF basis alone) or re-parse onto arrays the SCF
      !   pass left allocated (post-SCF basis after SCF).  imago.py
      !   refuses the combination first; this guard is for a
      !   hand-issued command line, and it sits here -- before any
      !   pass has run -- rather than in loen itself, which is
      !   reached only after the SCF and post-SCF blocks.
      if ((doSCF == 1) .or. (doPSCF == 1)) then
         write (20,*) 'loen is a standalone, geometry-only job and'
         write (20,*) 'cannot share a run with an SCF or post-SCF pass.'
         write (20,*) 'Give both basis codes as 0 (imago.py: -loen with'
         write (20,*) 'no -scf/-pscf, or -scf no -pscf no) and run the'
         write (20,*) 'other job separately.'
         call flush (20)
         stop 'readJobID: SCF or PSCF pass requested in a loen run'
      endif
   endif

end subroutine readJobID


subroutine initCLP

   ! Make sure that nothing funny is declared.
   implicit none

   ! Make sure that the first command line parameter to be read is #1.
   nextArg = 1

   ! Initialize all command line parameters to negative one.
   basisCode_SCF   = -1
   basisCode_PSCF  = -1
   doSCF           = -1
   doPSCF          = -1
   excitedQN_n     = -1
   excitedQN_l     = -1
   doDOS_SCF       = -1
   doBond_SCF      = -1
   doDIMO_SCF      = -1
   doOptc_SCF      = -1
   doSYBD_SCF      = -1
   doForce_SCF     = -1
   doField_SCF     = -1
   doMTOP_SCF      = -1
   doDOS_PSCF      = -1
   doBond_PSCF     = -1
   doDIMO_PSCF     = -1
   doOptc_PSCF     = -1
   doSYBD_PSCF     = -1
   doForce_PSCF    = -1
   doField_PSCF    = -1
   doMTOP_PSCF     = -1
   serialXYZ       = -1
   doLoEn          = -1

end subroutine initCLP


end module O_CommandLine

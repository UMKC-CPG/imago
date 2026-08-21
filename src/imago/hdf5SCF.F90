!! SPDX-License-Identifier: ECL-2.0
!! Copyright (c) 2026 Paul Rulis

module O_SCFHDF5

   ! Use the HDF5 module for HDF5 defined types (e.g. size_t and hid_t).
   use HDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define access
   public

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module data.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Declare the file ID.
   integer(hid_t) :: scf_fid

   ! Declare the property list for the scf file and its associated parameters.
   integer(hid_t)  :: scf_plid
   integer         :: mdc_nelmts  ! Meta-data cache num elements.
   integer(size_t) :: rdcc_nelmts ! Raw-data chunk cache num elements.
   integer(size_t) :: rdcc_nbytes ! Raw-data chunk cache size in bytes.
   real            :: rdcc_w0     ! Raw-data chunk cache weighting parameter.

   ! Declare the kPoint group ID.
   integer(hid_t) :: kPoint_gid

   ! Declare shared variables that are used for creating attributes that will
   !   track the completion of all datasets.
   integer(hid_t) :: attribInt_dsid ! Attribute dataspace.
   integer(hsize_t), dimension (1) :: attribIntDims ! Dataspace dimensionality

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Begin list of module subroutines.!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   contains

subroutine initHDF5_SCF (maxNumRayPoints, numStates)

   ! Use necessary modules.
   use O_TimeStamps

   ! Use the HDF5 module.
   use HDF5

   ! Use the subsection object modules for scf.
   use O_SCFIntegralsHDF5
   use O_SCFElecStatHDF5
   use O_SCFExchCorrHDF5
   use O_SCFEigValHDF5
   use O_SCFEigVecHDF5
   use O_SCFPotRhoHDF5
   use O_CommandLine, only: doSYBD_SCF, doMTOP_SCF
   use O_KPoints, only: kPointGroupName

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed dummy variables.
   integer, intent(in) :: maxNumRayPoints
   integer, intent(in) :: numStates

   ! Declare local variables.
   integer :: hdferr
   logical :: fileExists
   logical :: groupExists
   character*14 :: fileName
   character*17 :: kPointName

   ! Log the time we start to setup the SCF HDF5 files.
   call timeStampStart(6)

   ! Initialize the Fortran 90 HDF5 interface.
   call h5open_f(hdferr)
   if (hdferr < 0) stop 'Failed to open HDF library'

   ! Identify the file name of the hdf5 file we need to create/open.
   !   The construction lives in scfFileName (below) because the
   !   distributed integral stage reopens this file BY NAME on ranks
   !   that never execute this routine; one shared builder keeps the
   !   two spellings from drifting apart.
   call scfFileName (fileName)


   ! Create the property list for the scf hdf5 file and turn off
   !   chunk caching.
   call h5pcreate_f (H5P_FILE_ACCESS_F,scf_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create scf plid.'
   call h5pget_cache_f (scf_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,rdcc_w0,&
         & hdferr)
   if (hdferr /= 0) stop 'Failed to get scf plid cache settings.'
   call h5pset_cache_f (scf_plid,mdc_nelmts,0_size_t,0_size_t,rdcc_w0,hdferr)
   if (hdferr /= 0) stop 'Failed to set scf plid cache settings.'

   ! Create the name for the kPoint dependent data group.  The name
   !   records which k-point set (mesh, MTOP mesh, or band-structure
   !   path) produced the data; the one builder in O_KPoints serves
   !   this file and the PSCF file so the two forms cannot drift.
   call kPointGroupName (doSYBD_SCF, doMTOP_SCF, kPointName)

   ! Determine if an HDF5 file already exists for this calculation.
   inquire (file=fileName, exist=fileExists)

   ! If it does, then access the existing file. If not, then create one.
   if (fileExists .eqv. .true.) then
      ! We are continuing a previous calculation.

      ! Open the HDF5 file for reading / writing.
      call h5fopen_f (fileName,H5F_ACC_RDWR_F,scf_fid,hdferr,&
            & scf_plid)
      if (hdferr /= 0) stop 'Failed to open scf hdf5 file.'

      ! Check if a top-level group for the current kPoint set exists. If so,
      !   then we are continuing that kPoint set. Otherwise, we are starting
      !   a new kPoint set and will need to initialize it.
      call h5lexists_f(scf_fid,kPointName,groupExists,hdferr)

      ! If the group exists, then access the kPoint dependent data.
      if (groupExists .eqv. .true.) then

         ! Open the group.
         call h5gopen_f(scf_fid,kPointName,kPoint_gid,hdferr)

         ! Access the groups of the HDF5 file.
         call accessSCFIntegralHDF5 (kPoint_gid)
         call accessSCFEigVecHDF5 (kPoint_gid,numStates)
         call accessSCFEigValHDF5 (kPoint_gid,numStates)
      else

         ! All datasets will have an attached attribute logging that the
         !   calculation has successfully completed. (Checkpointing.) Thus,
         !   we need to create the shared attribute dataspace.
         attribIntDims(1) = 1
         call h5screate_simple_f (1,attribIntDims(1),attribInt_dsid,hdferr)
         if (hdferr /= 0) stop 'Failed to create the attribInt_dsid'

         ! Create the kPoint group that will hold all kPoint dependent results.
         call h5gcreate_f(scf_fid,kPointName,kPoint_gid,hdferr)

         ! Create the subgroups of the pscf hdf5 file. This must be done in this
         !   order due to dependencies on potPot_dsid and others.
         call initSCFIntegralHDF5 (kPoint_gid,attribInt_dsid,attribIntDims)
         call initSCFEigVecHDF5 (kPoint_gid,attribInt_dsid,attribIntDims,&
               & numStates)
         call initSCFEigValHDF5 (kPoint_gid,numStates)
      endif

      ! Access the kPoint independent groups of the HDF5 file.
      call accessSCFElecStatHDF5 (scf_fid)
      call accessSCFExchCorrHDF5 (scf_fid)
      call accessSCFPotRhoHDF5 (scf_fid)

   else
      ! We are starting a new calculation.

      ! Create the HDF5 file that will hold all the computed results. This
      !   uses the default file creation and file access properties.
      call h5fcreate_f (fileName,H5F_ACC_EXCL_F,scf_fid,hdferr,&
            & H5P_DEFAULT_F,scf_plid)
      if (hdferr /= 0) stop 'Failed to create scf hdf5 file.'

      ! All datasets will have an attached attribute logging that the
      !   calculation has successfully completed. (Checkpointing.) Thus,
      !   we need to create the shared attribute dataspace.
      attribIntDims(1) = 1
      call h5screate_simple_f (1,attribIntDims(1),attribInt_dsid,hdferr)
      if (hdferr /= 0) stop 'Failed to create the attribInt_dsid'

      ! Create the kPoint group that will hold all kPoint dependent results.
      call h5gcreate_f(scf_fid,kPointName,kPoint_gid,hdferr)

      ! Create the subgroups of the scf hdf5 file. This must be done in this
      !   order due to dependencies on potPot_dsid and others.
      call initSCFIntegralHDF5 (kPoint_gid,attribInt_dsid,attribIntDims)
      call initSCFEigVecHDF5 (kPoint_gid,attribInt_dsid,attribIntDims,numStates)
      call initSCFEigValHDF5 (kPoint_gid,numStates)
      call initSCFElecStatHDF5 (scf_fid,attribInt_dsid)
      call initSCFExchCorrHDF5 (scf_fid,attribInt_dsid,attribIntDims,&
            & maxNumRayPoints)
      call initSCFPotRhoHDF5 (scf_fid)
   endif


   ! Log the time we finish setting up the SCF HDF5 files.
   call timeStampEnd(6)

end subroutine initHDF5_SCF


subroutine closeHDF5_SCF

   ! Use the HDF5 module.
   use HDF5

   ! Use the subsection object modules for scf.
   use O_SCFIntegralsHDF5, only: closeSCFIntegralHDF5
   use O_SCFExchCorrHDF5,  only: closeSCFExchCorrHDF5
   use O_SCFElecStatHDF5,  only: closeSCFElecStatHDF5
   use O_SCFEigVecHDF5,  only: closeSCFEigVecHDF5
   use O_SCFEigValHDF5,  only: closeSCFEigValHDF5
   use O_SCFPotRhoHDF5,  only: closeSCFPotRhoHDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare local variables.
   integer :: hdferr

   ! Close access to all the subgroups and their parts.
   call closeSCFIntegralHDF5
   call closeSCFExchCorrHDF5
   call closeSCFElecStatHDF5
   call closeSCFEigVecHDF5
   call closeSCFEigValHDF5
   call closeSCFPotRhoHDF5

   ! Close the property list.
   call h5pclose_f (scf_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close scf_plid.'

   ! Close the file.
   call h5fclose_f (scf_fid,hdferr)
   if (hdferr /= 0) stop 'Failed to close scf_fid.'

   ! Close access to the HDF5 interface.
   call h5close_f (hdferr)
   if (hdferr /= 0) stop 'Failed to close the HDF5 interface.'

end subroutine closeHDF5_SCF


subroutine scfFileName (fileName)

   ! Construct the name of the SCF HDF5 file from the command-line
   !   state. This is a pure function of inputs every rank of a
   !   parallel run has parsed, so every rank can call it -- which is
   !   the point: writePotTermHDF5 below opens the file BY NAME on
   !   ranks that never ran initHDF5_SCF, and initHDF5_SCF itself
   !   calls this so the two spellings cannot drift apart.

   ! Use necessary modules.
   use O_CommandLine, only: excitedQN_n, excitedQN_l, basisCode_SCF

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed dummy variables.
   character*14, intent (out) :: fileName

   ! Declare local variables.
   character*2 :: edge

   ! Identify the edge name: the ground state, or the excited
   !   quantum-number pair (e.g. 1s, 2p) this calculation targets.
   if (excitedQN_n == 0) then
      write(edge,fmt="(a)") "gs"
   else
      if (excitedQN_l == 0) then
         write(edge,fmt="(i1,a1)") excitedQN_n, "s"
      elseif (excitedQN_l == 1) then
         write(edge,fmt="(i1,a1)") excitedQN_n, "p"
      elseif (excitedQN_l == 2) then
         write(edge,fmt="(i1,a1)") excitedQN_n, "d"
      elseif (excitedQN_l == 3) then
         write(edge,fmt="(i1,a1)") excitedQN_n, "f"
      elseif (excitedQN_l == 4) then
         write(edge,fmt="(i1,a1)") excitedQN_n, "g"
      endif
   endif

   ! Attach the basis-code suffix (minimal, full, or extended basis).
   if (basisCode_SCF == 1) then
      write(fileName,fmt="(a2,a12)") edge,"_scf-mb.hdf5"
   elseif (basisCode_SCF == 2) then
      write(fileName,fmt="(a2,a12)") edge,"_scf-fb.hdf5"
   elseif (basisCode_SCF == 3) then
      write(fileName,fmt="(a2,a12)") edge,"_scf-eb.hdf5"
   endif

end subroutine scfFileName


subroutine suspendHDF5_SCF

   ! Close the SCF file COMPLETELY so that the distributed term stage
   !   can use HDF5's own file lock as its write mutex (DESIGN 9.5;
   !   PSEUDOCODE 25.3): while the stage runs, the file is open by
   !   nobody between individual per-term writes, and any rank's
   !   openFileWithRetry can take its turn. Called on root only,
   !   after the root-only setup writes and the done-mask read.
   !
   ! The sections' close routines tear down the handles they track,
   !   but the handle population at this moment is configuration
   !   dependent (a stage skipped by restart or by rel == 0 leaves
   !   its attribute open where a stage that ran closed it), and
   !   HDF5 defers the REAL file close until the last handle dies --
   !   a single leaked handle would keep the lock held and starve
   !   every rank's retry loop. The closeAllFileObjects sweep is what
   !   turns this close into a guarantee. The HDF5 library interface
   !   itself stays open: resumeHDF5_SCF needs it moments later.

   ! Use necessary modules.
   use O_HDF5Utils, only: closeAllFileObjects

   ! Use the subsection object modules for scf.
   use O_SCFIntegralsHDF5, only: closeSCFIntegralHDF5
   use O_SCFExchCorrHDF5,  only: closeSCFExchCorrHDF5
   use O_SCFElecStatHDF5,  only: closeSCFElecStatHDF5
   use O_SCFEigVecHDF5,  only: closeSCFEigVecHDF5
   use O_SCFEigValHDF5,  only: closeSCFEigValHDF5
   use O_SCFPotRhoHDF5,  only: closeSCFPotRhoHDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Declare local variables.
   integer :: hdferr

   ! Close access to all the subgroups and their parts, exactly as
   !   the end-of-run teardown does.
   call closeSCFIntegralHDF5
   call closeSCFExchCorrHDF5
   call closeSCFElecStatHDF5
   call closeSCFEigVecHDF5
   call closeSCFEigValHDF5
   call closeSCFPotRhoHDF5

   ! Sweep whatever the section teardowns did not know about (the
   !   k-point group itself, and any stage-skipped status attribute
   !   that is still open).
   call closeAllFileObjects (scf_fid)

   ! Close the property list and the file. With the sweep done, this
   !   close truly releases the file and its lock.
   call h5pclose_f (scf_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to close scf_plid at suspend.'
   call h5fclose_f (scf_fid,hdferr)
   if (hdferr /= 0) stop 'Failed to close scf_fid at suspend.'

end subroutine suspendHDF5_SCF


subroutine resumeHDF5_SCF (numStates)

   ! Reopen the SCF file and rebuild every handle after the
   !   distributed term stage (PSEUDOCODE 25.3). Called on root only,
   !   and only after the stage-end barrier: an earlier reopen would
   !   hold the file lock against every rank still writing its terms.
   !   This is the SAME path the restart branch of initHDF5_SCF
   !   exercises when it finds an existing file and k-point group --
   !   the access routines re-allocate the handle arrays their close
   !   counterparts deallocated -- so nothing here is new machinery,
   !   only the proven close/access pair driven within one run.

   ! Use necessary modules.
   use O_CommandLine, only: doSYBD_SCF, doMTOP_SCF
   use O_KPoints, only: kPointGroupName

   ! Use the subsection object modules for scf.
   use O_SCFIntegralsHDF5, only: accessSCFIntegralHDF5
   use O_SCFExchCorrHDF5,  only: accessSCFExchCorrHDF5
   use O_SCFElecStatHDF5,  only: accessSCFElecStatHDF5
   use O_SCFEigVecHDF5,  only: accessSCFEigVecHDF5
   use O_SCFEigValHDF5,  only: accessSCFEigValHDF5
   use O_SCFPotRhoHDF5,  only: accessSCFPotRhoHDF5

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed dummy variables.
   integer, intent(in) :: numStates

   ! Declare local variables.
   integer :: hdferr
   character*14 :: fileName
   character*17 :: kPointName

   ! Recreate the cache-tuned property list exactly as initHDF5_SCF
   !   built it, then reopen the file and the k-point group by name.
   call h5pcreate_f (H5P_FILE_ACCESS_F,scf_plid,hdferr)
   if (hdferr /= 0) stop 'Failed to create scf plid at resume.'
   call h5pget_cache_f (scf_plid,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,&
         & rdcc_w0,hdferr)
   if (hdferr /= 0) stop 'Failed to get scf plid cache at resume.'
   call h5pset_cache_f (scf_plid,mdc_nelmts,0_size_t,0_size_t,&
         & rdcc_w0,hdferr)
   if (hdferr /= 0) stop 'Failed to set scf plid cache at resume.'

   call scfFileName (fileName)
   call h5fopen_f (fileName,H5F_ACC_RDWR_F,scf_fid,hdferr,scf_plid)
   if (hdferr /= 0) stop 'Failed to reopen scf hdf5 file at resume.'

   call kPointGroupName (doSYBD_SCF, doMTOP_SCF, kPointName)
   call h5gopen_f (scf_fid,kPointName,kPoint_gid,hdferr)
   if (hdferr /= 0) stop 'Failed to reopen kpoint group at resume.'

   ! Access the k-point dependent groups, then the file-level ones.
   call accessSCFIntegralHDF5 (kPoint_gid)
   call accessSCFEigVecHDF5 (kPoint_gid,numStates)
   call accessSCFEigValHDF5 (kPoint_gid,numStates)
   call accessSCFElecStatHDF5 (scf_fid)
   call accessSCFExchCorrHDF5 (scf_fid)
   call accessSCFPotRhoHDF5 (scf_fid)

end subroutine resumeHDF5_SCF


subroutine writePotTermHDF5 (termIndex, packedBuffer, packedVVDims,&
      & lockWaitSeconds, writeSeconds)

   ! Write ONE finished three-centre potential term -- every
   !   k-point's packed, orthogonalized matrix plus the completion
   !   attribute -- under the lock discipline of DESIGN 9.5 and
   !   PSEUDOCODE 25.5. The file is acquired through HDF5's own lock
   !   (openFileWithRetry), the term's objects are opened BY NAME
   !   using the parameters O_SCFIntegralsHDF5 exports beside their
   !   creation code, and everything is closed again so the next
   !   rank can take the lock. On one rank the open never contends
   !   and this degenerates to a plain open-write-close.
   !
   ! This routine lives HERE, with the other file-level operations,
   !   because it OPENS THE FILE: in the whole program only this
   !   module does that, and keeping that authority in one place is
   !   the module family's organizing rule (ARCHITECTURE 2). Every
   !   rank calls it -- root included, whose handle-holding view of
   !   the file is suspended for the whole stage -- and no rank may
   !   hold the file open between calls, or the other ranks' retry
   !   loops would starve.

   ! Use necessary modules.
   use O_Kinds
   use O_KPoints, only: numKPoints, kPointGroupName
   use O_CommandLine, only: doSYBD_SCF, doMTOP_SCF
   use O_HDF5Utils, only: openFileWithRetry
   use O_MPI, only: stopMPI
   use O_SCFIntegralsHDF5, only: intgGroupName, potOverlapGroupName,&
         & statusAttribName

   ! Make sure that no funny variables are defined.
   implicit none

   ! Define passed dummy variables.
   integer, intent (in) :: termIndex
   real (kind=double), dimension (:,:,:), intent (in) :: packedBuffer
   integer(hsize_t), dimension (2), intent (in) :: packedVVDims
   real (kind=double), intent (inout) :: lockWaitSeconds
   real (kind=double), intent (inout) :: writeSeconds

   ! Declare local variables.
   integer :: i
   integer :: hdferr
   integer(hid_t) :: fid            ! This open of the file.
   integer(hid_t) :: kPointOpen_gid ! The k-point set group.
   integer(hid_t) :: intg_gid       ! atomIntgGroup within it.
   integer(hid_t) :: potOL_gid      ! atomPotOverlap within that.
   integer(hid_t) :: term_gid       ! This term's own group.
   integer(hid_t) :: dsid           ! One k-point's dataset.
   integer(hid_t) :: aid            ! The term's status attribute.
   integer(hsize_t), dimension (1) :: attribIntDims
   character*17 :: kPointName
   character*14 :: fileName
   character*30 :: objectName
   integer :: clockGrant, clockDone, clockRate
   real (kind=double) :: waited

   ! The library interface may not be initialized on a worker rank,
   !   which never ran initHDF5_SCF. h5open_f is reference counted,
   !   so calling it again where it is already open is harmless.
   call h5open_f (hdferr)
   if (hdferr < 0) call stopMPI ('Failed to open HDF library.')

   ! Reconstruct the two names that locate this run's data, then
   !   take our turn at the file.
   call scfFileName (fileName)
   call kPointGroupName (doSYBD_SCF, doMTOP_SCF, kPointName)
   call openFileWithRetry (fileName, fid, waited)
   lockWaitSeconds = lockWaitSeconds + waited
   call system_clock (clockGrant, clockRate)

   ! Open the term's objects by name:
   !   /<kPointName>/atomIntgGroup/atomPotOverlap/<term>/<kpoint>
   !   with the term and k-point spelled as seven-digit numbers,
   !   exactly as initSCFIntegralHDF5 created them.
   call h5gopen_f (fid, kPointName, kPointOpen_gid, hdferr)
   if (hdferr /= 0) call stopMPI ('Failed to open kpoint group.')
   call h5gopen_f (kPointOpen_gid, intgGroupName, intg_gid, hdferr)
   if (hdferr /= 0) call stopMPI ('Failed to open atomIntgGroup.')
   call h5gopen_f (intg_gid, potOverlapGroupName, potOL_gid, hdferr)
   if (hdferr /= 0) call stopMPI ('Failed to open atomPotOverlap.')
   write (objectName,fmt="(i7.7)") termIndex
   call h5gopen_f (potOL_gid, trim (objectName), term_gid, hdferr)
   if (hdferr /= 0) call stopMPI ('Failed to open the term group.')

   ! Write every k-point's packed matrix for this term.
   do i = 1, numKPoints
      write (objectName,fmt="(i7.7)") i
      call h5dopen_f (term_gid, trim (objectName), dsid, hdferr)
      if (hdferr /= 0) then
         call stopMPI ('Failed to open a term dataset.')
      endif
      call h5dwrite_f (dsid, H5T_NATIVE_DOUBLE, packedBuffer(:,:,i),&
            & packedVVDims, hdferr)
      if (hdferr /= 0) then
         call stopMPI ('Failed to write a term dataset.')
      endif
      call h5dclose_f (dsid, hdferr)
      if (hdferr /= 0) then
         call stopMPI ('Failed to close a term dataset.')
      endif
   enddo

   ! Mark the term complete. The attribute is the unit of both
   !   restart and distribution (DESIGN 9.5), and it is written ONLY
   !   here, after the data: a rank killed at any earlier moment
   !   leaves the term undone, and the next run deals it again.
   attribIntDims(1) = 1
   call h5aopen_f (term_gid, statusAttribName, aid, hdferr)
   if (hdferr /= 0) call stopMPI ('Failed to open a term attribute.')
   call h5awrite_f (aid, H5T_NATIVE_INTEGER, 1, attribIntDims, hdferr)
   if (hdferr /= 0) call stopMPI ('Failed to write a term attribute.')
   call h5aclose_f (aid, hdferr)
   if (hdferr /= 0) call stopMPI ('Failed to close a term attribute.')

   ! Release the lock: close every object, then the file.
   call h5gclose_f (term_gid, hdferr)
   if (hdferr /= 0) call stopMPI ('Failed to close the term group.')
   call h5gclose_f (potOL_gid, hdferr)
   if (hdferr /= 0) call stopMPI ('Failed to close atomPotOverlap.')
   call h5gclose_f (intg_gid, hdferr)
   if (hdferr /= 0) call stopMPI ('Failed to close atomIntgGroup.')
   call h5gclose_f (kPointOpen_gid, hdferr)
   if (hdferr /= 0) call stopMPI ('Failed to close kpoint group.')
   call h5fclose_f (fid, hdferr)
   if (hdferr /= 0) call stopMPI ('Failed to close the scf file.')

   ! Account the write time separately from the lock wait: together
   !   they are the PA3 acceptance numbers (PSEUDOCODE 25.6) that
   !   answer DESIGN 9.5's unmeasured-write-share question.
   call system_clock (clockDone)
   writeSeconds = writeSeconds&
         & + real (clockDone - clockGrant, double)&
         & / real (clockRate, double)

end subroutine writePotTermHDF5


end module O_SCFHDF5

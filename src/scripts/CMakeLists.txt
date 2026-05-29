# CMake file for Imago set of scripts.

# Document the current location
message("At the set of scripts.")

set (SCRIPTS
   ab2pf.py angle_utils.py
   bond_analysis.py build_initial_potentials.py
   condense.py
   deadmd.py dump2skl.py
   element_data.py element_pct.py
   graspElems.py
   hdf5Comp.py
   imago.py initial_potential_db.py interFit.py
   lmp2skl.py
   make_veusz_graph.py makeBDB.py makePDOS.py makePotDB.py make_reactions.py
   makeSGDB.py makeSYBD.py makeinput.py mod_struct.py
   nanoTube.py
   org_radwavefn.py osrecurintg.py osrecurintgana.py osrecurintgnum.py
   osrecurintglib.py osrecurintg_makenum.py
   pdb2skl.py plot_deadmd.py plotgraph.py pot2plot.py processPOPTC.py
   runIsoAtoms.py
   skl2isaacs.py skl2pdb.py skl2vasp.py struct2skl.py structure_control.py
   unpackImagoDB.py
   vasp2skl.py viewCell.py
   xyz2skl.py)

set (RC_FILES
   bond_analysisrc.py
   collectSYBDrc.py condenserc.py
   dump2sklrc.py
   esviewrc.py
   hdf5Comprc.py
   graspElemsrc.py
   imagorc.py
   insert3cborc.py interFitrc.py
   makeBDBrc.py makeBONDrc.py makeFittedRhoVrc.py
   makeinputrc.py makePDOSrc.py makePotDBrc.py make_reactionsrc.py
   makeSYBDrc.py mod_structrc.py
   osrecurintg_makenumrc.py osrecurintgrc.py
   pdb2sklrc.py plotgraphrc.py porosityrc.py pot2plotrc.py processPOPTCrc.py
   runIsoAtomsrc.py
   skl2isaacsrc.py skl2pdbrc.py skl2vasprc.py struct2sklrc.py
   vasp2sklrc.py viewCellrc.py
   xyz2sklrc.py)

install (PROGRAMS ${SCRIPTS} DESTINATION bin)

# Install the kaleidoscope flight-runner package (DESIGN 6.2)
# as a directory under bin/.  Unlike the flat scripts above it is
# a Python *package*, so the whole directory is copied; it becomes
# importable as `kaleidoscope` by any script run out of bin/
# (bin/ is the script's sys.path[0]) and by a Parsl worker whose
# PYTHONPATH includes bin/.  Only the .py sources are installed --
# FILES_MATCHING excludes __pycache__ and any other stray files.
install (DIRECTORY kaleidoscope
         DESTINATION bin
         FILES_MATCHING PATTERN "*.py")

# SPECIAL PROCESS FOR RC FILES: #################################
# Install RC files only if they do not already exist, so that user
# edits made after the initial install are preserved across future
# "make install" runs.  To force an overwrite (e.g. when upstream
# defaults change), reconfigure with:
#   cmake -DIMAGO_FORCE_RC_INSTALL=ON ...
# and re-run "make install".
option(IMAGO_FORCE_RC_INSTALL
    "Overwrite existing RC files on install even if they already exist"
    OFF)

foreach(RC_FILE ${RC_FILES})
    if(IMAGO_FORCE_RC_INSTALL)
        install(PROGRAMS ${RC_FILE} DESTINATION ".imago")
    else()
        install(CODE "
            set(_src \"${CMAKE_CURRENT_SOURCE_DIR}/${RC_FILE}\")
            set(_dest
                \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/.imago\")
            set(_dst \"\${_dest}/${RC_FILE}\")
            if(NOT EXISTS \"\${_dst}\")
                file(INSTALL DESTINATION \"\${_dest}\"
                    TYPE PROGRAM
                    FILES \"\${_src}\")
            else()
                message(STATUS
                    \"Skipping existing RC file: \${_dst}\")
            endif()
        ")
    endif()
endforeach()

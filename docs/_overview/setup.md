---
title: Setup
order: 2
---
# Setting up Imago

Imago is intended to be cloned from github, built on the local machine, and ran from source. Currently there are no downloadable release builds availble. This page is intended to take you through the entire setup process. On this page is the general instructions for setting up and installing imago, but there are some scenario-specific instructions available under **Tutorials** if you so desire.

The following specific scenarios are available:
- [Compile and install for system-wide use](link to tutorial page)
- [Compile for a group](link to tutorial page)
- [Compile and install for a single user](link to tutorial page)
- [Blend group and single user options](link to tutorial page)

If you are following one of the above guides, this page serves to describe the big-picture idea core to each tutorial and may still be helpful.

A general overview of the following install process is:
- Setting up Compilers
- Cloning from Github
- User environment setup
- Compiling and install Imago
- Validating install
- Updating Imago

## Pre-requisites

The full requirements of the OLCAO package are as follows:

##### Compilers and build environment:
 
* A C compiler
* A Fortran compiler capable of compiling Fortran 90/95 (e.g., gfortran or ifort).
* The cmake (version >= 3.1.0) and GNU make programs.
    * The cmake executable may be `cmake` or it may be `cmake3` depending on your machine.

##### Perl and Perl libraries:

* Perl (version 5.8.8 or later)
* Not needed anymore: Inline::C
* Not needed anymore: Math::MatrixReal

##### Python and Python libraries

* Python (version 3)
* sympy
* numpy
* scipy
* pandas
* vedo
* matplotlib
* veusz (optional)
* perl-app-cpanminus (optional)

##### Other requirements:

* Linear algebra libraries (lapack and blas, mkl, atlas, etc.)
* HDF5 (including the Fortran interface and deflate filter)

Additional guidance on setting these up can be found [here](link to prerequisites page under tutorials) and a basic approach can be found in the next section.

***IMPORTANT:*** Make sure that the fortran compiler used to compile your HDF5 install is the same one you'll be using. Different compilers are known to cause problems

## Setting up Compilers

We will assume that you have access to a Fortran and C compiler. Perhaps you are using a system wide install of gfortran and gcc or you have installed a compiler using conda. Regardless, you need to make sure that you have HDF5 correctly wrapping your desired compiler. Use `h5fc -showconfig` to see the HDF5 configuration. In the output, look to ensure that a "deflate" filter was turned on. Then, identify which compilers (Fortran and C) were used to compile HDF5. It is vital that the compiler that was used to build HDF5 be the same kind that you have on your system. E.g., if HDF5 was built with gfortan, but your system has ifort (or vice-versa), then you will probably have a problem because these two Fortran compilers do not seem to play nice with each other's module files.

If you don't have HDF5 installed, then you will need to do that either by asking your system administrators to do it, or by using conda, or by doing a self install. In any case, you will need to make sure as above that the "deflate" filter is on and that the Fortran interface is included. See [here](https://github.com/UMKC-CPG/olcao/wiki/HDF5-Installation-Guide) to learn how to install HDF5 yourself in your $HOME.

For Perl, the version does not matter too much as long as it isn't really old. The hard part may be setting up the necessary Perl modules. Look in the above scenario links to learn how to install the necessary Perl modules for that scenario.

For Python, use version 3. The Python modules are usually easy to install. Look in the above scenario links to learn how to install the necessary Python modules for that scenario.

## Cloning from github

Once you have all required dependencies it's time to begin installation. If you already know how, use your preferred method to clone the Imago repository and move on to the next section. Otherwise, keep reading.

Steps:
- use git to clone the repository to your desired location (Ideally: $HOME). You may use either https:// (first) methods or ssh (second) methodsto do so.
   - `git clone https://github.com/UMKC-CPG/imago.git`
   - `git clone git@github.com:UMKC-CPG/imago.git`
- By default this will create an olcao/ directory wherever you run the clone command. If you wish to have a custom-named directiory, instead do:
   - `git clone https://github.com/UMKC-CPG/imago.git <directoryName>` or 
   - `git clone git@github.com:UMKC-CPG/imago.git <directoryName>`
   - Where `<directoryName>` is the desired path to the install directory. Usually it is:
      - $HOME/my/olacao/install/dir
      - $HOME/olcao (default)

## Setting up User Environment

Regardless of your unique scenario, it is highly recommended for each user to set up an "Imago Friendly" environment. This makes setup and use much easier on the users.

To set up the user environment perform the following steps:
- cd into the hidden .imago directory in the Repository
- Edit the imagorc (NOT imagorc.py) file and edit the environment variables such that they correctly point to the neccesary directories.
   - Make sure that IMAGO\_BRANCH is the path to the repository you used when cloning imago
   - Make IMAGO\_DIR the path from system root to the start of your IMAGO\_BRANCH path. $HOME/$IMAGO\_BRANCH by default
   - Note: IMAGO\_RC can be different from IMAGO\_DIR and will be where Imago will look for the imagorc.py file for defaults loading
- Edit your .bashrc or bash profile file to include the following:
   - Ask it to source the `imagorc` file you just edited.
   - Include the following aliases:  
     + `alias cmakedebug='cmake ../.. -DCMAKE_BUILD_TYPE=DEBUG'`  
     + `alias cmakerelease='cmake ../.. -DCMAKE_BUILD_TYPE=RELEASE'`  
     + `alias cdrelease='cd $IMAGO_DIR/build/release/'`  
     + `alias cddebug='cd $IMAGO_DIR/build/debug/'`  
     + `alias cdsrc='cd $IMAGO_DIR/src/olcao/'`  
     + `alias cdjobs='cd $IMAGO_DIR/jobs/'`
     + `alias interactive='srun -p interactive --mem=8G --pty /bin/bash'`
     + `alias interactive_login='srun -p interactive --mem=20G --pty /bin/bash --login'`
     + `alias interactive_x11='srun -p interactive --x11=first --mem=20G --pty /bin/bash --login'`
     + `alias la='ls -latr'`
     + `alias shist='sacct --starttime $(date -d yesterday +%D-%R) --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss,MaxVMSize,nnodes,ncpus,nodelist'`

After you finish editing your bashrc or bash profile, either source it or restart the terminal. This will set up your environment.

## Compiling and Installing Imago

The first time you compile Imago:
  * You need to make a compilation directory first:
    - `cd $IMAGO_DIR`
    - `mkdir build`
    - `cd build`
    - `mkdir release`
    - `cd release`
  * Then run cmake, compile/install, and unpack Imago databases:
    - `cmakerelease`
    - `make install`
    - `unpackImagoDB.py`

Later, if you are just editing existing files/scripts:
  * Go to the compile directory.
    - `cdrelease`
  * Compile and install
    - `make install`

Obviously, if you do more than just edit existing scripts and source files, then you will need to re-run cmake and may need to reinstall the Imago databases.

This is the basic outline of what you need to do to install the Imago package. See the above scenarios and machine specific instructions for more details. Happy computing!

## Validating the Imago Installation

* If all of the scripts were installed correctly, then the last step of running `unpackImagoDB.py` should spend a bit of time (5-15 seconds or so) and then return with no errors or other output.
* Go to the check directory `cd $IMAGO_DIR/check` and run the `runCheck` script. If all goes well, you should see a report printed to the screen indicating success of various test calculations. Success in this case means that a calculation you are running now produces output that matches with previously computed data.
* If things are not working well at this point you should probably contact your system administrator to save yourself time. Consider providing your administrator with your environment `env`, a copy-paste of any errors that were dumped to the screen, the `runtime` file if present, the `intermediate/fort.20` file if present, and the directory where you are doing this work when you ask for help.

### Updating OLCAO

If you have used git to clone the olcao files directly from this repository, the task of keeping your install up to date is very simple if you are not doing any development yourself. Do the following steps:

  * Change to your OLCAO root directory with `cd $IMAGO_DIR`.
  * Make sure that your remote directory is the UMKC-CPG repository (this github repository) by typing `git remote -v` after you . This should show something like `origin git@github.com/UMKC-CPG/imago.git`. If it does not, then you will need to change your remote directory using something like `git remote add imago github.com/UMKC-CPG/imago.git`.
  * An assumption is that you are updating your local master branch which you can get to with `git checkout master`.
  * If you know for certain that you want to merge the github version into your current branch you can use `git pull`. Proceed to the last step.
  * If you just want to see the changes you can use `git fetch origin`. This should show you all the changes (if any).
  * If you like what you see, then you can then merge the origin/master into your own local master with `git merge master origin/master`.
  * If there are changes and you wish to install them, go to the "build/release" directory with `cdrelease`, and run `make install`. This will run cmake again if needed and then compile and install.

You are now up to date!


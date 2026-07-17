---
layout: default
title: Imago Workflow
order: 3
---

# Imago Workflow

In general performing calculations in Imago follows these steps:
- Source atomic structure file (cif, pbd, skl, etc)
- Convert structure file to [skeleton](link to skeleton file page) file (skl) if it isn't one
   - cif2skl.py, pdb2skl.py, etc
- Create a job directory for the structure and place the skeleton file in the directory
- Make physics-informed decisions about what settings to use for the calculations
- Run makeinput.py with your desired settings on the skeleton file
- Run the imago.py script either by:
   - Editing then submitting the generated slurm file. Or
   - Running your desired calculations in the commandline using imago.py
- Post processing and analysis

Additionally, while currently in development, kaleidoscope will allow you to automate this process. Once ready, you can create a "builder" to create a large set of different inputfiles (the makeinput.py step) and then run them while keeping track of any that fail.

## Sourcing Atomic Structure Files

There are several external programs and databases that Imago has built-in support for using the conversion scripts. Supported external filetypes and their usual sources are given in a table at the bottom of this section. All conversion scripts use the naming convention [extension]2skl.py. Where [extension] is replaced by a filetype identifier (often the file's extenstion)

|Filetype                                 |Most common srouces             | [extension]  |
|:----------------------------------------|:-------------------------------|:------------:|
|Crystallographic Information Framework   |Crystallography Open Database   |cif           |
|LAMMPS dump file                         |LAMMPS                          |dump or lmp   |
|Protein DataBank                         |RCSB                            |pdb           |
|VASP position file                       |VASP                            |vasp          |
|Imago Structure data file                |Imago                           |struct        |
|xyz direct space file                    |N/A                             |xyz           |

### Crystallography Open Database

A very comprehensive source of material atomic structures is the Crystallography Open Database. Within Imago there is a script called cod\_fish.py. It can make api calls and download cif files directly for you. Additionally, you can define search parameters. See the script documentation for usage instructions [here](link to the cod_fish.py script documentation):

## Converting to Skeleton Files

The table in the section on sourcing atomic structure files has a column for [extension]. Simply replace [extension] with the value in that column in [extension]2skl.py to call the correct script. The cannonical command line usage is:
`[extension]2skl.py -i path/to/file -o path/to/output/location`

NOTE: there are 2 exceptions -- cif2skl and dump2skl. Their unique commandline usages are listed below.
- cif2skl.py:
   - `cif2skl.py [path/to/input.cif] [path/to/output.skl]`
- dump2skl.py:
   - `dump2skl.py -d [path/to/LAMMPS/dump/file] -a [path/to/LAMMPS/data/file] <-t [timestep] | -f [framenumber]>`

Other options and detailed documentation can be found in each scripts individual documentation page under Reference Manual and Script Catalog.

## Making the Inputs

Use the makeinput.py documentation to familiarize yourself with the available options; and use the skeleton file documentation page to familiarize yourself with the structure of skeleton files.

Once you have an idea what your doing (even if it's just running defaults), it is recommended you make a jobs/ directory in which to store all your input files. In this directory--or wherever you put your inputs--you will create a directory for each skeleton file you wish to use as input. Enter this directory and in it place your skeleton file and name it imago.skl for ease. Run makeinput.py with the options you wish to generate all the necessary inputs. This will include a default slurm file which you can edit if your machine uses slurm submissions.

## Running Imago

Now that you have your inputs ready, it's time to run the calculations you want. What calculations are available and a description of each and their outputs can be found on the imago.py documentation page.

Either edit the slurm to inlcude, line by line, all the calculations you wish to run one at a time. imago.py can only run one type of calculation at a time. Or, you can run jobs using the commandline one at a time. Again, one run of imago.py can only do one calculation. It will place all output files in the job directory you run it from. This is another reason why it's good to have one directory per skeleton file for organization.

## Post-processing

While some files are intended to be human readable, many require some amount of post-processing to make any useful products. There are several scripts that come with Imago to help you do this and quickly produce usable items. Covering them all is out of the scope of this page, thus if you are interested, direct yourself to the [post-processing page](link to post-processing page) to find a list of scripts and example use cases. Additionally, some information about the post-processing nature of some of the calculations can be found on the imago.py documentation page.

After this, Imago's role in your current analysis is over. Take what you get and do good science.


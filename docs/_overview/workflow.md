---
layout: default
title: Imago Workflow
nav_order: 3
---

# Imago Workflow

In general, performing a calculation in Imago follows these steps:

1. Obtain an atomic structure file (CIF, PDB, a LAMMPS dump, a VASP
   file, or a skeleton written by hand).
2. Convert it to a **skeleton** file (`imago.skl`) if it is not one
   already, with the matching `<extension>2skl.py` script.
3. Create a job directory for the structure and place the skeleton
   in it.
4. Make the physics-informed decisions about the calculation: the
   k-point sampling, the basis, the exchange-correlation functional,
   whether a band path is wanted.
5. Run `makeinput.py` with those settings in the job directory.
6. Run `imago.py`, either by submitting the generated batch script
   or directly on the command line.
7. Post-process and analyse.

The [first calculation tutorial](../tutorials/first_calculation.html)
performs all seven steps on silicon with every command shown and
every output file explained. This page is the map; that page is the
walk.

## Sourcing Atomic Structure Files

Imago has built-in support for several external file types through
its conversion scripts. All of them follow the naming convention
`<extension>2skl.py`, where `<extension>` identifies the file type,
usually by its file extension.

|Filetype                                 |Most common sources             | [extension]  | Script |
|:----------------------------------------|:-------------------------------|:------------:|:------:|
|Crystallographic Information Framework   |Crystallography Open Database   |cif           |[cif2skl.py](../reference_manual/script_cat/cif2skl.html)|
|LAMMPS dump file                         |LAMMPS                          |dump          |[dump2skl.py](../reference_manual/script_cat/dump2skl.html)|
|Protein DataBank                         |RCSB                            |pdb           |[pdb2skl.py](../reference_manual/script_cat/pdb2skl.html)|
|VASP position file                       |VASP                            |vasp          |[vasp2skl.py](../reference_manual/script_cat/vasp2skl.html)|
|Imago Structure data file                |Imago                           |struct        |[struct2skl.py](../reference_manual/script_cat/struct2skl.html)|
|xyz direct space file                    |N/A                             |xyz           |[xyz2skl.py](../reference_manual/script_cat/xyz2skl.html)|

### Crystallography Open Database

A very comprehensive source of material atomic structures is the
Crystallography Open Database. Imago's
[cod_fish.py](../reference_manual/script_cat/cod_fish.html) searches
it by composition and downloads CIF files directly, pinning the
database revision so that a structure can be fetched again exactly.

## Converting to Skeleton Files

For the PDB, VASP, Imago-structure and XYZ converters the canonical
command line is:

```
<extension>2skl.py -i path/to/file -o path/to/output.skl
```

Two converters take their inputs differently:

- `cif2skl.py path/to/input.cif [path/to/output.skl]`, which also
  resolves and preserves the crystal's space group;
- `dump2skl.py -d path/to/dump -a path/to/data <-t timestep | -f frame>`,
  which needs the LAMMPS data file to identify elements and a frame
  to extract.

Other options and detailed documentation are on each script's page
in the Script Catalog.

The skeleton format itself is small: a title block, the cell
parameters, the atoms of the asymmetric unit in fractional or
Cartesian coordinates, a space-group token, a supercell line, and
`prim` or `full`. The
[first calculation tutorial](../tutorials/first_calculation.html#read-the-skeleton-file)
explains each line using the shipped silicon skeleton, and
`$IMAGO_DIR/skl/` holds four worked examples: silicon, diamond,
graphite and quartz.

## Making the Inputs

Use the [makeinput.py](../reference_manual/script_cat/makeinput.html)
page to familiarize yourself with the available options. Once you
have an idea what you are doing (even if it is just running the
defaults), it is recommended that you make a `jobs/` directory in
which to store all your calculations. In it, create a directory for
each skeleton file you wish to run, place the skeleton there named
`imago.skl`, and run `makeinput.py` with the options you wish. This
generates the `inputs/` directory the engine reads, a `summary` of
the system, and a batch submission script you can edit if your
machine uses SLURM.

Two options deserve attention on every run: the k-point sampling
(`-kp` or `-kpd`), whose default of the Gamma point alone suits a
large amorphous cell but not a small crystal; and `-slurm`, whose
default partition and account are the authors' own.

## Running Imago

Now that your inputs are ready, run the calculations you want. What
calculations are available, and a description of each with its
outputs, is on the [imago.py](../reference_manual/script_cat/imago.html)
page.

One run of `imago.py` performs one calculation type, so either edit
the batch script to run the calculations you want one after another,
or run them from the command line one at a time. All output files
land in the job directory you run from, which is another reason to
keep one directory per skeleton. The SCF is performed once and
reused: after the first calculation in a directory, any later
calculation that needs the same basis skips straight to its own
property pass.

## Post-processing

Some output files are directly plottable; many require
post-processing to make useful products. The
[post-processing section of the Script Catalog](../reference_manual/script_cat/index.html#post-processing)
lists the scripts that do this, with
[plotgraph.py](../reference_manual/script_cat/plotgraph.html) as the
default plotter for densities of states, optical spectra and band
structures. Some information about the post-processing nature of
each calculation is also on the imago.py page.

After this, Imago's role in your current analysis is over. Take what
you get and do good science.

## Automating many calculations

The `kaleidoscope` package dispatches sweeps of calculations across
a cluster, with run reuse and a historical guidance database that
predicts the k-point density a new structure needs. Today it is
driven through
[build_initial_potentials.py](../reference_manual/script_cat/build_initial_potentials.html)
for the curated reference-structure database and is not yet a
general command-line tool for user campaigns.

---
layout: default
title: skl2vasp.py
script: skl2vasp.py
description: 'Converts an imago.skl skeleton file into VASP input. At present only the POSCAR file is written; the options for pseudopotential type, GW, Gamma-point sampling, job type and cell type are accepted and stored for the POTCAR, KPOINTS and INCAR generation that is not yet ported.'
parent: Script Catalog
---

# skl2vasp

Takes an Imago skeleton file and produces input files for a VASP
calculation, with a few intelligent options for choosing parameters.
The usual purpose is to relax a structure in VASP before bringing it
back into Imago with [vasp2skl.py](vasp2skl.html).

**Currently only the POSCAR file is generated.** The original Perl
script also produced POTCAR, KPOINTS and INCAR files, but those parts
are not yet ported to Python. The command-line options that governed
them (potential type and sub-type, the GW flag, the Gamma k-point
flag, the job type and the cell type) are parsed and stored so that
the interface is already in place for when those files are
implemented, but today they have no effect on the output.

**Defaults.** Every option's default value is read from the resource
control file `skl2vasprc.py` in the directory named by `$IMAGO_RC`.
Edit that file to change the defaults for every run; pass options on
the command line to override them for one run.

## Usage

**Commandline:**

```
skl2vasp.py [-j JOBTYPE] [-p POTTYPE [-s SUBTYPE]] [--gw] [--gamma]
            [-t CELLTYPE] [-i INPUT] [--help]
```

**Variable Table:**

| Variable   | Accepted Values                                         | Description |
|:----------:|:-------------------------------------------------------:|:------------|
|JOBTYPE     |relaxfull, relaxion1, relaxion2, relaxvol, cij, phonon, gpt | Type of VASP calculation to set up. Chooses optimal parameters for some of the input files. |
|POTTYPE     |potLDA, potGGA, potpawLDA, potpawGGA, potpawPBE, potpawLDA5x, potpawPBE5x | Base type of pseudopotential. See the table below. Do not replace the `x` with a number. |
|SUBTYPE     |s, h, pv, sv, d                                          | Sub-type of the potential. Many are valid only for particular elements or base types; it is your responsibility to know what you are asking for. |
|CELLTYPE    |tric, mono, ortho, tet, trig, hex, sc, bcc, fcc          | Bravais lattice type, overriding the auto-detection from the space group. |
|INPUT       |String; file path                                        | Imago skeleton input file. |

**Option Table:**

| Option          | Outputs      | Defaults                        | Description |
|:---------------:|:------------:|:-------------------------------:|:------------|
|-j / --job       |N/A           |relaxfull                        | Type of calculation to set up. Stored for future INCAR generation. |
|-p / --pot       |N/A           |potpawPBE                        | Base pseudopotential type. Stored for future POTCAR generation. |
|-s / --subpot    |N/A           |None (standard potential)        | Potential sub-type. Stored for future POTCAR generation. |
|--gw             |N/A           |Off                              | Append the GW sub-sub-type to the potential sub-type. Stored. |
|--gamma          |N/A           |Off                              | Use a Gamma-centred k-point instead of a general Monkhorst-Pack mesh. Stored for future KPOINTS generation. |
|-t / --celltype  |N/A           |Auto-detect from the skeleton    | Override the Bravais lattice cell type. |
|-i / --input     |N/A           |imago.skl                        | Imago skeleton input file. |
|--help           |N/A           |N/A                              | Print usage information and exit. |

## Pseudopotential Types

| POTTYPE       | Meaning |
|:-------------:|:--------|
|potLDA         | Ultrasoft LDA |
|potGGA         | Ultrasoft GGA |
|potpawLDA      | PAW LDA (VASP 4.x) |
|potpawGGA      | PAW GGA (VASP 4.x) |
|potpawPBE      | PAW PBE (VASP 4.x) |
|potpawLDA5x    | PAW LDA (VASP 5.x) |
|potpawPBE5x    | PAW PBE (VASP 5.x) |

| SUBTYPE  | Meaning |
|:--------:|:--------|
|s         | Soft |
|h         | Hard |
|pv        | Semi-core p states treated as valence |
|sv        | Semi-core s and p states treated as valence |
|d         | Semi-core d states treated as valence |

## Output

A `POSCAR` file in the current directory, with the atoms of each
element grouped together as VASP requires.

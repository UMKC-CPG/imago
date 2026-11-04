---
layout: default
title: dump2skl.py
script: dump2skl.py
description: 'Converts one frame of a LAMMPS dump file into an imago.skl skeleton file, using the LAMMPS data file to identify elements by their masses. The frame is chosen by frame number or by timestep.'
parent: Script Catalog
---

# dump2skl

Takes a LAMMPS coordinate dump file and its associated data file and
writes an `imago.skl` skeleton file for one chosen timestep. This is
the usual route from a molecular-dynamics trajectory (for example an
amorphous model produced by melt-and-quench) into an Imago
calculation.

**Choosing the frame.** You must specify either a timestep number
(`-t`) or a frame number (`-f`). Timesteps are the values listed in
the dump file under `ITEM: TIMESTEP`. Frame numbers are integers from
1 through the number of frames printed in the file. A negative frame
number counts backwards from the final frame, so `-f -1` means the
last frame.

**Identifying elements.** A LAMMPS dump file labels atoms only by a
numbered atom type. The script reads the `Masses` section of the data
file and compares each type's mass against the atomic masses in
Imago's element database; the element whose mass is closest is
assigned to that type.

**Defaults.** Every option's default value is read from the resource
control file `dump2sklrc.py` in the directory named by `$IMAGO_RC`.
Edit that file to change the defaults for every run; pass options on
the command line to override them for one run.

## Usage

**Commandline:**

```
dump2skl.py -d DUMPFILE -a DATAFILE (-f FRAME | -t TIMESTEP)
            [-n NAME] [--np] [--help]
```

**Variable Table:**

| Variable   | Accepted Values         | Description |
|:----------:|:-----------------------:|:------------|
|DUMPFILE    |String; file path        | Path to the LAMMPS dump file. |
|DATAFILE    |String; file path        | Path to the LAMMPS data file. |
|FRAME       |integer; 1 to N or negative | Frame number to extract. Negative values count from the end. |
|TIMESTEP    |integer                  | Timestep number matching an `ITEM: TIMESTEP` value in the dump file. |
|NAME        |String                   | System name written to the skeleton file header. |

**Option Table:**

| Option                | Outputs        | Defaults       | Description |
|:---------------------:|:--------------:|:--------------:|:------------|
|-d / --dump            |N/A             |None; required  | Path to the LAMMPS dump file. |
|-a / --data            |N/A             |None; required  | Path to the LAMMPS data file. |
|-f / --frame           |N/A             |Not set         | Frame number (1 to N, or negative to count from the end). Give either this or `-t`. |
|-t / --timestep        |N/A             |Not set         | Timestep number matching an `ITEM: TIMESTEP` value in the dump file. Give either this or `-f`. |
|-n / --name            |N/A             |name            | Name for the system, written to the skeleton file header. |
|--np / --nonPeriodic   |N/A             |Off             | Add 100 Angstroms of padding on each side of the cell so that the system behaves as non-periodic. See below. |
|--help                 |N/A             |N/A             | Print usage information and exit. |

## Option Details

### Non-periodic Systems

> --np / --nonPeriodic

Imago always applies periodic boundary conditions. A LAMMPS dump of a
bulk amorphous model is already periodic, so no padding is wanted and
the default is correct. For an isolated cluster or molecule that
happened to be simulated in a large LAMMPS box, `--np` adds 100
Angstroms of empty space on every side of the cell so that periodic
images cannot interact. Note that the short form is `--np` with two
dashes, because a single-dash `-np` would collide with `-n`.

## What the Script Does

1. Reads the LAMMPS data file to determine which chemical element
   corresponds to each LAMMPS atom type, by comparing masses.
2. Scans the dump file to count the frames and record their timestep
   values.
3. Resolves the frame or timestep you asked for.
4. Reads the box geometry and atom coordinates for that frame.
5. Processes the coordinates: atoms that drifted outside the box are
   wrapped back in, the optional padding is added, everything is
   shifted so the box starts at the origin, and the Cartesian
   coordinates are converted to fractional.
6. Writes the `imago.skl` file.

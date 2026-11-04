---
layout: default
title: unpackImagoDB.py
script: unpackImagoDB.py
description: 'One-time post-install step. Unpacks the five database tarballs shipped with Imago (atomic basis sets, atomic potentials, atomic SCF precursor data, space groups, and symmetric band k-point paths) in $IMAGO_DATA and regenerates the space group files.'
parent: Script Catalog
---

# unpackImagoDB

Unpacks the database files that come with the Imago package. The
databases are stored as compressed tarballs in the directory named
by `$IMAGO_DATA` and must be unpacked once before Imago can be used.
The build system installs this script to `bin/`; running it is the
last step of the installation described on the
[setup page](../../overview/setup.html).

After unpacking the space group database the script also runs the
`remake` script inside the `spaceDB/` directory, which regenerates
the individual space group files and their numeric symbolic links
from the packed representation (see [makeSGDB.py](makeSGDB.html)).

## Usage

**Commandline:**

```
unpackImagoDB.py [--help]
```

The script takes no arguments. It relies entirely on `$IMAGO_DATA`
to locate the tarballs, and it stops with an error if that variable
is not set.

## Databases Unpacked

| Tarball            | Directory      | Contents |
|:-------------------|:---------------|:---------|
|atomicBDB.tgz       |atomicBDB/      | Atomic basis set database, one directory per element (built by [makeBDB.py](makeBDB.html)). |
|atomicPDB.tgz       |atomicPDB/      | Atomic potential database (built by [makePotDB.py](makePotDB.html)). |
|precursorDB.tgz     |precursorDB/    | Precursor data for the atomic SCF program. |
|spaceDB.tgz         |spaceDB/        | Space group symmetry operations, one file per setting. |
|sybdDB.tgz          |sybdDB/         | Symmetric band structure k-point path definitions, one file per lattice type, used by `makeinput.py -sybdpath`. |

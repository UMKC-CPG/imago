---
layout: default
title: osrecurintg_makenum.py
script: osrecurintg_makenum.py
description: 'Derives the integrand equations used for numerical, brute-force evaluation of the multicentre Gaussian integrals, and prints them as Fortran code chunks for pasting into the osrecurintg numerical module. Separated from osrecurintg.py because the derivation is slow and need not be repeated on every run.'
parent: Script Catalog
---

# osrecurintg_makenum

Produces the equations used to compute **numerical** solutions of the
multicentre Gaussian-type-orbital integrals: the explicit integrands
that a brute-force quadrature on a mesh evaluates, against which the
analytic Obara-Saika solutions from [osrecurintg.py](osrecurintg.html)
are tested. It prints chunks of Fortran code that are copied into
`osrecurintgnum.py`, the module that builds the test program.

The derivation is separated from the program generation because it is
somewhat time-consuming, and once the equations have been determined
there is no reason to recreate them on every test or production run.
Like osrecurintg, this is a developer's tool run only when an
integral type is added.

**A caution from the script itself.** The output equations may be
copied directly, but that is not necessarily the best way to use
them. They can be separated into subcomponents organised into nested
loops to increase efficiency, and depending on the integration mesh
the numerical computation can be slow, so that reorganisation is
worthwhile. It is essential that no errors are introduced when doing
so; if a test can be performed to check that, perform it.

**Defaults.** Every integral type's default switch is read from
`osrecurintg_makenumrc.py` in the directory named by `$IMAGO_RC`.

## Usage

**Commandline:**

```
osrecurintg_makenum.py [-o] [-k] [-n] [-e] [-m] [-mv] [-d] [-dk]
                       [-dn] [-de] [-ko] [--help]
```

**Option Table:**

| Option               | Defaults | Description |
|:--------------------:|:--------:|:------------|
|-o / --overlap        |Off       | Two-centre overlap. |
|-k / --kinetic        |Off       | Two-centre kinetic energy. |
|-n / --nuclear        |Off       | Three-centre nuclear attraction. |
|-e / --electron       |Off       | Three-centre overlap with an s-type third centre. |
|-m / --momentum       |Off       | Two-centre momentum. |
|-mv / --massvel       |Off       | Mass velocity. |
|-d / --dipole         |Off       | Dipole moment. |
|-dk / --dkinetic      |Off       | Kinetic energy derivative. |
|-dn / --dnuclear      |Off       | Nuclear potential derivative. |
|-de / --delectron     |Off       | Electronic potential derivative. |
|-ko / --KOverlap      |Off       | Two-centre overlap with a plane wave. |
|--help                |N/A       | Print usage information and exit. |

## Output

`temp_Fortran_code.f90` in the current directory, holding the
generated code chunks for the requested integral types.

---
layout: default
title: osrecurintg.py
script: osrecurintg.py
description: 'Code generator for the Obara-Saika recursive evaluation of one-electron Gaussian integrals. Writes the Fortran source (gaussIntegrals.f90, optionally vectorised) that the Imago engine compiles for overlap, kinetic, nuclear, electronic, momentum, mass-velocity, dipole, derivative and plane-wave integrals, or a self-checking test program comparing analytic against numerical integration.'
parent: Script Catalog
---

# osrecurintg

A **code generator**. The integrals between Gaussian-type orbitals
that Imago needs, overlap, kinetic energy, nuclear attraction and the
rest, have closed-form solutions obtained by the Obara-Saika
recursion, but writing those solutions out by hand for every
combination of angular momenta up to f or g is error-prone and
enormous. This program performs the recursion symbolically and
writes the result as Fortran source. Developers run it once when an
integral type or the maximum angular momentum changes, and commit the
generated file into the engine; users of Imago never need it.

It has two modes:

* **Test mode** (the default) writes `gaussIntgTest.f90`, a
  standalone program that evaluates each requested integral both
  analytically, from the recursion, and numerically, by brute-force
  quadrature on a mesh, and compares them. This is how a new integral
  type is validated.
* **Production mode** (`-p`) writes `gaussIntegrals.f90`, or
  `gaussIntegrals.vec.f90` with `-v`, containing only the analytic
  solutions in the form the engine compiles.

Three sibling modules do the work and are imported, never run:
`osrecurintgana.py` produces the analytic solutions and their Fortran
rendering, `osrecurintgnum.py` the numerical ones for the test
program, and `osrecurintglib.py` the shared Fortran-printing
routines. [osrecurintg_makenum.py](osrecurintg_makenum.html) is the
separate, slow step that derives the numerical-integrand equations
those modules embed.

**Defaults.** Every option's default value is read from
`osrecurintgrc.py` in the directory named by `$IMAGO_RC`, where each
integral type has its own switch.

## Usage

**Commandline:**

```
osrecurintg.py [-p] [-v] [-l LMAX]
               [-o] [-k] [-n] [-e] [-m] [-mv] [-d] [-dk]
               [-dncb] [-dnbb] [-dnbc] [-decb] [-debb] [-debc] [-ko]
               [-a] [--help]
```

**Option Table:**

| Option                  | Defaults | Description |
|:-----------------------:|:--------:|:------------|
|-p / --production        |Off (test code) | Produce code for the production engine instead of the test program. |
|-v / --vectorize         |Off       | Vectorise the integrals. Ignored unless `-p` is also given. |
|-l / --max_lam           |3 (f)     | Highest angular momentum to produce equations for: 0 = s, 1 = p, 2 = d, 3 = f, 4 = g. |
|-o / --overlap           |Off       | Two-centre overlap integral. |
|-k / --kinetic           |Off       | Two-centre kinetic energy integral. |
|-n / --nuclear           |Off       | Three-centre nuclear attraction integral. |
|-e / --electron          |Off       | Three-centre overlap integral with an s-type third centre, used for the electronic potential. |
|-m / --momentum          |Off       | Two-centre momentum integral. The factor of minus i h-bar is implicit, so this equals the Pulay correction integral. |
|-mv / --massvel          |Off       | Two-centre mass-velocity integral. |
|-d / --dipole            |Off       | Two-centre dipole moment integral. |
|-dk / --dkinetic         |Off       | Two-centre derivative of the kinetic energy. |
|-dncb / --dnuclearcb     |Off       | Three-centre derivative nuclear attraction, C distinct from B: `<A|1/C C|B>`. |
|-dnbb / --dnuclearbb     |Off       | Three-centre derivative nuclear attraction, C equal to B: `<A|1/B B|B>`. |
|-dnbc / --dnuclearbc     |Off       | Three-centre derivative nuclear attraction, C distinct from B: `<A|1/B B|C>`. |
|-decb / --delectroncb    |Off       | Three-centre derivative electronic potential, C distinct from B: `<A|C|B>`. |
|-debb / --delectronbb    |Off       | Three-centre derivative electronic potential, C equal to B: `<A|B|B>`. |
|-debc / --delectronbc    |Off       | Three-centre derivative electronic potential, B distinct from C: `<A|B|C>`. |
|-ko / --Koverlap         |Off       | Two-centre overlap with a plane wave, used for the polarization calculation. |
|-a / --all               |Off       | Activate every integral type above. |
|--help                   |N/A       | Print usage information and exit. |

## Output

`gaussIntgTest.f90` in test mode, or `gaussIntegrals.f90`
(`gaussIntegrals.vec.f90` with `-v`) in production mode, written to
the current directory. Compile the test program with the same
compiler as the engine and run it; it reports the analytic and
numerical values side by side for each integral and angular-momentum
pair.

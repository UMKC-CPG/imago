---
layout: default
title: imago.py
script: imago.py
description: The main driver for running electronic structure calculation jobs. When run in a directory prepared by the makeinput.py script, running this will run an electronic structure calculation as specified by the command line arguments and place the output files in the current working directory.
order: 1
parent: script_catalog
---

# Imago 
> ### Imago driver for running jobs

Orchestrates self-consistent field (SCF) and post-SCF calculations by managing input
files, selecting the correct executable, and collecting
output. Supports checkpointing: completed calculations are
skipped on restart.

## Usage

```
USAGE: uolcao [-scf $basis] [-pscf $basis] [
              [-dos  [$edge]]  | [-scfdos  [$edge]] |
              [-bond [$edge]]  | [-scfbond [$edge]] |
              [-dimo [$edge]]  | [-scfdimo [$edge]] |
              [-mtop [$edge]]  | [-scfmtop [$edge]] |
              [-optc [$edge]]  | [-scfoptc [$edge]] |
              [-pacs [$edge]]  | [-scfpacs [$edge]] |
              [-nlop [$edge]]  | [-scfnlop [$edge]] |
              [-sige [$edge]]  | [-scfsige [$edge]] |
              [-sybd [$edge]]  | [-scfsybd [$edge]] |
              [-force [$edge]] | [-scfforce [$edge]] |
              [-field [$edge]] | [-scffield [$edge]] |
              [-loen] ]
              [-serialxyz]
              [-valgrind]
              -help
```

| Variable  | Accepted Values          | Description |
|:--------: |:---------------:         |:-----------:|
|edge       | gs,1s,2s,2p,3s,3p,3d,... |Electron Orbitals |
|basis      | MB,FB,EB,NO              |Minimal, full, or extended basis respectively|

| Option       | Outputs      | Defaults              |  Description |
|:---------:   |:------------:|:-------:              |:-------------|
| [-scf](#self-consistent-field-options)| N/A| FB     | Sets the basis for the scf portion of calculation | 
| [-pscf](#self-consistent-field-options)| N/A| Operation Dependant  | Sets the basis for the post scf portion of calculation|
| [-serialxyz](#-serialxyz)   | N/A          | Not Set               | Sets the x,y, and z component calculations to be run in serial |
| [-valgrind](#-valgrin)    | N/A          | Not Set               | Runs fortran scripts in a valgrind environment |
| [-dos/-scfdos](#density-of-states) | TDOS localization index plot and PDOS raw file | scf: FB; pscf: FB | Density of States Calculation |
|[-bond/-scfbond](#bond-order-and-q)| bond and Q\* raw files | scf: FB; pscf: MB | Bond order and Q\* calculation |
| [-dimo/-scfdimo](#dipole-moment)| total dipole moment file | scf: MB; pscf: MB | Dipole Moment Calculation |
|[-mtop/-scfmtop](#polarization-properties)| Polarization file | scf: FB; pscf: FB| Polarization Calculation|
|[-optc/-scfoptc](#optical-properties)| distinct files that contain x,y,z decompositions such as the optical conductivity, epsilon1, epsilon2, energy loss function (ELF), refractive index, absorption coefficient, reflectivity, and the imaginary epsilon1. Also included will be a file that contains the total epsilon1, epsilon2, and ELF (without x,y,z decomposition). | scf: FB; pscf: EB | Full Optical Properties calculation |
|[-pacs/-scfpacs](#photo-absorption-cross-section)| a single core-level photo-absorption spectrum file (energy vs. total and x,y,z components) | scf: FB; pscf: EB     |Will do a spectral calculation from an initial state to a final state. The type of calculation, and the potentials used depends on the case of the requested edge.  You MUST provide an edge for the -pacs option. PACS = Photo absorption cross section. | Photo-absorption cross section. Edge argument required |
|[-nlop/scfnlop](#non-linear-optical-properties) | chi1 and chi2 files (real and imaginary parts of the nonlinear susceptibility, each x,y,z decomposed) plus a combined total file | scf: FB; pscf: EB     | non-linear optical properties | 
|[-sige/-scfsige](#sigmae-curve)| a single conductivity curve file, sigma(E) vs. energy (total and x,y,z components) | scf: FB; pscf: FB     | Sigma(E) curve of optical transistion near the Fermi energy |
|[-sybd/-scfsybd](#symmetric-band-structure)| a plottable band-structure file (energies along the k-path) plus a valence-dimension (vdim) raw file | scf: FB; pscf: FB     | Symmetric Band Structure; should be a post-scf job |
|[-force/-scfforce](#inter-atomic-forces)| a force-integral data file (experimental; format preliminary) | scf: FB; pscf: FB     | Experimental; In development |
|[-field/-scffield](#quantum-and-electronic-properties)| 1-D field profiles along a, b, c; a charge-center file; and an XDMF + HDF5 volumetric dataset (wave function, charge density, potential) | scf: FB; pscf: FB     | Complex (or real for Gamma k-point) wave function calculation and electronic property suite |
|[-loen](#local-environment-analysis)| a per-site table of bispectrum components describing each atom's local environment | N/A                   | Local Environment analysis |

## Detailed Calculation Information

The -[option] and -scf[option] convention is consistent. Prepending scf to the calculation option will skip the scf portion of the calculation and instead use a pre-calculated scf for the property calculations. 

Most output files follow this naming convention, where the extension is option-specific:
**[$edge]-[option]-[$basis].[extension]**
All operations will output a log file detailing the calculation process and certain variables with the naming convention:
**[$edge]-pscf+[option]-[$basis].out**

#### Self-Consistent Field Options 

> -scf/-pscf
Question For Rulis: What exactly are the differences between the basis sets, why would someone pick them, and what are they good for?

#### Computation-specific Options

These options don't change the calculations performed, instead they affect HOW imago performs the calculations.

##### -serialxyz

This option will save memory by running the x, y, and z components one after another instead of in parallel. This will increase runtime.
 
##### -valgrind

Runs the program in a valgrind environment. This is for memory usage analysis and debugging. Not necessary for most scientific applications. Use when trying to optimize the program on your machine.

#### Density of States

> -dos/-scfdos

Used to calculate the number of states available to elctrons in the chosen system. If requesting an edge, it will use the potential for that edge, be default the edge is gs. Upon successful completion it will yield two major output types: TDOS localization index plot files, and PDOS raw files.

**TDOS Index Plot Files:**
- Directly plottable. For quick analysis, you can use the the [plotting file]. 
- Follows the [$edge]-dos-[$basis].t.plot naming convention.
- Contains two human-readable columns seperated by a tab. One for energy in Ev, and one for the Total Density of States (TDOS).
- Question for Rulis: what is the .loci.plot file?

**PDOS Raw File**:
- Must be postprocessed by the makePDOS.py script
- Follows the [$edge]-dos-[$basis].p.raw naming convention
- Information about the postprocessing to come

#### Bond Order and Q\*

> -bond/-scfbond

The bond order calculation covers all the atoms in the system and outputs a [$edge]-bond-[$basis].raw file. This must be post processed by the makeBOND.py file. The raw file is human readable containing entries for each atom with each entry having the following information:
- Atom Number (ATOM\_NUM)
   - NOTE: this is not the atomic number, it is an internal label used to identify the atoms in the system 
- System Number (SYSTEM\_NUM)
- Atom name, the abbreviation (ELEMENT\_NAME) 
- Element, species, and type id's (\*\_ID)
- Neutral Valence- (NEUT\_VALE...), atom- (ATOM...), and Orbital (ATOM\_ORBITAL...) charge (...\_CHARGE)
- The number of atoms bonded to the atom (NUM\_BONDED\_ATOMS) followed by a table identifing the bonded atoms
- The number of bond angles (NUM\_BOND\_ANGLES) followed by a table listing the bond angles

Question for Rulis: Where does the Q\* information live and how is it extracted?
Answer: the \*\_CHARGE fields

#### Dipole Moment

> -dimo/-scfdimo

The dipole moment calculation outputs a total dipole moment output file. Follows the default naming convention with the extenstion .t.plot. The file contains the following data:
- Electric Moment [Elec Mom (x,y,z)]
   - has three values for the x,y, and z components
- Electric Polarization [Elec. Polarization (x,y,z) [C/m^2]]
   - Three values for the x, y, and z components
- Nuclear Moment (x,y,z): 
   - Three values for the x, y, and z components
- Dipole Moment (a.u.)
   - Three values for the x, y, and z components
- Dipole Moment (Debya)
   - Three values for the x, y, and z components
- Total Dipole Moment (Debye)
   - Single Value
- Polarization (x,y,z)
   - Three values for the x, y, and z components
- Total Polarization [C/m^2]
   - Single value

#### Polarization Properties

> -mtop/-scfmtop
Output: polarization output file

#### Optical Properties

> -optc/-scfoptc

Computes the linear optical response of the system over an energy
range. The work is done in two stages:

1. **Optical conductivity and epsilon2.** The interband optical
   conductivity and the imaginary dielectric function (epsilon2) are
   evaluated directly from the transitions between occupied and
   unoccupied states. Following Cohen and Chelikowsky (section 4.1,
   eq. 4.10, after Ehrenreich and Cohen, Phys. Rev. 115, 786
   (1959)), each transition is broadened by a Gaussian and summed
   over the Brillouin zone. The energy axis runs from near zero up to
   the largest transition energy and is reported in eV.
2. **Kramers-Kronig conversion.** The Kramers-Kronig relation is
   applied to the epsilon2 spectrum to recover the real dielectric
   function (epsilon1), and from epsilon1 and epsilon2 the remaining
   quantities are derived: the energy loss function, refractive
   index, extinction coefficient, reflectivity, and absorption
   coefficient. The relevant closed forms are
   `ELF = eps2 / (eps1^2 + eps2^2)`,
   `n = sqrt[(sqrt(eps1^2 + eps2^2) + eps1)/2]`, and
   `k = sqrt[(sqrt(eps1^2 + eps2^2) - eps1)/2]`.

Every quantity is resolved into x, y, and z Cartesian components
plus a "total" column (the average of the three, i.e. the isotropic
polycrystalline value). One additional combined file collects the
total epsilon1, epsilon2, and ELF without the x/y/z decomposition.

**Output files.** All files follow the naming pattern
`[\$edge]_optc-[\$basis].[extension]`, where the leading `t` in the
extension marks the total (whole-system) result:

| Extension     | Description                                       |
|:------------- |:------------------------------------------------- |
| `t.cond.plot` | Optical conductivity, sigma(E)                    |
| `t.eps2.plot` | Imaginary dielectric function, epsilon2           |
| `t.eps1.plot` | Real dielectric function, epsilon1                |
| `t.eps1i.plot`| Alternate imaginary-integrand epsilon1            |
| `t.elf.plot`  | Energy loss function, ELF                         |
| `t.nref.plot` | Refractive index, n                               |
| `t.kext.plot` | Extinction coefficient, k (kappa)                 |
| `t.Rref.plot` | Reflectivity, R                                   |
| `t.aabs.plot` | Absorption coefficient, alpha                     |
| `t.plot`      | Combined total epsilon1, epsilon2, ELF (no x/y/z) |

**File format.** Each decomposed `.plot` file is plain text: a
one-line column header followed by rows of five fixed-width,
scientific-notation columns. The five columns are:

```
Energy(eV)   total   x   y   z
```

The header labels the total and component columns per quantity --
for example `Energy totalCond xCond yCond zCond` for conductivity,
`Energy totaln xn yn zn` for the refractive index, `... totalk ...`
for the extinction coefficient, `... totalReflectivity ...` for
reflectivity, and `... totalAbsorpCoeff xAlpha ...` for absorption.

The combined `.t.plot` file instead has just four columns --
energy, epsilon1, epsilon2, and ELF -- giving a quick, plot-ready
summary of the dielectric response without the Cartesian breakdown.

**Spin-polarized runs** additionally emit `.up` and `.dn` variants
of every file. **Partial optical properties** (POPTC) add a parallel
set of `.p.*` files decomposed by atom type, by atom, or by
angular-momentum channel; these are post-processed by the
`processPOPTC.py` script.

#### Non-Linear Optical Properties

> -nlop/-scfnlop

Computes the nonlinear optical susceptibility of the system over an
energy range. It reuses the linear optical-properties machinery, so
the calculation again proceeds in two stages:

1. **Imaginary susceptibility.** The imaginary part of the nonlinear
   susceptibility (chi2) is computed from the transition sum, driven
   by its own energy cutoff, maximum transition energy, energy step,
   and Gaussian broadening width set in the input file.
2. **Kramers-Kronig conversion.** The Kramers-Kronig relation is
   applied to chi2, yielding the real part (chi1) and a combined
   summary file, exactly as for `-optc`.

As with the linear case, each component is resolved into x, y, and z
plus a "total" (three-component average) column.

**Output files.** These follow the same naming pattern as the linear
optical properties, `[\$edge]_optc-[\$basis].[extension]`. Note that
they reuse the `optc` name tag (they are *not* tagged `nlop`), and
the extension distinguishes them:

| Extension    | Description                                          |
|:------------ |:---------------------------------------------------- |
| `chi2.plot`  | Imaginary part of the nonlinear susceptibility, chi2 |
| `chi1.plot`  | Real part, chi1, from Kramers-Kronig conversion      |
| `plot`       | Combined total chi1, chi2, and loss function (no x/y/z) |
| `chi2.50`    | Raw archival copy of the imaginary part (chi2)       |

**File format.** The `chi2.plot` and `chi1.plot` files use the same
plain-text layout as the linear optical files: a one-line header
followed by five fixed-width, scientific-notation columns
(`Energy(eV) total x y z`). The combined `plot` file has four
columns.

A note for readers: the column headers in `chi1.plot` and the
combined `plot` file read `Epsilon1`, `Epsilon2`, and `ELF`,
inherited from the linear optical calculation. For a nonlinear run
these labels stand in for the corresponding chi1 / chi2 quantities;
the data are the nonlinear susceptibility, not the linear dielectric
function.

#### Photo-Absorption Cross Section

> -pacs/-scfpacs

Computes a core-level photo-absorption (PACS) spectrum: the
transition probability from a deep core state of a target atom into
the conduction band, as a function of photon energy. This is the
same class of calculation as XANES/ELNES, and you **must** supply an
edge (e.g. `1s`, `2p`) to select which core state is excited; the
edge also determines the excited-state potential that is used.

The spectrum is broadened by a Gaussian in fixed energy steps. What
makes a PACS run distinctive is its energy onset: instead of
starting near zero, the energy scale begins at the total-energy
difference between the ground state and the core-excited state
(rounded down to a multiple of five), because the transitions
originate from a deep core level rather than the valence band.
Unlike the linear and nonlinear optical properties, PACS performs no
Kramers-Kronig conversion, so it emits a single spectrum file.

**Output files.** Naming pattern `[\$edge]_pacs-[\$basis].[extension]`:

| Extension | Description                                             |
|:--------- |:------------------------------------------------------- |
| `plot`    | Photo-absorption spectrum vs. energy                    |
| `plot.50` | Raw archival copy of the spectrum (Fortran unit suffix) |

**File format.** The `plot` file is plain text: a one-line header
followed by five fixed-width, scientific-notation columns --
`Energy(eV) total x y z` -- where `total` is the average of the
three Cartesian components. (Spin-polarized runs emit `.up` and
`.dn` variants of each file.)

#### Sigma(E) Curve

> -sige/-scfsige

Computes the energy-resolved electrical conductivity, sigma(E),
arising from optical transitions near the Fermi energy. Where the
`-optc` calculation targets the full interband optical response,
`-sige` focuses on the low-energy transitions that govern electronic
transport, and it accumulates a single frequency-dependent
conductivity curve.

The transition sum is broadened and scaled so that the reported
conductivity is expressed in units of inverse micro-ohm-centimeters,
`(micro-ohm-cm)^-1`. The run also reports a single scalar total
conductivity ("The total sigma is: ...") in the calculation log
(`.out`).

**Output files.** Naming pattern `[\$edge]_sige-[\$basis].[extension]`:

| Extension   | Description                                          |
|:----------- |:---------------------------------------------------- |
| `cond.plot` | Conductivity sigma(E) vs. energy                     |
| `cond.50`   | Raw archival copy of the curve (Fortran unit suffix) |

**File format.** The `cond.plot` file is plain text with **no
header line** -- it is pure numeric data in five fixed-width,
scientific-notation columns. The columns are `Energy(eV) total x y
z`, where `total` is the average of the three Cartesian components.
(Spin-polarized runs emit `.up` and `.dn` variants of each file.)

#### Symmetric Band Structure

> -sybd/-scfsybd

Computes the electronic band structure along a high-symmetry path
through the Brillouin zone. The path is chosen by using the
crystal's symmetry to connect the special k-points, which keeps the
calculation small, so this is best suited to highly symmetric
systems such as crystals. The path itself is selected when the job
is prepared (the `-sybdpath` option of `makeinput.py`), and this
should be run as a post-SCF job. The resulting band data are
formatted for plotting by the `makeSYBD.py` script, which `imago.py`
runs automatically.

**Output files.** Naming pattern `[\$edge]_sybd-[\$basis].[extension]`,
with one companion file that carries the `vdim` (valence dimension)
tag:

| File                     | Description |
|:------------------------ |:----------- |
| `[\$edge]_sybd-[\$basis].plot` | Plottable band structure: band energies along the k-path |
| `vdim-[\$basis].raw`       | Listing of the valence basis orbitals (see below) |

Two additional raw copies, `[\$edge]_sybd-[\$basis].raw.31` and
`.raw.33`, are retained in the working directory. (Spin-polarized
runs emit `.up` and `.dn` band-structure plots.)

**Band-structure `.plot` file.** Plain text with one row per k-point
sampled along the symmetry path. The first column is the cumulative
distance along the path (the horizontal axis of a band-structure
plot); every remaining column is one band's energy, in eV, at that
k-point. To plot the band structure, draw each energy column against
the first column. For the kyanite example run, this file holds 300
path points and 318 band columns.

**Valence-dimension `vdim` file.** Plain text with one row per
valence basis orbital, enumerating the basis that sets the number of
bands. Each row gives a running state index, the atom number, the
element name, the atomic number, the element / species / type
identifiers, the orbital quantum numbers (n, l, m), and the atom
position. This is useful for relating individual bands back to the
orbitals and atoms they come from.

#### Quantum and Electronic Properties

> -field/-scffield

Evaluates the system's quantum and electronic fields on a real-space
grid, producing a suite of spatially resolved quantities. For each
of the following the calculation samples the value throughout the
cell:

- the wave function, split into its real and imaginary parts (the
  imaginary part is zero for a Gamma-point calculation, where the
  wave function is real);
- the electron charge density;
- the potential.

Each of these is reported in several forms: the actual
self-consistent system ("live"), a spin form (spin-up minus
spin-down, meaningful for spin-polarized runs), the difference from
a neutral-atom reference (which highlights how bonding rearranges
charge), and the neutral-atom reference itself.

**Output files.** Naming pattern `[\$edge]_field-[\$basis].[extension]`:

| Extension                          | Description |
|:---------------------------------- |:----------- |
| `prof-a.dat`, `prof-b.dat`, `prof-c.dat` | One-dimensional line profiles of the fields along the a, b, and c cell axes |
| `rho.dat`                          | Charge-center data for the cell |
| `xdmf3`                            | An XDMF descriptor for the full three-dimensional volumetric data (see below) |

**Profile `.dat` files.** Plain text with a one-line column header
followed by one row per sampled point along the chosen axis. The
columns are the position along the axis followed by the wave
function real part, wave function imaginary part, wave function
value, charge density, and potential (all for the live spin-up +
spin-down form). The three files differ only in which cell axis (a,
b, or c) is traversed.

**Charge-center `rho.dat` file.** A small plain-text block of
charge-center values for the cell.

**Volumetric `xdmf3` file.** An XDMF (XML) file that, together with a
companion HDF5 file of the same base name, describes the full
three-dimensional field data for a visualization tool such as
ParaView or VisIt. It defines a uniform structured grid (a 10x10x10
mesh in the kyanite example) and one node-centered scalar attribute
for every field and form listed above -- the wave function (real and
imaginary parts), charge density, and potential, each in its live,
spin, difference, and neutral variants. The bulk of the data lives
in the paired HDF5 file; the `xdmf3` file is the small text index
that points into it.

#### Local Environment Analysis

> -loen

Quantifies the local environment around each atom -- the geometric
arrangement of its neighbors -- and expresses it as a set of
rotationally invariant *bispectrum components*. Neighbors are
gathered out to a radial cutoff and smoothly down-weighted near the
cutoff edge, and the bispectrum then encodes the shape of that
neighbor density in a way that does not depend on how the cell is
oriented. The result is a compact numerical fingerprint of each
site's surroundings, useful for comparing and classifying atomic
environments. Unlike the other calculations there is no `-scfloen`
form; `-loen` runs as its own post-SCF-style pass.

**Output files.** Naming pattern `[\$edge]_loen-[\$basis].[extension]`:

| Extension | Description                                              |
|:--------- |:-------------------------------------------------------- |
| `plot`    | Per-site table of bispectrum components (see below)      |
| `plot.21` | Raw archival copy of the table (retained in working dir) |

**File format.** The `plot` file is plain text with a one-line
header, then one row per site. Each row begins with identity columns
-- `site#`, `element`, `species` (the per-element species), `type_sp`
(the per-element-species type), and `type_flat` (the flat global
type index) -- followed by one column per bispectrum coupling
channel (headed `2j_NNN`), and finally a `total` column that sums the
components for that site. The number of component columns depends on
the coupling parameters chosen for the analysis.

#### Inter-Atomic Forces

> -force/-scfforce

Computes the forces acting between the atoms in the system, obtained
from the derivatives of the Hamiltonian and overlap integrals with
respect to atomic position.

> **Experimental / in development.** This calculation is not yet
> complete, and both its results and its output format should be
> treated as preliminary and subject to change. The current output
> is a diagnostic dump of the underlying force integrals rather than
> a finished per-atom force table.

**Output files.** Naming pattern `[\$edge]_force-[\$basis].[extension]`:

| Extension | Description                                              |
|:--------- |:-------------------------------------------------------- |
| `dat`     | Force-integral data for the system (see the note above)  |
| `dat.98`  | Raw archival copy of the data (retained in working dir)  |

Spin-polarized runs emit `.up.dat` and `.dn.dat` for the two spin
channels.

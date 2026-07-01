---
script: imago.py
description: The main driver for running electronic structure calculation jobs. When run in a directory prepared by the makeinput.py script, running this will run an electronic structure calculation as specified by the command line arguments and place the output files in the current working directory.
order: 1
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
|[-optc/-scfoptc](#optical-properties)| distinct files that contain x,y,z decompositions such as the optical conductivity, epsilon1, epsilon2, energy loss function (ELF), refractive index, absorption coefficient, reflectivity, and the imaginary epsilon1. Also included will be a file that contains the total epsilon1, epsilon2, and ELF (without x,y,z decomposition). | scf: FB; pscf: FB | Full Optical Properties calculation |
|[-pacs/-scfpacs](#photo-absorption-cross-section)| tbd          | scf: FB; pscf: EB     |Will do a spectral calculation from an initial state to a final state. The type of calculation, and the potentials used depends on the case of the requested edge.  You MUST provide an edge for the -pacs option. PACS = Photo absorption cross section. | Photo-absorption cross section. Edge argument required |
|[-nlop/scfnlop](#non-linear-optical-properties) | tbd          | scf: FB; pscf: EB     | non-linear optical properties | 
|[-sige/-scfsige](#sigmae-curve)| tbd          | scf: FB; pscf: FB     | Sigma(E) curve of optical transistion near the Fermi energy |
|[-sybd/-scfsybd](#symmetric-band-structure)| tbd          | scf: FB; pscf: FB     | Symmetric Band Structure; should be a post-scf job |
|[-force/-scfforce](#inter-atomic-forces)| tbd        | scf: FB; pscf: FB     | Experimental; In development |
|[-field/-scffield](#quantum-and-electronix-properties)| field data | scf: FB; pscf: FB     | Complex (or real for Gamma k-point) wave function calculation and electronic property suite |
|[-loen](#local-environment-analysis)| tbd          | N/A                   | Local Environment analysis |

## Detailed Calculation Information

The -[option] and -scf[option] convention is consistent. Prepending scf to the calculation option will skip the scf portion of the calculation and instead use a pre-calculated scf for the property calculations. 

Most output files follow this naming convention, where the extension is option-specific:
**[$edge]-[option]-[$basis].[extension]**
All operations will output a log file detailing the calculation process and certain variables with the naming convention:
**[$edge]-pscf+[option]-[$basis].out**

#### Self-Consistent Field Options 

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

Output: polarization output file

#### Optical Properties

#### Non-Linear Optical Properties

#### Photo-Absorption Cross Section

#### Sigma(E) Curve

#### Symmetric Band Structure

Simplifies the Band Structure calculation by using symmetry to simplify the math. Best used for highly symmetric systems like crystals. Outputs a 
#### Quantum and Electronic Properties

#### Local Environment Analysis

#### Inter-Atomic Forces


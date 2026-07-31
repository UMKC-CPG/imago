---
layout: default
title: makeinput.py
script: makeinput.py
description: 'Takes a skeleton file input (default: imago.skl) and prepares all the files imago.py requires for input. Including a default slurm file for submitting to a system utilizing slurm submission.'
parent: Script Catalog
---

# MakeInput

Takes a skeleton file input (default: imago.skl) and prepares all the files imago.py requires for input. Including a default slurm file for submitting to a system utilizing slurm submission.

**Suggested usage:** within your imago directory, make a jobs/ folder. Organize subdirectories as you will, but each skeleton file should have its own subdirectory. Example: imago/jobs/c/diamond/ for a diamond structure. Place the corresponding skeleton file within the directory, and run makeinput.py from that directory.

## Usage

**Commandline:**

```
makeinput [-basisdb \$atomicBDB] [-potdb \$atomicPDB]
          [-modpot \$modElementName \$minModTerm \$maxModTerm
                   \$numModTerms]
          [-subbasis \$basisSubOut \$basisSubIn [-subbasis ...]]
          [-subpot \$potSubOut \$potSubIn [-subpot ...]]
          [[[-scfkp \$a \$b \$c] [-pscfkp \$a \$b \$c]
           | [-kp \$a \$b \$c]]
          [-kpshift \$a \$b \$c]
          [-printbz \$bz \$scaleFactor]
          [-xccode \$xcCode]
          [-xcmesh [-numvect \$numSampVectors]
                   [-weight \$xcInWeight \$xcOutWeight]
                   [-samp \$xcInSamp \$xcOutSamp \$xcSpacingSamp]]
          [-target <[-atom \$targetAtom] || [-atxyz \$x \$y \$z] ||
                    [-atabc \$a \$b \$c]> [-sphere \$targetRadius]
                    [-zone \$targetZone] [-operand \$targetOp]
                    [-relate \$targetRelation]]
          [-block <-abc \$froma \$toa \$fromb \$tob \$fromz \$toc>
                    [-zone \$blockZone] [-operand \$blockOp]
                    [-relate \$blockRelation]]
          [-reduce [-level \$reduceLevel] [-thick \$reduceThick]
                   [-cutoff \$reduceCutoff] [-operand \$reduceOp]
                   [-tolerance \$reduceTolerance]
                   [-selection \$reduceSelection]]
          [-xanes [-sphere \$xanesRadius]
                  [-atom \$xanesAtom1 [\$xanesAtom2 [...]]]]
          [-sybdpath \$sybdPath]
          [-rel]
          [-statefactor \$factor]
          [-pdb] [-cif] [-basisVis]
          [-emu]
          [-nocore]
          [-slurm [-p \$partition] [-a \$account] [-t \$time]
                  [-m \$memory] [-n \$cpus] [-N \$nodes]]
          [-help]
```

**Variable Table:**

| Variable        | Accepted Values                | Description |
|:---------:      |:---------------------:         |:-----------:|
|atomicBDB        |String; Directory Path          | Path to a different atomic basis function database |
|atmoicPDB        |String; Directory Path          | Path to a different potential function database |
|modElementName   |Element; lowercase              | What element in the system to adjust the potential for |
|minModTerm       |float                           | Minimum value of the adjusted potential |
|maxModTerm       |float                           | Maximum value of the adjusted potential |
|numModTerms      |Integer                         | Number of terms in the adjusted potential |
|basisSubOut      |Element or species; lowercase   | What element or species to make a basis definition substitution | 
|basisSubIn       |integer                         | What basis definition to substitute in |
|potSubOut        |Element or species; lowercase   | What element or species to make a potential definition substitution |
|potSubIn         |integer                         | What potential definition to substitute in |
|a                |integer (float)                 | The size of on side of the kpoint mesh (integer) or one component of the k-point mesh shift (float) |
|b                |integer (float)                 | The size of on side of the kpoint mesh (integer) or one component of the k-point mesh shift (float) |
|c                |integer (float)                 | The size of on side of the kpoint mesh (integer) or one component of the k-point mesh shift (float) |
|bz               |integer                         | Which Brillouin Zone to print |
|scaleFactor      |float                           | Scale factor of the Brillouin Zone scaling |
|xcCode           |integer matching an xc code     | Exchange Correlation Code. Accepted integers and their corresponding exchange codes listed in xc\_code.dat |
|numSampVectors   | integer                        | the number of radial sampling vectors |
|xcInWeight       |float                           | weight of the points within the radial cutoff|
|xcOutWeight      |float                           | weight of the points outside the radial cutoff|
|xcInSamp         |float                           | Used to compute the minimum radius a sampling point should be |
|xcOutSamp        |float                           | Used to compute the maximum radius a sampling point can be |
|xcSpacingSamp    |float                           | Computes the spacing between points |
|targetAtom       | atomic number                  | Which atom to make the center of the target sphere|
|x                |float                           | x position of the center of the target sphere |
|y                |float                           | y position of the center of the target sphere |
|z                |float                           | z position of the center of the target sphere |
|targetRadius     |float                           | Radius of the target sphere |
|targetZone       | in or out                      | Whether to consider the atoms inside or outside the zone sphere |
|targetOp         |Species, type, reduce           | What category to group the atoms (species or type   
|targetRelation   | relation                       | Groups atoms by their relation| *Question:* What are the relations and how are they used?
|froma            | float                          | Start of the target zone along the a axis|
|toa              | float or a                     | End of the target zone along the a axis 'a' selects the entire cell axis|
|fromb            | float                          | Start of the target zone along the b axis|
|tob              | float or b                     | End of the target zone along the b axis. 'b' selects the entire cell axis|
|fromz            | float                          | Start of the target zone along the c axis|
|toc              | float or c                     | End of the target zone along the c axis. 'c' selects the entire cell axis|
|blockZone        | in or out                      | Whether to consider the atoms inside or outside the zone block
|blockOp          | Species, type, reduce          | How to consider the atoms in the target block |
|blockRelation    | alike or diff                  | What relation to group by |
|reduceLevel      | Integer                        | the number of shells to use when reducing |
|reduceThick      | float                          | thickness of the shells |
|reduceCutoff     | float                          | maximum cutoff radius of shells |
|reduceOp         | Species, type                  | Operand for the reduce command. How to label the groups |
|reduceTolerance  | float                          | distance threshold for comparing the same shell number from different atoms|
|reduceSelection  |                                | Experimental | 
|xanesRadius      | float                          | Radius of the sphere around the xanes atom |
|xanesAtom[#]     |list(atomic numbers) or N..M ranges | what atoms to target |
|sybdPath         | path matching: '\$IMAGO\_DATA/sybdDB/\*'| path to the type of cell for symmetric band structure calculation|
|factor           | integer                        | The number of states to calculate. Must be a multiple of the number of valence elecrons |
|partition        |string                          |Partition name requested in the slurm file|
|account          |string                          |Account name; what account to charge|
|time             |HH:MM:SS                        |How much time to request in the slurm file
|memory           |slurm memory                    | [number]G format specifying that amount of memory in the slurm file|
|cpus             |integer                         | Number of cpus requested in the slurm file|
|nodes            |integer                         | Number of nodes requested in the slurm file|

**Option Table:**

| Option    | Outputs      | Defaults              | Descripition |
|:---------:|:------------:|:---------------------:|:-------------|
|-basisdb   |N/A           |$IMAGO\_DATA/atomicBDB | Give an alternate location for basis function database|
|-potdb     |N/A           |$IMAGO\_DATA/atomicPDB | Give an alternate location for potential function database| 
|-subbasis  |N/A           |No Substitutions       | Subsitute a basis set for an element or species; repeatable for different elements/species|
|-subpot    |N/A           |No Substitutions       | Substitute a potential for an element or species; repeatable|
|-pot       |N/A           |Default tagged entry   | Override for the augmented per-element potential.|
|-modpot    |N/A           |No Modifications       | Modify potential to have explicitly set minimum and maximum values; not repeatable|
|-scfkp     |N/A           |Gamma kpoint/1 1 1     | Defines the kpoint mesh for scf calculations|
|-pscfkp    |N/A           |Gamma kpoint/1 1 1     | Defines the kpoint mesh for post-scf calculations|
|-kp        |N/A           |Gamma kpoint/1 1 1     | Apply the same kpoint selection to both SCF and post-SCF|
|-scfkpd    |N/A           |Not used               | Specify the k-point volume density for SCF calculations|
|-pscfkpd   |N/A           |Not used               | Specify the k-point volume density for post-SCF calculations|
|-kpd       |N/A           |Not used               | Use the same k-point volumen density for both SCF and post-SCF calculations|
|-scfkpint  |KPOINT\_INTG\_CODE in input/kp-scf.dat|0                      | Select, by integer ID, the k-point integration method for scf calculation|
|-pscfkpint |KPOINT\_INTG\_CODE in input/kp-pscf.dat|0|Select, by integer ID, the k-point integration method for post-scf calculations|
|-kpint     |N/A           | 0                     | Use the same integration method for scf and pscf calculations|
|-kpshift   |KP\_SHIFT\_A\_B\_C in input/kp-scf.dat|Lattice type dependant |Force a specific shift to the kpoint mesh by fractional amounts a, b, and c|
|-printbz   |Additional file describing the Brilloiun zone|First zone, no scaling |Selects and prints a brillouin zone based on the first integer and scaled by the second.|
|-converg   |N/A           |Sourced from makeinputrc (conver\_main)|Set the scf convergence limit for the run|
|-xccode    |N/A           | 100 (Wigner interpolation)   | set the exchange correlation code to be used.|
|-xcmesh    |N/A           | see dedicated section | Define the real space mesh for exchange correlation code. See dedicated section for suboptions |
|-target    |N/A           | Off. See dedicated section for supoption defaults.   | Define a central point and spherical zone for a reduce call. See dedicated section for suboptions |
|-block     |N/A           | Off. See Dedicated section for suboption defaults| Define a block or slab zone for a reduce call. See dedicated section for suboptions  |
|-reduce    |N/A           | Off. See section for suboption defaults| Collect information about each atom and group all atoms based on similarity. |
|-xanes     |N/A           | Off. See section for suboption defaults | Make one input file set for each atomic species or listed xanes atoms. See dedicated section for suboptions|
|-sybdpath  |N/A           |Not used               | Specify cell type and path to be used in sybd calculation. |
|-rel       |N/A           |Non-relativistic       | Prepare Imago input files for relativistic calculation. **NOT YET FUNCTIONAL**|
|-statefactor| N/A         |2.5                    | Select the number of states to calculate. |
|-pdb       |N/A           |Not generated          | Generate a Protein Data Bank crystal structure file|
|-cif       |N/A           |Not generated          | Generate a Crystallographic Interchange Format (CIF) structure file|
|-basisVis  | Atomic Orbital Basis visualization data files | Not generated| Produce radial portion data files of the basis functions and a set of POVRay scene files   |
|-emu       |EMU Configuration files   | Not generated| Genrate EMU configuration files  |
|-slurm     | custom slurm script   | On. See section for suboption defaults| Generate a slurm file with the given parameters. See dedicated section for suboptions   |

## Option Details

### Database Substitutions

As an electronic structure calculation, Imago needs various information about the potential and basis functions for every element. This information is stored in a database (a directory). While Imago comes with default values, you may wish to use your own, non-default potential and basis functions. If so, when running makeinput.py, you can define database substitution paths.

When doing this, Imago expects certain conventions to work properly. Ensure that your database follow these for it to work properly. Basis and Potential databases require different conventions so see each section on its own.

#### **Basis Databases**

> -basisdb

Within the database directory each element has its own subdirectory named by the atomic symbol in all lowercase (e.g. al for aluminum). Inside you need

#### **Potential Databases**

> -potdb

Within the database directory each element has its own subdirectory named by the atomic symbol in all lowercase (e.g. al for aluminum). Inide each subdirectory you need the following plaintext files:

* coeff1
   - The first line is the number of coeffecients in the file
   - Each coeffecient has its own line and five columns, each with a value in scientific notation.
   - Each column is seperated by a single space
* coeff1.isolated
   - This has the same format as coeff1 except:
   - The coeffecients are those which apply to the atom in vacuum.
* pot1
   - A series of lines with the first one being a variable and the next its value. These are, in order:
   - NUCLEAR\_CHARGE\_\_ALPHA
    + The entry is two columns of numbers separated by a single space
   - COVALENT\_RADIUS
    + a single number
   - NUM\_ALPHAS
    + The number of coeffecients. This should be the same as the first link if the coeff\* files.
   - ALPHAS
    + Two columns of numbers in scientific notation seperated by a space.

### Per-element Substitutions

While you can substitute out entire databases, you can also make substitutions for basis sets and potential functions on a per-element basis. The basis and potential substitutions have different requirements so see the details for each.

Additionally, if species are defined within the structure you can make the substitutions on a species-by-species definition rather than just by element.

#### **Basis Substitutions**

> -subbasis

Within each element and species entry, the basis sets are defined in contract\*.dat files where \* is an integer value. You can create additional basis by creating one of these files with a different number. Then, when you run -subbasis you provide that number after the species/element target. Then all atoms matching the selected target or within the selected species will use the contract\*.dat file where \* matches the given integer.

Example:

al contains: contract1.dat and contract2.dat 
when running makeinput, if no -subbasis is given then Imago will use contract1.dat. However if `-subbasis al 2` is provided, then all aluminum atoms will use contract2.dat instead of contract1.dat but everything else will still use contract1.dat.

#### **Potential Substitutions**

> -subpot

Within each element and species entry, the potential functions sets are defined in pot\* files where \* is an integer value. You can create additional potentials by creating one of these files with a different number. Then, when you run -subpot you provide that number after the species/element target. Then all atoms matching the selected target or within the selected species will use the pot\* file where \* matches the given integer.

Example:

al contains: pot1 and pot2 
when running makeinput, if no -subpot is given then Imago will use pot1. However if `-subpot al 2` is provided, then all aluminum atoms will use pot2 instead of pot1 but everything else will still use pot1.

### Modifying Potentials

In additions to full substitutions, you can modify potentials using coeff1.LABEL files within the database. Each element has a default LABEL which is used in most cases. Furthermore, you can adjust potentials by defining a minimum and maximum value to use.

#### **-pot**

Manual override for the augmented per-element potential database (DESIGN 5.6).  When an element carries an s\_gaussian\_pot.toml database, LABEL names which entry to apply uniformly across the structure (e.g. "default\_solid" or "isolated").  Without -pot, each database's default-tagged entry is used.  LABEL must exist in every augmented database it applies to; a missing label is a hard error rather than a silent fall-back, because a manual override expresses deliberate intent.  Elements that have no augmented database fall back to the legacy pot1/coeff1 files (with a warning when -pot was requested).

#### **-modpot**

Modify the potential for the element specified by ELEM so that the minimum value is MIN, the maximum value is MAX, and the number of terms is NTERMS.  Note that this option can only be applied to one element as the program is presently written.

### K-point Options

K-points are used for sampling and integration. Thus, they can have a large effect on both runtime and output. There are several ways to define the k-point mesh, but all have three variations:
- -scf\* which will make the definitions for scf calculations
- -pscf\* which will make the definitions for post-scf calculations
- -\* which will make the same definitions for both types of calculations

Additionally, you can selected from a few predefined integration methods using the kpint options. They follow the same scf\*,pscf\*,\* convention as defining the mesh.

#### **K-point mesh by dimensions**

> -scfkp,-pscfkp,-kp

Specify the mesh of k-points to be used for the SCF calculation. The three parameters are the size of the mesh in the a, b, and c directions.  (e.g. 2 2 2 would produce 8 kpoints before symmetry reductions, and 2 3 4 would produce 24 kpoints before symmetry reductions.)  This applies to the primary SCF programs (setup and main) and also to any subroutines run within the SCF stage (e.g. dos or bond).  Please note that most of the time things like dos or bond are automatically run in the post-SCF stage where a larger number of kpoints are typically used.  One further important note:  If the kpoint scheme given is 1 1 1 then it is assumed that the one point is at the gamma site.  If no kpoint definition is given then 1 general kpoint is used.  This is important since the program will run differently (faster) if the gamma kpoint is used since all the integral matrices will be real with no imaginary component.

**NOTE:** for the sybd post-scf calculation, the used k-point mesh is instead given by the path associated with the crystalline cell or -sybdpath parameter if given.

#### **K-point mesh by density**

> -scfkpd,-pscfkpd,-kpd

Specify the k-point volume density for the SCF calculation as a single number (kpoints per unit reciprocal-space volume, in Bohr^-3).  The total kpoint count will be density * V\_BZ, distributed as uniformly as possible across the three axes. Instead of giving explicit mesh dimensions, the program writes a style-code-2 k-point file and lets Imago compute the per-axis mesh counts from the reciprocal cell geometry at runtime.

#### **Shifting the K-point mesh**

> -kpshift

This option does not differentiate between pscf and scf calculations.

Force a specific shift to the kpoint mesh by the given fractional a, b, c amounts.

#### **K-point Integration Method**

Select the k-point integration method for the SCF calculation.  M is an integer: 0 = Gaussian/histogram (default), 1 = LAT (linear analytic tetrahedron). Higher integers are reserved for future methods. This value is written into the KPOINT\_INTG\_CODE field of the k-point input file.  Works with both mesh mode (-kp) and density mode (-kpd).

**NOTE:** Per-group flages (-scfkpint,-pscfkpint) will override this

### Targeting and Reduction

Making when you make species, you allow yourself to compute the individual contributions of different groups of atoms by species. Targeting allows you to select atoms based on regions and the reduce call makes the species.

#### **Making species out of elements**

> -reduce

**Suboptions Table**

|Option     | Default      | Description|
|:---------:|:------------:|:-----------|
|-level     | 2            |The number of shells |
|-thick     | 0.1 A        |Thickness of each shell |
|-cutoff    | 4.0 A        |Maximum cutoff range |
|-operand   | species      |Currently only with species |
|-tolerance | 0.05         |Distance threshold for comparing shells with the same number from different atoms|
|-selection | 0 (all atoms)|Used to apply to different groups of atoms; usually used with target or block options|

Collect information about each atom in order to group all the atoms based on the similarity of the information.  A group is defined by a series of spherical shells (levels), what neighbor atoms are in the shells, and what the distance is to each shell from the center.  Parameters can be adjusted to determine the number of shells (-level), the thickness of the shells (-thick), a maximum cutoff radius (-cutoff), and a distance threshold for comparing the same shell number from different atoms (-tolerance). The basic idea is to find the nearest atom to the current atom, define a shell of a given thickness and record what atoms are in that shell.  Then repeat (reduceLevel) times with the next nearest atom outside the shell.  Once all the shells have been defined for each atom, the results are compared in terms of the atoms in each shell and the shell distances.  This grouping method is not at all defined by boundaries and so does not have in and out -zone parameters.  It will also not have the ability to make grouping dissimilar (e.g. it will not find a group of atoms that are similar and then make them all different).  This method will only work to make species out of elements, and will not make types out of species.  It can be applied differently to different groups of atoms via -selection.

#### **Spherical Targeting**

> -target

**Suboptions Table**

|Suboption  |Default    |Description|
|:---------:|:---------:|:----------|
|-atom      |undefined  |Define an atom to make the center of the spherical target|
|-atxyz     |undefined  |Define the center of the spherical target by xyz coordinates|
|-atabc     |undefined  |Define the center of the spherical target by abc coordinates|
|-sphere    |3.50 A     |Radius of the sphere|
|-zone      |in         |Whether to consider atoms inside the zone or outside the zone|
|-operand   |species    |What information to consider about the atoms|
|-relate    |diff       |Group atoms by whether they are alike or different|

Consider a point given either by an x,y,z location (-atxyz), an a,b,c location (-atabc), or an atom number (-atom).  Then the atoms either in or out of the zone (sphere of given -sphere radius) will be considered in terms of their -operand (species, type, reduce) and will be grouped by their -relate relation (alike or diff).  If the reduce operand is given then this target is used as a selecting tool for a reduce call.

#### **Block or Slab targeting**

> -block

**Suboptions Table**

|Suboption  |Default       |Description|
|:-------:  |:-----------: |:----------|
|-abc       |undefined     | The block to consider; defined by from coord-\>to coord paradigm|
|-zone      |in            | Whether to consider atoms inside or outside of the zone|
|-operand   |species       | What information to consider about the atoms (element,species,type)|
|-relate    |diff          | Group atoms by whether they are alike or different|

Consider a zone defined like a slab or block of the system using -abc to specify from/to values for each lattice direction.  Note that it is possible to use the letters 'a', 'b', or 'c' in the place of actual numbers when one wants the block "To" value to be the maximum.  (e.g. "-abc 0 a 0 b 0 c" would include the entire cell.)

### Exchange Correlation

Imago uses exchange correlation methods to account for relativistic effects that classical potentials ignore. This is used to increase accuracy of calculations. Several methods are included with Imago, and there is no easy way to define extra methods at the moment.

#### **Selecting Exchange Correlation Code**

> -xccode

Several exchange correlations methods exist, and several are provided with Imago. These are identified internally by an ID number which you use to select the code. The ID numbers (XC\_Code) and their corresponding methods (Functional) are given below.

|Functional             |XC\_Code      |Spin  |Rel | GGA|
|:----------------------|:-----------:|:----:|:--:|:--:|
|Wigner(LDA)            |  100        | 1    | 0  |  0 |
|Ceperley-Alder(LDA)    |  101        | 1    | 0  |  0 |
|Hedin-Lundqvist(LDA)   |  102        | 1    | 0  |  0 |
|Ceperley-Alder(LSDA)   |  150        | 2    | 0  |  0 |
|von\_Barth-Hedin(LSDA)  |  151        | 2    | 0  |  0 |
|unknown\_old\_vB\_H(LSDA) |  152        | 2    | 0  |  0 |
|PBE(GGA)               |  200        | 1    | 0  |  1 | 
|Wigner(rel)            |  300        | 1    | 1  |  0 |
|Ceperley-Alder(rel)    |  301        | 1    | 1  |  0 |
|Hedin-Lundqvist(rel)   |  302        | 1    | 1  |  0 |
|Ceperley-Alder(rel)    |  350        | 2    | 1  |  0 |
|von\_Barth-Hedin(rel)   |  351        | 2    | 1  |  0 |
|Old\_vBarth-Hedin(rel)  |  352        | 2    | 1  |  0 |
|PBE(rel)               |  400        | 1    | 1  |  1 |

#### **Exchange Correlation Mesh**

> -xcmesh

**Suboption Table**

|Suboption     | Default   | Description |
|:------------:|:---------:|:------------|
|-numvect      | 100       | The number of radial vectors used for sampling|
|-weight       | 0.5/0.5   | 
|-samp         | 0.1 3.5 0.8| the distribution of points along the vector directions: in the sphere, out of the sphere, and their spacing respectively |

Exchange correlation methods have to sample the charge density, and the mesh used to do so is speherical and centered on atoms. The suboptions define how this mesh is constructed. Each atom has its own mesh.

### Extra File Creation

There are several options which can be turned on that will have makeinput create extra files. These are most often used when you plan to use some other program to assisst in your postprocessing. While some are used to generate files needed for certain calculations.

#### **Xanes**

> -xanes

**Suboption Table**

|Suboption  | Default      | Description|
|:---------:|:------------:|:-----------|
|-sphere    |3.5 A         | Radius around the target atom(s)|
|-atom      |1 from each species| What atoms to include by atom number |

Instead of generating one set of input files, make one set for each of the atomic species or listed xanes atoms in the system. Each input file set will be different in that a xanes atom of the current input file set will have the core orbitals included in the calculation, and the atoms within a given radius (default 3.50 A) of the xanes atom will all have different types.
 
#### **Brillouin zone**

> -printbz

Create an additional file as part of the input containing a description of the Brillouin zone suitable for importing in a Python script.  The two required options are an integer identifying which Brillouin zone to print and a scale factor that causes the size of the Brillouin zone to be scaled.

#### **Crystallographic Interchange Format**

> -cif

Will generate an atomic structure file using the Crystallographic Interchange Format. Allows for a greater number of species for each atom. Preferred when doing so.

#### **Protein Data Bank**

> -pdb

Generates a Protein Data Bank (PDB) crystal structure file.

#### **Basis Visualization**

> - basisVis

Ask the contract program to produce files in the .inputTemp directory for visualizing the complete atomic orbital basis.  The files are numerical data of the radial part of the basis functions, and a set of POVRay scene files (one for each type of atomic orbital in the system).

#### **EMU Config**

> -emu

Generates EMU configuration files.

#### **SLURM**

**Suboption Table**

|Suboption  | Default      | Description|
|:---------:|:------------:|:-----------|
|-p         |rulisp-lab    | Partition to charge|
|-a         |rulisp-lab    | Account to charge|
|-t         |00:60:00      | Wall time (HH:MM:SS)|
|-m         |10G           | Memory|
|-n         |1             | Number of cpus|
|-N         |1             | Number of nodes|

This is the only extra file option on by default. This will generate a SLURM file that you can use to submit jobs easily if your machine uses that paradigm.

### Convergence Limit

> -converg

The self-consisitent field (scf) calculation can take a very long time depending on the system. If you wish to, you can define a limit after which it will **Just use what it has**?

### SYBD Path

> -sybdpath

Specify the type of cell and particular path to be used in the symmetric band calculation (SYBD).  You should be consistent with the actual cell type or you will get a warning (even though it will still allow you to create the input files).  The valid options are present in the \$IMAGO\_DATA/sybdDB directory.  Simply specify one of those names as the sybdPath.

### States to Calculate

> -statefactor

Select the number of states to calculate as a multiple of the number of valence electrons. By default, 2.5.

### Including Core Electrons

> -nocore

The main advantage of using an OLCAO method, which Imago does, instead of an LCAO method is that we orthoganalize the core electrons out of the equation to increase calculation speed.

Setting this option will include core orbitals in the valence (no orthogonalization). This causes the basis functions that are normally orthogonalized out of the problem to instead be included as part of the list of valence basis functions so that they are also part of the SCF process.  This allows one to compute the core energy eigenvalues. If this option is not given, then basis functions that are designated as "core" will be orthogonalized out of the eigenvalue problem.

### Relativistic Calculations

> -rel

While it is not yet operational, setting -rel will prepare the input files for relativistic calculations.



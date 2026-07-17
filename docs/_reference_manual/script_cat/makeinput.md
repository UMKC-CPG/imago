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
| Option    | Outputs      | Defaults     | Descripition |
| 

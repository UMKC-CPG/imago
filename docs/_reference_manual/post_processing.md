---
title: Post Processing
---

# Post Processing

Each run of Imago will generate several output files for each type of calculation. However, none of these files are plots in and of themselves. Thus, Imago contains several post-processing scripts that take the data and turn them into usable figures, or other file formats used by external tools.

**Post Processing Scripts**

|Script              | Accepted inputs    | Generate Outputs   | Description |
|:------------------:|:------------------:|:------------------:|:------------|
|makeBOND.py         |Bond order raw file | Files useable by external programs | Takes bond order raw file and creates many files useable in openDX, Paraview/VisIT, Origin, or other similar tools|
|makeFittedRhoV.py   |                    |                    | Generates files to visualize various charge density and potential functions on a 3D real-space mesh.|
|makePDOS.py         |
|plot\_deadmd.py     |
|plotgraph.py        |
|pot2plot.py         |Imago scfV and structure.dat| Files containing plotable data columns| |
|processPOPTC.py     |Raw Epsilon 2 pOptc data| Control file for makePDOS| An intermediary step to make using makePDOS.py easier|

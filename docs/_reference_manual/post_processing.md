---
title: Post Processing
layout: default
---

# Post Processing

> ***This page is a stub and incomplete, it is a work in progress***

Each run of Imago will generate several output files for each type of calculation. However, none of these files are plots in and of themselves. Thus, Imago contains several post-processing scripts that take the data and turn them into usable figures, or other file formats used by external tools.

***Note: This table is incomplete but contains the most commonly used scripts***

**Post Processing Scripts**

|Script              | Accepted inputs    | Generate Outputs   | Description |
|:------------------:|:------------------:|:------------------:|:------------|
|makeBOND.py         |Bond order raw file | Files useable by external programs | Takes bond order raw file and creates many files useable in openDX, Paraview/VisIT, Origin, or other similar tools|
|makeFittedRhoV.py   |SCF potential and structure files; also needs potential database | input file for imagoRhoV fortran program and openDX plottable data file | Generates files to visualize various charge density and potential functions on a 3D real-space mesh.|
|makePDOS.py         |Control and raw PDOS data files|PDOS.plot file|Uses the commandline and control file to collect data from raw PDOS files to produce a .plot file useable by plotgraph.py
|plot\_deadmd.py     |dat files produced by deadmd.py|PNG plot| Accepts data files made by deadmd.py and saves 8 plots about the genetic algorithm runs.|
|plotgraph.py        |Various files       | Matplotlib script  | Accepts various files produced by Imago calculations and creates, then runs, a matplotlib script to plot the data.|
|pot2plot.py         |Imago scfV and structure.dat| Files containing plotable data columns| |
|processPOPTC.py     |Raw Epsilon 2 pOptc data| Control file for makePDOS| An intermediary step to make using makePDOS.py easier|

## plotgraph.py

In most cases plotgraph.py is sufficient to graph and analyze calculation outputs. It supports four types of figures and should be the default for them: general, density of states (dos), optical properties (optc), and symmetric band structure (sybd). If your data type does not fit within these categories or the program otherwise fails, look into some of the other post-processing scripts available.




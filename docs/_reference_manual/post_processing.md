---
title: Post Processing
layout: default
---

# Post Processing

> ***This page is a stub and incomplete, it is a work in progress***

Each run of Imago will generate several output files for each type of calculation. However, none of these files are plots in and of themselves. Thus, Imago contains several post-processing scripts that take the data and turn them into usable figures, or other file formats used by external tools.

***Note: This table is incomplete but contains the most commonly used scripts.*** Each script name links to its page in the Script Catalog.

**Post Processing Scripts**

|Script              | Accepted inputs    | Generate Outputs   | Description |
|:------------------:|:------------------:|:------------------:|:------------|
|[makeBOND.py](script_cat/makeBOND.html) |Bond order raw file | Files useable by external programs | Takes bond order raw file and creates many files useable in openDX, Paraview/VisIT, Origin, or other similar tools|
|[makeFittedRhoV.py](script_cat/makeFittedRhoV.html) |SCF potential and structure files; also needs potential database | input file for imagoRhoV fortran program and openDX plottable data file | Generates files to visualize various charge density and potential functions on a 3D real-space mesh.|
|[makePDOS.py](script_cat/makePDOS.html) |Control and raw PDOS data files|PDOS.plot file|Uses the commandline and control file to collect data from raw PDOS files to produce a .plot file useable by plotgraph.py
|[makeSYBD.py](script_cat/makeSYBD.html) |SYBD input, output and raw eigenvalue files|gs_sybd-fb.plot|Merges a symmetric band structure run's files into one plot file annotated with high-symmetry k-point labels. Run automatically by imago.py at the end of a -sybd job.|
|[plot\_deadmd.py](script_cat/plot_deadmd.html) |dat files produced by deadmd.py|PNG plot| Accepts data files made by deadmd.py and saves plots about the genetic algorithm runs.|
|[plotgraph.py](script_cat/plotgraph.html) |Various files       | Matplotlib, plotly or veusz script  | Accepts various files produced by Imago calculations and creates, then runs, a plotting script to display the data.|
|[pot2plot.py](script_cat/pot2plot.html) |Imago scfV and structure.dat| Files containing plotable data columns| Writes, per element, the Gaussian potential coefficients of every potential type and the potential evaluated on a radial grid.|
|[processPOPTC.py](script_cat/processPOPTC.html) |Raw Epsilon 2 pOptc data| Raw partial optical files readable by makePDOS| Runs imagoKKc on each partial epsilon 2 pair to produce partial epsilon 1, ELF, refractive index, reflectivity and absorption.|

## plotgraph.py

In most cases plotgraph.py is sufficient to graph and analyze calculation outputs. It supports four types of figures and should be the default for them: general, density of states (dos), optical properties (optc), and symmetric band structure (sybd). If your data type does not fit within these categories or the program otherwise fails, look into some of the other post-processing scripts available.




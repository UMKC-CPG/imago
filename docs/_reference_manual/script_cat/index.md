---
layout: default
title: Script Catalog
nav_order: 2
has_toc: false
---

# Script Catalog

The Script Catalog is a reference page for every script that comes
with Imago. The first section groups them by the task you are doing,
so that you can find the right tool without knowing its name. The
second lists every script alphabetically with a one-line description,
which the site search also indexes. Every script prints its own
options with `--help`, and its page here explains them.

## Scripts by task

### Obtaining a structure

| Script | Purpose |
|:-------|:--------|
|[cod_fish.py](cod_fish.html)      | Search the Crystallography Open Database by composition and download CIF files by pinned revision. |
|[cif2skl.py](cif2skl.html)        | Convert a CIF to a skeleton, preserving the space group. |
|[pdb2skl.py](pdb2skl.html), [xyz2skl.py](xyz2skl.html), [vasp2skl.py](vasp2skl.html), [struct2skl.py](struct2skl.html), [dump2skl.py](dump2skl.html) | Convert PDB, XYZ, VASP, Imago structure.dat and LAMMPS dump files to skeletons. |
|[skl2pdb.py](skl2pdb.html), [skl2vasp.py](skl2vasp.html), [skl2isaacs.py](skl2isaacs.html) | Export a skeleton to PDB, VASP POSCAR or ISAACS input. |
|[mod_struct.py](mod_struct.html)  | Supercells, translations, rotations, surfaces, vacuum, cutting, perturbation. |
|[viewCell.py](viewCell.html)      | 3D view of the cell, the Brillouin zone and the k-points. |

### Preparing and running a calculation

| Script | Purpose |
|:-------|:--------|
|[makeinput.py](makeinput.html)    | Expand a skeleton into the engine's input set and a batch script. |
|[imago.py](imago.html)            | Run the SCF and any post-SCF property calculation. |
|[tidy_scratch.py](tidy_scratch.html) | Reclaim the scratch space finished runs leave behind. |

### Post-processing

Each run of Imago generates several output files for each type of
calculation. Some are directly plottable; others must be processed
first into figures, or into formats that external tools read. In
most cases `plotgraph.py` is sufficient to plot a result: it handles
general x-y data, densities of states, optical spectra and band
structures, and is the default for all of them. If your data do not
fit those categories, or the script fails on them, look at the other
scripts below.

| Script | Accepted inputs | Generated outputs | Description |
|:-------|:----------------|:------------------|:------------|
|[plotgraph.py](plotgraph.html)          | Plottable files from Imago or any post-processing script | A matplotlib, plotly or veusz script, run for you | Plots general, DOS, optical and band-structure data with control over subplots, stacking, axes and curve styling. |
|[makePDOS.py](makePDOS.html)            | Raw PDOS file and a control file | PDOS.plot | Combines the per-atom or per-type partial DOS into the curves you ask for, for plotgraph. |
|[makeSYBD.py](makeSYBD.html)            | SYBD input, output and raw eigenvalue files | gs_sybd-fb.plot | Merges a band-structure run's files into one plot file with high-symmetry labels. Run automatically by imago.py. |
|[makeBOND.py](makeBOND.html)            | Bond order raw file | Files for openDX, Paraview, VisIt, Origin | Scatter data, structural models, bond profiles and meshes from bond orders and Q*. |
|[processPOPTC.py](processPOPTC.html)    | Raw partial epsilon-2 file | Raw partial optical files readable by makePDOS | Runs imagoKKc on each partial pair for partial epsilon-1, ELF, refractive index, reflectivity and absorption. |
|[makeFittedRhoV.py](makeFittedRhoV.html)| SCF potential, structure.dat and the potential database | Profiles along each axis and openDX volumes | Evaluates the fitted charge density and potential on a 3D mesh. |
|[pot2plot.py](pot2plot.html)            | SCF potential and structure.dat | Per-element plottable files | The Gaussian coefficients and the evaluated potential of every potential type. |
|[bond_analysis.py](bond_analysis.html)  | A skeleton | One file per analysis mode | Geometric bond lengths, angles, coordination, Q^n, rings, and model files, with no calculation needed. |
|[plot_deadmd.py](plot_deadmd.html)      | The data files of a DEAD-MD run | PNG plots | Fitness, energy, density and composition versus generation. |

### Building and curating the databases

| Script | Purpose |
|:-------|:--------|
|[unpackImagoDB.py](unpackImagoDB.html) | Unpack the shipped databases after installing. |
|[makeBDB.py](makeBDB.html), [makePotDB.py](makePotDB.html), [makeSGDB.py](makeSGDB.html), [runIsoAtoms.py](runIsoAtoms.html) | Rebuild the basis, potential, space-group and isolated-atom databases. |
|[build_initial_potentials.py](build_initial_potentials.html), [expand_manifest.py](expand_manifest.html), [makegroups.py](makegroups.html) | The initial-potential database from a curated set of reference solids. |
|[guidance_harvest.py](guidance_harvest.html), [guidance_promote.py](guidance_promote.html) | The historical guidance database that predicts k-point densities. |
|[cluster_probe.py](cluster_probe.html)  | Write a starter cluster-dispatch settings file. |

### Amorphous models by molecular dynamics

| Script | Purpose |
|:-------|:--------|
|[condense.py](condense.html)            | Generate the LAMMPS inputs for a reactive condensation of precursor molecules. |
|[make_reactions.py](make_reactions.html) | Build the bond/react templates for one pair of molecules. |
|[deadmd.py](deadmd.html), [deadmd_validity.py](deadmd_validity.html), [pdf_neutron.py](pdf_neutron.html) | The genetic-algorithm search over condensation parameters, its validity gate, and the pair distribution function it scores against. |
|[porosity.py](porosity.html), [element_pct.py](element_pct.html) | Porosity and composition of a model. |

### Relativistic basis sets and developer tools

| Script | Purpose |
|:-------|:--------|
|[graspElems.py](graspElems.html), [org_radwavefn.py](org_radwavefn.html), [make_veusz_graph.py](make_veusz_graph.html), [interFit.py](interFit.html) | Run Grasp2K, read and plot its orbitals, and fit them to Gaussians. |
|[osrecurintg.py](osrecurintg.html), [osrecurintg_makenum.py](osrecurintg_makenum.html) | Generate the Fortran source for the Gaussian integrals. |
|[hdf5Comp.py](hdf5Comp.html)            | Compare matrices between Imago HDF5 files. |
|[insert3cbo.py](insert3cbo.html)        | Three-centre bond orders. Currently not functional. |

## All scripts

{% for page in site.reference_manual %}
   {% if page.script %}
   <h3><a href="{{ page.url }}"> {{ page.script }} </a></h3>

   {{ page.description }}
   {% endif %}
{% endfor %}

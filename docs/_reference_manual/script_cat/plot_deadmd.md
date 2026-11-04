---
layout: default
title: plot_deadmd.py
script: plot_deadmd.py
description: 'Plotting utility for DEAD-MD genetic algorithm runs. Reads the per-generation data files written by deadmd.py and saves PNG plots of fitness, energy, density, element percentage and composition versus generation, optionally displaying them over X11.'
parent: Script Catalog
---

# plot_deadmd

Reads the `.dat` files that `deadmd.py` writes during a genetic
algorithm run and saves PNG plots of the run's progress to the current
directory. By default every plot is generated; a plot whose input
file is missing, as happens for an older run that did not record it,
is skipped with a note. Run it from the same directory where
`deadmd.py` was executed.

## Usage

**Commandline:**

```
plot_deadmd.py [ELEMENT] [--plots PLOT [PLOT ...]] [--show] [--help]
```

**Examples:**

```
plot_deadmd.py h                                   # all plots, hydrogen system
plot_deadmd.py h --show                            # also display over X11
plot_deadmd.py h --plots best-density best-element --show
plot_deadmd.py h --plots original-fitness original-energy \
                         original-density original-element
```

**Variable Table:**

| Variable   | Accepted Values               | Description |
|:----------:|:-----------------------------:|:------------|
|ELEMENT     |Chemical symbol, e.g. h, c, si | The tracked element, used to label the element-percentage plots. Optional: when omitted it is read from `deadmd.in`, and a mismatch with `deadmd.in` is reported. |
|PLOT        |One of the names below         | A plot to generate. |

**Option Table:**

| Option    | Outputs           | Defaults   | Description |
|:---------:|:-----------------:|:----------:|:------------|
|--plots    |Selected PNGs      |All plots   | Space-separated list of the plots to generate. |
|--show     |N/A                |Off         | Also pop up each plot interactively. Requires an active X11 display, i.e. an SSH session with `-X` or `-Y` so that `DISPLAY` is set. If no display is detected a warning is printed and the PNGs are still saved. |
|--help     |N/A                |N/A         | Print usage information and exit. |

## Available Plots

| Name                       | Input file                            | Output PNG                        | Shows |
|:---------------------------|:--------------------------------------|:----------------------------------|:------|
|original-fitness            |averages_per_gen.dat                   |fitness_per_gen.png                | Fitness: best individual vs average. |
|original-energy             |averages_per_gen.dat                   |total_energy_per_gen.png           | Total energy: best individual vs average. |
|original-density            |averages_per_gen.dat                   |density_per_gen.png                | Density: best individual vs average. |
|original-element            |averages_per_gen.dat                   |element_pct_per_gen.png            | Element percentage: best individual vs average. |
|original-energy-per-atom    |averages_per_gen.dat                   |energy_per_atom_per_gen.png        | Energy per atom (E/N): best individual vs average. |
|best-fitness                |best_by_fitness_per_gen.dat            |best_by_fitness_per_gen.png        | Fitness of the best-fitness member per generation. |
|lowest-energy               |best_by_energy_per_gen.dat             |best_by_energy_per_gen.png         | Lowest total energy per generation. |
|best-density                |best_by_density_per_gen.dat            |best_by_density_per_gen.png        | Actual vs target density of the best-density-match member, plus deviation from target. |
|best-element                |best_by_element_per_gen.dat            |best_by_element_per_gen.png        | Actual vs target element percentage of the best-element-match member, plus deviation from target. |
|lowest-energy-per-atom      |best_by_energy_per_atom_per_gen.dat    |best_by_energy_per_atom_per_gen.png| Lowest energy per atom per generation. |
|best-composition            |best_by_composition_per_gen.dat        |best_by_composition_per_gen.png    | Atomic percent of each target element of the lowest-composition-metric member vs target, plus the metric. Only present for runs with a target composition. |

The `averages_per_gen.dat` file is written by `deadmd.py`'s
per-generation averaging; the `best_by_*` files are written by its
best-quantity recording.

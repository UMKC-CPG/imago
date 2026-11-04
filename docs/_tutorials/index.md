---
# Landing page of the Tutorials collection.  `nav_order: 0` pins it
# to the top of the collection's list in the sidebar; the tutorial
# pages that follow carry their own nav_order so that they appear
# in the order a new user should read them.
layout: default
title: Start Here
nav_order: 0
---

# Tutorials

These pages walk through complete tasks from start to finish, with
every command shown and every file that appears explained. They are
meant to be followed at a terminal. The reference pages under
[Script Catalog](../reference_manual/script_cat/index.html) describe
each option of each script; the tutorials show which options you
actually reach for and why.

Read them in order if you are new to Imago.

| Tutorial | What you will have at the end |
|:---------|:------------------------------|
|[Prerequisites](prerequisites.html)          | A compiler, HDF5 and Python environment that Imago can be built against. |
|[Installing Imago](installation.html)        | A working install under your home directory, with the databases unpacked and the environment sourced. |
|[Your first calculation](first_calculation.html) | The density of states of crystalline silicon, from the shipped skeleton file to a plot on screen. |
|[A band structure](band_structure.html)      | The band structure of the same silicon crystal along the high-symmetry path, reusing the first tutorial's SCF. |

Every tutorial assumes the one before it. If a command here does not
behave as described, the most likely cause is an environment
variable that is not set; the installation tutorial's final section
lists what to check.

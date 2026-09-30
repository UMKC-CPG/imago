#!/usr/bin/env python3
## SPDX-License-Identifier: ECL-2.0
## Copyright (c) 2026 Paul Rulis

"""
=============================================================================
PLOT_DEADMD: Plotting utility for DEAD-MD genetic algorithm runs
=============================================================================

Reads the dat files produced by deadmd.py and saves PNG plots to the
current directory. By default all plots are generated (a plot whose
input file is missing, e.g. an older run, is skipped with a note). Use the
optional --plots flag to select a subset, and --show to pop up the
plots interactively via X11.

Usage:
    plot_deadmd.py [element] [--plots PLOT [PLOT ...]] [--show]

Positional arguments:
    element     Chemical symbol of the target element (e.g. h, c, si).
                Used to label the element percentage plots.

Optional arguments:
    --plots     Space-separated list of plots to generate. If omitted,
                all plots are generated. Available names:

                  original-fitness  -- fitness: best individual vs average
                  original-energy   -- total energy: best individual vs avg
                  original-density  -- density: best individual vs average
                  original-element  -- element %: best individual vs average
                  original-energy-per-atom -- energy per atom (E/N):
                                       best individual vs average

                  best-fitness      -- fitness of the best-fitness member
                  lowest-energy     -- lowest total energy per generation
                  best-density      -- actual vs target density of the
                                       best-density-match member, plus
                                       deviation from target
                  best-element      -- actual vs target element % of the
                                       best-element-match member, plus
                                       deviation from target
                  lowest-energy-per-atom -- lowest energy per atom (E/N)
                                       per generation
                  best-composition  -- at% of each target element of the
                                       lowest-composition-metric member
                                       vs target, plus the metric

    --show      Pop up each plot interactively in addition to saving the
                PNG. Requires an active X11 display (i.e. SSH with -X or
                -Y forwarding, so that the DISPLAY environment variable
                is set). If --show is passed but no display is detected,
                a warning is printed and the flag is ignored -- plots are
                still saved as PNGs.

Examples:
    # Generate all eight plots for a hydrogen system:
    plot_deadmd.py h

    # Generate all eight plots and display them interactively via X11:
    plot_deadmd.py h --show

    # Generate only the two target-match plots and display them:
    plot_deadmd.py h --plots best-density best-element --show

    # Generate the original four summary plots only:
    plot_deadmd.py h --plots original-fitness original-energy \
                             original-density original-element

Input files (must exist in the current directory):
    averages_per_gen.dat    written by record_generation_averages()
    best_by_fitness_per_gen.dat             \\
    best_by_energy_per_gen.dat               |
    best_by_energy_per_atom_per_gen.dat      |  written by
    best_by_density_per_gen.dat              |  record_best_quantities()
    best_by_element_per_gen.dat              |
    best_by_composition_per_gen.dat         /   (target_composition only)

Output files (written to the current directory):
    fitness_per_gen.png
    total_energy_per_gen.png
    density_per_gen.png
    element_pct_per_gen.png
    energy_per_atom_per_gen.png
    best_by_fitness_per_gen.png
    best_by_energy_per_gen.png
    best_by_density_per_gen.png
    best_by_element_per_gen.png
    best_by_energy_per_atom_per_gen.png
    best_by_composition_per_gen.png

Run from the same directory where deadmd.py was executed.
=============================================================================
"""

import sys
import os
import argparse
from datetime import datetime
import numpy as np

# ---------------------------------------------------------------------------
# Backend selection must happen before pyplot is imported. We check for
# --show in sys.argv here so we can make the decision early. If --show is
# requested and a DISPLAY is available, we use the interactive
# TkAgg backend (X11). Otherwise we force
# the non-interactive Agg backend which is always safe on HPC nodes.
# ---------------------------------------------------------------------------
_SHOW_REQUESTED = "--show" in sys.argv
_DISPLAY_AVAILABLE = bool(os.environ.get("DISPLAY", ""))

if _SHOW_REQUESTED and not _DISPLAY_AVAILABLE:
    print(
        "Warning: --show was requested but no DISPLAY environment variable "
        "is set. Plots will be saved as PNGs only. To enable interactive "
        "display, connect via 'ssh -X' or 'ssh -Y'."
    )
    _SHOW_REQUESTED = False

import matplotlib
if not _SHOW_REQUESTED:
    matplotlib.use("Agg")  # non-interactive, safe for HPC
else:
    # Ask for TkAgg explicitly. Left to itself matplotlib prefers QtAgg
    # when PyQt is installed, and on the cluster nodes Qt's xcb plugin
    # cannot load (libxcb-cursor0 is missing), which aborts the script.
    # If Tk is unavailable, fall back to matplotlib's own choice.
    try:
        matplotlib.use("TkAgg")
    except ImportError:
        pass

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


# ---------------------------------------------------------------------------
# All recognised plot names, in display order.
# ---------------------------------------------------------------------------

ALL_PLOTS = [
    "original-fitness",
    "original-energy",
    "original-density",
    "original-element",
    "original-energy-per-atom",
    "best-fitness",
    "lowest-energy",
    "best-density",
    "best-element",
    "lowest-energy-per-atom",
    "best-composition",
]


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Plot DEAD-MD genetic algorithm convergence results.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Available plot names:\n"
            "  original-fitness  original-energy  original-density\n"
            "  original-element  original-energy-per-atom\n"
            "  best-fitness      lowest-energy\n"
            "  best-density      best-element\n"
            "  lowest-energy-per-atom  best-composition\n"
        )
    )
    parser.add_argument(
        "element",
        type=str,
        nargs="?",
        default=None,
        help=("Chemical symbol of the tracked element (e.g. h, c, si). "
              "Optional: read from deadmd.in when omitted, and a "
              "mismatch with deadmd.in is reported.")
    )
    parser.add_argument(
        "--plots",
        nargs="+",
        choices=ALL_PLOTS,
        default=ALL_PLOTS,
        metavar="PLOT",
        help=(
            "Plots to generate (default: all). Choose from: "
            + ", ".join(ALL_PLOTS)
        )
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help=(
            "Display each plot interactively via X11 in addition to "
            "saving the PNG. Requires SSH with -X or -Y and a valid "
            "DISPLAY environment variable."
        )
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# File loading
# ---------------------------------------------------------------------------

def load_dat_file(filename, required=True):
    """
    Load a whitespace-delimited dat file with a one-line header and return
    a dictionary mapping column names to numpy arrays of float values.

    Args:
        filename: path to the dat file
        required: if False, a missing file returns None (plot skipped)

    Returns:
        dict of {column_name: np.ndarray}, or None
    """
    if not os.path.exists(filename):
        if not required:
            print(f"Skipped: '{filename}' not found (older run, or the "
                  f"option that writes it is off).")
            return None
        print(f"Error: '{filename}' not found. Run deadmd.py first.")
        sys.exit(1)

    with open(filename, "r") as dat_file:
        lines = [line.strip() for line in dat_file if line.strip()]

    if len(lines) < 2:
        print(f"Error: '{filename}' contains no data rows.")
        sys.exit(1)

    headers = lines[0].split()
    data    = {header: [] for header in headers}

    for line in lines[1:]:
        values = line.split()
        try:
            parsed = [float(v) for v in values]
        except ValueError:
            # Skip repeated header lines that appear when the dat file
            # is appended to across multiple deadmd.py runs.
            continue
        # A generation where every member crashed is written as a short
        # row of nan; pad it so all columns keep the same length.
        parsed += [float("nan")] * (len(headers) - len(parsed))
        for header, value in zip(headers, parsed):
            data[header].append(value)

    return {key: np.array(vals) for key, vals in data.items()}


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def style_axes(axes, title, xlabel, ylabel):
    """Apply consistent title, axis labels, grid, and legend to an Axes."""
    axes.set_title(title, fontsize=13, fontweight="bold", pad=10)
    axes.set_xlabel(xlabel, fontsize=11)
    axes.set_ylabel(ylabel, fontsize=11)
    axes.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    axes.grid(True, linestyle="--", alpha=0.5)
    axes.legend(fontsize=10)
    axes.tick_params(labelsize=10)


def save_figure(fig, filename, show=False):
    """
    Save a figure to disk and optionally display it interactively.

    When show=True the window is opened non-blocking so that all plots
    appear at once. A single blocking pause is issued at the end of
    main() after all figures have been drawn, keeping every window open
    until the user dismisses them.

    Args:
        fig:      matplotlib Figure object
        filename: output PNG path
        show:     if True, pop up the window via X11 (non-blocking)
    """
    fig.savefig(filename, dpi=150, bbox_inches="tight")
    print(f"Saved: {filename}")
    if show:
        plt.show(block=False)
        plt.pause(0.1)  # let the window event loop initialise


# Minimum y-axis spans. Without them matplotlib zooms to the data, so a
# 0.02 g/cm3 or 1 at% difference fills the whole panel and looks huge.
DENSITY_SPAN = 0.2   # g/cm3
ELEMENT_SPAN = 5.0   # at%


def set_min_span(axes, span, from_zero=False):
    """Widen the y range to at least `span`: centred on the data, or
    starting at 0 (deviation panels) when from_zero is True."""
    low, high = axes.get_ylim()
    if from_zero:
        axes.set_ylim(0.0, max(high, span))
    elif high - low < span:
        middle = 0.5 * (low + high)
        axes.set_ylim(middle - 0.5 * span, middle + 0.5 * span)


def element_from_deadmd_in(filename="deadmd.in"):
    """
    Return the element deadmd.py tracks in element_evolve and the
    best_by_element / element_pct columns: the `target_element` line if
    present, else the first element of `target_composition`. None if
    deadmd.in is missing or has neither.
    """
    if not os.path.exists(filename):
        return None
    first_composition = None
    with open(filename) as in_file:
        for line in in_file:
            words = line.split("#")[0].split()
            if len(words) >= 2 and words[0] == "target_element":
                return words[1].lower()
            if (len(words) >= 2 and words[0] == "target_composition"
                    and first_composition is None):
                first_composition = words[1].lower()
    return first_composition


# ---------------------------------------------------------------------------
# Original four plots (best-by-fitness vs population average)
# ---------------------------------------------------------------------------

def plot_original_fitness(generations, best_fitness, avg_fitness, show):
    """
    Plot best-by-fitness individual and population average fitness
    per generation. Lower fitness is better (minimisation problem).
    """
    fig, axes = plt.subplots(figsize=(8, 5))

    axes.plot(generations, best_fitness, marker="o", linewidth=2,
              color="steelblue", label="Best individual")
    axes.plot(generations, avg_fitness, marker="s", linewidth=2,
              linestyle="--", color="coral", label="Population average")

    style_axes(axes,
               title="Fitness per Generation",
               xlabel="Generation",
               ylabel="Fitness (lower is better)")

    fig.tight_layout()
    save_figure(fig, "fitness_per_gen.png", show)


def plot_original_energy(generations, best_energy, avg_energy, show):
    """
    Plot total energy of the best-by-fitness individual and population
    average per generation.
    """
    fig, axes = plt.subplots(figsize=(8, 5))

    axes.plot(generations, best_energy, marker="o", linewidth=2,
              color="steelblue", label="Best individual")
    axes.plot(generations, avg_energy, marker="s", linewidth=2,
              linestyle="--", color="coral", label="Population average")

    style_axes(axes,
               title="Total Energy per Generation",
               xlabel="Generation",
               ylabel="Total Energy (kcal/mol)")

    fig.tight_layout()
    save_figure(fig, "total_energy_per_gen.png", show)


def average_energy_per_atom(avgs, n_gens):
    """
    Population-average E/N per generation. Runs from 2026-09-30 on write
    it as avg_energy_per_atom in averages_per_gen.dat; for older runs it
    is recomputed from generation_<g>/<m>/lammps/{totE,natoms}_evolve,
    skipping crashed members (sentinel >= 1e90) as deadmd.py does.
    """
    if "avg_energy_per_atom" in avgs:
        return avgs["avg_energy_per_atom"]
    averages = []
    for gen in range(1, n_gens + 1):
        values = []
        member = 1
        while os.path.isdir(f"generation_{gen}/{member}"):
            lammps_dir = f"generation_{gen}/{member}/lammps"
            try:
                with open(f"{lammps_dir}/totE_evolve") as e_file:
                    energy = float(e_file.read().split()[-1])
                with open(f"{lammps_dir}/natoms_evolve") as n_file:
                    natoms = float(n_file.read().split()[-1])
            except (OSError, ValueError, IndexError):
                energy = natoms = 1e99
            if energy < 1e90 and natoms < 1e90:
                values.append(energy / natoms)
            member += 1
        averages.append(np.mean(values) if values else float("nan"))
    return np.array(averages)


def plot_original_energy_per_atom(generations, best_epa, avg_epa, show):
    """
    Plot energy per atom (E/N, the quantity the fitness energy metric
    uses) of the best-by-fitness individual and population average per
    generation.
    """
    fig, axes = plt.subplots(figsize=(8, 5))

    axes.plot(generations, best_epa, marker="o", linewidth=2,
              color="steelblue", label="Best individual")
    axes.plot(generations, avg_epa, marker="s", linewidth=2,
              linestyle="--", color="coral", label="Population average")

    style_axes(axes,
               title="Energy per Atom per Generation",
               xlabel="Generation",
               ylabel="Energy per atom (kcal/mol/atom)")

    fig.tight_layout()
    save_figure(fig, "energy_per_atom_per_gen.png", show)


def plot_original_density(generations, best_density, avg_density, show):
    """
    Plot density of the best-by-fitness individual and population
    average per generation.
    """
    fig, axes = plt.subplots(figsize=(8, 5))

    axes.plot(generations, best_density, marker="o", linewidth=2,
              color="steelblue", label="Best individual")
    axes.plot(generations, avg_density, marker="s", linewidth=2,
              linestyle="--", color="coral", label="Population average")

    style_axes(axes,
               title="Density per Generation",
               xlabel="Generation",
               ylabel="Density (g/cm³)")
    set_min_span(axes, DENSITY_SPAN)

    fig.tight_layout()
    save_figure(fig, "density_per_gen.png", show)


def plot_original_element(generations, best_pct, avg_pct, element, show):
    """
    Plot element percentage of the best-by-fitness individual and
    population average per generation.
    """
    element_label = element.capitalize()
    fig, axes = plt.subplots(figsize=(8, 5))

    axes.plot(generations, best_pct, marker="o", linewidth=2,
              color="steelblue", label="Best individual")
    axes.plot(generations, avg_pct, marker="s", linewidth=2,
              linestyle="--", color="coral", label="Population average")

    style_axes(axes,
               title=f"{element_label} Percentage per Generation",
               xlabel="Generation",
               ylabel=f"{element_label} Percentage (%)")
    set_min_span(axes, ELEMENT_SPAN)

    fig.tight_layout()
    save_figure(fig, "element_pct_per_gen.png", show)


# ---------------------------------------------------------------------------
# New four plots (best by independent criteria)
# ---------------------------------------------------------------------------

def plot_best_fitness(generations, fitness, show):
    """
    Plot the fitness score of the best-fitness member per generation.
    This member is selected solely by lowest fitness, independent of
    how well it matches target density or element percentage.
    """
    fig, axes = plt.subplots(figsize=(8, 5))

    axes.plot(generations, fitness, marker="o", linewidth=2,
              color="steelblue", label="Best-fitness member")

    style_axes(axes,
               title="Best Fitness per Generation",
               xlabel="Generation",
               ylabel="Fitness (lower is better)")

    fig.tight_layout()
    save_figure(fig, "best_by_fitness_per_gen.png", show)


def plot_lowest_energy(generations, energy, show):
    """
    Plot the lowest total energy found in each generation. The member
    with the lowest energy may differ from the best-fitness member since
    fitness also includes density and element percentage contributions.
    """
    fig, axes = plt.subplots(figsize=(8, 5))

    axes.plot(generations, energy, marker="o", linewidth=2,
              color="darkgreen", label="Lowest-energy member")

    style_axes(axes,
               title="Lowest Total Energy per Generation",
               xlabel="Generation",
               ylabel="Total Energy (kcal/mol)")

    fig.tight_layout()
    save_figure(fig, "best_by_energy_per_gen.png", show)


def plot_best_density(generations, actual, target, deviation, show):
    """
    Two-panel plot for the member whose actual density is closest to its
    own genome target density each generation.

    Top panel:  actual density vs target density on the same axes,
                showing how closely the best-match member achieved its
                target.
    Bottom panel: absolute deviation (|actual - target|) vs generation,
                  showing whether the GA is converging toward target
                  density over time.
    """
    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(8, 8),
                                          sharex=True)

    ax_top.plot(generations, actual, marker="o", linewidth=2,
                color="steelblue", label="Actual density")
    ax_top.plot(generations, target, marker="s", linewidth=2,
                linestyle="--", color="coral", label="Target density")

    style_axes(ax_top,
               title="",
               xlabel="",
               ylabel="Density (g/cm³)")

    ax_bot.plot(generations, deviation, marker="^", linewidth=2,
                color="purple", label="|actual - target|")

    style_axes(ax_bot,
               title="",
               xlabel="Generation",
               ylabel="Density Deviation (g/cm³)")
    ax_bot.set_title("")
    set_min_span(ax_top, DENSITY_SPAN)
    set_min_span(ax_bot, DENSITY_SPAN / 2, from_zero=True)

    fig.suptitle("Best Density Match per Generation",
                 fontsize=14, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, "best_by_density_per_gen.png", show)


def plot_best_element(generations, actual, target, deviation, element, show):
    """
    Two-panel plot for the member whose actual element percentage is
    closest to its own genome target element percentage each generation.

    Top panel:  actual percentage vs target percentage on the same axes.
    Bottom panel: absolute deviation (|actual - target|) vs generation.
    """
    element_label = element.capitalize()
    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(8, 8),
                                          sharex=True)

    ax_top.plot(generations, actual, marker="o", linewidth=2,
                color="steelblue",
                label=f"Actual {element_label} %")
    ax_top.plot(generations, target, marker="s", linewidth=2,
                linestyle="--", color="coral",
                label=f"Target {element_label} %")

    style_axes(ax_top,
               title="",
               xlabel="",
               ylabel=f"{element_label} Percentage (%)")

    ax_bot.plot(generations, deviation, marker="^", linewidth=2,
                color="purple", label="|actual - target|")

    style_axes(ax_bot,
               title="",
               xlabel="Generation",
               ylabel=f"{element_label} % Deviation")
    ax_bot.set_title("")

    set_min_span(ax_top, ELEMENT_SPAN)
    set_min_span(ax_bot, ELEMENT_SPAN / 2, from_zero=True)

    fig.suptitle(
        f"Best {element_label} Percentage Match per Generation",
        fontsize=14, fontweight="bold"
    )
    fig.tight_layout()
    save_figure(fig, "best_by_element_per_gen.png", show)


def plot_lowest_energy_per_atom(generations, energy_per_atom, show):
    """
    Plot the lowest energy per atom (E/N) in each generation, the
    quantity the fitness energy metric uses. Its member can differ from
    the lowest-total-energy member because members keep different atom
    counts.
    """
    fig, axes = plt.subplots(figsize=(8, 5))

    axes.plot(generations, energy_per_atom, marker="o", linewidth=2,
              color="darkgreen", label="Lowest-E/N member")

    style_axes(axes,
               title="Lowest Energy per Atom per Generation",
               xlabel="Generation",
               ylabel="Energy per atom (kcal/mol/atom)")

    fig.tight_layout()
    save_figure(fig, "best_by_energy_per_atom_per_gen.png", show)


def plot_best_composition(data, show):
    """
    Two-panel plot for the member with the lowest composition metric
    each generation (best_by_composition_per_gen.dat).

    Top panel:  at% of each target element (solid) and its target
                (dashed, same colour; the target is in the column
                header, e.g. 'si(9.43)').
    Bottom panel: the composition metric
                  (1/n) sum_El |x_El - x_target,El| / sigma_El.
    """
    generations = data["gen"].astype(int)
    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(8, 8),
                                          sharex=True)

    element_columns = [c for c in data
                       if c not in ("gen", "member", "comp_metric")]
    for column in element_columns:
        element, _, target = column.partition("(")
        line, = ax_top.plot(generations, data[column], marker="o",
                            linewidth=2, label=f"{element.capitalize()} %")
        if target:
            ax_top.axhline(float(target.rstrip(")")), linestyle="--",
                           color=line.get_color(),
                           label=f"{element.capitalize()} target")

    style_axes(ax_top, title="", xlabel="",
               ylabel="Composition (at%)")
    ax_top.legend(fontsize=9, ncol=2)

    ax_bot.plot(generations, data["comp_metric"], marker="^",
                linewidth=2, color="purple", label="Composition metric")
    style_axes(ax_bot, title="", xlabel="Generation",
               ylabel="Composition metric (lower is better)")

    fig.suptitle("Best Composition Match per Generation",
                 fontsize=14, fontweight="bold")
    fig.tight_layout()
    save_figure(fig, "best_by_composition_per_gen.png", show)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args    = parse_args()
    element = args.element
    tracked = element_from_deadmd_in()
    if element is None:
        if tracked is None:
            print("Error: no element given and none found in deadmd.in "
                  "(target_element or target_composition).")
            sys.exit(1)
        element = tracked
    elif tracked is not None and element.lower() != tracked:
        print(f"Warning: you gave '{element}', but deadmd.in makes deadmd.py "
              f"track '{tracked}' (target_element, else the first element "
              f"of target_composition). Labelling as '{tracked}'.")
        element = tracked
    plots   = set(args.plots)
    show    = _SHOW_REQUESTED  # resolved at import time against DISPLAY

    # Load only the files needed for the requested plots.
    needs_original = any(p.startswith("original-") for p in plots)

    if needs_original:
        best = load_dat_file("best_by_fitness_per_gen.dat")
        avgs = load_dat_file("averages_per_gen.dat")
        if len(best["gen"]) != len(avgs["gen"]):
            print(
                f"Error: best_by_fitness_per_gen.dat has "
                f"{len(best['gen'])} rows but averages_per_gen.dat has "
                f"{len(avgs['gen'])} rows. The dat files are from "
                f"different runs -- delete them all and rerun deadmd.py."
            )
            sys.exit(1)
        generations_orig = best["gen"].astype(int)

    # --- Original plots ---
    if "original-fitness" in plots:
        plot_original_fitness(
            generations_orig,
            best_fitness=best["fitness"],
            avg_fitness=avgs["avg_fitness"],
            show=show
        )

    if "original-energy" in plots:
        plot_original_energy(
            generations_orig,
            best_energy=best["total_energy"],
            avg_energy=avgs["avg_total_energy"],
            show=show
        )

    if "original-density" in plots:
        plot_original_density(
            generations_orig,
            best_density=best["density"],
            avg_density=avgs["avg_density"],
            show=show
        )

    if "original-element" in plots:
        plot_original_element(
            generations_orig,
            best_pct=best["element_pct"],
            avg_pct=avgs["avg_element_pct"],
            element=element,
            show=show
        )

    if "original-energy-per-atom" in plots:
        plot_original_energy_per_atom(
            generations_orig,
            best_epa=best["energy_per_atom"],
            avg_epa=average_energy_per_atom(avgs, len(generations_orig)),
            show=show
        )

    # --- New plots ---
    if "best-fitness" in plots:
        bf = load_dat_file("best_by_fitness_per_gen.dat")
        plot_best_fitness(
            generations=bf["gen"].astype(int),
            fitness=bf["fitness"],
            show=show
        )

    if "lowest-energy" in plots:
        le = load_dat_file("best_by_energy_per_gen.dat")
        plot_lowest_energy(
            generations=le["gen"].astype(int),
            energy=le["total_energy"],
            show=show
        )

    if "best-density" in plots:
        bd = load_dat_file("best_by_density_per_gen.dat")
        plot_best_density(
            generations=bd["gen"].astype(int),
            actual=bd["actual_density"],
            target=bd["target_density"],
            deviation=bd["deviation"],
            show=show
        )

    if "best-element" in plots:
        be = load_dat_file("best_by_element_per_gen.dat")
        plot_best_element(
            generations=be["gen"].astype(int),
            actual=be["actual_pct"],
            target=be["target_pct"],
            deviation=be["deviation"],
            element=element,
            show=show
        )

    if "lowest-energy-per-atom" in plots:
        epa = load_dat_file("best_by_energy_per_atom_per_gen.dat",
                            required=False)
        if epa is not None:
            plot_lowest_energy_per_atom(
                generations=epa["gen"].astype(int),
                energy_per_atom=epa["energy_per_atom"],
                show=show
            )

    if "best-composition" in plots:
        bc = load_dat_file("best_by_composition_per_gen.dat",
                           required=False)
        if bc is not None:
            plot_best_composition(bc, show=show)

    print("Done.")
    if show:
        # Hold the figures on screen until the viewer closes them.
        #   Each figure was raised with plt.show(block=False) as it was
        #   drawn (see save_figure), which returns immediately, so
        #   without something blocking here the interpreter would exit
        #   and take every window with it.  A final blocking plt.show()
        #   is matplotlib's own idiom for that wait, and it makes
        #   closing the windows the way to end the session -- no
        #   terminal keystroke needed, and no reading from standard
        #   input anywhere in this script.
        print("Close the plot windows to exit.")
        plt.show()


def record_command():
    """Append the issued command line to a file named "command" in
    the current directory, so the exact invocation can be recovered
    later.  This is a standing project convention: each run appends
    a dated block, so the file builds up a history of how the script
    was called."""

    with open("command", "a") as cmd:
        now = datetime.now()
        stamp = now.strftime("%b. %d, %Y: %H:%M:%S")
        cmd.write(f"Date: {stamp}\n")
        cmd.write("Cmnd:")
        for argument in sys.argv:
            cmd.write(f" {argument}")
        cmd.write("\n\n")


if __name__ == "__main__":
    record_command()
    main()

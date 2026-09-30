#!/usr/bin/env python3
## SPDX-License-Identifier: ECL-2.0
## Copyright (c) 2026 Mohammed Belhadj Larbi

"""
=============================================================================
DEAD_MD: Dynamic Evolutionary Algorithm-Driven Molecular Dynamics
=============================================================================

Version history: see git log.

Author:      Mohammed Belhadj Larbi
Institution: University of Missouri-Kansas City (UMKC)

-----------------------------------------------------------------------------
Description
-----------------------------------------------------------------------------
DEAD-MD is a genetic algorithm (GA) that searches the parameters of a
condense.py + LAMMPS fix bond/react condensation (cell size, compression
rate, per-reaction probabilities, Rmax and caps) for amorphous models that
match a target: experimental composition and density (replica runs), or
low energy at a chosen composition range (exploratory runs).  Every model is
built through reaction chemistry, so every member is a physically reachable
network.

-----------------------------------------------------------------------------
Algorithm
-----------------------------------------------------------------------------
1. Draw an initial population of genomes (seeded by ga_seed).
2. For every member: write condense.in from its genome, run condense.py
   (seeded per member) and LAMMPS.  Members run CONCURRENTLY on num_cores
   cores, each in its own srun step.
3. Read each member's results (energy, atom count, density, composition,
   bond fractions, optional G(r) and validity check) and compute its fitness
   (lower is better).  A crashed or invalid member gets the sentinel 1e99.
4. Select parents (truncation / rank / tournament k), make children by 1- or
   2-point crossover, mutate (uniform redraw or Gaussian creep).
5. The next generation = elites (kept, NOT re-simulated: copied and
   re-scored) + children.
6. Repeat until gen_max.  Convergence criteria (stop when the best fitness
   stops improving by more than the seed noise, or when the target is
   reached) will be added in the near future; until then the run always
   goes to gen_max.

-----------------------------------------------------------------------------
Genome (one row per member)
-----------------------------------------------------------------------------
    gene 0              target element at%  (single-element metric only)
    gene 1              target density: condense.py stops compressing there.
                        In replica runs set lower = upper (the experimental
                        density): it is then a fixed stop value, not searched.
    gene 2              cell_size (A)
    gene 3              squish factor of the genome stage (condensation rate)
    genes 4 .. 3+n      Rmax of reaction r (n = reactions_num)
    genes 4+n .. 3+2n   probability of reaction r
    genes 4+2n ..       max_rxn cap of each reaction that has "max_rxn LO HI",
                        in line order (whole numbers)
Every gene is drawn and mutated inside its own range; lower = upper fixes it.

-----------------------------------------------------------------------------
Fitness (lower is better)
-----------------------------------------------------------------------------
    fitness = norm * sum_i (w_i * metric_i),  norm = n_active / sum_i w_i
Metrics (each dimensionless, 0 = on target):
    energy       (E/N - lowest E/N of generation 1) / energy_tolerance
                 (signed: members below the reference score negative)
    density      |rho - gene 1| / gene 1  (only checks that the density
                 stop was reached; ~0 whenever it was)
    element      |at% - gene 0| / gene 0  (legacy single-element metric;
                 off when only target_composition is used)
    composition  (1/n) sum_El |x - x_target| / sigma_El    (optional)
    bond         (1/n) sum max(0, f - max_pct) / tolerance (optional)
    pdf          scale-optimised R-factor vs exp_pdf_file  (optional)
E = minimized total energy, N = final atom count, rho = density, x = at%.
w_i are the weight_* keywords (default 1.0); norm keeps the fitness scale
the same whatever the number of active metrics.

-----------------------------------------------------------------------------
Output
-----------------------------------------------------------------------------
Run directory: run_stamp/ (deadmd.py, condense.py, deadmd.in, ga_seed);
averages_per_gen.dat; best_by_{fitness,energy,energy_per_atom,density,element,composition,
pdf}_per_gen.dat; generation_<g>/fitness_scores_gen.dat (all members, ranked).
Member directory generation_<g>/<m>/: condense.in, seeds, member.out,
ELITE_COPY_OF (elite copies only) and lammps/ with final.data and the
*_evolve files (totE, natoms, density, composition, rxn_counts, bond_limit,
fitness_scores_gen, ...).

Reproducibility: the same ga_seed + deadmd.in + LAMMPS binary + ranks per
member (cores_per_member) give the same run; a different binary or rank
count does not.

----------------------------------------------------------------------------
Example deadmd.in
----------------------------------------------------------------------------
The following shows a complete annotated input file.  Every keyword is
case-insensitive.  Lines are parsed by splitting on whitespace; order
within the file does not matter.

    num_cores          4

    # Members run concurrently on the num_cores cores (num_cores must not
    # exceed the job's #SBATCH -n). By default the cores are shared out
    # automatically: with at least num_cores members to simulate, each gets
    # 1 core and num_cores run at once; with fewer, all run at once and
    # share the cores (10 cores, 4 members -> 3 3 2 2). cores_per_member
    # fixes the MPI ranks of every member instead (optional); a-BC:H runs
    # with inter_unlinked need 1 (serial-exact only, see DIAGNOSIS pin 9).
    # Each member's screen output (condense, LAMMPS) goes to member.out in
    # its directory, not to deadmd.o.
    cores_per_member   1

    # lmp picks the LAMMPS binary: parallel (lmp_parallel, MPI; default) or
    # serial (lmp, one core per member; forces cores_per_member 1). Serial
    # is exact for bond/react inter_unlinked (DIAGNOSIS pin 9).
    lmp                parallel

    # ga_seed seeds the GA's own random numbers (initial genomes, selection,
    # crossover, mutation) and, through a per-member seed derived from
    # (ga_seed, generation, member), every packmol / bond/react / velocity
    # draw in condense.py. Same ga_seed + same deadmd.in + same LAMMPS
    # binary + same cores_per_member -> the same run. Optional; when omitted a seed is drawn and
    # printed in deadmd.o, and either way it is kept in run_stamp/ga_seed.
    ga_seed            12345

    gen_max            20
    population_size    10
    # elitism_rate, mutation_rate and kill_rate are PERCENTS: % of the
    # population kept as elites (at least 1 if > 0), % chance per gene of
    # a mutation, % of the population not allowed to be parents.
    elitism_rate       20
    mutation_rate      5
    # Optional mutation operator (default uniform = redraw over the range):
    # mutation_type creep moves a mutated gene by a Gaussian step of sd
    # creep_sigma x (upper - lower); creep_jump % of mutations are still a
    # uniform redraw.
    mutation_type      creep
    creep_sigma        0.1
    creep_jump         10

    # Optional bond metric: % of El1 atoms bonded to >= 1 El2 atom
    # must stay <= max_pct; metric = mean of max(0, f - max_pct)/tolerance
    # (tolerance default 1 %). E.g. the a-SiCN:H C-C ceiling:
    bond_limit         c c 9.0 tolerance 1.0
    weight_bond        1.0
    kill_rate          4
    crossover_type     1
    selection          tournament 3

    target_element     h
    target_element_lower  30.0
    target_element_upper  40.0

    target_density_lower  1.8
    target_density_upper  2.4

    # Optional composition metric (replica runs): the member's at% of each
    # listed element against the target, in units of a tolerance sigma,
    #     composition metric = (1/n) * sum over El of |x - x_target| / sigma
    # (n = number of listed elements; 1 = on average every element is off
    # by one tolerance). List any subset of elements (e.g. leave out N when
    # the precursor fixes it); the values must not sum to more than 100.
    # Default sigma: 1 at% per element, 2 at% for H; override any
    # with composition_tolerance. With target_composition the three
    # target_element lines above become optional; left out, the
    # single-element metric is off. weight_composition defaults to 1.0.
    target_composition     si 9.43 c 21.61 n 6.67 h 62.29
    composition_tolerance  h 2.0 si 1.0
    weight_composition     1.0

    composition_num    2
    c6h19nsi2_1  c-1  50
    c2h6_1       c-1  10

    cell_size_lower    20.0
    cell_size_upper    40.0

    # Optional range of gene 3, the condensation rate (squish factor: the
    # final box edge fraction of a genome stage).  Used for the initial draw
    # AND for mutation.  Default 0.15 0.33.  A member can reach its target
    # density only if (starting density) / squish^3 >= target.
    squish_lower       0.15
    squish_upper       0.33

    max_speed          5.0

    # The energy metric needs no reference energy: each member's energy per
    # atom (total energy / final atom count, since reactions delete atoms)
    # is compared with the LOWEST energy per atom of generation 1:
    #     energy metric = (E/N - best gen-1 E/N) / energy_tolerance
    # The reference is generation 1's lowest energy per atom, fixed for the
    # run: that member scores 0, one energy_tolerance higher per atom scores
    # +1, and later members that beat it score NEGATIVE (in proportion), so
    # lower energy always lowers the fitness. No absolute value is taken.
    # energy_tolerance is in the LAMMPS energy unit per atom (kcal/mol per
    # atom for units real); default 0.5.
    energy_tolerance   0.5
    # ref_energy is no longer used (a leftover line is accepted and ignored).

    # weight_energy, weight_density and weight_element (and weight_composition,
    # weight_bond, weight_pdf) scale the fitness metrics relative to each
    # other.  All default to 1.0.  Only their RATIOS change the ranking.
    # Set one to 0.0 to ignore that metric.
    weight_energy      1.0
    weight_density     2.0
    weight_element     1.5

    # exp_pdf_file provides the path to a reference reduced PDF G(r).  The
    # file must contain two whitespace-separated COLUMNS: the first column is
    # r in Angstrom (0.01 spacing, strictly increasing, reaching at least
    # 10.00), the second is G(r).  When this keyword is present and
    # weight_pdf is greater than zero, a scale-optimised R-factor between
    # each member's G(r) and this reference is added to the fitness as a metric.
    # Omitting either keyword disables the PDF metric entirely.
    # Each member's G(r) is computed by pdf_neutron.py (normalised,
    # G = 4 pi rho r (g - 1)), out to half the member's box width.  For a
    # CALCULATED reference, make it with pdf_neutron.py on the reference
    # structure with the same pdf_weighting and pdf_sigma, e.g.
    #     pdf_neutron.py -i ref.skl -o ref_G.dat -rmax 12 -full
    # (-full so the file reaches 10 A).  Old rpdf output is a different,
    # unnormalised quantity and must not be used as the reference.
    exp_pdf_file       ref_G.dat
    weight_pdf         1.0
    # pdf_weighting: "neutron" (bound coherent scattering lengths; H-X
    # peaks come out negative) or "none" (every pair counts 1, i.e. the
    # geometric PDF).  Default neutron.  pdf_sigma: Gaussian broadening in
    # Angstrom standing in for thermal motion and Q resolution; default 0.05.
    pdf_weighting      neutron
    pdf_sigma          0.05

    # Optional validity gate (deadmd_validity.py, run on each member's
    # final.data).  A member that fails is INVALID, not worse, and gets
    # the crash sentinel.  Hard checks: every atom's bond count is allowed
    # for its element; final atoms = start - deleted atoms per event x
    # events; every bond within +-max_strain of its r0.  valence lines
    # override the default allowed bond counts (H 1, C 4, N 3, O 2, Si 4,
    # B 3-6, ...), e.g. a-BC:H cage carbon.  Default: off.
    validity_gate      yes
    valence            C 4 6
    max_strain         0.30

    # One stage line: squish_step_size ensemble t_start t_end t_damp run_steps
    # The condensation rate (squish_factor) is taken from the member's genome;
    # compression stops when the member's target density (gene 1) is reached,
    # then the stage finishes its steps at the fixed box.  Use a small
    # squish_step_size (e.g. 10) so the density stop is precise.
    stage              100 nvt 300.0 300.0 100.0 50000

    # Additional stage lines (optional) include squish_factor explicitly as
    # the first value.  Use 1.0 for an equilibration step (no compression);
    # every input should end with such a hold stage.
    stage              1.0 100 nvt 300.0 1000.0 100.0 100000


    # One line per reaction:
    #   mol1 site1 mol2 site2 prob_lower prob_upper [rmax_lower rmax_upper]
    # Each reaction gets its own Rmax gene and its own probability gene,
    # drawn (and mutated) within that reaction's own ranges.  Rmax (A) is
    # fix bond/react's largest initiator distance; without the two optional
    # values it stays at 5.0.  Every line then ends with inter (the bond
    # forms between two molecules) or intra (within one, e.g. si-h),
    # optionally followed by rmin R (smallest initiator distance, A,
    # default 0.0).  These are not genes; they go to condense.in as is.
    # Optional "max_rxn LO HI" (whole numbers) makes fix bond/react's cap on
    # that reaction a gene: each member draws a cap in [LO, HI] and the
    # reaction stops after that many events (a cap, never a quota).
    # "max_rxn N" (or LO = HI) fixes it.  Without it the reaction is
    # uncapped.  The cap genes sit at the END of the genome (after the
    # probabilities), one per capped reaction, in line order.
    reactions_num 3
    c6h19nsi2_1 c-1 c6h19nsi2_1 c-1 0.2 0.8 3.5 5.0 inter max_rxn 2 10
    c6h19nsi2_1 c-1 c2h6_1 c-1 0.5 0.6 inter
    c6h19nsi2_1 si-1 c6h19nsi2_1 h-1 1e-4 1e-3 intra rmin 2.0

"""

import os
import re
import sys
from datetime import datetime
import numpy as np
import random
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor


# The run's GA seed (set in initialize_population); member_seed derives
# each member's condense.py seed from it.
GA_SEED = None

# Mutation settings (set in initialize_population from deadmd.in).
#   MUTATION_TYPE "uniform": a mutated gene is redrawn uniformly over its
#                 whole range (the original behaviour, default).
#   MUTATION_TYPE "creep":   a mutated gene is moved by a Gaussian step,
#                 x' = x + N(0, CREEP_SIGMA * (upper - lower)), reflected
#                 back into [lower, upper]; with probability CREEP_JUMP it
#                 is redrawn uniformly instead (keeps some exploration).
MUTATION_TYPE = "uniform"
CREEP_SIGMA   = 0.1
CREEP_JUMP    = 0.1

# max_rxn genes (set in initialize_population from deadmd.in). One entry
# (reaction index r, lower, upper) per reaction line carrying "max_rxn LO HI",
# in line order; entry j is gene 4 + 2 * reactions_num + j, a whole number
# written to that reaction's condense.in line as "max_rxn N". Reactions
# without it are uncapped and have no gene, so older inputs keep their
# genome (and their random draws) unchanged.
MAXRXN_GENES = []

# Bond-limit metric (set in initialize_population from deadmd.in).
# BOND_LIMITS: list of (el1, el2, max_pct, tolerance), one per
# "bond_limit El1 El2 max_pct [tolerance t]" line. For each, f = % of El1
# atoms bonded to at least one El2 atom (from final.data), and the member's
# bond metric = (1/n) sum max(0, f - max_pct) / tolerance. Empty = off.
BOND_LIMITS = []
WEIGHT_BOND = 1.0

# Element from a LAMMPS mass (final.data has no element names); same table
# as deadmd_validity.py.
MASS_TO_ELEMENT = {
    1.008: "h", 10.811: "b", 12.011: "c", 14.007: "n", 15.999: "o",
    18.998: "f", 28.086: "si", 30.974: "p", 32.06: "s", 35.45: "cl",
    72.63: "ge",
}


def bond_fractions(data_file, pairs):
    """
    For each (el1, el2) in pairs, the percentage of el1 atoms that have at
    least one bond to an el2 atom, read from a LAMMPS data file (Masses,
    Atoms, Bonds sections). For el1 == el2 (e.g. c c) this is the % of C
    atoms in at least one C-C bond, the reading settled for the C-C
    ceiling (jobs/a-SiCN:H/DIAGNOSIS.md, 2026-09-09).
    """
    mass, atype, bonds, section = {}, {}, [], None
    with open(data_file) as f:
        for line in f:
            w = line.split()
            if not w:
                continue
            if w[0] in ("Masses", "Atoms", "Bonds", "Velocities", "Angles",
                        "Dihedrals", "Impropers") or "Coeffs" in w[:2]:
                section = w[0]
                continue
            if not w[0].lstrip("-").isdigit():
                continue
            if section == "Masses":
                mass[int(w[0])] = float(w[1])
            elif section == "Atoms":          # atom_style full: id mol type ...
                atype[int(w[0])] = int(w[2])
            elif section == "Bonds":
                bonds.append((int(w[2]), int(w[3])))
    elem_of_type = {}
    for t, m in mass.items():
        best = min(MASS_TO_ELEMENT, key=lambda x: abs(x - m))
        if abs(best - m) > 0.1:
            raise ValueError(f"no element with mass {m} in {data_file}")
        elem_of_type[t] = MASS_TO_ELEMENT[best]
    elem = {a: elem_of_type[t] for a, t in atype.items()}
    partners = {a: set() for a in elem}
    for a, b in bonds:
        partners[a].add(elem[b])
        partners[b].add(elem[a])
    result = {}
    for el1, el2 in pairs:
        atoms = [a for a in elem if elem[a] == el1]
        hit = sum(1 for a in atoms if el2 in partners[a])
        result[(el1, el2)] = 100.0 * hit / len(atoms) if atoms else 0.0
    return result


def member_seed(gen_count, member):
    """
    Seed for one member's condense.py run, fixed by (ga_seed,
    generation, member) alone. It does not draw from the GA's random
    stream, so it cannot shift the GA's own choices, and two members never
    share packing / reaction / velocity draws.
    """
    return random.Random(f"{GA_SEED}-{gen_count}-{member}").randint(
        1, 2**31 - 1)


def stamp_run():
    """
    Copy the code and input that made this run into run_stamp/, with the
    GA seed, so the run can be read or repeated after the scripts change.
    """
    os.makedirs("run_stamp", exist_ok=True)
    for src in (os.path.abspath(__file__), shutil.which("condense.py"),
                "deadmd.in"):
        if src and os.path.exists(src):
            shutil.copy2(src, "run_stamp")
    with open("run_stamp/ga_seed", "w") as f:
        f.write(f"{GA_SEED}\n")


def print_banner():
        print("\n")
        print("==========================================================")
        print(" DEAD-MD: Dynamic Evolutionary Algorithm-Driven Molecular Dynamics")
        print(" Version: 1.0")
        print(" Author:  Mohammed Belhadj Larbi")
        print(" Institution: University of Missouri-Kansas City (UMKC)")
        print("==========================================================")
        print("\n")


def initialize_population():
    """ 
    Reads the input file deadmd.in, and then based on that it will create 
    an initial population of N LAMMPS simulations with their necessary 
    files (data files, input files, reaction templates).

    Args:
        None

    Returns:
        Population, population size, mutation rate, kill rate, elitism rate,
        target element, number of cores, selection method, maximum generation
        number, target_element_lower, target element lower and upper bounds,
        target density lower and upper bounds, cell size lower and upper bounds,
        bonding probabilities lower and upper bounds, crossover type used, max
        molecular speed, a list of simulation stages (each stage is a dict
        with keys squish_factor, squish_step_size, ensemble_type, t_start,
        t_end, t_damp, run_steps; squish_factor is None for the first stage
        so that the genome condensation rate is injected at run time),
        tournament size k, ref_energy (obsolete, ignored; None unless an old
        input still has the line), and the three fitness term weights (weight_energy, weight_density,
        weight_element) that scale the energy, density, and element-percentage
        contributions to the total fitness score. All three weights default to
        1.0 when the corresponding keyword is omitted from deadmd.in.
        exp_pdf is a numpy array (r, G(r)) of the reference PDF on a 0.01
        Angstrom grid, or None when exp_pdf_file is absent or weight_pdf is 0. weight_pdf
        scales the PDF R-factor term and defaults to 0.0 so the term is
        inactive unless the user explicitly enables it in deadmd.in.
    """
    population = []
    molecules = []
    molecule_family = []
    num_molecule = []
    rxns = []
    rxn_options = []   # per reaction: react options passed to condense
    rxn_caps = []      # per reaction: (lower, upper) max_rxn range, or None
    stages = []   # list of stage dicts; see stage parsing below
    k = 2 # default value if user does not provide it or choses a selection
          # method that does not require a k value

    # Fitness term weights default to 1.0 (equal weighting). The user may
    # override any or all of them with weight_energy / weight_density /
    # weight_element in deadmd.in without changing the other terms.
    # weight_pdf defaults to 0.0: the PDF term is inactive unless the
    # user supplies both exp_pdf_file and a positive weight_pdf.
    weight_energy  = 1.0
    weight_density = 1.0
    ref_energy     = None   # obsolete, ignored
    weight_element = 1.0
    weight_pdf     = 0.0
    exp_pdf_file   = None
    # How pdf_neutron.py computes each member's G(r) for the PDF term.
    pdf_weighting  = "neutron"
    pdf_sigma      = 0.05
    # Validity gate (deadmd_validity.py): off unless "validity_gate yes".
    # "valence El n1 n2 ..." lines override the allowed bond counts;
    # "max_strain" is the largest allowed |bond strain| (fraction).
    validity_gate  = False
    valence        = {}
    max_strain     = 0.30
    # Range of gene 3, the condensation rate (squish factor: the final box
    # edge fraction of a genome stage). Used for both the initial draw and
    # mutation. Optional; the defaults are the old hardcoded init range.
    squish_lower   = 0.15
    squish_upper   = 0.33
    ga_seed        = None   # drawn below when not given
    cores_per_member = None # None: plan_cores shares the cores automatically
    lmp_mode       = "parallel"   # lmp serial|parallel: which LAMMPS binary
    # Composition term (optional): target at% per element and the
    # tolerance each deviation is measured in. Defaults: 1 at% for every
    # element except H, 2 at%; "composition_tolerance El t ..." overrides.
    target_composition    = None
    composition_tolerance = {}
    weight_composition    = 1.0
    # Energy term scale (LAMMPS energy unit per atom, kcal/mol for units
    # real): energy_term = (E/N - E_ref) / energy_tolerance.
    energy_tolerance      = 0.5
    # Mutation operator (see MUTATION_TYPE at the top of the file):
    # "mutation_type uniform|creep", "creep_sigma" (fraction of each gene's
    # range), "creep_jump" (% of mutations that are a uniform redraw).
    mutation_type         = "uniform"
    creep_sigma           = 0.1
    creep_jump            = 0.1
    # Bond-limit metric: "bond_limit El1 El2 max_pct [tolerance t]" lines
    # (see BOND_LIMITS at the top of the file) and its weight.
    bond_limits           = []
    weight_bond           = 1.0

# Read the deadmd.in input file
    with open("deadmd.in", "r") as input_file:
        input_file_lines = [line.strip() for line in input_file if line.strip()]

        for line_num, line in enumerate(input_file_lines):
            words = line.split()
            if (words[0].lower() == "num_cores"):
                num_cores = int(words[1])
            elif (words[0].lower() == "lmp"):
                lmp_mode = words[1].lower()
            elif (words[0].lower() == "cores_per_member"):
                cores_per_member = int(words[1])
            elif (words[0].lower() == "ga_seed"):
                ga_seed = int(words[1])
            elif (words[0].lower() == "gen_max"):
                gen_max = int(words[1])
            elif (words[0].lower() == "population_size"):
                population_size = int(words[1])
            elif (words[0].lower() == "elitism_rate"):
                elitism_rate = int(words[1])
            elif (words[0].lower()) == "mutation_rate":
                mutation_rate = float(words[1])/100
            elif (words[0].lower() == "mutation_type"):
                mutation_type = words[1].lower()
            elif (words[0].lower() == "creep_sigma"):
                creep_sigma = float(words[1])
            elif (words[0].lower() == "creep_jump"):
                creep_jump = float(words[1])/100
            elif (words[0].lower() == "bond_limit"):
                # bond_limit El1 El2 max_pct [tolerance t]
                try:
                    tol = 1.0
                    if len(words) == 6 and words[4].lower() == "tolerance":
                        tol = float(words[5])
                    elif len(words) != 4:
                        raise ValueError
                    bond_limits.append((words[1].lower(), words[2].lower(),
                                        float(words[3]), tol))
                except ValueError:
                    sys.exit(f"[DEAD-MD] bond_limit needs: El1 El2 max_pct "
                             f"[tolerance t], got: {line}")
            elif (words[0].lower() == "weight_bond"):
                weight_bond = float(words[1])
            elif (words[0].lower() == "kill_rate"):
                kill_rate = int(words[1])
            elif (words[0].lower() == "crossover_type"):
                crossover_type = int(words[1])
            elif (words[0].lower() == "selection"):
                selection_method = words[1]
                if selection_method == "truncation":
                    selection = 1
                elif selection_method == "rank":
                    selection = 2
                elif selection_method == "tournament":
                    selection = 3
                    k = int(words[2])
            elif (words[0].lower() == "target_element"):
                target_element = words[1].lower()
            elif (words[0].lower() == "target_element_lower"):
                target_element_lower = float(words[1])
            elif (words[0].lower() == "target_element_upper"):
                target_element_upper = float(words[1])
            elif (words[0].lower() == "target_density_lower"):
                target_density_lower = float(words[1])
            elif (words[0].lower()) == "target_density_upper":
                target_density_upper = float(words[1])
            elif (words[0].lower() == "composition_num"):
                composition_num = int(words[1])

                for molecule_line_index in range(line_num + 1, line_num + 1 + 
                        composition_num):
                    parts = input_file_lines[molecule_line_index].split()
                    molecules.append(parts[0])
                    molecule_family.append(parts[1])
                    num_molecule.append(parts[2])
                print(molecules)
                print(molecule_family)
                print(num_molecule)
            elif (words[0].lower() == "cell_size_lower"):
                cell_size_lower = float(words[1])
            elif (words[0].lower() == "cell_size_upper"):
                cell_size_upper = float(words[1])
            elif (words[0].lower() == "squish_lower"):
                squish_lower = float(words[1])
            elif (words[0].lower() == "squish_upper"):
                squish_upper = float(words[1])
            elif (words[0].lower() == "max_speed"):
                max_speed = float(words[1])
            elif (words[0].lower() == "stage"):
                # A stage line with 6 values (7 tokens including "stage")
                # is a genome stage: the condensation rate (squish_factor)
                # is drawn from the member's genome at run time, so
                # squish_factor is stored as None here. A stage line with
                # 7 values (8 tokens) is a fixed stage where the user
                # explicitly provides squish_factor as the first value --
                # useful for equilibration steps that should not compress
                # (squish_factor 1.0) or for any post-condensation stage.
                if len(words) == 7:
                    stages.append({
                        "squish_factor":    None,
                        "squish_step_size": int(words[1]),
                        "ensemble_type":    words[2],
                        "t_start":          float(words[3]),
                        "t_end":            float(words[4]),
                        "t_damp":           float(words[5]),
                        "run_steps":        int(words[6]),
                    })
                elif len(words) == 8:
                    stages.append({
                        "squish_factor":    float(words[1]),
                        "squish_step_size": int(words[2]),
                        "ensemble_type":    words[3],
                        "t_start":          float(words[4]),
                        "t_end":            float(words[5]),
                        "t_damp":           float(words[6]),
                        "run_steps":        int(words[7]),
                    })
                else:
                    raise ValueError(
                        f"stage line must have 6 or 7 values "
                        f"(got {len(words) - 1}): {line}"
                    )
            elif (words[0].lower() == "ref_energy"):
                # Obsolete: the energy term is now relative to the best
                # member of each generation. Accepted so old inputs run.
                ref_energy = float(words[1])
                print("[DEAD-MD] ref_energy is no longer used and is "
                      "ignored (energy is ranked within each generation).",
                      flush=True)
            elif (words[0].lower() == "weight_energy"):
                weight_energy = float(words[1])
            elif (words[0].lower() == "weight_density"):
                weight_density = float(words[1])
            elif (words[0].lower() == "weight_element"):
                weight_element = float(words[1])
            elif (words[0].lower() == "target_composition"):
                # target_composition El1 x1 El2 x2 ... (at%, O-free etc.)
                vals = words[1:]
                if not vals or len(vals) % 2:
                    sys.exit("[DEAD-MD] target_composition needs element / "
                             f"at% pairs, got: {line}")
                try:
                    target_composition = {vals[j].lower(): float(vals[j + 1])
                                          for j in range(0, len(vals), 2)}
                except ValueError:
                    sys.exit("[DEAD-MD] target_composition needs element / "
                             f"at% pairs, got: {line}")
            elif (words[0].lower() == "composition_tolerance"):
                vals = words[1:]
                if not vals or len(vals) % 2:
                    sys.exit("[DEAD-MD] composition_tolerance needs element "
                             f"/ at% pairs, got: {line}")
                try:
                    for j in range(0, len(vals), 2):
                        composition_tolerance[vals[j].lower()] = float(
                            vals[j + 1])
                except ValueError:
                    sys.exit("[DEAD-MD] composition_tolerance needs element "
                             f"/ at% pairs, got: {line}")
            elif (words[0].lower() == "weight_composition"):
                weight_composition = float(words[1])
            elif (words[0].lower() == "energy_tolerance"):
                energy_tolerance = float(words[1])
                if energy_tolerance <= 0.0:
                    sys.exit(f"[DEAD-MD] energy_tolerance must be positive, "
                             f"got {energy_tolerance}")
            elif (words[0].lower() == "exp_pdf_file"):
                exp_pdf_file = words[1]
            elif (words[0].lower() == "pdf_weighting"):
                if words[1].lower() not in ("neutron", "none"):
                    sys.exit(f"[DEAD-MD] pdf_weighting must be neutron or "
                             f"none, got {words[1]}")
                pdf_weighting = words[1].lower()
            elif (words[0].lower() == "pdf_sigma"):
                pdf_sigma = float(words[1])
            elif (words[0].lower() == "validity_gate"):
                if words[1].lower() not in ("yes", "no"):
                    sys.exit(f"[DEAD-MD] validity_gate must be yes or no, "
                             f"got {words[1]}")
                validity_gate = words[1].lower() == "yes"
            elif (words[0].lower() == "valence"):
                valence[words[1].capitalize()] = [int(x) for x in words[2:]]
            elif (words[0].lower() == "max_strain"):
                max_strain = float(words[1])
            elif (words[0].lower() == "weight_pdf"):
                weight_pdf = float(words[1])
            elif (words[0].lower() == "reactions_num"):
                reactions_num = int(words[1])

                for rxn_line_index in range(line_num + 1, line_num +
                        1 + reactions_num):
                    words = input_file_lines[rxn_line_index].split()
                    # Trailing react options (intra/inter, rmin R) are not
                    # genes: split them off and pass them to condense as is.
                    k = 4
                    while k < len(words) and re.match(r'^[0-9.eE+-]+$',
                                                      words[k]):
                        k += 1
                    rxns.append(words[:k])
                    rxn_options.append(words[k:])
                # One Rmax range and one probability range per reaction,
                # from that reaction's own line. Gene 4 + r is reaction r's
                # Rmax, gene 4 + reactions_num + r its probability.
                for rxn, opts in zip(rxns, rxn_options):
                    if len(rxn) not in (6, 8):
                        sys.exit("[DEAD-MD] each reaction line needs 6 or 8 "
                                 "values (mol1 site1 mol2 site2 prob_lower "
                                 "prob_upper [rmax_lower rmax_upper]), then "
                                 "optionally intra|inter and rmin R, got: "
                                 f"{' '.join(rxn + opts)}")
                    o = [w.lower() for w in opts]
                    if not {"inter", "intra"} & set(o):
                        sys.exit("[DEAD-MD] each reaction line must end "
                                 "with inter (two molecules) or intra "
                                 "(within one), got: "
                                 f"{' '.join(rxn + opts)}")
                    kept, cap = [], None
                    while o:
                        w = o.pop(0)
                        if w in ("intra", "inter"):
                            kept.append(w)
                            continue
                        if w == "rmin" and o and re.match(
                                r'^[0-9.eE+-]+$', o[0]):
                            kept += [w, o.pop(0)]
                            continue
                        if w == "max_rxn" and o and o[0].isdigit():
                            # max_rxn LO HI (a gene) or max_rxn N (fixed).
                            lo = int(o.pop(0))
                            hi = int(o.pop(0)) if o and o[0].isdigit() else lo
                            if not 1 <= lo <= hi:
                                sys.exit("[DEAD-MD] max_rxn range must satisfy "
                                         "1 <= lower <= upper (whole numbers), "
                                         f"got: {' '.join(rxn + opts)}")
                            cap = (lo, hi)
                            continue
                        sys.exit(f"[DEAD-MD] unknown reaction option '{w}' "
                                 "(expected intra, inter, rmin R or max_rxn "
                                 f"LO HI) in: {' '.join(rxn + opts)}")
                    # max_rxn is a gene: its range is not passed through; the
                    # member's value is written to condense.in instead.
                    opts[:] = kept
                    rxn_caps.append(cap)
                bonding_probability_lower = [float(rxn[4]) for rxn in rxns]
                bonding_probability_upper = [float(rxn[5]) for rxn in rxns]
                # Without an Rmax range the reaction keeps condense's 5.0.
                rmax_lower = [float(rxn[6]) if len(rxn) == 8 else 5.0
                              for rxn in rxns]
                rmax_upper = [float(rxn[7]) if len(rxn) == 8 else 5.0
                              for rxn in rxns]
                for lo, hi, rxn in zip(rmax_lower, rmax_upper, rxns):
                    if not 0.0 < lo <= hi:
                        sys.exit("[DEAD-MD] reaction rmax range must satisfy "
                                 "0 < lower <= upper, got: "
                                 f"{' '.join(rxn)}")
                for lo, hi, rxn in zip(bonding_probability_lower,
                                       bonding_probability_upper, rxns):
                    if not 0.0 <= lo <= hi <= 1.0:
                        sys.exit("[DEAD-MD] reaction probability range must "
                                 "satisfy 0 <= lower <= upper <= 1, got: "
                                 f"{' '.join(rxn)}")

        # Check the input now, before any member is simulated, so that a
        # typo stops the run in seconds instead of after generation 1.
        keywords_given = {line.split()[0].lower() for line in input_file_lines}
        required = ["num_cores", "gen_max", "population_size", "elitism_rate",
                    "mutation_rate", "kill_rate", "crossover_type",
                    "selection", "target_density_lower",
                    "target_density_upper", "composition_num",
                    "cell_size_lower", "cell_size_upper", "max_speed",
                    "stage", "reactions_num"]
        # The single-element term is required unless target_composition
        # replaces it.
        element_kws = ["target_element", "target_element_lower",
                       "target_element_upper"]
        if target_composition is None:
            required += element_kws
        missing = [kw for kw in required if kw not in keywords_given]
        if missing:
            sys.exit(f"[DEAD-MD] deadmd.in is missing required keyword(s): "
                     f"{', '.join(missing)}")
        if target_composition is not None:
            for el, x in target_composition.items():
                if not 0.0 <= x <= 100.0:
                    sys.exit(f"[DEAD-MD] target_composition: {el} {x} is "
                             f"not an at% between 0 and 100")
            total = sum(target_composition.values())
            if total > 100.5:
                sys.exit(f"[DEAD-MD] target_composition sums to {total:.2f} "
                         f"at%, more than 100")
            extra = set(composition_tolerance) - set(target_composition)
            if extra:
                sys.exit(f"[DEAD-MD] composition_tolerance lists element(s) "
                         f"not in target_composition: {', '.join(sorted(extra))}")
            for el, t in composition_tolerance.items():
                if t <= 0.0:
                    sys.exit(f"[DEAD-MD] composition_tolerance for {el} must "
                             f"be positive, got {t}")
            tolerance = {el: composition_tolerance.get(
                             el, 2.0 if el == "h" else 1.0)
                         for el in target_composition}
            print("[DEAD-MD] composition term: target "
                  + " ".join(f"{el} {x}" for el, x in target_composition.items())
                  + " at%, tolerance "
                  + " ".join(f"{el} {t}" for el, t in tolerance.items())
                  + f", weight {weight_composition}", flush=True)
            if not set(element_kws) <= keywords_given:
                # Gene 0 still exists: fix it at the first element's target
                # and switch the single-element term off.
                if set(element_kws) & keywords_given:
                    sys.exit("[DEAD-MD] give all three target_element "
                             "keywords or none (target_composition replaces "
                             "them)")
                target_element = next(iter(target_composition))
                target_element_lower = target_element_upper = \
                    target_composition[target_element]
                weight_element = 0.0
            elif weight_element > 0.0:
                print("[DEAD-MD] note: both target_element and "
                      "target_composition are active (set weight_element 0 "
                      "to use the composition term alone)", flush=True)
            composition_options = (target_composition, tolerance,
                                   weight_composition)
        else:
            composition_options = None
        if mutation_type not in ("uniform", "creep"):
            sys.exit(f"[DEAD-MD] mutation_type must be uniform or creep, "
                     f"got {mutation_type}")
        if not 0.0 < creep_sigma <= 1.0:
            sys.exit(f"[DEAD-MD] creep_sigma is a fraction of the gene range "
                     f"(0 < creep_sigma <= 1), got {creep_sigma}")
        if not 0.0 <= creep_jump <= 1.0:
            sys.exit(f"[DEAD-MD] creep_jump is a percent (0-100), got "
                     f"{creep_jump * 100}")
        for el1, el2, mx, tol in bond_limits:
            if el1 not in MASS_TO_ELEMENT.values() or \
               el2 not in MASS_TO_ELEMENT.values():
                sys.exit(f"[DEAD-MD] bond_limit: unknown element in "
                         f"{el1} {el2}")
            if not 0.0 <= mx <= 100.0 or tol <= 0.0:
                sys.exit(f"[DEAD-MD] bond_limit {el1} {el2}: max_pct must be "
                         f"0-100 and tolerance > 0, got {mx}, {tol}")
        global BOND_LIMITS, WEIGHT_BOND
        BOND_LIMITS, WEIGHT_BOND = bond_limits, weight_bond
        if bond_limits:
            print("[DEAD-MD] bond metric: "
                  + "; ".join(f"% {a} bonded to {b} <= {mx} (tolerance {t})"
                              for a, b, mx, t in bond_limits)
                  + f", weight {weight_bond}", flush=True)
        global MUTATION_TYPE, CREEP_SIGMA, CREEP_JUMP
        MUTATION_TYPE, CREEP_SIGMA, CREEP_JUMP = (mutation_type, creep_sigma,
                                                  creep_jump)
        if mutation_type == "creep":
            print(f"[DEAD-MD] mutation: creep, sigma {creep_sigma} x gene "
                  f"range, {creep_jump * 100:g}% of mutations are a uniform "
                  f"redraw", flush=True)
        else:
            print("[DEAD-MD] mutation: uniform redraw over the gene range",
                  flush=True)
        if crossover_type not in (1, 2):
            sys.exit(f"[DEAD-MD] crossover_type must be 1 (one point) or "
                     f"2 (two points), got {crossover_type}")
        if selection_method not in ("truncation", "rank", "tournament"):
            sys.exit(f"[DEAD-MD] selection must be truncation, rank or "
                     f"tournament, got {selection_method}")
        if selection == 3 and not 1 <= k <= population_size:
            sys.exit(f"[DEAD-MD] tournament size k must be between 1 and "
                     f"population_size ({population_size}), got {k}")
        if lmp_mode not in ("serial", "parallel"):
            sys.exit(f"[DEAD-MD] lmp must be serial or parallel, got "
                     f"{lmp_mode}")
        if lmp_mode == "serial":
            # The serial binary runs one rank: one core per member.
            if cores_per_member not in (None, 1):
                sys.exit(f"[DEAD-MD] lmp serial runs each member on 1 core; "
                         f"cores_per_member {cores_per_member} needs lmp "
                         f"parallel")
            cores_per_member = 1
        # Find the executables now, not after the first member fails.
        lmp_bin = "lmp" if lmp_mode == "serial" else "lmp_parallel"
        for exe in (lmp_bin, "srun"):
            if shutil.which(exe) is None:
                sys.exit(f"[DEAD-MD] cannot find the executable '{exe}' in "
                         f"PATH (needed for lmp {lmp_mode}); install it or "
                         f"add its directory to PATH")
        if cores_per_member is not None and not 1 <= cores_per_member <= num_cores:
            sys.exit(f"[DEAD-MD] cores_per_member must be between 1 and "
                     f"num_cores ({num_cores}), got {cores_per_member}")
        slurm_tasks = os.environ.get("SLURM_NTASKS")
        if slurm_tasks and num_cores > int(slurm_tasks):
            sys.exit(f"[DEAD-MD] num_cores {num_cores} is more than the "
                     f"job's slurm tasks ({slurm_tasks}, #SBATCH -n)")
        if not 0.0 < squish_lower <= squish_upper:
            sys.exit(f"[DEAD-MD] squish range must satisfy 0 < squish_lower "
                     f"<= squish_upper, got {squish_lower} {squish_upper}")


        # Seed the GA before its first random draw (the initial genomes).
        global GA_SEED
        if ga_seed is None:
            ga_seed = random.SystemRandom().randint(1, 2**31 - 1)
        GA_SEED = ga_seed
        random.seed(ga_seed)
        np.random.seed(ga_seed % 2**32)
        print(f"[DEAD-MD] ga_seed {ga_seed}", flush=True)

        # Filling the population array. The structure is a 2d array. 
        # Each row is a member (individual or genome).
        # Each column = [target_element_percentage, target_density, 
        #  cell_size, condensation_rate, rmax_1, ..., rmax_n,
        #  bonding_probability_1, ..., bonding_probability_n]
        # There is one Rmax and one bonding probability per reaction, so the
        #  array is (population_size, 4 + 2 * reactions_num + number of max_rxn genes).
        # The sequence of probabilities in the array is in respect to the sequence
        #  of the reaction you list in the deadmd.in or lamps.in
        # Reactions with "max_rxn LO HI" add one whole-number gene each at
        #  the END (gene 4 + 2 * reactions_num + j), so the other genes keep
        #  their positions.
        global MAXRXN_GENES
        MAXRXN_GENES = [(r, cap[0], cap[1])
                        for r, cap in enumerate(rxn_caps) if cap is not None]
        if MAXRXN_GENES:
            print("[DEAD-MD] max_rxn genes: "
                  + "; ".join(f"reaction {r + 1} ({' '.join(rxns[r][:4])}) "
                              f"{lo}-{hi}" for r, lo, hi in MAXRXN_GENES),
                  flush=True)
        population = np.empty((population_size,
                               4 + 2 * reactions_num + len(MAXRXN_GENES)))
        for i in range(population_size):
            population[i][0] = random.uniform(target_element_lower, target_element_upper)
            population[i][1] = random.uniform(target_density_lower, target_density_upper)
            population[i][2] = round(random.uniform(cell_size_lower, cell_size_upper), 2)
            population[i][3] = round(random.uniform(squish_lower, squish_upper), 2)
            for r in range(reactions_num):
                population[i][4 + r] = round(
                    random.uniform(rmax_lower[r], rmax_upper[r]), 2)
                population[i][4 + reactions_num + r] = random.uniform(
                    bonding_probability_lower[r], bonding_probability_upper[r])
        # Cap genes are drawn after every other gene of every member, so a
        # run with max_rxn and the same run without it (same ga_seed) share
        # genes 0 .. 3 + 2n in generation 1.
        for i in range(population_size):
            for j, (r, lo, hi) in enumerate(MAXRXN_GENES):
                population[i][4 + 2 * reactions_num + j] = random.randint(lo, hi)
        print(population)

        # Load the experimental PDF if the user supplied a file and a
        # positive weight. The file is two-column whitespace-delimited:
        # column 0 is r (Angstrom), column 1 is G(r), on the 0.01 Angstrom
        # grid pdf_neutron.py writes.
        exp_pdf = None
        if exp_pdf_file is not None and weight_pdf > 0.0:
            # Store both columns so fitness_function can align the
            # experimental r-grid against the simulated r-grid at
            # run time. Column 0 is r (Angstrom), column 1 is G(r).
            exp_pdf = np.loadtxt(exp_pdf_file)
            # fitness_function aligns the two grids by slicing, so a
            # repeated or backwards r value, or a file that stops short
            # of the 10 Angstrom pdf_neutron.py can reach, would silently fail
            # every member.
            if np.any(np.diff(exp_pdf[:, 0]) <= 0):
                sys.exit(f"[DEAD-MD] {exp_pdf_file}: the r column must be "
                         f"strictly increasing (repeated or backwards rows)")
            if exp_pdf[-1, 0] < 10.0 - 1e-6:
                sys.exit(f"[DEAD-MD] {exp_pdf_file}: r must reach 10.00 "
                         f"Angstrom, file stops at {exp_pdf[-1, 0]}")
            print(f"[DEAD-MD] Experimental PDF loaded: "
                  f"{exp_pdf.shape[0]} points")

    return (
                population,
                population_size,
                mutation_rate,
                kill_rate,
                elitism_rate, 
                target_element,
                num_cores, 
                selection, 
                gen_max,
                target_element_lower,
                target_element_upper,
                target_density_lower,
                target_density_upper,
                cell_size_lower, 
                cell_size_upper,
                squish_lower,
                squish_upper,
                rmax_lower,
                rmax_upper,
                bonding_probability_lower, 
                bonding_probability_upper,
                crossover_type,
                max_speed,
                stages,
                k,
                molecules,
                num_molecule,
                composition_num,
                reactions_num,
                rxns,
                rxn_options,
                ref_energy,
                weight_energy,
                weight_density,
                weight_element,
                exp_pdf,
                weight_pdf,
                (pdf_weighting, pdf_sigma),
                (valence, max_strain) if validity_gate else None,
                cores_per_member,
                lmp_mode,
                composition_options,
                energy_tolerance
            )


def get_element_percentage(skl_file, target_element):
    """
    Obtain the atomic percentage of a given target element by calling the external
    'element_pct.py' script as a subprocess.
    """
    element_percentage = subprocess.run(["element_pct.py", skl_file, target_element],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, check=True)
    return float(element_percentage.stdout.strip())


def get_composition(skl_file):
    """
    Return {element: at%} for every element in an imago .skl file (the
    atom lines of its "fractional N" section, as element_pct.py reads them).
    """
    counts = {}
    total = 0
    reading = False
    with open(skl_file) as f:
        for line in f:
            s = line.strip().lower()
            if s.startswith("fractional"):
                parts = s.split()
                if len(parts) == 2 and parts[1].isdigit():
                    total = int(parts[1])
                    reading = True
                continue
            if reading and (s == "" or s.startswith("space")):
                break
            if reading:
                parts = s.split()
                if len(parts) >= 4:
                    el = re.sub(r"\d+$", "", parts[0])
                    counts[el] = counts.get(el, 0) + 1
    if total == 0:
        raise ValueError(f"{skl_file}: atom count not found or zero")
    return {el: 100.0 * n / total for el, n in counts.items()}


def read_last_value(filename):
    """
    Read the last non-empty line from a file and return it as a float.
    """
    with open(filename, "r") as f_values:
        lines = [line.strip() for line in f_values if line.strip()]
        return float (lines[-1])


def plan_cores(n_run, num_cores, cores_per_member=None):
    """
    Share num_cores among the n_run members simulated this generation.

    Returns (cores, slots): cores lists the MPI rank count of each member,
    in run order; slots is how many members run at the same time. A member
    starts as soon as a slot is free.

    cores_per_member given (keyword in deadmd.in): every member gets that
    many cores and num_cores // cores_per_member run at once.
    Automatic (default):
      n_run >= num_cores: 1 core each, num_cores at once. Many independent
        runs give more throughput than fewer MPI runs, because MPI scaling
        on a few hundred to a few thousand atoms is far from linear.
      n_run <  num_cores: all members at once, the cores shared out as
        evenly as possible (10 cores, 4 members -> 3 3 2 2), so no core
        sits idle.
    """
    if n_run == 0:
        return [], 0
    if cores_per_member:
        return ([cores_per_member] * n_run,
                min(n_run, max(1, num_cores // cores_per_member)))
    if n_run >= num_cores:
        return [1] * n_run, num_cores
    base, extra = divmod(num_cores, n_run)
    return [base + 1 if i < extra else base for i in range(n_run)], n_run


def srun_mem_option():
    """
    Memory share for one srun step. Without it each step claims the job's
    whole memory and Slurm runs the steps one after another (tested
    2026-09-25, jobs/features/concurrent/srun_test).
    """
    per_cpu = os.environ.get("SLURM_MEM_PER_CPU")
    if per_cpu:
        return f"--mem-per-cpu={per_cpu}M"
    per_node = os.environ.get("SLURM_MEM_PER_NODE")
    ntasks = os.environ.get("SLURM_NTASKS")
    if per_node and ntasks:
        nnodes = int(os.environ.get("SLURM_NNODES", "1"))
        return f"--mem-per-cpu={int(per_node) * nnodes // int(ntasks)}M"
    return ""


def run_member(member_dir, member, cores, gen_count, target_element,
               exp_pdf, pdf_options, gate_options, mem_opt,
               lmp_bin="lmp_parallel"):
    """
    Build and run one member: condense.py, LAMMPS on `cores` MPI
    ranks, then the per-member analysis. Several of these run at once in
    threads, so nothing here changes the working directory: every command
    gets cwd=, and all their screen output goes to member_dir/member.out.
    """
    lammps_dir = os.path.join(member_dir, "lammps")
    print(f"[DEAD-MD] Gen {gen_count} Member {member}: starting on "
          f"{cores} core(s)", file=sys.stderr, flush=True)
    with open(os.path.join(member_dir, "member.out"), "a") as out:

        def run(cmd, cwd):
            subprocess.run(cmd, shell=True, check=True, cwd=cwd,
                           stdout=out, stderr=subprocess.STDOUT)

        # Wrap each individual in a try/except so that one failed member
        # (e.g. condense dying on a missing angle type, or LAMMPS crashing
        # on a neighbor list overflow) does not kill the entire EA run.
        # Instead, the member is assigned worst-case fitness values (1e99)
        # so the EA treats it as a bad individual and discards it
        # naturally via selection/kill_rate.
        try:
            # No -s: the member seed is the "seed" line of its condense.in
            # (written below with the rest of the member's input).
            run("condense.py", member_dir)
            # lmp_bin is lmp_parallel (MPI build, openmpi/gcc libs baked in
            # as RPATH, so no module load is needed) or the serial lmp
            # (always 1 core). Both go through srun: --exact and a memory
            # share let several members' steps share the job's cores at
            # once, and Slurm keeps each on its own core(s).
            run(f"srun --exact -n {cores} {mem_opt} {lmp_bin} "
                f"-in lammps.in", lammps_dir)
            run("dump2skl.py -d dump.coarse -a lammps.dat -f -1", lammps_dir)
            # Calculate elemental percentage.
            pct = get_element_percentage(
                os.path.join(lammps_dir, "imago.skl"), target_element)
            with open(os.path.join(lammps_dir, "element_evolve"), "a") as f:
                f.write(f"{pct}\n")
            # Every element's at%, one "El x El x ..." line per run; the
            # composition term of fitness_function reads the last line.
            comp = get_composition(os.path.join(lammps_dir, "imago.skl"))
            with open(os.path.join(lammps_dir, "composition_evolve"),
                      "a") as f:
                f.write(" ".join(f"{el} {x:.4f}" for el, x in comp.items())
                        + "\n")
            # Compute the member's normalised G(r) when a reference PDF has
            # been loaded. pdf_neutron.py reads imago.skl and writes
            # pdf_neutron.plot inside lammps/, out to half the box width.
            # fitness_function reads pdf_neutron.plot later.
            if exp_pdf is not None:
                weighting, sigma = pdf_options
                run(f"pdf_neutron.py -weighting {weighting} -sigma {sigma}",
                    lammps_dir)
            # Final atom count, printed by the lammps.in tail that
            # condense.py writes. The fitness divides the total
            # energy by it, since each reaction removes atoms.
            # A missing count is written as the 1e99 sentinel, which
            # fitness_function treats like a crash.
            with open(os.path.join(lammps_dir, "log.lammps")) as f_log:
                counts = re.findall(r"Number of atoms is (\d+);",
                                    f_log.read())
            natoms = counts[-1] if counts else "1e99"
            with open(os.path.join(lammps_dir, "natoms_evolve"), "a") as f:
                f.write(f"{natoms}\n")
            # Validity gate: an invalid model (wrong coordination, broken
            # bookkeeping, absurd bond strain) is not scored at all.
            # deadmd_validity.py writes its reasons to validity.txt;
            # validity_evolve gets 1 (valid) or 0, and fitness_function
            # gives a 0 the crash sentinel.
            if gate_options is not None:
                gate_valence, gate_strain = gate_options
                cmd = f"deadmd_validity.py -max_strain {gate_strain}"
                for el, counts_allowed in gate_valence.items():
                    cmd += (f" -valence {el} "
                            + " ".join(map(str, counts_allowed)))
                valid = subprocess.run(cmd, shell=True, cwd=lammps_dir,
                                       stdout=subprocess.DEVNULL
                                       ).returncode == 0
                with open(os.path.join(lammps_dir, "validity_evolve"),
                          "a") as f:
                    f.write(f"{int(valid)}\n")
                if not valid:
                    print(f"[DEAD-MD] Gen {gen_count} Member {member}: "
                          f"INVALID (see lammps/validity.txt), assigning "
                          f"worst fitness.", file=sys.stderr, flush=True)

        except subprocess.CalledProcessError as e:
            print(f"[DEAD-MD] Gen {gen_count} Member {member}: CRASHED"
                  f" (command: {e.cmd}; output in {member_dir}/member.out),"
                  f" assigning worst fitness.", file=sys.stderr, flush=True)
            # fitness_function reads totE_evolve, density_evolve,
            # element_evolve and natoms_evolve from lammps/. Write
            # sentinels so the EA can continue rather than crashing on
            # missing files.
            os.makedirs(lammps_dir, exist_ok=True)
            for fname in ["totE_evolve", "density_evolve", "element_evolve",
                          "natoms_evolve"]:
                with open(os.path.join(lammps_dir, fname), "a") as f:
                    f.write("1e99\n")
            return
    print(f"[DEAD-MD] Gen {gen_count} Member {member}: finished",
          file=sys.stderr, flush=True)


def run_lammps_simulations(population, population_size, num_cores, target_element,
                           molecules, num_molecule, composition_num, reactions_num,
                           rxns, rxn_options, max_speed, stages, gen_count,
                           exp_pdf=None,
                           pdf_options=("neutron", 0.05), gate_options=None,
                           elite_sources=None, cores_per_member=None,
                           lmp_mode="parallel"):
    """
    This function will create the necessary input files for a lammps condensation
    and run lammps for each member of the Genetic Algorithm (GA)

    Members run concurrently: plan_cores shares num_cores among the members
    to simulate, and each one runs in its own thread (run_member) with its
    own srun step. condense.in files and elite copies are written first,
    in this thread.

    elite_sources maps a member number (1-based) to the absolute path of the
    previous-generation member it is an elite copy of. Elites are NOT rerun:
    that member's directory is copied (without the large dump.* files) and
    a note ELITE_COPY_OF is left pointing at the original, so fitness_function
    and the recorders read the elite's raw results (energy, atom count,
    density, composition, PDF) and re-score it with its new generation.
    """
    if elite_sources is None:
        elite_sources = {}
    # Creating seperate directories for generation and each member, each contains 
    # the lamps.in file. Lammps for each member is ran from the corresponding directory
    gen_dir = os.path.abspath(f"generation_{gen_count}")
    os.makedirs(gen_dir, exist_ok=True)
    to_run = []
    for member in range(1, population_size + 1):
        member_dir = os.path.join(gen_dir, str(member))
        if member in elite_sources:
            src = elite_sources[member]
            # An elite kept for several generations: point at the member
            # that was actually simulated (where the dumps are), not at
            # the previous copy.
            origin = src
            if os.path.exists(f"{src}/ELITE_COPY_OF"):
                with open(f"{src}/ELITE_COPY_OF") as f_prev:
                    origin = f_prev.readline().strip()
            shutil.copytree(src, member_dir, dirs_exist_ok=True,
                            ignore=shutil.ignore_patterns("dump.*"))
            with open(f"{member_dir}/ELITE_COPY_OF", "w") as f_note:
                f_note.write(f"{origin}\n(trajectory dumps are there; this "
                             f"member was not rerun)\n")
            print(f"[DEAD-MD] Gen {gen_count} Member {member}: elite, "
                  f"copied from {src}, originally {origin} (not rerun)",
                  file=sys.stderr, flush=True)
            continue
        os.makedirs(member_dir, exist_ok=True)
        with open(os.path.join(member_dir, "condense.in"), "w") as lamps_input_file:
            lamps_input_file.write(f"composition {composition_num}\n")
            for mol, num in zip(molecules, num_molecule):
                lamps_input_file.write(f"{mol} family1 {num}\n")
             
            lamps_input_file.write(f"\ncell_size {population[member - 1][2]}\n\n")
            lamps_input_file.write(f"max_speed {max_speed}\n")
            lamps_input_file.write(
                f"seed {member_seed(gen_count, member)}\n")
            # Write all simulation stages. For genome stages (squish_factor
            # is None) the condensation rate is drawn from the member's
            # genome so it varies across members. Fixed stages use the
            # squish_factor the user specified in deadmd.in directly.
            for stage in stages:
                sf = (population[member - 1][3]
                      if stage["squish_factor"] is None
                      else stage["squish_factor"])
                lamps_input_file.write(
                    f"stage {sf} {stage['squish_step_size']} "
                    f"{stage['ensemble_type']} {stage['t_start']} "
                    f"{stage['t_end']} {stage['t_damp']} "
                    f"{stage['run_steps']}\n"
                )
            lamps_input_file.write("\n")
            lamps_input_file.write(f"target_density {population[member - 1][1]}\n\n")
            
            lamps_input_file.write(f"reactions {reactions_num}\n")

            caps = {r: int(round(population[member - 1][4 + 2 * reactions_num + j]))
                    for j, (r, lo, hi) in enumerate(MAXRXN_GENES)}
            for rxn_index in range(reactions_num):
                lamps_input_file.write(f"{rxns[rxn_index][0]} {rxns[rxn_index][1]} "
                                       f"{rxns[rxn_index][2]} {rxns[rxn_index][3]} "
                                       f"{population[member - 1][4 + reactions_num + rxn_index]} "
                                       f"{population[member - 1][4 + rxn_index]}"
                                       f"{''.join(' ' + w for w in rxn_options[rxn_index])}"
                                       + (f" max_rxn {caps[rxn_index]}"
                                          if rxn_index in caps else "")
                                       + "\n")

        to_run.append(member)

    cores, slots = plan_cores(len(to_run), num_cores, cores_per_member)
    lmp_bin = "lmp" if lmp_mode == "serial" else "lmp_parallel"
    if to_run:
        per = (f"{cores[0]} core(s) each" if len(set(cores)) == 1
               else "cores " + " ".join(map(str, cores)))
        note = ""
        if slots < len(to_run) and len(to_run) % slots:
            lo = len(to_run) - len(to_run) % slots
            note = (f"; the last round uses {len(to_run) % slots} of {slots}"
                    f" slots -- {lo} or {lo + slots} members to simulate"
                    f" would leave none idle")
        print(f"[DEAD-MD] Gen {gen_count}: {len(to_run)} member(s) to "
              f"simulate on {num_cores} core(s) with {lmp_bin}: {per}, "
              f"{slots} at a time{note}", flush=True)
    mem_opt = srun_mem_option()
    with ThreadPoolExecutor(max_workers=max(slots, 1)) as pool:
        jobs = [pool.submit(run_member, os.path.join(gen_dir, str(member)),
                            member, k, gen_count, target_element, exp_pdf,
                            pdf_options, gate_options, mem_opt, lmp_bin)
                for member, k in zip(to_run, cores)]
        for job in jobs:
            job.result()


def _pdf_rfactor(G_sim, exp_pdf):
    """
    Compute the scale-optimised R-factor from two G(r) arrays.

    The optimal scalar s is determined analytically by minimising the
    sum of squared residuals between s*G_sim and G_exp:

        s = (G_sim . G_exp) / |G_sim|^2

    This makes the comparison insensitive to the overall amplitude (e.g. the
    choice of neutron normalisation) -- only peak positions, relative
    heights and signs matter.

    The R-factor is then:

        R = sqrt( sum((s*G_sim - G_exp)^2) / sum(G_exp^2) )

    Args:
        G_sim (np.ndarray): Simulated G(r) values on a common r-grid.
        exp_pdf (np.ndarray): Experimental G(r) values on the same
            r-grid as G_sim.

    Returns:
        float: The scale-optimised R-factor (lower is better).

    Raises:
        ValueError: If G_sim is all-zero (degenerate or crashed
            simulation).
    """
    if np.dot(G_sim, G_sim) < 1e-12:
        raise ValueError("all-zero G(r) -- degenerate or crashed simulation")
    s = np.dot(G_sim, exp_pdf) / np.dot(G_sim, G_sim)
    return float(np.sqrt(
        np.sum((s * G_sim - exp_pdf) ** 2)
        / (np.sum(exp_pdf ** 2) + 1e-12)
    ))


def compute_pdf_rfactor(sim_pdf_file, exp_pdf_file):
    """
    Compute the scale-optimised R-factor between a simulated G(r) file
    and an experimental G(r) file.

    Both files must have two columns: r (Angstrom) and G(r), with a
    0.01 Angstrom increment. The experimental file may extend beyond
    10 Angstrom; it is truncated to r <= 10 Angstrom to match the
    range of the simulated pdf_neutron.plot output.

    This is the public entry point intended for use outside deadmd.py.
    Internally, fitness_function calls _pdf_rfactor directly to avoid
    reloading the experimental file on every member every generation.

    Args:
        sim_pdf_file (str): Path to the simulated pdf_neutron.plot file
            produced by pdf_neutron.py.
        exp_pdf_file (str): Path to the experimental G(r) data file
            (two-column: r, G(r)).

    Returns:
        float: The scale-optimised R-factor (lower is better).

    Raises:
        FileNotFoundError: If either file does not exist.
        ValueError: If the simulated G(r) is all-zero.
        OSError: If either file cannot be read.
    """
    sim_data = np.loadtxt(sim_pdf_file)
    G_sim    = sim_data[:, 1]
    exp_data = np.loadtxt(exp_pdf_file)
    r_exp    = exp_data[:, 0]
    g_exp    = exp_data[:, 1]
    # Align by the start of the simulated r-grid. Find the first
    # experimental point at r >= r_sim[0], then take exactly
    # len(G_sim) consecutive points forward. Points before r_sim[0]
    # or after r_sim[-1] in the experimental file are silently ignored.
    r_sim     = sim_data[:, 0]
    start_idx = int(np.searchsorted(r_exp, r_sim[0] - 1e-9))
    exp_pdf   = g_exp[start_idx : start_idx + len(G_sim)]
    if len(exp_pdf) != len(G_sim):
        raise ValueError(
            f"Experimental G(r) does not have enough points from"
            f" r={r_sim[0]:.4f} to cover the simulated range"
            f" ({len(G_sim)} points needed, {len(exp_pdf)} available)."
            f" Ensure both files use the same r-increment."
        )
    return _pdf_rfactor(G_sim, exp_pdf)


def best_energy_per_atom(population_size, gen_count):
    """
    Lowest energy per atom (total energy / final atom count) among the
    members of generation gen_count that did not crash and were not failed
    by the validity gate. None if no member qualifies. main() takes the
    first generation's value as the fixed energy reference of the run.
    """
    SENTINEL = 1e99
    e_per_atom = []
    for i in range(population_size):
        d = f"generation_{gen_count}/{i + 1}/lammps"
        e, n = (read_last_value(f"{d}/totE_evolve"),
                read_last_value(f"{d}/natoms_evolve"))
        if (os.path.exists(f"{d}/validity_evolve")
                and read_last_value(f"{d}/validity_evolve") == 0):
            n = SENTINEL
        if e < SENTINEL and n < SENTINEL:
            e_per_atom.append(e / n)
    return min(e_per_atom) if e_per_atom else None


def fitness_function(population, population_size, gen_count, energy_ref,
                     weight_energy=1.0, weight_density=1.0,
                     weight_element=1.0, exp_pdf=None, weight_pdf=0.0,
                     composition_options=None, energy_tolerance=0.5):
    """
    Compute a generation-independent fitness score for each individual.

    Three or four terms are summed (lower total is better), depending on
    whether an experimental PDF has been supplied:

      1. Energy term -- energy per atom (total energy E over the final
         atom count N) relative to energy_ref, the LOWEST energy per atom
         of generation 1's scored members (crashed / invalid excluded),
         fixed for the whole run:
             energy_term = (E/N - energy_ref) / energy_tolerance
         Generation 1's best scores 0; one energy_tolerance higher per atom
         scores +1; a later member that beats it scores negative, in
         proportion (no absolute value, so lower energy is rewarded). Nothing has
         to be calibrated per material, and because the reference is fixed
         the fitness values ARE comparable across generations.

      2. Density term -- relative absolute deviation of the simulated
         density from each member's own target density (population[i][1]),
         which was the value given to LAMMPS as input:
             density_term = |actual - target| / target
         abs IS used here because density has a specific target value:
         overshooting and undershooting are equally wrong. Using the
         member's own target keeps this term generation-independent.

      3. Element percentage term -- same form as the density term, using
         the member's own target element percentage (population[i][0]):
             element_term = |actual - target| / target
         Same reasoning: there is a specific desired percentage, so
         both directions of deviation are penalized symmetrically.

      4. PDF R-factor term (optional) -- scale-optimised R-factor between
         the member's G(r) from pdf_neutron.py (to half its box width) and
         the reference G(r) over the same r range. The optimal scale factor s is found analytically:
             s        = (G_sim . G_exp) / |G_sim|^2
             pdf_term = sqrt( sum((s*G_sim - G_exp)^2) / sum(G_exp^2) )
         Using an optimal scale makes the comparison insensitive to the
         overall amplitude; only peak
         positions and relative heights drive the score. Active only when
         exp_pdf is not None and weight_pdf > 0.

      5. Composition term (optional, target_composition in deadmd.in) --
         for each target element, the member's at% x (from imago.skl)
         against the target x_target, in units of that element's
         tolerance sigma (default 1 at%, H 2 at%), squared and averaged:
             composition_term = (1/n) * sum over El of |x - x_target| / sigma
         n = number of listed elements: 1 means the elements are off by one
         tolerance on average.
         Active when composition_options is not None; weight
         weight_composition.

    Crashed members (sentinel value 1e99) bypass all computation and
    receive a sentinel fitness so selection discards them naturally.

    The three terms are combined as a weighted sum:
        fitness = weight_energy  * energy_term
                + weight_density * density_term
                + weight_element * element_term
    All weights default to 1.0 (equal emphasis). Raising a weight makes
    the GA optimise that quantity more aggressively; setting a weight to
    0.0 removes that term from the selection pressure entirely.

    Args:
        population:      numpy array of genomes, shape
                         (population_size, 4 + 2 * reactions_num + number of max_rxn genes)
        population_size: number of individuals in the current generation
        gen_count:       current generation index
        energy_ref:      fixed reference energy per atom (generation 1's
                         lowest; see best_energy_per_atom)
        weight_energy:   multiplicative weight for the energy term
                         (default 1.0; keyword weight_energy in deadmd.in).
        weight_density:  multiplicative weight for the density term
                         (default 1.0; keyword weight_density in deadmd.in).
        weight_element:  multiplicative weight for the element-percentage term
                         (default 1.0; keyword weight_element in deadmd.in).
        exp_pdf:         numpy array of experimental G(r), shape (1000,),
                         covering 0.01-10.00 Angstrom at 0.01 Angstrom
                         spacing. None when the PDF term is not requested.
        weight_pdf:      multiplicative weight for the PDF R-factor term
                         (default 0.0; keyword weight_pdf in deadmd.in).
                         The PDF term is active only when exp_pdf is not
                         None and weight_pdf > 0.

    Returns:
        fitnesses: list of floats, one per member (lower is better).
        As a side effect, writes generation_{gen_count}/fitness_scores_gen.dat
        containing one row per member with columns: member, rank, fitness.
        Rank 1 is the best (lowest fitness) and rank N is the worst -- the
        opposite of the rank-based selection method, where the best member
        receives the highest rank (N) to give it the greatest selection
        probability. The rank here is purely a human-readable label for
        inspecting generation output and has no effect on the algorithm.
    """
    SENTINEL = 1e99
    fitnesses = []
    use_pdf   = (exp_pdf is not None and weight_pdf > 0.0)
    use_comp  = composition_options is not None
    if use_comp:
        target_comp, comp_tol, weight_composition = composition_options
    else:
        weight_composition = 0.0

    for i in range(population_size):
        member_dir = str(i + 1)

        totE_file = (f"generation_{gen_count}/{member_dir}"
                     f"/lammps/totE_evolve")
        density_file = (f"generation_{gen_count}/{member_dir}"
                        f"/lammps/density_evolve")
        element_pct_file = (f"generation_{gen_count}/{member_dir}"
                            f"/lammps/element_evolve")
        natoms_file = (f"generation_{gen_count}/{member_dir}"
                       f"/lammps/natoms_evolve")

        total_energy    = read_last_value(totE_file)
        actual_density  = read_last_value(density_file)
        actual_elem_pct = read_last_value(element_pct_file)
        natoms          = read_last_value(natoms_file)
        # A member failed by the validity gate is treated like a crash.
        validity_file = (f"generation_{gen_count}/{member_dir}"
                         f"/lammps/validity_evolve")
        if os.path.exists(validity_file) and read_last_value(validity_file) == 0:
            natoms = SENTINEL

        # Pass crashed members straight through without computing terms.
        if total_energy >= SENTINEL or natoms >= SENTINEL:
            fitnesses.append(SENTINEL)
            with open(f"generation_{gen_count}/{member_dir}"
                      f"/lammps/fitness_scores_gen", "a") as fitness_file:
                fitness_file.write(f"{SENTINEL}\n")
            continue

        # Each member carries its own target density and element percentage
        # in its genome. Using these as per-member references keeps the
        # density and element terms generation-independent.
        target_density  = population[i][1]
        target_elem_pct = population[i][0]

        # Per-atom energy: each reaction removes atoms, so members with
        # different reaction counts have different atom counts and their
        # total energies are not comparable. Relative to energy_ref, the
        # lowest energy per atom of generation 1 (fixed for the run).
        energy_term  = ((total_energy / natoms - energy_ref)
                        / energy_tolerance)
        density_term = (abs(actual_density - target_density)
                        / (target_density + 1e-12))
        element_term = (abs(actual_elem_pct - target_elem_pct)
                        / (target_elem_pct + 1e-12))

        # --- Optional composition term (chi^2 / n) -------------------------
        comp_term = 0.0
        if use_comp:
            comp_file = (f"generation_{gen_count}/{member_dir}"
                         f"/lammps/composition_evolve")
            try:
                with open(comp_file) as f_comp:
                    words = [l.split() for l in f_comp if l.strip()][-1]
                comp = {words[j]: float(words[j + 1])
                        for j in range(0, len(words), 2)}
            except (FileNotFoundError, IndexError, ValueError) as comp_err:
                print(f"[DEAD-MD] composition error gen {gen_count}"
                      f" member {i + 1}: {comp_err}",
                      file=sys.stderr, flush=True)
                fitnesses.append(SENTINEL)
                with open(f"generation_{gen_count}/{member_dir}"
                          f"/lammps/fitness_scores_gen", "a") as ff:
                    ff.write(f"{SENTINEL}\n")
                continue
            comp_term = (sum(abs(comp.get(el, 0.0) - x) / comp_tol[el]
                             for el, x in target_comp.items())
                         / len(target_comp))
            with open(f"generation_{gen_count}/{member_dir}"
                      f"/lammps/composition_term_evolve", "a") as ff:
                ff.write(f"{comp_term}\n")

        # --- Optional bond-limit metric ----------------------------------
        # f = % of El1 atoms bonded to >= 1 El2 atom (final.data); metric =
        # mean over the limits of max(0, f - max_pct) / tolerance, so 0 for
        # a member inside every limit.
        bond_term = 0.0
        if BOND_LIMITS:
            lammps_dir = f"generation_{gen_count}/{member_dir}/lammps"
            try:
                fr = bond_fractions(f"{lammps_dir}/final.data",
                                    [(a, b) for a, b, _, _ in BOND_LIMITS])
            except (FileNotFoundError, ValueError, KeyError) as bond_err:
                print(f"[DEAD-MD] bond metric error gen {gen_count}"
                      f" member {i + 1}: {bond_err}",
                      file=sys.stderr, flush=True)
                fitnesses.append(SENTINEL)
                with open(f"{lammps_dir}/fitness_scores_gen", "a") as ff:
                    ff.write(f"{SENTINEL}\n")
                continue
            bond_term = (sum(max(0.0, fr[(a, b)] - mx) / t
                             for a, b, mx, t in BOND_LIMITS)
                         / len(BOND_LIMITS))
            with open(f"{lammps_dir}/bond_limit_evolve", "a") as ff:
                ff.write(" ".join(f"{a}-{b} {fr[(a, b)]:.4f}"
                                  for a, b, _, _ in BOND_LIMITS) + "\n")
            with open(f"{lammps_dir}/bond_term_evolve", "a") as ff:
                ff.write(f"{bond_term}\n")

        # --- Optional PDF R-factor term ----------------------------------
        # Read pdf_neutron.plot (inside lammps/) and compute a
        # scale-optimised R-factor against the experimental G(r). The
        # optimal scalar s is determined analytically:
        #   s = (G_sim . G_exp) / |G_sim|^2
        # so the comparison is insensitive to the overall amplitude -- only
        # peak positions, relative heights and signs matter.
        # If pdf_neutron.plot is missing or degenerate the member is failed with
        # SENTINEL so the EA discards it via selection/kill_rate.
        pdf_term = 0.0
        if use_pdf:
            rpdf_file = (f"generation_{gen_count}/{member_dir}"
                         f"/lammps/pdf_neutron.plot")
            try:
                sim_data  = np.loadtxt(rpdf_file)
                r_sim     = sim_data[:, 0]
                G_sim     = sim_data[:, 1]
                start_idx = int(np.searchsorted(
                    exp_pdf[:, 0], r_sim[0] - 1e-9))
                exp_g     = exp_pdf[start_idx:
                                    start_idx + len(G_sim), 1]
                pdf_term  = _pdf_rfactor(G_sim, exp_g)
            except (FileNotFoundError, ValueError, OSError) as pdf_err:
                print(f"[DEAD-MD] PDF fitness error gen {gen_count}"
                      f" member {i + 1}: {pdf_err}",
                      file=sys.stderr, flush=True)
                fitnesses.append(SENTINEL)
                with open(f"generation_{gen_count}/{member_dir}"
                          f"/lammps/fitness_scores_gen", "a") as ff:
                    ff.write(f"{SENTINEL}\n")
                with open(f"generation_{gen_count}/{member_dir}"
                          f"/lammps/pdf_rfactor_evolve", "a") as ff:
                    ff.write(f"{SENTINEL}\n")
                continue
            # Write the R-factor so record_* functions can read it
            # without re-loading and re-computing.
            with open(f"generation_{gen_count}/{member_dir}"
                      f"/lammps/pdf_rfactor_evolve", "a") as ff:
                ff.write(f"{pdf_term}\n")

        # Normalize weights so their sum always equals the number of
        # active terms. At equal weights the factor is 1, reproducing
        # the original three-term behaviour exactly. When the PDF term
        # is active num_terms becomes 4 and the total_weight grows
        # accordingly, keeping the fitness magnitude stable.
        # A term switched off by target_composition (weight_element 0 with
        # no target_element lines) is not counted as active.
        use_bond     = bool(BOND_LIMITS)
        num_terms    = (2 + (1 if weight_element > 0.0 or not use_comp else 0)
                        + (1 if use_pdf else 0) + (1 if use_comp else 0)
                        + (1 if use_bond else 0))
        total_weight = (weight_energy + weight_density + weight_element
                        + (weight_pdf if use_pdf else 0.0)
                        + (weight_composition if use_comp else 0.0)
                        + (WEIGHT_BOND if use_bond else 0.0))
        norm         = num_terms / (total_weight + 1e-12)

        fitness = norm * (weight_energy  * energy_term
                          + weight_density * density_term
                          + weight_element * element_term
                          + (weight_pdf * pdf_term if use_pdf else 0.0)
                          + (weight_composition * comp_term
                             if use_comp else 0.0)
                          + (WEIGHT_BOND * bond_term if use_bond else 0.0))
        fitnesses.append(fitness)

        with open(f"generation_{gen_count}/{member_dir}"
                  f"/lammps/fitness_scores_gen", "a") as fitness_file:
            fitness_file.write(f"{fitness}\n")

    # Write a single generation-level summary so all member scores and
    # their ranks are visible in one place without opening individual
    # member directories. Rank 1 is the best (lowest fitness). Crashed
    # members carry the sentinel value 1e99 and naturally rank last.
    sorted_indices = sorted(range(population_size), key=lambda i: fitnesses[i])
    ranks = [0] * population_size
    for rank, idx in enumerate(sorted_indices, start=1):
        ranks[idx] = rank

    fitness_scores_gen = f"generation_{gen_count}/fitness_scores_gen.dat"
    with open(fitness_scores_gen, "w") as summary_file:
        summary_file.write(f"{'member':<10} {'rank':<8} {'fitness':<20}\n")
        for i, score in enumerate(fitnesses):
            summary_file.write(
                f"{i + 1:<10} {ranks[i]:<8} {score:<20.6f}\n")

    return fitnesses


def record_generation_averages(population_size, gen_count, fitness_scores,
                               exp_pdf=None):
    """
    Compute the average fitness, total energy, energy per atom (E/N),
    density, element percentage,
    and (when PDF is active) PDF R-factor across all members of the current
    generation and append one line to averages_per_gen.dat.

    Crashed members are assigned sentinel values of 1e99 by the simulation
    runner and are excluded from all averages so that a single crash does
    not corrupt the generation statistics. If every member crashed, 'nan'
    is written so matplotlib renders a visible gap rather than a spike.

    Args:
        population_size: number of members in the population
        gen_count:       current generation number
        fitness_scores:  list of fitness values for all members
        exp_pdf:         experimental G(r) array or None; controls whether
                         the pdf_rfactor column is collected and written.
    """
    # Any value at or above this threshold is a crash sentinel.
    SENTINEL = 1e90

    total_energies = []
    natoms_list    = []
    densities      = []
    element_pcts   = []
    pdf_rfactors   = []

    for member in range(1, population_size + 1):
        gen_dir = f"generation_{gen_count}/{member}/lammps"
        total_energies.append(
            read_last_value(f"{gen_dir}/totE_evolve"))
        natoms_list.append(
            read_last_value(f"{gen_dir}/natoms_evolve"))
        densities.append(
            read_last_value(f"{gen_dir}/density_evolve"))
        element_pcts.append(
            read_last_value(f"{gen_dir}/element_evolve"))
        if exp_pdf is not None:
            try:
                pdf_rfactors.append(
                    read_last_value(f"{gen_dir}/pdf_rfactor_evolve"))
            except (FileNotFoundError, ValueError):
                pdf_rfactors.append(float("nan"))

    # Filter crashed members out of each quantity independently so a
    # partial crash still produces a meaningful average from survivors.
    valid_fitness   = [v for v in fitness_scores  if v < SENTINEL]
    valid_energies  = [v for v in total_energies  if v < SENTINEL]
    valid_epa       = [e / n for e, n in zip(total_energies, natoms_list)
                       if e < SENTINEL and n < SENTINEL]
    valid_densities = [v for v in densities       if v < SENTINEL]
    valid_elements  = [v for v in element_pcts    if v < SENTINEL]
    valid_rfactors  = ([v for v in pdf_rfactors
                        if v < SENTINEL and not np.isnan(v)]
                       if exp_pdf is not None else [])

    avg_fitness      = (sum(valid_fitness)   / len(valid_fitness)
                        if valid_fitness   else float("nan"))
    avg_total_energy = (sum(valid_energies)  / len(valid_energies)
                        if valid_energies  else float("nan"))
    avg_epa          = (sum(valid_epa)       / len(valid_epa)
                        if valid_epa       else float("nan"))
    avg_density      = (sum(valid_densities) / len(valid_densities)
                        if valid_densities else float("nan"))
    avg_element_pct  = (sum(valid_elements)  / len(valid_elements)
                        if valid_elements  else float("nan"))
    avg_pdf_rfactor  = (sum(valid_rfactors)  / len(valid_rfactors)
                        if valid_rfactors  else float("nan"))

    with open("averages_per_gen.dat", "a") as out_file:
        if gen_count == 1:
            header = (f"{'gen':<6} {'avg_fitness':<14}"
                      f" {'avg_total_energy':<18} {'avg_density':<13}"
                      f" avg_element_pct  avg_energy_per_atom")
            if exp_pdf is not None:
                header += "  avg_pdf_rfactor"
            out_file.write(header + "\n")
        line = (f"{gen_count:<6} {avg_fitness:<14.6f}"
                f" {avg_total_energy:<18.6f} {avg_density:<13.6f}"
                f" {avg_element_pct:<16.6f} {avg_epa:.6f}")
        if exp_pdf is not None:
            line += f"  {avg_pdf_rfactor:.6f}"
        out_file.write(line + "\n")


def record_best_quantities(population, population_size, gen_count,
                           fitness_scores, exp_pdf=None,
                           composition_options=None):
    """
    For each generation, find the member that is best by four independent
    criteria and append one line per criterion to its own dat file in the
    top-level run directory. The four criteria are:

      1. best_by_fitness_per_gen.dat -- lowest overall fitness score
      2. best_by_energy_per_gen.dat  -- lowest total energy (LJ potential)
      2b. best_by_energy_per_atom_per_gen.dat -- lowest energy per atom
      3. best_by_density_per_gen.dat -- actual density closest to the
                                        member's own genome target density
      4. best_by_element_per_gen.dat -- actual element percentage closest
                                        to the member's own genome target
                                        element percentage
      4b. best_by_composition_per_gen.dat -- lowest composition metric
                                        (only with target_composition)

    Each criterion may identify a different member. Recording them
    separately supports the planned fitness function redesign
    (DESIGN.md section 5) where density and element percentage are
    penalised by deviation from target rather than by raw value.

    Args:
        population:      2D numpy array of genome values
                           col 0: target element percentage
                           col 1: target density
                           col 2: cell size
                           col 3: condensation rate
                           col 4: bonding probability
        population_size: number of members in the population
        gen_count:       current generation number
        fitness_scores:  list of fitness values (one per member)
    """
    # Any value at or above this threshold is a crash sentinel, not a
    # real simulation result. If all members crashed, nan is written so
    # the plot shows a gap rather than a sentinel spike.
    SENTINEL = 1e90
    NAN      = float("nan")

    # Read actual simulation results for every member.
    total_energies = []
    densities      = []
    element_pcts   = []
    natoms_list    = []
    pdf_rfactors   = []

    for member in range(1, population_size + 1):
        gen_dir = f"generation_{gen_count}/{member}/lammps"
        total_energies.append(
            read_last_value(f"{gen_dir}/totE_evolve"))
        densities.append(
            read_last_value(f"{gen_dir}/density_evolve"))
        element_pcts.append(
            read_last_value(f"{gen_dir}/element_evolve"))
        natoms_list.append(
            read_last_value(f"{gen_dir}/natoms_evolve"))
        if exp_pdf is not None:
            try:
                pdf_rfactors.append(
                    read_last_value(f"{gen_dir}/pdf_rfactor_evolve"))
            except (FileNotFoundError, ValueError):
                pdf_rfactors.append(SENTINEL)

    # --- Criterion 1: lowest fitness ---
    best_fit_idx    = int(np.argmin(fitness_scores))
    best_fit_member = best_fit_idx + 1
    all_crashed     = fitness_scores[best_fit_idx] >= SENTINEL

    # One line per generation: the lowest-fitness member and its raw
    # quantities (energy_per_atom = total_energy / final atom count, the
    # quantity the energy metric uses; pdf_rfactor only with a PDF file).
    with open("best_by_fitness_per_gen.dat", "a") as out_file:
        if gen_count == 1:
            header = (f"{'gen':<6} {'member':<8} {'fitness':<14}"
                      f" {'total_energy':<14} {'energy_per_atom':<16}"
                      f" {'density':<12} element_pct")
            if exp_pdf is not None:
                header += "  pdf_rfactor"
            out_file.write(header + "\n")
        if all_crashed:
            line = (f"{gen_count:<6} {'nan':<8} {'nan':<14}"
                    f" {'nan':<14} {'nan':<16} {'nan':<12} nan")
            if exp_pdf is not None:
                line += "  nan"
        else:
            i = best_fit_idx
            line = (f"{gen_count:<6} {best_fit_member:<8}"
                    f" {fitness_scores[i]:<14.6f}"
                    f" {total_energies[i]:<14.6f}"
                    f" {total_energies[i] / natoms_list[i]:<16.6f}"
                    f" {densities[i]:<12.6f}"
                    f" {element_pcts[i]:.6f}")
            if exp_pdf is not None:
                line += f"  {pdf_rfactors[i]:.6f}"
        out_file.write(line + "\n")

    # --- Criterion 2: lowest total energy ---
    best_energy_idx    = int(np.argmin(total_energies))
    best_energy_member = best_energy_idx + 1

    with open("best_by_energy_per_gen.dat", "a") as out_file:
        if gen_count == 1:
            out_file.write(
                f"{'gen':<6} {'member':<8} total_energy\n"
            )
        if total_energies[best_energy_idx] >= SENTINEL:
            out_file.write(f"{gen_count:<6} {'nan':<8} nan\n")
        else:
            out_file.write(
                f"{gen_count:<6} {best_energy_member:<8}"
                f" {total_energies[best_energy_idx]:.6f}\n"
            )

    # --- Criterion 2b: lowest energy per atom (E/N) ---
    # Its own file because the lowest-E/N member can differ from the
    # lowest-total-E member (members keep different atom counts). E/N is
    # the quantity the energy metric uses. Crashed members (E or N at the
    # sentinel) are skipped.
    energies_per_atom = [
        total_energies[i] / natoms_list[i]
        if total_energies[i] < SENTINEL and natoms_list[i] < SENTINEL
        else SENTINEL
        for i in range(population_size)
    ]
    best_epa_idx = int(np.argmin(energies_per_atom))

    with open("best_by_energy_per_atom_per_gen.dat", "a") as out_file:
        if gen_count == 1:
            out_file.write(
                f"{'gen':<6} {'member':<8} {'energy_per_atom':<16}"
                f" {'total_energy':<14} natoms\n"
            )
        if energies_per_atom[best_epa_idx] >= SENTINEL:
            out_file.write(
                f"{gen_count:<6} {'nan':<8} {'nan':<16} {'nan':<14} nan\n"
            )
        else:
            out_file.write(
                f"{gen_count:<6} {best_epa_idx + 1:<8}"
                f" {energies_per_atom[best_epa_idx]:<16.6f}"
                f" {total_energies[best_epa_idx]:<14.6f}"
                f" {int(natoms_list[best_epa_idx])}\n"
            )

    # --- Criterion 3: actual density closest to genome target density ---
    density_deviations = [
        abs(densities[i] - population[i][1])
        for i in range(population_size)
    ]
    best_density_idx    = int(np.argmin(density_deviations))
    best_density_member = best_density_idx + 1

    with open("best_by_density_per_gen.dat", "a") as out_file:
        if gen_count == 1:
            out_file.write(
                f"{'gen':<6} {'member':<8} {'actual_density':<16}"
                f" {'target_density':<16} deviation\n"
            )
        if densities[best_density_idx] >= SENTINEL:
            out_file.write(
                f"{gen_count:<6} {'nan':<8} {'nan':<16}"
                f" {'nan':<16} nan\n"
            )
        else:
            out_file.write(
                f"{gen_count:<6} {best_density_member:<8}"
                f" {densities[best_density_idx]:<16.6f}"
                f" {population[best_density_idx][1]:<16.6f}"
                f" {density_deviations[best_density_idx]:.6f}\n"
            )

    # --- Criterion 4: actual element pct closest to genome target ---
    element_deviations = [
        abs(element_pcts[i] - population[i][0])
        for i in range(population_size)
    ]
    best_element_idx    = int(np.argmin(element_deviations))
    best_element_member = best_element_idx + 1

    with open("best_by_element_per_gen.dat", "a") as out_file:
        if gen_count == 1:
            out_file.write(
                f"{'gen':<6} {'member':<8} {'actual_pct':<12}"
                f" {'target_pct':<12} deviation\n"
            )
        if element_pcts[best_element_idx] >= SENTINEL:
            out_file.write(
                f"{gen_count:<6} {'nan':<8} {'nan':<12}"
                f" {'nan':<12} nan\n"
            )
        else:
            out_file.write(
                f"{gen_count:<6} {best_element_member:<8}"
                f" {element_pcts[best_element_idx]:<12.6f}"
                f" {population[best_element_idx][0]:<12.6f}"
                f" {element_deviations[best_element_idx]:.6f}\n"
            )

    # --- Criterion 4b: lowest composition metric (target_composition) ---
    # composition metric = (1/n) sum_El |x_El - x_target,El| / sigma_El,
    # written by fitness_function to composition_term_evolve; the member's
    # at% come from composition_evolve. Crashed members (fitness sentinel)
    # are skipped. Columns: the metric, then actual at% of each target
    # element (target in the header).
    if composition_options is not None:
        target_comp = composition_options[0]
        best_c = None
        for i in range(population_size):
            if fitness_scores[i] >= SENTINEL:
                continue
            gen_dir = f"generation_{gen_count}/{i + 1}/lammps"
            try:
                metric = read_last_value(f"{gen_dir}/composition_term_evolve")
            except (FileNotFoundError, ValueError):
                continue
            if best_c is None or metric < best_c[1]:
                best_c = (i, metric)
        with open("best_by_composition_per_gen.dat", "a") as out_file:
            if gen_count == 1:
                out_file.write(
                    f"{'gen':<6} {'member':<8} {'comp_metric':<12} "
                    + " ".join(f"{el + '(' + format(x, 'g') + ')':<12}"
                               for el, x in target_comp.items()) + "\n")
            if best_c is None:
                out_file.write(f"{gen_count:<6} {'nan':<8} nan\n")
            else:
                i, metric = best_c
                with open(f"generation_{gen_count}/{i + 1}/lammps"
                          f"/composition_evolve") as f_comp:
                    w = [l.split() for l in f_comp if l.strip()][-1]
                comp = {w[j]: float(w[j + 1]) for j in range(0, len(w), 2)}
                out_file.write(
                    f"{gen_count:<6} {i + 1:<8} {metric:<12.6f} "
                    + " ".join(f"{comp.get(el, 0.0):<12.4f}"
                               for el in target_comp) + "\n")

    # --- Criterion 5: lowest PDF R-factor (best structural match) ------
    # Active only when experimental PDF data was provided. Lower R-factor
    # means the simulated G(r) most closely matches the experimental one.
    if exp_pdf is not None:

        # Rank all members by pdf_rfactor (rank 1 = lowest = best).
        # Crashed/missing members carry the sentinel and rank last.
        # The ranking is computed first so the winner's rank is available
        # when writing both the per-generation and simulation-wide files.
        sorted_rf_idx = sorted(
            range(population_size),
            key=lambda i: pdf_rfactors[i]
        )
        rf_ranks = [0] * population_size
        for rank, idx in enumerate(sorted_rf_idx, start=1):
            rf_ranks[idx] = rank

        valid_rf = [(v, i) for i, v in enumerate(pdf_rfactors)
                    if v < SENTINEL and not np.isnan(v)]
        if valid_rf:
            best_pdf_val, best_pdf_idx = min(valid_rf, key=lambda x: x[0])
            best_pdf_member = best_pdf_idx + 1
            best_pdf_rank   = rf_ranks[best_pdf_idx]
        else:
            best_pdf_val    = None
            best_pdf_member = None
            best_pdf_rank   = None

        # Simulation-wide file: one row per generation, best member only.
        # After appending this generation's result, re-read all rows and
        # rewrite the file with cross-generation ranks so rank 1 always
        # identifies the generation with the lowest pdf_rfactor overall.
        sim_pdf_file = "best_by_pdf_per_gen.dat"
        new_row = (
            gen_count,
            best_pdf_member if best_pdf_member is not None else float("nan"),
            best_pdf_val    if best_pdf_val    is not None else float("nan"),
        )

        # Load all previously recorded rows (skip header line).
        prior_rows = []
        try:
            with open(sim_pdf_file, "r") as rf:
                lines = rf.readlines()
            for line in lines[1:]:          # lines[0] is the header
                parts = line.split()
                prior_rows.append((
                    int(parts[0]),
                    parts[1],               # member (may be "nan")
                    float(parts[3]),        # pdf_rfactor (col 3, skip rank)
                ))
        except FileNotFoundError:
            prior_rows = []

        all_rows = prior_rows + [new_row]

        # Rank generations by pdf_rfactor; nan/crashed entries rank last.
        def _sort_key(row):
            v = row[2] if not isinstance(row[2], float) else row[2]
            return (np.isnan(v), v)

        sorted_gens = sorted(range(len(all_rows)), key=lambda i: _sort_key(all_rows[i]))
        gen_ranks   = [0] * len(all_rows)
        for rank, idx in enumerate(sorted_gens, start=1):
            gen_ranks[idx] = rank

        with open(sim_pdf_file, "w") as out_file:
            out_file.write(
                f"{'gen':<6} {'member':<8} {'rank':<6} pdf_rfactor\n"
            )
            for idx, row in enumerate(all_rows):
                g, m, v = row
                if np.isnan(v):
                    out_file.write(
                        f"{g:<6} {'nan':<8} {'nan':<6} nan\n"
                    )
                else:
                    out_file.write(
                        f"{g:<6} {str(m):<8} {gen_ranks[idx]:<6}"
                        f" {v:.6f}\n"
                    )

        # Per-generation file: all members ranked by pdf_rfactor,
        # mirroring fitness_scores_gen.dat. Rank 1 = best (lowest).
        pdf_rfactor_gen = f"generation_{gen_count}/pdf_rfactor_gen.dat"
        with open(pdf_rfactor_gen, "w") as rf_file:
            rf_file.write(f"{'member':<10} {'rank':<8} pdf_rfactor\n")
            for i, rfactor in enumerate(pdf_rfactors):
                rf_file.write(
                    f"{i + 1:<10} {rf_ranks[i]:<8} {rfactor:.6f}\n"
                )


def parents_selection(selection, fitness_scores, kill_rate, population_size, k):
    """
    Select the top performing individuals based on fitness (lower fitness is better
    and means a fitter memberand the chosen selection method.
    Returns indices of slected parents.
    """
    num_parents = int(population_size - (population_size * kill_rate / 100))

    # Make sure number of parents is at least 2 and an even numebr to enable crossover.
    if num_parents < 2:
        num_parents = 2
    if num_parents % 2 != 0:
        num_parents += 1

    # Selection option 1 is Truncation selection
    if selection == 1:
        # Pair index with fitness, sort by fitness ascending
        ranked = sorted(
                enumerate(fitness_scores),
                key=lambda x: x[1])
    

        # Keep the top performers according to the kill_rate

        parents_indices = [indx for indx, score in ranked[:num_parents]]

    # Selection option 2 is rank selection
    elif selection == 2:
        # Rank individuals in an ascending order
        ranked = sorted(
                range(population_size),
                key=lambda i: fitness_scores[i]
        )        

        # assign linear ranks (best gets higher rank)
        ranks = np.array([population_size -i for i in range(population_size)], dtype=float)

        # Normalize to probabilities
        probabilities = ranks / np.sum(ranks)

        # Sample parents positions (with replacement)
        chosen_positions = np.random.choice(
                population_size,
                size=num_parents,
                replace=True,
                p=probabilities
        )

        # Map back to actual population indices
        parents_indices = [ranked[pos] for pos in chosen_positions]

    # Selection option 3 is tournament selection
    elif selection == 3:
        parents_indices = []

        for _ in range(num_parents):

            # Pick k random individuals to compete
            competitors = random.sample(range(population_size), k)

            # Select the fittest of that small tournament (lowest fitness)
            best = min(competitors, key=lambda i: fitness_scores[i])

            # add to parents list
            parents_indices.append(best)

    else:
        raise ValueError("selection must be truncation, rank or tournament")

    return parents_indices


def crossover_method(parent1, parent2, crossover_type):
    """
    Perform a crossover between two parents.

    Args:
        parent1 (list or np.ndarray): Genome of parent 1
        parent2 (list or np.ndarray): Genome of parent 2
    Returns:
    child1, child2 (same type as parents)
    """
    genome_length = len(parent1)

    # One point crossover
    if crossover_type == 1:

        crossover_point = random.randint(1, genome_length -1)

        child1 = np.concatenate((parent1[:crossover_point], parent2[crossover_point:]))
        child2 = np.concatenate((parent2[:crossover_point], parent1[crossover_point:]))
    
    # Two points crossover
    elif crossover_type == 2:
        crossover_point1, crossover_point2 = sorted(random.sample(range(1, genome_length), 2))
        child1 = np.concatenate(
                (parent1[:crossover_point1], parent2[crossover_point1:crossover_point2],
                 parent1[crossover_point2:])
                ) 
        child2 = np.concatenate(
                (parent2[:crossover_point1], parent1[crossover_point1:crossover_point2],
                 parent2[crossover_point2:])
                ) 



    return child1, child2


def mutation(child, mutation_rate, target_element_lower, target_element_upper,
             target_density_lower, target_density_upper, cell_size_lower,
             cell_size_upper, squish_lower, squish_upper,
             rmax_lower, rmax_upper,
             bonding_probability_lower, bonding_probability_upper):

    """
    Mutate a child genome by randomly changing genes with a giving mutation 
    probability.

    Args:
        child (list or np.ndarray): Genome to mutate
        mutation_rate (float): Probability of mutating each gene

    Returns:
        np.ndarray: Mutated genome
    """

    mutated_child = child.copy()

    def new_value(x, lo, hi):
        # One mutated gene: uniform redraw, or a Gaussian creep step
        # reflected back into [lo, hi]. A fixed gene (lo == hi) stays fixed.
        if MUTATION_TYPE == "uniform" or hi <= lo or random.random() < CREEP_JUMP:
            return random.uniform(lo, hi)
        y = x + random.gauss(0.0, CREEP_SIGMA * (hi - lo))
        if y < lo:
            y = 2 * lo - y
        elif y > hi:
            y = 2 * hi - y
        return min(max(y, lo), hi)

    # Every gene is mutated independently with the same probability
    # mutation_rate, each within its own range.
    # gene 0: target element percentage
    if random.random() < mutation_rate:
        mutated_child[0] = new_value(child[0], target_element_lower,
                                     target_element_upper)

    # gene 1: target density
    if random.random() < mutation_rate:
        mutated_child[1] = new_value(child[1], target_density_lower,
                                     target_density_upper)

    # gene 2: cell_size
    if random.random() < mutation_rate:
        mutated_child[2] = round(new_value(child[2], cell_size_lower,
                                           cell_size_upper), 2)

    # gene 3: condensation rate (squish factor)
    if random.random() < mutation_rate:
        mutated_child[3] = round(new_value(child[3], squish_lower,
                                           squish_upper), 2)

    # genes 4 + r: Rmax of reaction r; genes 4 + n + r: its bonding
    # probability, each within that reaction's own range.
    n = len(bonding_probability_lower)
    for r in range(n):
        if random.random() < mutation_rate:
            mutated_child[4 + r] = round(new_value(child[4 + r], rmax_lower[r],
                                                   rmax_upper[r]), 2)
    for r in range(n):
        if random.random() < mutation_rate:
            mutated_child[4 + n + r] = new_value(
                child[4 + n + r], bonding_probability_lower[r],
                bonding_probability_upper[r])

    # genes 4 + 2n + j: max_rxn caps (whole numbers). A uniform redraw picks
    # any whole number in [lo, hi] with equal chance; a creep step is rounded.
    for j, (r, lo, hi) in enumerate(MAXRXN_GENES):
        if random.random() < mutation_rate:
            g = 4 + 2 * n + j
            if MUTATION_TYPE == "uniform" or hi <= lo or \
               random.random() < CREEP_JUMP:
                mutated_child[g] = random.randint(lo, hi)
            else:
                y = child[g] + random.gauss(0.0, CREEP_SIGMA * (hi - lo))
                if y < lo:
                    y = 2 * lo - y
                elif y > hi:
                    y = 2 * hi - y
                mutated_child[g] = min(max(int(round(y)), lo), hi)

    return mutated_child


def create_offsprings(population, parents_indices, mutation_rate,
                      target_element_lower, target_element_upper,
                      target_density_lower, target_density_upper,
                      cell_size_lower, cell_size_upper,
                      squish_lower, squish_upper,
                      rmax_lower, rmax_upper,
                      bonding_probability_lower, bonding_probability_upper,
                      crossover_type, num_offspring):
    """
    Generate exactly num_offspring children from selected parents.

    Two distinct parents are picked at random from parents_indices each
    iteration, crossed over, and mutated. The loop repeats until at least
    num_offspring children exist, then the list is sliced to exactly that
    count. This guarantees the caller always receives a full generation
    worth of children regardless of kill_rate and elitism_rate settings.

    Picking randomly (rather than sequentially) gives every fit parent an
    equal mating probability and avoids pairing bias.

    Args:
        population (list): Current population, where each element is a
            genome.
        parents_indices (list): Indices of selected parents in the
            population.
        num_offspring (int): Exact number of children to produce. Should
            be population_size - elitism_count so that elites + children
            always reconstitute a full-size generation.

    Returns:
        list: Exactly num_offspring child genomes produced via crossover
              and mutation.
    """
    children = []

    while len(children) < num_offspring:
        # Pick two distinct parents at random from the fit pool.
        idx1, idx2 = random.sample(range(len(parents_indices)), 2)
        p1 = population[parents_indices[idx1]]
        p2 = population[parents_indices[idx2]]

        child1, child2 = crossover_method(p1, p2, crossover_type)

        child1 = mutation(child1, mutation_rate, target_element_lower,
                          target_element_upper, target_density_lower,
                          target_density_upper, cell_size_lower,
                          cell_size_upper, squish_lower, squish_upper,
                          rmax_lower, rmax_upper,
                          bonding_probability_lower,
                          bonding_probability_upper)
        child2 = mutation(child2, mutation_rate, target_element_lower,
                          target_element_upper, target_density_lower,
                          target_density_upper, cell_size_lower,
                          cell_size_upper, squish_lower, squish_upper,
                          rmax_lower, rmax_upper,
                          bonding_probability_lower,
                          bonding_probability_upper)

        children.append(child1)
        children.append(child2)

    return children[:num_offspring]


def compute_elitism_count(population_size, elitism_rate):
    """
    Return the number of elite individuals to carry into the next generation.

    Centralises the elitism logic so that create_offsprings (which needs
    to know how many slots are left for children) and new_generation (which
    fills those slots) always agree on the count.

    Args:
        population_size (int): Total number of individuals in the population.
        elitism_rate (float): Percentage of population preserved as elites.

    Returns:
        int: Number of elites (>= 1 when elitism_rate > 0, else 0).
    """
    count = int(elitism_rate * population_size / 100)
    # Ensure at least one elite survives when elitism is enabled.
    if elitism_rate > 0 and count == 0:
        count = 1
    return min(count, population_size)


def new_generation(population, fitness_scores, children, elitism_count):
    """
    Create the next generation population.

    Args:
        population (list): Current population (list of genomes)
        fitness_scores (list): Fitness score per individual
        children (list): Offspring genomes
        elitism_count (int): Number of top individuals to carry over
            unchanged. Compute this with compute_elitism_count() before
            calling.

    Returns:
        list: New generation
    """
    population_size = len(population)

    ranked_indices = sorted(
        range(population_size),
        key=lambda i: fitness_scores[i]
    )

    elites = [population[i] for i in ranked_indices[:elitism_count]]

    # ---- Fill the rest with children ----
    new_gen = elites.copy()

    for child in children:
        if len(new_gen) < population_size:
            new_gen.append(child)
        else:
            break

    return new_gen


def main():

    print_banner()

    (
        population, 
        population_size, 
        mutation_rate, 
        kill_rate, 
        elitism_rate,
        target_element, 
        num_cores, 
        selection, 
        gen_max,
        target_element_lower,
        target_element_upper,
        target_density_lower,
        target_density_upper,
        cell_size_lower, 
        cell_size_upper,
        squish_lower,
        squish_upper,
        rmax_lower,
        rmax_upper,
        bonding_probability_lower, 
        bonding_probability_upper,
        crossover_type,
        max_speed,
        stages,
        k,
        molecules,
        num_molecule,
        composition_num,
        reactions_num,
        rxns,
        rxn_options,
        ref_energy,
        weight_energy,
        weight_density,
        weight_element,
        exp_pdf,
        weight_pdf,
        pdf_options,
        gate_options,
        cores_per_member,
        lmp_mode,
        composition_options,
        energy_tolerance
    ) = initialize_population()

    stamp_run()

    gen_count = 1
    elite_sources = {}
    energy_ref = None   # generation 1's lowest energy per atom, then fixed

    # -----GA loop excution -----
    while gen_count <= gen_max:
        print(f"---- Current generation is: {gen_count}/{gen_max} -----", flush=True)
        print(f"---- Population # in current gen is: {len(population)} ---", flush=True)
        run_lammps_simulations(population, population_size, num_cores,
                               target_element, molecules, num_molecule,
                               composition_num, reactions_num, rxns,
                               rxn_options, max_speed, stages, gen_count, exp_pdf,
                               pdf_options, gate_options, elite_sources,
                               cores_per_member, lmp_mode)
        # Fix the energy reference at the first generation that has a
        # scored member (normally generation 1).
        if energy_ref is None:
            energy_ref = best_energy_per_atom(population_size, gen_count)
            if energy_ref is not None:
                print(f"[DEAD-MD] Energy reference fixed at generation "
                      f"{gen_count}'s lowest energy per atom: {energy_ref}"
                      f" (energy_tolerance {energy_tolerance})",
                      flush=True)
        fitness_scores = fitness_function(population, population_size,
                                          gen_count,
                                          energy_ref if energy_ref is not None
                                          else 0.0,
                                          weight_energy, weight_density,
                                          weight_element, exp_pdf,
                                          weight_pdf, composition_options,
                                          energy_tolerance)
        record_generation_averages(population_size, gen_count,
                                   fitness_scores, exp_pdf)
        record_best_quantities(population, population_size, gen_count,
                               fitness_scores, exp_pdf, composition_options)
        print (fitness_scores)
        parents_indices = parents_selection(
            selection, fitness_scores, kill_rate, population_size, k)
        print(parents_indices)

        elitism_count = compute_elitism_count(population_size, elitism_rate)
        num_offspring = population_size - elitism_count

        children = create_offsprings(population, parents_indices,
                                     mutation_rate,
                                     target_element_lower,
                                     target_element_upper,
                                     target_density_lower,
                                     target_density_upper,
                                     cell_size_lower, cell_size_upper,
                                     squish_lower, squish_upper,
                                     rmax_lower, rmax_upper,
                                     bonding_probability_lower,
                                     bonding_probability_upper,
                                     crossover_type, num_offspring)
        # print(f"Children are {children}\n")
        population = new_generation(population, fitness_scores, children,
                                    elitism_count)
        # new_generation puts the elites first, in fitness order. Record
        # where each one's results are so the next generation copies them
        # instead of rerunning LAMMPS.
        ranked = sorted(range(population_size),
                        key=lambda i: fitness_scores[i])
        elite_sources = {
            pos + 1: os.path.abspath(f"generation_{gen_count}/{src + 1}")
            for pos, src in enumerate(ranked[:elitism_count])}
        # print(f"new population is {population}\n")

        gen_count += 1

        
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

---
layout: default
title: deadmd.py
script: deadmd.py
description: 'DEAD-MD, a genetic algorithm that searches the parameters of a condense.py plus LAMMPS bond/react condensation (cell size, compression rate, per-reaction probabilities, Rmax and caps) for amorphous models matching a target composition and density, or minimizing energy. Runs members concurrently as srun steps inside one SLURM job, driven by a deadmd.in keyword file.'
parent: Script Catalog
---

# deadmd

DEAD-MD (Dynamic Evolutionary Algorithm-Driven Molecular Dynamics) is
a genetic algorithm (GA) that searches the parameters of a
[condense.py](condense.html) plus LAMMPS `fix bond/react` condensation
for amorphous models that match a target. The searched parameters are
the cell size, the compression rate, and each reaction's probability,
Rmax and event cap. The target is either an experimental composition
and density (replica runs) or low energy at a chosen composition
range (exploratory runs). Every model is built through reaction
chemistry, so every member of the population is a physically
reachable network.

Author: Mohammed Belhadj Larbi, University of Missouri-Kansas City.

## Algorithm

1. Draw an initial population of genomes, seeded by `ga_seed`.
2. For every member, write `condense.in` from its genome, run
   `condense.py` (seeded per member) and LAMMPS. Members run
   **concurrently** on `num_cores` cores, each in its own `srun`
   step.
3. Read each member's results (energy, atom count, density,
   composition, bond fractions, optional G(r) and validity check)
   and compute its fitness. Lower is better. A crashed or invalid
   member gets the sentinel 1e99.
4. Select parents (truncation, rank, or tournament of k), make
   children by one- or two-point crossover, and mutate (uniform
   redraw or Gaussian creep).
5. The next generation is the elites, kept and *not* re-simulated
   but copied and re-scored, plus the children.
6. Repeat until `gen_max`. Convergence criteria are planned; until
   then the run always goes to `gen_max`.

## Usage

**Commandline:**

```
deadmd.py
```

There are no command-line options. The script reads `deadmd.in` from
the current directory and must be launched **inside a SLURM job**,
since it starts every member as an `srun` step and reads the job's
task count from the environment. A typical batch script requests
`num_cores` tasks and runs `deadmd.py` as its only command. Progress
goes to the job's standard output; each member's own condense and
LAMMPS output goes to `member.out` in its directory.

## The Genome

One row per member:

| Gene         | Meaning |
|:------------:|:--------|
|0             | Target element atomic percent (single-element metric only). |
|1             | Target density: condense.py stops compressing there. In replica runs set lower = upper to the experimental density; it is then a fixed stop value, not searched. |
|2             | Cell size in Angstroms. |
|3             | Squish factor of the genome stage, i.e. the condensation rate. |
|4 .. 3+n      | Rmax of reaction r, for n = `reactions_num` reactions. |
|4+n .. 3+2n   | Probability of reaction r. |
|4+2n ..       | The `max_rxn` cap of each reaction given with `max_rxn LO HI`, in line order, as whole numbers. |

Every gene is drawn and mutated inside its own range; setting
lower = upper fixes it.

## Fitness

```
fitness = norm * sum_i (w_i * metric_i),   norm = n_active / sum_i w_i
```

Each metric is dimensionless with 0 meaning on target. The weights
are the `weight_*` keywords (default 1.0), and `norm` keeps the
fitness scale the same whatever the number of active metrics. Set a
weight to 0.0 to ignore its metric.

| Metric       | Definition |
|:------------:|:-----------|
|energy        | (E/N minus the lowest E/N of generation 1) / `energy_tolerance`. Signed: members below the reference score negative, so lower energy always lowers the fitness. E is the minimized total energy and N the final atom count, since reactions delete atoms. |
|density       | abs(rho minus gene 1) / gene 1. Only checks that the density stop was reached; near 0 whenever it was. |
|element       | abs(at% minus gene 0) / gene 0. Legacy single-element metric; off when only `target_composition` is used. |
|composition   | (1/n) sum over listed elements of abs(x minus x_target) / sigma. Optional. |
|bond          | (1/n) sum of max(0, f minus max_pct) / tolerance. Optional. |
|pdf           | Scale-optimised R-factor against `exp_pdf_file`. Optional. |

## The Input File

Every keyword in `deadmd.in` is case-insensitive, lines are split on
whitespace, and the order of lines does not matter. The groups below
follow the annotated example in the script's own docstring, which is
the authoritative reference.

**Resources and seed**

| Keyword             | Meaning |
|:--------------------|:--------|
|`num_cores N`        | Members run concurrently on N cores. N must not exceed the job's task count. By default the cores are shared out automatically: with at least N members each gets one core and N run at once; with fewer, all run at once and share the cores. |
|`cores_per_member N` | Fix the MPI ranks of every member instead. Optional. Runs with `inter_unlinked` reactions need 1. |
|`lmp parallel|serial`| How the one LAMMPS binary is launched: several MPI ranks per member (default) or one rank, which forces `cores_per_member 1`. |
|`ga_seed N`          | Seeds the GA's own random numbers and, through a per-member seed derived from (seed, generation, member), every packmol, bond/react and velocity draw in condense.py. Optional; when omitted a seed is drawn, printed, and kept in `run_stamp/ga_seed`. |

**GA control**

| Keyword                                  | Meaning |
|:-----------------------------------------|:--------|
|`gen_max N`, `population_size N`          | Generations and members per generation. |
|`elitism_rate P`, `mutation_rate P`, `kill_rate P` | Percents: of the population kept as elites (at least 1 if above 0), chance per gene of a mutation, and of the population not allowed to be parents. |
|`mutation_type uniform|creep`, `creep_sigma S`, `creep_jump P` | Uniform redraw over the range (default), or a Gaussian step of standard deviation S times the range, with P percent of mutations still a uniform redraw. |
|`crossover_type 1|2`                      | One- or two-point crossover. |
|`selection truncation|rank|tournament K`  | Parent selection scheme. |

**Targets and metrics**

| Keyword                                              | Meaning |
|:-----------------------------------------------------|:--------|
|`target_element EL`, `target_element_lower`, `_upper` | Range of gene 0 for the single-element metric. Optional once `target_composition` is given. |
|`target_density_lower`, `target_density_upper`        | Range of gene 1. |
|`target_composition EL PCT EL PCT ...`                | Composition metric: the member's atomic percent of each listed element against the target, in units of a tolerance sigma. Any subset of elements may be listed; the values must not sum to more than 100. |
|`composition_tolerance EL S ...`                      | Override the default sigma of 1 at% per element and 2 at% for H. |
|`bond_limit EL1 EL2 MAX_PCT [tolerance T]`            | The percent of EL1 atoms bonded to at least one EL2 atom must stay at or below MAX_PCT. |
|`energy_tolerance T`                                  | Denominator of the energy metric, in the LAMMPS energy unit per atom. Default 0.5. |
|`weight_energy`, `weight_density`, `weight_element`, `weight_composition`, `weight_bond`, `weight_pdf` | Metric weights. All default to 1.0; only their ratios change the ranking. |
|`exp_pdf_file PATH`, `pdf_weighting neutron|none`, `pdf_sigma S` | Reference reduced PDF G(r) for the pdf metric: two columns, r in Angstrom at 0.01 spacing reaching at least 10.00, and G(r). Each member's G(r) is computed by [pdf_neutron.py](pdf_neutron.html); make a calculated reference with the same weighting and sigma. |
|`validity_gate yes`, `valence EL N ...`, `max_strain S` | Run [deadmd_validity.py](deadmd_validity.html) on each member's final structure. A member that fails is invalid, not worse, and gets the crash sentinel. |

**Condensation**

| Keyword                              | Meaning |
|:-------------------------------------|:--------|
|`composition_num N` then N lines `MOL FAMILY COUNT` | The molecules to pack, as in condense.in. |
|`cell_size_lower`, `cell_size_upper`  | Range of gene 2. |
|`squish_lower`, `squish_upper`        | Range of gene 3. Default 0.15 to 0.33. A member can reach its target density only if (starting density) / squish^3 is at least the target. |
|`max_speed V`                         | Maximum initial molecule speed. |
|`stage STEP ENSEMBLE T0 T1 TDAMP STEPS` | The genome stage: its squish factor comes from gene 3 and compression stops when gene 1's density is reached. Use a small STEP (e.g. 10) so the stop is precise. |
|`stage SQUISH STEP ENSEMBLE T0 T1 TDAMP STEPS` | Additional stages with an explicit squish factor. Use 1.0 for an equilibration hold; every input should end with one. |
|`reactions_num N` then N lines        | `mol1 site1 mol2 site2 prob_lo prob_hi [rmax_lo rmax_hi] inter|intra [rmin R] [max_rxn LO HI]`. Each reaction gets its own probability and Rmax genes drawn within its ranges; `inter`/`intra` and `rmin` are passed to condense.in as is; `max_rxn LO HI` makes the event cap a gene. |

## Output

| Location | Contents |
|:---------|:---------|
|`run_stamp/`                       | Copies of deadmd.py, condense.py and deadmd.in that made the run, plus `ga_seed`, so the run can be read or repeated after the scripts change. |
|`averages_per_gen.dat`             | Per-generation averages. |
|`best_by_{fitness,energy,energy_per_atom,density,element,composition,pdf}_per_gen.dat` | The best member of each generation by each criterion. |
|`generation_<g>/fitness_scores_gen.dat` | All members of generation g, ranked. |
|`generation_<g>/<m>/`              | One member: `condense.in`, `seeds`, `member.out`, `ELITE_COPY_OF` for elite copies, and `lammps/` with `final.data` and the `*_evolve` files (total energy, atom count, density, composition, reaction counts, bond limit). |

Plot the run with [plot_deadmd.py](plot_deadmd.html). The same
`ga_seed`, `deadmd.in`, LAMMPS binary and ranks per member reproduce
a run; a different binary or rank count does not.

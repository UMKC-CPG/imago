# Vision

## Purpose

Imago is an all-electron electronic structure code that
implements the OLCAO (Orthogonalized Linear Combination of
Atomic Orbitals) method.  It uses periodic boundary conditions
and an LCAO basis to compute electronic properties of a broad
range of material systems: crystals, amorphous solids,
nanoparticles, molecules, interfaces, and grain boundaries.
It is an academic code used for research and student training.

The current development focus comprises four concurrent prongs.
Two of them address Brillouin-zone integration: implementing the
Linear Analytic Tetrahedral (LAT) method, and correcting the
treatment of eigenvector-dependent quantities under IBZ symmetry
reduction. The third addresses SCF startup: augmenting today's
isolated-atom initial potential database with potentials drawn
from converged solid-state calculations, opening a path toward
future environment-aware prediction of the starting potential.
The fourth is infrastructural: a shared system for submitting,
tracking, and harvesting batches of Imago calculations on an
HPC cluster, supporting both near-term needs (the database
build, convergence sweeps, validation harnesses) and longer-term
goals (ab-initio molecular dynamics, high-throughput screening).

## Goals

1. **Implement LAT k-point integration.** Replace the current
   Gaussian-broadening approach for Brillouin-zone integration with
   the Linear Analytic Tetrahedral method (Bloechl, Jepsen, &
   Andersen, PRB 49, 16223, 1994). This provides exact analytic
   integration within each tetrahedron, eliminating the arbitrary
   broadening parameter and improving accuracy at lower k-point
   densities.
2. **Correct partial properties under IBZ reduction.** Ensure that
   eigenvector-dependent quantities (PDOS, bond order, effective
   charge) are computed correctly when the k-point mesh is reduced
   to the irreducible Brillouin zone. The current approach of
   multiplying IBZ contributions by star multiplicity is incorrect
   for these quantities.
3. **Augment the initial SCF potential database.** Extend today's
   isolated-atom potential database with potentials taken from
   converged Imago calculations on a curated set of real solids
   (drawn from the Crystallography Open Database and/or the
   Materials Project). The original single-atom potentials are
   preserved -- both as a fallback when no improved potential is
   available, and as the baseline against which improvements are
   measured. The first deliverable is one improved potential per
   element. The database grows progressively, not exhaustively:
   new (element, environment) entries are added as evidence
   warrants, with no attempt to cover every combination. The
   schema and lookup interface must be designed from the outset
   to permit future keying by additional attributes -- local
   environment, expected charge state, coordination -- without
   disturbing callers. Two payoffs: faster SCF convergence on
   typical systems, and the ability to run useful non-SCF
   (single-pass) calculations whose results remain informative
   for many purposes.
4. **Enable high-throughput Imago calculation campaigns.**
   Provide shared infrastructure for submitting, tracking, and
   harvesting batches of Imago calculations on an HPC cluster,
   so that scripts whose purpose is *what* to compute and *what*
   to learn -- initial-potential database build, convergence
   sweeps, validation harnesses, future ab-initio molecular
   dynamics, future high-throughput screening -- do not each
   reinvent *how* to submit and watch jobs. The infrastructure
   rests on three named layers. `imago.py` continues to run one
   Imago calculation per invocation, refactored to expose a
   callable Python API alongside its existing CLI entry point.
   An ASE Calculator wraps `imago.py` and adapts imago to the
   broader computational-materials community; committing to ASE
   yields LAMMPS interoperability, ab-initio molecular dynamics
   through ASE's MD integrators, trivial convergence-test loops,
   and acceptance by every downstream workflow tool that
   consumes ASE -- essentially for free. Parsl orchestrates
   calculations on SLURM as the v1 dispatch backend; its
   Python-as-workflow idiom handles both embarrassingly parallel
   sweeps (thousands of independent jobs) and tightly iterative
   inner loops (AIMD, adaptive convergence) under one mental
   model. **Snakemake may join later** as an outer orchestration
   layer if DAG-of-codes workflows -- e.g., SCF -> bands ->
   optics chains spanning many structures -- become a regular
   need; in that scenario the inner Parsl logic is unchanged
   and Snakemake is additive rather than a replacement.

## Design Principles

1. **Eigenvalue symmetry != eigenvector symmetry.** Point-group
   operations guarantee e_n(k) = e_n(Rk), but eigenvectors at Rk
   are related to those at k by a basis transformation that mixes
   orbital coefficients. Any integration scheme that weights
   eigenvector-dependent quantities by IBZ star multiplicity alone
   is incorrect. This was confirmed empirically: bond orders for
   symmetry-equivalent O atoms in KNbO3 disagreed under IBZ
   reduction but agreed with the full mesh.
2. **Full mesh for post-processing.** Use the IBZ for the expensive
   SCF diagonalization, but reconstruct the full BZ mesh before
   computing partial properties. For eigenvalue-only quantities
   (TDOS), unfolding eigenvalues by symmetry suffices. For
   eigenvector-dependent quantities, either use the full mesh
   directly or apply atom-permutation corrections.
3. **Phased implementation.** Implement LAT TDOS first (eigenvalues
   only), validate against Gaussian-broadened results at high
   k-point density, then extend to integrated partial properties
   via electronPopulation_LAT, and finally to energy-resolved
   PDOS with cornerIntgWt_LAT. Each
   phase must be validated before proceeding to the next.
4. **Stable potential representation, flexible container.** The
   per-potential numerical representation (the existing radial /
   Gaussian-fit form) is preserved. How potentials are stored,
   labeled, and combined into a starting potential is open to
   redesign -- one file per (element, label) pair, one file per
   element containing many labeled potentials, or other
   arrangements as needed. Downstream consumers see only the
   assembled starting potential, not how it was retrieved or
   composed, so the labeling scheme and retrieval logic can grow
   more sophisticated (e.g., interpolation across nearby labels
   when no exact match exists) without disturbing callers.
5. **Curation and regeneration are both part of the deliverable.**
   The value of the improved-potential prong depends directly on
   the breadth and quality of the reference solids used to
   generate the database. A documented curation strategy -- which
   structures, why, and how to extend the set as new chemistries
   are encountered -- is required from the start, not added as an
   afterthought. The database itself must be regeneratable from
   the curated set, not a hand-edited artifact, so that adding
   entries or rebuilding after a code change is a scriptable
   operation rather than a manual exercise.
6. **Cost discipline at SCF setup.** Any prediction or lookup
   performed at SCF initialization must remain cheap relative to
   the SCF iterations it is meant to shorten or replace. This
   bounds how elaborate the environment-keyed prediction step can
   grow, and guards against over-engineering the prediction layer
   for marginal accuracy gains.
7. **Iteration count is the headline metric.** Every reference
   calculation that contributes to the database records its SCF
   iteration count. Across a representative benchmark, runs that
   start from improved initial potentials should take at least
   20% fewer iterations on average than runs starting from the
   original isolated-atom potentials. The original potentials are
   kept in the database precisely so this comparison can be made
   on the same code path, with only the starting potential
   changed.
8. **Decouple adapter from orchestrator.** The ASE Calculator
   adapts imago to the materials-simulation community; Parsl
   orchestrates calculations on the cluster. These two layers are
   independent. Replacing either should not require changing the
   other, and new adapters (other ASE-compatible calculators) or
   new orchestrators (Snakemake on top, or a different dispatch
   backend) can slot in along well-defined seams.
9. **General-purpose orchestration.** Workflow tools are chosen
   so that the skills students acquire while developing or running
   Imago transfer beyond computational materials. Domain-specific
   machinery lives at the adapter layer (ASE); the campaign layer
   (Parsl) is ordinary scientific Python. A student who leaves
   the field carries portable workflow skills rather than
   Imago-specific muscle memory.
10. **Complete-and-report at the campaign level.** A single
    failed calculation never fails an entire campaign by default.
    Per-job failures are recorded and surfaced; the campaign
    driver script decides whether the aggregate outcome is
    acceptable for its scientific purpose. This applies whether
    the cause is convergence non-attainment, cluster-side job
    loss, or post-processing error.

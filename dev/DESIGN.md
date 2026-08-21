# Design

> **Document hierarchy:** VISION -> ARCHITECTURE -> **DESIGN**
> -> PSEUDOCODE -> Code. For goals and principles, see
> `VISION.md`. For repository layout and module map, see
> `ARCHITECTURE.md`.

---

## 1. LAT K-Point Integration

### 1.1 Overview and Motivation

The current Brillouin-zone integration uses Gaussian broadening:
each eigenvalue is smeared by a Gaussian of width sigma, and the
contributions are summed with k-point weights. This introduces
an arbitrary broadening parameter that affects the shape of the
DOS and requires dense k-point meshes for convergence.

The Linear Analytic Tetrahedral (LAT) method (Bloechl, Jepsen,
& Andersen, PRB 49, 16223, 1994) decomposes the BZ into
tetrahedra and integrates analytically within each one. This
eliminates the broadening parameter and provides better accuracy
at lower k-point densities.

### 1.2 Tetrahedra Generation

The uniform Monkhorst-Pack mesh defines a grid of nA x nB x nC
parallelepipeds. Each parallelepiped has 8 corners. Six
tetrahedra sharing one of the box's four long diagonals tile
the box without overlap, and which diagonal is chosen is free.
Taking the M1-M8 diagonal gives the standard set (Bloechl
1994):

```
Parallelepiped corners at grid position (a, b, c):
  M1 = (a,   b,   c  )    M5 = (a+1, b+1, c  )
  M2 = (a+1, b,   c  )    M6 = (a+1, b,   c+1)
  M3 = (a,   b+1, c  )    M7 = (a,   b+1, c+1)
  M4 = (a,   b,   c+1)    M8 = (a+1, b+1, c+1)

Six tetrahedra sharing diagonal M1-M8:
  T1: M1, M2, M5, M8      T4: M1, M4, M7, M8
  T2: M1, M3, M5, M8      T5: M1, M4, M6, M8
  T3: M1, M3, M7, M8      T6: M1, M2, M6, M8
```

The rule behind that table generalizes to any diagonal, and
the implementation needs the general form. Walk from one end
of the chosen diagonal to the other changing one coordinate
at a time; the four corners visited are a tetrahedron, and
the six orders in which the three coordinates can be changed
give the six tetrahedra, each containing both ends. Starting
from M1 this reproduces the table above exactly.

The mesh is periodic (the BZ is a torus), so indices wrap with
modular arithmetic: `mod(a, nA) + 1` etc.

**Every box is cut all four ways.** Choosing one diagonal for
the whole mesh, as the standard presentation does, produces a
decomposition that the crystal point group does not carry
onto itself: all six tetrahedra of a box share one diagonal,
a box has four, and an operation carrying the chosen diagonal
onto a different one produces a decomposition with nothing in
common with the original. Imago therefore cuts each box once
per diagonal and gives each of the resulting 24 tetrahedra a
quarter of the weight. Each of the four cuts tiles the box on
its own, so their average is a valid quadrature, and an
operation that maps one diagonal to another now finds that
cut already present. The count of diagonals used is an input,
defaulting to all four; see the end of this section.

  numTetrahedra = 6 * numTetraDiagonals * nA * nB * nC

All tetrahedra span an equal fraction of the BZ:

  tetraVol = 1 / numTetrahedra

**The quarter weight needs no separate factor.** With four
times as many tetrahedra, `tetraVol` is a quarter of what it
was, and each of the four cuts contributes 6*nA*nB*nC
tetrahedra at that weight -- one quarter of the zone each,
summing to one. Every consumer loops `do t = 1, numTetrahedra`
and scales by `tetraVol`, so the change is confined to
`generateTetrahedra` and no consumer is touched. The same
sentence covers the one-diagonal setting: the weight follows
the count without anyone having to know the count.

**What this fixes, measured.** On a 4x4x4 cubic mesh, all 48
operations leave the four-diagonal decomposition invariant,
against 12 of 48 for the single-diagonal form, whose worst
case maps NONE of its 384 tetrahedra onto a tetrahedron in
the original set. On a hexagonal mesh it improves matters
without curing them: 8 of 24 against 4 of 24. The argument
needs an operation to carry a grid box onto a grid box, which
holds when the operations permute the mesh axes up to sign
(cubic, tetragonal, orthorhombic) and fails for hexagonal and
rhombohedral, where a six-fold rotation sends one axis onto
the SUM of two, so the image of a box is a sheared
parallelepiped that is not a box of the grid. No choice of
diagonal repairs a decomposition whose boxes are not
preserved. Section 1.7 covers the remainder.

**Why not the shortest diagonal.** The literature recommends
choosing the shortest diagonal per box, and that advice is
about interpolation accuracy -- more compact tetrahedra
represent the bands better -- not about symmetry. On a cubic
mesh all four diagonals are the same length, so the rule
gives no guidance; the most even-handed per-box rule
available measured 2 of 48, WORSE than the global choice, and
left mesh points belonging to unequal numbers of tetrahedra
(16, 24 and 40). Averaging all four includes the shortest
alongside the others, which gives up some interpolation
accuracy on a strongly skewed mesh in exchange for a
decomposition with no preferred direction.

**Cost.** Four times the tetrahedron loop in every consumer,
and four times the `tetrahedra` array: 4 x 24 x numFullMeshKP
integers. Totals shift slightly, because this is a different
quadrature, so stored tetrahedron baselines move; the shift
is a change of quadrature and not a regression.

**The single-diagonal form stays available, per phase.** Four
times a per-iteration cost is worth avoiding where it buys
nothing, and there are two cases where it provably buys
nothing.

  1. A run that asks only for totals. A total sums over
     exactly the set the asymmetry permutes, so it cancels
     (section 1.7); the asymmetry has nothing to act on.
  2. The SCF occupation path, whatever the run asks for.
     `computeElectronPopulation_LAT` pools every corner's
     weight onto that corner's IBZ representative, so every
     consumer reading it per k-point already sees the star
     sum, which is symmetric. This is the same structural
     argument that makes effective charge and bond order
     immune (section 1.7), and it matters here because those
     weights are rebuilt every SCF iteration.

So the number of diagonals is an input, read from the k-point
file beside `KPOINT_INTG_CODE`:

```
NUM_TETRA_DIAGONALS
4                        ! 4 = all four (default), 1 = one
```

It belongs there rather than in the shared control section
for two reasons. It is a property of the k-point integration,
which is what that file describes. And the SCF and post-SCF
phases read separate k-point files, so putting it there gives
per-phase control for nothing: a run may take one diagonal
for the SCF occupation and all four for the post-SCF
spectra. `makeinput.py` exposes it as a sibling of the
existing `-scfkpint` and `-pscfkpint` options.

The default is 4 in BOTH phases. The economy above is real
but it rests on an argument about the consumers that exist
today, and a later consumer of the LAT occupations would
inherit the bias silently. A user who wants the saving asks
for it and can be pointed at the reasoning; a user who does
nothing gets the isotropic quadrature.

Adding the tag is an input FORMAT change, so `makeinput.py`
and the `skl/` examples move with it, and any other producer
of these files must be checked.

**What the two settings do and do not guarantee.** Do not
read `NUM_TETRA_DIAGONALS 1` as "the broken one". The
equality of symmetry-equivalent atoms is delivered by the
symmetrization of section 1.7, which is on by default and
works with either setting. What four diagonals adds is
isotropy of the quadrature ITSELF, which shows up where there
is no orbit to average over -- low-symmetry cells -- and in
the shape of a spectrum rather than in the equality of atoms.
One diagonal with symmetrization is a coherent choice; one
diagonal without it is the configuration C149 was opened
about.

The tetrahedra reference the FULL uniform mesh, not the
IBZ-reduced kpoints. Under Option A, eigenvalues are
computed at IBZ points during SCF, then unfolded to the
full mesh for post-processing via a mapping array.

The `generateTetrahedra` call must happen AFTER the axial
kpoint counts are known (after `computeAxialKPoints` for
style code 2, or after reading for style codes 0 and 1).
It is called from `initializeKPoints`, not `readKPoints`.

`tetraVol` is computed in `initializeKPoints` after the
tetrahedra are generated.

**New module-level data in O_KPoints:**
```fortran
integer :: numTetrahedra
integer :: numFullMeshKP
integer :: numTetraDiagonals
    ! How many of the box's four long diagonals are cut
    ! along: 4 by default, 1 for the cheaper single-cut
    ! decomposition. Read from the k-point file, so the
    ! SCF and post-SCF phases carry their own value.
real(kind=double) :: tetraVol
integer, allocatable, dimension(:,:) :: tetrahedra
    ! (4, numTetrahedra) -- indices into the full mesh
integer, allocatable, dimension(:) :: fullKPToIBZKPMap
    ! (numFullMeshKP) -- maps each full mesh kpoint
    ! to its IBZ representative index
```

The `fullKPToIBZKPMap` is produced by the IBZ folding in
`initializeKPointMesh` (from the `kPointTracker` array).
For each full-mesh index i, `fullKPToIBZKPMap(i)` gives
the IBZ kpoint index. Eigenvalue lookup for the TDOS
becomes:
  `eigenValues(band, fullKPToIBZKPMap(fullMeshIdx), spin)`

### 1.3 LAT TDOS (Eigenvalues Only)

The LAT TDOS loop structure inverts relative to Gaussian
broadening:

- **Gaussian (current):** outer = k-points, inner = states,
  innermost = energy bins (Gaussian smear)
- **LAT:** outer = bands, middle = tetrahedra, innermost =
  energy bins (analytic formula)

For each band n and each tetrahedron, sort the 4 corner
eigenvalues: e1 <= e2 <= e3 <= e4. The DOS contribution g(E)
per tetrahedron (Lehmann-Taut / Bloechl):

| Range          | Formula                                    |
|----------------|--------------------------------------------|
| E < e1         | 0                                          |
| e1 <= E < e2   | 3(E-e1)^2 / [(e2-e1)(e3-e1)(e4-e1)]       |
| e2 <= E < e3   | Bloechl eqs. 14-16 (see reference below)   |
| e3 <= E < e4   | 3(e4-E)^2 / [(e4-e1)(e4-e2)(e4-e3)]       |
| E >= e4        | 0                                          |

The middle range (e2 <= E < e3) has a more complex formula
involving cross-terms between all four eigenvalues. The exact
expressions are in Bloechl 1994, equations 14-16.

The Bloechl correction terms improve accuracy at lower k-point
densities, and an earlier draft of this section asked for them
here. They do not belong here: the correction applies to the
cumulative INTEGRATION weights, not to the DOS curve these
formulas produce. It is specified in 1.3.1 and consumed by 1.5
and 1.6.

**Degenerate-corner guards:** When two or more corner energies
coincide, denominators in the analytic formulas go to zero. All
formulas must include guards: `if (abs(e2-e1) < eps) then ...`.

**Units:** Eigenvalues are in Hartree. The DOS output is in
states/eV. The Bloechl formulas give DOS in units of
1/(energy units of e), so a Hartree-to-eV conversion factor
is needed in the output, consistent with the current code's
use of `sigmaSqrtPi / hartree`.

**Unified corner DOS subroutine.** The per-tetrahedron TDOS
g(E) is the *sum* of four per-corner density weights
dw_c/dE (the energy derivatives of the Bloechl cumulative
corner weights, eqs. 18-21). These per-corner derivatives
are the more fundamental quantity: the TDOS needs only
their sum, but the PDOS (section 1.4) needs them
individually to weight each corner's Mulliken projection.

A single subroutine
`bloechlCornerDOSWt(E, eps, cornerDOSWt_LAT)` computes
the four dw_c/dE values. It serves both paths:

  TDOS:  dosContrib = sum(cornerDOSWt_LAT(1:4))
  PDOS:  pdosComp(alpha) += cornerDOSWt_LAT(c)
         * V_T * proj(alpha, n, k_c)

This eliminates duplicated case logic. The inline
dosContrib formulas in `computeTDOS_LAT` are replaced by
a call to `bloechlCornerDOSWt` + sum. The subroutine
follows the same case structure (Cases 0-3) as
`bloechlCornerWeights` but returns the energy derivative
of each corner weight (the per-corner spectral density)
rather than the cumulative value.

The identity
`sum(cornerDOSWt_LAT) == dosContrib` provides a built-in
self-consistency check.

**Distinction from `bloechlCornerWeights`.** The
cumulative corner integration weights `cornerIntgWt_LAT`
(returned by `bloechlCornerWeights`) give the fraction of
the tetrahedron's occupation attributed to each corner up
to energy E. They are dimensionless and are the correct
quantity for integrated properties like
`electronPopulation_LAT` (section 1.5), which evaluates
them at a single energy (the Fermi level). The corner DOS
weights `cornerDOSWt_LAT` (returned by
`bloechlCornerDOSWt`) have units of 1/energy and are the
correct quantity for the energy-resolved DOS density.
Both subroutines share the same case structure and
intermediate variables; only the final expressions
differ.

**Diagnostic integral unit fix.** The trapezoidal integral
that computes "Spin States Calculated" (the integrated
area under the TDOS) multiplies the TDOS (in states/eV)
by deltaDOS (stored in Hartree). This produces
states/27.211 instead of states. The integral must include
a hartree conversion factor: `deltaDOS * hartree`. This
same fix applies to the integrated-area diagnostic in
both `computeTDOS_LAT` and `computeDOS`.

**BZ weight normalization.** The Gaussian DOS path uses
`kPointWeight` as its BZ integration weight. By
convention, `kPointWeight` sums to 2.0 (set as
`weightSum` in kpoints.f90). This factor of 2 accounts
for the two electron spin states per band in a
spin-unpolarized calculation: each band contributes
`sum(kPointWeight)/spin` = 2/1 = 2 spin states. The
tetrahedron BZ integration uses `tetraVol`, which sums
to 1.0 (just the geometric BZ fraction). To produce
DOS values on the same scale as the Gaussian path, the
LAT accumulation must include `sum(kPointWeight)` as a
multiplicative factor alongside `tetraVol`. Both
`computeTDOS_LAT` and `integratePDOS_LAT` require this
factor. For spin-polarized calculations (spin=2), the
factor becomes 2/2 = 1 spin state per band, which is
also correct.

#### 1.3.1 The Bloechl correction terms

Linear interpolation of eps(k) inside a tetrahedron misplaces
the iso-energy surface, and the resulting error falls off slowly
with mesh density. Bloechl's correction (PRB 49, 16223 (1994),
eq. 22) compensates it by adjusting each corner's cumulative
integration weight:

!! dw_i = (1/40) * D_T(E_F) * sum_{j=1..4} (eps_j - eps_i)
!!      = (1/10) * D_T(E_F) * (epsBar - eps_i)

where `D_T(E_F)` is THIS tetrahedron's density of states at the
Fermi level and `epsBar` is the mean of its four corner
eigenvalues. The second form follows from
`sum_j (eps_j - eps_i) = 4*(epsBar - eps_i)` and is the one to
implement: it makes the three properties below visible instead
of leaving them to be derived.

**It sums to zero over the four corners.**
`sum_i (epsBar - eps_i) = 0` identically, so the correction
moves weight BETWEEN corners and never changes a tetrahedron's
total. Three things follow, and they are what make this
correction safe to add to a working SCF path: the electron count
is unchanged, therefore the Fermi search is unchanged, therefore
the `numElectrons * spin/2` calibration of 1.6d cannot be broken
by adding it. What changes is WHICH k-points hold the
occupation -- hence the charge density and the band energy, and
nothing else.

**It vanishes wherever `D_T(E_F) = 0`.** Every wholly occupied
and wholly empty tetrahedron contributes nothing, so insulators
are untouched and only straddling tetrahedra move. The si_fd-3m
agreement that validated the uncorrected implementation
therefore survives unchanged, and is not evidence about the
correction either way.

**Its sign is predictable, which makes it testable.** `dw_i` is
positive for corners below the tetrahedron mean, so weight
shifts toward lower-energy corners and the band energy falls.
The uncorrected si_cmce ladder sits ABOVE its Gaussian
counterpart at every mesh and is still descending at 720
k-points (TODO C138(a)), so the correction pushes LAT in exactly
the direction that would close that gap. Whether it closes
enough of it is the measurement.

**Which consumers.** The two that evaluate cumulative
integration weights: `electronPopulation_LAT` for the
bond/effective-charge path (1.5) and the SCF occupation path
(1.6). NOT the DOS path (1.3): that uses `cornerDOSWt_LAT`, the
energy derivative dw/dE, and eq. 22 corrects w rather than
dw/dE. A DOS correction is a separate question this section does
not answer.

**A separate routine, not a change to `bloechlCornerWeights`.**
`bloechlCornerCorrection(energy, sortedEps, cornerCorrWt)`
returns the four `dw_i`; the two consumers above add them to the
weights they already fetch. Three reasons, in order of weight:
`bloechlCornerWeights` stays the literal transcription of the
paper's uncorrected weight expressions, so a reader can check it
against the reference without mentally subtracting a correction;
the correction gets its own built-in self-consistency check,
`sum(cornerCorrWt) == 0`, as strong as the existing
`sum(cornerDOSWt_LAT) == dosContrib` identity and testable
without reference to any other quantity; and the two can be
measured apart, which folding them together forecloses.

**The input it needs is already computed.** `D_T(E_F)` is
`sum(cornerDOSWt)` for that tetrahedron before the `tetraVol`
factor, and `latElectronCount` already calls
`bloechlCornerDOSWt` at the same energy in the same loop (it
needs the total as the derivative for its Newton step). The
correction routine takes the same `(energy, sortedEps)` pair as
its siblings and forms the quantity itself, so it stays a pure
function of one tetrahedron's corners and shares their
degenerate-corner guards.

**Always on, with no new option.** A switch would be an
option-contract change (6.2.10) and would reach the cache key
through makeinput's outputs, which is a large cost for a
validation convenience. It is not needed: the uncorrected ladder
is preserved under
`jobs/si_fingerprint/seed/ladders/linear-tetrahedral/`, so the
comparison is a re-run against a saved baseline rather than two
live configurations. Note that such a re-run needs `--force` or
cleared run directories -- the cache asks whether this is the
same calculation, not whether the engine has since improved
(6.2.5), so a corrected binary alone will not miss.

**Equation 22 is the whole of it.** An earlier reading of the
reference list here cited "eqs. 22-24" as the correction, which
invited the belief that two further terms were owed. They are
not: 23 and 24 quantify the discrepancy between the true Fermi
surface and the linearly interpolated polyhedral one, a
comparison the paper makes rather than a formula anything
computes. Nothing in this design or its implementation uses
them, and a later reader should not go looking for the missing
two-thirds of a correction that is already complete.

### 1.4 LAT PDOS (Energy-Resolved Partial DOS)

The corner DOS weights `cornerDOSWt_LAT` returned by
`bloechlCornerDOSWt` (section 1.3) determine how much of
each corner's partial-DOS projection to include at energy
E. These weights are occupation-independent: they
distribute spectral density among the four tetrahedron
corners at any energy, whether occupied or unoccupied.
They are the LAT replacement for the Gaussian broadening
+ `kPointWeight` mechanism used in the current PDOS path.

The PDOS contribution from tetrahedron T, band n is:

  dPDOS_alpha(E) = (V_T/V_BZ) * sum_{c=1..4}
      cornerDOSWt_LAT(c) * p_alpha(k_sigma(c), n)

where `cornerDOSWt_LAT(c)` is the Bloechl corner DOS
weight from `bloechlCornerDOSWt`, p_alpha(k, n) is the
Mulliken projection of channel alpha onto band n at
k-point k (`oneValeRealAccum` in the current code), and
sigma is the permutation that sorts eigenvalues. The
`cornerDOSWt_LAT` values depend on E, the four sorted
corner eigenvalues, and the corner index c -- they are
recomputed for every (tetrahedron, band, energy point)
combination and never stored as a module-level array.

**The fundamental constraint:** p_alpha is needed at all 4
corner k-points of each tetrahedron simultaneously. The
current code computes projections one k-point at a time and
immediately discards them. This forces a two-pass design:

**Pass 1 -- k-streaming (compute projections):**
  For each k-point, read eigenvectors from HDF5, compute
  Mulliken projections for all bands and channels, store in
  P(alpha, band, kpoint).

**Pass 2 -- tetrahedron integration:**
  For each band and tetrahedron, sort corner eigenvalues,
  compute `bloechlCornerDOSWt`, accumulate weighted
  projections into PDOS.

**Memory strategy (resolved D1):** Store projections
only at IBZ k-points: P(alpha, band, k_IBZ). When
assembling tetrahedron corners, look up the IBZ
representative via fullKPToIBZKPMap and apply atom
permutation on-the-fly to map the channel index:

  P(alpha, n, k_full) =
      P(permuted_alpha, n, fullKPToIBZKPMap(k_full))

where the channel permutation follows from atomPerm
and the operation stored in fullKPToIBZOpMap(k_full).
For mode 0 (per atom-type, per l-shell), R carries each
atom onto an atom of the same type, so the type-level sum
maps onto itself and no channel permutation is needed.
That closure is enforced by `buildAtomPerm` rather than
assumed from the meaning of a type; see section 2.3.

This reduces memory by the IBZ reduction factor
(typically 4-48x depending on point group symmetry):
- Moderate system (200 channels, 500 states, 500 IBZ
  kpts from 2000 full): ~0.4 GB
- Large system (5000 channels, 500 states, 1000 IBZ
  kpts from 5000 full): ~20 GB

The atom permutation infrastructure (atomPerm,
fullKPToIBZOpMap) is shared with the Q* and bond order
fix (section 2.4). One additional array is introduced:
invAtomPerm (section 2.4, item 4), which provides the
inverse mapping R^{-1}(A) for backward channel
permutation during tetrahedron corner assembly.

**Inverse atom permutation.** The channel permutation
for modes 1-2 requires R^{-1}(A) where R is the forward
operation (k_IBZ → k_full) stored in fullKPToIBZOpMap.
Since atomPerm stores R(A), we build invAtomPerm(R, B)
= A where atomPerm(R, A) = B, giving R^{-1}(B)
directly. It is built in O_AtomicSites alongside
atomPerm, with array shape (numPointOps, numAtomSites).

**Per-mode channel permutation rules:**

  Mode  Channel           Permutation rule
  ───────────────────────────────────────────────────
  0     per-type, per-l   None: type-level sum is
                          invariant under R (closure
                          enforced by buildAtomPerm)
  1     per-atom total    invAtomPerm(R, atomIdx)
  2     per-atom, per-l   invAtomPerm remaps atom;
                          l-shell offset unchanged
                          (a type shares one basis,
                          so the offset carries over)
  3     per-atom, per-lm  Not supported: requires
                          D^l(R) rotation matrices
  ───────────────────────────────────────────────────

For modes 1-2 a precomputed channelPermTable(R, alpha)
avoids repeated index decode/encode in the inner loop.
Mode 0 needs no table (identity mapping). Mode 2
decodes alpha into (atom, l-offset), permutes the atom
via invAtomPerm, and re-encodes using the permuted
atom's cumulative offset in cumulNumDOS.

**Mode 3 restriction.** When kPointIntgCode == 1 and
detailCodePDOS == 3, the program stops with a clear
error message. Individual Cartesian Gaussian projections
(px, py, pz separately) mix under rotation via D^l(R)
(section 2.3). Atom permutation alone does not suffice
and the full rotation matrices are not available.

**Refactored computeDOS structure.** computeDOS gains
an internal branch on kPointIntgCode inside the spin
loop. The setup phase (pdosIndex construction,
cumulNumDOS, allocations) and output phase (file
writing, normalization) are shared between Gaussian and
LAT. Only the computation phase -- filling
pdosComplete -- differs:

  computeDOS(inSCF):
    Setup (shared, unchanged)
    do h = 1, spin
      if kPointIntgCode == 1:
        LAT two-pass → fills pdosComplete
      else:
        Gaussian single-pass → fills pdosComplete
      endif
      Output (shared): normalization, file writing
    enddo

The Gaussian single-pass is unchanged, preserving its
memory efficiency (no projection array needed).

**imago.F90 dispatch.** The calling sequence becomes:

  if (kPointIntgCode == 1) then
    call computeTDOS_LAT
  endif
  call computeDOS(inSCF)

computeTDOS_LAT remains the validated eigenvalue-only
TDOS writer (fort.60/61). Inside computeDOS, the LAT
branch writes PDOS (fort.70/71) and localization index
(fort.80/81) but skips TDOS output (already written).
The Gaussian branch writes all three as before.

**Normalization.** The Gaussian path normalizes
pdosComplete by electronFactor (ratio of exact electron
count to Gaussian-broadened count) to correct for
broadening-tail truncation. For LAT, tetrahedron corner
weights provide exact BZ integration, so electronFactor
≈ 1.0. The LAT branch computes and logs this ratio as a
diagnostic (to fort.20) but does not apply it to
pdosComplete.

### 1.5 electronPopulation_LAT for Integrated Properties

For integrated (energy-summed) partial properties -- effective
charge, bond order -- the existing k-point loops barely need
to change if we precompute the LAT analog of
`electronPopulation`.

Define `electronPopulation_LAT(n, k, spin)`: the fractional
electron occupation of state (n, k) as determined by
tetrahedron integration. It answers the same physical
question as `electronPopulation` -- "how occupied is this
state, weighted for BZ integration?" -- but computed via the
LAT method instead of Gaussian broadening + Fermi filling.

  electronPopulation_LAT(n, k) =
      sum over all tetrahedra T containing corner k {
          (V_T / V_BZ) * w_c(E_Fermi)
      }

where w_c(E_Fermi) is the cumulative Bloechl corner
weight from `bloechlCornerWeights` (eqs. 18-21 evaluated
at the Fermi energy). Size: numStates x numKPoints x spin --
for 500 bands, 2000 kpts: ~8 MB.

**Naming convention:** The `_LAT` suffix identifies the
integration method. If future methods fill the same role
(a per-state occupation weight for BZ integration), they
follow the pattern `electronPopulation_XXX`.

**Usage:** In `computeBond` and the effective-charge branch
of `computeDOS`, replace:
  `electronPopulation(stateSpinKPointIndex)`
with:
  `electronPopulation_LAT(j, i, h)` (band j, kpoint i,
  spin h)

This unification means:
- Bond order loop: structure unchanged, swap weight source
- Effective charge loop: same
- Energy-resolved PDOS: `bloechlCornerDOSWt` (two-pass,
  occupation-independent; see section 1.4)
- TDOS: simplest case, eigenvalues only

The integration method (Gaussian vs. LAT) becomes a
pluggable parameter in a common projection-then-weight
framework.

**Distinction from `bloechlCornerDOSWt`.** The corner
DOS weights `cornerDOSWt_LAT` (sections 1.3-1.4) are
energy-resolved, occupation-independent, and transient
(recomputed per tetrahedron per energy point).
`electronPopulation_LAT` uses the cumulative corner
integration weights `cornerIntgWt_LAT` from
`bloechlCornerWeights`, evaluated once at the Fermi
energy; it is occupation-integrated, stored as a
module-level array, and has the same lifecycle as
`electronPopulation`.

**Resolved (D2):** Replacing `electronPopulation` with
`electronPopulation_LAT` provides the correct occupation
weight, but does not by itself cover bond order
accumulation correctly. The accumulation loop must also
apply atom permutation to distribute each IBZ k-point's
contribution correctly. See section 2.5 for the full
analysis.

### 1.6 LAT in the SCF Occupation Path

Sections 1.3-1.5 place LAT entirely behind the *post-SCF*
properties: TDOS reads eigenvalues, PDOS reads the
energy-resolved corner weights, and `electronPopulation_LAT`
serves bond order and effective charge. The SCF itself uses
none of them. `valeCharge` accumulates the valence charge
from `electronPopulation` -- the Gaussian-broadened,
thermally-filled array `populateStates` builds -- with no
branch on `kPointIntgCode`. A ground-state SCF therefore
returns a bit-identical total energy under either
integration code, which is confirmed empirically: the same
solid run at ten meshes under code 0 and code 1 agreed to
all printed digits.

**Why this matters.** The quantity a k-point convergence
study reads is the SCF total energy versus mesh. A metal's
energy oscillates there rather than settling, because the
occupation function has a step at the Fermi surface and a
uniform quadrature over a discontinuous integrand converges
slowly and non-monotonically. Tetrahedron integration is
the standard remedy that introduces no smearing parameter
-- the converged answer stays a single number rather than a
surface over (sigma, N_k) -- but it cannot help while the
SCF does not consume it. This section specifies that
consumption.

**What does NOT arise: the symmetry permutation.** Section
2.5 requires atom permutation for bond order because that
quantity is decomposed per site, which splits a symmetry
orbit into its members. The SCF's `potRho` is not: its
index runs over potential TYPES (`potDim` is the cumulative
alpha sum over `numPotTypes`), so each component is already
summed over a whole orbit of equivalent atoms. An orbit sum
is invariant under the operations that reduce the mesh, so
the weighted IBZ accumulation equals the full-BZ value
exactly, with no permutation applied or needed. LAT changes
only the scalar occupation `f(n,k)` multiplying the density
matrix and so cannot disturb this argument. The precondition
is that the reduction group is a genuine symmetry of the
Hamiltonian, which is section 2's invariant and is enforced
elsewhere.

**(a) The Fermi level is the substantive decision.** DESIGN
1.5 evaluates the corner weights "once at the Fermi energy",
which suffices for a property computed after the SCF has
finished. Inside the SCF the Fermi level is re-determined
each iteration from the current eigenvalues. Two choices:

  1. Keep `populateStates`' existing Gaussian/thermal search
     for `occupiedEnergy` and use LAT only to compute the
     weights at that energy.
  2. Determine the Fermi level from the tetrahedron
     integral itself, by requiring the LAT-integrated
     electron count to equal `numElectrons`.

**Choice 2 is required, not preferred.** Under choice 1 the
electron count implied by the LAT weights does not equal
`numElectrons`, because the two integration schemes place
the Fermi level differently for the same spectrum. The
resulting charge is wrong by a factor the SCF will partly
absorb into the potential, which makes the error hard to
see and easy to mistake for slow convergence. The
determination and the weights must come from one scheme.

The root find is cheap if it exploits what the LAT
machinery already provides. `N(E)` is monotone
non-decreasing, so bracketing is safe; `dN/dE` is the
density of states, which is exactly `cornerDOSWt_LAT`
(section 1.3 notes it is the energy derivative of
`cornerIntgWt_LAT`), so Newton converges in a few steps
rather than the thirty-odd a bisection to machine
precision would take; most (tetrahedron, band) pairs are
fully occupied or fully empty at a trial energy and
contribute a constant, so the working set after the first
pass is only the straddling pairs; and the previous
iteration's Fermi level is a good starting guess. The
existing Gaussian search bounded by `fermiSearchLimit` is
the structural model to follow.

**(b) Lifecycle: per iteration, not once.** The weights
depend on the Fermi level, which moves every cycle, so
`computeElectronPopulation_LAT` runs once per SCF
iteration rather than once per run. The cost is
`numTetrahedra x numStates` = `24 x numFullMeshKP x bands`
per evaluation, each item a four-corner sort and a few
polynomial evaluations -- still well under one percent of
the iteration's diagonalizations. The factor is 24 rather
than the 6 of a single-diagonal decomposition because
section 1.2 cuts every box four ways. Note it scales with the FULL
mesh while the diagonalizations scale with the IBZ, so the
margin narrows for a low-symmetry cell where the IBZ saves
little.

**(c) The substitution point is the unpack, not the
accumulation.** `valeCharge` reads `electronPopulation` as
a flat one-dimensional array walked in (kpoint, spin,
state) order -- explicitly documented there as not matching
the eigenvalue sort -- and unpacks it into
`structuredElectronPopulation(numStates, numKPoints, spin)`.
`electronPopulation_LAT` is already in that target shape,
so the LAT path replaces the unpack loop rather than
adding a reordering step. The accumulation below it is
unchanged. Routing the LAT array through the flat-index
loop is the one mistake available here and would scramble
the occupations silently.

**(d) The normalization convention must be stated, not
assumed.** `electronPopulation` has the k-point weight
folded in: `populateStates` multiplies by `kPointWeight`
when filling it. The LAT expression of section 1.5 carries
its own Brillouin-zone fraction, `(V_T / V_BZ) * w_c(E_F)`.
Both are therefore weight-inclusive, but by different
routes, and the two must be shown to agree in convention
before the substitution is trusted. A mismatch is a
constant factor on the valence charge, which an SCF partly
absorbs -- the same failure mode as (a), and as hard to
attribute.

**(e) Thermal smearing and LAT are alternatives, not
layers.** `thermalSigma` broadens the occupations in the
Gaussian path. Tetrahedron integration determines
occupations geometrically and needs no broadening, so a run
that sets both is asking for two different answers. The
LAT path ignores `thermalSigma` for the SCF occupation and
says so where the option is documented, rather than
silently applying one and discarding the other. (Note
`subroutine bond` already zeroes `thermalSigma` around its
own populate call for a related reason.)

**(f) A band-structure path is not a zone integral.**  A
symmetric band structure (`-sybd`, `-scfsybd`) replaces the
loaded k-point set with a 1-D path (`makePathKPoints`, equal
weights `2/N`), and its Fermi estimate is the equal-weight
fill of the path's own eigenvalues by `populateStandard` --
a 1-D sample along high-symmetry lines, useful for placing
the zero of the plot and nothing more.  Tetrahedron
integration has nothing to act on there: no uniform mesh, no
IBZ fold, no tetrahedra.  Yet the LAT request travels with
the k-point FILE (`kPointIntgCode`, read by `readKPoints`),
which is written per deck and does not know which job will
read it, so a deck prepared with `-pscfkpint 1` (or
`-scfkpint 1`) hands the SYBD job a code of 1.  Two consumers
downstream of the SYBD branch would then act on it: the
tetrahedra block at the end of `initializeKPoints`
(tetrahedra of the FILE's mesh counts, a grid nobody
computed) and `populateStates`, whose LAT branch reads the
IBZ map that only a mesh build allocates.

The rule: **on the SYBD path the LAT request is ignored**, and
it is ignored at ONE switch.  The SYBD branches of
`initializeKPoints` set `kPointIntgCode = 0` before any
consumer runs and write a note to the output saying so and
why; the tetrahedra block, `populateStates`, and every other
consumer keyed on the code then take the standard path
without a guard of their own.  This is a property of the
band-structure job, not of the k-point file, so it lives in
imago rather than in makeinput.  MTOP is different and is
NOT switched: its full, unreduced mesh receives identity IBZ
maps and genuine tetrahedra from `initializeKPointMesh(0)`,
so LAT on the MTOP mesh is coherent.

**Spin: one Fermi level, matching the Gaussian path.** That
path merges every (state, spin, kpoint) triplet into a
single sorted list and fills it in ascending energy until
the cumulative charge reaches `numElectrons`, with each
state receiving `kPointWeight / spin`. There is no
per-channel constraint anywhere, so the magnetic moment is
an OUTCOME of which channel's states sort lower rather than
an imposed input -- which is the physically right choice,
since both channels of a collinear system in equilibrium
share one chemical potential. LAT adopts it unchanged. Two
per-channel Fermi levels would be a physics change wearing
an integration-scheme costume, and would make LAT and
Gaussian runs of one system differ for a reason having
nothing to do with integration.

Non-collinear magnetism, when it comes, generalizes this
reasonably: spinor states replace the two-channel index but
one chemical potential and one sorted fill survive. What
does not survive is `spin` as a simple divisor and
`potRho`'s total/difference decomposition. That is a larger
change than the Fermi level and argues for building LAT
inside the present collinear structure rather than
anticipating it.

**Open questions.**

  1. Whether a run may switch integration scheme between
     the SCF and the post-SCF properties -- makeinput
     already exposes `-scfkpint` and `-pscfkpint`
     separately, so the input format permits it, but
     whether it is meaningful is not settled.
  2. The convergence-tolerance interaction: a tetrahedron
     Fermi level found to a loose tolerance introduces
     noise into the total energy that could masquerade as
     unconverged k-sampling. The root-find tolerance
     should be tied to the SCF convergence criterion
     rather than fixed.

---

### 1.7 Symmetry of Atom-Resolved Tetrahedron Results

Section 1.2 makes the decomposition point-group invariant for
lattices whose operations permute the mesh axes up to sign,
and leaves hexagonal and rhombohedral cells partly exposed.
This section closes the remainder by symmetrizing the result
rather than the geometry. The two are not alternatives: 1.2
earns the symmetry where it can, and this imposes it where
1.2 cannot reach.

**What is already immune, and must not be "repaired".** The
integrated partial properties -- effective charge and bond
order -- are unaffected by any of this, for a structural
reason worth stating so that a later reader does not add a
correction that would double-count.
`computeElectronPopulation_LAT` pools every tetrahedron
corner's weight onto the corner's IBZ representative, and
`computeBond` then spreads that pooled weight back across the
star, dividing by the star size. Pooling followed by even
redistribution IS the average over the star, so equivalent
atoms receive equal weight whatever the decomposition does.

The exposed quantities are precisely those that attach a
weight to an INDIVIDUAL full-mesh corner rather than pooling
it: the energy-resolved ones. That is `integratePDOS_LAT`
(section 1.4) and the two optical accumulators (section 12).

**The remedy: average the finished result over the group.**
For a quantity resolved onto atoms, replace each value by its
average over the point group, using the permutation table
that the unfolding already builds:

```
pdos_sym(alpha, E) = (1 / numPointOps)
    * sum over R of pdos(channelPermTbl(R, alpha), E)

poptc_sym(a, b, c, E) = (1 / numPointOps)
    * sum over R of poptc(partialPerm(R, a),
                          partialPerm(R, b), c, E)
```

No orbit has to be enumerated. Summing over every operation
of the group IS the average over each atom's orbit, because
the operations that carry an atom to a given orbit member are
a coset and contribute that member equally often. For the
same reason the DIRECTION of the permutation does not matter:
`channelPermTbl` is built from `invAtomPerm` and `partialPerm`
from `atomPerm`, and either gives the same average, since R
runs over the whole group in both cases.

**This is not a cosmetic patch, and the equivalence is worth
recording.** Averaging an atom-resolved result over its orbit
is exactly equal to replacing each integration weight by its
average over the star of its k-point. Writing p for a
projection and w for the weight, and using the transformation
law of section 2.3,

```
sum_k w_bar(k) p_A(k)
    = (1/|G|) sum_R sum_k w(Rk) p_A(k)
    = (1/|G|) sum_R sum_k w(k) p_{RA}(k)
    = (1/|G|) sum_R (unsymmetrized result for atom RA)
```

So this is the projection of the quadrature onto the
symmetric subspace, not an adjustment applied to a finished
number. Two consequences follow and both matter in practice.
Totals are preserved EXACTLY, because summing the averaged
weights over k gives back the original sum, so no total
spectrum and no TDOS moves. And the operation is exact for
every lattice, including the hexagonal and rhombohedral cases
1.2 cannot reach, because it averages over the star that
actually exists rather than over an assumed geometry.

What it does NOT do is remove the tiling bias itself. It
removes the part of the bias that breaks the symmetry. In a
cell where an atom is alone in its orbit there is nothing to
average and this does nothing at all, which is the reason 1.2
is worth its cost rather than being superseded here.

**Which quantities, and which modes.** Only the atom-grouped
decompositions need it; the type-grouped ones are already
invariant because every operation carries an atom onto an
atom of the same type, so a type-level sum maps onto itself
(section 2.5).

```
  quantity        applies to          skip for
  ------------------------------------------------------
  LAT PDOS        detailCodePDOS      mode 0 (per type);
                  1 and 2             mode 3 is refused
                                      on this path (1.4)
  LAT optical     detailCodePOPTC     codes 1 and 2
                  3 and 4             (type grouped)
```

The Gaussian pathway needs none of this. Its star average
distributes each IBZ point's contribution evenly over the
star by construction, so it is already symmetric -- where it
is applied at all, which is the separate defect recorded as
TODO C148.

**Default, and the opt-out.** Symmetrization is ON by
default. The reason is not tidiness: on a hexagonal or
rhombohedral cell it is the only thing standing between the
user and unequal values for atoms that must be equivalent,
and a default that ships those silently is the situation this
work exists to end.

A new tagged line in the k-point file turns it off, beside
`KPOINT_INTG_CODE` and `NUM_TETRA_DIAGONALS`:

```
SYMMETRIZE_LAT_PARTIALS
1                        ! 1 = on (default), 0 = off
```

The k-point file is the right home even though this reads
like an output setting, and the equivalence derived above is
why: averaging the result over the group IS averaging the
integration weights over the star, so this is a property of
the integration and belongs with the other two. Keeping all
three together also means a reader looking for "what did this
run do about symmetry" finds one place rather than two, and
it inherits the same per-phase independence for free. Adding
it is an input FORMAT change and therefore also touches
`makeinput.py` and the `skl/` examples.

The opt-out exists for diagnosis. Turning it off is how the
residual asymmetry gets measured, which is how the size of
C149 was established in the first place and how any later
claim about it must be checked.

**Report the spread, always.** Before averaging, compute the
largest deviation within each group of channels that the
averaging will merge, and write it to fort.20. An imposed
equality that leaves no trace is indistinguishable from an
earned one, and a reader who cannot see which they have is
being misled by a number that looks like a result. This also
gives every run a free measurement of the residual, which is
what makes the diagnostic switch above a rarely-needed tool
rather than the only route to the number.

**Seam inventory.** Every quantity the new code consumes,
where it comes from, who allocates it, and when it is
released -- including the release, which is the item whose
omission produced the defect described in PSEUDOCODE 19.2.1.

```
quantity           supplied by              lifetime
---------------------------------------------------------
numPointOps        read with the k-point    set in
                     file; O_KPoints          initializeKPoints;
                                              NOT set on the
                                              SYBD path
atomPerm,          buildAtomPerm and        allocated in
  invAtomPerm        buildInvAtomPerm,        setupSCF (SCF)
                     O_AtomicSites            and intgPSCF
                                              (PSCF); absent
                                              for k-point
                                              style code 0
channelPermTbl     buildChannelPermTable,   local to
                     from invAtomPerm         computeDOS; the
                                              symmetrization
                                              must run before
                                              it is released
                                              and before the
                                              output phase
pdosComplete       allocated in             filled by
                     computeDOS's setup       integratePDOS_LAT,
                                              then symmetrized,
                                              then written
partialPerm        buildPOPTCIndex, from    released by
                     atomPerm; built only     cleanUpPOPTCIndex
                     when detailCodePOPTC     from subroutine
                     >= 3 AND atomPerm is     optc, AFTER the
                     allocated                spectra are
                                              written
optcCondPOPTC      allocated in             filled by the
                     computeOptcSpectra,      accumulators, then
                     O_OptcSpectra            symmetrized, then
                                              printed by
                                              printOptcSpectra
symmetrizeLAT      readKPoints, from the    read once per
  Partials           k-point file;            k-point set, so
                     O_KPoints                the SCF and PSCF
                                              phases carry
                                              their own value
```

Two consequences of that table are binding on the code.
`atomPerm` is absent for k-point style code 0, so both
permutation tables may be unallocated; the symmetrization is
guarded on the table's allocation and says in fort.20 that it
was skipped, rather than failing or silently doing nothing.
And the optical symmetrization must sit between the
accumulators and `printOptcSpectra`, inside the window that
today's `cleanUpPOPTCIndex` placement keeps `partialPerm`
alive -- which is exactly the window that did not exist
before that placement was corrected.

**Memory.** The averaged values must be built from the
unsymmetrized ones, so the sum cannot be accumulated in
place. For the PDOS a single channel-indexed scratch vector
per energy point suffices. For the optical partials, work one
energy point at a time: the slab is `sumNumPartials` squared
by `dim3`, which is small, whereas a full copy of
`optcCondPOPTC` is not.

---

## 2. IBZ Correctness for Eigenvector-Dependent Quantities

### 2.1 Physical Basis for Symmetry in K-Space

Understanding why the IBZ reduction works for some quantities
and fails for others requires tracing the symmetry argument
from its origin in atomic geometry through to its consequences
for computed properties.

**From atomic symmetry to Hamiltonian symmetry.** A crystal's
point group is defined by the nuclear positions: a symmetry
operation R maps every nucleus onto another nucleus of the
same species. The electronic Hamiltonian -- kinetic energy,
electron-nuclear attraction, and electron-electron interaction
-- inherits this symmetry because R leaves the potential
energy landscape unchanged. Formally, [H, R] = 0, so the
eigenstates of H must transform as irreducible
representations of the point group. This is what forces
e_n(k) = e_n(Rk) for any point-group operation R.

The key point: although the quantity of interest is the
electronic wavefunction, we are relying on the atomic
geometry to set our expectations for electronic symmetry.
This reliance is rigorously justified (within the
Born-Oppenheimer approximation) because H is entirely
determined by the nuclear configuration. The one caveat is
symmetry-broken electronic phases (magnetic ordering, charge
ordering, Jahn-Teller distortions) where the electronic
ground state spontaneously adopts lower symmetry than the
nuclear geometry -- but these are special cases outside the
scope of a standard DFT/HF calculation.

**Basis closure under symmetry operations.** Does the choice
of an LCAO basis with Cartesian Gaussian atomic orbitals
undermine the symmetry argument? No -- provided the basis is
properly constructed.

The eigenvalues of H are basis-independent; the symmetry
relation e_n(k) = e_n(Rk) is a statement about the physics,
not the representation. The practical requirement is that the
basis set must be *closed* under the point-group operations.
For atom-centered functions this means: every symmetry-
equivalent atom (related by some operation R) must carry the
same set of basis functions. When R maps atom A to atom B, it
also maps the basis functions on A to the corresponding
functions on B. If the basis is identical on both atoms, the
full basis set is invariant under R, and the Hamiltonian
matrix respects the symmetry.

For Cartesian Gaussians there is one subtlety: the angular
parts (x^a * y^b * z^c) do not individually transform as
irreducible representations. A rotation mixes them. For
example, the six degree-2 Cartesian functions (xx, yy, zz,
xy, xz, yz) span a reducible representation (5 d-type + 1
s-type). Under R, a single Cartesian Gaussian maps to a
linear combination of Cartesian Gaussians of the same degree
on the rotated atom. As long as the complete set of
Cartesians for each degree is included on each equivalent
atom, the basis remains closed and the symmetry is preserved.

**Eigenvalues vs. eigenvectors.** The symmetry argument
guarantees that eigenvalues are invariant: e_n(k) = e_n(Rk).
Eigenvalues are scalar quantities unchanged by unitary
transformation, so the IBZ reduction is safe for any property
that depends only on eigenvalues -- total DOS, total energy,
band structure.

Eigenvectors, however, are *not* invariant. The eigenvectors
at Rk are related to those at k by a unitary transformation
that mixes the orbital expansion coefficients according to R.
They are different vectors in the basis representation. This
distinction is the root cause of the IBZ problem described in
the following subsections: any property that depends on the
wavefunction expansion coefficients (PDOS, bond order,
effective charge) cannot be correctly computed by simply
scaling the IBZ representative's contribution by the star
multiplicity.

**Time-reversal symmetry.** In addition to the crystallographic
point group, time-reversal symmetry (e_n(k) = e_n(-k)) often
doubles the effective symmetry, even when -k is not related
to k by any point-group operation alone. Most codes fold this
in when constructing the IBZ.

### 2.2 The Problem

Observed in KNbO3: bond orders for all O atoms are identical
(within machine precision) when using the full k-point mesh,
but differ for each O atom when using IBZ-reduced k-points.
This is incorrect -- the crystal symmetry requires them to be
identical.

### 2.3 Root Cause: Eigenvector Transformation Law

The root cause of the IBZ bug is that eigenvector-dependent
quantities do not simply scale with the star multiplicity.
This subsection develops the transformation rules needed to
determine which quantities are safe and which are not.

**LCAO wavefunctions and Mulliken projections.** In LCAO
the wavefunction for band n at k-point k is:

  psi_n(k) = sum_mu  c_mu(n,k) * phi_mu(k)

where mu indexes every orbital on every atom in the unit
cell. The coefficients c_mu(n,k) are the eigenvector
components (`valeVale` in Imago). The Mulliken population
of orbital mu in state (n,k) partitions the state's
electron count among basis functions:

  p_mu(n,k) = Re[ c_mu*(n,k)
                * sum_nu c_nu(n,k) * S_{mu,nu}(k) ]

where S is the overlap matrix. These satisfy the partition
sum_mu p_mu(n,k) = 1. Every property of interest --
effective charge Q*, bond order, PDOS -- is built from
sums of p_mu over specific index ranges (per atom, per
atom pair, per angular momentum shell). The Mulliken
overlap (bond order contribution) between orbitals mu on
atom A and nu on atom B is similarly:

  b_{mu,nu}(n,k) = Re[ c_mu*(n,k)
                      * c_nu(n,k) * S_{mu,nu}(k) ]

**How a symmetry operation transforms the eigenvector.**
A point-group operation R acts on three things at once:

  (a) Atoms: R maps atom A to atom R(A), preserving
      species and basis function types.
  (b) Orbitals: angular parts transform by the
      representation matrix D^l(R) for each l-shell.
      Example -- 90-degree rotation around z on p-orbs:

              [ 0  -1   0 ]
     D(R)  = [ 1   0   0 ]   px -> -py
              [ 0   0   1 ]   py -> px, pz -> pz

  (c) K-points: k maps to Rk in the Brillouin zone.

The complete eigenvector transformation is:

  C_n(Rk) = D(R) * C_n(k)                          (1)

where D(R) is block-diagonal: one block per atom, each
block being D^l(R) for that atom's orbital set. The
overlap matrix transforms consistently:

  S(Rk) = D(R) * S(k) * D(R)^dagger                (2)

**The Mulliken vector argument.** Define the Mulliken
vector M(n,k) = S(k) * C_n(k), so that the population
is p_mu(n,k) = Re[ c_mu*(n,k) * M_mu(n,k) ].

Under R, both the eigenvector and the Mulliken vector
transform by the same matrix:

  C_n(Rk) = D(R) * C_n(k)                  [from (1)]

  M(n,Rk) = S(Rk) * C_n(Rk)
           = [D(R) S(k) D(R)^dag] [D(R) C_n(k)]
           = D(R) * S(k) * C_n(k)
           = D(R) * M(n,k)                         (3)

The cancellation D(R)^dag D(R) = I in the middle step
is the key identity -- it ensures that eigenvectors and
their overlap-weighted counterparts transform in lockstep.

**Shell-sum invariance (proof).** The Mulliken population
of orbital mu at the symmetry-related point Rk is:

  p_mu(n,Rk) = Re[ sum_{a,b}  D(R)*_{mu,a} c_a*
                             * D(R)_{mu,b}  M_b    ]

(suppressing (n,k) arguments for brevity). Now sum over
all mu in a complete l-shell on atom R(A). These mu are
the images of the l-shell on atom A under D(R), so the
sum invokes the unitarity of D(R) within that subspace:

  sum_mu  D(R)*_{mu,a} * D(R)_{mu,b} = delta_{a,b}

The off-diagonal terms (a != b) vanish, giving:

  sum_{mu in l-shell on R(A)}  p_mu(n, Rk)

      = Re[ sum_a  c_a*(n,k) * M_a(n,k) ]

      = sum_{a in l-shell on A}  p_a(n, k)         (4)

The l-shell-summed Mulliken population on atom R(A) at
Rk equals the l-shell-summed population on atom A at k.
Atom totals obey the same relation by summing all shells.
Bond order between atom pairs follows by the same proof
applied to the two-atom coefficient product:

  B(R(A), R(B), n, Rk) = B(A, B, n, k)             (5)

In every case, no explicit rotation matrices are needed
-- only the atom relabeling A -> R(A).

**Why a type-level sum needs no relabeling at all.**
Equations (4) and (5) relate a quantity on atom A to the
same quantity on atom R(A). A sum taken over *every* atom
of one type needs one further property: that R carries
each atom of the type onto an atom of the same type. The
set being summed then maps onto itself, and the sum is
unchanged without any relabeling being applied.

It is tempting to justify that property by saying point
group operations permute atoms only within a species.
**That premise is not safe in Imago**, because a type is
a user-assigned grouping and need not be a symmetry
orbit. Two ordinary cases break it. In an amorphous cell
a type is a bin of locally similar atomic environments
and carries no symmetry content whatever. In a point
defect supercell, types are deliberately assigned from
the *pre-defect* symmetry -- to keep the number of
distinct potentials down and the cost with it -- even
though the defect has destroyed that symmetry.

What actually guarantees the property is a runtime check.
`buildAtomPerm` (O_AtomicSites) accepts a partner for
R(A) only among atoms sharing A's `atomTypeAssn`, and
stops with a fatal error when no such partner exists. So
whatever the types happen to mean physically, the
operations Imago actually reduces by are closed on each
type -- verified at startup rather than assumed.

The direction of a mis-typing decides whether it is
caught or is harmless:

- Typing **coarser** than the orbits (one type is a union
  of several orbits) is still closed, because a union of
  orbits of a group remains closed under any subgroup of
  it. This is where the defect supercell lands: types
  built from the original symmetry are unions of the
  reduced group's true orbits. Harmless, and it passes.
- Typing **finer** than the orbits (a single orbit split
  between two types) is not closed, and this is exactly
  the case `buildAtomPerm` rejects.

An amorphous cell reaches the same safety from the other
direction: with no symmetry to declare it carries only
the identity operation, so the permutation is trivial.

**The exception is style code 0.** With an explicit
k-point list Imago cannot build the symmetry maps, so
neither the unfolding nor the `buildAtomPerm` check runs
-- only the warning described in section 2.6. A hand
supplied, already-reduced k-point list is therefore the
one configuration in which a type partition finer than
the orbits can pass unnoticed and corrupt a type-level
decomposition.

**Why individual orbitals break the pattern.** For a
single orbital mu (not summed over a shell), the Mulliken
population at Rk is:

  p_mu(n,Rk) = Re[ sum_{a,b}  D(R)*_{mu,a} D(R)_{mu,b}
                             * c_a* * M_b              ]

Without a sum over mu, unitarity cannot be invoked and
the cross terms (a != b) persist. The individual-orbital
projection depends on the full D(R) matrix, not simply
on which atom mu belongs to.

Example: p-orbitals under 90-degree rotation around z.
Suppose the coefficients at k for one band are c_px=0.5,
c_py=0.3, c_pz=0.1. Applying D(R) gives at Rk:
c_px = -0.3, c_py = 0.5, c_pz = 0.1.

Simplified Mulliken projections (|c|^2):

  At  k:  |c_px|^2 = 0.25   |c_py|^2 = 0.09
  At Rk:  |c_px|^2 = 0.09   |c_py|^2 = 0.25

The px and py projections swap -- they are NOT related
by a simple atom permutation. But the p-shell sum is
0.35 at both k and Rk, preserved by ||c||^2 = ||Dc||^2.
This is why PDOS modes that sum over complete l-shells
(detailCodePDOS 0-2) work with atom permutation, while
individual-orbital PDOS (detailCodePDOS 3) does not.


### 2.4 The Atom Permutation Fix

From equation (4), the fix for any quantity that sums
Mulliken projections over a complete l-shell or over an
entire atom does not require rotating eigenvectors -- it
requires only knowing which atom maps to which under each
symmetry operation: the atom permutation table.

**Corrected accumulation.** For each IBZ k-point k_i,
instead of multiplying the projection by the star
multiplicity, loop over each operation R_s in the star
of k_i and accumulate into the permuted atom indices:

  Effective charge:
    chargeContrib(R_s(A)) += p_A(n,k_i) * f(n,k_i)

  Bond order:
    bondContrib(R_s(A), R_s(B)) += b(A,B,n,k_i)
                                 * f(n,k_i)

where f(n,k_i) is the occupation weight.

**Example 1: charge in a mirror-symmetric 1D chain.**
Two atoms per cell (A, B), related by mirror m that
swaps A with B. Each has one s-orbital. The IBZ contains
one k-point; the star is {k, mk} with size 2.

Eigenvectors at k for the bonding band:

  c_A = 0.8,  c_B = 0.6

At the mirror image mk, D(m) swaps coefficients:

  c_A = 0.6,  c_B = 0.8   (same eigenvalue)

Simplified Mulliken projections (|c|^2):

  At  k:   p_A = 0.64    p_B = 0.36
  At mk:   p_A = 0.36    p_B = 0.64

Note the atom permutation: p_A(mk) = p_B(k).

Naive IBZ weighting (WRONG -- this is the KNbO3 bug):

  charge(A) = 2 * 0.64 = 1.28
  charge(B) = 2 * 0.36 = 0.72    symmetry violated

Atom permutation (CORRECT):

  From identity: charge(A) += 0.64  charge(B) += 0.36
  From mirror:   charge(B) += 0.64  charge(A) += 0.36

  Result:  charge(A) = 1.00,  charge(B) = 1.00

**Example 2: bond order with C3 rotation symmetry.**
Three atoms A, B, C with 120-degree rotation symmetry.
R1 (120 deg): A->B, B->C, C->A. R2 (240 deg): A->C,
B->A, C->B. The IBZ k-point has star size 3.

Mulliken overlaps at k_1 for one occupied band:

  b(A,B) = 0.15    b(A,C) = 0.10    b(B,C) = 0.20

Naive star-weight multiplication (WRONG):

  BO(A,B) = 3*0.15 = 0.45
  BO(A,C) = 3*0.10 = 0.30
  BO(B,C) = 3*0.20 = 0.60          C3 violated

Atom permutation (CORRECT):

  R0: BO(A,B)+=0.15  BO(A,C)+=0.10  BO(B,C)+=0.20
  R1: BO(B,C)+=0.15  BO(B,A)+=0.10  BO(C,A)+=0.20
  R2: BO(C,A)+=0.15  BO(C,B)+=0.10  BO(A,B)+=0.20

Collecting (BO is symmetric):

  BO(A,B) = 0.15 + 0.10 + 0.20 = 0.45
  BO(B,C) = 0.20 + 0.15 + 0.10 = 0.45
  BO(A,C) = 0.10 + 0.20 + 0.15 = 0.45     C3 OK

**Required infrastructure:**

1. `atomPerm(numOps, numAtomSites)`: atom permutation
   table. atomPerm(R, A) = B means operation R maps atom
   A to atom B.
2. `fullKPToIBZOpMap(numFullMeshKP)`: for each full-
   mesh k-point i, the index of the symmetry operation
   R such that R(k_IBZ) = k_full(i) -- i.e., the
   operation that maps the IBZ representative to the
   full-mesh k-point. Currently `fullKPToIBZKPMap`
   gives only the IBZ index, not which operation did
   the mapping. For an IBZ k-point itself, the stored
   operation is the identity.
3. Star decomposition: for each IBZ k-point, the set
   of operations in its star. This falls naturally out
   of `fullKPToIBZOpMap` by collecting all full-mesh
   k-points that share the same IBZ representative.
4. `invAtomPerm(numOps, numAtomSites)`: inverse atom
   permutation. invAtomPerm(R, B) = A where
   atomPerm(R, A) = B, i.e., R^{-1}(B) = A. Built
   in O_AtomicSites alongside atomPerm. Required for
   LAT PDOS channel unfolding (section 1.4): when
   assembling projections at full-mesh corner k_f
   mapped from IBZ point k_i by forward operation R,
   the channel index must be transformed by R^{-1}
   to reference the stored IBZ-kpoint projection.


### 2.5 Per-Quantity Implications

The table below summarizes what each computed quantity
requires for correct IBZ unfolding, based on the shell-
sum invariance (4) and bond order invariance (5):

  Quantity                  Needed            Why
  ────────────────────────────────────────────────────
  TDOS (eigenvalues)        fullKPToIBZKPMap  e(Rk)=e(k)
  Q* (effective charge)     atom perm         eq. (4)
  Bond order                atom perm         eq. (5)
  PDOS mode 0 (type, l)    nothing extra      *
  PDOS mode 1 (atom total) atom perm         eq. (4)
  PDOS mode 2 (atom, l)    atom perm         eq. (4)
  PDOS mode 3 (atom, lm)   D^l(R) matrices   **
  POPTC codes 1, 2 (type)  nothing extra      * ***
  POPTC codes 3, 4 (atom)  atom perm twice   eq. (5) ***
  ────────────────────────────────────────────────────

*Mode 0 sums projections over all atoms of the same
type. Every operation carries each atom onto an atom of
the same type, so the set being summed maps onto itself
and the type-level sum is invariant -- no correction
needed. Note that this closure is enforced at startup by
`buildAtomPerm`, not inferred from a type being a
symmetry orbit; in an amorphous cell or a defect
supercell it is not one. See section 2.3.

**Mode 3 resolves individual Cartesian Gaussian
components (px, py, pz separately). Under rotation these
mix via D^l(R), so atom permutation alone does not
suffice. Correct unfolding requires the full
representation matrices for each l-shell. This is
deferred -- mode 3 is rarely used in practice.

***The partial optical properties (POPTC) decompose a
momentum matrix element between a *pair* of groups, so
each row above carries over from its PDOS analogue with
one change and one caveat.

The change: an operation acts on both group indices at
once. The pair matrix is conjugated by the permutation
rather than re-indexed along one axis,

    M(a,b) at k_full = M(invAtomPerm(R,a),
                         invAtomPerm(R,b)) at k_IBZ

which is the two-atom invariance (5) already used for
bond order, applied to a different quantity. No new table
is needed; invAtomPerm is simply used twice.

Which codes need it follows from the grouping alone, and the
numbering is arranged so that it can be read off the code
number. Section 11.3 orders the codes with grouping as the
major key: codes 1 and 2 group by TYPE, so the type-level
sum argument (*) covers both and neither needs a correction;
codes 3 and 4 group by ATOM and need the conjugation. The
test is a threshold, `detailCodePOPTC >= 3`, not a set of
cases.

The resolution axis does not enter, because both offered
resolutions -- whole group, and QN_nl summed over its m
components -- are complete-shell sums, which is exactly the
condition equation (4) needs. An nlm resolution would break
that and require D^l(R), which is one of the two reasons
section 11.2 gives for not offering it. The consequence is
worth stating plainly: **no offered POPTC decomposition
depends on the deferred D^l(R) matrices.**

The caveat, and it limits what codes 3 and 4 can claim: the
momentum operator is a vector, so an operation both
relabels the atoms and mixes the Cartesian components,
P_i(Rk) = sum_j R_ij P_j(k). The atom permutation handles
the first and not the second. The isotropic column that
`printSpectrumPOPTC` writes is the sum of the three
components divided by three -- a trace, invariant under
the rotation -- so for that column the atom permutation
is sufficient and the result is correct. The separate x,
y and z columns need the rotation as well and remain
unverified until it is added. This is the same defect
that affects the per-axis columns of the TOTAL spectra,
independently of any decomposition, so it is one fix for
both rather than something peculiar to POPTC.

Do not attempt to verify any of this with the identity
that the partials sum to the total. That identity holds
under ANY permutation of the two indices, because
permuting summands does not change their sum, and it is
therefore exactly blind to a missing or wrong unfolding.
A passing sum rule proves nothing here. Verify instead by
running a structure whose symmetry-equivalent atoms are
inequivalently oriented on a full mesh and on a reduced
mesh, and requiring the two to agree per atom.

One cost note, because the obvious implementation is
needlessly expensive. The star sum is a fixed linear map
on the pair matrix and does not depend on energy, and
broadening is linear. Symmetrize the pair matrix once per
k-point per transition, before broadening, and leave the
weighted accumulation over IBZ points as it stands.
Placing the star loop inside the energy loop instead
would multiply the innermost work -- a partials by
partials by three slab per transition per energy point,
with the partial count taken from the table in section
11.3 -- by the IBZ reduction factor of 4 to 48, for the
same answer.

**This correction belongs to the Gaussian integration
pathway.** Everything above assumes the spectrum is
accumulated by visiting irreducible k-points and weighting
them, which is what makes a star average necessary in the
first place. Section 12 adds a tetrahedron pathway that
visits full-mesh corners directly and applies the operation
once per corner as it fetches each quantity. On that
pathway the star average is not performed at all -- applying
both would count the symmetry twice. Section 12.6 gives the
argument. Both pathways are correct under reduction; they
arrive there differently.

**Resolution of D2.** The open question asked whether
replacing electronPopulation with electronPopulation_LAT
covers bond order accumulation correctly. The answer is
no. The occupation weight from electronPopulation_LAT is
correct (it properly sums tetrahedron contributions for
each IBZ k-point), but the Mulliken projection is
computed only at the IBZ representative. Multiplying a
correct weight by a non-permuted projection distributes
charge incorrectly among symmetry-equivalent atoms.

The fix has two parts that work together:

  (a) electronPopulation_LAT provides the correct
      occupation weight per (band, kpoint, spin).
  (b) The accumulation loop distributes each IBZ
      k-point's contribution across atom pairs using
      the atom permutation table (section 2.4).

Together, (a) and (b) give correct per-atom Q* and per-
pair bond order. Either alone is insufficient. Note that
total charge summed over all atoms IS correct with (a)
alone -- the error is only in the per-atom distribution.

### 2.6 Relation to LAT and Pragmatic Options

Both the LAT and Gaussian integration paths share the same
IBZ symmetry issue. LAT does not bypass it -- eigenvectors
still come only from IBZ k-points regardless of how the
occupation weights are computed. The distinction between
the two paths is purely in the weight calculation, not in
how projections must be unfolded.

**LAT-specific note on PDOS.** For LAT PDOS (section
1.4), the tetrahedron integration needs Mulliken
projections at all four corners of each tetrahedron
simultaneously. The corners are full-mesh k-point indices.
If we only diagonalize at IBZ k-points, each full-mesh
corner's projection must be "unfolded" from its IBZ
representative via atom permutation. Concretely, if
full-mesh corner k_f maps to IBZ point k_i via operation
R (stored in fullKPToIBZOpMap). Because R maps k_IBZ
to k_f (forward direction: R(k_IBZ) = k_f), the
inverse R^{-1} maps atoms at k_f back to those at
k_IBZ:

  p_{atom A, l-shell}(k_f, n) =
      p_{atom R^{-1}(A), l-shell}(k_i, n)

This is equation (4) applied to corner unfolding.

**Decided approach (Option A):** Use IBZ for SCF
diagonalization, unfold eigenvalues (and eventually
eigenvectors) to the full BZ for post-processing. The SCF
benefits from fewer diagonalizations, while tetrahedra
reference the full mesh for integration.

**Workflow scenarios:**

Imago reads separate kpoint files for SCF (`fort.15` =
`kp-scf.dat`) and PSCF (`fort.16` = `kp-pscf.dat`). Each
kpoint set goes through `readKPoints` + `initializeKPoints`
independently. The SCF and PSCF meshes can differ in
density, style code, and integration code.

Supported combinations involving LAT:
1. SCF with IBZ → DOS/bond with LAT (within SCF phase,
   same kpoints): The SCF diagonalizes at IBZ points. When
   `doDOS_SCF=1` and `kPointIntgCode=1`, the DOS routine
   unfolds eigenvalues to the full mesh via `fullKPToIBZKPMap`
   and uses tetrahedra for integration.
2. SCF with IBZ → PSCF with LAT (different, typically
   denser mesh): The PSCF reads its own kpoint file, builds
   its own mesh, IBZ reduces it, builds its own tetrahedra
   and `fullKPToIBZKPMap`. The PSCF diagonalizes at its IBZ
   points, then LAT DOS unfolds to its full mesh.
3. SCF with IBZ → PSCF with Gaussian (standard current
   behavior): No tetrahedra needed. IBZ kpoints and weights
   used directly.

The `fullKPToIBZKPMap`, `tetrahedra`, `numTetrahedra`, and
`tetraVol` are per-kpoint-set data stored in `O_KPoints`.
They are rebuilt each time `initializeKPoints` runs (once
for SCF, once for PSCF). When `kPointIntgCode == 0`
(Gaussian), they are not allocated.
`electronPopulation_LAT` lives in `O_Populate` alongside
its sibling `electronPopulation`; tetrahedra data is
passed in from O_KPoints as arguments.

**Implementation strategy:**

- The IBZ folding in `initializeKPointMesh` saves the
  full-to-IBZ mapping as
  `fullKPToIBZKPMap(numFullMeshKP)`. For each full-mesh
  index i, `fullKPToIBZKPMap(i)` gives the IBZ kpoint
  index whose eigenvalues are identical.
- `fullKPToIBZOpMap(numFullMeshKP)` (new): for each
  full-mesh index i, the index of the point group
  operation R such that R(k_IBZ) = k_full(i). The
  folding loop applies each operation to the IBZ
  representative and checks for matches among full-mesh
  k-points; the matching operation is stored. For an
  IBZ k-point itself, the identity operation is stored.
  Saved alongside `fullKPToIBZKPMap` during IBZ folding
  in `initializeKPointMesh`.
- `atomPerm(numPointOps, numAtomSites)` (new): for each
  point-group operation and atom, the index of the image
  atom. Built once during `initializeKPoints` from the
  point-group operations and atomic positions.
- `generateTetrahedra` uses `numAxialKPoints` to build
  tetrahedra referencing full-mesh indices (1 to
  `numFullMeshKP`). Called from `initializeKPoints`
  after the mesh is constructed.
- `tetraVol = 1 / numTetrahedra`, computed in
  `initializeKPoints` after tetrahedra are generated.
- For TDOS, eigenvalue lookup at full-mesh corner k is:
    `eigenValues(band, fullKPToIBZKPMap(k), spin)`
- For PDOS, projection at full-mesh corner k for an
  l-shell channel on atom A is:
    `p(R^{-1}(A), l, band, fullKPToIBZKPMap(k))`
  where R = operation(fullKPToIBZOpMap(k)) maps the
  IBZ representative to k (forward direction). See
  section 2.3 equation (4).
- **Style code 0 warning (resolved D3):** When
  `kPointStyleCode == 0` (explicit k-point list), Imago
  does not build the full mesh internally and therefore
  cannot construct `fullKPToIBZKPMap`, `fullKPToIBZOpMap`,
  or `atomPerm`. Emit a prominent warning at initialization
  that decomposition properties (effective charge, bond
  order, partial DOS, partial optical properties) will not
  be correct unless the user has taken extreme care to
  provide a symmetric k-point mesh.
  For style codes 1 and 2, Imago builds the full mesh
  and all symmetry maps internally, so the atom permutation
  fix works for both Gaussian and LAT integration paths.
  makeinput.py no longer produces style code 0 files;
  mesh mode (`-kp`) now writes style code 1.

**SYBD and MTOP paths bypass atomPerm.**  Symmetric band
structure (`-sybd`, `-scfsybd`) replaces the loaded k-point
set with a 1-D path between user-specified high-symmetry
vertices.  On that path every k-point is its own end product
-- band-structure output is per-k-point eigenvalues, and the
planned partial decomposition (future work) is a direct
per-atom projection at the very k-point being plotted, with
no star to unfold.  Modern-theory-of-polarization (`-mtop`,
`-scfmtop`) likewise replaces the loaded set: it builds its
own FULL, unreduced mesh from the `MTOP_INPUT_DATA` counts,
because the Berry phase is accumulated along strings of
k-points that must all be present -- there is no IBZ, and no
star to unfold.  In neither case are there shell-summed
quantities for atomPerm to reconstruct, so the table is not
needed.

This matters in practice because the SYBD branches of
`initializeKPoints` call `makePathKPoints`, and the MTOP
branches call `initializeKPointMesh(0)` on the MTOP counts,
and BOTH skip all of the point-ops setup that the style-code
0/1/2 branches do (`numPointOps` assignment, the point-op
array allocations, `computeRealPointOps`,
`computeAxisClasses`).  `abcRealPointOps` and
`abcRealFracTrans` therefore stay unallocated, as will
`xyzRealPointOps` when section 13.4 adds it, since it is
built on the same branches.
Consequently, the calls to `buildAtomPerm` and
`buildInvAtomPerm` in `setupSCF` (SCF path) and `intgPSCF`
(PSCF path) are guarded with `if ((doSYBD_SCF /= 1) .and.
(doMTOP_SCF /= 1))` and `if ((doSYBD_PSCF /= 1) .and.
(doMTOP_PSCF /= 1))` respectively.  For the same reason the
`RESOLVED_KP_*` records of section 3.9 / PSEUDOCODE 4d.5 are
written only when a style-code 1 or 2 branch actually built
and resolved the mesh; the SYBD and MTOP branches emit
none, exactly as style code 0 emits none, because the axis
classes those records report were never computed there.

**The HDF5 k-point group name is built in ONE place.**  The
SCF and PSCF HDF5 files each keep their k-point-dependent data
under a group whose 17-character name says which k-point set
produced it, so a re-run on a different set does not silently
read another set's matrices.  The name has three forms, one per
branch of `initializeKPoints`: a mesh gives
`nnnnn_nnnnn_nnnnn` (the three axial counts); an MTOP run gives
`nnnn_nnnn_nnnn_MP` (the MTOP counts, tagged); a band-structure
path gives `nnnn_` (the number of path k-points) followed by
the high-symmetry-point letters in path order, cut off at the
17th character.  Because the letters are appended one at a
time, they are placed by substring assignment; a formatted
internal write that names the destination among its own items
is not legal Fortran and overflows the record on the first
append.  Both `initHDF5_SCF` and `initHDF5_PSCF` obtain the
name from one O_KPoints routine, `kPointGroupName`, passing
their own SYBD/MTOP flags; the inputs it reads (`numPathKP`,
`numPaths`, `numHighSymKP`, `highSymKPChar`, `numAxialKPoints`)
are all loaded by `parseInput` (SYBD and MTOP control data are
read unconditionally) and, for the axial counts, set by
`initializeKPoints`, both of which run before either HDF5
initializer on their respective passes.  Nothing else -- no
script, no other Fortran -- reads a group name back by format,
so the two files share the form and cannot drift apart.

The downstream consumers of `atomPerm` and `invAtomPerm`
are `computeBond` (effective charge / bond order star
distribution) and the LAT PDOS channel-permutation step
in `dos.F90`.  Both are themselves gated by their own
`doBond_*` / `doDOS_*` flags, so a pure `-sybd`, `-scfsybd`,
`-mtop` or `-scfmtop` run with no decomposition flag never
reaches them.  Combining those with `-bond` or `-dos` is not
physically meaningful (you cannot integrate Q* or PDOS
over a 1-D path) and is left as an unguarded combination
for now -- if it occurs, those consumers will trip on the
unallocated `atomPerm` and fail loudly rather than emit
silent wrong answers.  Adding an explicit early refusal
for that combination is a worthwhile follow-up.

Note: `kPointWeight` is irrelevant for LAT integration --
the DOS contribution is determined entirely by the
analytic formula over each tetrahedron. The weight array
still matters for SCF electron counting (charge density
integration).


### 2.7 Conv-abc On-Disk Symmetry Operations with Cell-Mode Flag

`buildAtomPerm` (PSEUDOCODE 4) consumes point group
operations to permute atoms across the star of each IBZ
k-point.  The operations must act on atom positions in
the basis of the lattice currently held in `O_Lattice`
-- otherwise the matrix-vector product mixes
inconsistent quantities and no image atom matches.  The
skeleton's `full` / `prim` flag controls which lattice
ends up in `O_Lattice`: the conventional cell in `full`
mode, a primitive reduction in `prim` mode.  This
subsection pins down the basis convention used to
thread operations from the space-group database through
the kp file and into the runtime so the result is
correct for both modes and for every cell type (cubic,
hexagonal, monoclinic, ...) the code supports.

**The basis-mismatch issue.** Operations in
`share/spaceDB/<sg>` are stored in the natural axes of
the conventional crystallographic setting: rotation
matrix entries act on conventional-cell-abc fractional
vectors and per-operation translation vectors are
conventional-cell fractions.  When `apply_space_group()`
reduces a centered conventional cell to its primitive
form, the in-memory lattice in `O_Lattice` is
overwritten and is no longer the cell those operations
were written for.  Using the operations as-is against
primitive-cell atom positions mixes bases and produces
wrong images.  The same mismatch arises in `full` mode
for non-orthogonal-conventional cells (hex, trigonal-
hex setting, monoclinic, triclinic), where the
conventional-abc operation matrix differs from the form
needed to act on the loaded cell even when the loaded
cell happens to coincide with the conventional one.

**Two flavors of invariance.**  Two on-disk forms each
solve the cell-choice (full vs prim) part of the
problem:

- *Cartesian xyz Bohr* is *cell-choice-invariant*: a
  single similarity transform on the producer side
  (`M_conv * R_conv_abc * M_conv^{-1}` for rotations,
  `M_conv * t_conv_abc` for translations) converts the
  spaceDB operation to a Cartesian form that no longer
  references any particular cell choice.  Numerical
  values still depend on the lattice parameter,
  however, since the translation scales with `a_conv`:
  diamond's `(0.25, 0.25, 0.25)` d-glide becomes
  `(1.685..., 1.685..., 1.685...) Bohr` on disk, and a
  different Fd-3m material at a different lattice
  parameter writes a different number.
- *Conventional-abc fractional* is *also* cell-choice-
  invariant -- the spaceDB entries are already in this
  form -- and additionally *lattice-parameter-
  invariant*: every Fd-3m crystal stores the same
  `(0.25, 0.25, 0.25)` for its d-glide regardless of
  `a_conv`.  The space-group geometry is factored
  cleanly out of the material-specific geometry.

This design uses the conventional-abc fractional form.
Both options produce identical physics; the conv-abc
form keeps the on-disk values readable (they trace
directly back to `share/spaceDB/<sg>`), lets the
producer emit spaceDB entries with no transform, and
pushes the entire basis change onto the consumer where
it can be specialized cleanly for the two cell modes.

**On-disk format.**  The kp file carries the spaceDB
operations as-is alongside two small metadata blocks:

- `CONV_LATTICE`: a 3x3 block giving the conventional
  cell in Bohr (one row per lattice vector, columns
  are xyz).  Required because the consumer no longer
  carries implicit knowledge of which cell the
  operations live in.
- `CELL_MODE`: a string flag, `full` or `prim`,
  matching the skeleton's choice and therefore
  matching whatever ended up in `O_Lattice` at
  consumer time.

The `POINT_OPS` matrices and per-operation translation
vectors are written exactly as they appear in
`share/spaceDB/<sg>` -- no producer-side similarity
transform.  `makeinput.py`'s only added responsibility
is emitting the `CONV_LATTICE` block from
`sc.full_cell_real_lattice`, the conventional-cell
snapshot taken inside `apply_space_group()` before the
primitive reduction may overwrite the in-memory
lattice, and writing the `CELL_MODE` flag straight
from the skeleton.

**Consumer-side math.**  At runtime, imago reads the
spaceDB operations into `convAbcPointOps` and
`convAbcFracTrans`.  Each `POINT_OPS` file row is read
into the matching array *row*, so `convAbcPointOps`
holds the standard crystallographic rotation `R`
acting as `r' = R*r` -- not its transpose.  Let `L`
(`= realVectors`) and `Lc` be the loaded and
conventional lattices with their vectors as *columns*.
The `CONV_LATTICE` block stores `Lc`'s vectors as
*rows*, so `Lc = transpose(convLattice)`.  A Cartesian
point `x` has loaded fractional coordinates
`L^{-1} x` and conventional fractional coordinates
`Lc^{-1} x`, so the change of basis carrying loaded
fractional to conventional fractional is

```
T = Lc^{-1} * L          (r_conv = T * r_loaded)
```

formed once before the per-operation loop.

`computeRealPointOps` conjugates each operation into
the loaded direct-space basis.  Direct-space
fractional coordinates are *covariant*, so an
operation acting as `r' = R*r` transforms by the
ordinary similarity

```
R_loaded = T^{-1} * R_conv * T
t_loaded = T^{-1} * t_conv
```

producing `abcRealPointOps` and `abcRealFracTrans` in
the basis of whichever cell `O_Lattice` holds -- which
is exactly the basis `buildAtomPerm` uses for atom
positions (themselves built as `L^{-1} x`, i.e. a dot
with the *columns* of `invRealVectors`).  In `full`
mode the `CELL_MODE` flag tells the routine that
`L == Lc`, so `T = I` and the loop collapses to a copy
of the on-disk arrays.  In `prim` mode the flag
selects the full conjugation path.  No cell-shape,
centering, or symmetry assumption is built into the
math.

`computeRecipPointOps` then builds the reciprocal-
space twin.  k-point (reciprocal) fractional
coordinates are *contravariant*, the dual of the
covariant direct coordinates, so the reciprocal
representation of an operation is the **inverse
transpose** of the direct one:

```
R_recip = (R_loaded)^{-T}
```

It is computed directly from `abcRealPointOps` (so
`computeRealPointOps` runs first), and needs no lattice
matrices of its own.  kpoint folding via
`abcRecipPointOps` and atom permutation via
`abcRealPointOps` thus descend from a single source of
operations and a single change-of-basis step `T`.  No
full/prim branching exists anywhere outside the two
transform routines.

**Generality.**  `R_loaded` and `R_recip` come out as
integer matrices with determinant +/-1 that close the
point group -- the defining property of valid lattice
automorphisms -- for every cell type imago describes
(verified numerically for cubic, hexagonal, and
trigonal cells).  The earlier `C = M_loaded^{-1} *
M_conv = invRealVectors^T * M_conv` form is **not**
this `T`: it transposes the cell matrix incorrectly
and was correct only for cubic-like cells, whose abc
rotations are signed permutation matrices for which
the transpose error and the inverse-transpose
distinction both vanish.  The reciprocal side
compounded the problem by reusing `C` *and* assuming
the direct and reciprocal reps were identical.

**Diagnostic history.**  The bug surfaced when a
kaleidoscope SLURM flight ran graphite (SG 186,
hexagonal) and silica (SG 152, trigonal) alongside
silicon and diamond (cubic).  The cubic cells
converged; both 120-degree cells stopped at
`buildAtomPerm: no atom match found`.  For silica,
operation 2 (the 3-fold screw) sent Si site 1 at
`(0.465, 0, 1/3)` to `(0.633, 0.537, 2/3)` -- a
vacancy -- when the true image is the existing atom at
`(0, 0.465, 2/3)`.  Reproducing the reported value to
1e-9 pinned down two compounding transpose errors that
cancel only for cubic: the atom Cartesian-to-
fractional conversion dotted the *row* of
`invRealVectors` instead of the column, and the
rotation was applied as `R^T` because the operations
were stored column-major.  The fix un-flips the op
storage, replaces `C` with the correct `T` conjugation
(real) and its inverse transpose (reciprocal), and
dots the `invRealVectors` column in `buildAtomPerm`.
Validation: `R_loaded`/`R_recip` integer + group-
closed for all four structures, the atom permutation
closes for every operation and atom, and all four runs
converge through `buildAtomPerm` and the SCF.

**Input-precision gate.**  The conjugation is now exact,
but a *second*, independent precision hazard remains: a
fractional coordinate written as a short truncation of
a repeating decimal (e.g. `0.66667` for 2/3) is not
symmetry-exact, so a rotation maps it ~1e-5 off its true
image -- right at `buildAtomPerm`'s match tolerance.
This is an input-quality problem, not a math problem, so
it is caught at the source: `StructureControl.
_check_repeating_fraction` (called from `read_imago_skl`)
rejects any fractional coordinate that matches a third
or sixth (1/3, 2/3, 1/6, 5/6 -- the only non-terminating
decimals reachable by crystallographic 3- and 6-fold
axes) to fewer than 8 decimal places, naming the atom
and the suspected fraction.  Only thirds and sixths are
checked; denominators with no crystallographic basis
(7, 9, 11, ...) are excluded so a general-position
coordinate near, say, 2/9 is never falsely rejected.

---

## 3. K-Point Mesh: Density Input, Selection, and Reduction

### 3.1 Motivation

The existing workflow requires the user to specify explicit
axial k-point counts (e.g. `-kp 4 4 4`). This forces the user
to know the lattice geometry in advance and manually choose
counts that give adequate sampling density. A density-based
option (`-kpd D`) lets the user specify a single number -- the
minimum k-point volume density (kpoints per unit reciprocal-
space volume, Bohr^-3) -- and the program determines the per-
axis counts automatically. The total kpoint count is at least
`D * recipCellVolume`, distributed as uniformly as possible
across the three axes. This gives users a geometry-independent
knob: the same density value produces the same sampling
quality regardless of cell size.

This mechanism already exists in `imago` (`kPointStyleCode=2`,
`computeAxialKPoints`). The goal here is to expose it through
`makeinput.py` so that the normal user workflow supports it.

### 3.2 Pipeline When Using Density Mode

When the user passes `-kpd`, `-scfkpd`, or `-pscfkpd`:

1. `makeinput.py` extracts the point group operations from
   the space group database -- through the shared
   `symmetry.read_conv_abc_point_ops` reader (ARCHITECTURE
   2, 7), the single parser both the kp writer and the
   initial-potential producer read that file through -- and
   writes kpoint files (`kp-scf.dat` and/or `kp-pscf.dat`)
   directly, bypassing the `makeKPoints` executable.
2. Each file uses `kPointStyleCode=2` and contains:
   - `KPOINT_STYLE_CODE` = 2
   - `KPOINT_INTG_CODE` = 0 (default; histogram)
   - `MIN_KP_LINE_DENSITY` = the user's volume density
     (label is historical; the value is a volume density)
   - `KP_SHIFT_A_B_C` = the shift REQUEST -- the user's
     explicit shift, or an AUTO sentinel that imago
     resolves once the resolved counts are known (3.9).
     makeinput no longer selects a shift by crystal
     system, because the choice depends on the counts'
     parity, and in density mode the counts are resolved
     inside imago
   - `NUM_POINT_OPS` = number of point group operations
   - `POINT_OPS` = the 3x3 rotation matrices (abc coords)
3. At runtime, `imago` reads this file in `readKPoints`
   (including the point group operations and the shift
   request), then in `initializeKPoints`:
   - forms the reciprocal point operations and their axis
     classes (3.8), then the symmetry-compatible axial
     counts from the density (3.7)
   - resolves the shift: an AUTO request is selected from
     the counts and operations (3.9); an explicit request
     is honored, with any single-point axis zeroed (3.6)
   - builds the full Monkhorst-Pack mesh (3.9), then
     either folds it to the IBZ with reciprocal-lattice-
     periodic matching (3.10) or keeps the full mesh for
     tetrahedron integration
   - `convertKPointsToXYZ` — transform to Cartesian

### 3.3 Interaction with Existing Options

- `-kpd D` sets both SCF and PSCF density to D.
- `-scfkpd D` sets only the SCF density.
- `-pscfkpd D` sets only the PSCF density.
- Density mode is all-or-nothing: if any density option is
  given, both SCF and PSCF use density mode. An unset group
  defaults to a density of 1 (equivalent to a single
  k-point per axis). Explicit-mesh options (`-kp`, `-scfkp`,
  `-pscfkp`) are mutually exclusive with density options;
  if both are given, density wins and a warning is printed.
- `-kpshift` applies to both pipelines.
- `-printbz` is skipped when density mode is active. A note
  is printed informing the user that BZ visualization only
  works with an explicit k-point mesh.

### 3.4 File Format for Style Code 2

The file written by `makeinput.py` must match what `readKPoints`
in `kpoints.f90` expects. Specifically:

```
KPOINT_STYLE_CODE
2
KPOINT_INTG_CODE
0
MIN_KP_LINE_DENSITY
<density_value>
KP_SHIFT_A_B_C
<shift_a> <shift_b> <shift_c>
NUM_POINT_OPS
<n>
POINT_OPS
<3x3 matrix for op 1, one row per line>

<3x3 matrix for op 2, one row per line>
...
```

The point group operations are the rotational parts of the
space group symmetry operations with translations stripped.
They are given in fractional (abc) coordinates. Blank lines
between operations are optional (for readability) and are
skipped by the Fortran reader.

`<shift_a> <shift_b> <shift_c>` is the shift request. It is
either three explicit fractional offsets, or the AUTO
sentinel `-1 -1 -1`, which directs imago to select the shift
per 3.9 once the counts are resolved. A whole-Gamma request
is written as the explicit `0 0 0` (3.6), not the sentinel:
Gamma is a fully resolved single point at the origin, not an
auto-selection.

### 3.5 Impact on Memory Estimation and Summary

When using density mode, the exact number of k-points and the
axial mesh dimensions are not known at `makeinput` time (they
depend on the reciprocal lattice geometry, which is computed
inside `imago`). Two consequences:

- **Memory estimation:** Skipped entirely in density mode.
  To be revisited in the future.
- **Summary output:** The k-point count and mesh dimensions
  cannot be printed. The summary should indicate that
  density mode is active and print the density value instead
  of the count and mesh array.

### 3.6 Gamma-Point Requests and the Single-K-Point Shift Rule

A k-point shift is meaningful only on an axis sampled by more
than one k-point. On an axis with a single point, the "shift"
becomes that lone point's absolute coordinate, so a nonzero
shift silently moves the sample off the origin (Gamma). The
rule is therefore:

**A shift is applied only when more than one k-point is used.
The Gamma point -- a single k-point at the origin (0,0,0) with
no shift -- is requested with an explicit `0` sentinel, not
inferred from a `1 1 1` mesh.**

The sentinels, per group (SCF / post-SCF), are:

- Mesh mode: `-kp 0 0 0`, `-scfkp 0 0 0`, `-pscfkp 0 0 0`.
- Density mode: `-kpd 0`, `-scfkpd 0`, `-pscfkpd 0`.

Any of these marks that group as Gamma. A Gamma group is
written *canonically* -- always a single 1x1x1 style-code-1
mesh with shift `0 0 0`, regardless of which flag requested it
-- so there is exactly one on-disk representation of Gamma.
`imago.py`'s `check_gamma_kp` recognizes that form (1x1x1
mesh, zero shift) and selects the gamma-specialized real
executable `imagoG`, whose integral matrices are real (faster,
roughly half the memory). No change to `check_gamma_kp` is
needed: its existing style-code-1 detection already covers the
canonical Gamma file.

Two job families override that selection, and `init_exes` is
the ONE place the executable is decided.  A symmetric band
structure (`-sybd`, `-scfsybd`) and a modern-theory-of-
polarization run (`-mtop`, `-scfmtop`) always run on the
general executable `imago`, whatever the k-point files say,
because both replace the loaded k-point set inside imago
(section 2.6): the band path and the polarization strings
are sets of general k-points that only the complex arithmetic
can evaluate, and the polarization routine itself exists only
in the multi-k build.  The MTOP mesh, moreover, is written by
makeinput from the post-SCF mesh request, so on a deck whose
post-SCF group is Gamma the `MTOP_INPUT_DATA` counts are the
`0 0 0` sentinel -- a mesh with no strings to walk.  `imago.py`
refuses `-mtop` and `-scfmtop` on such a deck with a message
saying to regenerate it with a post-SCF mesh, and imago itself
stops on any non-positive MTOP count, so a hand-edited input
fails loudly rather than dividing by zero.

By contrast, a `1 1 1` mesh is a *general* single-point request,
not the Gamma sentinel: it is written with whatever shift was
requested (the AUTO sentinel by default, or `-kpshift`), and it
runs on the general complex executable `imago`. Where that lone
point lands is decided by imago at run time under the single-
point rule above: the shift on a single-point axis is dropped,
so with the mesh convention of 3.9 the point sits at the origin.
The two requests thus sample the same k-point but differ in
routing and in what the file promises -- the `0` sentinel is
the only form that guarantees Gamma on disk and selects the
real-arithmetic executable. (There is no "shifted single
point": a shift is meaningful only where an axis carries more
than one point, and the origin-plus-shift reading of a lone
point is exactly what the single-point rule forbids.) The
no-option default is a `1 1 1` request, so prior behavior is
unchanged; only the explicit `0` sentinel is new. A density of
1.0 stays a genuine density request, distinct from the `0`
Gamma sentinel.

A zero mixed with positive mesh counts (e.g. `0 0 1`) is
rejected as a fatal typo for the all-zero sentinel.

The rule lives in `makeinput.reconcile` (per-group resolution
of `kp_gamma` and the display `kp_note`) and `_make_kp` (the
canonical Gamma write). It supersedes an earlier behavior in
which a `1 1 1` mesh was labeled `(Gamma)` while a nonzero
auto-shift was still written -- a mismatch that routed the job
to the complex executable despite the label.

### 3.7 Mesh Selection: What the Density-to-Mesh Map Must Do

Sections 3.1 through 3.6 describe how a requested volume density
travels from the command line to the `imago` reader. This
section specifies the step that happens next and inside the
Fortran: turning that one density number into three integer
axial counts `(n_a, n_b, n_c)` and a shift. The density is the
knob the user and the convergence ladder (7.7, 7.8) turn; this
map is what the knob drives.

**The density itself is a sound rung parameter.** The quantity
that sets sampling quality is the full uniform mesh, of size
`n_a * n_b * n_c`, before any symmetry reduction. That product
is monotone in the requested density -- a larger density asks
for and receives a finer mesh -- which is exactly the geometry-
independent guarantee 3.1 promises. The convergence ladder may
therefore keep indexing its rungs by density. What must be
fixed is not the knob but the map: the rule that chooses the
integer mesh *shape* for a given density.

**Two requirements define a good mesh shape, and they do not
conflict.**

1. *Isotropy.* The spacing between adjacent k-points along each
   reciprocal axis is `|b_i| / n_i`, where `|b_i|` is the
   reciprocal lattice vector length. A good mesh makes these
   three spacings as nearly equal as the integer counts allow,
   so the Brillouin zone is sampled evenly in every direction
   rather than finely along one axis and coarsely along
   another.

2. *Symmetry compatibility.* The mesh should be carried onto
   itself by every operation of the cell's point group (3.8).
   A mesh that is not symmetry-compatible samples the zone in a
   way the crystal's own symmetry says is redundant on one side
   and sparse on another, and it reduces poorly because rotated
   mesh points miss their partners.

These two pull in the same direction, not against each other,
because reciprocal axes that a symmetry operation exchanges
necessarily have equal length. Isotropy then already wants
their counts equal, which is precisely what symmetry
compatibility requires. The isotropy target and the symmetry
constraint agree on every cell.

**The selection this section mandates.** Given the target
volume density `D` and the reciprocal cell:

1. Compute the reciprocal axis lengths `|b_i|` and the point
   group's axis classes (3.8).
2. Compute the continuous isotropic counts. With a common
   spacing `h = (prod|b_i| / (recipCellVolume * D))^(1/3)`, the
   real per-axis count is `x_i = |b_i| / h`. Axes in one class
   have equal `|b_i|` and so equal `x_i`; assign the class a
   single shared real count.
3. Round each class's shared count to the nearest positive
   integer. This is the isotropic choice, and because it is
   applied per class the result is symmetry-compatible by
   construction.
4. If the resulting full-mesh product falls below the density
   floor `D * recipCellVolume`, raise counts a whole class at a
   time -- never a single axis within a multi-axis class --
   choosing at each step the class whose increment most improves
   uniformity, until the floor is met.

The floor semantics of 3.1 are unchanged: the delivered mesh is
still at least `D * recipCellVolume` points. What changes is
that every candidate mesh along the way respects the symmetry
classes and stays as isotropic as the integers permit.

Steps 2 and 3 -- turning a target spacing into a shared,
per-class, rounded integer mesh -- are the reusable core of this
map. The crystalline mesh climb's opening floor (3.12.4) applies
exactly that core to a *fixed* spacing (the one that puts a small
cap of points on the densest axis) rather than one derived from a
density, and it skips step 4's density-floor raise, so the two
paths distribute and round a spacing identically and differ only
in the spacing they start from.

### 3.8 Symmetry-Compatible Axial Counts

A uniform mesh with counts `(n_a, n_b, n_c)` and a shift `s`,
built as in 3.2, is invariant under a point group operation `M`
(written in the reciprocal abc basis) precisely when `M` maps
every mesh point onto another mesh point. Writing `D` for the
diagonal matrix `diag(n_a, n_b, n_c)`, the mesh points are the
lattice `D^{-1} (Z^3 + s)` inside the cell (3.9), and `M`
carries this lattice onto itself exactly when

    D M D^{-1}   is an integer matrix,

and, for the shifted mesh, additionally when

    (D M D^{-1} - I) s   is an integer vector.

The first condition constrains the counts; the second
constrains the shift (3.9). The counts condition has a simple
combinatorial reading. For the signed-permutation matrices that
make up crystallographic point groups in the abc basis, a
nonzero off-diagonal entry `M[i, j]` means operation `M` sends
axis `j` onto axis `i`. `D M D^{-1}` is then integral in that
entry only if `n_i / n_j` is an integer; since the group also
contains the inverse operation, `n_j / n_i` must be integral
too, which forces `n_i = n_j`. Hence:

**Axis-class rule.** Define two axes to be *coupled* when some
group operation has a nonzero entry connecting them
(`M[i, j] != 0` or `M[j, i] != 0`, `i != j`). Take the
transitive closure to get axis *classes*. A mesh is symmetry-
compatible if and only if all axes in the same class share one
count.

The classes follow the crystal family. A cubic group couples
all three axes into one class, so `n_a = n_b = n_c`. A
hexagonal or tetragonal group couples the two in-plane axes and
leaves the principal axis free: two classes, `n_a = n_b` with
`n_c` independent. An orthorhombic group couples nothing: three
free classes. A triclinic group likewise leaves all three free,
its only operations being the identity and inversion, which are
diagonal. This is why the selection in 3.7 assigns counts by
class: the class structure is exactly the set of equalities the
group imposes, and it is derived from the operations themselves
rather than from a hard-coded table of crystal systems, so it
is correct for every cell `imago` can describe.

### 3.9 Shift Selection

The uniform mesh of 3.2 is a Monkhorst-Pack mesh: a regular grid
over the reciprocal cell, folded to the irreducible zone by the
point group with equal base weights. Its shift `s` (the
`KP_SHIFT_A_B_C` value of 3.4) is the Monkhorst-Pack grid offset,
expressed in units of the grid spacing. In fractional reciprocal
abc coordinates the mesh on axis `i` with count `n_i` is

    k_i = (m + s_i) / n_i,    m = 0, 1, ..., n_i - 1

so the mesh is the lattice `D^{-1} (Z^3 + s)` taken inside the
cell (3.8). Points are compared modulo a reciprocal lattice
vector throughout (3.10), so any point may equally be quoted in
`[-1/2, 1/2)`; the formula above fixes only which set of points
is meant, and it is the same set the classic Monkhorst-Pack
prescription produces. The shift has exactly two principled
values per axis, and their meaning does not depend on whether
the count is even or odd:

- `s = 0` places a sample on the origin: a Gamma-centered mesh,
  for every count. A single point (`n = 1`) with `s = 0` IS the
  Gamma point.
- `s = 1/2` centers the samples between grid nodes, so Gamma is
  absent: the Monkhorst-Pack off-Gamma mesh, for every count.

That parity independence is the reason the formula carries no
`-1/2` offset. A grid written `(m + s)/n - 1/2` describes the
same two families of point sets, but with the labels swapped
for odd counts (`s = 0` would then miss Gamma whenever `n` is
odd, and a lone point would sit at the zone corner), and every
statement in this section and in 3.6 that reasons about "Gamma
present" versus "Gamma absent" would be silently false for odd
counts. The classic Monkhorst-Pack prescription itself,
`(2m - n + 1)/(2n)`, is exactly `(m + 1/2)/n - 1/2` -- the
offset form with the half-shift built in -- and the standalone
`makeKPoints` program (`src/makeKPoints`) still builds that
grid and then subtracts its own shift from it. That program is
no longer on the makeinput path, and its grid convention is
not this one: here the shift is measured from the origin, and
the classic MP grid is the `s = 1/2` member of the family.

An intermediate offset (a quarter of the spacing, say) is neither
of these and has no place here; it displaces the samples so that
the point group no longer carries the mesh onto itself. Two
considerations select between the two principled values.

**Invariance.** A shift preserves symmetry only when it satisfies
the second condition of 3.8, `(D M D^{-1} - I) s` integer for
every operation `M`. A shift that violates it moves the samples
so the group no longer maps the mesh onto itself, which both
biases the sampling and defeats the reduction. The origin
`s = 0` satisfies the condition for every point group, since the
origin is a fixed point of every operation -- so a Gamma-centered
mesh is symmetry-compatible for any cell whatsoever. The half-
shift is admissible only for some groups: a cubic group admits
the body-half `(1/2, 1/2, 1/2)`; a hexagonal group forbids the
in-plane half-shift and admits only `(0, 0, 0)` and
`(0, 0, 1/2)`. The invariant set is group-dependent and must be
computed from the operations, not assumed from the crystal
family name.

**Reduction and convergence quality.** Among the invariant
shifts, a half-shift on an axis with an even count generally
folds the mesh furthest and tends to converge faster, because it
centers the samples symmetrically in the zone interior and keeps
them off the zone-boundary planes where fewer operations act. The
effect on reduction is large: a cubic `4x4x4` mesh folds to 4
irreducible points under `(1/2, 1/2, 1/2)` but to 20 under a
quarter-spacing offset that breaks the mesh's inversion symmetry
-- five times the eigenvalue problems for identical sampling.

**The rule.** Default to the Gamma-centered mesh `s = 0`. Because
it is invariant under every point group, it is the universally
correct choice and the required one wherever no half-shift is
admissible (face-centered cubic, the in-plane axes of hexagonal
and trigonal cells, and the other cases where the half-shift
fails the invariance test). Where the point group does admit a
half-shift on an even-count axis -- a body-centered cubic or
simple cubic cell, the principal axis of a hexagonal cell -- use
it, for its faster convergence and deeper reduction. In all
cases an axis carrying a single k-point takes no shift, per the
single-point rule of 3.6, since a shift would move that lone
sample off Gamma.

This selection is made where the counts are known. In density
mode the counts are resolved inside imago (3.7), so the
automatic shift is resolved there too, after the counts and
before the mesh is built; makeinput records only the request
(3.2, 3.4). An explicit user shift is carried through
unchanged; imago honors it verbatim and only warns if it is
not invariant, since an override is the user's prerogative.

This is a Monkhorst-Pack mesh throughout. It is distinct from the
special-point method, which builds a small, individually placed,
unequally weighted set of k-points rather than offsetting a
uniform grid; that method is not what this pipeline implements,
and its point coordinates are not mesh offsets. The Gamma-
centered default also ties to 3.10: an even Gamma-centered mesh
places samples on the zone-boundary planes, which are reciprocal-
lattice-periodic partners, so the periodic comparison of 3.10 is
what lets that default reduce fully.

### 3.10 IBZ Reduction with Reciprocal-Lattice Periodicity

Folding the uniform mesh to the irreducible zone (the loop of
3.2's `initializeKPointMesh`) applies each operation to each
mesh point and merges the image onto whatever mesh point it
coincides with. The merge is exact and loses nothing: when a
rotated point lands on a partner, `epsilon(M k) = epsilon(k)`
holds identically because `M` is a true symmetry of the crystal,
so transferring the partner's weight to the representative
leaves the Brillouin-zone integral unchanged. Every mesh point
is assigned to exactly one representative and contributes its
weight once, so the weights still sum to the full value.

**The periodicity gap.** Reciprocal space is periodic: two
k-points that differ by a reciprocal lattice vector are the same
physical point. A rotated mesh point may equal a partner only
*after* translation by such a vector -- in abc coordinates, when
the fractional difference is an integer rather than zero. A
merge test that compares raw fractional differences against a
small threshold misses these wrapped coincidences and leaves the
two points unmerged. The reduction is still *correct* -- a
missed merge only means the irreducible set carries an extra
point whose weight is accounted for -- but it is *incomplete*,
producing a larger IBZ and more eigenvalue problems than
symmetry allows.

**The periodic comparison.** The mesh must be reduced with a
reciprocal-lattice-periodic match: two points coincide when
their fractional difference is integral to within the threshold.
The basis-independent form compares each component of the
difference against its nearest integer,

    delta_i - nint(delta_i)   negligible for i = a, b, c,

which needs no assumption about the coordinate interval and is
correct for orthogonal and non-orthogonal cells alike. This is
strictly a superset of the raw comparison: it finds every merge
the raw test finds, plus the wrapped ones. It can therefore only
shrink the irreducible set, never enlarge it, and never creates
a false merge, since an integer fractional difference is a true
reciprocal-lattice translation.

**When it matters.** For orthogonal cells (cubic, tetragonal,
orthorhombic) the abc-basis operations are signed permutations
that preserve each coordinate's magnitude, and a mesh shifted
into the half-open interval never has a rotated image escape it,
so the wrapped and raw comparisons coincide and the periodicity
check is a no-op. For non-orthogonal cells (hexagonal, trigonal,
triclinic, and centered monoclinic) the operations mix
components, rotated images routinely fall outside the interval,
and the periodic comparison is what recovers the full reduction.
A hexagonal `6x6x6` mesh, for instance, folds to 28 irreducible
points with the periodic test but only 40 without -- a 43 %
inflation carried silently into every downstream cost. Because
`imago` is intended for all cell types, including the
hexagonal and trigonal cells the space-group machinery already
supports (2.7), the periodic comparison is required, not
optional.

**Postconditions.** Whatever the cell, the completed reduction
must satisfy three invariants, which serve as checkable
correctness conditions:

- The star sizes (the number of full-mesh points folded onto
  each representative) sum to the full-mesh point count.
- Each star size divides the point group order, since a star is
  an orbit under the group.
- The irreducible weights sum to the full weight value
  (`weightSum`, 3.2).

### 3.11 Density-Indexed Ladder and Duplicate Rungs

The convergence ladder (7.7) requests a sequence of increasing
densities and reads back the resulting energies. Because 3.7
delivers a mesh that is a *step function* of density -- a
well-behaved map holds the same integer mesh over a range of
densities, then jumps to the next -- two adjacent rungs will
sometimes resolve to the identical mesh. This is expected and
correct behavior of a good map, not a fault: distinct density
requests that fall in the same step are genuinely the same
calculation.

The consequence for the acceptance rule (7.8) is that two rungs
can produce a byte-identical energy, giving an energy delta of
exactly zero. The two-sided flatness test must not read such a
zero as physical convergence, because it is an artifact of the
map, not evidence that the energy has stopped moving. The guard
belongs in 7.8 and is small: before applying the flatness test,
recognize when consecutive rungs resolved to the same mesh --
comparing the resolved k-point sets, which 8.2 already surfaces
as `kpoint_count` -- and collapse such duplicates so the test
only ever compares genuinely distinct calculations. With the
mesh map of 3.7 producing smoothly varying shapes, true
duplicates become the *only* source of zero deltas, so this
guard has nothing to catch but them.

This keeps the ladder indexed by density, as 3.1 intends, while
making its rungs trustworthy. It is deliberately less invasive
than re-indexing the ladder by resolved mesh: the density knob
and the `makeinput.py` pipeline of 3.2 are preserved, and only
the density-to-mesh map (3.7) and the duplicate guard (7.8)
change. The guard stays valid, but for the convergence *search*
itself, 3.12 revisits this choice: once 3.7 makes the mesh coarse
for high-symmetry cells, a fixed density grid cannot resolve
convergence, and the search climbs meshes directly.

### 3.12 Adaptive Mesh Climb for Convergence

Section 3.11 keeps the convergence search on a *fixed* grid of
requested densities. That works when the density-to-mesh map is
fine-grained, but 3.7 makes the map coarse for a high-symmetry
cell: a cubic group forces `m x m x m`, so a wide spread of
densities collapses onto a handful of meshes, and a fixed grid
cannot place enough *distinct* rungs to see the energy flatten.
A cubic cell whose grid tops out at `[5,5,5]` may still be moving
by more than a meV per atom there, with no finer mesh available
below the grid ceiling. This section replaces the fixed grid
with an adaptive climb through meshes, and in doing so answers
the ladder-axis question 3.11 set aside.

#### 3.12.1 Density is the currency, the mesh is the step

Two roles were conflated in the density ladder, and separating
them resolves the tension:

- **Density is the currency.** It is cell-size-independent and
  therefore comparable across materials (3.1), which is exactly
  what the historical guidance dataspace needs as its key: the
  predictor stores and returns a `kpoint_density` (7.2, 7.6), so
  a converged density learned from one material can seed another.
- **The mesh is the step.** It is what imago actually integrates
  over, it is discrete, and for a high-symmetry cell only a few
  exist. It is the natural unit for *searching* for convergence,
  because each distinct mesh is a genuinely different calculation
  and consecutive meshes differ by the least the symmetry allows.

So the producer **predicts and records in density, but searches
in mesh space**: a predicted density seeds a starting mesh, the
search climbs meshes until the energy is flat, and the converged
result is recorded back as the density that reaches that mesh.
Each half uses the unit it is good at.

**There is a third role, and neither of the two above fills it:
the order rungs are COMPARED in.** The two-sided test asks
whether a rung's energy sits between its neighbours', which
presumes each later rung is a *better sampling* than the one
before. Neither unit delivers that. Density is monotone in the
mesh and still fails to order the energies; the mesh is monotone
and fails too. Measured on fcc at the standard shift, the
irreducible point count -- what the calculation actually
evaluates -- swings by a factor of two to three between
consecutive rungs:

```
odd  n:  19  44  85 146 231 344 489 670
even n:   8  16  29  47  72 104 145 195 256
```

so `[20,20,20]` (256 points) is a step *backwards* from
`[19,19,19]` (670). Both schemes show it, Gaussian worse than
LAT, so it is a property of the mesh sequence and not of the
integration (TODO D21). The comparison order therefore wants the
irreducible count, which is neither the currency nor the step:
the count is not comparable across materials, so it cannot serve
as the predictor key, but it is what orders a single climb's
rungs by the information each actually adds. Reading the parity
split as an fcc quirk misses the general point -- a low-symmetry
cell has no parity, and bcc's counts are already monotone.

#### 3.12.2 The rung rule and the stride

The climb reuses machinery 3.7 already defines. Selecting axial
counts for a density is a loop that repeatedly increments the
axis class whose bump most evens the three reciprocal-axis
spacings (`spacingSpread`), stopping when it reaches the density
floor. A convergence climb is the *same* increment with a
different stop condition -- keep going until the energy is flat
rather than until a density is met. One rung to the next:

> Increment the axis class that minimizes `spacingSpread` (3.7).

Because the class structure comes from the point group (3.8),
this is one rule for every crystal system, with no per-material
table. The sequences it produces (real reciprocal lattices):

```
cubic  Si:  [2,2,2] -> [3,3,3] -> [4,4,4] -> [5,5,5] -> [6,6,6]
hex  graphite:  [2,2,1] -> [3,3,1] -> [4,4,1] -> [4,4,2] -> [5,5,2]
ortho  Si:  [2,2,2] -> [2,2,3] -> [2,3,3] -> [2,3,4] -> [2,4,4]
```

A cubic cell moves in lockstep -- one degree of freedom, so the
steps are coarse, but each is the *finest* mesh the symmetry
permits (nothing exists between `[4,4,4]` and `[5,5,5]`).

*Finest is not the same as monotone, and only the first was ever
established.* The rule guarantees no mesh is skipped; it does not
guarantee that rung k+1 samples better than rung k, and for fcc
it does not (3.12.1). Every consumer of this ladder that treats
adjacency as refinement inherits that gap -- the two-sided test
above all. A


hexagonal cell climbs its in-plane pair several steps for each
principal-axis bump, holding the spacing even while the counts
stay anisotropic (`na == nb` is preserved because they are one
class). A triclinic or orthorhombic cell moves one axis at a
time -- the finest granularity, and the best convergence
resolution. Every rung is symmetry-compatible and distinct from
the last by construction, so the climb produces no duplicate
meshes and the duplicate guard of 3.11 becomes a safety net
rather than the mechanism.

The rung rule fixes the *direction* of a step; *how far* to
step is separate. The climb advances by a **stride** -- a
number of consecutive rung-rule increments -- and computes an
energy only at the stride's endpoints, skipping the meshes
between. A stride of one is the fine climb just described,
which computes every ladder position. A larger stride crosses
many ladder positions for the cost of a single calculation,
which is what lets the bracket search (3.12.3) cover an unknown
convergence distance cheaply; the skipped meshes are not wasted
work, they simply go uncomputed unless the refine step later
asks for one. Because a stride is only the rung rule applied
repeatedly, every mesh it lands on is still symmetry-compatible
and distinct -- the stride changes the spacing of the
*computed* points along the one ladder, never the ladder
itself.

The shift rides along for free: each rung's shift is chosen by
`selectShift` (3.9) from that rung's counts, so it tracks the
count parity as the climb proceeds (a cubic climb alternates the
body-half shift on even meshes with Gamma-centered on odd).

#### 3.12.3 The stop condition: bracket, then refine

A rung is the converged one when its per-atom energy is within
the k-point threshold of BOTH its neighbours on the ladder --
the same two-sided test as 7.8 step 3c, over consecutive
*distinct* meshes. This is the *authoritative* convergence
call; everything below only decides which rungs to compute so
the call can be made cheaply.

**The threshold MUST sit above the ladder's own rung-to-rung
scatter.** A bar beneath the noise cannot be cleared by
convergence, only by two coincidences in a row -- and a search
that demands persistence will therefore reject correct answers
and occasionally accept lucky ones. Measured scatter at the top
of these ladders, eV per atom: fcc Al 0.0008 under LAT and 0.0047
under unsmeared Gaussian, bcc Fe ~0.0015, the transition metals
~0.001. The threshold this design carried was 5e-4 -- below all
of them. That single mismatch produced every symptom the seed
runs showed: a Gaussian ladder declaring convergence on scatter,
ladders stopping one rung short of confirming what they had
found, and four of five elemental metals climbing to the ceiling.
The bar was therefore raised to 2e-3 eV per atom, which converged
all thirteen seed solids where 5e-4 converged two, while moving
the insulators only from `[12,12,12]` to `[10,10,10]`.

**The bar is now 1e-3 eV per atom, because the population it
judges has changed.** The floor principle above is unaltered and
is not a taste: a threshold must sit above the scatter of the
ladders it is applied to. What changed is *which ladders those
are*. Every scatter figure quoted above -- Al at 0.0047 under
unsmeared Gaussian, Fe at ~0.0015, the transition metals at
~0.001 -- was measured on a METAL. That is not incidental. A
metal's energy oscillates as the mesh crosses the Fermi surface,
which is the whole reason its ladder is noisy, and it is also
the reason a metal cannot converge in k-points at all. Those
ladders no longer reach the flatness test: the gap test stops a
metal on every search shape before any convergence work is done.
So 2e-3 was a floor set by ladders the test now never sees, and
"all thirteen converge" counted five solids that today stop for
a different reason entirely.

What remains for the threshold to judge is insulators, whose
ladders settle rather than oscillate. For the six ordinary
`si_fd-3m` seeds, 1e-3 and 2e-3 pick the identical mesh,
`[10,10,10]`, so the tightening costs them nothing measurable.

The reason to spend anything at all here is the **gap**, not the
energy. `gap_ev` is read off whichever rung the climb stopped
on, and it is a predictor key (7.6) that nothing downstream
re-converges -- unlike the potential, which the consuming SCF
re-converges by construction. A looser bar stops earlier and
therefore records a coarser-mesh gap. That is the defect carried
as D22, and while this threshold does not fix it, a bar no
looser than the ladders need stops widening it.

**The case to watch is `si_ia-3`,** a narrow-gap insulator that
sits closest to the metals in behaviour. The recorded evidence
disagrees with itself about which meshes it reaches at which
bar, so it must be re-measured live rather than read off either
table.

This suits the deliverable rather than merely rescuing it. The
initial-potential database wants a rough good starting point for
a later self-consistent calculation, and that calculation
re-converges what it is given; the guidance dataspace is advice a
curator weighs, not an answer. 5e-4 eV per atom is 0.5 meV per
atom -- a publication-grade bar applied to a starting guess, and
that is the sense in which the original number was wrong, quite
apart from sitting under the noise. The same reasoning already
justifies the metal short-circuit below; it was simply never
applied to the threshold. Note where it stops: the argument
licenses *not* being tighter than the work needs. It does not
license being looser, which is why the bar is set by measurement
from below rather than by how rough a starting guess may be.

Two cautions the seed runs earned. The threshold is a *harvest*
setting and is not part of the run-reuse cache key (6.2.5), so
re-judging a workspace at a new threshold costs nothing -- but it
also means a stored result carries no memory of the bar it was
judged against, which is why an entry records its
`metric_threshold` (7.2). And relaxing the threshold also loosens
the bracket, since that phase reads
`threshold * stride_flatness_multiple`: a looser bar changes
which rungs get *computed*, not merely how they are judged. Two
climbs at different thresholds are not the same ladder scored
twice. The test needs the candidate rung
plus one below and one above, so the search must place three
computed points around a candidate before it can declare
convergence there -- and one further point above for each extra
consecutive flat rung a persistent search demands (below).

The default climb does not walk the ladder rung by rung. It
**brackets, then refines**. The *bracket* phase strides from the
seed by growing steps -- one ladder position, then two, four,
eight (3.12.2) -- computing an energy at each stride's endpoint
until a stride comes back flat: its two endpoints within the
bracket phase's per-atom flatness threshold, which is
deliberately *looser* than the convergence test's (below). A flat
stride is a small energy change spread over *many* added
k-points, so it is strong evidence the energy has settled. And
because the energy is steep below convergence and flat above,
the converged rung lies in the LAST
*non-flat* stride interval -- the one across which the energy
went from moving to settled -- so that interval is the bracket,
not the first flat stride above it. The *refine* phase then
**fills the bracket from the bottom up** -- computing ladder
positions lowest-first -- and re-applies the two-sided test to
the growing consecutive block after each new rung, stopping at
the FIRST mesh that passes. Because the fill climbs from the
bottom and the test runs as it goes, that first pass is the
SMALLEST converged mesh, and the search computes only up to it
and the neighbours the test needs -- never the rungs above. A
convergence low in the bracket therefore costs only the fill
beneath it, not the whole interval: a cubic cell that settles at
its lower edge is confirmed without computing the wide meshes
near the top. Filling rather than bisecting is deliberate: the
bracket is small by construction, its width being the stride of
that last non-flat interval, which the max-stride cap keeps
short; and the two-sided test needs a short run of *consecutive*
rungs -- three to confirm a single settled rung, one more for
each extra rung the persistence rule demands (below) -- which
filling lowest-first supplies directly. Filling low-to-high is
no more calculations than a sparse bisection's three-point
probes, and it is simpler, finds the smallest converged rung by
construction, and reuses the stop test unchanged.
The flatness trace the material carries -- the ladder
the harvest re-judges -- is the *consecutive filled bracket*
around the converged rung, the rungs the two-sided test actually
compared; the sparse bracket endpoints are search scaffolding,
not part of it, so the harvest re-judges a clean consecutive
ladder exactly as a fixed grid's (7.8, 3.12.4).

A flat stride can lie: an oscillating energy -- the near-metal
case -- can dip and return across a stride and read flat by
coincidence. Refine catches it, because no rung in a falsely
bracketed interval passes the two-sided test, and the search
then resumes striding from the top of the bracket. The coarse
phase only *proposes* a bracket; the two-sided test *disposes*.
So the search keeps against oscillation exactly the robustness
the fine climb has by construction -- it never trusts a stride,
only a verified interior rung.

That *propose-versus-dispose* split is why the bracket phase can
afford a looser flatness threshold than the convergence test --
and why it should. The bracket threshold is a fixed multiple of
the convergence threshold (3.12.6). A stride whose endpoints have
settled to within that looser bound has *nearly* converged, so
bracketing there is right: it stops striding one geometric step
sooner and never computes the next, far larger endpoint the
strict threshold would have demanded -- the endpoint that
dominates the cost, since a stride can double the mesh and a mesh
costs about its point count. The refine still judges with the
strict threshold, so a stride that reads loosely flat but has not
truly converged is caught exactly as an oscillation is: no
interior rung passes, and the search resumes striding upward. The
looser bracket can therefore only ever move where the search
*looks*, never where it *converges*. And the trade is asymmetric
in its favour: an over-eager bracket wastes a few *cheap*
low-mesh fills before it resumes, while the overshoot it avoids
is a single *expensive* high-mesh calculation. The looser
threshold shaves that overshoot whether the convergence sits low
or high on the ladder, because it always brackets at an
equal-or-lower endpoint than the strict test would.

The fine **unit-step climb** is the degenerate stride-of-one
search: it skips the bracket and computes every ladder position,
paying the most calculations for the finest resolution and the
most conservative reading of a noisy energy. It is kept as an
explicit option (3.12.5) for a material a curator would rather
not have searched by stride.

How many flat interior rungs the refine phase requires **scales
with the prediction's confidence**. A confident prediction (7.6)
is an independent statement of where convergence lies, so a
single two-sided-flat rung near the predicted mesh is
corroboration from two directions and the search stops there. A
cold or bootstrap search has no such corroboration and is
establishing ground truth the seed database will be trusted on,
so it requires the flatness to *persist* -- a second consecutive
flat interior rung -- before stopping. The extra rung is cheap
insurance against a coincidental plateau. Because that test
judges that many consecutive interior rungs, the refine fills one
rung MORE than that above the bottom of the flat stride, not just
one. The extra rung matters because the bottom rung of the flat
stride need not itself be settled -- its own lower neighbour may
still be moving -- so the first rung the two-sided test can
confirm is often one step higher, and confirming the flatness
persists from there needs computed neighbours higher still.
Filling too few above would expose too few interior candidates,
and a persistent search could never confirm.
Confidence also sets where the climb is seeded (3.12.4): a warm
seed starts near the answer, so the first stride is already flat
and the bracket phase all but vanishes; a cold seed starts low,
and the geometric stride is what keeps the long bracket cheap.

Because a genuinely hard or ill-posed material could stride
without ever flattening, the climb carries a **ceiling** that
bounds the bracket phase, and two kinds compose. A fixed maximum
per-axis count is the always-present hard backstop -- it
prevents a runaway today and needs nothing else built. A cost
ceiling drawn from the resource dataspace (8) is the
operationally meaningful bound, since a count is a poor proxy for
effort (a large cell at a modest mesh can dwarf a small cell at a
fine one); it layers on once DESIGN 8 predicts run cost. The
climb stops at whichever bites first. When the ceiling is what
the *bracket* phase reaches -- the next stride would step past
it before any stride read flat -- the search does not give up at
the cap: it fills and refines the final interval, from the
highest computed endpoint up to the ceiling, so a convergence a
geometric stride jumped over just below the cap is still found.
Only a material with no converged rung even in that last
interval -- still visibly steep at the ceiling -- is reported
as non-converged (7.8 step 3d / 7.9), exactly as a fixed grid
that never flattened is today, and it carries the energy trace
(`grid_energies`) so a curator can read why it stopped. A
monotone, narrowing trace at the ceiling is a slow but genuine
insulator that simply wants a finer mesh -- raise the ceiling. An
*oscillating* trace should almost never reach this point, because
an oscillating energy is the signature of a metal, and a metal is
recognised and settled far sooner, by its gap, as the next
paragraphs describe.

Metals are what would otherwise be dragged to that ceiling, and
they are recognised directly, from the quantity that defines
them. A metal has no band gap; its energy oscillates as the mesh
crosses the Fermi surface and never settles, which is precisely
why chasing its convergence is futile. The gap is already
computed at every rung -- part of each rung's result -- so the
search reads it for free rather than inferring metallicity from a
side effect. This is the direct signal an earlier design lacked.
That design watched instead for a finer mesh that *raised* the
energy past a margin, treating a large upward step as the
fingerprint of a Fermi-surface oscillation. The proxy fails for
the common near-metal whose oscillation is small: its energy
wanders by more than the convergence threshold yet by far less
than any margin set safely above real convergence wobble, so it
is never flagged and the climb is dragged to the ceiling after
all. si_cmce did exactly this -- twenty-plus rungs, its energy
reversing by a few times the threshold, never once by the
fifty-fold margin the proxy needed. Reading the gap sidesteps the
proxy entirely: a vanishing gap *is* metallicity, not a symptom
of it to be read out of the energy's behaviour.

So the stop logic gains a single rule, checked as each rung is
computed: the first rung whose gap is essentially zero -- below a
small threshold (3.12.6) -- declares the material a metal, and
the search stops. It settles at that rung, or one step above for
a slightly denser sampling, and records the result as a metal: a
deliberately rough starting potential. The test is *live*, not a
scan of a fixed mesh. A near-metal can show a small but non-zero
gap on a coarse mesh and have it close only as the mesh refines,
so the trigger is the first rung that actually reads zero,
wherever on the climb that falls -- a material that opens at a
gap of, say, 0.15 eV and reads zero two rungs later is declared a
metal there, and settles there.

**Which rung, when more than one is in hand.** A confident opening
grid resolves several rungs at once, and a refine fill lands rungs
*below* ones already computed, so at the moment the test fires the
gapless rung need not be the densest mesh on the ladder. It is the
coarsest gapless rung that is settled on. The two readings
coincide on a plain upward walk, which is how the distinction
comes to be lost. They part exactly where it costs something:
settling on the densest mesh can settle on a rung that read a
*gap*, and anything that later re-reads that single rung -- the
harvest does (7.8) -- then sees an insulator where the climb saw a
metal, and records a k-point convergence claim that was never
made. The coarsest gapless rung is also the cheapest, and
roughness is the whole intent. The justification is the
deliverable itself (6.2, 3.12). The initial-potential database
wants a *rough* good starting point for a later self-consistent
calculation, not a converged energy. A metal's energy cannot be
converged in k-points by any mesh worth paying for -- that is the
whole difficulty -- but its potential at a modest, floor-level
mesh is a perfectly serviceable starting guess, and far better
than the isolated-atom potential it replaces. Refining a metal's
mesh beyond that buys nothing the deliverable needs.

The gap test leans on the opening floor (3.12.4), for the same
reason the retired proxy did but more simply. The single Gamma
point (`[1,1,1]`) samples the whole zone at one place and can
misreport a gap in either direction; a rung that coarse is not a
mesh whose gap can be believed. The crystalline climb never opens
there -- its first rung is floored at a physically meaningful
resolution (3.12.4) -- so every gap the test reads comes from a
mesh already dense enough to trust, and no coarse-mesh guard is
needed to keep a Gamma-point artifact from being misread.
Non-crystalline system_types seed at or near Gamma by convention
(7.9) and converge without a bracket climb, so this metal
short-circuit is a crystalline-path concern only. An insulator
never triggers it: at floor-and-above meshes its gap reads
clearly non-zero and stays there, so it converges by the
two-sided test above, its gap never nearing the trigger. The seed
solids bear this out -- the one metal reads essentially zero
throughout, while every insulator, down to the narrowest-gapped,
sits several times the threshold above it.

Like the bracket phase's looser flatness threshold, the metal
test is a heuristic with a deliberately cheap failure mode. Its
one dial is the gap threshold (3.12.6), set low enough that no
real insulator crosses it yet high enough to catch a true metal's
essentially-zero reading.

*The cheap-failure claim was one-sided, and the seed runs found
the other side.* This design argued that a wrong firing costs at
worst a re-run, "never a wrong recorded energy". That covers the
test firing when it should not. It does not cover the test
FAILING to fire, which is what happened: fcc Al and fcc Cu were
recorded as gapped insulators, with `gap_ev` of 0.124 and 0.185
eV against a 0.05 eV cutoff and `gap_kind` of "indirect".
Aluminium is the textbook free-electron metal. The readings are
the finite-mesh artifact of 1.6 -- a level spacing in the
globally sorted spectrum rather than a gap at the Fermi level --
and the cost is not a re-run: `gap_ev` is a predictor key (7.6),
so such an entry teaches the dataspace a gap that does not exist.
The root is that the test reads the gap of ONE rung, the rung the
climb stopped on, which makes the classification depend on the
convergence threshold. That is a defect in what an entry
*measures* rather than in when the climb stops, and it is carried
as TODO D22 together with the matching degradation of insulator
gaps. Until it is resolved, a metal reaching the dataspace as an
insulator is a live possibility, not a hypothetical.

And it
*retires* machinery rather than adding it -- with metallicity
read straight from the gap, the rising-stride proxy and the
coarse-mesh guard it needed are both gone, and the search wants
no separate oscillation or stall test on top. It branches once,
on the gap, and gives each material the treatment that suits it:
the insulator its two-sided convergence, the metal a prompt,
rough, floor-level stop.

**The test belongs to every search shape**, and the branch is
taken before the shape is chosen rather than inside any one of
them. This design once scoped it to the automatic bracket-refine
climb, on the reasoning that the fine unit-step walk is the
conservative shape a curator pins on purpose and should be left
to compute every rung. That confused two different things. A
*stopping rule* is properly a matter of search shape, because the
shapes exist precisely to disagree about which rungs are worth
computing. A *classification* is not. Recognising that a material
is a metal is the second kind, and a shape that declines to make
the classification is not being conservative. It is blind, and it
is blind in the direction that costs the most.

Follow what the unit-step climb actually does with a metal when
it cannot make the call. It walks every rung to the ceiling
looking for a flat interior. A metal has none to find -- that is
the whole point of the classification -- so it exhausts the
ladder and stops non-converged. And a non-converged stop harvests
*nothing*: no potential is extracted, no entry is staged. So
withholding the test does not buy the same rough potential more
slowly. It converts a serviceable rough potential into no
potential at all, after paying for the entire ladder. It also
leaves the harvest's single-rung gap reading as the only metal
judgment anywhere in the system, which is exactly the reading
that recorded Al and Cu as insulators above.

The diagnostic use the old scoping was protecting survives, by an
existing dial rather than a special case. A band gap cannot be
negative, so a `metal_gap_threshold` (3.12.6) set below zero is a
test no rung can ever trigger. A curator who wants every rung of
a known metal computed -- to read the ladder itself, as the
seed-run gap ladders did -- sets it negative and gets the
walk-to-the-ceiling behaviour back, with the reason written down
in the manifest rather than implied by which climb shape was
chosen.

#### 3.12.4 Seeding the climb from a prediction

The guidance predictor (7.6) returns a `kpoint_density` for the
new structure (or, when under-trained, signals the bootstrap
path). `selectAxialCounts` (3.7) converts that density to a
starting mesh on the ladder, and the climb begins a rung or more
below it -- so the predicted rung has a lower neighbour and the
climb only has to move upward to acquire its upper one. How far
below to begin scales with confidence (below), so it is not a
fixed offset. The effect of the dataspace is to shorten the
climb, not to change its path:

- **Cold** (under-trained, or a novel material): start low and
  climb far. This is the bootstrap regime, where a wide search
  is correct because nothing is known about where convergence
  lies.
- **Warm** (a populated dataspace): start near the answer and
  confirm a flat interior in a rung or two. In the limit the
  climb is just "run the predicted mesh and its neighbours,
  confirm flat."

For crystalline system_types the cold start also carries a
**floor**. No crystalline material reaches the k-point accuracy the
harvest demands on a mesh as coarse as a single Gamma point or a
handful of points per axis, so beginning the climb there only
spends rungs in a regime the search must leave anyway -- and worse,
those coarsest meshes report their gap unreliably enough (3.12.3)
to mislead the metal test. So the crystalline climb opens no lower
than a
floor rung, defined as a **cap of a few points per axis**: the
densest reciprocal axis (the largest `|b_i|`) gets the cap count,
and every other axis, being coarser in reciprocal space, is sampled
to that same k-spacing and so gets fewer -- never more -- down to a
single point. A cubic cell floors at `[4,4,4]`; an anisotropic one
at `[4,4,2]`, `[4,3,2]`, `[4,1,1]`, and the like, capped at the
count on every axis. The cap is a per-axis maximum, not a total
point count, so unlike a fixed-density floor it never forces extra
k-points onto a strongly anisotropic cell whose long real axis
needs only a point or two. The opening rung is then the higher of
the seed-derived rung and this floor, so a confident warm seed --
already above the floor for any real convergence -- is untouched,
and only the cold bootstrap is lifted out of the coarse regime.
Non-crystalline system_types are exempt: a molecule converges at
Gamma-only (7.9) and must not be floored up off it. This floor is
what lets the metal test read a trustworthy gap at every rung
(3.12.3): the climb no longer visits a mesh too coarse to judge a
gap on, so no separate coarse-mesh gate is needed.

`predictor_confidence` (7.6) is the natural dial for how much
confirmation to demand, how far below the prediction to begin,
and which dispatch mode to use (3.12.5) -- a confident
prediction warrants a short, tight search; a weak one warrants a
wider one. Whatever the search finds, the converged rung is
recorded back as a density (7.8), so the dataspace stays keyed on
the transferable currency and every converged material sharpens
the next prediction. The density a converged mesh represents is
its full-mesh volume density -- the product of its axial counts
divided by the reciprocal cell volume -- which is self-consistent
with 3.7, so a future prediction of that density reproduces that
mesh in that cell. Its resolved mesh is recorded alongside the
density (7.2), so the exact calculation is auditable -- a density
round-trips through 3.7 only up to the rounding the map applies,
whereas the mesh plus the cell is exact.

#### 3.12.5 Orchestration

A fixed grid is dispatched once and harvested after (6.2, 7.7);
an adaptive climb is a control loop -- submit a rung, wait,
judge, decide the next -- so the producer gains an outer loop
around the dispatch layer rather than a single fan-out. This is
exactly the shape Principle 12 reserves for client Python: the
dependent, per-unit iteration lives in the producer, not in the
dispatcher. (6.2.3 frames the alternative -- folding the
iteration inside one unit's wingbeat -- as the default for
adaptive convergence; this climb deliberately takes the
producer-loop road instead, for the reason just given.) The
dispatch core stays domain-ignorant (Principles 9 and 12): it
runs whatever units it is handed and reports outcomes; the climb
logic that reads energies and chooses the next mesh lives in the
producer.

The parallelism inverts cleanly. Within one material the rungs
are serial -- rung N+1's mesh is not known until rung N's energy
is judged -- but across materials the climbs are independent, and
that independence is total: a chain's next rung depends only on
its own last energy, never on any other chain's. So no chain need
ever wait for another. The producer keeps one rung of every
active chain in the air at once, waits for whichever rung lands
first, judges that one chain, and -- the moment its energy is in
hand -- either sends that chain's next rung or retires the chain
because it converged or hit the ceiling. A chain that finishes a
rung early climbs on immediately rather than idling until some
slower chain's rung completes. Cluster throughput is preserved --
many chains climb at once -- and no chain is ever paced by an
unrelated one. The confident mode lays down its whole small grid
at once (below); those grid rungs simply enter the in-flight set
together and are judged as a group once all have landed.

This is a heavier interaction with the dispatch layer than the
one-shot flight of 6.2: the producer sends rungs and collects
them one at a time rather than as a single fan-out, so the
dispatch core exposes its send-off and its collect as separately
callable steps (6.2.3). The producer/dispatch wiring for the
climb is elaborated where predict-then-verify is specified (7.7).

As a chain retires, the workers that were running its rungs fall
idle, and because retirement is final -- a converged or
ceiling-stopped chain never re-enters the active set -- the pool
of idle workers only ever grows over a climb's life. Handing
those freed workers to the chains still climbing, so a late,
expensive chain finishes sooner, is a forward extension: it waits
on a parallel imago and is designed, but not built, in 6.2.11.

The serial dependency is worth paying only when the search is
genuinely uncertain, so `predictor_confidence` (7.6) gates two
dispatch modes. A **confident** prediction has already narrowed
convergence to a small mesh neighbourhood, so the producer
dispatches a small *fixed* mesh grid -- the predicted rung and a
rung or two on each side -- in a single parallel round, then
judges it exactly as a climb round. This trades the climb's
minimal calculation count for lower wall-clock latency, the
better bargain when the prediction is trustworthy. The adaptive
serial **climb** is reserved for the cold and moderate-confidence
cases, where the search does not know in advance how far it must
go; it is the bracket-then-refine search of 3.12.3 by default,
with the fine unit-step climb available as an explicit option for
the most conservative reading. All three shapes -- grid,
bracket-refine, unit climb -- share the rung rule (3.12.2), the
two-sided stop test and the metal test (3.12.3), and the harvest
(7.8); they differ
only in how the ladder is sampled: the grid lays a fixed
neighbourhood all at once, the bracket-refine climb strides then
fills the bracket it lands in, and the unit climb walks every
rung.

#### 3.12.6 Relationship to the rest of the chain

This section supersedes the *fixed verify grid* of predict-then-
verify (7.7); the predictor (7.6) and the harvest (7.8) are
unchanged, since they already speak in density and only ever
needed a converged rung to record. The dataspace schema (7.2)
gains one field -- the converged rung's resolved mesh, recorded
in the verification block beside its density (3.12.4). The
duplicate-rung guard (3.11 / 7.8) remains valid as a safety net.
The mesh selection (3.7) and shift selection (3.9) are reused
verbatim as the climb's rung and shift rules.

The design choices are settled above; what remains is numeric
tuning, best fixed by experiment on the seed set. These are
empirically-justified knobs, so by Principle 11 they live in an
auditable config -- the manifest `[harvest]` block (5.7), where
`kpoint_convergence_threshold` (the metric threshold) already
lives, or the site rc for the cost budget -- never as hardcoded
script constants:

- The `predictor_confidence` thresholds that gate one-versus-two
  flat rungs (3.12.3) and the grid-versus-climb dispatch mode
  (3.12.5).
- The bracket phase's geometric stride growth and any cap on the
  largest stride, plus the choice of climb shape -- bracket-
  refine (the default) or the fine unit-step climb (3.12.2 /
  3.12.3 / 3.12.5).
- The multiple by which the bracket phase's flatness threshold is
  looser than the convergence threshold (3.12.3), which sets how
  eagerly a nearly-settled stride is bracketed and so how much of
  the top-end overshoot is shaved.
- The gap threshold below which a rung's computed band gap counts
  as essentially zero, so the climb calls the material a metal and
  settles at once on a rough, floor-level mesh (3.12.3); low enough
  that no real insulator crosses it, high enough to catch a true
  metal's near-zero reading.
- The per-axis cap that sets the crystalline climb's opening floor
  (3.12.4) -- the most points any axis of the coarsest starting
  mesh may carry, `4` by default (so `[4,4,4]` for a cubic cell,
  fewer per axis for an anisotropic one). Keeping the climb above
  this floor is what lets the gap test (3.12.3) read a
  trustworthy gap at the opening rung, with no coarse-mesh guard.
- The value of the fixed per-axis count ceiling, and the cost
  budget once the resource dataspace (8) supplies one (3.12.3).
- The width of the confident-mode fixed mesh grid -- how many
  rungs to each side of the predicted one (3.12.5), and how far
  below the prediction the climb mode begins (3.12.4).

---

## 4. UFF Bond Parameter Database

### 4.1 Motivation

The current `bonds.dat` enumerates specific element pairs
with hand-tuned force constants and rest bond lengths.  It
covers only six elements (H, B, C, N, O, Si) across 14 bond
types.  Every force constant is a uniform 5000.0 kcal/mol/A^2
-- roughly 15 times stiffer than physically realistic values
for typical covalent bonds.  Expanding coverage to 36 or more
elements in the pair-listing format is impractical: 36
elements yield up to 666 unique pairs, and 54 elements
yield 1,485.

A scalable alternative stores per-element parameters and
computes the equilibrium bond length and harmonic force
constant for any pair on the fly, using the Universal Force
Field (UFF) formulas (Rappe et al. 1992).  The new file is
named `bond_parameters.dat` to distinguish it from the legacy
`bonds.dat` format.

### 4.2 UFF Bond Stretching Model

The UFF describes bond stretching with a harmonic potential:

  E = K_ij * (r - r_ij)^2

where K_ij is the force constant (kcal/mol/A^2) and r_ij
is the equilibrium bond length (Angstroms).  The LAMMPS
`bond_style harmonic` uses this same convention: the K
parameter absorbs the factor of 1/2 that appears in the
physics textbook form E = (1/2) k x^2.

For any element pair (i, j), UFF defines:

**Equilibrium bond length:**

  r_ij = r_i + r_j - r_EN

  r_EN = r_i * r_j * (sqrt(chi_i) - sqrt(chi_j))^2
         / (chi_i * r_i + chi_j * r_j)

where r_i and r_j are single-bond covalent radii
(Angstroms) and r_EN is an electronegativity correction
that shortens bonds between elements of unequal
electronegativity.  For homonuclear bonds (same element),
chi_i = chi_j so r_EN = 0 and the equilibrium length
reduces to r_ij = 2 * r_i.

**Force constant (LAMMPS harmonic convention):**

  K_ij = 332.06 * Zstar_i * Zstar_j / r_ij^3

The prefactor 332.06 = 664.12 / 2 absorbs the 1/2 that
converts from the UFF convention E = (1/2) k (r-r0)^2 to
the LAMMPS convention E = K (r-r0)^2.  Zstar_i and Zstar_j
are the UFF effective charges (dimensionless).  The 664.12
constant carries units of kcal*A/mol and encodes the
fundamental relationship between bond stiffness, effective
nuclear charges, and bond length.

Each element requires only three tabulated parameters:

  Parameter  Meaning                       Units
  --------------------------------------------------
  r_i        Single-bond covalent radius   Angstroms
  Zstar_i    Effective charge              (none)
  chi_i      GMP electronegativity         eV

The formula is symmetric in i and j, so the computed values
are independent of element ordering.  However, the calling
code continues to enforce the Z_1 <= Z_2 convention used
throughout the codebase (bond analysis output, tag
construction, pair matching).  The `get_bond_params` method
accepts Z arguments in any order but the callers in
`create_lammps_files` canonicalize to Z_1 <= Z_2 before
constructing bond tags, exactly as the current code does.

**Validation against established force fields.**  UFF-derived
values agree with AMBER within 10-20 %, which is typical
inter-force-field variation for single bonds:

  Bond   UFF K       UFF r0   AMBER K     AMBER r0
         (kcal/      (A)      (kcal/      (A)
         mol/A^2)             mol/A^2)
  ----------------------------------------------------
  C-H    ~331        1.11     ~340        1.09
  C-C    ~350        1.51     ~310        1.53
  C-N    ~360        1.44     ~337        1.47
  O-H    ~540        0.98     ~553        0.96
  Si-O   ~285        1.63     --          1.61

### 4.3 Element Coverage

The UFF provides parameters for every element from Z = 1
through Z = 103.  The initial database covers Z = 1 through
Z = 54 (hydrogen through xenon), spanning:

- All main-group elements through the 5th period
- All 3d transition metals (Sc through Zn)
- All 4d transition metals (Y through Cd)
- Halogens, chalcogens, and pnictogens
- Noble gases (He, Ne, Ar, Kr, Xe)

Noble gases are included for table completeness.  Their very
small effective charges yield negligible force constants, so
they will not produce meaningful bonds in practice.

**Contiguity requirement.**  The table must contain every
element from Z = 1 through `NUM_UFF_ELEMENTS` with no gaps.
The validation check in `get_bond_params` tests
`z > num_uff_elements`, so a gap would leave an uninitialized
slot that could silently produce wrong results.  Extension
beyond Z = 54 requires appending rows for every Z up to the
new maximum -- no code changes are needed.

### 4.4 New bond_parameters.dat File Format

The new file `bond_parameters.dat` replaces `bonds.dat`.  Its
format changes from pair-enumeration to a per-element parameter
table.  Comment lines (beginning with `#`) are permitted and
skipped by the reader.  The tagged-section structure is
preserved for consistency with other Imago data files:

```
# UFF bond stretching parameters.
#
# Source: Rappe, A. K.; Casewit, C. J.; Colwell, K. S.;
#   Goddard, W. A., III; Skiff, W. M.
#   J. Am. Chem. Soc. 1992, 114, 10024-10035.
#   DOI: 10.1021/ja00051a040
#
# For any element pair (i, j):
#   r_ij = r_i + r_j - r_EN
#   r_EN = r_i r_j (sqrt(chi_i) - sqrt(chi_j))^2
#          / (chi_i r_i + chi_j r_j)
#   K_ij = 332.06 Zstar_i Zstar_j / r_ij^3
#          (kcal/mol/A^2, LAMMPS harmonic convention)
NUM_UFF_ELEMENTS
54
UFF_BOND_PARAMS
#  Z    r_i     Zstar_i  chi_i    Element
   1   0.3540  0.7120   4.5280  # H
   2   0.8490  0.0980   9.6600  # He
   3   1.3360  1.0260   3.0060  # Li
  ...
  54   ...                       # Xe
```

Each data line provides: the atomic number Z, the covalent
radius r_i (Angstroms), the effective charge Zstar_i, and the
GMP electronegativity chi_i (eV).  The element comment at the
end of each line is optional but aids readability.

**Reader semantics.**  The Z column on each row is used as
the array index: the reader stores the parameters at
position Z in the arrays, not at the sequential row number.
This makes the file order-independent and robust against
accidental reordering.  If a Z value appears that is outside
the range 1..`NUM_UFF_ELEMENTS`, the reader exits with an
error.

### 4.5 Bond Scale Factor

A new parameter `bond_parameter_scale` provides a global
multiplier for all bond force constants.  The default value
(1.0) is hardcoded in `Condense.__init__` and can be
overridden by a `bond_parameter_scale` keyword in the
condense.in input file (read by `parse_input_file()`):

  bond_parameter_scale 0.9

Force-field parameters are deliberately kept out of
`condenserc.py` and `ScriptSettings` to avoid cluttering
the CLI with rarely-touched calibration knobs.  See 4.8.8
item 5 for the broader rationale and for the parallel
treatment of the angle parameters.

The value is dimensionless and multiplies every computed K_ij
before writing the LAMMPS Bond Coeffs section.  Values below
1.0 loosen all bonds; values above 1.0 stiffen them.

Only K_ij is scaled -- the equilibrium bond length r_ij is
left unchanged.  This lets the user tune overall bond rigidity
while preserving the equilibrium geometry of the system.  The
motivation is that UFF force constants are approximate by
nature (10-20 % inter-force-field variation is typical), so a
global rescaling provides a simple empirical knob for tuning
the dynamic behavior of the condensation simulation without
modifying the underlying database.

### 4.6 Impact on condense.py

The `BondData` class is restructured:

1. **Reading.**  `init_bond_data()` reads the per-element
   parameter table (r_i, Zstar_i, chi_i) and stores it in
   arrays indexed by atomic number Z.  This replaces the
   pair-enumerated `hooke_bond_coeffs` list.

2. **Querying.**  A new method `get_bond_params(z1, z2)`
   computes and returns (K_ij, r_ij) for any element pair
   using the UFF formulas from section 4.2.  Argument order
   does not matter (the formula is symmetric).

3. **Bond lookup in create_lammps_files.**  The existing
   linear scan over `hooke_bond_coeffs` is replaced by a
   single call to `get_bond_params(z1, z2)` -- O(1) per
   lookup instead of O(n).

4. **Bond lookup in normalize_types.**  The same linear
   scan appears a second time in `normalize_types()`,
   where unique bond types are matched to coefficients
   for rewriting LAMMPS files with unified type indices.
   This scan is replaced by `get_bond_params(z1, z2)` in
   exactly the same way as item 3.

5. **Scale factor plumbing.**  The default value of
   `bond_parameter_scale` (1.0) is defined in
   `condenserc.py` and loaded by `assign_rc_defaults()`.
   It can then be overridden by the `bond_parameter_scale`
   keyword in condense.in (read by `parse_input_file()`).
   The multiplier is applied to K_ij in **both** output
   paths: `create_lammps_files` (initial LAMMPS data file)
   and `normalize_types` (rewritten LAMMPS data file with
   unified type indices).  Both paths must produce the
   same scaled force constants.

6. **Error handling.**  If either Z exceeds the range of
   the parameter table, print the element symbol and Z
   number and exit with a clear message directing the user
   to extend the bond_parameters.dat table.

### 4.7 Backward Compatibility

The new `bond_parameters.dat` is not readable by the Perl
`BondData.pm` module, which expects the old `bonds.dat`
pair-listing format.  Since `condense.py` is the active
development path and `BondData.pm` belongs to the deprecated
Perl toolchain, this is accepted.  Users who still need the
Perl `condense` script can retain a local copy of the old
`bonds.dat` file.

**Build system.**  `src/data/CMakeLists.txt` must be updated
to install `bond_parameters.dat` instead of `bonds.dat` in
the DATABASES list.  The old `bonds.dat` is removed from the
install set.

### 4.8 Geometry-Derived Angle Parameters

#### 4.8.1 Motivation

The old `angles.dat` file listed 56 explicit triplet entries
covering only seven elements (H, B, C, N, O, Si).  Every
entry used a uniform spring constant k = 500.0 kcal/mol/rad^2
regardless of the element triplet.  Adding a new element
required manually enumerating every triplet and rest angle
it participates in -- an unsustainable maintenance burden and
a frequent source of "Cannot find angle in the database"
failures when the system contains elements outside the seven.

Section 4.8 of the previous design (preserved below in
section 4.8.2) analyzed why the per-element UFF strategy
that worked for bonds does not transfer directly to angles:
the same element triplet can have multiple physically
distinct rest angles (e.g., C-C-C at 60, 109.5, 120, and
180 degrees), and the UFF angle potential is a cosine
Fourier series rather than a simple harmonic.

The key insight is that the Imago bond analysis already
computes the actual bond angles for every atom in every
molecule.  These observed angles encode the real electronic
structure -- hybridization, strain, ring membership, and
neighbor effects -- for the specific system at hand.  They
are more accurate than any generic lookup table, and they
are already available at runtime.  The design below uses
these angles directly as equilibrium values, eliminating
the need for an external angle database entirely.

#### 4.8.2 Prior Analysis (retained for reference)

**Why per-element UFF does not transfer directly to angles.**

1. **Multiple rest angles per triplet.**  The same element
   triplet (e.g., C-C-C) appears in angles.dat with several
   distinct equilibrium angles (60 deg for cyclopropane,
   108 deg for cyclopentane, 180 deg for linear chains).
   A per-element UFF lookup gives one natural angle per
   vertex atom type (e.g., C_3 = 109.47 deg), which cannot
   distinguish these chemical environments.

2. **UFF angle potential form mismatch.**  The UFF angle
   bending potential is a cosine Fourier series (Rappe
   eq. 8), not a simple harmonic.  LAMMPS `angle_style
   harmonic` uses E = K (theta - theta_0)^2.  Adopting
   UFF angle parameters would require either approximating
   the cosine form as harmonic near theta_0 or switching
   to a different LAMMPS angle style -- both are significant
   changes beyond the data file.

3. **Complex force constant formula.**  The UFF angle force
   constant K_IJK (Rappe eq. 13) depends on the bond lengths
   of both arms (r_IJ and r_JK), all three effective charges,
   and the vertex atom's natural angle.  This is considerably
   more involved than the clean two-element bond formula
   K_ij = 332.06 * Zstar_i * Zstar_j / r_ij^3.

The geometry-derived approach (sections 4.8.3-4.8.9) resolves
all three issues: it uses observed angles (solving 1), feeds
them into a LAMMPS harmonic potential (avoiding 2), and uses
a simplified force constant formula based on bond stiffnesses
already computed by `get_bond_params()` (simplifying 3).

#### 4.8.3 Approach: Cluster Observed Angles by Triplet

For each molecule in the system, the Imago bond analysis
(via `bond_analysis.py`) already computes every bond angle.
The `create_lammps_files` method already iterates over these
angles and constructs triplet tags of the form
(Z_end1, Z_vertex, Z_end2).  Currently it searches
`angles.dat` for a matching entry.  The new approach replaces
that database lookup with the following procedure:

1. **Collect.**  For each angle instance, extract the
   full triplet (Z1, Zv, Z2) with Z1 <= Z2, and the
   observed angle theta_obs.

2. **Cluster by triplet.**  Group all angle instances
   that share the same (Z1, Zv, Z2) triplet.  Within
   each triplet group, sort the observed angles and
   greedy-merge values into a cluster while two
   conditions hold: (a) the candidate is within +/-
   `angle_cluster_tolerance` degrees of the running
   mean, and (b) the resulting cluster span (max - min)
   remains within `2 * angle_cluster_tolerance`.  The
   spread cap prevents a long chain of closely-spaced
   observations from silently sweeping values from
   opposite ends of a wide distribution into a single
   cluster.  When either condition fails, the current
   cluster is finalized and a new one begins at the
   candidate.  The cluster's rest angle theta_0 is the
   mean of its members.  The same spread cap is applied
   when `normalize_types()` re-clusters across sources
   (see 4.8.8 item 4a), so local and cross-source
   clustering use consistent semantics.

3. **Assign types.**  Each cluster becomes one LAMMPS
   angle type.  Every angle instance is assigned to the
   cluster whose mean it contributed to.

**Example.**  Suppose carbon vertex atoms yield observed
angles of 108.3, 109.1, 109.8, 120.2, 119.7, and 60.1
degrees, all for the C-C-C triplet, with
`angle_cluster_tolerance = 5.0`:

- Cluster 1: {60.1} -> theta_0 = 60.1 (ring)
- Cluster 2: {108.3, 109.1, 109.8} -> theta_0 = 109.1 (sp3)
- Cluster 3: {119.7, 120.2} -> theta_0 = 120.0 (sp2)

This produces three angle types instead of six individual
entries.

#### 4.8.4 Force Constant Formula

The angular spring constant K is computed from the UFF
per-element parameters already stored in
`bond_parameters.dat` (section 4.4).  The formula uses
the bond stiffnesses of the two arms:

  K_angle = C_angle * sqrt(K_bond_IJ * K_bond_JK)

where K_bond_IJ is the UFF harmonic bond force constant
for the (Z1, Zv) pair and K_bond_JK is for the (Zv, Z2)
pair, both obtained from `get_bond_params()`.  The
geometric mean captures the essential physics: stiffer
bonds produce stiffer angles.  The calibration constant
C_angle is dimensionless and converts bond stiffness
(kcal/mol/A^2) into an angular stiffness scale
(kcal/mol/rad^2).

Unlike the UFF bond constant (332.06, well-established
from the Rappe paper), C_angle is a project-specific
heuristic with no published source (see Provenance below).
It is therefore exposed as a user-tunable keyword in
condense.in:

  angle_stiffness_coeff 0.15

The default value (0.15) is defined in `condenserc.py`.
Together with `angle_parameter_scale` (section 4.8.5),
the user has two complementary controls:
`angle_stiffness_coeff` sets the base conversion from
bond stiffness to angle stiffness, while
`angle_parameter_scale` applies a uniform global
multiplier on top.  The final force constant written to
LAMMPS is:

  K_final = angle_stiffness_coeff
            * sqrt(K_bond_IJ * K_bond_JK)
            * angle_parameter_scale

**Provenance.**  This formula is a project-specific
heuristic, not drawn from a published force field.  The
full UFF angle bending force constant K_IJK (Rappe et al.
eq. 13) depends on the bond lengths of both arms, all
three effective charges, the equilibrium angle, and uses
a cosine Fourier expansion rather than a harmonic
potential.  Adopting the full UFF angle treatment would
require either switching LAMMPS to a cosine angle style
or performing a non-trivial harmonic approximation of the
Fourier series near each equilibrium angle.  The geometric
mean heuristic sidesteps both issues by staying within
the LAMMPS `angle_style harmonic` framework (Thompson et
al. 2022; E = K (theta - theta_0)^2) while still
producing element-dependent K values that track the
underlying bond stiffnesses.

**Calibration.**  Published harmonic angle force constants
for small organic molecules typically fall in the range
30-100 kcal/mol/rad^2.  For reference, the AMBER ff94
force field (Cornell et al. 1995) assigns C-C-C angles
K ~ 40 kcal/mol/rad^2 and H-C-H angles K ~ 35
kcal/mol/rad^2; the OPLS-AA force field (Jorgensen et al.
1996) gives similar values.  These are considerably softer
than the uniform k = 500 used in the old `angles.dat`.

Typical UFF bond force constants from `get_bond_params()`
are 200-700 kcal/mol/A^2.  For a C-C-C angle, both arms
give K_bond ~ 470 kcal/mol/A^2, so sqrt(470 * 470) = 470.
The default `angle_stiffness_coeff` of 0.15 yields
K_angle ~ 70 kcal/mol/rad^2, which is within the range
of published values.  Users should calibrate this value
against a known system (e.g., a small organic molecule
with published force field parameters).

**Note on the uniform k = 500 in the old database.**  The
old `angles.dat` used k = 500 for every entry.  This is
extremely stiff -- roughly 5-10x typical literature values.
It is not physically motivated; it appears to have been
chosen as a "rigid enough" default.  The computed K values
from the formula above will be significantly softer and
more physically realistic.  If the user needs the old
stiff behavior, `angle_parameter_scale` can be set to a
large value (e.g., 5.0-7.0).

#### 4.8.5 Angle Scale Factor

A new parameter `angle_parameter_scale` provides a global
multiplier for all computed angle force constants, following
the same pattern as `bond_parameter_scale` (section 4.5).
The default value (1.0) is hardcoded in `Condense.__init__`
and can be overridden by a keyword in condense.in (read by
`parse_input_file()`):

  angle_parameter_scale 0.8

The value is dimensionless and multiplies every computed
K_angle before writing the LAMMPS Angle Coeffs section.
Values below 1.0 loosen all angular springs; values above
1.0 stiffen them.  Only K_angle is scaled -- the rest angle
theta_0 is left unchanged.

#### 4.8.6 Angle Cluster Tolerance

A new parameter `angle_cluster_tolerance` controls how
aggressively observed angles are merged into shared types.
The default value (5.0 degrees) is hardcoded in
`Condense.__init__` and can be overridden in condense.in
(read by `parse_input_file()`):

  angle_cluster_tolerance 3.0

**Scope.**  This parameter and the clustering algorithm
described below apply only to condense.py's two call
sites -- local clustering in `create_lammps_files()` and
cross-source re-clustering in `normalize_types()`.  The
template producer `make_reactions.py` uses a fixed
tolerance of 0.0 (identity-only merge) regardless of
this parameter's value, which is what makes reaction
templates reusable across any condense.py simulation;
see section 4.8.10 for the full rationale.

**Clustering algorithm.**  Within each (Z1, Zv, Z2)
triplet group, the observed angles are sorted in
ascending order and then merged greedily: the first angle
starts a new cluster; each subsequent angle is added to
the current cluster if it falls within
`angle_cluster_tolerance` of the cluster's running mean,
otherwise it starts a new cluster.  This greedy approach
is simple, deterministic, and keeps the type count low.

Note that a chain of angles spaced just under the
tolerance apart (e.g., 105, 108, 111, 114 with a 5-degree
tolerance) will merge into one cluster because each new
member is compared to the running mean, not to the first
member.  This is the intended behavior: it favors fewer,
broader clusters, which reduces the angle type count and
lowers the risk of bond/react type-mismatch failures.

**Interaction with bond/react type count.**  LAMMPS
bond/react requires that the atom, bond, and angle types
in the pre- and post-reaction templates match the types
in the main data file.  Every distinct angle type in the
system increases the combinatorial space that must be
consistent across all files.  A larger tolerance produces
fewer, coarser angle types, which reduces the risk of
type-mismatch failures in bond/react.  A smaller tolerance
preserves finer geometric detail but creates more types.

The default of 5.0 degrees is a practical compromise.
For systems with many distinct molecular species or
complex reaction networks, increasing the tolerance to
8-10 degrees may be necessary to keep the type count
manageable.

#### 4.8.7 Look-Ahead Angles for Bond/React Products (deferred)

**The problem.**  The clustering procedure in section 4.8.3
discovers angle types from the initial molecular geometries.
But LAMMPS bond/react creates new bonds between molecules,
and those new bonds produce new angles that did not exist
in any isolated molecule.

Consider B12H12 and CH4.  In the isolated molecules, no
C-B bond exists, so no C-B-H or C-B-B angle is ever
observed.  After bond/react fires and creates a C-B bond,
the post-reaction template would need angle types for
every triplet that includes the new bond.  If those types
are not present in the LAMMPS data file's Angle Coeffs
section, bond/react would fail.

**Current state of the code.**  The Perl `makeReactions`
script (and its Python port `make_reactions.py`) adds one
new *bond* between the trigger atoms in the post-reaction
template but does **not** add any new *angles*.  The Perl
source (lines 2508-2513) states this explicitly:

  "In the future it might be necessary to *add* bond
   angles through the bonding atoms (after the S are
   removed), but presently we do not do that.  (Hence
   the bonded molecules may be too floppy.)"

A commented-out `addBondAngle` subroutine (Perl lines
2990-3141) shows that an attempt was started: it computes
angles via the law of cosines from post-reaction
coordinates, builds angle tags, and registers new angle
types.  But the subroutine was never activated.

Because the post-reaction templates do not currently
contain any new angles, there are no novel angle types
for condense.py to "look ahead" to.  The look-ahead
mechanism in condense.py and the angle-creation logic
in makeReactions are two halves of the same problem.

**Empirical confirmation (2026-04-25, 60-mol B12H12 run).**
Inspection of a representative `make_reactions.py`-generated
postRxn template (`postRxn.b12h12_1_b-1_b12h12_1_b-1.data`)
makes the floppy-joint behavior precise.  preRxn carries 240
angles, postRxn carries 230 -- a net deficit of exactly ten.
The bond count changes correctly (62 -> 61: minus two B-H
bonds for the deleted hydrogens, plus one new B-B bond), but
the angle count is wrong by ten.

The ten missing angles decompose symmetrically.  Around the
side-1 initiator (template atom 1), preRxn has 15 angles
centred on it: ten B-B-B angles among the five B neighbours
{2,3,4,5,6} and five H-B-B angles using the H atom (atom 12)
to be deleted.  postRxn keeps the ten B-B-B angles and drops
the five H-B-B angles -- correct as far as it goes, since
atom 12 no longer exists.  But the five new B-B-B angles
that the new bond 1-18 should produce (`18-1-2`, `18-1-3`,
`18-1-4`, `18-1-5`, `18-1-6`) are simply absent.  The
identical pattern holds on side 2: five missing angles
`1-18-19` through `1-18-23` involving the new bond on the
atom-18 end.  Five plus five equals the observed deficit.

**Physical consequence.**  The new B-B bond carries only a
bond-stretch potential; the angle term contributes nothing
along the new bond axis.  As a result the two cages can
rotate freely about the new B-B bond axis with no restoring
torque, and the bond axis itself can swing relative to the
local cage symmetry without any angular penalty.  In a single
isolated dimer this is already a structural inaccuracy, but
in a chain of N inter-cage reactions the deficit grows
linearly: each joint contributes ten missing angle terms, so
an N-mer chain has 10N unconstrained angular degrees of
freedom at its joints.  This becomes more serious than the
isolated-dimer case once chains are long enough to coil back
on themselves, because unrestricted joint rotation lets
chains reach geometries that the bond/angle topology would
otherwise forbid -- and once the bond/react fix accepts a
reaction, the resulting topology persists for the rest of
the run.

**Phased approach.**  This work is split into two steps
to keep each change testable:

1. **Step 1 (this design, sections 4.8.3-4.8.6 and
   4.8.8):**  Replace the angles.dat database with
   geometry-derived clustering and computed force
   constants.  The system works for all intra-molecular
   angles that already exist in the pre- and post-
   reaction templates.  Post-reaction bonding sites
   remain "floppy" (same as today) because no new
   angles are created by makeReactions.

2. **Step 2 (future work):**  Activate angle creation in
   `make_reactions.py` (porting and completing the
   commented-out Perl subroutine).  Once post-reaction
   templates carry the new angles, condense.py's
   `normalize_types()` will automatically pick them up
   during its template-scanning pass.  The rest angle
   theta_0 can be computed from the post-reaction atom
   coordinates (law of cosines, as the Perl prototype
   does), and K_angle can be computed from the same
   formula (section 4.8.4).  The clustering tolerance
   should be applied when deciding whether a novel
   post-reaction angle merges with an existing type or
   creates a new one.

Step 2 will also need to address whether new *bond* types
(not just angles) can appear in post-reaction templates.
The bond case is simpler because `get_bond_params()` can
always compute K and r0 for any element pair, but the
type must still be registered in the unified type list.

  >> DESIGN QUESTION D5a (deferred to step 2): When
  >> post-reaction angles are added to the templates,
  >> should novel angles always be added as distinct
  >> types to ensure the post-reaction geometry is
  >> exactly preserved, or should they be merged into
  >> the existing cluster list if they fall within
  >> `angle_cluster_tolerance` of an existing type?
  >> Merging keeps the type count down (reducing
  >> bond/react fragility), but exact preservation may
  >> matter for the product geometry.

#### 4.8.8 Impact on condense.py (step 1)

The `AngleData` class is eliminated.  No external data
file is read for angles.  The changes are:

1. **AngleData class: remove.**  The class and its
   `init_angle_data()` method are deleted from every
   source file that currently carries them -- both the
   copy in `condense.py` and the duplicated copy in
   `make_reactions.py` (lines ~113-193) along with its
   `self.angle_data` instantiation (line ~568).  All
   references to `self.angle_data` and
   `ad.hooke_angle_coeffs` are removed from
   `create_lammps_files()`, `normalize_types()`, and
   `make_reactions.py`'s template-emission path (the
   hooke_angle_coeffs scan near line 2518).

2. **Clustering step: add (shared helper).**  A new
   helper routine implements the cluster-by-triplet
   algorithm from section 4.8.3.  Input: the list of
   all angle instances with their (Z1, Zv, Z2,
   theta_obs) tuples from one producer's scope.
   Output: a list of local angle types, each with
   (Z1, Zv, Z2, theta_0) and an observation count,
   plus a mapping from each instance to its local
   type index.  This helper is shared by both
   producers -- `create_lammps_files()` in
   condense.py and the template-emission path in
   `make_reactions.py` -- so that local clustering
   semantics are byte-identical across sources.  The
   cross-source unification described in item 4 then
   operates on the per-source outputs.

3. **Force constant computation: add (condense.py
   only).**  For each angle type in the lammps.dat
   produced by `create_lammps_files()`, compute K_angle
   using the formula from section 4.8.4, calling
   `get_bond_params()` for the two arm bond stiffnesses.
   Apply `angle_stiffness_coeff` and
   `angle_parameter_scale` (both factor into the
   K_angle formula).  `normalize_types()` recomputes the
   same K values in item 4b, which is safe because
   K_angle depends only on the triplet and not on
   theta_0.

   The template producer `make_reactions.py` does **not**
   compute K_angle locally.  Reaction template files
   (pre-, post-, and map-) carry only connectivity,
   per-atom angle entries, and the tag tail
   "{theta_0_local} {t}" -- no K value is ever written
   into a template.  Since `normalize_types()` in item
   4b recomputes K authoritatively from the triplet
   for every final cluster, a redundant producer-side K
   in `make_reactions.py` would be neither written nor
   consumed.  Skipping it keeps `make_reactions.py`
   independent of `BondData` and avoids plumbing
   `angle_stiffness_coeff` and `angle_parameter_scale`
   into a script that never writes their effect to disk.
   See section 4.8.10 for the full rationale -- the
   same decision is what makes reaction templates
   reusable across any condense.py simulation.

4. **normalize_types(): cross-source unification.**
   `normalize_types()` is the authoritative
   cross-source clusterer for angle types, not a
   passive consumer of producer-emitted theta_0
   values.  Each producer (`create_lammps_files()`
   and `make_reactions.py`) locally clusters its
   own observations and emits one type per local
   cluster, with the cluster-mean theta_0 carried
   in the tag tail "{theta_0_local} {t}".  Because
   the two producers see different observation
   populations of the same physical angle, their
   local theta_0 for a chemically identical triplet
   can differ by a few tenths of a degree.  A plain
   string comparison on the tag tail would split
   such cases into distinct types and break
   bond/react type-ID matching across lammps.dat
   and the reaction templates.  `normalize_types()`
   absorbs this drift in four steps:

   a. **Cross-source clustering.**  Collect every
      angle type emitted by every source -- the
      lammps.dat produced by `create_lammps_files()`
      and every reaction template produced by
      `make_reactions.py` -- each carrying
      (z1, zv, z2, theta_0_local, obs_count,
      source, local_type_id).  Group by canonical
      triplet (z1 <= z2).  Within each group, sort
      by theta_0_local and greedy-merge while the
      candidate is within `angle_cluster_tolerance`
      of the running cluster mean.  Apply a spread
      cap (max cluster span <= 2 *
      `angle_cluster_tolerance`) to prevent greedy
      chaining across a wide distribution.  The
      running mean, weighted by obs_count, is the
      final canonical theta_0 for the merged
      cluster.

      **Associativity of obs_count weighting.**  The
      weighted-mean formula used here is associative
      under pre-merging: clustering identical
      observations into one local record with
      obs_count > 1 before cross-source merging
      produces the same final theta_0 as passing
      each observation through as a separate
      weight-1 record.  This property is what makes
      `make_reactions.py`'s local clustering step
      safe at any tolerance that only merges
      bit-identical theta values (see section 4.8.10
      and PSEUDOCODE 10d): the template producer can
      collapse duplicate observations to keep the
      template file compact without perturbing what
      `normalize_types()` computes.  Merging
      non-identical observations at the producer
      would not preserve associativity and would
      create the T_m > T_c hazard described in
      section 4.8.10.

   b. **Force constant computation.**  For each
      final cluster, compute K_angle via
      `get_bond_params()` using the formula from
      section 4.8.4, with `angle_stiffness_coeff`
      and `angle_parameter_scale` applied (same
      formula as item 3).  Because K_angle depends
      only on the arm bond stiffnesses, which are a
      function of the triplet (z1, zv, z2) alone
      and not of theta_0, the K_angle computed here
      matches whatever any producer locally
      computed for a member cluster -- cross-source
      merging does not alter K_angle, only theta_0.

   c. **Tag rewrite and type-ID remap.**  Each
      (source, local_type_id) pair maps to exactly
      one final cluster id.  Walk the Angles
      section of lammps.dat and every angle
      reference in every reaction template, and
      rewrite the per-angle type id to the global
      cluster id.  Rewrite the tag tail
      "{theta_0_local} {t}" to
      "{theta_0_final} {global_t}" so any tool that
      later inspects the tag sees a consistent
      value.  The rewrite is deterministic given
      the cluster map, so repeated runs on
      identical inputs produce byte-identical
      output.

   d. **Cluster-map diagnostic.**  Emit a log file
      or log section listing, for every final
      cluster: global id, canonical theta_0,
      (z1, zv, z2), and every contributing
      (source, local_theta_0, obs_count) tuple.
      Students debugging a bond/react type mismatch
      should be able to open this file and see at a
      glance which observations were merged, which
      were split, and why.  This diagnostic is the
      main debuggability payback for routing all
      clustering through `normalize_types()` rather
      than accepting per-producer string tags.

5. **Parameter plumbing: add.**  Introduce three new
   force-field parameters on the `Condense` class with
   hardcoded defaults in `Condense.__init__`:
   `angle_stiffness_coeff` (default 0.15),
   `angle_parameter_scale` (default 1.0), and
   `angle_cluster_tolerance` (default 5.0).  These
   follow the `bond_parameter_scale` precedent:
   force-field parameters are deliberately kept out of
   `condenserc.py` and `ScriptSettings` to avoid
   cluttering the CLI with rarely-touched knobs.  User
   overrides are accepted through matching keywords in
   `condense.in`, parsed by `parse_input_file()`.

   All three parameters are **condense.py-scoped**:
   they govern the clustering and K computation done
   inside `create_lammps_files()` and
   `normalize_types()`, not anywhere in
   `make_reactions.py`.  `angle_stiffness_coeff` and
   `angle_parameter_scale` appear only where K_angle is
   computed, which per item 3 is condense.py only.
   `angle_cluster_tolerance` governs the two clustering
   sites that live inside condense.py -- local
   clustering in `create_lammps_files()` and the
   cross-source re-clustering in `normalize_types()` --
   along with the matching spread-cap policy
   `2 * tolerance` at each.  The template producer
   `make_reactions.py` uses a fixed tolerance of 0.0 at
   its local clustering site (identity-only merge over
   the 0.5-degree-quantized observed angles from
   `bondAnalysis.ba`), so no user-tunable tolerance is
   exposed there and no parameter coordination is
   required between the two scripts.  Section 4.8.10
   explains why this asymmetry is what makes reaction
   templates reusable across any condense.py simulation
   regardless of its `angle_cluster_tolerance` value.

6. **angles.dat: retire.**  Remove from
   `src/data/CMakeLists.txt` DATABASES list.  Remove
   from `share/` install target.  The file may be kept
   in the repository for historical reference but is no
   longer read by any code path.

Note: the look-ahead pass described in section 4.8.7 is
deferred to step 2.  Step 1 handles all angles that
already exist in the pre- and post-reaction templates
(which is the same set that the old angles.dat handled).
The bonding-site "floppiness" is unchanged from the
current behavior.

#### 4.8.9 Backward Compatibility

The new approach is not backward-compatible with the old
`angles.dat` format.  Since `condense.py` is the active
development path and the Perl toolchain is deprecated,
this is accepted (same reasoning as section 4.7 for bonds).

The behavioral difference is that angle rest values will
now come from the system's own geometry rather than a
curated database.  For well-prepared input structures
(which is the expected use case), this produces identical
or better rest angles.  For poorly prepared structures,
the rest angles will reflect the input geometry -- which
is arguably more honest than imposing idealized angles
that the structure does not actually have.

Users who relied on the old uniform k = 500 behavior can
approximate it by setting `angle_parameter_scale` to a
large value, though the per-triplet variation in K will
still be present.

#### 4.8.10 Template Reusability (make_reactions.py tolerance = 0)

The reaction template files produced by
`make_reactions.py` (pre-, post-, and map-files) are
intended to be reusable across any `condense.py`
simulation, regardless of the `angle_cluster_tolerance`
value that the downstream simulation happens to use.
This reusability is not automatic -- it depends on the
template producer being careful about what it does with
its own local clustering step.  This section explains
the constraint, the failure mode if it is violated, and
why the chosen design (tolerance = 0 at
`make_reactions.py`) makes the guarantee hold.

**Producer and consumer tolerances.**  Let T_m be the
tolerance `make_reactions.py` uses for its local
clustering and T_c be the tolerance `condense.py` uses
inside `normalize_types()` for cross-source clustering.
The two scripts are independent: a given template is
generated once and may later be fed to many different
`condense.py` runs that each set T_c however the user
pleases.

**The asymmetric hazard: T_m > T_c is dangerous, T_m <=
T_c is safe.**  Local clustering in either producer is
in principle an optimization: it compresses the record
stream that `normalize_types()` consumes.  But that
optimization is only semantically neutral while it
merges bit-identical observations.  Once T_m is large
enough to merge non-identical observations, the
producer has collapsed information that the consumer
cannot recover.

- **T_m < T_c (stricter local, looser global).**  The
  producer emits finer local types.  The consumer's
  looser cross-source merge can still fold those finer
  splits together correctly, because the weighted-mean
  formula in item 4a is associative (see the
  associativity note there).  No physics error; at
  worst, slightly more records pass through the cross-
  source step.
- **T_m > T_c (looser local, stricter global).**  The
  producer has already fused physically distinct angles
  into a single local type with a single theta_0 and a
  single local_type_id.  The consumer sees one record
  for what the user wanted to treat as two types.  All
  occurrences in the template share one global type
  after the cross-source step, so LAMMPS applies the
  wrong equilibrium angle and wrong force constant to
  some fraction of the angles, off by up to 2 * T_m
  degrees.  This is not a crash -- the simulation runs
  -- but the results are wrong, and debugging the
  discrepancy requires tracing local_type_id lineage
  through both producers.

**The chosen design: T_m = 0 at `make_reactions.py`.**
Fixing T_m to 0 makes the guarantee symmetric and
absolute: any T_c >= 0 is safe, and T_c >= 0 covers
every valid configuration of `condense.py`.  At
T_m = 0, `cluster_angles()` (PSEUDOCODE 10a) merges
only bit-identical theta_obs values.  Observed angles
in `bondAnalysis.ba` are already quantized to 0.5
degrees by the reader (see `make_reactions.py
_read_angle_data`), so bit-identity is a well-defined
operation and collapses duplicates cleanly -- e.g. a
benzene ring's six geometrically identical C-C-C angles
at 120.0 degrees become a single local record with
`obs_count = 6` rather than six separate records.

**What does not happen under T_m = 0.**  Because no
non-identical angles are ever fused at the producer,
no interpretive merging is performed before the
consumer sees the data.  `normalize_types()` is
therefore free to apply any T_c it likes -- the raw
0.5-degree-quantized local_theta_0 values are still
present in every record, and the weighted-mean cross-
source cluster in item 4a produces the same final
theta_0 for any T_c that it would have produced if the
template had listed each observation separately.  The
associativity property in item 4a is what formally
guarantees this.

**Why not enforce T_m via a manifest instead.**  An
earlier design sketch considered having
`make_reactions.py` write a manifest recording the T_m
it used, and having `condense.py` verify on read.  Two
reasons that approach was rejected:

1. **Needless coordination.**  The manifest introduces
   a hand-off protocol between two independent scripts
   that otherwise share nothing but the template file
   format.  Fixing T_m = 0 removes the coordination
   entirely -- there is no parameter to record, check,
   or propagate.

2. **User experience under mismatch.**  A manifest
   check would force users to regenerate templates any
   time they tuned `angle_cluster_tolerance` in
   `condense.in`.  Since templates are often produced
   by one researcher and consumed by many, regenerating
   them is not cheap.  Fixed T_m = 0 means the
   condense.in tolerance can be re-tuned freely
   without touching the template files.

The cost of T_m = 0 is modestly larger template files
(one local record per distinct 0.5-degree-quantized
observed angle rather than one per coarse cluster),
but the observation data per template is small in
absolute terms and the reusability guarantee is worth
the trade.

---

## 5. Initial SCF Potential Database

### 5.1 Overview

This section pins down the algorithms, data structures,
and exact TOML schema for the augmented initial-potential
database introduced in VISION Goal 3 and architected in
ARCHITECTURE Section 8.  The augmented per-element file
`share/atomicPDB/<elem>/s_gaussian_pot.toml` carries
multiple labeled potential entries; `makeinput.py` reads
the file at input-generation time and emits the chosen
entry's numerical content into the Imago input file in
today's on-the-wire format.  The Fortran side does not
change.

### 5.2 TOML Schema (version 2)

Schema v2 extends the Phase-1 v1 schema with two additions
needed to support fingerprint-driven manifest selection
(5.6) without churning the on-the-wire Imago input format:

- A required per-entry `default` boolean.  Exactly one
  `[[potential]]` per file carries `default = true`; this
  is the entry picked when no scheme matches.
- An optional `[[potential.fingerprint]]` inner array on
  every entry.  Each fingerprint record carries a
  `method` name, a method-specific `sub_spec` inline
  table, and a method-specific payload.  The schema
  validates structural presence and `(method, sub_spec)`
  uniqueness; payload shape is the responsibility of the
  matcher that consumes the record (ARCHITECTURE 8.9).

There is no on-disk compatibility with v1 files: the
reader rejects any `schema_version != 2`.  Existing v1
files (none in production at the time of the bump) are
regenerated by the producer (5.7), which adds the
`default` tag and the fingerprint records the curator
declares.  That reject-and-regenerate stance was free only
because nothing was in the field to lose; 5.2.5 states the
versioning rule and the reader's policy for the bumps that
follow, once a curated database is worth carrying forward.

**Top-level keys (required):**

  Field            Type    Description
  --------------------------------------------------------
  schema_version   int     Currently 2.  A file at any
                           other version is migrated or
                           refused by the version gate
                           (5.2.5); the reader never
                           reads one as-is.
  element_symbol   string  E.g., "Au".  Must match the
                           parent directory name.
  nuclear_z        real    Atomic number Z.  Nominally an
                           integer, but stored and emitted as a
                           real: Imago consumes Z as a real
                           number, and the legacy pot1 file
                           already records it as one.
  nuclear_alpha    float   Alpha in the nuclear potential
                           form Z * exp(-alpha * r^2).
                           Per-element constant; the same
                           value applies to every entry
                           in the file.
  covalent_radius  float   Reserved for future use.
                           atomSCF currently writes 1.0.

**Per-entry keys, under each `[[potential]]` (required):**

  Field            Type    Description
  --------------------------------------------------------
  label            string  Unique within file.  The label
                           "isolated" is reserved for the
                           atomSCF-derived entry; the
                           label "default_solid" is the
                           single-bulk improved entry
                           (Phase 1 deliverable).
  default          bool    Selection hint.  Exactly one
                           entry per file carries true
                           (rule 7).  The entry used when
                           no scheme matches and no
                           `-pot LABEL` override is given.
                           Independent of `"isolated"`:
                           the curator may mark either
                           the isolated baseline or a
                           curated improved entry as the
                           default.
  description      string  Human-readable note.
  num_gaussians    int     Count of Gaussian basis
                           functions in the electronic
                           potential expansion.
  alpha_min        float   Smallest alpha in the original
                           geometric series from which
                           alphas[] was generated.
                           Informational.
  alpha_max        float   Largest alpha in the original
                           geometric series.
                           Informational.
  coefficients     array   Length = num_gaussians.  The
                           entry's potential: the
                           harvested potential of the
                           first representative atom that
                           mapped to this environment,
                           stored verbatim (5.2.3).  A
                           statistical mean over every
                           mapped atom -- with its
                           per-coefficient spread, an atom
                           multiplicity, and a model count
                           -- is a deferred refinement
                           (5.2.3, "Deferred").
  alphas           array   Length = num_gaussians.
                           Authoritative explicit alpha
                           values per basis function
                           (col 2 of coeff1).  May
                           diverge from the geometric
                           series implied by alpha_min
                           and alpha_max if a future
                           entry uses a non-geometric
                           layout.  In the present regime
                           every atom of an element shares
                           one set of alphas, so all
                           entries in a file carry
                           identical alphas.  (The deferred
                           statistical merge of 5.2.3 would
                           assert this equality before
                           averaging coefficients.)

**Per-entry fingerprint records, under
`[[potential.fingerprint]]` (optional, may repeat):**

Each fingerprint record describes the local environment
of the reference atom site at the moment that entry's
numerical potential was harvested.  The matcher (8.9)
named by `method` interprets the `sub_spec` and the
record's remaining fields.

  Field      Type          Description
  ------------------------------------------------------
  method     string        Matcher name.  Currently one
                           of `"reduce"` or
                           `"bispectrum"`.  The matcher
                           defines the meaning of
                           `sub_spec` and the names and
                           shapes of additional fields
                           on the record.
  sub_spec   inline table  Method-specific parameters
                           that fully qualify the
                           fingerprint.  For
                           `"bispectrum"`, e.g.,
                           `{ twoj1 = 8, twoj2 = 8 }`.
                           For `"reduce"`, e.g.,
                           `{ level = 2, thick = 0.5,
                              cutoff = 5.0,
                              tolerance = 0.05 }`.
                           Two fingerprint records with
                           the same `method` but
                           different `sub_spec` keys or
                           values are non-comparable and
                           coexist on the same entry
                           (rule 8).
  preferred  bool          Optional, default false.  Marks
                           this entry's canonical record for
                           a matcher family: the one whose
                           `sub_spec` names the settings the
                           consumer's file-dictated
                           (crystalline) match computes its
                           query with (5.6.5 step 2), and the
                           one the dedup keys on (5.2.3).
                           Scoped to the ENTRY: exactly one
                           record per family present on that
                           entry carries `preferred = true`
                           (rule 10), so every harvested entry
                           flags its own.  The preferred
                           `sub_spec`
                           for a family is uniform across the
                           whole database (set once in the
                           curation manifest, 5.7).  Storing
                           extra non-preferred sub_specs of
                           the same family is always allowed;
                           only the divergent *preferred
                           flag* is forbidden.

Additional fields on the record are matcher-specific and
not validated by the schema.  As examples for the two
matchers Phase 2 ships with:

- `"bispectrum"` records carry `values` (array of reals,
  length `twoj2 + 1` -- the count of coupling channels `j`
  in `|j1 - j2| <= j <= j1 + j2`, with `twoj1 >= twoj2`).
- `"reduce"` records carry a `shell_code` inline table:
  the central atom's `element` symbol plus a `levels`
  array, one entry per reduction level holding that
  shell's `distance` and a `neighbors` list of neighbor
  element symbols (the neighbor count is implicit in the
  list length).  The shells are built over a periodic
  neighbour list, so a neighbour is counted once per
  periodic image and the counts are those of the physical
  environment rather than of the chosen cell -- which is
  what lets one structure's record be matched against
  another's (5.11).  The neighbor multiset is element-only --
  *not* element/species.  Species numbering is local to a
  single structure (one structure's "species 2" has no
  relation to another's), so it would not transfer to the
  query structures this stored fingerprint is later
  matched against (5.6.5); element symbols are global and
  transferable, so the cross-structure descriptor keeps
  only them.  Within `group_reduce` the species component
  still distinguishes atoms, but that comparison never
  leaves one structure.  All symbols are lowercased to
  match the CLI element/species token convention.

Fingerprint records inherit provenance from their parent
`[potential.provenance]` block: the same reference run
that produced the numerical potential also produced the
fingerprint.

**Per-entry provenance, under `[potential.provenance]`
(required):**

  Field           Type    Description
  --------------------------------------------------------
  source          string  "atomSCF" or "Imago".
  commit          string  Git SHA of the generating tool
                          at write time.
  generated_at    string  ISO-8601 UTC timestamp.

  Additional fields, required when source == "Imago":

  Field                  Type    Description
  --------------------------------------------------------
  reference_id           string  E.g., "COD-XXXXXXX",
                                 "MP-XXXXXX", or a local
                                 tag for in-house
                                 reference structures.
  atom_site              int     Site index in the
                                 reference structure
                                 whose converged
                                 potential is captured.
  kpoint_spec            string  Free-form record of the
                                 k-point specification
                                 used in the reference
                                 run (axial counts +
                                 shifts, density value,
                                 etc.).  Recorded for
                                 provenance only.
  type_assignment        string  Scheme that assigned this
                                 entry's species/type, e.g.
                                 "symmetry", "reduce", or
                                 "bispectrum".  Derives
                                 each fingerprint's native
                                 vs witness role (5.2.2):
                                 method M is native iff
                                 M == type_assignment.
  scf_threshold          float   SCF convergence
                                 threshold of the
                                 reference run.
  scf_iterations         int     Iteration count of the
                                 run that produced this
                                 potential.  Feeds the
                                 validation harness
                                 (5.8).

**Validation rules** (enforced at load time):

1. `schema_version` must equal 2 by the time the rest
   of these rules run.  A file that arrives at another
   version passes through the version gate (5.2.5)
   first: an older one is migrated forward when every
   field the intervening bumps added is derivable, and
   refused with an actionable message when one is not;
   a newer one is always refused.  The gate runs before
   rule 3 so that an out-of-date file is diagnosed as
   out of date rather than as missing a field.
2. `element_symbol` must match the parent directory
   name (case-insensitive).
3. Every required field must be present -- both the
   top-level keys (`schema_version`, `element_symbol`,
   `nuclear_z`, `nuclear_alpha`, `covalent_radius`)
   and the per-entry keys inside each `[[potential]]`,
   including `default`.  A missing field is a hard
   error with the file path, the entry label (when
   the missing field is per-entry), and the field
   name in the message.
4. `len(coefficients) == len(alphas) == num_gaussians`.
5. Labels must be unique within the file.
6. At least one entry with `label = "isolated"` is
   required.  The legacy baseline must always be
   present so that the validation harness (5.8) can
   compare against it on the same code path.  This
   rule is independent of rule 7: the `"isolated"`
   entry need not be the default entry.
7. Exactly one entry per file must carry
   `default = true`.  Zero or multiple defaults is a
   hard error.  The producer guarantees this: it
   tags the customized environment when the
   manifest declares one, otherwise it falls back to
   the `"isolated"` baseline (5.7 rule 7), so a
   hand-free element still loads with exactly one
   default.
8. Within any one `[[potential]]` entry's fingerprint
   array, the pair `(method, sub_spec)` must be unique.
   Two fingerprint records with the same `method` and
   the same `sub_spec` keys-and-values are a hard
   error.  Records with the same `method` but
   differing `sub_spec` are explicitly allowed (e.g.,
   one bispectrum fingerprint at `twoj1=8, twoj2=8`
   alongside another at `twoj1=6, twoj2=4`).
9. `method` must be one of the matcher names
   registered in `matchers.py` at load time (8.9).
   Unknown methods are a hard error rather than a
   silent skip, so that a typo in the manifest fails
   loudly rather than quietly omitting the fingerprint
   from the lookup.
10. The preferred flag is scoped to the **entry**, not to
    the file.  For each `[[potential]]` and each `method`
    appearing among *that entry's* fingerprint records,
    exactly one of them carries `preferred = true`.  Zero
    (a family present on the entry but none preferred) or
    two-or-more preferred for one `method` on one entry is
    a hard error.  An entry carrying no fingerprints at all
    -- the `"isolated"` baseline -- is exempt, having no
    family present to prefer.

    Additionally, and across the whole file: every preferred
    record of a given `method` must share one `sub_spec`.
    This is what makes the flag *mean* something -- it names
    the canonical `sub_spec` for that family, the settings
    the consumer computes its query with (5.6.5 step 2) --
    and two flagged records of one family disagreeing on
    settings would leave no canonical answer.

    The per-entry scope is forced by what the flag is for.
    Fingerprints are comparable only when computed at the
    same `sub_spec`, so the consumer reads any flagged
    record of a family to learn *which settings to use*
    (its payload is never read) and then compares its query
    against every entry.  The dedup (5.2.3) goes further:
    it keys on "the transferable descriptor every harvested
    entry shares," so it asks *each* entry for that entry's
    canonical bispectrum.  Both readings require every
    harvested entry to flag its own record, which is exactly
    what the producer writes (5.7: "the producer stamps the
    preferred flag onto each matching record").

    A file-wide "exactly one flagged record" rule would say
    the opposite, and would also make one arbitrary entry
    load-bearing: these files grow incrementally, entries
    are merged by the dedup and may one day be removed, and
    a canonical `sub_spec` recorded on a single entry would
    vanish with it, leaving every other entry's fingerprints
    uninterpretable.  Per-entry flagging is self-describing:
    read one entry alone and you know which of its records
    is canonical.

    The database-wide constraint that *all elements'*
    preferred records of a family share one `sub_spec`
    (manifest rule 11) is cross-file and cannot be checked
    from a single per-element file; the loader trusts the
    producer to have written a consistent database, in
    keeping with VISION Principle 5 (the database is
    produced from the manifest(s) by the producer, never
    hand-edited).

#### 5.2.1 Label naming

The `label` of every producer-harvested `[[potential]]`
entry is assembled mechanically from the OLCAO identity
coordinate of the harvested atom site, prefixed by the
reference solid it came from.  The canonical form is:

    <reference_id>-<element><species>-t<type>-a<site>

all lowercase.  Example: a single Si site in a diamond-Si
reference run yields `si_diamond-si1-t1-a1`.

The five components and their sources:

  Component     Source
  --------------------------------------------------------
  reference_id  The manifest `reference_id` of the
                reference solid (5.7).  The only
                human-minted part of the label.
  element       Periodic-table symbol of the site,
                lowercased (makeinput's element/species/
                type model: "elements are defined by the
                periodic table").
  species       OLCAO species number of the site
                (`atom_species_id`).  Species are defined
                by the structure of the system (the si1,
                si2 tags carried in the skeleton file).
  type          OLCAO potential-type number of the site
                (`atom_type_id`).  Types are defined by
                the needs of the calculation and assigned
                by the grouping pass (crystallographic
                equivalency, reduce, target, or block).
  site          The 1-based `atom_site` index of the
                harvested atom in the reference structure;
                equals `provenance.atom_site`.  The `a`
                prefix keeps it from colliding with the
                species token.

`element` and `species` fuse into the single OLCAO token
the CLI already speaks (`si1`), so a label reads straight
into a `-pot LABEL scope=si1` override (5.6.1) with no
translation step.

**Why these five and not fewer.**  The label is not the
environment encoder (the fingerprint is, 5.6.5) and not
the prose (the `description` is).  It is a typeable handle
plus an exact back-pointer: from the label alone a reader
reconstructs which run, which species, which potential
type, and which atom produced the entry.  Every component
but `reference_id` is a value the run already holds, so
the scheme invents no names and cannot drift -- the type
integer is the realized grouping verdict, not a separate
interpretation of the environment that could disagree
with the fingerprint.

**Type integers are per-run.**  The grouping pass relabels
and compresses type numbers on each run (the species/type
relax-and-renumber step in makeinput), so `t3` is
meaningful only within its `reference_id`.  The prefix
scopes it, so this is not a defect: `si_vac-si1-t2-a17`
reads as "type 2 as the si_vac run assigned it."

**Uniqueness is automatic.**  Two entries in one element
file share an element by construction; they differ in the
`-a<site>` component when they come from the same
reference solid (atom_site is unique per structure) and
in the `<reference_id>` prefix otherwise.  So rule 6
((element, label) unique) holds by construction provided
`reference_id` is manifest-unique (rule 5) and label-safe
(lowercase letters, digits, `-`, `_`; no spaces, since
the whole label is typed into `-pot`).

**Assembled at harvest, not authored in the manifest.**
The `type` number is not known until the grouping pass
runs, so the label cannot be written into the manifest
ahead of time.  The ordering of the producer pipeline
makes this a non-issue: every reference run passes through
makeinput, which assigns species (from the skeleton tags)
and type (from the grouping pass) while building the run's
Imago input, so by the time any run has finished executing
the `(species, type)` of every site is a settled fact.
The harvest stage runs strictly after all executions and
mints the storage label only then.

To hand the harvester those numbers without re-parsing the
input, the producer co-opts a file makeinput already
writes: `datSkl.map`.  Today that file records the mapping
between the sorted `imago.dat` atom numbering and the
original skeleton numbering (two columns, DAT# and
SKELETON#).  It is emitted at the exact point in the
output path where the sorted per-atom element, species,
and type arrays are in hand, so it gains three further
columns -- the site's element, `atom_species_id`, and
`atom_type_id`.  The `atom_site` of a manifest entry is a
skeleton-numbering index, so the harvester looks up the
SKELETON# row and reads that site's `(species, type)`
straight off the map.  The assigner records its own
verdict where the run lives; the harvester reads it back
and assembles the label.  Consequently the manifest entry's
`label` field (rule 3) becomes an optional curator override
of the derived default rather than a required field.  See
TODO C87 for the producer and makeinput changes.

**Reserved labels are unaffected.**  `"isolated"` (the
atomSCF baseline, rule 6) and `"default_solid"` (the
Phase-1 single-bulk improved entry) keep their fixed
names; the mechanical scheme governs only the
producer-harvested solid entries that Phase 2 adds.

#### 5.2.2 Native and witness fingerprints

A fingerprint plays one of two roles for the entry that
carries it, and the role is *derived*, not stored -- it
follows from how the reference run assigned its types.

- **Native.**  The run's species/type partition was
  computed *from* this fingerprint method; its grouping
  decided which atoms share the harvested potential (the
  environment-based species pass of 5.6.4).
- **Witness.**  The method did *not* drive the partition.
  It was computed only to record the geometry under a
  second descriptor, while the types were assigned by
  something else -- crystallographic symmetry, a
  position-based flag, or the *other* fingerprint method.

Method `M`'s role on an entry is exactly
`M == provenance.type_assignment`: native when the names
match, witness otherwise.  Deriving the role rather than
storing a per-record flag removes a field that could drift
out of agreement with `type_assignment`.

**Symmetric dual harvesting.**  Every harvested atom
records *both* registered fingerprint methods -- the
native one and the witness one -- regardless of which
assigned the types.  The witness computation is cheap next
to any SCF: the Python-side reduce shells are free given
the structure the run already wrote, and a bispectrum
(loen) pass is a single non-self-consistent Imago
invocation, small relative to the converged run that
produced the potential.  So there is no eager/deferred
asymmetry: a reduce-assigned run still computes the
bispectrum witness, and a bispectrum-assigned run still
computes the reduce witness.  Every stored environment is
thus dual-indexed, so a later run that assigns types by
*either* method can find a match.

**A witness is valid but approximate.**  The witness
fingerprint is faithful to its own atom -- it is that
atom's true descriptor under the second method.  What it
inherits from the assigning method is the *potential* it
points at: that potential was converged under the
assigning method's grouping, which may be coarser than the
witness method's own grouping would have been.  A future
query that matches a witness therefore imports a potential
built under a different partition.  This is safe because
the imported potential is only the *starting guess* for
the new run's self-consistent iteration (5.6.5, 5.6.6):
the new run sets its own partition and relaxes the guess
to the truth, so a witness import can cost convergence
speed, never correctness.

Symmetry-assigned entries are the limiting good case: when
crystallographic symmetry assigns the types
(`type_assignment = "symmetry"`), *both* methods are
witnesses, but exact ones -- symmetry-equivalent atoms
share one environment, so their reduce and their
bispectrum fingerprints are identical across the type and
the shared potential genuinely belongs to each.  The
witness coarseness above appears only in the disordered,
cross-method case.

**Producer convention: reduce is never the assigner.**
The native/witness machinery above is symmetric in
principle -- either method can assign -- but the producer
that populates this database (5.7) uses only two assigning
methods: crystallographic symmetry for ordered references,
bispectrum for disordered ones.  Reduce is always
harvested as a witness, never as the method that draws
species boundaries.  Two consequences follow.  The
bispectrum index is always complete: every environment
bispectrum can distinguish becomes its own entry.  And
reduce never shapes the *entry set*, so it is a pure
add-on column -- were the reduce scheme ever retired,
removing it is a clean deletion of the reduce fingerprint
records, with no entry to re-derive and no reference run
to repeat.  A database populated with reduce as an
assigner would not collapse so cleanly, which is why the
convention is fixed here rather than left to each run.

#### 5.2.3 Environment storage model: dedup and weights

The database stores **distinct environments, not atoms.**
This is the rule that keeps it from exploding as model
sizes grow into the tens of thousands of atoms, and it is
the substrate both the present nearest-neighbour lookup
and the future learned predictor (5.2.4) want.

**Why not one entry per atom.**  A naive harvest that
emitted one entry per atom would grow the database with
the *total atom count ever harvested*: a single
10,000-atom amorphous model would add 10,000 entries, and
the consumer's per-query cost -- (new-system atoms) x
(database entries) -- would grow with it.  But the atoms
of such a model are not 10,000 distinct environments; they
cluster heavily (that clustering is the whole premise of
grouping atoms into types).  Two atoms with near-identical
fingerprints also carry near-identical potentials, so the
second adds almost no information.  The quantity of value
is *coverage of environment space*, which is far
lower-dimensional than the atom count and **saturates**:
once the chemistry is covered, a new model of the same
material contributes almost no new environments.

**Dedup on insert, skip duplicates.**  Harvesting is
therefore an *insert-or-skip*.  When a harvested
environment matches one already stored, the producer does
not append a second entry and does not alter the existing
one: the first representative's potential stands and the
duplicate is dropped.  When the environment is new to the
database, it is appended.  The database therefore grows
with environment *diversity* -- sub-linear in atoms, and
saturating once the chemistry is covered -- rather than
with the atom count.

**Skip-on-match makes the build idempotent.**  Because a
match is skipped, re-running an unchanged manifest moves
nothing: every environment it would harvest is already
present.  The producer grows the database incrementally --
load the existing per-element file, harvest only this
batch's solids, append the environments new to it, and
write the file back -- so several manifests can accrete
into one database without any count drifting, because
nothing is accumulated on a match (5.7).  What this leaner
model gives up -- corroboration counts, a statistical mean
across duplicates, removal of a solid, and exact wholesale
rebuild -- is recorded under "Deferred" below.
The dedup tolerance is the producer-side mirror of the
consumer-side similarity floor (5.6.5, TODO C61): the
consumer asks "is this query close enough to a stored
environment to *import* its potential?", and the producer
asks "is this new environment close enough to a stored one
that storing it *adds nothing*?" -- one tolerance, two
uses.

**Dedup subsumes the symmetry gate.**  No special case is
needed for crystals.  A symmetry-assigned type's atoms are
identical, so they dedup to a single entry automatically --
exactly "one representative per type."  A disordered type's
atoms vary, so they dedup to however many genuinely
distinct environments exist, generally far fewer than the
atom count.  The crystalline case is just the degenerate
collapse of the same rule.

**Dedup keys on bispectrum, not all methods.**  Merging
happens at two scales.  *Within* one reference run the
assigning method -- symmetry for an ordered reference,
bispectrum for a disordered one -- defines the distinct
environments, so symmetry-equivalent or bispectrum-
equivalent atoms already collapse to one representative.
*Across* runs, where the per-element database actually
accumulates, the dedup keys on the **bispectrum**
descriptor at the preferred `sub_spec`: it is the
transferable one every entry carries (native for a
disordered reference, witness for an ordered one, 5.2.2),
whereas symmetry equivalence has no meaning between two
different structures.  Bispectrum subsumes the within-run
symmetry collapse too, since symmetry-equivalent atoms
carry identical bispectrum fingerprints.  Reduce never
gates the dedup -- it rides along only as a tolerated-
imprecise witness; keying on all methods would let reduce
split entries bispectrum considers identical, constraining
the entry count, the opposite of treating reduce as a
droppable column.  So the dedup asks one question, "are
these the same environment under bispectrum?", and the
surviving entry's reduce witness is one representative's,
accepted as approximate (5.2.2).

**The representative's potential is kept as-is.**  Two
atoms that dedup into one environment can still carry
*slightly different* harvested coefficients -- the
fingerprint characterizes geometry, but the converged
potential can differ in detail.  The leaner model does not
reconcile them: the entry's `coefficients` are the
harvested potential of the *first* atom that mapped to the
environment, stored verbatim.  This matches the consumer,
which already imports a single representative's potential
(5.6.5), so the lookup needs no averaging to work.
Averaging every mapped atom's potential into a mean, and
reporting the spread of those potentials as a confidence,
is a worthwhile refinement -- but it requires keeping
per-solid statistics, and is deferred (see "Deferred"
below).

**Search cost.**  Controlling the entry count via dedup is
the first-order fix.  Beyond it, the search is already
partitioned -- per element, and per `(method, sub_spec)`
-- so no query ever compares across elements or
incomparable descriptors, and fingerprint vectors admit
spatial indexing (a k-d / ball tree, or approximate
nearest neighbour) for sub-linear lookup once an index is
large.  Redundant per-method copies that the bundled shape
still carries (5.2.4) can additionally be collapsed into a
per-method index *in memory at load time*, so search
efficiency does not wait on any on-disk change.

#### 5.2.4 Forward compatibility and the learned predictor

The schema today is **bundled**: one entry carries one
potential, the reduce record, and the bispectrum record
together.  At scale a **normalized** shape is more
efficient: a pool of distinct potentials, plus a separate
deduplicated index per method that references them (so a
potential is stored once, and the coarser method's
fingerprints are not duplicated).  The migration from the
first to the second is **lossless for everything the
database is for**, provided two constraints -- already
satisfied -- are honoured:

1. **Retain type/observation identity.**  The label
   (`...-t<type>-a<site>`, 5.2.1) and provenance carry the
   type a potential belongs to, so the potential pool can
   be reformed by grouping on that identity.  This must
   not be dropped as an optimization.
2. **Dedup conservatively** (5.2.3).  Because the bundled
   shape keeps an entry whenever *any* method finds it
   novel, every method's coverage is present; normalizing
   merely collapses the redundancy each method carries.
   Had dedup discarded a method's novelty, normalization
   could not recover it.

Under those, the bundled form is a lossless *superset* of
the normalized form, and the conversion is a contained
schema-version bump confined to `initial_potential_db.py`
(ARCHITECTURE 8.7); producers and consumers above the
library barely change.  Deferral costs only a constant
factor of disk (duplicated potentials and coarser-method
fingerprints, bounded by environment diversity, not atom
count) -- not data, and not a scaling regression.

**The one thing normalization drops** is the *per-atom
pairing*: that one atom carried reduce fingerprint A *and*
bispectrum fingerprint B.  After the split you still see
that A and B map to the same potential, but within a
multi-atom type you cannot always tell which A paired with
which B.  That pairing is useful only for cross-method
correlation or a learned reduce/bispectrum translator,
both deprioritized as fragile, approximate paths.  If we
ever want it, carrying the originating entry's identity
(its label) onto each split record at normalization
preserves it at no schema cost -- so even this is
recoverable by choice, not lost by default.

**The learned predictor (a future consumer).**  A natural
endpoint is a neural network trained on the database to
predict a converged-quality potential directly from an
atom's environment, so that the starting guess becomes the
*answer* and the SCF iteration is shortened or skipped.
This is a *different* consumer from the nearest-neighbour
lookup, and it reshapes none of the above -- it *wants*
exactly the deduplicated, coverage-oriented corpus the
dedup model produces (and, once 5.2.3's deferred
statistical merge lands, the per-environment weights it
adds):

- At inference the network pays no per-query search cost;
  the lookup is amortized into its weights.  So the search
  concern of 5.2.3 does not apply to it.
- Its *training* is harmed by raw near-duplicates: they
  teach nothing new and *bias* the model toward
  over-represented environments (a network trained on raw
  per-atom data would be dominated by bulk-like sites and
  under-learn the rare defect environments that matter
  most).  A deduplicated set -- with the atom multiplicity
  of 5.2.3's deferred statistical merge available as a
  *sample weight* rather than baked in as duplication -- is
  the correct corpus.
- It may eventually want a *richer* descriptor than the
  coarse lookup key (possibly the raw local environment),
  so the fidelity retained per stored environment is worth
  keeping in mind -- a "do not paint ourselves into a
  corner" note, not a present decision.

Designing for dedup, coverage, and weights therefore
serves the nearest-neighbour present and the learned-
predictor future with one structure.  Should the predictor
become a committed goal rather than a possibility, it
should be promoted to a VISION-level objective; it is
recorded here as the dataspace's intended trajectory.

**Deferred: the statistical merge and exact rebuild.**
The leaner skip-on-match model keeps one representative per
environment and nothing more.  A later revisit (TODO C103)
may upgrade it to a *statistical merge* that folds every
atom mapping to an environment into the stored entry,
adding:

- a per-coefficient **mean** in place of the lone
  representative (averaging out fit noise across
  corroborating atoms) and a per-coefficient **standard
  deviation** for the spread of potentials that collapsed
  into the environment.  That spread is the empirical test
  of the database's founding premise that near-identical
  fingerprints carry near-identical potentials: a small
  spread says the assigning fingerprint is a faithful
  proxy for the potential and the mean is trustworthy; a
  large spread says the descriptor, or its tolerance, is
  collapsing potentials that genuinely differ, and is the
  signal to tighten the tolerance or distrust the entry.
  The same number gives the consumer a per-entry
  confidence and the learned predictor (5.2.4) a variance
  to weight by;
- an atom **multiplicity** (how many atoms, across all
  models, mapped here -- the natural sample weight a
  learned predictor wants) and a **model_count** (how many
  distinct reference solids contributed), which measures
  independent corroboration: an environment seen in five
  models is more transferable than one seen five-thousand
  times in a single model;
- the machinery that keeps both order-free and
  re-runnable: one small **contribution record per
  reference solid** (its `reference_id`, its `atom_count`,
  and the per-coefficient `coeff_sum` and `coeff_sumsq` of
  its atoms' potentials), from which the mean, spread, and
  counts are *derived*.  Keying the fold on the
  contributing solid is what lets the upgrade stay
  idempotent while also gaining the two properties
  skip-on-match gives up -- **removal** of a solid (drop
  its record, re-derive) and an **exact wholesale rebuild**
  to a byte-identical file.  Because coefficient-by-
  coefficient averaging is well-defined only on one shared
  set of alphas, the merge would assert the two entries'
  alpha sets are equal before combining (in the present
  regime the alpha set is a per-element constant, so the
  assertion only guards a future where alphas vary per
  environment).

**Implementation status.**  The minimal build (TODO C88)
loads the existing per-element file, harvests one
representative per distinct environment keyed on the
bispectrum descriptor at the preferred `sub_spec` (5.2.2),
and inserts-or-skips into it; it records the
`type_assignment` provenance field (5.2) so each
fingerprint's native/witness role is known, but none of
the statistical fields above.  The present interim
producer -- the C60 harvest -- predates even that: it
writes one entry per declared site and does not yet dedup.
The schema in 5.2 specifies the minimal target; the
statistical merge is the C103 follow-up.

#### 5.2.5 Schema versioning and migration

**The required-field set *is* the version.**  Adding a
required field, dropping one, or changing the type of one
is a schema version bump, and the bump is not optional
housekeeping -- it is what makes every other guard in this
section able to work.  Rule 1 (5.2) is the only check that
runs before any field-level check, so it is the only place
a reader can say something useful about an out-of-date
file.  It can only fire if the version actually moved.

The cost of skipping the bump is on record.  The
`type_assignment` field was added to the required Imago
provenance set (5.2.2) while `schema_version` stayed at 2.
A file written before that addition was therefore a
legal-looking v2 that failed a v2 required-field check,
and the reader could only report `missing required field:
type_assignment` -- naming a symptom, with no version
context and no recovery but hand-deletion of the file.
The rule above is what prevents a repeat.

**Four outcomes for the reader.**  Given a file's
`schema_version` and the version this build writes:

1. **Equal** -- load and validate as 5.2 describes.
2. **Older, and every field the intervening bumps added
   is derivable** -- migrate the parsed data in memory,
   stamp it with the current version, and load it.  The
   next `save` (5.5) writes the file forward, so the
   migration is paid once and silently: because the
   producer refreshes the `isolated` baseline on every run
   (5.7), any file it touches is rewritten current.
3. **Older, and some added field is not derivable** --
   refuse.
4. **Newer than this build writes** -- refuse, and say so
   in those terms.  The recovery here is the opposite of
   case 3's: the file is fine and the *code* is behind, so
   the user should update Imago rather than regenerate a
   database that a newer build wrote correctly.

**The honesty test for a migration.**  Each bump declares,
for every field it newly requires, either a *derivation*
from what an older file already carries, or that the field
is not derivable.  A default value standing in for an
unknown is not a derivation, and the distinction is the
whole point: the database's value rests on its provenance
being true, so a migrator that invents provenance destroys
what it was meant to preserve.

`type_assignment` is the worked example of not-derivable.
It names the scheme that drew the run's type partition,
and the native/witness role of every fingerprint in the
entry follows from it (5.2.2).  Nothing else in an older
file records which scheme ran.  Filling it with the
commonest value would leave a file written under `reduce`
or `bispectrum` silently mislabelled -- an error no
downstream check could catch, because the file would then
be perfectly well-formed.  So it is refused, not guessed.

The v1 -> v2 `default` flag is the contrasting case, and a
genuinely derivable one: rule 7 already documents the
fallback (tag the `isolated` baseline when no curator
choice is recorded), so a v1 file can be migrated to a v2
file that says exactly what a v1 file meant.

**The error contract.**  Every refusal names the file
path, the version the file carries, the version this build
writes, and a concrete recovery -- regenerate with
`build_initial_potentials.py`, or update Imago for case 4
-- and case 3 additionally names the field that blocked
the migration.  A bare missing-field message is never the
right report for a version problem.

**Where the machinery lives.**  The version table and its
derivations sit in `initial_potential_db.py` beside the
reader (ARCHITECTURE 8.7), which keeps a bump contained in
the one library, exactly as 5.2.4 assumes when it calls
the eventual normalization "a contained schema-version
bump".  A bulk `potential_migrate.py` -- rewriting a whole
installed tree in place, mirroring `guidance_migrate`
(ARCHITECTURE 10) and `resource_migrate` (ARCHITECTURE
11) -- is the natural companion, but it is a convenience
rather than a correctness requirement, because in-memory
migration plus the next save already carries forward every
file the producer touches.  It earns its keep only for
files no producer run will revisit.

No retroactive bump is proposed for the `type_assignment`
addition: the databases that predated it have since been
regenerated, so there is nothing left in the field for a
v2 -> v3 migration to rescue.  The rule starts from the
next change.

### 5.3 Sketch (gold, two entries with fingerprints)

The sketch below uses simplified float notation for
readability; the actual emitter writes 16 significant
digits per 5.5.  The Au example carries one curated
improved entry tagged `default = true`; the isolated
baseline carries `default = false`.  The
`default_solid` entry exposes two bispectrum
fingerprints at different `(twoj1, twoj2)` settings so
the same database file can serve calculations that
request either parameter pair.

```toml
schema_version  = 2
element_symbol  = "Au"
nuclear_z       = 7.9e+01
nuclear_alpha   = 4.0e-01
covalent_radius = 1.0e+00

[[potential]]
label         = "isolated"
default       = false
description   = "Single isolated Au atom (from atomSCF)."
num_gaussians = 32
alpha_min     = 1.0e-03
alpha_max     = 1.0e+02
coefficients = [
   1.2345678901234567e-03,
   ...  (32 entries, one per line)  ...
]
alphas = [
   1.0000000000000000e-03,
   ...  (32 entries, one per line)  ...
]

[potential.provenance]
source       = "atomSCF"
commit       = "abcdef1"
generated_at = "2026-05-08T14:00:00Z"

[[potential]]
label         = "default_solid"
default       = true
description   = "Au in fcc bulk (Fm-3m)."
num_gaussians = 32
alpha_min     = 1.0e-03
alpha_max     = 1.0e+02
coefficients  = [ ... ]
alphas        = [ ... ]

[potential.provenance]
source         = "Imago"
commit         = "fedcba2"
generated_at   = "2026-05-08T14:30:00Z"
reference_id   = "COD-1011098"
atom_site      = 1
kpoint_spec    = "12 12 12 0 0 0"
scf_threshold  = 1.0e-6
scf_iterations = 28

# The preferred = true record is the one the consumer uses
# for a file-dictated (crystalline) match (5.6.5 step 2); the
# preferred sub_spec for a family is uniform database-wide.
[[potential.fingerprint]]
method    = "bispectrum"
sub_spec  = { twoj1 = 8, twoj2 = 8 }
preferred = true
values    = [
   1.2345678901234567e-01,
   ...  (9 entries; length = twoj2 + 1 = 9)  ...
]

# A second sub_spec is freely stored (e.g. for validation
# comparison); it just may not also be preferred.
[[potential.fingerprint]]
method   = "bispectrum"
sub_spec = { twoj1 = 6, twoj2 = 4 }
values   = [
   ...  (5 entries)  ...
]

[[potential.fingerprint]]
method    = "reduce"
sub_spec  = { level = 2, thick = 5.0e-01, cutoff = 5.0e+00 }
preferred = true
shell_code.element = "au"
shell_code.levels  = [
   { distance = 2.88e+00, neighbors = ["au", "au", "au"] },
   { distance = 4.07e+00, neighbors = ["au", "au"] },
]
```

Reduce `shell_code` records the central atom's element
and, per reduction level, the shell distance and the
neighbor element symbols -- element-only, so the
descriptor transfers across structures (5.2).

### 5.4 In-Memory Representation

**Purpose of `initial_potential_db.py`.**  This is
the file-format **library**: a small, passive helper
module that knows exactly one thing -- how to read,
validate, look up entries in, and write per-element
`s_gaussian_pot.toml` files.  It contains no
orchestration, no SCF runs, and no curation logic.
Its only external dependency is `tomllib` (Python
stdlib).

It is imported by every other script in the chain:
`makeinput.py` for the runtime lookup,
`build_initial_potentials.py` for `save()` calls,
`bench_initial_potential.py` indirectly via Imago
runs that go through `makeinput.py`.  The
library / script split -- this module is the
**library**, `build_initial_potentials.py` is the
**producer**, `makeinput.py` is the **consumer** --
keeps read-only callers from pulling in manifest
handling or SCF-wingbeat code they don't use,
isolates any future schema-version bump or format
swap to one file (per ARCHITECTURE 8.7), and lets
unit tests cover the file format with synthetic
byte strings rather than real Imago runs.

The module's module-level docstring must capture
this purpose and role explicitly so future readers
do not have to reverse-engineer it from the call
sites.

The module exposes a small public surface:

```python
@dataclass
class FingerprintRecord:
    method:   str               # matcher name, e.g.
                                #   "bispectrum",
                                #   "reduce"
    sub_spec: dict[str, Any]    # method-specific
                                #   parameters (TOML
                                #   inline table)
    payload:  dict[str, Any]    # remaining record
                                #   fields; structure
                                #   is matcher-defined
                                #   (e.g., values for
                                #   bispectrum,
                                #   shell_code for
                                #   reduce)
    preferred: bool             # the family's consumer-
                                #   chosen record (5.6.5
                                #   step 2); <=1 true per
                                #   method per file

@dataclass
class PotentialEntry:
    label:        str
    default:      bool                       # rule 7
    description:  str
    num_gaussians: int
    alpha_min:    float
    alpha_max:    float
    coefficients: list[float]
    alphas:       list[float]
    provenance:   dict[str, Any]
    fingerprints: list[FingerprintRecord]    # may be
                                             # empty

@dataclass
class ElementDatabase:
    schema_version: int
    element_symbol: str
    nuclear_z:      float    # real: Z used as a real number
    nuclear_alpha:  float
    covalent_radius: float
    potentials:     list[PotentialEntry]
```

Public functions:

  Function            Behavior
  -----------------------------------------------------
  load(path)          Reads a TOML file via tomllib;
                      validates per 5.2; returns an
                      ElementDatabase; raises a clear
                      error on any rule violation.
  lookup(db, lbl)     Returns the PotentialEntry whose
                      label == lbl, or raises KeyError.
  baseline(db)        Returns the entry with label ==
                      "isolated"; guaranteed to succeed
                      by rule 6.  Used by the
                      validation harness (5.8), which
                      always compares against the
                      isolated-atom starting point on
                      the same code path.
  default_entry(db)   Returns the entry with default ==
                      true; guaranteed to succeed by
                      rule 7.  Used by makeinput.py
                      (5.6) when no scheme matches and
                      no `-pot LABEL` override is
                      given.  Distinct from baseline:
                      a curator may mark the isolated
                      entry default (then the two
                      functions return the same
                      object) or mark a curated
                      improved entry default (then the
                      two functions return different
                      objects).
  find_fingerprint(   Returns the FingerprintRecord on
    entry, method,    `entry` whose `(method,
    sub_spec)         sub_spec)` matches; raises
                      KeyError if absent.  Sub-spec
                      comparison is by canonical
                      key-and-value equality.  Used by
                      the matcher dispatch in
                      makeinput.py (5.6, 8.9).
  save(db, path)      Writes the file via the hand-
                      formatted emitter in 5.5.

`load` is read-only; `save` is the single point through
which the build pipeline (5.7) writes the database.
Other scripts (`makeinput.py`, the validation harness)
only call `load`, `lookup`, `baseline`, `default_entry`,
and `find_fingerprint`.

### 5.5 Hand-Formatted TOML Emitter

The emitter is deterministic at the bit level: given
the same `ElementDatabase` as input, it produces
byte-identical output bytes.  Emitter determinism is
non-negotiable because it cleanly separates
"formatting changed" from "numbers changed" in any
diff of a produced file.  It does **not** imply
file-level byte-identity across pipeline runs: the
build pipeline (5.7) refreshes provenance
timestamps every run, and SCF / fit numerical drift
(floating-point accumulation order, threading,
external library versions, and development changes
to the solver) can perturb the numerical content as
well.  The strict bit-level guarantee lives at the
emitter, not at the pipeline.

**Layout:**
- One blank line between the top-level keys block and
  the first `[[potential]]` block.
- One blank line between each entry's numerical body
  and its `[potential.provenance]` sub-block.
- One blank line between consecutive `[[potential]]`
  blocks.
- File ends with exactly one trailing newline.

**Scalar formatting:**
- Integers: bare digits; no underscores or signs.
- Floats: `"%.16e"` -- 16 significant digits in
  scientific notation.  This guarantees round-trip
  safety for IEEE-754 doubles.
- Strings: double-quoted; backslash-escape `"`, `\`,
  and control characters per TOML 1.0.

**Array formatting:**
- `coefficients` and `alphas` are emitted as multi-line
  arrays: opening `[`, one value per line indented by
  3 spaces, then closing `]`.  Trailing comma on every
  element (TOML 1.0 allows it) keeps diffs clean when
  arrays are extended.

**Key alignment within a block:**
- All `=` signs are vertically aligned at one space
  past the longest key name in that block.  The exact
  alignment column depends on which fields are present,
  but is fully deterministic for any given entry's
  data, satisfying idempotency.

### 5.6 Selection Algorithm

The Phase-2 selection algorithm unifies *species grouping*
and *manifest-entry pick* under one CLI: each
environment-based grouping flag drives both, with the
same parameters.  Position-based flags assign groupings
spatially; the manifest pick for those atoms falls
through to a manual `-pot LABEL` override or to the
default-tagged entry.  Types inherit potentials from
their parent species and are perturbed only by
electronic-state flags (XANES today).

#### 5.6.1 CLI surface

```
-pot LABEL   scope=SPEC (optional)
             scope=~SPEC (optional)
                        Manual override for the augmented
                        potential database.  Apply the
                        entry named LABEL to the atoms
                        selected by SPEC.  Without scope=,
                        LABEL applies across the whole
                        structure (the global form).  SPEC
                        is an element (scope=si -> all
                        silicon) or a species (scope=si1
                        -> species si1 only); scope=~SPEC
                        excludes that element or species.
                        Repeatable: different scopes may
                        carry different labels.  Optional.
                        See precedence in 5.6.3.

-nofingerprint          Disable the environment fingerprint
                        pick.  Force every species to its
                        default-tagged (isolated-atom) entry,
                        skipping the fingerprint match in
                        5.6.5.  A manual `-pot` override still
                        applies.  This is the opt-out for a
                        run that wants the plain isolated-atom
                        potential regardless of what
                        fingerprints the database carries --
                        the way `-reduce` ran before Phase 2.
                        Optional; the default is fingerprint
                        matching ON whenever the database
                        carries comparable records.

-target name=NAME ...   Position-based species/type
-block  name=NAME ...   grouping.  Existing flag
                        families gain a `name=NAME`
                        keyword (NAME matches
                        `[A-Za-z0-9_-]+` and is unique
                        across all spatial flags in
                        the run).  The name is the
                        handle environment-based
                        flags refer to via `scope=`.

-reduce ...  scope=NAME (optional)
             scope=~NAME (optional)
-bispec ...  scope=NAME (optional)
             scope=~NAME (optional)
                        Environment-based species
                        grouping AND manifest-entry
                        pick.  The same parameters
                        drive both: shell-code
                        comparison for `-reduce`,
                        bispectrum-vector comparison
                        for `-bispec`.  Optional
                        `scope=NAME` restricts the
                        scheme to atoms inside the
                        named spatial region;
                        `scope=~NAME` restricts to
                        atoms outside.  Without
                        `scope=`, the scheme applies
                        to the whole structure.

-xanes ...              Electronic-state flag.
                        Creates new types within
                        affected species (the
                        core-hole atom and its
                        in-sphere neighbors).
                        Unchanged from current
                        behavior; layered on top of
                        the species pass.
```

A note on `scope=`: the keyword is deliberately reused
across two option families, but its argument differs by
the pipeline stage the option acts in.  For the *grouping*
ops (`-reduce`, `-bispec`) `scope=` names a spatial region
declared by a `-target`/`-block name=NAME` (5.6.4) -- it
selects *which atoms to group*.  For the *assignment* op
(`-pot`) `scope=` names an already-resolved element or
species -- it selects *which assigned atoms receive the
label*.  Grouping precedes assignment, so by the time
`-pot` runs every atom already has an `(element, species)`,
which is exactly what its `scope=` refers to.  A future
relaxation may let `-pot scope=` also accept a named
spatial region for symmetry with the grouping ops; until a
concrete need appears it is element/species only (TODO).

**Retirement of legacy potential/basis substitution.**
The historical `-subpot` / `-subbasis` options substituted
an alternate *numbered* legacy file (`pot<N>` / `coeff<N>`,
`contract<N>.dat`) for a targeted element or species.  The
potential half is being removed: an audit of the installed
database (all 103 element directories) found only `pot1` /
`coeff1` present -- no element ever shipped a `pot2` or
higher, so `-subpot` had nothing to substitute and was
never exercised in practice.  Its capability is fully
subsumed by the augmented database plus the now-scoped
`-pot`: an alternate potential is a *labeled entry* the
curator harvests (5.7), selected per element or species by
`-pot LABEL scope=SPEC` rather than by a magic file index.
The basis half (`-subbasis`) follows the same trajectory
once the basis-set database gains an augmented, labeled
form analogous to the potential database; the two
deprecate together, preserving their long-standing
symmetry.

#### 5.6.2 Mutual exclusion

`-reduce` and `-bispec` are mutually exclusive: at most
one environment-based scheme per run, regardless of
scoping.  Both present is a hard error at CLI parse
time.  (A future relaxation may permit disjoint scopes
to host different environment schemes; for now the rule
is global.)

Multiple `-target` and `-block` flags compose freely
in their existing order-dependent way, with `name=` as
a new but additive field.  Multiple environment-scheme
flags of the same type are also disallowed: at most one
`-reduce` *or* one `-bispec`.  Multiple invocations
would compete for the same species partition.

#### 5.6.3 Per-element preflight

Before the per-atom work, for each unique element symbol
in the parsed structure:

1. Construct the path
   `share/atomicPDB/<elem>/s_gaussian_pot.toml`.
2. If the file exists, call
   `initial_potential_db.load(path)`.  On any validation
   error, abort with a message naming the file and the
   failing rule.
3. If the file does not exist, fall back to the legacy
   `pot1` / `coeff1` reader path for this element.
   Emit an info-level message that the augmented
   database is not yet populated for this element.  The
   element's atoms will participate in species grouping
   normally but cannot match any fingerprint scheme;
   they receive the legacy isolated-atom potential.
4. **Coverage note** (whenever fingerprint matching is
   enabled -- the default; `-nofingerprint` turns it off).
   Check whether the loaded database carries at least one
   fingerprint record for *any* registered matcher.  If none
   does, do **not** abort: emit an info-level message naming
   the element, and let its atoms fall through to the
   default-tagged entry (the isolated-atom potential) at the
   per-species pick (5.6.5 step 3).  The fingerprint pick is
   *decoupled* from how species are grouped (5.6.4, 5.6.5):
   it runs for every species regardless of the grouping
   scheme, but it is still a bonus layered on top of grouping
   that the default entry always backstops, so an element
   with no stored fingerprints simply takes the default
   potential -- exactly as a run did before Phase 2.  Aborting
   here would be wrong: it would break an otherwise-valid run
   merely because the database is not yet populated for that
   element, which is the *current* state of every shipped
   element (none carries fingerprint records yet -- see the
   producer side-quest in TODO C91), and which is never the
   user's error.  The info note keeps the fall-through visible
   -- a user who expected a fingerprint match is told why
   every species received the same default potential -- without
   making it fatal.  The note is suppressed under
   `-nofingerprint`, where the default is the deliberate
   choice.  (The bispectrum *grouping* path still reports its
   own loen-coverage condition inside `makegroups.py` (5.10);
   that is a separate report about grouping, not about this
   potential-pick coverage.)

#### 5.6.4 Species pass

Walk `settings.methods` in CLI order (same dispatch as
today's `assign_group`).  For each flag whose `op` is
`"species"`:

- **Position-based** (`-target name=N`, `-block name=N`).
  Apply spatial grouping as today: atoms passing the
  spatial test receive new species IDs; atoms outside
  keep their current assignment.  The named region is
  remembered for any subsequent environment-scheme
  scope reference.
- **Environment-based, reduce only** (`-reduce
  scope=N`).  Determine the atom set:
    - No `scope=`: all atoms of the active element.
    - `scope=N`: only atoms inside the spatial region
      named `N`.
    - `scope=~N`: only atoms outside the named region.
  Compute the per-atom reduce fingerprints over that set
  in-Python, then bucket atoms whose vectors agree within
  tolerance into species (the matcher's `distance` and
  `default_similarity_floor`, ARCHITECTURE 8.9).  Atoms
  outside the scope keep whatever assignment earlier flags
  produced.

  **Bispectrum is not grouped here.**  A Fortran-side
  descriptor can only come from a completed Imago run, so
  bispectrum grouping happens *before* makeinput, in
  `makegroups.py` (5.10): it runs the loen sequence and
  rewrites the skeleton with explicit per-element species
  tags.  By the time this species pass runs, those atoms
  already carry their bispectrum-derived species (read from
  the skeleton, structure_control `use_file_species=True`),
  so makeinput does nothing special for them -- it is a
  plain input-writer.

After the species pass, every atom has a final
`atom_species_id[atom]`.

**Grouping and the potential pick are independent.**  The
species pass decides only *how atoms are partitioned* into
species.  Which database entry each species then receives is a
separate step (5.6.5) that runs for *every* species, no matter
how it was grouped -- crystallographic default, position-based
(`-target`/`-block`), or environment-based (`-reduce`,
`-bispec`).  In particular, the per-atom reduce descriptors
computed here for *bucketing* are not the only road to a
fingerprint match: a crystallographically grouped species,
which the species pass never touched with a matcher, is still
eligible for the fingerprint pick, because 5.6.5 computes
whatever query descriptor it needs at pick time from the
species' own atoms.  Phase 2's first cut wired the pick to the
reduce grouping's descriptors and so matched only when
`-reduce` was active; C93 decoupled the two, making the pick
the default for any grouping (with `-nofingerprint` the
opt-out).

**Where the representative comes from depends on the system.**
The order of operations is always the same -- *group first,
then summarize each group into one representative, then match
the representative against the database* (5.6.5 step 2) -- but
what performs the grouping differs by system, and that decides
whether the representative is fresh work or reused:

- *Crystalline or position-grouped systems.*  Crystallography
  (or a `-target`/`-block` flag) defines the species without
  any environment matcher.  No per-atom environment descriptor
  was computed during grouping, so the pick computes a fresh
  representative for each species from its member atoms -- one
  query per matcher family it tries.

- *Amorphous and other non-crystalline systems.*  There is no
  symmetry to define the species, so an environment scheme
  (`-reduce` or `-bispec`) must run *first*, precisely to do
  the grouping.  Only after that grouping exists can a
  representative be formed for each group.  The per-atom
  descriptors the scheme computed for bucketing are exactly
  what that representative is built from, so for the family
  that did the grouping the descriptors are computed once and
  reused -- the grouping pass and the database match share one
  descriptor computation rather than repeating it.  When the
  user invoked the scheme explicitly, that choice is *honored*:
  the user's `sub_spec` drives both the grouping and the match,
  and if the database happens to carry no fingerprint at that
  `sub_spec` the species simply misses and takes the default
  entry (5.6.5 step 3).  This is best-effort, never an error --
  the database does not overrule the user's choice, and a user
  who deviates from the curation convention (5.7) accepts that
  the stored fingerprints may not match.  In normal use the
  user runs at the convention `sub_spec`, so the reuse is exact
  and matches are dense.

Which descriptor family the match uses depends on the regime,
and the bispectrum-then-reduce priority applies only to the
file-dictated case -- it never overrules an explicit user
choice:

- *User chose a scheme* (`-reduce` or `-bispec`).  The match
  uses *that* family at the user's `sub_spec`; the other family
  is never run behind the user's back.  So `-reduce` on an
  amorphous model matches against reduce fingerprints only,
  even when the database also carries bispectrum ones.

- *Species are file-dictated* (crystalline or pre-assigned, no
  environment flag).  The consumer has no user `sub_spec`, so
  the *database* decides: it reads the `preferred` record per
  family (the 5.7 convention -- one preferred `sub_spec` per
  family, uniform database-wide) and uses bispectrum when the
  database carries a preferred bispectrum record, otherwise
  reduce.  Exactly one family, one `sub_spec`, one query -- the
  match is a simple best-effort lookup, never a search across
  multiple sub_specs or a per-element tuning exercise.

#### 5.6.5 Manifest-entry pick per species

For each `(element, species)` pair appearing in the
final assignment, pick exactly one `PotentialEntry`
from that element's database.  The pick runs for *every*
species, independent of how the species was grouped
(5.6.4) -- the C93 decoupling, so a crystallographically
grouped species is as eligible for a fingerprint match as
a `-reduce`-grouped one.  Precedence, top to bottom:

1. **`-pot` (manual override).**  Each `-pot LABEL`
   applies its label to the species in its `scope=`
   (an unscoped `-pot LABEL` applies to every species in
   every element; a `scope=SPEC` form applies only to the
   named element or species).  When several `-pot` flags
   are given, a more specific scope wins over a broader one
   for the species they share (species beats element beats
   global); two equally-specific scopes naming the same
   species are a hard error at parse time.  `KeyError` for
   any in-scope element that lacks the label is fatal --
   `-pot` is a deliberate override and a silent fallback
   would mask user intent.
2. **Fingerprint match** (enabled by default; disabled for
   the whole run by `-nofingerprint`, and skipped for a
   loen-descriptor build -- both send every species straight to
   step 3, the loen-build case for the reason below).  Match is
   a single best-effort lookup,
   not a search: it picks exactly one descriptor family and one
   `sub_spec`, computes one query, and accepts a miss.  Which
   family and `sub_spec` are used depends on the regime (5.6.4):

   - **User chose a scheme** (`-reduce` or `-bispec`).  Match
     uses *that* family at the user's `sub_spec`.  The per-atom
     descriptors already computed for grouping are reused for
     the representative.  The database is not allowed to
     overrule the choice: if it carries no fingerprint at that
     `(method, sub_spec)`, the species misses and falls to step
     3 -- silently, since a deviation from the convention (5.7)
     is the user's call, not an error.

   - **Species are file-dictated** (crystalline or
     pre-assigned).  There is no user `sub_spec`, so the
     database decides: read the `preferred` record per family
     (the 5.7 convention -- one preferred `sub_spec` per family,
     uniform database-wide), and use bispectrum when the
     element's database carries a preferred bispectrum record,
     otherwise reduce.  Compute that one query -- the bispectrum
     query Fortran-side via the loen path (fast), the reduce
     query in-process from geometry.  No loen run is triggered
     for an element whose database carries no preferred
     bispectrum record.

   A **loen-descriptor build skips this match entirely.**
   makeinput builds the loen input itself -- the provisional
   `imago.dat` a loen run reads to compute a bispectrum
   descriptor (5.10.2) -- by running makeinput with `-loeninput`.
   That build must not perform the match: in the file-dictated
   regime the match may need a bispectrum descriptor, which is
   computed by a loen run, whose own first step is exactly this
   build.  A match here would therefore invoke the build within
   itself, without end.  So a `-loeninput` build always takes the
   default-tagged entry (step 3).  This is not a compromise --
   the bispectrum is geometric, so the potential the loen run
   sees is irrelevant to the descriptor it produces (5.10.2) --
   and it closes the one path by which potential resolution could
   reach back into itself.  This skip is the third condition that
   sends a species to step 3, beside `-nofingerprint` and a
   `-pot` override.

   In both regimes, once the family and `sub_spec` are fixed:
   ask the matcher to summarize the species' atoms into one
   representative fingerprint via its `representative` method
   (8.9); among the element's entries carrying a comparable
   `(method, sub_spec)` fingerprint, pick the one minimizing
   `distance(representative, entry_fingerprint)`; accept it
   only if that distance is within the matcher's
   `default_similarity_floor`, otherwise warn (naming the
   species and the best-but-rejected entry) and fall to step 3.
   The matcher chooses representative semantics appropriate to
   its descriptor space: `BispecMatcher.representative` returns
   the element-wise mean of the member vectors, `ReduceMatcher`
   returns the first member's shell-code (symmetry-equivalent
   or reduce-grouped atoms carry identical shell-codes by
   construction, so any member is exact); future matchers may
   use a medoid or another scheme -- the protocol does not pin
   the choice.

   Each matcher carries a heuristic default for its similarity
   floor; the concrete numbers (e.g., 0.05 for `ReduceMatcher`,
   0.10 for `BispecMatcher`) are starting values intended to be
   tuned during the Phase-2 validation pass (TODO C61) against
   the benchmark systems' actual fingerprint-distance
   distributions, and users may override them per scheme on the
   CLI when a particular system warrants a tighter or looser
   tolerance.
3. **Default tag.**  For any species not assigned by steps
   1-2 — `-nofingerprint` in force, no fingerprint family
   matched within its floor, or the database carried no
   comparable fingerprints — use `default_entry(db)` for that
   element.  Guaranteed to succeed by rule 7.

**Why one representative, not a per-atom pick.**  A species
groups atoms whose environments are interchangeable to within the
matcher's tolerance, yet every atom of the species must receive
the *same* initial potential.  The tempting question -- "matched
against *which* atom?" -- is a trap.  Singling out the first atom
(or any fixed atom) makes the chosen potential depend on the order
the atoms happen to appear in, which for a non-crystalline cell is
arbitrary; renumbering the atoms would then change the result,
which is unacceptable.  Collapsing the whole group into one
order-independent representative removes the question entirely --
no atom is privileged, and the answer is stable under renumbering.

The two descriptor families sit in genuinely different positions
here.  Reduce fingerprints match as a hard yes/no, so every atom
in a reduce species carries an *identical* shell-code by
construction; "which member speaks for the group" is moot and the
first member is exact.  Bispectrum fingerprints match by
closeness, so the members of a bispectrum species scatter around a
centroid; the element-wise mean is the order-independent summary
that speaks for all of them.  Should that scatter ever grow wide
enough that the choice of representative changes which database
entry is matched, the correct response is to tighten the scheme's
similarity tolerance and split the species -- not to smooth over a
loose cluster.

We deliberately do *not* search for the single database entry that
best fits *all* of a species' atoms at once (an all-to-all
comparison).  Such a search is affordable -- a species holds tens
of atoms, not the whole system -- but it is the wrong tool on two
counts.  First, it hides exactly the loose-cluster signal above,
quietly returning a least-bad compromise where a split was called
for.  Second, the initial potential is only a *starting guess*:
the self-consistent iteration relaxes it to the true potential
regardless of where it began, so optimizing a launch point the
calculation is about to iterate away is poor return on complexity.
A representative of a tight cluster is good enough by the only
standard that matters here -- fast, reliable convergence.  The
medoid (the member collectively closest to the others) is the
robustness fallback should a family's clusters ever prove skewed
by an outlier; the protocol leaves the choice to each matcher
(8.9) rather than pinning it in this section.

The chosen entry is attached to the species and
flows through to every type born from it (5.6.6).

#### 5.6.6 Type pass and electronic-state perturbation

Types are subdivisions of a species made on electronic
grounds, not geometric grounds.  They are not driven by
the matcher.  Algorithm:

1. Start the type pass with one type per species: every
   species' single child type is type 1.  Each type
   inherits its parent species' chosen
   `PotentialEntry`.
2. Apply each electronic-state flag in CLI order.
   Today this is just `-xanes`, which:
   - Splits affected species into a new type for the
     core-hole atom plus a new type for its
     in-sphere neighbors.
   - Assigns those new types a XANES-specific
     potential (the existing core-hole machinery,
     unchanged in Phase 2).
3. Atoms unaffected by any electronic-state flag keep
   their species' inherited potential.

From Imago's perspective only the total type count
matters: the per-type potential list emitted into the
Imago input file is the flattened
`(element, species, type) -> chosen_potential` mapping.

#### 5.6.7 Emit

For each Imago-level type, format and write into the
Imago input file, in today's on-the-wire format: the
element file's top-level `nuclear_z`, `nuclear_alpha`,
`covalent_radius`, plus the chosen entry's
`num_gaussians`, `alpha_min`, `alpha_max`,
`coefficients`, and `alphas`.  Imago itself is unaware
of the manifest, the fingerprint records, or the
matcher — it sees only the resolved per-type numbers.

### 5.7 Build Pipeline Algorithm

**Purpose of `build_initial_potentials.py`.**  This
is the script that **produces** the augmented
per-element database files
(`share/atomicPDB/<elem>/s_gaussian_pot.toml`).  It
takes a curated set of reference solids as input,
**delegates each solid's SCF run to kaleidoscope**
(section 6.2) as a small *predict-then-verify*
flight, harvests the converged potential at each
named atom site, and writes the results into
per-element database files via
`initial_potential_db.save()`.  It is the
**producer** half of the library / producer /
consumer split documented in 5.4, and a *client* of
kaleidoscope in the sense of ARCHITECTURE 9.7: it
decides *what* to compute, while kaleidoscope owns
*running*, *caching*, and *tracking* the batch.

**The curation manifest -- what it is, what it
contains, why.**  The pipeline's primary input is
a human-readable **curation manifest**.  It
declares *which solids, under which labels, with
which SCF settings*.  A build may take more than
one manifest -- each is a batch of reference
solids, and the curated set is their union (the
database grows incrementally, 5.2.3) -- so no
single file need ever list every solid.  A
manifest's job is threefold:

1. **Declare the curation set.**  The manifest
   *is* the curation strategy made explicit.
   Adding a reference solid means adding a
   manifest entry and running the producer;
   reviewing the curation set means reading the
   manifest(s).  (Removing a solid from the
   database is a deferred capability, 5.2.3 --
   dropping its manifest entry stops re-adding
   it but does not retract what it already
   contributed.)
2. **Tell the pipeline what to harvest.**  The
   harvest is automatic -- one representative per
   distinct environment (5.2.2, 5.2.3) -- so the
   manifest need not enumerate atom sites.  An
   entry, when present, only *annotates* an
   auto-discovered environment (its default tag,
   its description, or a pinned representative).
3. **Record the SCF settings used.**  k-points,
   convergence threshold, etc., recorded into
   the provenance fields of 5.2 so every
   database entry carries the conditions of its
   reference run.

VISION Principle 5 ("the database must be
regeneratable from the curated set, not a
hand-edited artifact") rules out hardcoding the
curation set inside the pipeline script and rules
out folder-of-files conventions that lose
metadata.  The manifest is the smallest piece of
structured data that closes the gap: every
curation choice captured in version-controlled
files alongside the structure files they point
at, so growing the database stays a scripted
operation on (manifests, structure files, Imago
build), not a hand-edit.  (An *exact* wholesale
rebuild from those inputs is a deferred
capability, 5.2.3; the leaner build grows the
database in place and treats the database file
itself as the artifact to keep.)

The producer grows the database with an *idempotent*
insert-or-skip (5.2.3): it loads each affected element
file that already exists (seeding a fresh `"isolated"`
baseline from the current `pot1`/`coeff1` only when no
file is there yet), harvests this run's solids, appends
the environments new to the file, and skips those already
present.  Re-running an unchanged manifest therefore moves
nothing, and a manifest of new solids folds in only their
new environments -- so the build is safe to repeat and
safe to extend, one batch at a time, without any count
drifting.  What it does *not* do is reconcile a duplicate
(the first representative stands) or retract a solid; a
statistical merge across duplicates, removal, and an exact
wholesale rebuild are the deferred upgrade of 5.2.3.

**What it contains** (per the schema sketched
below): a database-wide `[characterization]`
recipe (the fingerprint sub_specs, 5.2 rule 11);
per-solid fields the SCF run needs
(`structure_path`, `kpoint_spec`,
`scf_threshold`) plus a stable `reference_id` and
an optional `source_description`; and optional
per-entry customizations (`default`, `description`,
`label`, or a pinning `atom_site`) that annotate
auto-discovered environments -- the harvest takes
one representative per environment (5.2.2) and
needs no enumerated sites.

**What it deliberately omits**: the numerical
potentials themselves (those are SCF outputs);
iteration counts and convergence metrics (those go
in `share/curation/run_log.toml`); any executable
logic.  The manifest is data driving a script,
not behavior.

**Build analogy.**  The manifest is the **build
configuration** for the database.  The structure
files are the source.  The Imago build is the
toolchain.  `build_initial_potentials.py` is the
build script.  The augmented database is the
compiled output.  Same role `pyproject.toml`
plays for a Python package, or a Makefile plays
for a binary.

The script's module-level docstring must carry
this purpose, manifest rationale, and build
analogy explicitly so future students reading the
source can build the mental model from the file
itself.

**Where the schema lives, and how a manifest is
authored.**  The manifest schema -- the dataclasses
a parsed manifest becomes, the strict reader
`load_manifest_v2`, the relaxed structure-only
reader `load_structure_sources`, and the writer
`format_manifest` -- lives in its own neutral leaf
library, `curation_manifest.py`, so the producer and
the authoring tools share one definition: the
producer imports it to *run* a manifest, while
`cod_fish.py` and `expand_manifest.py` import it to
*write* one (and to source the default recipe and run
settings), and the library itself depends only on the
lower libraries it validates against
(`initial_potential_db`, `guidance_db`).  A curator
does not hand-write a complete manifest from nothing:
`cod_fish.py` discovers and pins structures from COD
and, by default, prints a *complete, runnable*
manifest -- a `schema_version` header, the shared
`[characterization]` recipe and `[defaults]` run
settings (their values taken from the shared library,
not hardcoded in cod_fish), and one
`[[reference_solid]]` stub per structure, each with a
`reference_id` auto-derived from the CIF metadata
(`<formula>_<symbol>_<number>_<year>`).  Its
`--sketch-only` mode instead prints the bare stubs
(no recipe, no defaults), for the cases that still
want `expand_manifest.py`: a manifest of local
(non-COD) structures, a non-default recipe, or
interactive per-structure customization.
`expand_manifest.py` reads such a sketch -- from a
file, or from standard input when none is named -- and
fills in the `[characterization]` recipe and the
`[defaults]` run settings (the same shared-library
values); entry customizations are optional, added
interactively or left for the curator to write by
hand.  Each sketch stub also carries discovery hints
cod_fish read from the CIF.  The composition
(``elements``) is a non-schema hint the producer
ignores and the finished manifest omits; the
``source_description`` is persisted as a
`reference_solid` field, from which the harvest
composes each environment's description (5.2.1).  The
writer emits human-readable TOML -- shortest
round-trippable floats, inline `sub_spec` tables in
their authored order, `label` only when present and
`preferred` only when true -- and its output
round-trips through `load_manifest_v2`.

**Reproducibility (layered).**
`build_initial_potentials.py` reproducibility is
layered, with each layer's guarantee matched to what
the inputs can actually support:

- **Emitter determinism (bit-level, strict).**  The
  TOML emitter (5.5) writes byte-identical bytes for
  a fixed in-memory `ElementDatabase`.  Formatting,
  key order, and numeric printing never themselves
  introduce diff churn.
- **Pipeline numerical output (precision-level,
  loose).**  Given the same manifest, same `pot1` /
  `coeff1` files, and same Imago build, the
  pipeline's numerical outputs (`coefficients`,
  `alphas`, `num_gaussians`, `alpha_min`,
  `alpha_max`) should agree across runs at the
  precision the SCF / Gaussian-fit chain can reach.
  Bit-identity is not promised: floating-point
  accumulation order, threading, and external
  library versions can shift the last few bits.
  Development changes to the SCF or fit code that
  legitimately produce better-converged potentials
  are expected and welcome -- the database is a
  build artifact the producer grows, not an
  archival hand-curated table.
- **Provenance metadata (free).**  Timestamps,
  commits, and similar fields refresh on every run
  and are exempt from any reproducibility
  guarantee.

**Inputs:**

- A curation manifest (TOML; schema v2 specified below).  Default
  location is `share/atomicBDB/manifest.toml`; alternate via the
  `--manifest` flag.
- The historical guidance dataspace (`share/historicalGuidanceDB/`,
  DESIGN 7).  Consulted per reference solid to predict the
  converged k-point density and to size the verification grid
  around that prediction.  When the dataspace cannot predict for a
  solid (too few neighbors), the producer falls back to the wide
  default grid defined in the flight builder (DESIGN 6.2.8).
- The existing `share/atomicPDB/` tree -- the `pot1` and `coeff1`
  files (read when refreshing each element's `"isolated"` baseline
  entry) and any existing `s_gaussian_pot.toml` (loaded so the
  build grows it incrementally rather than replacing it, 5.2.3).
- An Imago build location.  The producer does not invoke it
  directly; kaleidoscope drives it through the wingbeat seam when
  dispatching the verification flights (and the follow-on
  `imago.py -loen -scf no` runs for Fortran-side fingerprint
  matchers).
- Network access to the Crystallography Open Database (COD).  Used
  only by the structure-materialization step (below) for
  `[[reference_solid]]` entries that declare a `cod_id`, to fetch
  the pinned revision once to a local file.  `structure_path`
  entries need no network.

**Manifest schema (version 2).**

The manifest is TOML.  Schema v2 adds two fields per
`[[reference_solid.entry]]`: a required `default`
boolean (which entry becomes the database file's
default-tagged entry per per-element-database rule 7),
and an optional list of `[[reference_solid.entry.fingerprint]]`
records declaring which fingerprints to harvest
during the reference run.  The Phase-1 v1 manifest is
not loaded by the v2 reader; curators must add the
`default` tag explicitly during the bump.

```toml
schema_version = 2

# The characterization recipe is a database-wide constant: one
# preferred sub_spec per method, applied to every harvested
# environment (5.2 rule 11).  Declared once here, never repeated
# on an entry.
[characterization]
  [[characterization.fingerprint]]
  method   = "bispectrum"
  sub_spec = { twoj1 = 8, twoj2 = 8 }
  [[characterization.fingerprint]]
  method   = "reduce"
  sub_spec = { level = 2, thick = 0.5, cutoff = 5.0,
               tolerance = 0.05 }

# The run settings shared by every reference solid -- declared
# once here and inherited by each [[reference_solid]] that does
# not override them (validation rule 2).
[defaults]
kpoint_spec        = { density = 60.0, shift = [0.0, 0.0, 0.0] }
scf_threshold      = 1.0e-6
basis              = "fb"
functional         = "wigner"
kpoint_integration = "linear-tetrahedral"

# Harvest settings shared by every reference solid.  Unlike the
# [defaults] run settings above, these do not feed any single
# calculation -- they govern how the finished runs are read back
# into the database.  Optional block: a solid inherits these when
# it does not name its own, and an omitted setting falls back to
# the producer's built-in value (validation rule 2).
[harvest]
kpoint_convergence_threshold = 1.0e-3   # eV/atom; 1 meV/atom

[[reference_solid]]
reference_id          = "au_fcc"
system_type           = "crystalline"
# Exactly one of (cod_id, cod_revision) or structure_path is set
# per [[reference_solid]] (validation rule 4 below).
cod_id                = 9008463
cod_revision          = "2023-04-12"
# structure_path      = "au_fcc.skel"     # alternative form
# Run settings inherit from [defaults]; name any of kpoint_spec /
# scf_threshold / basis / functional / kpoint_integration here
# only to override this one solid's value.
source_description    = "Au in fcc bulk (Fm-3m), COD 9008463."

  # Entries are OPTIONAL customizations on the auto-harvested
  # environments, not the harvest list.  This one designates the
  # element default and gives a hand-written description; an
  # environment no entry mentions is still harvested, taking
  # default = false and an auto-composed description.
  [[reference_solid.entry]]
  element     = "Au"          # optional; cross-checked at harvest
  atom_site   = 1             # optional handle: a representative
                              #   atom of the target environment
  label       = "default_solid"
  default     = true
  description = "Au in fcc bulk (Fm-3m)."

    # A per-entry fingerprint is a RARE override: an extra,
    # non-preferred sub_spec to harvest for this environment
    # alongside the database-wide preferred ones above.  The
    # common entry carries none.
    [[reference_solid.entry.fingerprint]]
    method   = "bispectrum"
    sub_spec = { twoj1 = 6, twoj2 = 4 }

[[reference_solid]]
# ... another solid ...
```

**Per-solid fields.**  The five run settings -- `kpoint_spec`,
`scf_threshold`, `basis`, `functional`, `kpoint_integration` --
may be given per solid or omitted to inherit the `[defaults]`
block (below); each must be *resolvable* one way or the other
(rule 2).  The one harvest setting,
`kpoint_convergence_threshold`, may likewise be given per solid
or omitted to inherit the `[harvest]` block -- but it carries a
built-in default, so unlike the run settings it need not be
resolvable from the manifest at all.  The rest are per solid.

- `reference_id` (string): stable, human-readable identifier for
  the reference solid.  Used as the kaleidoscope flight/run stable
  id (DESIGN 6.2.4) and as the filename stem of the materialized
  local structure; must match `[A-Za-z0-9_-]+` and be unique
  across the manifest.
- `system_type` (string): one of `"crystalline"`, `"amorphous"`,
  `"nanostructure"`, `"molecular"`.  Required field: the producer
  must declare this so the guidance-dataspace predictor (DESIGN
  7) can switch to the correct sub-model when called with the
  reference solid's feature vector.  The bulk of curation
  entries are crystalline (the canonical reference-solid case);
  amorphous and molecular references are rare in this manifest
  but supported for completeness.
- `source_description` (string, optional): a one-line, structure-
  level description of the reference solid, persisted from the
  hint cod_fish reads from the CIF.  The harvest composes each
  environment's auto-description from it, qualified by the species
  and site it discovers (5.2.1); a curator entry customization can still
  replace any environment's description by hand.
- `cod_id` (positive integer, optional iff `structure_path` is
  set): Crystallography Open Database entry ID.  The
  structure-materialization step fetches the structure once, at
  the pinned `cod_revision`, to a local file before the solid's
  flight is dispatched.
- `cod_revision` (non-empty string, required iff `cod_id` is set):
  pinned COD revision token (an ISO date or COD-supplied revision
  identifier).  Pinning keeps the build deterministic against
  upstream COD edits: re-running months later fetches the same
  bytes by construction.
- `structure_path` (string, optional iff `cod_id` is set):
  relative path under the manifest's directory.  Hand-authored
  escape hatch for materials not in COD (a new polymorph, a
  hypothetical structure, an unpublished result).
- `kpoint_spec` (inline table): k-point mesh specification using
  Imago's k-point style code 2 (minimum-density mode; see
  `src/imago/kpoints.f90:60-67`).  Fields:
    - `density` (real, units Bohr^-3): `minKPointDensity`, the
      minimum number of k-points per unit reciprocal-space
      volume.  Imago picks the per-axis counts needed to meet
      this density given the lattice.  Style 2 over style 1
      (explicit axial counts) because style 2 stays correct when
      the lattice is tweaked; explicit counts would silently
      under- or over-sample.
    - `shift` (array of three reals): fractional shift along the
      a, b, c reciprocal axes.
- `scf_threshold` (real): SCF convergence threshold for
  the reference run.  Recorded in provenance.
- `basis` (string): the basis set the reference SCF run uses --
  one of `"mb"` (minimal), `"fb"` (full), `"eb"` (extended).
  The initial-potential producer uses the full basis (`"fb"`)
  for reference-quality potentials.  Together with `functional`
  and `kpoint_integration` it selects the guidance predictor's
  sub-model (DESIGN 7.6) and is recorded on the produced
  guidance entry's context.  At present the basis differs from
  `functional` and `kpoint_integration` in *where* it takes
  effect: it is not currently a makeinput setting -- makeinput
  writes all three basis sets into `imago.dat` and the basis is
  chosen at the imago run itself (`scf_basis`, coded `fb -> 2`).
  This may change: a future makeinput could own the basis
  selection directly.  The seam does not depend on the present
  arrangement -- the producer translates each manifest field
  into the tools' own settings and the wingbeat routes each to
  whichever tool recognises it (6.2.10), so the basis can move
  to makeinput later without reworking the seam.
- `functional` (string): the exchange-correlation functional
  token.  `"wigner"` is the Imago default (the Wigner
  interpolation method, makeinput `-xccode` 100).  A predictor
  sub-model selector: a prediction never mixes data converged
  under different functionals.
- `kpoint_integration` (string): the Brillouin-zone integration
  method, as a sub-model-selecting token.
  `"linear-tetrahedral"` (the producer's default; makeinput
  `-scfkpint` 1, the linear analytic tetrahedron) is
  parameter-free, while a Gaussian-smeared method carries its
  smearing width in the token (e.g. `"gaussian-0.1"`, makeinput
  `-scfkpint` 0).  The producer maps this token to makeinput's
  integer integration code and, when the token names a width,
  forwards that width as makeinput's thermal smearing sigma (the
  `-thermsmear` option, written into `THERMAL_SMEARING_SIGMA`, in
  eV).  A bare `"gaussian"` names no width, so makeinput keeps its
  rc-sourced default (no smearing).

  **The tetrahedron method is this path's default, and the reason
  is metals.**  The producer cannot know whether a reference solid
  is a metal before it runs one: that is what the k-point ladder
  discovers, and often not until several rungs up.  So the
  integration scheme has to be chosen while the answer is still
  unknown, which means choosing the one that is safe under both
  outcomes rather than the one that suits the commoner case.

  Unsmeared Gaussian integration assigns each state to the
  occupied set or the empty set by where its eigenvalue falls
  relative to the Fermi level, one mesh point at a time.  In a
  metal the Fermi level cuts through a partly filled band, so a
  small refinement of the mesh can move a state across it and
  change the total energy by a full state's worth.  Refining the
  mesh therefore does not settle the energy smoothly; it rattles
  it, and the rattle does not shrink in proportion to the mesh
  spacing.  That is the mechanism behind the noise floor the
  convergence threshold has to sit above (3.12.3).  The
  tetrahedron method instead interpolates the bands linearly
  between mesh points and integrates over the occupied volume, so
  the occupied fraction of a tetrahedron straddling the Fermi
  level varies continuously as the mesh moves, and the Blochl
  correction removes the leading curvature error of that linear
  interpolation (1.6).

  **The choice costs an insulator nothing.**  Where a real gap
  sits at the Fermi level every tetrahedron is wholly occupied or
  wholly empty, the Blochl weights reduce to one quarter per
  corner, summing over the tetrahedra that share a k-point returns
  the uniform mesh weight, and the tetrahedron method reproduces
  the Gaussian answer exactly.  This is measured, not argued:
  `si_fd-3m_227_2001` at mesh 6-6-6 returns the same total energy
  under both schemes.  So the default is not a trade -- it is free
  where it does not matter and load-bearing where it does, which
  is what lets it be set before the metal question is answered.

  **The default belongs to the ground state, not to every
  calculation.**  A core-level spectroscopy run (XANES/ELNES)
  names `"gaussian"` and gets it: `populateLAT` stops with a
  message naming Gaussian integration as the supported path,
  because the core-hole correction is written against the flat
  sorted occupation array and its band arithmetic does not carry
  over to the `(band, kpoint, spin)` array unexamined.  A wrong
  correction there misplaces exactly one electron, which presents
  as a convergence problem and is never traced back.  That
  refusal is deliberate and is tracked with the rest of the
  occupation work; it does not qualify the ground-state default,
  it scopes it.
- `cell` (string): which cell the reference run computes in --
  `"prim"` (the primitive reduction of the structure's space
  group, the default) or `"full"` (its conventional cell).  It
  becomes the
  `full` / `prim` token in the materialized `imago.skl`, which
  `structure_control` reads to decide whether to reduce; the
  reduction is performed by the space-group operations already
  in hand, not by re-deriving symmetry, and runs one way only
  (`full -> prim`, since the atoms a full cell adds come from
  those operations).

  This setting changes **cost, not physics**.  The harvested
  quantity is a per-species, per-environment potential, and a
  fingerprint is a local descriptor evaluated under periodic
  boundary conditions, so both are cell-invariant: the same
  structure yields the same entries either way.  What changes is
  the size of the problem.  A primitive cell of an n-fold centred
  lattice holds n times fewer atoms and, at a fixed k-point
  density, takes n times more k-points.

  **The primitive cell is the default, on measured evidence.**  A
  scaling argument suggests the saving should be large -- both
  dominant cost terms carry atom count cubed against one factor
  of k-point count, so `1/n^2`, or sixteen-fold for an F-centred
  cubic cell.  The measurement does not bear that out: at these
  system sizes the cubic term does not dominate, and the same
  structure converged in both cells costs

      diamond Si (8 atoms -> 2)   5.1 s -> 2.7 s   1.9x
      Si III BC8 (16 atoms -> 8)  43.6 s -> 22.8 s  1.9x

  at its converged mesh.  So the honest figure is about **twice**,
  not sixteen times, and it should improve for larger cells where
  the cubic term does begin to tell.  Across a whole seeding
  campaign the advantage is smaller again -- about 1.33x -- because
  the primitive cell's larger reciprocal cell makes the mesh climb
  walk further before the energy flattens (eleven rungs against
  seven for diamond).  The campaign figure is the pessimistic one:
  a climb is a one-time cost per material, while a production run
  reads its k-point density from the guidance dataspace and pays
  only the single converged calculation.

  Two independent checks say the choice costs nothing in
  correctness.  Converged energies agree to **0.002 meV/atom** on
  every insulator tested (0.8 meV/atom on the one metal, which
  settles at a deliberately rough mesh anyway).  And the axis
  classes the climb depends on -- which reciprocal axes must share
  a k-point count -- agree between this design's Python port and
  imago's own runtime computation for every cell and centring
  tested, including the case where the reduction genuinely changes
  them: a C-centred orthorhombic cell has three independent axes
  conventionally and only two primitively, because its primitive
  vectors `(a +/- b)/2` have equal length.  That is the conjugation
  of 2.7 being exercised and confirmed, not merely unexercised and
  assumed.

  Two consequences follow from *cost, not physics*.  First,
  `cell` is **not** a predictor sub-model selector, unlike
  `basis` / `functional` / `kpoint_integration`: a prediction may
  freely mix entries converged in either cell, because the
  quantity being predicted -- a k-point *density*, defined per
  unit reciprocal volume (3.7) -- is itself cell-invariant.
  Second, the `converged_mesh` an entry records is **not**
  invariant: it counts points along that cell's own reciprocal
  axes, so a `prim` mesh and a `full` mesh are not comparable as
  bare triples and must be read against the cell that produced
  them.

  `cell` governs only structures the producer *materializes*.  A
  curator-supplied `structure_path` skeleton already carries its
  own `full` / `prim` token, and that token stands: the producer
  writes a curator's file, never rewrites it.

  Alone among the run settings, `cell` is **exempt from the
  resolvability rule** (rule 2) and carries a built-in default of
  `"full"`.  The rule exists so that nothing the producer
  *emits* rides on an implicit default (VISION Principle 11), and
  the other five are each recorded on the entries a run produces
  -- they select the predictor's sub-model, or they land in
  provenance.  `cell` is recorded nowhere: it selects no
  sub-model, and no harvested value depends on it.  A manifest
  that never mentions a cell is therefore not leaving a
  provenance gap, only accepting the conventional cell.

  That exemption has a precise expiry.  The moment `cell` is
  recorded on an entry -- the natural next step, since a database
  otherwise cannot say which cell produced a given potential --
  it becomes emitted knowledge, Principle 11 applies in full, and
  it must join the other five as a required, resolvable setting.
  Because recording it also adds a required provenance field,
  that change is a schema version bump (5.2.5) carrying an
  honestly derivable migration: every file written before the
  default moved to `"prim"` was produced from a conventional
  cell, so the missing field fills with `"full"`.

  Moving the default gives recording it a concrete purpose, though
  a narrower one than it might first appear.  **It does not affect
  dedup, and must never be allowed to.**  The dedup keys on the
  preferred bispectrum descriptor (5.2.3), which the engine
  computes from a periodic neighbour list at seven significant
  figures and which is identical in both cells -- measured, not
  assumed -- so two harvests of one environment in different cells
  collapse to one entry exactly as they should.  Were `cell` ever
  folded into the dedup or the match, it would manufacture
  distinctions the physics does not have: the same environment
  stored twice because a curator drew a different cell, inflating
  the database with redundancy and teaching a learned predictor
  (5.2.4) that cell choice is a feature of an environment.  It is
  not.  `cell` belongs beside `commit` and `generated_at` -- it
  says where a number came from, never what the number means.

  What it does buy is **reconstruction**.  The Imago provenance
  fields exist so the originating SCF can be identified and re-run
  on demand, and the cell is now part of what defines that run: it
  fixes the atom count, and through the reciprocal cell the mesh
  the climb converges on.  Given a database entry and no recorded
  cell, that run cannot be reproduced.  The same fact serves the
  resource-and-cost dataspace (section 8), whose size signature is
  built on `atom_count` and `secular_dimension` -- both halved or
  doubled by this one setting, so a cost model fitted across mixed
  cells without knowing which would be fitting the cell choice as
  noise.

The one *harvest* setting a solid may carry is not a run setting
and resolves against the `[harvest]` block, not `[defaults]`:

- `kpoint_convergence_threshold` (real, units eV/atom): the
  per-atom energy-flatness tolerance the solid's k-point-density
  ladder must reach before its converged rung is harvested
  (DESIGN 7.8).  Given per solid to override just that solid's
  value, else inherited from the top-level `[harvest]` block,
  else the producer's built-in default of 1e-3 (1 meV/atom).
  It never enters a calculation -- it is consulted only when the
  finished runs are read back -- and the resolved value is
  recorded on each guidance entry the solid contributes.

The adaptive mesh climb's tuning knobs (3.12.6) live in an
optional `[harvest.kpoint_climb]` sub-table.  Unlike
`kpoint_convergence_threshold`, they are **database-wide** -- they
tune how the whole seed campaign searches, not one solid's physics
-- so they resolve `[harvest.kpoint_climb]` -> built-in default,
with no per-solid override.  Each knob is optional; an omitted one
takes the producer's provisional default (3.12.6), so a manifest
that names neither the sub-table nor any knob searches with the
built-in policy.  The knobs are:

- `max_count` (positive int): the fixed per-axis count ceiling the
  climb never exceeds, the backstop that guarantees termination
  when the energy never quite goes flat (3.12.3).
- `confidence_high` (real, [0,1]): the predictor-confidence at or
  above which a prediction is trusted enough to dispatch a small
  fixed mesh grid in one round rather than climb (3.12.5).
- `grid_width` (int >= 0): how many rungs to each side of the
  predicted rung the confident mode's fixed grid spans (3.12.5).
- `start_offset_moderate`, `start_offset_cold` (int >= 0): how
  many rungs below the predicted seed the climb mode begins -- the
  moderate offset for a low-confidence prediction, the larger cold
  offset for an under-trained one, since a weaker guide starts
  lower (3.12.4).
- `flat_needed_confident`, `flat_needed_cold` (positive int): how
  many consecutive flat interior rungs the stop test must see
  before it accepts convergence -- one when confident, two when
  cold, so a single lucky flat step cannot end a cold climb early
  (3.12.3).
- `max_stride` (positive int): the cap on the bracket phase's
  geometric stride growth, which keeps the interval a refine must
  then fill short (3.12.3).
- `climb_shape` (one of `bracket-refine`, `unit-step`): which
  search the low and under-trained cases run -- the bracket-then-
  refine climb by default, or the fine unit-step walk when a
  curator pins the most conservative reading (3.12.3 / 3.12.5).
- `stride_flatness_multiple` (real >= 1): the multiple by which
  the bracket phase's flatness threshold is looser than the
  convergence threshold, setting how eagerly a nearly-settled
  stride is bracketed and so how much top-end overshoot is shaved
  (3.12.3).
- `metal_gap_threshold` (real, eV): the band-gap value below
  which a rung is read as gapless, so the climb declares the
  material a metal and settles at once on a rough, floor-level
  mesh; low enough that no real insulator crosses it, high enough
  to catch a true metal's near-zero gap (3.12.3).  Applies to
  every climb shape.  A NEGATIVE value turns the test off, since
  no band gap can be negative -- the way a curator asks for every
  rung of a known metal to be computed, for a diagnostic ladder
  (3.12.3).  This needs no special case in the rule itself, which
  is why it is stated here as a use rather than coded as one.
- `crystalline_floor_axis_count` (positive int): the per-axis cap
  on a crystalline climb's opening floor mesh -- the most points
  any one axis of the coarsest starting mesh may carry, `4` by
  default, so a cubic cell floors at `[4,4,4]` and an anisotropic
  one lower per axis (3.12.4).  Holding the climb above this floor
  is what lets the metal test read a trustworthy gap at every
  rung, with no separate coarse-mesh guard (3.12.3).

**Per-entry fields (`[[reference_solid.entry]]`).**  An entry
is an *optional customization* on an auto-discovered environment, not
a harvest instruction (5.2.2, 5.2.3).  Every field is optional;
an entry exists only to override what the harvest would
otherwise produce, and an environment no entry mentions is
harvested all the same.

- `atom_site` (positive integer, *optional*): the handle that
  addresses the target environment -- a representative atom of
  it, 1-based into the structure (Imago's site-indexing
  convention).  The species number cannot be the handle because
  it is unknown until the run's grouping pass; a representative
  atom is known up front, or chosen after an exploratory run.
  When given, it also *pins* that exact atom as the
  environment's representative instead of the computed
  order-independent one.
- `element` (string, *optional*): element symbol the customization
  targets.  Cross-checked against the species at `atom_site`
  after Imago loads the structure.
- `label` (string, *optional*): the label this environment is
  written under in the element's `s_gaussian_pot.toml`.  When
  omitted, the producer derives it at harvest per 5.2.1
  (`<reference_id>-<element><species>-t<type>-a<site>`); when
  present it overrides that derived default.  An explicit
  `(element, label)` must be unique across the entire manifest;
  a derived label's `(reference_id, element, atom_site)` must be
  unique (rule 6).
- `default` (bool, *optional*): marks this environment the
  element's default.  At most one environment per element
  carries `default = true` (rule 7); when no customization does, the
  producer falls back to the element's `"isolated"` baseline,
  so each per-element file still gets its required-single
  default entry per 5.2 rule 7.  An environment with no customization
  is `default = false`.
- `description` (string, *optional*): one-sentence prose
  explanation of the chemical environment, copied verbatim into
  the database entry's `description` field (5.2).  When omitted,
  the producer auto-composes one from the reference solid's
  `source_description`, qualified by the species and site
  (5.2.1).

**The `[characterization]` block (database-wide recipe).**  The
fingerprint recipe -- one preferred `sub_spec` per method,
applied to every harvested environment -- is a constant of the
whole database, so it is declared once in a top-level
`[characterization]` block rather than repeated on every entry.
Each `[[characterization.fingerprint]]` carries a `method` and
its `sub_spec`.  Lifting the recipe here makes rule 11
structural: a single declaration per method cannot diverge
between elements.  Every harvested environment computes these
fingerprints, and the producer stamps the preferred flag onto
each matching record in the per-element files.

**The `[defaults]` block (shared run settings).**  The five run
settings -- `kpoint_spec`, `scf_threshold`, `basis`,
`functional`, `kpoint_integration` -- are typically identical
across a manifest's solids, so they may be declared once in a
top-level `[defaults]` block and inherited by every
`[[reference_solid]]`.  A solid may still name any of them
itself to override its own value; an omitted setting resolves to
the `[defaults]` value.  The block is *optional*: a manifest that
spells out all five on every solid needs none.  Each setting must
be *resolvable* for every solid -- present on the solid or in
`[defaults]` -- or the load is refused (rule 2).  Nothing the
producer emits depends on an implicit default (VISION Principles
5 and 11): the value is always written down, just possibly once
for the whole manifest rather than once per solid.

**The `[harvest]` block (shared harvest settings).**  A separate
optional block for the settings that govern how finished runs
are *read back* into the database, as opposed to how the
calculations are *run*.  Keeping the two apart is deliberate:
the `[defaults]` settings each feed an individual calculation
(they become makeinput/imago options baked into the run),
whereas a harvest setting is consulted only afterward, when the
producer interprets results that already exist on disk.  In v1
the block holds a single setting, `kpoint_convergence_threshold`
(real, units eV/atom): the per-atom energy flatness a solid's
k-point-density ladder must reach before its converged rung is
harvested (DESIGN 7.8).  Like the run settings it may be named
per solid to override just that solid's value, resolving to the
`[harvest]` block when the solid is silent; unlike them it also
carries a built-in producer default (1e-3 eV/atom = 1 meV per
atom), so a manifest may omit the block entirely.  The
*resolved* value -- whether authored or defaulted -- is recorded
on every guidance entry the run contributes (its
`metric_threshold`, 7.8), so the threshold actually used stays
recoverable even when it was never written in the manifest.  The
`[characterization]` block above is a harvest-configuration
block in the same spirit (it declares which fingerprints to
compute at harvest); it stays a sibling of `[harvest]` rather
than folding in, because it is database-wide and already fully
specified on its own terms.

**Per-entry fingerprint declarations
(`[[reference_solid.entry.fingerprint]]`).**  A per-entry
declaration is a *rare override*: it asks the producer to
harvest one *extra*, non-preferred fingerprint for this
environment alongside the database-wide preferred ones above.
The common entry carries none.  The record is written into the
database entry's `[[potential.fingerprint]]` array per 5.2.

- `method` (string): matcher name.  Must be a method known
  to the matcher registry in `matchers.py` (8.9).
- `sub_spec` (inline table): method-specific parameters.
  Two declarations on the same entry with the same `method`
  and the same `sub_spec` keys-and-values are a hard error
  (rule 8); same `method` with differing `sub_spec` is
  explicitly allowed -- the database stores as many sub_specs
  per family as the curator wants.
- A per-entry declaration may not be marked `preferred`: the
  preferred record of each family is fixed database-wide by the
  `[characterization]` block (rule 11), so per-entry overrides
  are always non-preferred alternates, stored for validation or
  comparison.

**Canonical-sub_spec convention (the `preferred` record).**  A
file-dictated match (5.6.5 step 2) needs one unambiguous
`sub_spec` per family to query at, so that a multi-element
structure costs a single loen run rather than one per distinct
sub_spec.  The `preferred` flag supplies it: across the whole
database, exactly one bispectrum `sub_spec` and one reduce
`sub_spec` are marked preferred, and every element's preferred
record sits at that same sub_spec.  The consumer reads the
preferred record straight off each element's database and never
searches across sub_specs.  Non-preferred records (alternate
sub_specs harvested for validation or comparison) ride along in
the same files without affecting the consumer -- they are
simply never the one it picks.  This keeps the runtime analysis
a simple best-effort lookup that never steps outside the
convention, while leaving the database free to accumulate
alternate descriptors over time.

Producing a fingerprint requires the matcher's compute step.
For Python-side matchers (`reduce`), this runs in-process
from the reference structure.  For Fortran-side matchers
(`bispectrum`), the producer runs `imago.py -loen -scf no`
on the reference structure with `method`/`sub_spec` mapped
into the loen-side input parameters; the output (`fort.21`)
is parsed and the row for `atom_site` becomes the
fingerprint record's payload.

**Validation rules.**  The manifest loader refuses to proceed if
any rule below fails — strict refusal, no last-wins fallback, no
warning-and-continue.  Behavior mirrors the per-element database
file (5.2):

1. `schema_version == 2`.
2. Every `[[reference_solid]]` carries `reference_id`,
   `system_type`, and *exactly one* of
   `{(cod_id, cod_revision), structure_path}` (see rule 4 for
   details).  `system_type` must be one of the four allowed
   values `{"crystalline", "amorphous", "nanostructure",
   "molecular"}`; any other value is a hard error, since the
   guidance predictor (DESIGN 7) switches its sub-model on it and
   the produced entry records it for forensics.  The five run
   settings -- `basis`, `functional`, `kpoint_integration`,
   `kpoint_spec`, `scf_threshold` -- must each be *resolvable*
   for the solid: present on the solid itself or supplied by the
   top-level `[defaults]` block.  They select the predictor
   sub-model (DESIGN 7.6) and are recorded on every produced
   entry's context, so nothing the producer emits depends on an
   implicit default (VISION Principles 5 and 11) -- the value is
   always written down, once per solid or once in `[defaults]`.
   The harvest setting `kpoint_convergence_threshold` is exempt
   from this resolvability requirement: it carries a built-in
   producer default (1e-3 eV/atom), so a solid that names neither
   it nor a `[harvest]` block is not refused -- the default
   applies and its resolved value is still recorded in provenance
   (the guidance entry's `metric_threshold`, 7.8).
   A top-level `[characterization]` block
   declaring at least one fingerprint is required (rule 10): it
   sets the database-wide preferred recipe, so a manifest without
   one is refused rather than silently producing a database with
   no preferred descriptors for the consumer to match against.
   The relaxed `--materialize-only` reader does not apply this --
   it only materializes structures and never harvests.
3. `[[reference_solid.entry]]` blocks are *optional* customizations
   (5.2.2): a reference solid may carry none.  When an entry is
   present, all its fields are optional; `atom_site`, when given,
   is the representative-atom handle that addresses the target
   environment, and the species and type numbers (unknown until
   the grouping pass runs) are never authored.  `label`, when
   present, is an explicit curator override (5.2.1); when absent
   the producer derives it at harvest from the run's site
   identity.
4. Exactly one of `cod_id` or `structure_path` is set on each
   `[[reference_solid]]`.  If `structure_path`, it resolves to
   an existing file under the manifest's directory.  If `cod_id`,
   it parses as a positive integer *and* `cod_revision` is
   present as a non-empty string.
5. `reference_id` is unique across the manifest *and* label-safe
   (lowercase letters, digits, `-`, `_`): it is embedded verbatim
   in every derived entry label and typed into `-pot`, so a
   non-conforming id is a hard error (5.2.1).
6. No two database entries may collide on a label.  A
   customization with an explicit `label` must be unique on `(element, label)`
   across the manifest — two solids cannot both produce, e.g.,
   `("Au", "default_solid")`.  Derived labels are unique by
   construction: each auto-discovered environment mints
   `<reference_id>-<element><species>-t<type>-a<site>` (5.2.1)
   from its own run identity, so two environments differ in
   species, type, or site within a solid and in the
   `reference_id` prefix across solids.  A silent overwrite that
   would mask a curation mistake is refused.
7. At most one harvested environment per element carries
   `default = true`, set by a `default = true` customization.  At
   load time no element may carry two such customizations.  At
   harvest, an element whose harvested environments include
   no default customization is not an error: the producer tags that
   element's `"isolated"` baseline as the default
   (`is_isolated_default_for`, PSEUDOCODE 11.4), so a
   hand-free element still yields exactly one default entry
   (5.2 rule 7).  When a customization does designate a
   default the manifest is the single source of truth; the
   isolated
   baseline is only the no-designation fallback.
8. Within any one `[[reference_solid.entry]]`'s fingerprint
   declarations, `(method, sub_spec)` is unique.  Same
   `(method, sub_spec)` declared twice on the same entry is
   a hard error.  Same `method` with differing `sub_spec`
   coexists on the same entry by design.  A per-entry
   declaration may not carry `preferred = true` (rule 10).
9. Every `method`, in `[characterization]` or in a
   `[[reference_solid.entry.fingerprint]]` declaration, must be
   a name known to the matcher registry in `makeinput.py` (8.9).
   Unknown methods are a hard error rather than a silent skip.
10. The `[characterization]` block declares *at most one*
    fingerprint per `method`, and that single declaration is the
    family's database-wide preferred record (its `sub_spec` is
    the one the consumer queries, 5.6.5 step 2).  A `method`
    named twice in `[characterization]`, or a per-entry
    declaration marked `preferred = true`, is a hard error: the
    preferred recipe has exactly one home.
11. Rule 11 -- the preferred `sub_spec` for a family is uniform
    across the whole database -- holds *structurally*: because
    the preferred recipe is the single `[characterization]`
    declaration per method, it cannot diverge between elements.
    Non-preferred per-entry overrides may use any `sub_spec`;
    only the global preferred one is constrained, and it is
    constrained by construction.

**Structure materialization.**

Before a reference solid can be handed to kaleidoscope its
structure must exist as a local file, because kaleidoscope keys
its run-reuse cache on the structure file's contents (DESIGN
6.2.5).  `materialize_structure(ref)` produces that file and
returns its path:

- For a `structure_path` entry it reads the on-disk file named by
  the manifest (resolved under the manifest's directory).  No
  network.  The file is an `imago.skl` (the run consumes skl);
  for a crystal the curator authors or converts it with the
  space group preserved (see `cif2skl`, ARCHITECTURE 9.5).
- For a `cod_id` entry it fetches the structure once from the
  Crystallography Open Database at the pinned `cod_revision`
  (`cod_fish.py get`, the canonical fetch this step imports) and
  converts the fetched CIF to an `imago.skl` with `cif2skl`,
  which preserves the CIF's space group -- recovering the
  asymmetric unit and the `spaceDB` setting rather than
  flattening to P1, because the Brillouin-zone integration
  samples the irreducible wedge using that space group
  (ARCHITECTURE 9.5).  A CIF whose space group cannot be
  resolved to a `spaceDB` setting is a hard error (no silent P1
  fallback for a crystal); the curator then supplies a
  pre-converted `structure_path` skl instead.

Both artifacts -- the fetched CIF and the converted skeleton --
are cached on disk under `share/curation/structures/`, beside
the producer's other reconstructible working artifacts
(ARCHITECTURE 8.1), and a converted skeleton is reused when it
is already present rather than rewritten.  **The cached
skeleton's
name must therefore carry every manifest setting that changes
what the conversion writes**, or a later run under different
settings will silently be handed the earlier run's file.  Today
that is exactly one setting, `cell`, so the skeleton is named
`<reference_id>-<cell>.skl`; the fetched CIF needs no such
qualifier, since `cod_id` and `cod_revision` already pin its
bytes and no manifest setting alters them.

The rule matters more than the present instance.  A cache keyed
too coarsely does not fail loudly -- it returns a stale artifact
that is perfectly well-formed, so the run succeeds and reports
the wrong thing.  A comparison of two settings would then find
them identical *because* the second never happened, which reads
as a confirmation rather than an error.  Any future setting that
reaches the conversion joins the name for that reason.

This step is the producer's only network access and is
**deliberately decoupled from any run cache**: its sole job is to
guarantee a local structure file and hand back its path.  It
carries no SCF results, no convergence state, and no
hit/miss comparison logic — those now belong to kaleidoscope's
run-reuse cache.  Pinning `cod_revision` keeps the build
deterministic against upstream COD edits: re-running months later
fetches the same bytes by construction.

*COD-fetch is strict.*  Fetch failures (network down, COD outage,
the pinned revision missing) error out, name the failing fetch,
and refuse to fall back to any other revision.  A silent fallback
would produce a structure inconsistent with the pinned manifest —
exactly the failure mode pinning exists to prevent.

**Run-reuse caching is kaleidoscope's job, not the producer's.**

The producer keeps no SCF cache of its own.  Because it
delegates every SCF run to kaleidoscope, the avoid-recompute
responsibility sits with kaleidoscope's run-reuse cache (DESIGN
6.2.5), which keys each run on the structure file's contents
together with the makeinput options and the Imago build
identity.  A single mechanism supplies every property a
producer-side cache would have had to provide:

- *Edits to harvest declarations stay cheap.*  Adding or changing
  a `[[reference_solid.entry]]` changes neither the structure nor
  the flight options, so the run is a cache hit and only the
  harvest step (`extract_potential`, PSEUDOCODE 11.4) re-runs.
- *Content, not path.*  The cache compares structure file
  *contents*, so renaming a `structure_path` file on disk does
  not force a re-run.
- *Threshold-sensitive, build-insensitive.*  Changing the SCF
  convergence limit or the k-point grid invalidates the cached run,
  as it must.  Rebuilding the engine does not: the build identity is
  recorded beside each run and never compared (6.2.5, VISION 16), so
  a re-run after a rebuild is something a curator asks for with
  `--force` rather than something the cache decides for them.

The reason the old design preferred direct comparison over a
content hash still holds — at ~100 reference solids run by hand,
naming the exact field that changed beats a bare "different hash"
— and DESIGN 6.2.5 carries that reasoning for the kaleidoscope
cache.

**Procedure.**

The pipeline runs in three phases — *build*, *converge*,
*harvest*.  The build phase predicts and seeds every solid; the
converge phase walks each solid's k-point mesh upward in adaptive
parallel rounds until its energy flattens; the harvest phase
reads each converged run.  Convergence is no longer a static grid
guessed up front and dispatched in one batch: it is the adaptive
mesh climb (DESIGN 3.12), which searches in *mesh* space -- so
every rung is a genuinely distinct calculation, never several
densities collapsing onto one symmetric mesh -- while it *keys*
and *records* its result in the transferable k-point *density*.
Predict in density, search in mesh, record in density.

1. Load and validate the manifest (the rules above).  Resolve the
   climb policy once for the run: merge the manifest's optional
   `[harvest.kpoint_climb]` knobs over the provisional defaults
   into the confidence thresholds and the per-axis `max_count`
   ceiling the climb reads (3.12.6).
2. For every element with a directory in `share/atomicPDB/`,
   load the element's existing `s_gaussian_pot.toml` if it is
   present (preserving every environment harvested by earlier
   runs), or start a new in-memory database if it is not; then
   refresh (or create) its `"isolated"` entry directly from the
   current `pot1` and `coeff1` files.  Loading rather than
   resetting is what makes the build incremental (5.2.3);
   refreshing the baseline guarantees it is always present
   (rule 6 of 5.2) and tracks any changes in atomSCF output.
3. **Build.**  For each `[[reference_solid]]` in the manifest:
   a. `materialize_structure(ref)` → a local structure file
      (read from disk for `structure_path`, or fetched once from
      COD at the pinned `cod_revision` for `cod_id`).
   b. Predict and seed.  `options = make_producer_options(ref)`
      translates the manifest physics into each tool's coded
      settings (`functional` → `xccode`, `kpoint_integration` →
      `scfkpint`, `basis` → `scf_basis`, `scf_threshold` →
      `converg`, `kpoint_spec.shift` → `kpshift`; 6.2.10), while
      `submodel = {basis, functional, kpoint_integration}` keeps
      the human names the predictor and record speak.  Then
      `predict_kpoint_density(structure, options, dataspace,
      system_type, submodel, center)` (DESIGN 7.7) consults the
      guidance dataspace and returns the **seed density**, the
      prediction's `confidence` and under-trained flag, and the
      `PredictionRecord` the harvest later recovers.  It lays *no*
      grid: the climb (step 4) searches outward from the seed.  A
      manifest `kpoint_spec.density` is passed as `center` -- a
      curator override that pins the seed and bypasses the
      predictor.  From that confidence and the run's resolved
      climb policy, build the solid's `ClimbConfig`: the dispatch
      mode (`PARALLEL_GRID` when confident, `CLIMB` when cold or
      under-trained), the flatness persistence `flat_needed`, and
      the reciprocal-cell geometry the rung mechanics read
      (3.12.4 / 3.12.5): the reciprocal magnitudes and cell
      volume from the loaded cell, and the axis classes
      recomputed in Python from the cell's conventional-abc
      space-group operations -- read through the same shared
      `symmetry.read_conv_abc_point_ops` the kp writer uses, so
      the classes the climb seeds from are derived from the very
      operations imago will run under (2.7).
   c. Emit one structure-only `imago.py -loen -scf no` unit per
      declared Fortran-side fingerprint (the bispectrum
      fingerprint depends on geometry alone, not on SCF
      convergence, so it need not wait for a converged mesh).
      These loen units are geometry-only and mesh-independent, so
      they belong to no climb round: they are dispatched once,
      together, in a small pre-flight batch before the climb, and
      their run dirs persist for the harvest's fingerprint step.
      The *convergence* units are not built here -- the climb
      builds each round's meshes as it runs (step 4).
4. **Converge.**  Drive every solid through the climb
   (`converge_by_climb`, DESIGN 3.12.5): serial within a solid,
   concurrent across, and no solid waiting on another.  Every
   solid's seed mesh(es) launch at once (`initial_meshes` from the
   seed density and mode -- a small grid when confident, one
   starting rung when cold).  Thereafter the producer collects
   rungs as they land, one at a time: each landing advances the
   single solid it belongs to -- reading that solid's rungs so far
   and either accepting a converged rung, stopping at the
   `max_count` ceiling, or launching one more mesh -- and leaves
   every other solid untouched.  Each mesh rides an explicit
   axial-count k-point file (`build_mesh_unit`, DESIGN 7.7)
   launched through kaleidoscope's wingbeat seam and run-reuse
   cache (DESIGN 6.2.5), so a mesh re-run later is a cache hit.  A
   solid leaves the climb the moment it converges or ceilings, and
   because a solid climbs on the instant its own rung lands, a
   late, expensive solid never holds back the ones already done.
   The producer runs no SCF itself.  Each solid's
   outcome is either its **settled rung** (the mesh, its
   energy, and the ascending distinct-mesh ladder below it) or
   `NON_CONVERGED` -- a ceiling, or a rung that failed to run.  A
   non-converged solid is flagged, never retried here: retries are
   the runner's job (Principle 12).

   **A rung whose SCF did not converge never joins the ladder.**
   Finishing and converging are different things, answered in
   different places: the flight entry's status says the job
   completed, while `converged` in the run's own result.toml says
   the SCF reached its fixed point.  A `NOT_CONVERGED` run exits
   cleanly and writes a total energy (6.2's run statuses already
   say its outputs "must not be harvested as a reference
   potential"), so nothing marks that energy as unfinished once it
   is a number on a ladder.  It is wherever the iteration
   happened to stop.  Handing it to a flatness test inverts what
   the test is asking: the test wants to know whether the energy
   has stopped moving with the MESH, and this one stopped moving
   for a reason that has nothing to do with the mesh.  It can read
   flat by coincidence, and it can break a plateau that was real.

   Such a rung is treated as a rung that did not run, which stops
   that material.  The alternative -- drop the rung and keep
   climbing -- was rejected: the climb chooses its next mesh FROM
   the ladder, so a ladder that does not grow yields the same
   request again, forever.  Stopping surfaces the problem to the
   curator, which is what a run that would not converge deserves.

   **The reason for the stop travels with the rung.**  Two
   different reasons yield a settled rung: the energy went flat
   (a k-converged insulator) or the material read gapless (a
   metal settled on a deliberately rough potential, 3.12.3).
   They are not interchangeable downstream, so the climb reports
   which one it was rather than leaving each later stage to work
   it out again from whatever evidence it happens to hold. Three
   verdicts, then: `converged`, `metal`, `not_converged`.

   This matters because a *classification* the system already
   made should not be re-derived from weaker evidence later. The
   climb reads the whole ladder; a later stage looking at one
   rung sees an artificial gap whose size depends on where the
   mesh points fell (1.6), and can reach the opposite conclusion
   with nothing to warn it. Carrying the verdict makes the metal
   path a decision taken once and remembered, which is also what
   lets the run log say what actually happened.
5. **Harvest.**  For each `[[reference_solid]]`:
   a. If the solid is `NON_CONVERGED`, flag it in the run log and
      harvest no potential (a non-converged sweep earns neither a
      potential nor a guidance entry).  Otherwise its converged
      rung names the mesh whose `kpt-mesh-<a>-<b>-<c>` run dir --
      the run the climb already dispatched -- carries the
      converged potential, and the steps below read that run.
   b. Record that run's SCF iteration count and convergence
      metrics in the run log.
   c. Discover the run's environments: every atom carries a
      species after the grouping pass (5.6.4), and the assigning
      method's partition defines the distinct environments
      (5.2.2).  For each distinct environment take one
      order-independent representative (5.6.5), unless a curator
      customization pins an `atom_site` for it, in which case that atom
      is the representative:
      i.   `extract_potential` for the representative's
           `atom_site` from the converged run.  The converged
           `scfV` output lists every potential type (a
           `NUM_TYPES` header + per-type Gaussian blocks under a
           `TOTAL__OR__SPIN_UP` channel); the harvest selects the
           site's type block -- its type number read from
           `datSkl.map` (ARCH 9.7) -- and takes each term's
           coefficient and alpha (columns 1 and 2) together.
      ii.  Compute the `[characterization]` fingerprints (and any
           per-entry override) at the representative.  Python-side
           matchers (e.g., `reduce`) compute in-process from the
           run's *expanded* structure -- `outputs["structure"]`,
           makeinput's `imago.fract-mi` full cell -- not the
           materialized source file, which for a space-grouped
           reference is only the asymmetric unit and carries
           neither the full-cell geometry the shells need nor the
           run's numbering.  Because that expanded skeleton is
           ordered by the run's sorted (dat) numbering while
           `atom_site` is a skeleton index, the harvest maps
           `atom_site` to the structure row through the same
           `datSkl.map` used in step i.  Fortran-side matchers
           (e.g., `bispectrum`) read the matching `-loen` unit
           that kaleidoscope already dispatched in step 3c and
           parse the row for `atom_site` from `fort.21`.  Build a
           `FingerprintRecord` (5.4) per method; the record at a
           `[characterization]` method's `sub_spec` is marked
           `preferred = true`.
      iii. Build the entry.  Its `label`, `default`, and
           `description` come from a curator customization when one
           addresses this environment; otherwise `label` is
           derived and `description` auto-composed (5.2.1), and
           `default = false`.  Carry the run-supplied numerical
           fields (the harvested coefficients, 5.2.3), the
           provenance (recording the solid's `system_type` for
           forensics, per rule 2), and the `FingerprintRecord`
           list from step ii.
      iv.  Insert-or-skip the entry into the in-memory
           `ElementDatabase` for its element (5.2.3).  If the
           environment duplicates one already present under the
           dedup rule, skip it -- the stored representative
           stands; otherwise append it.  An explicit
           customization `label` naming an existing entry
           replaces it.
   d. **Guidance contribution.**  The climb already holds the
      solid's converged rung and its ladder in memory, so the
      guidance entry is built *in place* -- no re-read of the
      workspace.  `record_converged(rung, rungs, config)` (DESIGN
      3.12.4) turns the converged rung into the entry's k-point
      density and the ascending flatness ladder; those chosen
      facts, together with the converged run's `result.toml` (gap,
      magnetization, SCF threshold, exact mesh, commit), feed the
      one entry builder `build_entry` -- the *same* builder the
      standalone density harvest uses, so the two paths cannot
      diverge -- and `save_entry` stages the result under
      `share/historicalGuidanceDB/staging/<system_type>/`.
      `build_entry` reads the exact converged mesh from
      `result.toml`, the same source in both harvest paths, so the
      mesh stored beside the density cannot differ between them.
      Every reference solid the producer converges thus becomes
      training data that sharpens the predictor for the next
      solid.  A climb that stops on its *two-sided* test carries
      at least the three distinct rungs that test required
      (3.12.3), so the stored flatness ladder is long enough for
      the curator's `auto_promote_ok` to re-judge, and a
      curator-pinned seed is verified by the climb like any other
      (there is no unverified single-point "trust" harvest here,
      unlike the density-era helper, DESIGN 7).  The metal
      short-circuit is the exception and the reason this is not
      stated as "every converged solid contributes an entry": it
      stops at the *first* gapless rung, so its ladder can be a
      single point and its settled density is not a convergence
      claim at all.  Metals therefore stage no guidance entry
      (7.8); their potential is still harvested, which is what
      this database wanted from them.
6. Save each affected `ElementDatabase` to disk via
   `initial_potential_db.save()` (5.5).
7. Write `share/curation/run_log.toml` capturing the manifest
   snapshot, per-run iteration counts, the settled mesh and its
   k-density chosen for each solid, and the Imago commit.  The
   validation harness (5.8) reads this log.

   Each row carries the climb's `verdict` verbatim -- one of
   `converged`, `metal`, `not_converged` -- and a `converged`
   boolean that means what it says: k-point converged, and so
   FALSE for a metal.  The two fields are not redundant, because
   they answer different questions.  `converged` answers "is
   this energy a converged one?", which the validation harness
   needs before it treats a row as a reproducible target.
   `verdict` answers "why did it stop?", which is what
   distinguishes the two false cases from each other: a metal
   row names a mesh and DID yield a potential, a not_converged
   row yielded nothing at all.  Collapsing them would force a
   reader to guess which kind of false it was looking at.

**Reporting: outcomes and problems, not progress.**  A run says
what it achieved and what went wrong with it.  It does not
narrate what it is doing.  This is VISION Principle 10 read at
the level of a single run: failures are to be *surfaced*, and a
failure printed among two hundred lines reporting that nothing
went wrong has not been surfaced in any sense that matters.

A pre-flight in which every reference structure arrives intact
prints nothing at all, because there is nothing the curator has
to act on; a structure that fails to materialize is named, with
its reason, whether or not anyone asked for detail.  The same
split governs the in-flight tidying of 6.2.12: a scratch tree
that will not remove is a broken assumption and is always
reported, while a tree pruned or refused according to policy is
the mechanism working and is not.

The narration is worth keeping for the occasions when a fetch or
a prune misbehaves and the question is *which* one, so
`--verbose` restores it rather than deleting it.  This is the
split `tidy_scratch.py` already draws between its summary and
its per-item listing, and the reason is the same: output printed
on every successful run stops being read, and real failures then
hide inside it.

Verbosity is deliberately **not** threaded through the build's
signatures.  It is a property of how the process talks to its
user, not of how a flight converges, and passing it down through
`converge_by_climb` into the prune hook would put a reporting
concern into four functions that are otherwise about physics and
dispatch.  It is instead one module-level setting, established
once from the parsed flags before any work begins and read by
the reporting helpers alone.

**Flags:**

- `--verbose`: print the per-item narration that is otherwise
  suppressed -- each reference solid as it is fetched and
  converted, and each scratch tree as it is pruned or refused.
  Failures and warnings appear with or without it, as does the
  closing summary of a completed build.  The pre-flight's own
  tally is not a closing summary and follows the rule above: it
  is printed when something failed, because it then says how much
  of the set survived, and when this flag asks for it.
- `--force`: forward a cache-bypass to kaleidoscope so every
  dispatched unit re-runs from scratch instead of reusing the
  run-reuse cache (DESIGN 6.2.5).  Fresh results are still written
  into that cache afterwards, so the next ordinary run is a
  warm-cache hit.
- `--manifest PATH`: alternate manifest location (default:
  `share/atomicBDB/manifest.toml`).
- `--element ELEM`: restrict the build to a single element's
  `s_gaussian_pot.toml`.  Reference solids whose entries
  contribute only to other elements are still dispatched (so the
  run-reuse cache is warmed and a follow-up run without
  `--element` benefits) but are skipped at the harvest step.

### 5.8 Validation Harness Algorithm

`bench_initial_potential.py` implements the headline
metric from VISION Principle 7: a >= 20% reduction in
average SCF iteration count when starting from improved
initial potentials versus the isolated-atom baseline.

**Inputs:**
- A benchmark manifest (TOML) listing test systems and
  their SCF settings.
- The augmented database from 5.7.

**Procedure:**

1. Load and validate the benchmark manifest.
2. For each test system:
   a. Run Imago with `-pot isolated`; record
      `iter_isolated`.
   b. Run Imago with `-pot default_solid`; record
      `iter_default_solid`.
   c. Compute per-system reduction:

      pct = (iter_isolated - iter_default_solid)
            / iter_isolated * 100

3. Aggregate:
   - `mean_pct` over all test systems.
   - `held_out_mean_pct` over only the test systems
     whose `reference_id` does not appear in the
     curation manifest.  The benchmark manifest must
     include at least one held-out system, or the
     harness aborts with a configuration error.
4. Emit a comparison report
   (`share/curation/bench_report.md`) listing per-
   system counts, per-system reductions, the overall
   mean, and the held-out mean.  The report ends with
   a pass/fail line: PASS if `mean_pct >= 20`, FAIL
   otherwise.  The held-out mean is reported as a
   sanity indicator but does not gate PASS/FAIL on
   its own.
5. Exit 0 on PASS, 1 on FAIL, so CI can gate on the
   harness if desired.

### 5.9 Open Design Questions

Two of the three Phase-2 questions from the original draft
of this section are resolved by the Phase-2 selection
algorithm (5.6), the manifest schema bump (5.2), and the
matcher protocol (ARCHITECTURE 8.9).  The remaining
question is parked for Phase 3; one new follow-up is
filed under Phase 2.

**Resolved by Phase 2:**

- **Per-site label selection** -- resolved.  The
  fingerprint-driven matcher dispatch in 5.6.5 picks a
  manifest entry per species (and therefore per atom)
  from the species-centroid fingerprint compared
  against curated reference fingerprints.  No global
  one-label-fits-all assumption remains.
- **Descriptor computation** -- resolved.  Each
  matcher (`reduce`, `bispectrum`, future
  `bispectrum-by-element`, future `soap`) defines its
  own descriptor in its `sub_spec` block plus the
  payload it stores in the manifest's fingerprint
  records.  Adding a new descriptor is a new matcher
  in the registry (8.9) plus a new
  `[[potential.fingerprint]]` shape -- no schema
  rewrite.

**Phase 2 follow-up:**

- **Element-aware bispectrum.**  The current
  `computeBispectrumComponent` in `loen.f90` does not
  account for neighbor element identity: a C atom
  surrounded by six O neighbors produces the same
  fingerprint as one surrounded by six C neighbors.
  Making it element-aware is real new Fortran work
  -- a new input parameter in `O_Input` (e.g.,
  `bispecByElement`), a per-neighbor-element
  accumulation in `computeBispectrumComponent`, an
  extended `fort.21` output format that labels the
  per-element vector slices, and a matcher-distance
  update that zips by neighbor-element symbol on the
  Python side.  Until then, the element-agnostic
  bispectrum is the only `bispectrum` matcher
  variant.  Scheduled in TODO under the Phase-2
  follow-up chain.

**Parked for Phase 3:**

- **Interpolation when no fingerprint match clears
  the floor.**  Phase 2 falls back to the
  default-tagged entry with a warning when the best
  fingerprint match exceeds the matcher's
  similarity floor (5.6.5 step 3).  A Phase-3 design
  could replace that fallback with a numerical
  blend: e.g., distance-weighted average across the
  K-nearest manifest entries in descriptor space,
  with weights derived from the matcher's distance
  metric, or a learned predictor over a corpus of
  reference fingerprints.  Postponed because the
  default-tag fallback is acceptable for the
  initial Phase-2 deliverable, and a sensible
  blending scheme requires accumulated experience
  with how often the floor is exceeded in practice.

### 5.10 Sequential loen for Fortran-side matchers

Matchers split into two families by where the descriptor
computation lives:

- **Python-side** (`reduce`).  Computes from
  `StructureControl` in-process during the species pass
  of 5.6.4.  No external program; no intermediate files.
- **Fortran-side** (`bispectrum`).  Computes inside the
  Imago engine's loen path, which needs a populated
  `imago.dat` to read the structure, the k-points, and
  the `(twoj1, twoj2)` parameters.  The run writes
  `fort.21`: one row per potential site (the columns are
  given in 5.10.3).

Because a Fortran-side descriptor can only come from a
*completed* engine run, it is obtained by a short,
explicit **sequence** of ordinary program runs
orchestrated from *outside* `makeinput.py` -- never by
`makeinput.py` launching a copy of itself.  `makeinput.py`
stays a plain input-writer: it reads whatever type
assignment the skeleton (`imago.skl`) already gives it
(the per-element species tags `si1`, `si2`, ... that
`structure_control` honours under `use_file_species=True`)
and writes the matching `imago.dat`.  All bispectrum
reasoning lives in the orchestrator and reuses the
`BispecMatcher` methods of 8.9 (`parse_loen_output`,
`distance`, `representative`).

The orchestrator is `build_initial_potentials.py` for the
producer (5.7), and an outside script or a human for any
other caller.  The producer already runs full Imago jobs
through kaleidoscope (DESIGN 6), so each loen run is just
one more dispatched unit -- and kaleidoscope's run-reuse
cache (6.2.5) already avoids recomputing it.  No special
machinery in `makeinput.py` is required.

#### 5.10.1 The two situations

A Fortran-side descriptor is needed in exactly two cases,
and only one of them rewrites the skeleton:

1. **Assign types by bispectrum** -- defect structures,
   amorphous materials, nanostructures, and the like.
   These are non-crystalline: the `imago.skl` is already
   in space group P1, or `makeinput.py` effectively makes
   it so, so every atom is explicit.  The bispectrum
   fingerprints decide which atoms share a type.
2. **Witness fingerprint for an already-typed crystal**
   -- a crystalline reference whose types come from
   symmetry.  Here the bispectrum is purely
   *informational* (a witness record for the database,
   5.2.2); the types are not derived from it.
   **Crystalline systems are never grouped by
   bispectrum.**  The reason is sharper than "the types
   are already known."  Grouping works in P1: it
   reassigns every atom a type from its fingerprint,
   which means it must first drop the space group and
   treat the cell as having no symmetry.  But the
   k-point machinery folds the full Brillouin-zone mesh
   onto the irreducible wedge *using that space group*,
   so regrouping a crystal in P1 would silently discard
   the symmetry the k-point sampling depends on and
   corrupt the band structure.  The witness path
   therefore leaves both the symmetry and the skeleton
   untouched.

Case 1 is the only one that writes a new skeleton; case 2
harvests one fingerprint per existing type and leaves the
skeleton untouched.

#### 5.10.2 The sequence

1. **First makeinput.**  Run `makeinput.py` on the
   skeleton with no environment grouping; the types are
   whatever the skeleton already declares (for case 1,
   typically every atom its own type in P1).  This writes
   a provisional `imago.dat`.  Because this is a
   `-loeninput` build, the fingerprint match is skipped and
   each atom takes the per-element default-tagged entry
   (5.6.5 step 2's loen-build skip -> step 3); the potential
   is irrelevant anyway, the bispectrum being geometric.
   The skip is what makes no `-pot` necessary -- and, more
   than a convenience, is what keeps this build from
   invoking itself (below).  The `LOEN_INPUT_DATA` block of
   this `imago.dat` carries the matcher's
   `(twoj1, twoj2, ...)` parameters via
   `BispecMatcher.to_loen_input` (5.10.5).
2. **Run loen.**  Invoke `imago.py -loen` against that
   `imago.dat`, producing `fort.21` (5.10.3).  **loen is a
   standalone, geometry-only job**: it reads the structure
   and the `LOEN_INPUT_DATA` parameters, needs no basis and
   no potential, and never shares a run with an SCF or
   post-SCF pass.  `imago.py`'s job table therefore defaults
   BOTH bases to `no` for `-loen` (so the explicit `-scf no`
   the producer passes is redundant but harmless), and both
   `imago.py` and imago's command-line parser (before any pass
   runs) refuse a run that combines loen with any basis,
   naming the reason.  The
   refusal is not mere tidiness: imago's input readers index
   per-type arrays by the basis code of the pass being
   parsed, and `loen` parses as a post-SCF pass, so an SCF
   basis with no post-SCF basis hands the reader an index of
   zero; and a post-SCF basis after an SCF pass parses the
   input a second time onto arrays the SCF pass left
   allocated.  A standalone run has neither problem: with
   both codes zero the command-line reader substitutes a
   valid slot index for parsing only, and neither pass runs.
3. **Orchestrate on `fort.21`.**
   - *Case 1 (grouping).*  Requires a P1 skeleton --
     space group 1 (the mandatory `space` line reads
     `1_a`) and a unit (`1 1 1`) supercell; the
     orchestrator refuses a symmetry-bearing skeleton
     here rather than silently drop its space group
     (5.10.1, 5.10.4).  Bucket the atoms by
     `BispecMatcher` fingerprint distance (the species
     logic of 5.6.4, run in the orchestrator rather than
     inside makeinput), then **rewrite the skeleton** with
     explicit per-element species tags (5.10.4).
   - *Case 2 (witness).*  Take one fingerprint per
     existing type (any member -- symmetry makes them
     identical) and attach it as a witness record (5.2.2).
     No skeleton rewrite.
4. **Second makeinput (case 1 only).**  Run `makeinput.py`
   again on the rewritten skeleton; it reads the new
   per-element species tags and writes the final, grouped
   `imago.dat`.  The run then proceeds (SCF and harvest,
   for the producer).

The two makeinput runs the orchestrator sequences (steps 1
and 4) are separate, ordinary processes passing only the
skeleton between them, so the *orchestration* does not
recurse.  There is a subtler self-invocation to guard,
though, and step 1's loen-build skip (5.6.5, above) is the
guard.  Step 1's build could otherwise reach the fingerprint
match of 5.6.5, and in the file-dictated regime that match
runs a loen descriptor computation whose own first step is a
makeinput build exactly like step 1 -- a genuine self-call.
The skip closes it at the source: a `-loeninput` build never
performs the match, so it cannot reach back into itself.
Without the skip the build recurses without end the moment
the database carries a preferred bispectrum record for the
element -- each level nesting a scratch directory inside the
last until the path length overflows.

#### 5.10.3 fort.21 carries its own identity

So that the orchestrator can map a `fort.21` row to the
atom and type it describes without a separate
cross-reference file, `fort.21` is **self-describing**:
each row leads with the site's identity and then the
descriptor.  The columns are:

  Column           Meaning
  -----------------------------------------------------
  site#            Engine potential-site index
                   (element-sorted "dat" order).
  element          Element symbol of the site.
  species          Per-element species index (the
                   `si1` / `si2` number).
  type_in_species  Per-element-species type index.
  type_flat        The single global type index
                   (1, 2, 3, ...) the engine derives by
                   element-sorted expansion.
  components       The `twoj2 + 1` real bispectrum values
                   (coupling channels j in the triangle
                   range; twoj1 >= twoj2).
  sum              Trailing sum column (ignored by the
                   matcher).

With these columns the mapping the orchestrator needs --
row -> (element, species, type) -- is read straight off
the file, which is more robust than re-deriving it through
`datSkl.map`; this codebase has a history of cross-file
numbering bugs, and a self-describing output removes that
class of error here.

*Implementation note.*  Today's `fort.21` carries only
`site#`, the `components`, and `sum` (loen.f90, the
write block opening `open(unit=21,...)`); it also writes
a header line.  Adding the `element` / `species` /
`type_in_species` / `type_flat` columns is a Fortran
change to the loen writer, which already has each site's
identity in hand.  `BispecMatcher.parse_loen_output` (8.9)
must then skip the header line and read the identity
columns -- the C55 first cut assumed a bare
components-plus-sum row with no header and must be revised
to the real format before any live use.

#### 5.10.4 Writing the new skeleton (case 1)

This rewrite runs only for a P1 skeleton -- no space
group line and a unit (`1 1 1`) supercell -- because case
1 is by construction non-crystalline (5.10.1).  The
orchestrator enforces this as a hard precondition: handed
a symmetry-bearing skeleton it refuses to rewrite rather
than discard the space group, which would corrupt
Brillouin-zone k-point folding.  A crystalline reference
belongs on the witness path (case 2), which leaves the
skeleton untouched.

When the orchestrator rewrites the skeleton with its
bispectrum grouping, it must follow the **per-element
species numbering** convention the engine expects: each
element's species restart at 1 -- e.g. `Si1, Si2, Si3`
then `O1, O2, O3` -- *not* a single run-on sequence like
`O4, O5, O6`.  The engine builds its flat type list
(1, 2, 3, ...) by element-sorted expansion of these
per-element species, so a run-on numbering would corrupt
the derived types.  A round-trip test -- group, write the
skeleton, reread it, and recover the same per-element
grouping -- guards this.

#### 5.10.5 Producer vs consumer parameter mapping

Both the producer and the consumer must map a
`(method, sub_spec)` pair into the loen-side input
parameters that `O_Input::readLoEnControl` reads.  The
mapping table lives in the matcher protocol (8.9):
each matcher exposes a `to_loen_input(sub_spec)`
method returning the parameter dict the
`LOEN_INPUT_DATA` block of `imago.dat` expects.

**Bispectrum parameter contract.**  The
`LOEN_INPUT_DATA` block as read by `O_Input` carries
the following parameters; `BispecMatcher.to_loen_input`
must populate all of them:

  Parameter      Source                Notes
  -----------------------------------------------------
  loenCode       Matcher constant      Selects the
                                       descriptor
                                       algorithm.
                                       `BispecMatcher`
                                       returns `1`
                                       (bispectrum
                                       component path
                                       in `loen.f90`).
  twoj1          `sub_spec["twoj1"]`   Integer; the
                                       larger of the
                                       two angular-
                                       momentum
                                       parameters.
                                       `O_Input` swaps
                                       internally to
                                       ensure twoj1 >=
                                       twoj2.
  twoj2          `sub_spec["twoj2"]`   Integer.  The
                                       output vector
                                       length is
                                       `twoj2 + 1` (the
                                       coupling channels
                                       j in the triangle
                                       range).
  max_neigh      `sub_spec.get(        Integer cap on
                  "max_neigh", 50)`    the per-site
                                       neighbor list
                                       length; sized
                                       for the cutoff
                                       reach below.
  cutoff         `sub_spec.get(        Real radial
                  "cutoff", 9.0)`      cutoff in Bohr
                                       on the neighbor
                                       list; wide
                                       enough to
                                       enclose every
                                       atom's first
                                       shell (see
                                       below).
  angleSqueeze   `sub_spec.get(        Real angular
                  "angle_squeeze",     compression
                  0.85)`               factor (see
                                       `loen.f90`
                                       notes on
                                       `angleSqueeze`).

The required `sub_spec` keys are `twoj1` and `twoj2`.
The remaining three are optional, with the defaults
shown -- the database-wide values `makeinput.py`
emits today.  The `cutoff` default of 9.0 Bohr (about
4.76 Angstrom) is chosen to enclose the first
coordination shell of *every* atom, including large,
loosely bonded cations whose first shell sits farther
out.  A cutoff too small to reach an atom's first
shell leaves that atom with an empty neighbor list and
an all-zero descriptor that carries no information and
cannot match anything in the database; 9.0 Bohr avoids
that failure mode across atom sizes.  The `max_neigh`
default of 50 caps the per-site neighbor list and must
be large enough for that reach -- a dense first shell
within 9.0 Bohr can hold a few dozen neighbors, and
the loen neighbor list has no internal bound check, so
an undersized cap would overrun its arrays.  (The
principled long-term answer to the single-global-cutoff
limitation is the element-aware cutoff, TODO C62 / D10.)
Two fingerprints whose `sub_spec` differs in *any* of
these five values produce different bispectrum vectors
and must coexist as separate fingerprint records per
DESIGN 5.2 rule 8.
When comparing `sub_spec`s for rule-8 uniqueness, the
canonical form is the post-default-resolution dict
(omitting a key is equivalent to specifying the
matcher's documented default for that key).

When element-aware bispectrum lands (TODO D10 / C62),
`BispecMatcher.to_loen_input` gains a `bispecByElement`
flag from `sub_spec.get("by_element", False)`.

Centralizing the mapping in the matcher keeps producer
and consumer aligned by construction.

#### 5.10.6 One declaration set, checked before dispatch

The producer touches the fingerprint declarations twice,
at opposite ends of a run.  Early, the **build** side asks
which Fortran-side declarations exist, because each one
needs its own dispatched `-loen -scf no` unit (5.10).
Late, the **harvest** side asks the same question again,
because each declaration becomes a stored fingerprint
record and a Fortran-side one is read from that unit's
descriptor.  The two answers must be the same answer.  If
the build's set is ever smaller than the harvest's, the
harvest reads a descriptor no run produced.

That failure has happened, and its shape is worth keeping
on record because it is the shape any future divergence
will take.  The harvest read its set as the database-wide
`[characterization]` recipe plus the environment's own
overrides; the build read only the entries' overrides.
While every manifest carried its bispectrum as a per-entry
declaration the two agreed by accident.  When the
characterization block became the normal home for the
recipe (5.7), a manifest with a recipe and no overrides
left the build with nothing to dispatch, and the producer
ran every solid's full convergence sweep before dying at
harvest on a descriptor that was never computed.

**One rule, two consumers.**  The composition of a
declaration set -- the recipe, then that environment's
overrides -- is defined once and used by both sides.  The
harvest applies it to one environment.  The build cannot
know the environments yet, since they are discovered from
the converged run (5.7), so it applies the same rule to
every entry the manifest declares *and* to the
override-less case, and takes the union.  The build's set
is therefore a superset of any set the harvest can later
present, by construction rather than by inspection.  This
is what makes the agreement structural: a future third
source, or a change in precedence, is written in one place
and both sides move together.

Two consequences are worth naming.  The union may build a
unit for a declaration no environment turns out to claim
-- a site-less customization, for instance, which the
harvest's environment discovery cannot yet match (5.7).
That is a wasted geometry-only run, which is cheap, and it
is the right direction to err: an extra descriptor costs
one short run, a missing one costs the whole flight.  And
the union deduplicates by calc tag, so one run still
serves every environment sharing a `(method, sub_spec)`.

**Checked before anything is dispatched.**  Structural
agreement is an argument, and an argument can be wrong.
The producer therefore also asserts the invariant directly:
once the units are assembled and before any is sent, every
Fortran-side declaration the harvest could read must have a
matching loen unit in the flight.  The check is a
set comparison over calc tags, costing nothing, and it
converts the entire class of build/harvest drift from a
failure discovered after minutes of cluster SCF time into
one raised before a single job is submitted, naming the
solid and the `sub_spec` that has no run.

This backstop keeps its value even though the single rule
above should make it unreachable.  It guards what the rule
does not: a unit that was composed correctly and then
dropped during flight assembly, or a filter that removes
units for reasons of its own.  A cheap invariant that can
only fire when something else is already broken is exactly
the invariant worth asserting -- it costs nothing when the
code is right and saves a whole run when it is not.

### 5.11 The reduce descriptor

Where 5.10 defines the path a Fortran-side descriptor
takes, this section defines what the Python-side one *is*.
The reduce descriptor answers "what does this atom see
around it?" with concentric spherical shells: for each
level, the distance out to that shell and the multiset of
elements sitting in it.

**The geometry it reads is a periodic neighbour list.**
For an atom, every periodic *image* of every atom that
falls within the `cutoff` is a neighbour, counted once per
image.  This includes images of the central atom itself,
which are ordinary neighbours in space -- in an
face-centred cubic lattice the entire second shell of a
site consists of images of that site.  Only the atom at
distance zero is excluded, since an atom is not its own
neighbour.

That definition is what makes the descriptor
**transferable**, which is the property 5.2 relies on when
it stores a shell code in one structure and matches it
against another (5.6.5).  A count of neighbours is a
property of the environment; it must not change because a
curator chose a different but equivalent cell for the same
material.

**The walk.**  Shells are built outward, one level at a
time:

1. Seed the level at the closest neighbour not yet
   assigned to a shell.
2. Sweep every neighbour whose distance falls in
   `[seed, seed + thick]` and within `cutoff` into this
   level.
3. Repeat for `level` levels.

Each shell records the seed distance and its neighbours.
Within one structure the neighbour multiset carries
`(element, species)`, because species distinguishes atoms
there; the multiset *stored* in the database carries
element symbols only, since species numbering is local to
a structure and would not transfer (5.2).

**Exhaustion is refused, not padded.**  A cutoff too small
to reach the requested number of levels leaves a level with
no neighbour to seed it.  That is refused, naming the level
and the cutoff, because an empty shell is a value the walk
did not find and inventing one would store a descriptor no
structure produced.

#### 5.11.1 Why this replaced a cell-atom walk

The descriptor was originally computed over the *atoms of
the cell*, using a minimum-image distance matrix: one
entry per central-cell atom, holding the shortest distance
to any of its images.  Distances were therefore
periodic-correct, but *multiplicity* was capped at the
number of atoms the cell happened to contain.  Diamond
silicon has four nearest and twelve second neighbours; its
eight-atom conventional cell reported four and three,
because the twelve second neighbours are images of only
three distinct cell atoms.  Its two-atom primitive cell
reported one neighbour and could not form a second shell
at all.

Two atoms of the same material in different cells thus
received different descriptors, which is precisely the
property a transferable fingerprint must not have.  The
defect was invisible for years because of where the
descriptor came from: it began as a *grouping* tool for
large disordered models, and in a cell big enough that
every neighbour is a distinct atom, the two definitions
agree exactly.  It was only when the same code was asked
to fingerprint small crystalline reference cells that the
cap began to bite.

**Grouping is unaffected by the change**, which was
measured rather than assumed before making it.  Grouping
compares atoms *within* one structure, and the old cap
truncated every symmetry-equivalent atom identically, so
the relative comparison survived even where the absolute
counts were wrong.  Across a 1296-atom amorphous silica
model and four small crystals, both walks produce
identical species partitions -- on the glass they agree
atom-for-atom on the shells themselves, while on an
eight-atom diamond cell the shells differ completely
(second-level counts of three against twelve) and the
partition is *still* identical.  So the correction changes
what the descriptor says without changing what grouping
does with it.

---

## 6. High-Throughput Calculation Flights

This section holds the algorithm- and contract-level
designs for VISION Goal 4, whose architecture is laid
out in ARCHITECTURE 9 (the layering, the four
scripts/packages, and the load-bearing VISION
principles).  Where ARCHITECTURE 9 says *what* each
layer is and *why* the boundaries fall where they do,
this section says *how* each one behaves in enough
detail to implement.

The designs land incrementally, one subsection per
TODO item, in dependency order: 6.1 is the `imago.py`
callable API (TODO D11), the foundation every higher
layer reaches through; the kaleidoscope dispatcher (D13),
the ASE adapter (D12), and structure acquisition (D14)
follow in later subsections.  Note the section-number
offset from ARCHITECTURE: DESIGN 6 corresponds to
ARCHITECTURE 9, just as DESIGN 5 corresponds to
ARCHITECTURE 8.  The mapping is by name and
cross-reference, not by number.

### 6.1 imago.py callable API

This subsection designs the refactor of ARCHITECTURE
9.2: turning today's command-line-only `imago.py`
driver into a callable Python API, with the CLI reduced
to a thin wrapper over it.  The API is the single seam
every higher layer (the ASE adapter, kaleidoscope, and
through it the database producer and the bench harness)
reaches through, so its contract is designed first and
deliberately Imago-native and dependency-free.

#### 6.1.1 What the first client needs

The API's first real client is the initial-potential
database producer (`build_initial_potentials.py`, C48,
running *through* kaleidoscope, ARCHITECTURE 9.7).
Designing the result object against that client's
concrete needs keeps the contract shaped by use rather
than by guesswork.  The producer, per the harvest
contract settled in 5.7 and ARCHITECTURE 9.7, needs the
API to tell it:

1. **Did the SCF converge?**  A clear converged /
   not-converged / failed verdict.  A non-converged or
   crashed reference run must never be harvested into
   the database, and per VISION Principle 10 the
   flight must learn this as *data* (a result it can
   record and skip), not as an exception that aborts
   the whole batch.
2. **Where is the converged potential?**  An absolute
   path to the converged `scfV` output file (the
   `<edge>_scfV-<basis>.dat` that today's `manage_output`
   writes from `fort.8`).  That file carries the potential
   for *every* OLCAO potential type in the material, not one
   bare coefficient block: a `NUM_TYPES` header, then a
   `TOTAL__OR__SPIN_UP` channel listing each type as a count
   line plus that many Gaussian-term lines (a redundant
   `SPIN_DN` channel follows -- the producer runs non-spin,
   so the `TOTAL__OR__SPIN_UP` channel *is* the total
   potential and the `SPIN_DN` copy is ignored; spin handling
   is deferred).  The harvest selects the named site's type
   block -- its type number read from `datSkl.map` (9.7) --
   and takes columns 1 and 2 of each term as the coefficient
   and its alpha.  Those alphas equal the basis input the
   producer fed makeinput -- the consistency "converged
   `scfV` matches input `scfV`" (5.7) names -- so coefficients
   and alphas are read together from the one converged block.
3. **Under what conditions did it run?**  The SCF
   settings actually used -- basis, k-point spec,
   convergence threshold, Imago build commit -- so the
   producer can fill the provenance fields of 5.2, and so
   kaleidoscope can form its run-reuse cache key (the SCF
   convergence limit plus the resolved structure bytes;
   ARCHITECTURE 9.6).  The build commit is wanted here for
   the provenance record only -- it is written down beside
   each run and never compared (6.2.5).
4. **How much work did it take?**  The SCF iteration
   count, both for the producer's run log and for the
   20%-iteration-reduction validation harness (5.8).

These four needs map directly onto the result object's
fields (6.1.2).  Nothing here is producer-specific in a
way that pollutes the contract: every field is a
plain fact about an Imago run that any client (a
convergence sweep, the ASE adapter reporting `energy`,
a future AIMD step) would also want.

#### 6.1.2 The result object

A single small, immutable result object is returned by
both entry modes (6.1.3).  It is a plain dataclass with
no ASE, Parsl, or makeinput dependency -- an
Imago-native record of one run's outcome.

```
ImagoResult
  status            RunStatus enum (see below)
  success           bool: True iff status is CONVERGED
  run_dir           absolute path to the run directory
                      (the project home, where named
                       outputs are written)
  temp_dir          absolute path to the intermediate
                      (IMAGO_TEMP-mirrored) working dir
  scf_iterations    int | None: SCF cycles to reach the
                      convergence threshold; None when
                      no SCF ran (e.g. -scf no) or when
                      the count could not be parsed
  converged         bool: SCF reached its threshold
                      (distinct from "ran without
                       crashing"; see status)
  reused_checkpoint bool: True if within-run-dir
                      checkpointing short-circuited some
                      or all of the work (6.1.5)
  total_energy      float | None: harvested total energy
                      in Hartree, when available; the
                      ASE adapter (D12) converts to eV
  measured          MeasuredQuantities | None: scalar
                      electronic-structure quantities plus
                      the resolved k-point mesh, harvested
                      from the converged SCF output (see
                      below).  None for
                      runs that did not converge or for
                      job types that do not compute
                      them.  Used by the guidance-
                      dataspace harvest (DESIGN 7.8)
  outputs           dict[str, str]: logical name ->
                      absolute path for each output file
                      produced (e.g. "scfV", "energy",
                      "iteration", plus property-specific
                      keys like "tdos", "bond"); the
                      producer reads outputs["scfV"].  Also
                      "structure" -> the run's expanded
                      full-cell skeleton (makeinput's
                      imago.fract-mi: every atom explicit,
                      space group 1, at the run's sorted
                      numbering), and "datSkl_map" -> the
                      sorted<->skeleton atom map.  The
                      reduce fingerprint harvest (5.7) reads
                      this pair: the expanded skeleton gives
                      the geometry its shells need, and the
                      map turns a manifest atom_site
                      (skeleton numbering) into the row of
                      that structure.  Both are present only
                      when the run went through makeinput
  job               echo of identity: edge, job_name,
                      basis_scf, basis_pscf
  runtime_seconds   float: wall-clock time of the run,
                      for kaleidoscope's status.toml
  message           human-readable summary or error text

MeasuredQuantities
  gap_ev               float | None: band gap in eV,
                         read from the eigenvalue
                         spectrum; 0.0 for metals; None
                         when not computed (closed-shell
                         molecular runs that skipped the
                         analysis, unsupported job types)
  gap_kind             str | None: "direct" / "indirect" /
                         "none"; "none" iff gap_ev == 0.0
  spin_polarization    float | None: fractional spin
                         polarization at the Fermi
                         level for metals; 0.0 for
                         closed-shell non-magnetic
                         systems
  total_magnetization  float | None: total magnetic
                         moment per formula unit in Bohr
                         magnetons; 0.0 for non-magnetic
  kpoint_mesh          list[int] | None: the resolved
                         axial k-point counts
                         [n_a, n_b, n_c] that imago
                         selected for this run
                         (PSEUDOCODE 4c.6).  This is the
                         mesh's identity: two runs with an
                         equal kpoint_mesh integrated over
                         the same k-points and are the
                         same calculation.  The k-density
                         ladder guard (7.8 step 3c) keys
                         on it.  None for job types that
                         build no mesh
  kpoint_count         int | None: the number of k-points
                         actually computed -- the IBZ size
                         for a symmetry-reduced run, the
                         full-mesh size otherwise
                         (PSEUDOCODE 4c.6).  The resource
                         dataspace's size signature (8.2).
                         None when no mesh was built
```

The `measured` block is the C76 follow-up: a small
extension to imago.py's post-SCF analysis path that
exposes these quantities through the callable API.
Their primary consumer is the guidance-dataspace harvest
(DESIGN 7.8), but they are also of general interest to
any client that wants to inspect what the calculation
produced without re-parsing Imago's native output files.
Every field is optional so a job type that does not
compute it (e.g. a band-structure-only run) cleanly
omits it.

`RunStatus` is an enum with four members, chosen so the
flight layer can branch on outcome without parsing
`message`:

- `CONVERGED` -- the run completed and the SCF reached
  its convergence threshold.  The only status for which
  `success` is True and `outputs["scfV"]` is safe to
  harvest.
- `NOT_CONVERGED` -- the run completed (the Fortran
  binary exited cleanly, the `fort.2` success file was
  written) but the SCF hit its iteration ceiling
  without converging.  Outputs exist but must not be
  harvested as a reference potential.
- `FAILED` -- the run did not complete: the Fortran
  binary aborted, the `fort.2` success file was absent,
  or a required input was missing at run time.  This is
  an *expected* run-level failure, returned (not
  raised) so the flight can record-and-continue.
- `SKIPPED` -- there was nothing to do because
  within-run-dir checkpointing found the requested work
  already complete (6.1.5).  `success` is True;
  `outputs` point at the pre-existing files.

The `CONVERGED`-vs-`NOT_CONVERGED` verdict needs a signal
the current driver does not produce: today it checks only
the `fort.2` success file, which certifies the binary ran
without an abortive error, *not* that the SCF converged.
The P6 pseudocode pass surfaced this gap and resolved it
with no new Fortran signal (PSEUDOCODE 12.5): read the
iteration file's last data row and compare its
convergence metric (column 4) against the
`CONVERGENCE_TEST` criterion in the run's own `imago.dat`;
converged iff it is below the criterion.  The same row
also yields the last iteration's total energy (column 5)
and -- because column 1 is a per-run cycle counter that
resets each SCF invocation -- the `scf_iterations` count
robust to the file's append-on-rerun behavior.

C76 widens this same row to a fixed eight columns so a
plain SCF run also surfaces the electronic-structure
signal the guidance harvest (C72) needs: column 6 is the
magnetization (0.0 and always emitted for a non-spin run,
the total magnetic moment otherwise), column 7 the raw
band gap in Hartree, and column 8 an integer gap-kind
code (0 = metal/no gap, 1 = direct, 2 = indirect) that
the parser maps via `GAP_KIND_BY_CODE`.  The writer
prefixes every field with an explicit blank so adjacent
values can never abut even when one fills its full field
width, keeping the columns whitespace-delimited for the
`split()`-based reader.  Metal detection uses a dedicated
cutoff `metalGapThresh = 1.0e-3` a.u. (about 0.027 eV, of
order room-temperature kT), deliberately far larger than
the 1e-8 numerical-degeneracy threshold: a true metal
sampled on a discrete k-point mesh shows a small finite
gap on the order of the level spacing at the Fermi energy
(~1e-4 to 1e-2 a.u.), and the kT-scale cutoff collapses
these mesh artifacts to a zero-gap metal while staying
well below any genuine semiconductor gap.  A coarse-mesh
metal whose artificial gap exceeds the cutoff is a
k-point convergence problem to cure with a denser mesh,
not a reason to raise the threshold.

The boundary on error handling is deliberate and
important for Principle 10.  *Run-level* failures
(non-convergence, a Fortran abort, a missing input
file) are reported as a returned `ImagoResult` with the
appropriate status -- they are normal outcomes of
running real calculations and must not abort a
flight.  *Contract* failures (the environment is not
configured: `$IMAGO_RC`/`$IMAGO_TEMP`/`$IMAGO_BIN`
unset; the named run directory does not exist or holds
no inputs; the lock file is already held by another
process, 6.1.4) raise an `ImagoError`.  These are
programmer or environment errors that no per-job retry
can fix, so they propagate.  Today's `imago_exit`
(which prints to the runtime log and calls `sys.exit`)
is replaced inside the API path: it must never call
`sys.exit`, because that would kill the long-lived
kaleidoscope worker driving many runs.  The thin CLI
wrapper (6.1.3) is the only place a process actually
exits.

#### 6.1.3 The two entry granularities

The API offers two entry points, so a caller joins at
whichever level it already has inputs for (ARCHITECTURE
9.2).  Both funnel into one private core (6.1.4) and
both return an `ImagoResult`.

- **`run_prepared(run_dir, *, settings=None) ->
  ImagoResult`** -- *prepared-directory mode.*  The
  given `run_dir` already holds the staged Imago inputs
  (`imago.dat`, `structure.dat`, `scfV.dat`, kp files),
  produced by makeinput or a prior step.  No makeinput
  call is made.  `settings` carries the run options
  (job type, bases, edge) that today come off the
  command line; when omitted, the same resource-control
  defaults apply as for a bare CLI `imago` invocation.
  This is the mode kaleidoscope's default wingbeat uses,
  because kaleidoscope (or makeinput, dispatched by it)
  has already built the directory.
- **`run_structure(structure, options, run_dir, *,
  settings=None) -> ImagoResult`** --
  *structure-and-options mode.*  Given a structure and
  a set of makeinput options, the API drives
  `makeinput.py` to build `run_dir` first, then calls
  `run_prepared` on it.  Input *preparation* still
  lives in makeinput; this mode simply calls it on the
  caller's behalf.  `structure` is, at this design
  stage, a path to an `imago.skl`; whether it may also
  be an in-memory `StructureControl` is deferred to
  D12/C64 (the ASE-free factory), so this contract does
  not yet commit to it -- see 6.1.6.

The **CLI wrapper** is the third, outermost layer and
the only one that touches `sys.argv` or exits the
process.  Today's `main()` (and the argv-bound
`ScriptSettings` it constructs) is split into three
responsibilities:

1. Parse `sys.argv` into run options (the existing
   argparse surface and `reconcile` logic, unchanged in
   meaning).
2. Decide the entry mode.  A bare `imago ...` on a
   directory that already has inputs is
   prepared-directory mode -- the overwhelmingly common
   case and today's only behavior -- so the CLI calls
   `run_prepared` on the current working directory.
   (A future CLI surface for structure-and-options mode
   is possible but is not required by D11; the CLI's
   job here is simply to keep doing what it does today,
   now through the API.)
3. Translate the returned `ImagoResult` into a process
   exit code: `CONVERGED`/`SKIPPED` -> 0;
   `NOT_CONVERGED`/`FAILED` -> non-zero, with `message`
   written to the runtime log; an uncaught `ImagoError`
   -> a non-zero exit with its message.

This split is what lets `ScriptSettings` stop being
constructed from `sys.argv` unconditionally.  Its
`reconcile` method already takes an `args` namespace
and contains all the job-type/edge/basis resolution;
the refactor separates *building* that namespace (from
argv, in the CLI; or from a plain options mapping, in
the API) from *reconciling* it into a settings object.
The argv-only side effects in today's constructor --
`recordCLP`, which appends the literal `sys.argv` to a
`command` file -- become CLI-only: in API mode there is
no meaningful argv to record, so the API instead
records the equivalent call provenance (entry mode,
run_dir, options) or skips the `command` file
entirely.  This is flagged as an open detail in 6.1.6.

#### 6.1.4 The private run core, and cwd discipline

Both entry modes converge on one private core that
performs the sequence today's `main()` runs inline:
resolve directories, acquire the lock, stage inputs,
build the job command line, execute the Fortran binary
(plus any immediate secondary jobs -- SYBD post-pass,
Kramers-Kronig for optical properties), collect and
rename outputs, parse the result fields, release the
lock, and return the `ImagoResult`.  The behavior is
identical to today's flow; the change is that it
returns a value instead of falling off the end of
`main()`, and reports failure by status instead of
`sys.exit`.

One genuine behavioral difference from the CLI must be
designed in: **current-working-directory discipline.**
Today `imago.py` is a one-shot process: `main()` does
`os.chdir(temp)` and never restores the cwd, which is
harmless because the process exits immediately after.
A kaleidoscope worker, by contrast, is a long-lived
process that may drive many run directories in
sequence.  The API core therefore must treat cwd as a
resource to acquire and release: it takes `run_dir`
explicitly (rather than implicitly trusting the
caller's cwd), changes into the working directory for
the duration of the run, and **restores the original
cwd on exit, including on failure** (a `try/finally` or
context manager).  Without this, one failed run would
leave a flight worker stranded in a stale temp
directory and corrupt every subsequent run's relative
path resolution.  This is the single most important
correctness difference between the CLI's one-shot
assumption and the API's reentrant requirement.

#### 6.1.5 Lock-file and checkpoint behavior preserved

Both existing robustness mechanisms carry over
unchanged in *meaning*; the design only clarifies how
they behave under concurrent, in-process use.

**The lock file is already per-run-directory, so
flight concurrency is safe by construction.**  Today
the lock (`imagoLock`) lives in the `temp` directory,
which `get_temp_dir` derives by mirroring the run
directory's path under `$IMAGO_TEMP`.  Two different run
directories therefore mirror to two different temp
directories and two different lock files.  A
kaleidoscope flight running thousands of independent
SCFs in parallel -- each in its own run directory --
takes thousands of independent locks that never
collide.  The API keeps the exact same acquire / mark /
release lifecycle (create on entry, stamp with the run
label, remove in the cleanup step).  The one contract
change: encountering an *already-held* lock is a
contract failure in API mode and raises `ImagoError`
(another process owns this run directory -- the caller
must not have dispatched two runs into the same
directory), whereas the CLI prints its existing
"Is another imago script running?" message and exits
non-zero.  The lock guards a single run directory; it
is never a process-global or flight-global lock.

**Checkpointing stays within the run directory and is
orthogonal to kaleidoscope's coarser cache.**  Imago's
internal checkpointing -- skipping completed SCF
integrals on restart, skipping a basis SCF that is
already complete when another job needs the same basis
-- is driven by `manage_input`'s staging logic plus the
Fortran binaries, and is untouched by this refactor.
The API runs the same `manage_input`, so a re-entered
run directory resumes exactly as a re-run CLI invocation
would.  The result object surfaces this with
`reused_checkpoint` (some work was short-circuited) and,
in the limiting case where *all* requested work was
already complete, `status = SKIPPED`.  This is a
deliberately different and finer-grained thing from
kaleidoscope's run-reuse cache (ARCHITECTURE 9.6): the
API/`imago.py` decides whether to redo work *within* a
run directory; kaleidoscope decides whether to *launch*
the run directory at all.  The clean statement of the
boundary: `imago.py` resumes within a run; kaleidoscope
launches or skips runs.  Designing kaleidoscope's
side of that boundary is D13's job, not D11's.

#### 6.1.6 Open details (for PSEUDOCODE / implementation)

These are deliberately deferred to the PSEUDOCODE pass
for D11 or to C63 implementation; none of them changes
the contract above.

- **Output-key enumeration.**  `manage_output` renames
  `fort.*` files by a job-type-specific (`jobID % 100`)
  scheme.  The exact set of logical keys in
  `outputs{}` per job type must be enumerated in
  pseudocode, factored out of the existing per-property
  `_manage_*_output` helpers so the API and the file
  layout cannot drift apart.
- **Parsing `scf_iterations` and `total_energy`.**
  Whether to count iterations and read the energy from
  the named output files (`<edge>_iter<basis>.dat`,
  `<edge>_energy<basis>.dat`) after the run, or to
  capture them from the Fortran stdout already written
  to the runtime log, is an implementation choice.
  Reading the settled output files is the more robust
  default and is the working assumption.
- **`run_structure` structure type.**  Whether
  `structure` may be an in-memory `StructureControl` in
  addition to an skl path depends on the ASE-free
  factory of D12/C64; D11 commits only to the skl-path
  form and leaves the richer signature to land with
  that work.
- **Call-provenance recording in API mode.**  What
  replaces `recordCLP`'s `command`-file append when
  there is no `sys.argv` (record the API call shape, or
  skip the file) is an implementation detail with no
  bearing on the returned contract.

### 6.2 kaleidoscope flight dispatcher

This subsection designs ARCHITECTURE 9.4 and 9.6: the
Parsl-based package that drives a *set* of Imago
calculations, tracks their outcomes, and resumes over
work already done.  It builds directly on 6.1 -- the
default unit of work is an `imago.py` callable-API call,
and the result it persists is the 6.1.2 `ImagoResult`.
It also **resolves the workspace-scheme open question of
ARCHITECTURE 9.8** (the stable-id convention, the
`<calc>` tag format, and the `status.toml` schema are
pinned in 6.2.4).

The governing constraint is VISION Principle 9:
kaleidoscope is *ordinary scientific Python* and stays
free of materials-specific coupling.  It dispatches,
tracks, and caches; it does not know what an SCF or a
potential is.  Everything domain-specific lives either
below it (the wingbeat, 6.2.2) or above it (the client,
6.2.6).  Three other load-bearing principles shape this
design.  Principle 8 keeps the wingbeat seam independent of
the execution adapter.  Principle 10 (complete-and-report)
ensures one failed unit never aborts the flight.  And
Principle 12 (the flight layer stays dumb; flight
description lives in Python) is the choice that
kaleidoscope never grows a flight description language
-- no DSL, no workflow grammar, no DAG engine.  The
`Flight` data model (6.2.1) is a flat list of
independent units; higher-order flight shape (multi-axis
sweeps, dependent phases, per-unit iteration) is composed
in client Python that builds the flat list, or absorbed
inside a custom wingbeat that owns one unit's internal
iteration (6.2.2).

A practical corollary of Principle 12 is the
**flight-builder split**.  Domain-agnostic *option-axis*
sweeps (sweep k-density values, sweep target atoms for
XANES, sweep basis sizes) live as helpers inside
kaleidoscope -- 6.2.8 is the first such helper, the
k-point-density convergence constructor for DESIGN 7
(predict-then-verify is its default strategy).  These
builders live in a `kaleidoscope/builders/` subpackage (one
module per builder, each named for the axis it sweeps: the
first is `builders/kpoint_convergence.py`), imported by a
client so the dumb core's import graph never pulls the
physics layer (`guidance_db`, `structure_control`) that a
builder depends on.  Domain-aware
*structure-axis* sweeps (supercell expansion, LAMMPS-
snapshot per-frame splitting, defect-site enumeration)
generate skl files and therefore live in
`structure_control` and acquisition; kaleidoscope only
consumes the resulting structure paths.  Both halves
ultimately produce flat lists of `CalcUnit`s; neither
grows a DSL.

#### 6.2.1 The unit of work and the flight

Kaleidoscope's data model is two plain, domain-agnostic
records.

```
CalcUnit
  id            stable per-structure key (6.2.4); the
                  curation reference_id for the producer,
                  a COD id for an acquisition flight
  calc          tuple[str, ...] of per-axis directory
                  components (6.2.4).  The empty tuple
                  means "no second level"; a single-
                  element tuple is one calc tag; a multi-
                  element tuple is a nested-axis sweep
                  (one element per varied axis, in
                  Flight.sweep.varied_axes order)
  structure     path to an imago.skl (or a structure
                  handle the chosen wingbeat understands)
  options       makeinput options for this unit
  wingbeat      which Wingbeat executes it (6.2.2);
                  defaults to the flight default
  kind          run role (6.2.9): a short label the core
                  stores and round-trips but never
                  interprets.  Default "convergence"; each
                  harvester reads only the kinds it knows
                  (e.g. "fingerprint" for loen runs)
  key_fields    client-declared cache identity (6.2.5):
                  scalar fields + names of key files

Flight
  root          workspace root directory (6.2.4)
  units         list[CalcUnit]
  default_wingbeat  Wingbeat used when a unit names none
  parsl_config  the Parsl Config (deployment, 6.2.3)
  sweep         SweepRecord | None: records varied_axes
                  order + fixed_axes when the flight
                  was built by the predict-then-verify
                  helper (6.2.8); None for hand-built
                  flights that did not declare a sweep
  on_outcome    optional per-unit callback (6.2.6)
  metadata      dict[str, dict]: opaque per-key tables the
                  dispatch core round-trips verbatim into
                  flight.toml as [flight.<key>] blocks but
                  never interprets (Principle 9).  Default
                  empty; the k-point convergence builder
                  (6.2.8/6.2.9) stashes its per-structure
                  PredictionRecords as
                  metadata["predictions"][<id>]

SweepRecord
  varied_axes   tuple[str, ...]: axis names in the
                  canonical order they appear at each
                  level of CalcUnit.calc
  fixed_axes    dict[str, str]: axis -> value for axes
                  that take the same value across every
                  unit in the flight (recorded as
                  context, not as a calc-tag level)
```

The `CalcUnit.calc` tuple is the on-disk path
representation of the unit's sweep position.  Each
element is a `<axis>-<value>` directory component per the
6.2.4 naming rules; iterating the tuple gives the levels
top-to-bottom.  `pathlib.Path(unit.id, *unit.calc)` builds
the unit's run-dir relative to `root/wingbeats/`; reading
an existing path back, split on `/` to recover the tuple.
The shape is deliberately a tuple of strings rather than
a dict of (axis -> value): the flight's `sweep` field
already records the axis order canonically, and
duplicating the axis names per unit would invite drift
between them.  Future extensions (per-axis annotations,
floats with units) become a `tuple[CalcAxis, ...]` swap
without changing the path-building code.

A client builds a `Flight` in process -- kaleidoscope
is a library first (Principle 9), not a CLI -- and hands
it to the dispatch entry point.  Kaleidoscope serializes
the flight to `flight.toml` in the workspace root so
a flight is inspectable and a resume has an
authoritative record of *what was asked for*, separate
from `status.toml`'s record of *what happened* (6.2.4).
Whether `flight.toml` may also be hand-authored as the
primary surface, rather than always generated from the
in-process `Flight`, is left open (6.2.7).

The producer (C48.3) is the worked example throughout the
rest of 6.2.  Under DESIGN 7's predict-then-verify
workflow, the producer's relationship to kaleidoscope
changes shape: rather than launching one `CalcUnit` per
curated reference solid, the producer asks the
flight-builder helper (6.2.8) to expand each reference
solid into a small **verification sub-grid** of
`CalcUnit`s -- one per k-density value chosen by the
predict-then-verify algorithm of 7.7.  Every unit in that
sub-grid shares `id = reference_id`, the curated skl as
`structure`, the default (Imago) wingbeat, and the same
`key_fields` (the SCF convergence limit as the one scalar,
the structure file as a key file -- 6.2.5); they
differ in `calc` (the per-grid-point tag per 6.2.4) and in
the swept k-density value carried in `options`.

The harvest step (6.2.6) then walks each reference solid's
sub-grid, picks the converged grid point per 7.8 (the
smallest k-density at which consecutive grid points'
energy delta falls below a threshold), and reads its
converged potential from that one run dir.  A reference
solid whose sub-grid fails to converge at the top of the
range is skipped -- no potential harvested, no guidance
entry staged -- per 7.9's non-convergence recovery.  This
is the shape that replaces the single-`CalcUnit`-per-solid
sketch the earlier draft of 6.2 carried: the producer is
now a **predict-then-verify client** of both kaleidoscope
and the historical-guidance DB.

**Trust mode for nearly-identical families.**  Not every
flight warrants the verification sub-grid.  When the
curator already knows -- from prior work, from a recent
seed flight on a sibling solid, or from a high-similarity
match they trust -- that a particular k-density value is
the right operating point for an entire family of
nearly-identical reference solids, requiring every solid
to re-verify it is wasted compute.  The flight-builder
helper (6.2.8) therefore exposes a `verify=False` mode
that collapses the sub-grid to a single `CalcUnit` at the
trusted predicted k-density.  In that mode the producer
still harvests the converged potential from the one run
dir, but does *not* auto-stage a new guidance entry: a
single converged calculation is weaker evidence than a
converged grid (it confirms the value works, but does not
demonstrate that smaller values would not have worked
too), so the harvest path treats trusted runs as
contributing potentials but not new guidance.  A curator
who wants the trusted value reinforced into the DB can
stage it manually via a `source = "manual"` entry per
7.4 / 7.8.  This is a deliberate asymmetry: trust mode
*consumes* the guidance DB without *amending* it.

#### 6.2.2 The pluggable wingbeat seam

A *wingbeat* is the seam (Principle 8) that isolates
kaleidoscope's dispatch core from how a unit actually
executes.  It is a small protocol:

```
Wingbeat.run(unit, wingbeat_dir) -> WingbeatOutcome
```

The wingbeat receives a unit and the prepared run
directory, executes the calculation however it likes,
and returns a **domain-agnostic** `WingbeatOutcome`:

```
WingbeatOutcome
  ok        bool: did the unit complete (not "succeed
              scientifically" -- see detail)
  detail    short opaque string the wingbeat chooses and
              kaleidoscope records but never interprets
              (e.g. "converged", "not_converged")
  runtime_seconds  float
  message   human-readable text
```

The crucial layering decision: kaleidoscope tracks a
generic lifecycle status (6.2.4) and stores `detail`
verbatim *without interpreting it*.  ARCHITECTURE 9.4's
list of surfaced outcomes ("converged, non-converged,
cluster-side loss, post-processing error") is therefore
a deliberate split -- cluster-side loss is
kaleidoscope's own (a Parsl task that vanished, 6.2.3),
while converged / non-converged are wingbeat-supplied
`detail` strings.  This is what lets kaleidoscope
surface convergence in `status.toml` *and* stay ignorant
of what convergence means.

- The **default wingbeat** (`ImagoWingbeat`) drives the 6.1
  API.  It first **stages the run directory**, which has three
  cases because the wingbeat is generic and only *some* clients
  prepare.  When the unit names a staging area -- the producer,
  whose driver-side prepare step (6.2.5) already built the
  inputs there -- the wingbeat does not build: it **commits the
  staged inputs** into the run directory.  When the run
  directory already holds staged inputs (a re-launch of a
  directory a prior run built) there is nothing to do.  Only
  when neither holds -- a client that never prepared -- does the
  wingbeat **build** the inputs itself with makeinput, from the
  unit's structure and its makeinput-side options.  In every
  case it then runs the prepared directory with
  `imago.run_prepared`.  It
  STILL **re-applies the unit's imago-side settings** on that
  run -- it partitions the unit's options (6.2.10) and passes
  the imago-side ones (`job` / `edge` / `scf_basis`) to
  `run_prepared`, because those are imago *runtime* options,
  not baked into the staged `imago.dat`.  A pre-staged
  directory carries no memory of the job it was built for, so
  the settings must travel with every invocation.  The job
  type and the SCF suppression live only in these settings, so
  if they are dropped imago no longer sees the unit's
  `-loen -scf no` request and falls back to its default job, a
  ground-state SCF.  (`-loen -scf no` never runs an SCF itself;
  the unwanted SCF is purely the dropped-settings fallback.)
  It maps the returned `ImagoResult`
  into a `WingbeatOutcome`: `ok =
  status in {CONVERGED, NOT_CONVERGED, SKIPPED}` (the binary
  *ran*), `detail = status.name.lower()`.  It also **persists
  the full `ImagoResult` into the run directory** as
  `result.toml` (6.2.6), so the Imago-native detail
  survives for the client to reload without
  kaleidoscope ever parsing it.  Alongside the measured
  fields it echoes one recorded fact into that file: the build
  identity from the unit's `record` mapping (6.2.4), written as
  `imago_commit`.  A guidance entry's provenance reads it there
  (7.8), which keeps that harvest on the one per-run file it
  already reads instead of opening the core's `status.toml`.
  The echo is where TODO C84 lands: once the binary reports its
  own build, the wingbeat prefers the engine's word and the
  recorded value becomes the fallback.
- An **ASE wingbeat** wraps `ImagoCalculator` (D12) for
  units that need ASE-MD or ASE-relaxation semantics; it
  too ultimately calls the 6.1 API underneath.
- A single flight may **blend wingbeats** per unit, so
  plain SCFs and adapter-wrapped calculations dispatch
  under one flight (ARCHITECTURE 9.4).  New adapters
  slot in by implementing the protocol; the dispatch
  core never changes (Principle 8).

#### 6.2.3 Parsl dispatch and complete-and-report

Each unit becomes one Parsl app: a `python_app` that
runs `unit.wingbeat.run(unit, wingbeat_dir)` on a worker.
Kaleidoscope's `parsl_config` (a Parsl `Config`, supplied
by the client/deployment) maps those apps onto SLURM via
a `HighThroughputExecutor` and a SLURM provider, so the
same dispatch code serves a laptop, an interactive node,
and a batch allocation -- only the `Config` changes.

That same `Config` also selects the *cluster topology*,
and both shapes are supported (VISION Goals 4, 7,
ARCHITECTURE 9.4): a single shared allocation whose
workers stream many units (pooled), or one scheduler job
per unit (for large or heterogeneous runs).  These are
`Config` shapes, not changes to the dispatch core.  Who
assembles that `Config` for the producer, and right-sizing
heterogeneous parallel units, are open (6.2.7,
ARCHITECTURE 9.8); the per-unit size such sizing needs is
what section 8 predicts.

The two workload shapes ARCHITECTURE 9.4 calls out are
both expressed in this one model:

- **Embarrassingly parallel sweeps** (thousands of
  independent SCFs): each unit is an independent app
  future; Parsl schedules them across the executor's
  workers.
- **Tightly iterative inner loops** (adaptive
  convergence, future AIMD): the *iteration* lives
  inside the unit's wingbeat (it calls the 6.1 API in a
  loop, or drives ASE's optimizer), so kaleidoscope still
  dispatches one unit; it does not need to model the
  inner loop as a DAG.  If a future client genuinely
  needs cross-unit data flow, Parsl's own futures compose
  -- but that is not required by D13.  The k-point
  convergence climb takes the *other* shape this one
  model allows: a producer-side control loop that sends
  and collects one rung at a time (3.12.5), not a wingbeat
  inner loop.  That is a deliberate choice -- Principle
  12 keeps the energy-reading, next-mesh logic in client
  Python -- and 3.12.5 records the reasoning; both shapes
  are valid, and the dispatch core is neutral between
  them.

**Complete-and-report (Principle 10)** is the dispatch
core's contract.  Kaleidoscope gathers all futures and
**catches exceptions per future** rather than letting
one propagate: a unit whose app raised, or whose Parsl
task was lost cluster-side, is recorded with the
appropriate status (6.2.4) and the flight continues.
No single unit failure aborts the batch.  When all
futures have resolved, kaleidoscope returns a
`FlightReport` (6.2.6); deciding whether the aggregate
is scientifically acceptable is the client's job, never
kaleidoscope's.

**Two phases, separately callable.**  Dispatch is two
steps in sequence: a *send-off* that walks the units,
reports the cache hits, and hands each miss to the
executor -- returning one future per unit -- and a
*collect* that resolves a single future, writes its
terminal status, and builds its report entry.  The
ordinary `dispatch` runs send-off and then collects every
future in unit order, which is the right shape for a
one-shot fan-out.  But a control-loop client cannot use
that shape: the climb (3.12.5) must react to *whichever*
rung lands first, not wait out a whole batch in unit
order, so it needs to send one rung, wait on that rung
alone, and send its successor the instant it lands.  Both
phases are therefore public.  The client sends the units
it has decided so far, then repeatedly collects the next
future to finish and decides what to send next; a future
exposes `done()` -- true once its result is ready -- so
the client can find which of several outstanding rungs has
landed without blocking on a particular one.  `dispatch`
itself is nothing more than send-off followed by
collecting all, in order: the convenience form of the two
public steps, kept identical so every existing one-shot
caller is unchanged.  The domain stays out of the core
throughout (Principles 9, 12): kaleidoscope moves units
and reports outcomes; which rung to send after an energy
lands is the producer's decision alone.

**No error correction in the core -- the custodian
boundary.**  A natural question is whether kaleidoscope
should auto-correct a failing run the way custodian (the
pymatgen job-babysitter for VASP and friends) does:
detect a known error signature, edit the inputs, and
rerun the job in place until it succeeds or hits a retry
ceiling.  It must not.  custodian's value *is* its
embedded domain knowledge -- it knows what a given solver
error means and how to repair it -- which is exactly the
coupling Principle 9 keeps out of the dispatch core (the
same reason `detail` is opaque, 6.2.2).  So a failure
here is terminal (`failed` or `lost`, 6.2.4): kaleidoscope
records it and the client decides acceptance.  Where
custodian-style correction *does* belong is one layer
down -- inside imago.py's own iterate-to-convergence
logic, or inside a smarter wingbeat that loops on the 6.1
API, precisely the "tightly iterative inner loop" shape
above.  In that arrangement custodian's true analog is
the wingbeat plus imago.py's intra-run resume, not
kaleidoscope: the layering mirrors VASP's own
`FireWorks/jobflow -> custodian -> VASP` as `kaleidoscope
-> wingbeat -> imago.py`.  kaleidoscope still dispatches one
unit and sees one outcome; recovery within a run lives
below it, and whole-run reuse on resume lives in the
cache (6.2.5).

#### 6.2.4 Workspace layout (resolves ARCHITECTURE 9.8)

This pins the strawman of ARCHITECTURE 9.6 into a
committed scheme.

```
<root>/
  flight.toml             generated from the Flight
                          (6.2.1): what to run.
  structures/<id>/        acquired/curated inputs.
  wingbeats/<id>[/<calc>]/      one working dir per calc:
      <staged makeinput inputs + run outputs>
      cache_key.toml      identity snapshot (6.2.5).
      result.toml         wingbeat-persisted native result
                            (6.2.6); Imago for ImagoWingbeat.
      status.toml         lifecycle + outcome (below).
  results/                client-aggregated outputs.
  logs/
```

**Stable-id convention.**  `<id>` is the client-supplied
stable per-structure key.  Kaleidoscope requires it to
be filesystem-safe and unique within the flight:
lowercased, restricted to `[a-z0-9_-]`, with any other
character rejected at `Flight` build time (not
silently rewritten -- a surprising rewrite would break
the cache hit-test, 6.2.5).  The producer uses the
curation `reference_id`; an acquisition flight uses the
COD id.  Uniqueness collisions abort the build with the
two offending units named.

**`<calc>` slot format.**  The optional second level
exists only when one structure hosts more than one
calculation (different bases, a property run vs. its SCF,
a sweep over a varied axis).  A unit with no calc tags
runs directly in `wingbeats/<id>/` with no second level.  When
present, every directory component obeys the
`[a-z0-9_-]` rule and must be unique among the calcs
sharing an `id`.  For the legacy single-calc-per-id case,
kaleidoscope derives a default tag from the wingbeat's job
identity (for the Imago wingbeat,
`"<job_name>-<basis_scf>"`, e.g. `"scf-mb"`), and errors
only if that derived tag still collides.

**Sweep flights: one directory level per varied axis.**
A *sweep* flight -- the Principle-12 shape where one
structure hosts a list of units that differ only by one
or more swept option values -- needs more than the
single-string tag above.  The convention is a **directory
tree, one level per varied axis**, in a stable order the
flight-builder helper (6.2.8) declares.  Concrete
examples for a graphite host with successive levels of
sweep complexity:

```
Single-axis (k-density sweep over 3 values):
  wingbeats/graphite/kpt-density-100/
  wingbeats/graphite/kpt-density-150/
  wingbeats/graphite/kpt-density-200/

Two-axis (cell x k-density):
  wingbeats/graphite/cell-2x2x2/kpt-density-100/
  wingbeats/graphite/cell-2x2x2/kpt-density-150/
  wingbeats/graphite/cell-2x2x2/kpt-density-200/
  wingbeats/graphite/cell-3x3x3/kpt-density-100/
  wingbeats/graphite/cell-3x3x3/kpt-density-150/
  wingbeats/graphite/cell-3x3x3/kpt-density-200/

Three-axis (basis x cell x k-density):
  wingbeats/graphite/basis-mb/cell-2x2x2/kpt-density-100/
  ...
```

The flight's `flight.toml` records the axis ordering
and the *fixed* axes (axes that take the same value for
every unit), so harvest and re-judging do not have to
recover them from the path:

```toml
[flight.sweep]
varied_axes = ["basis", "cell", "kpt-density"]
fixed_axes  = { functional = "lda" }
```

**Naming rules** that keep the convention bidirectional
(tags are recoverable into `(axis, value)` pairs by
parsing the path):

1. Every axis name and every value uses only
   `[a-z0-9-]` (lower-case, digits, hyphen); underscores
   are reserved for future use and are not permitted
   inside a level.
2. The first hyphen of a level splits axis from value:
   `kpt-density-200` parses as
   `("kpt-density", "200")`.  Multi-token axis names
   (`kpt-density`) are therefore allowed; multi-token
   values are not.
3. Decimal points in numeric values become `p`:
   `200.5` is recorded as `200p5`, parsed back as the
   real number `200.5`.  Negative numbers prefix `m`:
   `-0.3` becomes `m0p3`.
4. The flight's chosen `varied_axes` order is
   authoritative -- all units in the flight produce
   the same tree shape (no missing levels mid-tree),
   even when one unit happens to share a value with
   another at some level.

**Why a tree rather than one flat tag.**  Two reasons.
First, flat tags balloon: a 4-axis sweep at typical value
widths produces 60-80 character names that wrap and
break tab-completion.  Second, the tree mirrors how
humans actually navigate the data -- `ls
wingbeats/graphite/cell-3x3x3/` shows every k-density value
swept at that cell size, a natural slice.  `find wingbeats
-name 'kpt-density-200' -type d` finds the same
swept-value across cell sizes, the orthogonal slice.
Neither slice is convenient under a flat string.

**Kaleidoscope itself stays domain-agnostic about the
convention.**  Its dispatch core stores whatever string
(or sequence of strings, walked as nested directories) the
client set on the unit and only validates against the
`[a-z0-9_-]` rule.  The convention lives in the
flight-builder helper (6.2.8) so domain knowledge stays
out of kaleidoscope (Principle 9).  A sweep client that
bypasses the helper is responsible for setting the
per-axis directory components per these rules; the
legacy single-string fallback will collide and abort if
it does not, surfacing the mistake loudly.

**`status.toml` schema.**  One file per run directory,
rewritten as the unit moves through its lifecycle:

```
id               = "<id>"
calc             = "<calc>"     # omitted when None
status           = "queued" | "running" | "done"
                   | "failed" | "lost"
detail           = "<wingbeat string>"  # e.g. "converged"
wingbeat         = "imago" | "ase" | ...
submitted_at     = <iso8601>
started_at       = <iso8601>    # omitted until running
finished_at      = <iso8601>    # omitted until terminal
runtime_seconds  = <float>      # omitted until terminal
message          = "<text>"

[record]                        # omitted when unit.record is empty
<name>           = "<value>"    # free-form; see below
```

**The `[record]` table: written down, never interpreted.**  A client
may hang facts on a `CalcUnit`'s `record` mapping that belong to
the run but are not inputs to it -- the engine build identity is
the standing case (6.2.5).  Kaleidoscope copies the mapping
verbatim into this table when the unit launches, and the dispatch
core never interprets, compares, or acts on a single value in it.
It is deliberately not part of the cache key and deliberately not
part of the `FlightReport` (6.2.6): it exists so a curator reading
a directory months later, or a reuse plan naming what produced a
stored result, has something to read.  Because it is written at
launch it stays put across a cache hit, describing the run that
produced the result rather than the flight that reused it.

A *wingbeat* may read a value it recognises and copy it into its
own result file, and the Imago wingbeat does exactly that with the
build identity, echoing it into `result.toml` (6.2.2) so a guidance
entry's provenance can record it (7.8).  This keeps each reader on
the file that is its own: the core writes and reads `status.toml`
and never opens a wingbeat's result, while a domain harvest reads
`result.toml` and never opens the core's lifecycle file.  One fact
therefore lands in two files, which cannot disagree -- both are
copied from the same `record` mapping at launch.

The five `status` values are kaleidoscope-owned and
generic.  `queued` / `running` are lifecycle;
`done` / `failed` are terminal wingbeat outcomes
(`done` iff `WingbeatOutcome.ok`); `lost` is the
kaleidoscope-only category for a Parsl-side
disappearance (worker died, allocation expired) where no
`WingbeatOutcome` ever came back.  Convergence does **not**
appear as a status -- it rides in `detail`, per 6.2.2.

#### 6.2.5 The run-reuse cache

The cache is the general kaleidoscope mechanism of
ARCHITECTURE 9.6, split into mechanism (kaleidoscope) and
policy (client) so generality does not cost correctness.

**Mechanism (kaleidoscope).**  Before launching a unit,
kaleidoscope resolves its
`wingbeat_dir = wingbeats/<id>[/<calc>]/` and performs
the hit-test.  Under the 6.2.4 sweep-tag convention,
`<calc>` may expand to multiple nested directory levels
(`<axis1>-<value1>/<axis2>-<value2>/...`); the hit-test
keys on the full leaf path and otherwise behaves
identically -- the cache mechanism is oblivious to how
deep the per-unit tree is.

1. If the directory exists, holds a `cache_key.toml`
   that matches the unit's *current* key (below), and
   its `status.toml` reads `status = "done"`, the unit
   is a **hit**: skip the launch, and report the existing
   outcome straight from `status.toml` / `result.toml`.
2. Otherwise (no directory, key mismatch, or a
   non-`done` status) it is a **miss**: write a fresh
   `cache_key.toml`, set `status = "queued"`, dispatch,
   and update `status.toml` through the lifecycle.

**A key file that cannot be read is a miss, not an
error.**  The byte-comparison has two sides -- the
freshly built source and the copy staged in the run
directory -- and either may be absent: a prepare
directory reclaimed as scratch, a structure cache that
moved between runs, a run directory left half-written by
a job that died.  Every one of those means the same
thing, that this unit's identity cannot be established,
and the answer is always the same: re-run the unit rather
than trust it.  Raising instead lets a single unreadable
file abort a campaign that has already paid for hours of
converged rungs -- Principle 10 inverted, a per-unit
doubt failing the whole flight.

*A key file must be one every unit has, whatever its job
does.*  A miss is the right answer to a doubt, and the
wrong answer to a file that was never going to be there:
a name absent from every run of a given job kind makes
that kind permanently uncacheable, missing forever at
full price and never saying why.  So the declared names
are drawn from `inputs/`, which makeinput writes for
every unit, and never from a surface whose contents
depend on what the unit's job reads.  TODO D23 records
the case that established this: `kp-scf.dat` was declared
against the run-directory root, which only carries the
names a given job reads, so every fingerprint unit --
which runs no SCF -- missed on every campaign from the
day the name was added.

Resuming a flight is therefore *nothing more than
re-running it*: the hit-test over every unit naturally
skips the completed ones and re-dispatches the rest.

**The key has two parts**, mirroring the producer's
existing `is_cached_v2` (DESIGN 5.7) and generalizing it:

- **Scalar fields** -- written verbatim into `cache_key.toml` as
  TOML and compared field-by-field.  For the producer that is a
  single field, the SCF convergence limit (`converg`).  The swept
  k-density is not among them -- each k-point grid is its own unit
  with its own calc tag and run directory, 6.2.4, so the run-dir
  path already keeps them distinct.
- **Key files** -- each declared in `key_fields` as a path
  *relative to the unit's directory*, and compared by
  **byte-comparison of the freshly built copy against the
  staged copy at that same relative path**:
  `<prepare_dir>/<path>` against `<wingbeat_dir>/<path>`.
  For the producer these are `inputs/structure.dat` and
  `inputs/kp-scf.dat`, both the *resolved* files makeinput
  writes rather than the raw skeleton.  makeinput builds
  into `inputs/` in every directory it builds, so a key
  file named there exists for every unit whatever its job
  does; see "Both sides of the compare are `inputs/`"
  below for why that is the only surface the identity
  compare may use.
  Keying on makeinput's output is deliberate: those files
  bake in the inputs that change the result, so changing
  one misses the cache on its own, with no hand-maintained
  list of "options that matter" to fall out of date.  This
  keeps DESIGN 5.7's "byte-compared file copies, no
  hashing, for debuggability" property: a developer can
  diff the files to see *why* a cache missed, which a hash
  would hide.

  *Why two files, and the claim this corrects.*  This
  design previously named `structure.dat` alone and
  justified it by saying that file bakes in **every** input
  that changes the result.  That claim is false, and the
  exception is not obscure: the k-point integration scheme
  reaches `kp-scf.dat` as `KPOINT_INTG_CODE` (makeinput's
  `scfkpint` -> `kp_intg_code`) and appears nowhere in
  `structure.dat`.  One solid at one mesh under two
  different schemes therefore resolved to the same run
  directory with a matching key and a `done` status -- a
  hit that returns the other scheme's answer.  The same
  file also carries the point operations, so a run that
  suppresses the mesh reduction while keeping its atomic
  symmetry (TODO C136) is likewise indistinguishable from
  one that does not.
  This is a different class of fault from a stale hit and
  is worth separating.  The argument for keeping the engine
  build out of the key (below) rests on a false hit meaning
  *the same physics computed by an older binary*, with
  `--force` as the escape valve.  A false hit across
  integration schemes returns *different physics* under the
  name of the physics that was asked for, and prints
  nothing.  A later reader should not "simplify" this back
  to one key file.
  It has been harmless only because `kpoint_integration`
  has held one value in every run to date, and because LAT
  currently reaches the post-SCF properties but not the SCF
  occupation path (1.6), so the two schemes happen to give
  identical SCF energies.  Both of those are about to stop
  being true.

  *Why a second key FILE rather than a key scalar.*  The
  scheme could equally be added to the scalar list, and
  that is the worse choice.  Scalars are compared as a
  whole table, so introducing a name no stored
  `cache_key.toml` carries invalidates every cached unit in
  every surviving workspace at once -- a mass false miss,
  the failure this section otherwise works to avoid.  A key
  *file* costs nothing **provided the name is one every
  unit has**: the comparison walks the files the unit
  declares and byte-compares each against its staged copy,
  and every run directory stages `inputs/kp-scf.dat`
  whatever its job reads, so a same-scheme re-run still
  hits and only a genuine scheme change misses.  Read that
  proviso as binding rather than incidental -- naming this
  file against the run-directory root instead, where it is
  present only for jobs that run an SCF, is what D23
  records.  It also keeps the diagnosis
  legible -- two `kp-scf.dat` files diff to the single
  `KPOINT_INTG_CODE` line that differs.

**Policy (client).**  The client supplies the key fields
in `CalcUnit.key_fields`; only it knows which inputs
define identity for its calculations.  Kaleidoscope never
guesses -- a too-broad key risks false hits and wrong
science; a too-narrow key risks needless re-runs.  One cache
serves every client: the producer keeps no cache of its own.

**The key asks "is this the same calculation?", not "is this result
still good?" (VISION 16).**  The engine build is deliberately outside
it.  A rebuilt engine does not make a stored potential wrong -- that
potential is a starting point every later SCF re-converges -- and the
two ways of being wrong are not symmetric: a false *hit* has an escape
valve in `--force` (6.2.3), a false *miss* has none, so the hours are
simply spent again.  The build identity is the field that would make
the key miss most often, since ordinary development changes it
constantly and changes the physics almost never.

*What this gives up, deliberately.*  A directory produced by an older
engine is reused without comment, so a flight inherits that engine's
SCF bug for every unit that hits.  Accepted knowingly: the artifact is
a starting guess, not a published number, and a curator who knows a
fix landed re-runs with `--force`.  A later reader should not read the
absent comparison as an oversight and restore it -- least of all as a
whole-repository commit stamp, which is what this design carried
before and is wrong both ways at once: over-sensitive (a
documentation-only commit invalidates every cached run in the
workspace, measured in TODO C132) and under-protective (editing
Fortran and rebuilding without committing leaves the stamp unchanged,
so the stale-engine case it was meant to catch is the one it misses).

*Recorded, not compared; and a printed plan in place of the guard.*
The build identity is written into the `status.toml` the driver
already keeps per unit (6.2.4) -- no new file -- and stays put on a
hit, describing the run that produced the result rather than the
flight that reused it.  The driver then states its plan before
spending anything: how many units it will reuse and how many it will
run.  That count is the decision a curator makes, so it always
prints.  The per-unit lines behind it -- reuse or run for each, and on
a reuse the build and finish time read from that file -- are the
evidence for the count, and they follow the reporting rule of 5.7
(closed as TODO C131): a run says what it achieved, not what it is
doing, so per-item narration waits for `--verbose`.  The climb calls
the driver once per round (3.12.5), which is exactly why: an
unconditional line per unit would refill the screen C131 cleared, on
the same code path.  The plan is also available as a preview that
dispatches nothing, and there the per-unit lines print in full
without asking -- reading them one by one *is* the whole purpose of
a preview.  `--force` keeps its meaning --
driver-side, all-or-nothing, every hit becomes a miss -- and becomes
the ordinary way a re-run is requested.  TODO C84 makes the recorded
identity trustworthy by having the engine stamp its own; that work is
for the record, never for the cache.

**Prepare before the hit-test.**  Because the producer's
key file is makeinput's output, makeinput must run *before*
the hit-test can read it.  Once the structure is
*materialized* -- the network fetch that yields the local
skeleton (5.7) -- the producer runs makeinput itself, in the
driver, as a cheap *prepare* step: distinct from that
materialize step, and named to match the run directory it
builds (`build_run_dir` prepares, `run_prepared` runs,
6.2.2).  The prepare step is cheap because makeinput is pure
Python and, when a solid's species/type assignment needs it,
adds only a fast `imago -loen` run; the expensive
ground-state SCF is never part of it.  It produces
`structure.dat` and the staged inputs; the hit-test then
byte-compares this freshly built `inputs/structure.dat`
against the prior run's copy at the same relative path (the
prepare step must not clobber that reference before the test
-- a PSEUDOCODE detail) and, only on a miss, launches the
expensive imago SCF.

**Both sides of the compare are `inputs/`.**  makeinput
writes its outputs into an `inputs/` subdirectory of
whatever directory it builds, prepare and run alike, so
`inputs/` is the one surface both sides always have and
the identity compare uses it on both:
`<prepare_dir>/<path>` against `<wingbeat_dir>/<path>`,
the same relative path each side.  A *run* directory
additionally carries some of those files flattened at its
own root, because that is where imago reads them when the
unit actually runs -- but which names are flattened
depends on what that unit's job reads, so the root is a
job-dependent subset and is not a surface the identity
compare may use.  A *prepare* directory is never run and
so is never flattened at all.

**And a root copy that exists must agree with its
`inputs/` copy, or the unit misses.**  This is a second,
separate test, and what makes it worth stating is the pair
of cases it tells apart -- both of which one test asking
only "is the root copy there and equal?" would answer with
the same bare "miss":

- *root copy absent* -- this unit's job does not read that
  file.  Nothing is wrong and nothing is stale; the
  identity compare has already been made on `inputs/`.
- *root copy present and disagreeing* -- imago would read
  a file the key does not describe.  Miss, and re-run.

The second case is the failure "What a commit owes a
surviving run directory" (below) exists to prevent, and
that rule does prevent it at the source by clearing the
root copies a commit supersedes.  This test is the
backstop, held here rather than assumed away two
subsections down: it keeps the property that *the key
describes what the engine will read* local to the cache,
where a reader checking the cache can see it, and it holds
for a directory built before the clearing rule existed or
by any future path that writes a root copy without going
through a commit.  It costs one stat and one byte-compare
per key file, only when the root copy is there to check.

Deciding a hit in the driver, from local files, is what
keeps a re-run cheap: a hit never
reaches the scheduler, so the surviving misses are the only
units that cost a calculation.  Seating the key on
makeinput's output also closes a correctness gap a
raw-skeleton key leaves open -- a skeleton plus a
hand-listed set of scalars silently reuses a stale result
when an unlisted option changes (a functional, a potential),
whereas that change is already present in `structure.dat`
and misses the cache on its own.  Running makeinput in the
driver rather than the worker is also what lets a cached
unit skip the scheduler entirely; where the driver itself
runs (a login node or its own batch job) is settled in
6.2.11.

**What a commit owes a surviving run directory.**  The root copies the
compare above declines to read still have to be right, because they are
what the engine reads, and keeping them so is the write side's job.
The flattened root copies are not makeinput's doing: makeinput writes
`inputs/` and
nothing else.  It is `imago.py` that copies each file up to the
run-directory root on the first run, when the root copy is absent, and
that thereafter reads the root copy in preference to the staged one.
The root copy is a cache of `inputs/` -- populated once, invalidated
never.

That precedence is right for a directory a person runs by hand: build
the deck once, edit `kp-scf.dat` in place, re-run, and the edit
survives.  It is wrong the moment a driver commits authoritative staged
inputs into a directory that already holds root copies.  The commit
refreshes `inputs/`, the root copies keep the previous calculation's
contents, and the engine reads the root copies -- so the run executes
the OLD physics while the key file, the run's own `summary`, and the
flight report all describe the new.

This is not the false hit of the two-key-file correction above, and in
one respect it is worse.  The root-copy agreement test above catches
it, so the directory misses correctly, every time, forever: it pays the
full price of a calculation on every re-run and returns the old answer
on every one.  There is no `--force` escape, because `--force` only
turns hits into misses and this is already a miss.  The only recovery
is deleting the run directory -- which means the cache's whole purpose,
deciding reuse from local files, is what a curator must give up to get
a correct answer.

Note which way that dependency runs.  The agreement test makes the
stale directory *fail safe*; the rule below is what stops it arising,
so a curator never meets the unrecoverable miss at all.  Neither
substitutes for the other, and a later reader should not drop the test
on the grounds that the rule prevents the case, nor relax the rule on
the grounds that the test catches it.

**The rule: a commit removes the run directory's root copy of every
name it stages, and lets `imago.py` repopulate it.**  Delete, do not
overwrite.  Overwriting would require the commit step to know which
staged names get flattened and where; that knowledge already lives in
`imago.py`, and a second copy of it would drift.  Deleting needs no
such list: an absent root copy is unambiguous, `imago.py`'s existing
copy-up refills it from the staged file, and there stays exactly one
writer of the root copy and one source for its contents.

Inverting the precedence in `imago.py` -- reading `inputs/` first --
would fix the driver's case by silently discarding a hand edit in every
other case.  The precedence is not the defect.  What was missing is
anyone whose job it is to declare a root copy stale, and the commit is
the only step in a position to know.

The obligation belongs to the commit, not to makeinput.  makeinput
builds into `inputs/` and never owns the root, and on the producer path
it builds a *prepare* directory and never sees the run directory at
all, so no amount of cleaning on its side can reach these files.  A
`--reset` guaranteeing a clean build directory is a reasonable
convenience for a hand-run rebuild and for the wingbeat's
build-in-place branch (6.2.2), but it is a separate matter and must not
be mistaken for this fix.

*What a commit deliberately does not remove.*  The stored potential
(`gs_scfV-*.dat`) is an output of the previous run that the next run
picks up as a starting guess, and it stays.  That follows the position
this section already takes: the key asks whether this is the same
calculation, not whether the result is still good, and a potential is a
starting point every SCF re-converges.  Only the inputs decide what
physics runs.

*The case that established this.*  On 2026-07-31 the si_cmce metal was
re-run with `kpoint_integration` changed from `linear-tetrahedral` to
`gaussian`.  The cache behaved exactly as the two-key-file correction
intends -- all 27 rungs missed on `kp-scf.dat` and re-ran -- and every
root `kp-scf.dat` still read `KPOINT_INTG_CODE 1`.  All 27 energies
reproduced the tetrahedral ladder to the eighth decimal, the SCF output
still printed its tetrahedron population line, and the run's own
`summary` reported `SCF KP Integration = Gaussian`.  Eight minutes of
compute bought a verbatim copy of the answer it was meant to be
compared against.

The boundary with 6.1's checkpointing is clean and worth
restating: `imago.py` resumes *within* a run directory
(skip completed integrals, skip an already-done basis
SCF); kaleidoscope decides whether to *launch* the run
directory at all.  The two never overlap.

#### 6.2.6 Harvest handoff and the flight report

Kaleidoscope returns a `FlightReport`: one entry per
unit, each carrying `id`, `calc`, `status`, `detail`,
`wingbeat_dir`, `runtime_seconds`, and `message` -- exactly
the generic `status.toml` fields, nothing domain-specific.
An optional per-unit `on_outcome` callback (6.2.1) fires
as each unit reaches a terminal state, so a client can
stream-process rather than wait for the whole batch.

**Harvest stays on the client side.**  Kaleidoscope does
not read domain data out of run directories (Principle
9).  The handoff is the run directory itself: the wingbeat
persisted its native result there (`result.toml`), so the
client walks the report and, for each unit it deems
acceptable, opens `wingbeat_dir` and reads what it needs.  For
the producer that means: keep units whose `detail ==
"converged"`, reload the 6.1.2 `ImagoResult` from
`result.toml`, read the named site's type block from the
converged `scfV` at `result.outputs["scfV"]` (a multi-type
file; the site's type number comes from `datSkl.map`),
taking that block's coefficients and alphas together (5.7 /
ARCHITECTURE 9.7).  Non-converged or failed units are simply
skipped -- recorded in the report, never harvested.

This is the precise shape of the C48.3 producer-as-client
relationship: kaleidoscope runs and tracks the batch and
owns the cache; the producer declares the units and the
key, then harvests converged potentials from the run
dirs it is told about.

#### 6.2.7 Open details (for PSEUDOCODE / implementation)

Deferred to the PSEUDOCODE pass for D13 or to C68; none
changes the contracts above.

- **Parsl `Config` specifics.**  Executor type, SLURM
  provider parameters, worker counts, and walltime are
  deployment configuration the client supplies, not
  design.
- **`lost`-unit retry policy.**  Whether a `lost` unit is
  retried automatically (Parsl's own retry, or a
  kaleidoscope re-dispatch on the next flight run via
  its non-`done` status) versus left for the client to
  re-launch.  The cache mechanism already makes a plain
  re-run safe; the open question is only whether to
  retry *eagerly*.
- **`flight.toml` as an authoring surface.**  Whether
  it may be hand-written as the primary input rather than
  always generated from the in-process `Flight`.
- **Concurrency limits for tightly-iterative units.**
  Whether such units need a distinct executor or a
  resource cap so a few long inner loops do not starve a
  parallel sweep sharing the same allocation.
- **`result.toml` for non-Imago wingbeats.**  The Imago
  wingbeat persists an `ImagoResult`; what a future
  non-Imago wingbeat persists (and how a mixed-wingbeat
  client reads it back) is that wingbeat's contract, set
  when the wingbeat is added, not here.
- **Producer dispatch `Config` source -- RESOLVED
  (6.2.11).**  How the producer (and other clients) obtain
  a Parsl `Config` -- a tiered per-site resource-control
  file plus a few per-run CLI choices, assembled by a
  shared generator in kaleidoscope -- is settled in 6.2.11.
  Scheduler dispatch is the default; an in-process local
  run is the explicit opt-out.  The code is TODO C100.
- **Right-sizing heterogeneous parallel units.**  Giving
  each unit a block sized to its own predicted cost
  (per-size executors keyed on a section-8 hint) versus one
  uniform worker slice; deferred until imago is parallel
  and the cost predictor exists.

#### 6.2.8 Flight-builder helper for predict-then-verify

> **Where the producer's k-point convergence lives now.**  The
> producer drives k-point convergence through the adaptive mesh
> climb (DESIGN 3.12, producer flow in 5.7), which predicts in
> density but searches in mesh.  The one-shot density-grid builder
> `build_kpoint_convergence` this subsection designed is retired:
> DESIGN 7.7 documents its split into `predict_kpoint_density`
> (predict only) plus the climb machinery (`build_mesh_unit` per
> rung; `mesh_climb` drives the search).  The verification-grid
> design below is retained as the rationale the climb inherits --
> predict a converged operating point, then verify around it --
> now realized by a mesh search rather than a fixed grid of
> densities.

This subsection designs the **first option-axis builder
helper** living inside `src/scripts/kaleidoscope/`: a
small factory function that turns "a structure plus an
options dict plus a loaded `Dataspace` (DESIGN 7.4)" into
a flight of `CalcUnit`s laid out as a verification grid
around the predictor's predicted operating point.  It is
the corollary-of-Principle-12 builder split mentioned in
the 6.2 intro: option-axis sweeps live here, and this is
the first.  Structure-axis sweeps live elsewhere
(`structure_control`, acquisition) and are not in
kaleidoscope's scope.

**Why this helper is in kaleidoscope's package and not
upstream.**  Domain awareness creeps in along three axes
when building such a flight: which guidance entry to
consume, how wide to make the verification grid, and how
to spell the per-grid-point `<calc>` tag (6.2.4).  The
first two are the historical-guidance DB's contract
(DESIGN 7), the third is kaleidoscope's tag convention
(6.2.4).  Placing the helper inside
`src/scripts/kaleidoscope/` keeps Principle 12 honest:
the dispatch *core* (6.2.1-6.2.7) is dumb and stays so;
the *helper* is a separate optional convenience that
domain-aware clients can call.  A client building
sweeps without using the DB skips the helper and
constructs `CalcUnit`s directly.

**Inputs and outputs.**

```
build_kpoint_convergence(
    structure,                 # StructureControl or skl path
    options,                   # dict of non-swept RUN SETTINGS in
                               #   each tool's own coded vocabulary:
                               #   scf_basis, xccode, scfkpint,
                               #   converg, kpshift, ... held fixed
                               #   across the grid and copied
                               #   verbatim into every CalcUnit.
                               #   Carries NO physics-name keys and
                               #   no bookkeeping -- every key is a
                               #   real tool input, and the wingbeat
                               #   forwards it as-is (DESIGN 6.2.10).
    dataspace,                 # Dataspace loaded by
                               #   guidance_db.load() from
                               #   share/historicalGuidanceDB/
                               #   (DESIGN 7.4)
    system_type,               # "crystalline" / "amorphous"
                               #   / "nanostructure" /
                               #   "molecular" (DESIGN 7.2)
    submodel,                  # dict of the three physics-name
                               #   choices the predictor and the
                               #   PredictionRecord speak: "basis",
                               #   "functional", "kpoint_integration"
                               #   (DESIGN 7.6 step 2).  Kept OUT of
                               #   `options` because makeinput would
                               #   reject these names -- they are a
                               #   prediction input, not a tool
                               #   setting (DESIGN 6.2.10).
    verify       = True,       # False -> trust mode (see
                               #   6.2.1 trust-mode note)
    id           = None,       # flight-level unit id; if
                               #   None, derived from
                               #   structure path
    center       = None,       # curator-pinned k-density
                               #   (6.2.9): bypass the
                               #   predictor and build the
                               #   grid around this value
    root         = "",         # workspace root for the
                               #   returned flight; "" lets a
                               #   multi-structure producer
                               #   supply the shared root when
                               #   merging the per-structure
                               #   flights
)  ->  (Flight, PredictionRecord)
```

**Algorithm sketch.**

```
1.  Compute the query signature for the structure per
    DESIGN 7.4 (if `structure` is a path, load it into a
    StructureControl first):
        query_sig = compute_signature(structure,
                        system_type, dataspace.group_table)
    This builds the 13-d composition vector and (for
    crystalline) the lattice-family one-hot; system_type
    is fixed by argument.

2.  Query the predictor (skipped entirely when `center`
    is given -- the curator override consults no history).
    `predict` is the DESIGN 7.4 free function over the
    loaded Dataspace, not a method on a `db` object; the
    (basis, functional, kpoint_integration) sub-model is
    read from the `submodel` dict (DESIGN 7.6 step 2), NOT
    from `options` -- those physics names never enter the
    tool-facing options (DESIGN 6.2.10):
        result = predict(dataspace, query_sig,
                         submodel["basis"],
                         submodel["functional"],
                         submodel["kpoint_integration"])

    `result` is a PredictionResult (DESIGN 7.4) carrying
    the predicted converged k-density, a combined
    confidence score (stage-1 x stage-2 variance, DESIGN
    7.6), an `is_under_trained` flag, the intermediate
    `predicted_gap` / `predicted_magnetization`, and the
    neighbor entry_ids that produced the prediction.

3.  Decide the verification grid:

      if center is not None:
          # curator override (6.2.9): a tight verify grid
          #   centred on the pinned density, or a single
          #   point when not verifying.  No prediction.
          result      = None
          grid_values = (build_verification_grid(center, 1.0)
                         if verify else [center])
          policy      = "curator_override"

      elif not verify:
          # trust mode: one CalcUnit at the predicted
          #   value, no widening.
          grid_values = [result.predicted_kpoint_density]
          policy      = "trust_no_verify"

      elif result.is_under_trained:
          # no useful prior; wide-grid fallback
          #   (DESIGN 7.9).
          grid_values = default_wide_kpoint_density_grid()
          policy      = "wide_grid_no_prior"

      else:
          # predict-then-verify with variance-aware
          #   widening (DESIGN 7.7).
          grid_values = build_verification_grid(
                          result.predicted_kpoint_density,
                          result.confidence,
                       )
          policy      = "verify_around_prediction"

4.  Round and dedupe the grid to integer k-densities,
    then build one CalcUnit per value.  Each value is an
    integer in BOTH the makeinput `kpd` option and the
    `kpt-density-<int>` tag, so the on-disk tag parses
    back to exactly the swept value (6.2.4
    bidirectionality).  The wide-grid defaults (DESIGN
    7.9) are already integers; build_verification_grid's
    logspace floats round here.  Deduping collapses a
    degenerate grid where rounding merged two close
    logspace points:

      kpd_grid = sorted(set(round(v) for v in grid_values))
      units = []
      for kpd_int in kpd_grid:
          unit_options = dict(options)         # tool-only:
                                               #   forwarded as-is
          unit_options["kpd"] = kpd_int        # makeinput
                                               #   options key
          calc_axes = {"kpt-density": kpd_int} # tag-tree
                                               #   display name
          units.append(CalcUnit(
              id        = id,
              calc      = build_calc_tag(calc_axes),
              structure = structure,
              options   = unit_options,
              wingbeat  = "imago",
              kind      = "convergence",       # default role
              key_fields = standard_key_fields(),
          ))

    `build_calc_tag(calc_axes)` returns a
    `tuple[str, ...]` of `<axis>-<value>` directory
    components -- the on-disk representation of the
    sweep position per CalcUnit.calc (6.2.1) and the
    tag-tree convention (6.2.4).  In v1 it is the single
    `("kpt-density-<int>",)` component: one varied axis.

5.  Assemble the PredictionRecord -- the full field set,
    identical to DESIGN 7.7 step 5 and matching the
    [flight.predictions.<id>] fields the harvest hook
    recovers (DESIGN 7.8).  In predict mode `predict()`
    always returns a PredictionResult (never None, DESIGN
    7.4), so no None-guards are needed; the under-trained
    case is carried by `is_under_trained`:

      prediction_record = PredictionRecord(
          policy                   = policy,
          predicted_kpoint_density =
              result.predicted_kpoint_density,
          confidence               = result.confidence,
          is_under_trained         = result.is_under_trained,
          neighbor_entry_ids       = result.neighbor_entry_ids,
          predicted_gap            = result.predicted_gap,
          predicted_magnetization  =
              result.predicted_magnetization,
          system_type              = system_type,
          feature_vector           = query_sig,
          basis                    = submodel["basis"],
          functional               = submodel["functional"],
          kpoint_integration       =
              submodel["kpoint_integration"],
      )

    The (basis, functional, kpoint_integration) sub-model
    enters the helper as the separate `submodel` dict and is
    recorded ON the per-structure record and ONLY there -- it
    is neither copied into the flight-level `fixed_axes`
    (6.2.9) nor mixed into the tool-facing `options`, so the
    same fact never lives in two places and makeinput never
    sees a physics name it would reject (DESIGN 6.2.10).  This
    is what keeps a combined multi-structure flight whose
    structures use different sub-models harvestable: each
    structure's harvest reads its own sub-model from its own
    record (6.2.9; DESIGN 7.8 step 3f).

    In curator-override mode (`center` given) there is no
    `result`: the record instead documents the pinned value
    -- `predicted_kpoint_density = center`, `confidence =
    1.0`, `is_under_trained = False`, empty
    `neighbor_entry_ids`, and `None` predicted character --
    with `feature_vector = query_sig` and the three sub-model
    axes still recorded.

6.  Record the sweep shape so serialize_flight emits
    the [flight.sweep] block (PSEUDOCODE 13.1) and the
    harvest hook can recover the varied axis without
    re-deriving it from run-dir paths (DESIGN 7.8 step
    3a).  In v1 the single varied axis is k-density;
    `fixed_axes` is empty -- the sub-model the run used is
    carried on the per-structure record (step 5), not
    duplicated here (6.2.9):

      sweep = SweepRecord(
          varied_axes = ("kpt-density",),
          fixed_axes  = {},
      )

    Return (Flight(units=units, sweep=sweep, ...),
    prediction_record).  The helper attaches
    `prediction_record` to the Flight under
    `metadata["predictions"][id]` (a one-entry mapping for
    a single structure; 6.2.9) so the harvest hook (DESIGN
    7.8) recovers it from [flight.predictions.<id>] in
    flight.toml.  A multi-structure producer merges these
    one-entry mappings into a combined flight.
```

**Trust mode and the harvest contract.**  When
`verify=False` the helper builds a single-unit flight.
The producer (or any caller) still harvests the
converged potential or other deliverables from that one
run.  Per 6.2.1, trust-mode harvest does *not* auto-stage
a new guidance entry -- a single converged calculation
is weaker evidence than a converged grid, and the user
asked for trust, not for new evidence.  A curator who
wants the trusted value reinforced into the dataspace
can stage it manually per DESIGN 7.4 / 7.8.

**Cross-references.**  The pieces this helper coordinates:

- DESIGN 7.4 -- signature computation
  (`compute_signature`, returning a `Signature`) and the
  `predict(dataspace, query, basis, functional,
  kpoint_integration)` free function (not a method on a
  `db` object).
- DESIGN 7.6 -- the k-NN regression `predict()` runs over
  the dataspace, with the variance-aware confidence.
- DESIGN 7.7 -- `build_verification_grid` widening
  function.
- DESIGN 7.9 -- `default_wide_kpoint_density_grid` and
  the under-trained / no-prior fallback path.
- 6.2.4 -- the `<calc>` tag convention this helper's
  `build_calc_tag` emits.
- 6.2.6 -- the harvest handoff that consumes the
  `PredictionRecord` to write back into the dataspace.

**Single varied axis in v1.**  `build_kpoint_convergence`
sweeps exactly one axis -- k-density -- so `Flight.sweep`
always has `varied_axes = ("kpt-density",)` and every
`CalcUnit.calc` is a one-element tuple.  This matches
DESIGN 7.2's "exactly one verified target" and DESIGN
7.7's single-axis grid.  An earlier draft carried an
`extra_axes` parameter for additional swept axes, but its
sketch conflated a constant tag-level (dict `.update`)
with a true Cartesian sweep and was dropped; genuine
multi-axis sweeps belong in a future helper (next
paragraph), which can build the nested 6.2.4 tag tree and
a multi-entry `varied_axes`.

**Future option-axis builder helpers** will live alongside
this one (multi-axis sweeps, XANES-target sweeps, basis-
size sweeps), each following the same shape: an upstream
domain-aware library plus a small helper that turns its
query result into a `Flight` with the right tag
convention.  None of those helpers belongs inside the
dispatch core; all of them will share the path conventions
and the `PredictionRecord` mechanism this first helper
establishes.

#### 6.2.9 Multi-structure flights: per-structure prediction and run kinds

6.2.8 builds a flight for *one* structure.  A producer
that runs a *set* of structures -- the initial-potential
builder (5.7), and any later survey over many systems --
wants all their runs dispatched as one flat batch (a
single cluster submission, not one per structure).  Two
small additions to the 6.2.1 model let a single flight
carry many structures without losing per-structure
meaning.

**Per-structure prediction.**  A single structure has one
`PredictionRecord`.  For a set of structures a lone
record on the flight is wrong: each structure has its own
predicted operating point, its own confidence, and --
decisively -- its own `system_type`, which determines how
that structure's feature signature is computed.  A
multi-structure flight therefore carries a *mapping* from
structure id to prediction:

```
metadata["predictions"][<structure_id>] = PredictionRecord
```

The single-structure builder (6.2.8) produces a one-entry
mapping, so the single and multi cases share one shape and
the old singleton key is retired.  The harvester already
groups runs by structure id (6.2.6), so when it processes a
structure's group it looks that id up in the mapping; the
flight-wide trust check of 6.2.6 becomes a per-structure
check.

**Per-structure sub-model.**  The (basis, functional,
kpoint_integration) sub-model travels on the per-structure
`PredictionRecord` for the same reason `system_type` does:
in a combined multi-structure flight these axes can differ
per structure (5.7 allows a manifest to set them per
reference solid), so a single flight-level value is wrong.
The sub-model is recorded on the per-structure record and
*only* there -- it is deliberately not also written into
`sweep.fixed_axes`, even for a single-structure flight where
the three would be genuinely constant.  Storing the same
fact in two places invites both drift (the two copies
disagreeing after a later edit) and reader confusion (which
copy is authoritative?), so the design keeps exactly one
home.  The builder receives that sub-model as a dedicated
`submodel` dict, deliberately separate from the tool-facing
`options` it forwards to makeinput and imago (6.2.8 /
6.2.10): the three physics names are a prediction input, not
a tool setting.  The guidance harvest (7.8 step 3f) reads each
structure's sub-model back from *its own* record; a flight
with no record for a structure does not yield that
structure's sub-model and so cannot be harvested into the
dataspace (7.8).  (Three resolutions were weighed: require
one uniform sub-model per producer run; carry the three axes
per record; or a per-structure `fixed_axes` map.  The
per-record form was chosen -- the record already travels per
id and the prediction *was* made under that sub-model, so it
is the record's natural home, and it keeps the per-solid
override 5.7 grants.  `sweep.fixed_axes` remains a general
SweepRecord field for any future axis a flight truly holds
constant across every unit, but in v1 it has no occupant.)

**The mode rides on the prediction.**  A `PredictionRecord`
already names which grid path produced it in its `policy`
field.  That field *is* the per-structure convergence mode,
and a producer may choose it independently per structure:

- `wide_grid_no_prior` -- a broad sweep, used when guidance
  has no usable prior.
- `verify_around_prediction` -- a narrow sweep centred on
  the predicted value.
- `trust_no_verify` -- a single calculation at the predicted
  value; no sweep.
- `curator_override` -- a single point, or a tight sweep,
  centred on a value the curator pinned by hand (the 5.7
  `kpoint_spec` override; this mode bypasses the predictor).

This is the "broad / narrow / none / override" choice made
once per structure -- the nested-loop shape: one outer list
of structures, each with its own inner convergence strategy.

**Run kind.**  A flight may hold runs that are *not*
convergence-sweep points.  The initial-potential builder,
for one, also dispatches a structure-only
`imago -loen -scf no` run per declared fingerprint (5.7),
which the convergence harvester must not mistake for a grid
point.  So every `CalcUnit` carries a `kind` label -- a
short string the dispatch core stores and round-trips but
never interprets (Principle 9) -- and each harvester
consumes only the kinds it understands:

- `kind = "convergence"` -- a k-density grid point; read by
  the historical-guidance harvester (k-density + gap) and by
  the initial-potential harvester (converged potential).
- `kind = "fingerprint"` -- a structure-only loen run; read
  by the fingerprint harvester.

The default kind is `"convergence"`, so an ordinary
single-purpose sweep needs no annotation; only the extra
runs are tagged.  The convergence harvester selects
`kind == "convergence"` before grouping, which is why a
fingerprint run that shares its structure's id no longer
pollutes the grid.

*Why a label rather than a tag convention.*  One could
instead infer "is this a sweep point?" by parsing each
run's `<calc>` tag against the flight's varied axis (6.2.4)
-- a `kpt-density-200` tag is a sweep point, a
`loen-bispec-...` tag is not.  That works for today's two
kinds but is implicit string-matching; an explicit `kind`
reads plainly to a student and extends to a third kind
without new parsing rules.

**Harvest is general (6.2.6 restated).**  "Harvest" is the
general act of walking a finished flight and pulling out a
*specific* result; the swept k-density is only one target.
The same finished flight feeds several harvesters -- the
historical-guidance harvester (k-density), the
initial-potential harvester (potential), the fingerprint
harvester (descriptor) -- each selecting its `kind` and,
where it needs the prediction, indexing
`metadata["predictions"]` by structure id.  The dispatch
core stays domain-ignorant throughout.

#### 6.2.10 The makeinput/imago option-contract seam

*The collision this resolves.*  A unit carries a single
`options` dict (6.2.1), and the default wingbeat forwarded it
**whole** to both tools: `makeinput.build_run_dir(structure,
options, ...)` to build the input deck, and
`imago.ScriptSettings.from_options(options)` to drive the run.
That only works if every key is meaningful to both readers, and
it is not.  makeinput validates **strictly** -- an unrecognised
dest is a contract fault and raises (`unknown makeinput option:
'basis'`) -- while imago reads **leniently**, taking only the
keys it knows.  The two tools also use **disjoint vocabularies**
for the same physics: the basis choice is an imago run-time
selection (`scf_basis`, coded `fb -> 2` via `BASIS_CODE_MAP`)
and is never a makeinput dest at all, because makeinput writes
*all three* basis sets into `imago.dat` as an overlapping set
and the Fortran run selects one.  So the first imago-only key
the producer emitted (`basis`) aborted every unit before
makeinput could build anything.

*Two root causes, kept separate.*

1. **One options dictionary, two tools.**  The single `options`
   dictionary must serve a strict consumer (makeinput) and a
   lenient one (imago) whose recognised keys do not overlap.
2. **Wrong vocabulary.**  The producer emitted physics names
   (`functional`, `kpoint_integration`, `kpoint_shift`,
   `scf_threshold`) where the tools expect their own dest names
   and coded values (`xccode = 100`, `scfkpint = 1`, ...), and it
   emitted build-identity bookkeeping (`imago_commit`) that is
   not a tool input at all.

*Freedom to refactor the tools.*  `makeinput.py` and `imago.py`
are **ours to change** -- they are not fixed external programs.
Where a change to either makes the split cleaner (a new dest such
as the `-converg` option below, an exported `OPTION_KEYS` set, a
helper that reports which keys a tool recognises, or a tidier
`from_options` contract), prefer changing the tool over
contorting the wingbeat around its present shape.  The decisions
below assume that latitude.

*Decision 1 -- the wingbeat owns the split.*  Partitioning the
options dictionary is the **wingbeat's** responsibility, not
`run_structure`'s.  The wingbeat is the one component that
already runs *both* tools
and is, by construction, the imago-specific adapter (6.2.2);
letting it route options keeps the dispatch core domain-ignorant
(Principle 8) while placing the makeinput/imago knowledge exactly
where the rest of that knowledge already lives.  `run_structure`
is therefore no longer a splitter (this amends 6.1.3 / 6.3.6):
the wingbeat calls `makeinput.build_run_dir` and
`imago.run_prepared` itself, handing each its own separated set
of options.

*Decision 2 -- the producer speaks the tools' vocabulary; the
wingbeat is a pure router.*  `make_producer_options` (5.7)
translates the manifest's human-readable physics
(`functional = "wigner"`, `kpoint_integration =
"linear-tetrahedral"`, `basis = "fb"`) into the **dest-keyed,
coded** options the `from_options` APIs already require
(`xccode = 100`, `scfkpint = 1`, `scf_basis = "fb"`, ...).  The
wingbeat then routes purely by namespace and never value-codes:
the translation lives in the producer, where the physics intent
is known, mirroring how `-xccode` already takes the integer code
defined in `xc_code.dat`.

*The physics names still feed the builder -- through their own
channel.*  Those same three names also drive the flight builder's
predictor and its `PredictionRecord` (6.2.8), which need the human
words (`wigner`), not the codes (`100`).  They reach the builder as
a separate `submodel` dict, never through `options`: the coded
`options` stay strictly tool-facing, and the human sub-model keeps
its own home (6.2.8 / 6.2.9).  Only the basis appears in both
channels -- as `submodel["basis"]` and the run-time `scf_basis` --
because it is genuinely both a prediction input and an imago
setting; `functional` and `kpoint_integration` carry a different
value in each channel (`wigner` vs `100`) and so never collide.

*The routing rule (two buckets).*  The wingbeat splits a
unit's options into:

- **imago run options** -- the fixed, known imago key set
  `{job, edge, scf_basis, pscf_basis, serialxyz, valgrind}`
  (6.1), exported as `imago.OPTION_KEYS` so the wingbeat does not
  hard-code it -> handed to `from_options`.
- **makeinput build options** -- everything else, handed to
  `build_run_dir`.  makeinput **keeps its strict unknown-key
  check**, which now serves as the typo backstop: a key that is
  neither an imago key nor a real makeinput dest still raises, so
  the safety that strictness buys is preserved.

*Bookkeeping is not an option.*  `options` is a dictionary of
*tool inputs*, and a fact that reaches neither tool has no business
in it.  The build identity is the standing example: it is recorded
per run and never compared (6.2.5), so it is neither an imago
setting nor a makeinput dest nor part of the key.  It therefore
rides on the `CalcUnit` itself, in a small free-form `record`
mapping the driver copies verbatim into `status.toml` at launch
and never interprets.  That keeps a third "dropped before
forwarding" bucket out of the routing rule entirely, and it keeps
makeinput's strict check meaningful -- an unrecognised key in
`options` is now always a typo, never a deliberate passenger.

| Producer emits | Bucket | Tool dest / note |
| --- | --- | --- |
| `scf_basis` | imago | basis selection, `fb -> 2` |
| `kpd` | makeinput | k-point density (already correct) |
| `xccode` | makeinput | XC functional code (wigner = 100) |
| `scfkpint` | makeinput | k-point integration (LAT = 1) |
| `kpshift` | makeinput | gamma-centred mesh offset |
| `converg` | makeinput | SCF threshold (new dest; below) |
| build identity | not an option | rides on `unit.record` (above) |

*Designed for setting migration -- not rigid about today's
split.*  The routing is by each tool's recognised-key set, not a
hard-wired ownership map, so a setting can move between the tools
without reworking the seam.  Concretely, the basis lives on the
imago side **today** only because `scf_basis` is in
`imago.OPTION_KEYS` and no makeinput basis dest exists; if a
future makeinput grows a real basis option, dropping `scf_basis`
from `imago.OPTION_KEYS` lets the basis key fall through to
makeinput on its own.  The one case the present "imago set, else
makeinput" rule does not yet cover is a setting that must reach
**both** tools at once; that would be handled by routing on each
tool's explicit recognised-set and forwarding a shared key to
both -- the natural extension point, deferred until a real
both-tools setting exists.

*Decision 3 -- SCF convergence threads like `xccode`.*  makeinput
currently has **no** dest for the SCF convergence limit: it
sources `converg_main` from the rc file and writes it into
`imago.dat`.  To let the producer pin it per reference solid, add
a makeinput option `-converg` (dest `converg`, `type=float`) that
**overrides** the rc `converg_main` when supplied and falls back
to it when absent -- structurally identical to `-xccode`
defaulting to 100.  `make_producer_options` maps the manifest
`scf_threshold` onto `converg`.  The value stays a cache-key
scalar (6.2.5), so `standard_key_fields`' `_KEY_SCALAR_NAMES` is
`("converg",)` -- the cache keys on the dest the run actually
used, and on nothing else scalar.

*Related robustness fix (surfaced by the failed smoke run).*
When every unit fails this way, no `result.toml` is written, and
opening a missing one raised `FileNotFoundError`.  A failed or
result-less unit must be treated as **non-converged** (logged and
skipped, 5.7), never an uncaught crash: the report entry `collect`
builds from `status.toml` (6.2.3) carries the terminal status, and
the climb dispatcher's `next_rung` reads it before ever opening a
`result.toml` -- a unit whose status is not a completed run yields
the FAILED marker, so the material stops non-converged (4e.5) and no
missing `result.toml` is opened.

*Follow-on code (for TODO).*  (a) add makeinput `-converg`;
(b) rewrite `make_producer_options` to emit the dest-keyed, coded
vocabulary above; (c) export `imago.OPTION_KEYS` and retire the
third routing bucket -- delete `kaleidoscope.wingbeats`'
`CACHE_ONLY_KEYS` and the partition branch that reads it, and have
the producer's unit builders set `unit.record` instead, which is
where the build identity now travels.  Those two halves must land
*together*: with the bucket gone, a key that reaches neither tool
raises, so a `make_producer_options` still emitting `imago_commit`
into `options` would abort every unit with `unknown makeinput
option: 'imago_commit'` -- the very failure this seam was written
to fix.  (d) move the
partition into `ImagoWingbeat.run` and retire the single-shared-
options call to `run_structure`; (e) update `_KEY_SCALAR_NAMES`;
(f) the completion gate reads `status.toml` before any `result.toml`
-- `collect` records the terminal status in the report entry and the
climb dispatcher checks it before reading a result -- so a missing one
cannot crash the harvest.

#### 6.2.11 Cluster dispatch configuration

Section 6.2.3 establishes that a flight reaches SLURM purely
through its Parsl `Config`, and that *only the `Config`
changes* between a laptop, an interactive node, and a batch
allocation.  This subsection settles how a client -- the
producer first, other flights later -- actually obtains that
`Config`, so the hardcoded local executor of ARCHITECTURE 9.7
gives way to scheduler dispatch by default, with an in-process
local run kept as a deliberate opt-out.  It resolves the four
questions ARCHITECTURE 9.8 left open.

**The three configuration layers.**  A cluster submission is
assembled from three sources, each owned by whoever knows it
best (ARCHITECTURE 9.4):

1. *Site facts* -- the cluster and account, which do not
   change between runs and differ between users: queue names,
   the account string, cores (and accelerators) per node, and
   the commands a worker runs to bring up the imago
   environment.  These live in a dedicated settings file
   (decision 1).
2. *Per-run choices* -- what this particular flight wants: the
   dispatch shape, which queue, how many nodes, the time
   limit (decision 2).
3. *Per-unit size* -- how much one calculation needs.  For now
   a single uniform value; predicting it per calculation is
   deferred (decision 3).

**Decision 1 -- site facts live in a dedicated, tiered
settings file.**  Following the established `*rc.py`
convention (a module returning a `parameters_and_defaults()`
dictionary), cluster facts live in their own file rather than
mixed into `imagorc`.  The file is *tiered*: a newcomer fills
a small core and everything else takes a built-in default,
while a power user may supply as much detail as they wish.
Only the core is required; any omitted key uses its default.

*Getting-started core (enough to dispatch at all).*

- `partitions` -- the queue(s) available; the first is the
  default queue.
- `worker_init` -- the shell commands a worker runs before
  imago (activate the environment, load modules, set paths),
  so a worker can find imago.
- `account` -- the scheduler account string; required only
  where the cluster demands one, omitted otherwise.

*Performance tuning (optional; improves throughput).*

- `cores_per_node` -- lets the generator pack workers onto a
  node; defaults to one worker per node when absent.
- `workers_per_node` / `cores_per_worker` -- how many
  calculations run at once on a node and how many cores each
  gets (today, with serial imago, one core each and as many
  workers as cores).
- default `nodes`, `walltime`, and `default_topology`, so a
  common run needs no per-run options at all.
- `max_blocks` -- how many allocations the pooled shape may
  grow to when work backs up.
- `memory_per_node` / `memory_per_worker` -- two *distinct*
  memory concepts, deliberately not conflated.  A worker is one
  calculation, so `memory_per_worker` is the memory ONE
  calculation needs: it is the *request*, expressed in
  gigabytes, and it is what the generator spends.
  `memory_per_node` is a node's *physical capacity*, in
  megabytes -- what the hardware has, which is what
  `cluster_probe` can discover from the scheduler.  It is a
  *ceiling*, never a request: the generator does not spend it.
  Keeping them apart resolves a naming trap the single-field
  form fell into -- a field named for node capacity but spent
  as a per-job request, forcing every block to reserve a whole
  node's memory.  The `--mem` the generator writes is
  *derived*, not copied: SLURM's `--mem` is a per-node figure,
  and a node runs as many calculations at once as it packs
  workers, so the per-node request is
  `memory_per_worker x workers_on_the_node` -- one worker's
  worth under the per-job shape (one calculation per node),
  the packed worker count under the pooled shape.  This split
  also opens the forward path: a later automatic
  memory-estimator will predict a calculation's need per
  structure and feed `memory_per_worker` directly, while
  `memory_per_node` stands ready as the capacity ceiling to
  check that estimate against -- a limit, a warning threshold,
  or a packing constraint, as the estimator's design settles.

*Advanced and forward-looking (power users; future imago).*

- `launcher` -- how a single calculation is started across
  cores or ranks; trivial today (serial), the seam for MPI /
  GPU runs later.
- `ranks_per_worker` / `threads_per_rank` -- the hybrid
  parallel split once imago runs in parallel: how many MPI
  ranks one calculation spawns and how many OpenMP threads
  each rank drives.  Their product is the cores the
  calculation occupies, so the two together let a user trade
  message-passing breadth against shared-memory threading to
  match the machine and the problem.
- `binding` -- how ranks and threads are pinned to the
  hardware: to cores, to sockets, or to NUMA memory domains.
  Pinning keeps a rank's threads on cores that share a cache
  and a memory controller, which on a multi-socket node is
  often the difference between scaling and stalling on remote
  memory traffic.  Defaults to the scheduler's own placement
  when absent.
- `omp_places` / `omp_proc_bind` -- the finer OpenMP thread-
  placement controls (spread across sockets versus packed
  onto neighbouring cores), for a user who wants to tune
  thread locality beyond the coarse `binding` choice.
- accelerator facts (`gpus_per_node` and how to request them)
  for the future GPU path.
- per-queue overrides, so a setting may differ by queue.
- named profiles, so a user with access to several clusters
  selects one by name.
- `extra_scheduler_options` -- a raw passthrough of arbitrary
  scheduler directives, and a final escape hatch for settings
  the schema does not name, so a power user is never blocked.

The principle is that the *core is tiny and the rest is
invited*: approachable for someone bringing up their first
cluster, rewarding for someone who wants to tune it.

**How the layers resolve.**  Three of the settings above are
not values but *overlays* on the others, so the order they
apply in is part of the contract.  Most general first:

1. the built-in defaults of `parameters_and_defaults()`;
2. the named **profile**, when `--profile` selects one --
   a user with access to several clusters keeps one file;
3. the **per-queue override** for the queue this run will
   actually use, because a setting may legitimately differ
   by queue (a debug queue with a short walltime cap, a
   large-memory queue with a different per-node capacity);
4. the **per-run command-line flags** (decision 2), which
   are the most specific statement there is.

The queue overlay needs to know the queue, and the queue is
itself a per-run choice (`--partition`, defaulting to the
first entry of `partitions`).  So the queue is resolved
first, *from the profile-overlaid file*, then its override
is applied, and only then do the remaining per-run choices
take their defaults from the site the overlays produced.  An
override for a queue this run does not use is simply not
applied; a file may carry overrides for every queue on the
cluster.

**The loader owns every overlay.**  Reading the settings file
and overlaying it are one operation, not two, and the reader
takes the queue as an argument in order to perform both.  This
is a deliberate constraint rather than a convenience.  Were the
overlay a separate step, every piece of code that reads the
file would have to remember to take it, and a reader that
forgot would get settings that look complete and are quietly
wrong -- the cluster-wide walltime where the queue's cap
belongs, which a scheduler answers with a rejected or truncated
job rather than an error naming the cause.  There is more than
one such reader (the per-unit dispatch, and the driver's own
batch submission), and there will be more as the resource-cost
dataspace of section 8 grows its own.  Making un-overlaid
settings *unobtainable* is what keeps them honest: the mistake
is not merely avoided, it cannot be written.

**Every overlay merges per key, at every layer.**  Most
settings are a single value, and overlaying one simply
replaces it.  But a setting may itself be a *block* of
settings -- `orchestrator` and `md` are the two, holding
respectively the driver's cores, memory and walltime, and the
MD job's ranks, memory, walltime and bring-up -- and there an
overlay names only the keys it means to change.  The others
keep the value the layer beneath them gave.  A curator who
writes, of the debug queue, "the driver needs only two
gigabytes there" must not thereby lose the driver's core
count and time limit: they never mentioned those, and would
receive not an error but two plausible-looking fallbacks in
their place.

This is the same rule already stated for the per-run
`--orchestrator-*` flags below, and it is stated once here
because it governs all three overlays -- profile, queue, and
flags alike.  Whole-block replacement silently discards
facts the curator never meant to touch, and it does so at
whichever layer is allowed to do it; there is no layer at
which that is the desirable reading.

The merge goes exactly one level down.  A block of settings
holds plain values, not further blocks, so a deeper merge
would have nothing to descend into and would only obscure
what a given overlay can reach.

Two guards, in keeping with the strict-contract discipline
the rest of the settings file follows.  A key inside an
override that names no known setting is a configuration
error, not a silent no-op -- it is almost always a typo, and
a silently ignored typo in a resource request is exactly the
failure this file exists to prevent.  And an override may
not set `partitions` or `profiles`: those choose *which*
overlay applies, so letting an overlay rewrite them invites
a rule that refers to itself.

The typo guard descends exactly as far as the merge does --
one level, into a block.  The two must agree, and the reason
is that a block is where an unnoticed typo does its quietest
damage.  A misspelling at the top level leaves a stray
setting that nothing reads, which is bad enough; a
misspelling *inside* a block leaves the real key standing at
its old value beside the stray one, so the run proceeds with
the number the curator meant to change and no sign that
anything was ignored.  A queue override reading `rank` for
`ranks` would let a job run at the site's width while its
author believed they had widened it.  Checking only the
outer keys would leave the merge able to reach a place the
guard cannot see.

**Discovering site facts -- the `cluster_probe.py` tool.**
The settings file itself stays *pure data* -- like every other
`*rc.py`, it is just `parameters_and_defaults()`, with the two
required fields shipped as `None` and a `REQUIRED` comment.
The discovery logic is a *separate* program, `cluster_probe.py`,
because reading the machine is real work (subprocess queries,
parsing) that does not belong in a data file.  Much of the
settings file can be read straight off the machine, which makes
first-cluster bring-up far less daunting: `cluster_probe.py`
queries the scheduler and writes a *starter* copy of the
settings file with everything it can learn already filled in
and a brief plain-language note on every setting.  The tool is
*self-contained*: it carries its own copy of the schema and
only ever *writes* a `clusterrc.py`, never reads one, so it
needs no settings file to exist and does no directory lookup at
all -- a clean split where `cluster_config` reads the file and
`cluster_probe` creates it.  The cost is that the key list
lives in two places; a test keeps the tool's copy identical to
`clusterrc.parameters_and_defaults()` so they cannot drift.

Two honesty rules govern what it writes:

- *Only the scheduler is trusted, never the login node.*
  `sinfo` reports the *compute* nodes' queues, cores, memory,
  and accelerators (generic resources), and `sacctmgr` lists
  the accounts the user may charge.  The tool does **not** read
  the login node's own CPU layout (`lscpu` / `numactl`): that
  would describe the wrong machine, so no login-node fact --
  no socket/NUMA topology -- ever reaches the file.  (The
  `binding` / `omp_*` knobs those facts would inform are the
  deferred parallel seam anyway.)
- *A heterogeneous cluster is not guessed at.*  When the nodes
  disagree on a per-node number -- cores, memory, or GPUs --
  the tool does not silently pick one.  It leaves that setting
  blank, flags it `FILL IN`, and lists the distinct values it
  saw ("nodes vary -- core counts seen: 36, 48, 64, ...") so
  the user chooses deliberately.  Only when every node agrees
  is the value filled in.

What a query *cannot* supply is convention and policy:
`worker_init` -- the module loads and environment setup that
let a worker find imago -- is pure site convention, while the
*correct* `account` to charge and *which* of the listed queues
to prefer are policy.  Those stay blank.  The dividing line is
fact versus policy: the scheduler-known facts are filled (or,
when nodes disagree, offered as options), and the human
choices are left marked.  The tool is best-effort and
scheduler-specific (SLURM today); its output is a draft the
user reviews and edits, never an authority.

*Install relationship.*  The settings file ships as a
*template*: the install places `clusterrc.py` in `$IMAGO_RC`
(non-clobbering, so later edits survive reinstalls) with the
required core left as `None`.  It is deliberately *not* a
working configuration -- there are no universal defaults for a
site's queues, account, or worker bring-up -- so it fails loud
until populated.  Populating it (run `cluster_probe.py`, then
complete `worker_init` and `account`) is therefore a required
setup step for any user running on a cluster, exactly
analogous to `unpackImagoDB.py` for the databases.  A
local-only user never touches it: `--dispatch local`
(decision 2) short-circuits before any settings file is read.
Because `cluster_probe.py` only *writes* a `clusterrc.py` and
never reads one, it has no bootstrap dependency at all: it runs
before any settings file exists and needs neither `$IMAGO_RC`
nor a working-directory copy.

*Resolution precedence (the dispatch read).*  Only
`cluster_config` -- the dispatch read -- locates an existing
`clusterrc.py`, and the search order is deliberate: the working
directory first, then `$IMAGO_RC`.  So the convenient default
is the global copy in `$IMAGO_RC` -- populated once, picked up
by every run -- while a `clusterrc.py` dropped beside a
particular campaign overrides it for that run only, letting a
sweep pin different queues or walltime without disturbing the
global settings.

**Decision 2 -- per-run choices are command-line options,
optionally saved.**  The client exposes options --
`--dispatch local|slurm-pooled|slurm-per-job`, `--partition`,
`--nodes`, `--walltime` -- each defaulting from the site file
(the dispatch shape from `default_topology`, covered just
below), so a fully configured site needs no per-run options at
all.  These four size the *worker* job class.  A client that
submits its driver as a batch job exposes three more --
`--orchestrator-cores`, `--orchestrator-memory`,
`--orchestrator-walltime` -- which size the *orchestrator*
class and default from that block of the site file (the
override rule is set out with the block, below).

The everyday path is a single command (captured in the
`command` log the scripts already keep); for a reproducible
record, the client may also write the resolved *dispatch
choices* -- the shape, queue, node count, time limit, and
the profile that fed them -- as a human-readable file in the
run directory, beside the manifest.  The orchestrator's own
shape is not among them: the driver's `sbatch` script is
itself written to the data root, and it records that shape
exactly, in the form the scheduler received it.

The command-line default is `slurm-pooled`, because the whole
point of this work is that the producer and the database seed
*should* reach the scheduler rather than run on a login node
(ARCHITECTURE 9.7), and the seed and the convergence sweeps
are the many-small-similar-units workload one warm allocation
serves best: its packed workers stream the units with no
per-unit queue wait, and one allocation is simpler to reason
about and to release than many.  `slurm-per-job` stays one
flag -- or one site default -- away for large or heterogeneous
units.  The default is not hardcoded on the flag: `--dispatch`
takes its value from the site file's `default_topology` when
unset, exactly as `--partition`, `--nodes`, and `--walltime`
take theirs, so a site chooses its own default rather than
inheriting one the client baked in.  A run launched where no
site settings file is present is therefore a configuration
error, surfaced up front, rather than a quiet fall-back to a
serial local run -- the cluster behaviour is the default, and
a local run is the deliberate opt-out.  That opt-out is
`--dispatch local`, which needs no site facts and builds no
`Config` at all; the test suite, a
laptop session, and the materialize pre-flight all request it
explicitly, so they neither read a settings file nor touch the
scheduler.  The library entry point mirrors this: its
programmatic default is `local`, so in-process callers (tests
above all) opt *in* to a cluster, never out of one by
accident.

**Decision 3 -- one uniform per-unit size for now;
right-sizing deferred.**  Both cluster shapes (below) give
every calculation the *same* resource slice in this round.
Giving each calculation a slice matched to its own size
(right-sizing) needs both a parallel imago and a predictor of
per-calculation cost, neither of which exists yet, so it is
deferred.  The hook is already named: the per-unit size is
exactly what the resource-and-cost dataspace (section 8,
VISION Goal 6) predicts, and the provisioning consumer that
fills it in is TODO C81 -- a later layer that drops onto this
one without disturbing it.

**The two cluster shapes are two `Config` shapes.**  Both are
expressed entirely in the `Config` the generator builds; the
dispatch core (6.2.3) is unchanged.

- *Pooled* -- one allocation (optionally allowed to grow to
  `max_blocks`) whose workers stream many units.  The
  generator builds a high-throughput executor over a SLURM
  provider sized by the per-run nodes/walltime and the site's
  per-node packing.  Best for many small, similar units --
  the convergence sweeps and the database seed.
- *Per-job* -- one scheduler submission per unit.  The
  generator builds an executor in which each unit maps to its
  own one-worker block, so each calculation queues and runs
  independently.  Best for large or heterogeneous units.

**A block asks for its slice, not for the node.**  Under
either shape, a block's request names exactly what its own
workers need -- their cores and their memory -- and claims the
node no more broadly than that.  This follows from Decision 3
rather than adding to it: if every calculation is given one
uniform slice, then a block holding `w` workers is asking for
`w` slices, and the node's remaining cores are not its to
hold.  They belong to whoever asks next, be that another block
of the same run or another user entirely.

The reasoning bites hardest on the per-job shape, whose whole
promise -- each calculation queues and runs independently --
holds only if the calculations can be *scheduled*
independently.  A one-core calculation that claims a whole
node has not queued independently: it has queued for a node,
and each sibling behind it waits for a node of its own however
few cores it would use.  On a cluster that permits sharing,
one gets in and the rest stay pending while the cores they
asked for sit idle beside the one that ran.

The cores must be named explicitly, and this is the trap worth
recording.  A block that merely declines to claim the node,
without stating a core count, receives the scheduler's default
of a single core -- correct for a one-worker block by accident
and wrong for a packed one, whose workers would then contend
for one core between them.  So the cores are derived exactly
as the memory is, as `cores_per_worker x
workers_on_the_node`, and the two directives are written
together or not at all.

A site whose queue policy genuinely requires whole nodes is
served by `extra_scheduler_options`, which exists for exactly
this kind of local rule and costs the schema no new setting.

**Reclaiming a retired chain's workers (a forward note).**
This extension is designed here but *not built*; it waits on a
parallel imago, and (a)'s job now is only to keep the door to it
open (the three disciplines at the end).  The climb (3.12.5)
retires a chain the moment it converges or hits its ceiling, and
retirement is final -- a retired chain never climbs again.  Under
the pooled shape the run holds one allocation for its whole life,
so as chains retire the block goes on holding every core it asked
for while fewer and fewer of them have work to do.  Handing those
freed cores to the chains still climbing -- so one late,
expensive chain (the seed run's `si_cmce`, 28 rungs against the
others' eight or nine) can finish on several cores instead of
one -- is reclamation.  It is a *distinct* idea from the
right-sizing Decision 3 defers, and a cheaper one: right-sizing
predicts how much a given calculation will need, whereas
reclamation predicts nothing and merely hands out what is
demonstrably idle.

What makes it tractable is the same finality that 3.12.5 relies
on.  The hard half of dynamic resource-sharing is that resources
come *and go*: a general scheme must preempt, rebalance, and
guard against two jobs each holding half of what the other needs.
None of that arises here, because a retired chain never returns,
so the count of idle cores only ever grows over a climb's life.
An allocation of freed cores made now can never need to be taken
back, because nothing will arrive to reclaim it.  No preemption,
no rebalancing, no fairness policy, no deadlock -- the climb's
own shape deletes the whole difficult half of the problem.

The one part that repays a careful reading is the difference
between a *worker* and a *core*, because reclamation turns on it.
A core is one physical processor: real, finite, the thing a
calculation actually consumes.  A Parsl *worker* is a process
that runs one task at a time and then asks for another; it is a
consumer of work, not a slice of hardware.  Parsl's bookkeeping
counts *busy workers* -- worker 7 has a task, so give it no
other; worker 8 is idle, so it may take the next one -- and it
never counts cores at all.  The site's `cores_per_worker` is
arithmetic that *sizes the request* (a block of `w` workers asks
the scheduler for `w x cores_per_worker` cores, exactly as it
asks for `w x memory_per_worker` of memory); it is not a limit
Parsl imposes on a running task, and nothing stops a task from
using more cores than one.  Today every task runs imago on a
single core, so Parsl's count of idle workers *happens* to equal
the count of idle cores.  That equality is a coincidence of the
one-core-per-task regime, not a fact Parsl maintains.

Reclamation ends the coincidence deliberately.  When two chains
have retired, the producer sends a surviving chain's next rung as
one ordinary task whose wingbeat runs imago across, say, three
cores.  Parsl still sees one busy worker and reports the other
two idle -- and that report is *true*: those two worker processes
have no task.  What is no longer true is that idle workers imply
idle cores, because the three cores those workers would have used
are the ones the wide rung is running on.  Nothing has been
deceived; Parsl's number still means exactly what it always meant
(idle worker processes), it has simply stopped doubling as a core
count, which was never what it measured.  This is safe because
Parsl acts on its idle-worker count *only when something submits
a task*, and during a climb the sole submitter is the producer --
the very party that did the core arithmetic and knows those cores
are spoken for.  It will not submit against them.

That safety rests on one condition, and because breaking it fails
silently the condition must be stated: **during a climb the
producer must be the only submitter into its pool.**  If anything
else were to submit, Parsl would place that task on an
"idle" worker whose cores are in use, two imago processes would
contend for the same core, and the symptom would be a slowdown
with no error and no log line naming its cause.  Nothing today
submits into the climb's pool (the pre-flight is a separate
phase), so the condition holds; it is written down because a
future concurrent submitter would break it invisibly.

There is no way to have Parsl keep the core count honestly on our
behalf.  Its one mechanism for tasks of differing size,
`MPIExecutor`, partitions an allocation by *whole nodes* -- it is
built for a task that spans several machines -- whereas our pool
is many one-core workers packed onto a *single* node, a
granularity it cannot express.  So the core count has nowhere to
live but the producer, which is why reclamation is a producer
concern and not a Parsl configuration.

Three disciplines, cost-free today, keep (a) able to grow this
later; skipping them would force a retrofit:

- *A unit may carry a resource request the core round-trips but
  never interprets*, exactly as `kind` (6.2.9) is carried and not
  read.  The width a rung wants rides on the unit; kaleidoscope
  honours it and stays ignorant of what it means (Principle 12).
- *The producer tracks the set of outstanding chains, not a fixed
  worker-per-chain correspondence.*  What to send next, and how
  wide, is computed from what has retired at the moment of
  sending -- so the loop must not bake in "one chain, one
  worker."
- *The width is recorded but is never part of the cache key.*  A
  rung run on one core and the same rung run on three are the
  same physics and must remain one cache entry (6.2.5), or every
  change of width would silently re-run a ladder.  The width is a
  property of a *cost observation* (section 8), so it belongs in
  `status.toml`, never in `cache_key.toml`.

**Decision 4 -- the generator lives in the dispatcher
package.**  The helper that turns (site facts + per-run
choices) into a `Config` belongs in `kaleidoscope`, which
already owns dispatch and is imported by every flight client,
so the producer and future flights share one copy.  It reads
the site settings file and the per-run choices and returns
either a Parsl `Config` (for a cluster shape) or nothing (for
`local`).

**Changing the producer over.**  With the generator in place, the
producer stops forcing a local executor.  For a cluster shape
(now the default) it attaches the generated `Config` to the
flight (`flight.parsl_config`) and lets `dispatch` select the
Parsl path (6.2.3); for the `local` opt-out it attaches
nothing and dispatch runs in process.  This removes
the only reason the seed and database builds run on the login
node, and is the code tracked as TODO C100.  Because the
generator lives in kaleidoscope and `dispatch` makes the
local-versus-cluster choice itself, every client uses this
same change -- no client writes its own executor builder.  The
run-reuse cache-bypass the producer's executor helper used to
carry moves with it: it becomes a `dispatch` argument, since
it governs the cache the driver owns (6.2.5), not where a
unit runs.  The run-reuse cache (6.2.5) is unaffected:
workers execute each unit in its
own run directory on the shared filesystem exactly as the
local executor does, so a cluster run and a local run share
one cache.

**Driver location -- login node or its own batch job.**  The
decisions above settle where each *unit* runs -- the `Config`
and the dispatch shape.  They leave open where the *driver*
runs: the orchestrator process that reads the manifest,
prepares each unit (6.2.5), decides cache hits from local
files, and submits and awaits the rest.  That driver now does
real per-unit work before any SCF -- a makeinput build and,
when a solid's species/type assignment needs it, a fast
`imago -loen` run, once per unit including cache hits.  At seed
scale (a handful of solids) this is negligible and the driver
may run interactively.  At scale it is serial work that would
occupy a login node's terminal for the whole flight, so the
driver may itself be wrapped in a scheduler job.

**A separate orchestrator resource block.**  The driver's own
footprint is sized independently of the per-unit slice
(`memory_per_worker`, decision 1): the two describe different
things -- one calculation versus the orchestrator process --
and conflating them would missize both.  The orchestrator block
follows the dispatch shape.  Under `slurm-per-job` or `pooled`
the driver only prepares units and fans the calculations out to
worker jobs, so it asks for little -- modest memory, a core or
two, and a walltime long enough to outlast the flight it
supervises.  Under `local` the driver runs the SCFs itself, in
process, so its block must be compute-sized.  The block is a
new site/per-run setting alongside the worker sizing, not a
reuse of it.

**The orchestrator shape is overridable per run, key by key.**
The site block is a *default* shape, not a fixed one.  A run
overrides any of its three keys from the command line
(`--orchestrator-cores`, `--orchestrator-memory`,
`--orchestrator-walltime`, decision 2), and that is precisely
what keeps the settings file bounded: a second orchestrator
whose driver needs more memory -- a future builder, or this
producer under `--dispatch local`, where the driver runs the
SCFs itself -- says so for its own run instead of earning a
second block in the file.  Without the override, the file would
have to grow one block per orchestrator, which is the shape
ARCHITECTURE 9.4 rules out.

Three properties of the override, each a deliberate choice:

- **Per key, not whole block.**  Overriding the memory must
  leave the site's cores and walltime standing.  A
  whole-block replacement would silently discard the site
  facts the curator never meant to touch.  This is the
  general overlay rule of decision 1, applied at the last
  layer: the flags merge into the block the profile and
  queue overlays have already built, exactly as those two
  merged into the built-in defaults.
- **Walltime alone keeps a further fallback.**  An unset
  `cores` or `memory` simply goes unrequested and the
  scheduler applies its own default, which is harmless.  An
  unset walltime is not harmless -- a driver job with no time
  limit is a job that can hang a queue -- so when neither the
  flag nor the block names one, the driver's job takes the
  run's resolved `--walltime`.  It always carries a limit.
- **The worker flags do not reach the driver.**  `--walltime`
  and `--nodes` size the *worker* job class (decision 2).  A
  curator shortening `--walltime` to clear a short queue is
  speaking about the calculations, not about the process that
  submits them, and the orchestrator's own walltime continues
  to govern the driver's job unless
  `--orchestrator-walltime` says otherwise.

**A third job class -- the md block.**  `condense.py` writes a
LAMMPS submission file, and that job is neither a worker nor an
orchestrator.  A worker is one serial calculation; an
orchestrator is one driver process; an md job is many MPI ranks
of an external program filling a single node.  Sizing it as
either would missize it -- the same argument that separated the
orchestrator from the worker above, applied once more
(ARCHITECTURE 9.4).  The block holds `ranks`, `walltime`,
`memory`, and its own bring-up lines.

**Cluster facts come from the site file; everything else stays
with the script.**  `condense.py` already reads `condenserc.py`
for the chemistry and structure settings that are its own
business.  What it does not own is the cluster: the queue, the
account, and a node's core count are the same facts a flight
uses, and writing them down a second time in `condenserc.py`
would leave one truth in two files with nothing to mark which
had aged.  So the script reads cluster facts from `clusterrc.py`
and keeps the rest where they are.  Reading two settings files
is the price of not duplicating a fact, and it is the cheaper of
the two prices on offer.

**The md bring-up is its own, not `worker_init`.**  `worker_init`
starts imago; an external MD program does not need imago, and
does need something else entirely.  A separate `init` in the md
block is what lets a site record where it installed LAMMPS
without that location appearing anywhere in `src/` -- which is
the whole point of routing this through configuration instead of
a hard-coded template.

**Ranks derive from the node, and a missing node fact yields one
rank.**  A condensation run fills a node, so `ranks` defaults to
`cores_per_node` rather than to a number written into the
source.  A hard-coded count is exactly what left the previous
template asking for 125 tasks while also naming a partition
whose nodes have 48, so that the alternative it appeared to
offer could never have run.  `cores_per_node` is itself
optional, and where a site has not recorded it the generator
asks for a single rank and says so in a comment in the file it
writes.  A one-rank MD job is visibly wrong to whoever opens it,
which is the point: guessing a plausible core count instead
would produce a job that runs, and runs wrong, while never
announcing that the site was never configured.

**Sizing falls back; the bring-up does not.**  Where a site has
left the md block's *sizing* keys unset the generator supplies
them -- ranks from the node, walltime and memory from their
defaults -- the same latitude the orchestrator block has in
being allowed to stand empty.  The bring-up is different in
kind.  A rank count guessed wrong yields a job that runs
badly; an absent `init` yields a job that cannot start at all,
because nothing has put the MD program on the path.  So `init`
is marked in the settings file the way `worker_init` is:
shipped blank, flagged as required, and left blank by
`cluster_probe`, which cannot supply it -- where LAMMPS was
installed is site convention, exactly the fact-versus-policy
line the probe already draws.  The generator refuses to write
a submission file without it.

The *enforcement point* is chosen as deliberately as the
marking.  `init` is not added to the required core that the
loader checks on every read.  That check guards a dispatch,
and a site that flies calculations but never condenses has no
reason to record where an MD program lives; making it globally
required would refuse a flight over a setting the flight never
reads.  The requirement is real and it belongs to the md
generator, which is the only party that needs it.

**An unconfigured site is refused, not papered over.**
`condense.py` reads the site file through the same loader every
flight uses and inherits its verdict: where the required core
is unfilled, the run stops with an error naming what is
missing.  This is the rule decision 2 already states for a
dispatch, holding here for the same reason.  A file that looks
like a working submission while naming the wrong queue, the
wrong account, and no bring-up is worse than no file, because
it fails at the scheduler -- later, on another machine, with
nothing in the failure pointing back at the settings that were
never filled.  The case this protects is not the machine with
no configuration at all, which nobody condenses on; it is the
half-configured one, where enough is present that the omission
does not announce itself.

That refusal belongs at the *start* of the run, not at the
moment the file is written.  `condense.py` computes bonds,
angles and the whole LAMMPS input before it has any use for a
queue name, and an error surfacing after all of that would
throw the run away and be met a second time on the rerun.  So
the site file is read where the script reads its other
settings, before any work, and the submission file is later
written from settings already in hand.  This is the same "up
front, never deep in a run" discipline the loader already
states for a flight; it is worth restating here only because
the natural place to *use* the site file is nowhere near the
right place to *read* it.

Timing decides *which* file is read, besides deciding when the
refusal lands.  The loader searches the current directory
before `$IMAGO_RC`, so that a run may carry a settings file of
its own; and by the point the submission file is written the
script has moved into `lammps/`.  A read there would resolve
against a directory the run created rather than the one the
user launched from -- silently picking up a stray file, or
silently missing the intended one.  Reading at settings time
fixes the answer while the working directory still means what
the user meant by it.

**One thread per rank, written by the generator.**  The job
fills a node with MPI ranks, so each rank must confine itself
to one core.  A threaded BLAS left to its own devices assumes
it has the machine and starts a thread per core *in every
rank*, which on a forty-core node is sixteen hundred threads
contending for forty cores.  The generator therefore writes
`OMP_NUM_THREADS=1` itself rather than leaving it to the site's
`init`: it follows from the shape of the job -- ranks sized to
fill the node -- and that shape is the generator's own
arithmetic, not a site convention.  It is written after the
bring-up, so a module that sets a thread count of its own
cannot overwrite it.

**Every generated submission file brings up its own
environment.**  A bring-up is a sequence of shell commands, and
on a module-based cluster those are `module` commands, which
exist only in a shell whose profile has run.  A submission file
that omits the login shell appears to work whenever it happens
to be submitted from a shell already set up, and fails from
`cron`, from a workflow driver, or under `sbatch
--export=NONE` -- inheriting what it ought to establish.  Worse,
it fails *there* rather than where it was written, so the report
is of an unattended job dying rather than of a generator that
assumed its caller's environment.

Both generators therefore open with a login shell: the md job
because its `init` is module commands, and the orchestrator job
because `worker_init` may equally be.  The rule belongs to the
act of generating a submission file, not to either job class, so
neither generator is left as the one a site discovers the rule
through.

**Materialize on the login node, then submit.**  The one step
that needs the network is the structure fetch -- the
materialize pre-flight (`--materialize-only`, the `local`
opt-out of decision 2).  It runs on the login node first,
pinning and caching every structure, because compute nodes may
have no internet.  Only then is the driver's batch job
submitted.  Inside that job the prepare step (6.2.5) consumes
the already-fetched skeletons and touches no network, and
dispatch runs the calculations.  So the single
network-dependent step is isolated to the login node, up
front, and nothing downstream depends on a compute node
reaching COD.

**One flag, and one deferral.**  Which shape the driver's job
uses is the same per-run `--dispatch` choice (decision 2):
`local` inside the orchestrator job for seed scale now,
`slurm-per-job` or `pooled` later -- a flag, not a rewrite.
One thing is deferred with it.  At seed scale the driver
prepares every unit serially inside its job, which is exactly
what keeps a cache hit off the scheduler (6.2.5).  When that
serial prepare becomes the bottleneck, prepare-and-hit-test
can itself move onto dispatched worker units -- at the cost
that a hit then occupies a cheap worker slot rather than being
decided driver-local.  That transition is a later refinement,
turned on when the serial cost bites.

#### 6.2.12 Scratch reclamation

A finished run directory holds two very different kinds of
file, and the difference is almost entirely one of size.
The **kept** tier lives in the run directory itself: the
staged inputs, `result.toml`, `status.toml`,
`cache_key.toml`, the SCF potential, the descriptor, the
log.  The **scratch** tier lives behind the `intermediate`
symlink, which `imago.py` creates pointing at a temporary
area (typically a fast local or scratch filesystem), and
holds the engine's working files -- above all the HDF5 that
carries the wavefunctions.

Measured on a seed-scale producer run, the split is stark:
scratch was **99.7%** of the bytes (3.17 GB of 3.2 GB, about
25 MB per calculation, essentially all HDF5), against 222 KB
of kept files per calculation.  Everything else in scratch
was numbered copies of files the run had already written
home.  So reclaiming scratch is nearly all of the available
saving, and reclaiming anything else is nearly none of it.

**Why this is safe, and how we know.**  Two consumers could
in principle be broken by removing scratch, and neither is:

- The **harvest** reads only paths recorded in
  `result.toml`'s `outputs` table, and every one of them
  -- the energy, the iteration count, the SCF potential,
  the expanded structure, the `datSkl.map`, the log --
  resolves inside the run directory.  None points through
  `intermediate`.
- The **run-reuse cache** (6.2.5) decides a hit from
  `status.toml` and `cache_key.toml` alone.  Scratch is
  never consulted, so a pruned run still hits.

Those two facts are what make reclamation a tidying
operation rather than a destructive one: a pruned run
directory answers every question the producer ever asks of
it.  They are properties of the current design, though, not
laws, so a client that starts reading through
`intermediate` must say so (below) rather than discover the
loss later.

**Mechanism and policy, split as the cache is.**  Section
6.2.5 already divides the cache into a kaleidoscope
mechanism and a client-supplied key, because only the
client knows what defines identity for its calculations.
Reclamation divides the same way, and for the same reason:

- *Mechanism (the reclamation tool).*  Walk a root,
  resolve each run directory's `intermediate` target, and
  remove what a policy marks reclaimable -- with the
  dry-run preview, the reporting, and the refusal rules
  below.
- *Policy (client).*  Decide **when** a unit's scratch is
  reclaimable.  Only the client knows whether a finished
  run is finished *with*: the k-point producer is done
  with a rung the moment its `result.toml` lands, while a
  client that post-processes wavefunctions into a density
  of states needs the HDF5 until that step has run.

**Where the mechanism lives, and why not in the
dispatcher.**  The *split* above is the cache's, but the
*placement* is not: the cache mechanism sits inside
kaleidoscope, while this one sits beside it, in a tool of
its own (layer (c) below).  The reason is what reclamation
has to know.  It reads imago's own names -- the
`intermediate` link, the lock file, the wording of the
completion line -- out of `imago.py` rather than
re-spelling them, so that a rename in the engine cannot
silently desync the recognizer.  That is engine knowledge,
and within kaleidoscope there is exactly one place engine
knowledge is allowed: the **wingbeat**, the pluggable piece
that knows how to run one unit.  The dispatch core beneath
it -- the driver and the data model -- names no imago file
at all, which is what would let a different engine be flown
by the same dispatcher.  This is Principle 9 exactly:
domain-specific machinery lives at the adapter layer, and
the flight layer is ordinary scientific Python.  Principle
12 says the matching thing from the other side, naming the
wingbeat as where per-unit domain iteration belongs.
Reclamation is not
a way of running a unit, so it is not a wingbeat; and
putting it in the core would place engine knowledge in the
one layer that must not carry any.  It therefore lives
outside kaleidoscope altogether, and each of the three
layers below reaches it from the *client* side rather than
from within a flight.

A policy is handed `(run_dir, target)`: the run directory
and its already-resolved scratch.  The second argument is
part of the contract because one of the built-in policies
must look *inside* the scratch to decide (the job tree's,
below), and a client is no worse off for being given the
path it is judging.

Each kind of root supplies the default that fits it, and a
client overrides by passing its own.  For a workspace the
default is the conservative one -- reclaimable when the
status is `done` and `result.toml` exists.

**Two kinds of root, one contract per call.**  A
kaleidoscope workspace is not the only place imago scratch
accumulates, and it is not even the common one.  Every
ordinary `imago.py` run -- a student's job directory, a
hand-driven convergence test -- plants the same
`intermediate` link and leaves the same tens of megabytes
behind it.  Those runs sit in a **job tree**: run
directories with no `wingbeats/` and no flight above them.
Reclamation recognizes both roots and decides which it has
by looking rather than by being told:

- a root holding `wingbeats/` is a **workspace**;
- a root with no `wingbeats/` but `intermediate` links
  somewhere below it is a **job tree**;
- a root with neither is not a reclamation target, and is
  refused.

One call handles exactly one kind.  Mixing them in a single
report would gather two different safety contracts under
one set of totals, and the contracts really do differ --
which is the next point.

**The job tree has no `status.toml`, but it is not without
evidence.**  Workspace reclamation rests on a unit
declaring itself `done`.  A hand run writes no
`status.toml` and no `result.toml`, so that authority is
absent -- but the run is not silent, because a hand run is
not a separate code path.  The CLI is a thin wrapper over
the same callable core (6.1.3), so an `imago.py` run in a
job directory leaves the same two traces every kaleidoscope
unit does:

- **`imagoLock`, inside the scratch itself.**  It is
  created before any work begins and removed in the
  cleanup that always runs, so its presence means the run
  either owns the directory now or died without releasing
  it.  Either way the scratch is not ours to take.
- **`Program Sequence Complete.` in `runtime`.**  The
  driver appends it as it closes the log.

Both traces are `imago.py`'s to define, not the cleanup
tool's to guess.  The lock name, the log name, and the
completion marker are read from `imago.py` directly rather
than re-spelled as literals here, so that a change to how
the engine names or marks a run cannot silently desync the
recognizer.  The failure would be in the safe direction --
"absence of evidence is refusal" means an unrecognized
marker refuses rather than deletes -- but it would be
silent, and a silent refusal of everything is harder to
diagnose than a mismatch that fails loudly.

Requiring *both* is what makes the test sound.  The lock is
released a moment before the marker is written, so each
covers the other's blind spot, and together they close the
window a bare staleness test would leave open: a run
starting *now* takes the lock before it writes anything, so
a reclamation racing it is refused rather than mis-timed.
No arbitrary "old enough" threshold is needed anywhere.

Two details decide whether this works in practice.  The
`runtime` log is opened in **append** mode, so a directory
run four times holds four markers, and only the *last*
non-blank line describes the *current* state -- a run
interrupted after three good ones ends in something else
entirely.  The test therefore reads the tail, never
`grep`.  And the marker is written from a `finally`, so it
records that the driver reached cleanup, not that the
calculation succeeded: a run that failed but exited tidily
is reclaimable here, where the workspace contract would
have preserved it for the curator by finding no
`result.toml`.  That is the one place the two contracts
genuinely disagree, and a job tree has no success signal
with which to close the gap.  `--older-than` is the lever
for anyone who wants recent failures kept.

That age filter must be measured correctly, which is
subtler than it looks: the mtime of the scratch *directory*
moves only when entries are added to or removed from it, so
a job that has spent a week writing into an
already-created HDF5 still presents a week-old directory.
Age is therefore the newest mtime anywhere in the tree,
never the top directory's.  The size walk already visits
every file, so one pass yields both.

**Five refusals that hold for every root.**  Reclamation
deletes, so it is defined by what it will not do.  These
five apply to a workspace and a job tree alike; the job
tree then adds two more of its own:

1. **Never delete the run directory, only scratch.**  The
   kept tier is the record of the calculation and is two
   orders of magnitude smaller than the saving; there is
   no case for touching it.  The `intermediate` symlink
   itself is left in place, dangling, so the run directory
   still shows where its scratch was.
2. **Never delete a unit that is not finished.**  A
   `running` or `queued` status means the engine may still
   be writing, and a missing `result.toml` means the run
   produced nothing to keep.  Both are skipped and
   reported rather than pruned, since the second is
   usually the state a curator most wants to investigate.
3. **Never follow a link out of the scratch area.**  The
   target is resolved and checked to be under the scratch
   root before anything is removed.  An `intermediate`
   pointing somewhere unexpected -- the `FIXME` rename
   `imago.py` performs when a link is stale is one way it
   happens -- is reported and skipped.  A cleanup tool
   that can be redirected by a symlink is a hazard, not a
   convenience.
4. **Never descend through a symlink while walking.**  The
   third refusal guards where scratch may *be*; this one
   guards what is even *considered*.  `intermediate` is
   itself a symlink into the scratch area, so a walk that
   followed links would leave the root, find whatever lives
   over there, and could plan a removal outside the tree
   it was pointed at.  The walk therefore skips every
   symlinked subdirectory, which keeps the set of
   candidate run directories inside the root it was given
   by construction rather than by later checking.
5. **Never remove a scratch tree holding another run's
   scratch.**  Scratch mirrors the run directory's path,
   so a run nested inside another -- a `debug/`
   subdirectory of a job that was itself run in place --
   has its scratch nested inside the outer run's scratch
   too.  Removing the outer tree would take the inner one
   with it: at best double-counting the saving, at worst
   deleting the working files of a run the refusals above
   had just declined to touch, which turns refusal 2 into
   a formality.  The outer tree is skipped while any other
   run's scratch lies within it, and named as such.  Once
   the inner ones are reclaimed a second pass takes the
   outer, so nothing is lost -- only deferred.

   The containment test looks at *every* run directory the
   walk found, not merely the ones this call selected: a
   filter that excluded the inner run would otherwise let
   the outer removal delete it as collateral, which is
   precisely the case the refusal exists to stop.

**Two further refusals, for the job-tree contract.**

6. **Never reclaim a run that has not declared it
   finished.**  A hand run is reclaimable only when its
   scratch holds no `imagoLock` *and* the last non-blank
   line of its `runtime` log is the completion marker.
   Absence of evidence is refusal, not permission: a run
   directory with no `runtime` at all, or one whose log
   ends mid-stream, is reported and left alone.
7. **Never descend into a workspace from a job tree.**  A
   job tree may hold a workspace far below it, and walking
   in would apply the presumption-based contract to units
   that have a provable one.  Such directories are not
   descended, and every one is named in the report, so the
   bytes deliberately left behind stay visible instead of
   vanishing from the accounting.

**Dry run is the default.**  The tool previews what it
would remove, with sizes, and removes nothing until asked.
An operation whose entire purpose is deletion should make
the destructive path the one the user typed deliberately.

**Three layers, one code path.**  The logic lives in one
place and is reached three ways:

- (c) A **standalone tool**, the home of the logic:
  selective over a workspace by unit and by calc kind,
  over a job tree by path, and over either by age, with
  the preview and the refusals above.  This is what a user
  runs to reclaim a finished campaign, or to sweep up
  after a season of hand runs.
- (a) The producer's **`--clean-after`**, which calls that
  same logic once its harvest completes.  It supplies only
  the workspace and leans on the default policy to gate
  each unit, rather than naming the harvested units: the
  policy is the same "spent once `done` with a
  `result.toml`" rule the producer would apply anyway, so
  the two cannot diverge, and a unit the flight left
  unfinished is preserved rather than swept up as a side
  effect of the harvest.
- (b) **Prune-as-you-go**, which applies the client policy
  as the flight advances rather than after it ends,
  discarding a unit's superseded scratch while later units
  are still running.  It matters for a campaign large
  enough to exhaust scratch mid-flight.  There is no fixed
  headroom to name here: `$IMAGO_TEMP` has no per-user
  quota, so the ceiling is whatever a *shared* filesystem
  happens to have free -- unguaranteed, not the user's to
  count on, and reduced by everyone else's jobs at the
  same time.  That is a stronger reason to prune in flight
  than a generous private allowance would be, not a weaker
  one: a campaign cannot know in advance how much room it
  will actually get.  As in (a), the producer names no
  policy of its own and leans on the default -- a unit is
  spent once it is `done` with a `result.toml` -- so the
  in-flight prune, the post-harvest sweep, and the
  standalone tool all judge a workspace by one rule.  A
  client whose units are not spent that early passes its
  own, which is the whole reason the policy is an argument.
  (b) is also the layer that deletes while runs are still
  in progress, so it is the one whose policy and refusals
  have to be right; the two paragraphs after this list say
  where it hooks and what the mechanism had to grow to
  support it.

**How prune-as-you-go hooks.**  It adds nothing to
kaleidoscope, for the reason given above, and it does not
need to: a flight already tells its client when a unit
reaches a terminal state.  The `on_outcome` callback fires
once per unit, in landing order, carrying that unit's run
directory -- exactly the moment, and exactly the fact, a
prune needs.  So (b) is a *client* wiring: the client hands
the flight one callback that reports the landing and prunes
it.  The alternative considered was a reclamation policy
recorded on the flight itself, which the dispatcher would
act on; it was rejected because it buys uniformity across
clients at the price of teaching the dispatch core an
engine's file names.  A client wanting both progress
reporting and pruning composes them in its own callback,
which is a few lines it can read, rather than inheriting a
deletion it never asked for.

**A prune that fails is contained, but never hidden.**  This
is Principle 10 applied to housekeeping -- one failure never
fails the flight, and failures are recorded and surfaced --
so what follows is that principle worked out, not a new
rule invented here.  The
callback runs inside the dispatcher's collect step, which
does not guard it, so an exception escaping it would
propagate out through the climb and abandon the campaign --
trading hours of cluster time for a piece of housekeeping.
Nothing may escape, therefore.  But *contained* must not
become *ignored*, and the distinction that makes both
possible is between a refusal and a failure.  A refusal --
unfinished, nested, too recent -- is the mechanism working
exactly as the rules above say it should, and is reported
at the same quiet level the standalone tool uses.  A
failure -- a removal attempted and declined by the
filesystem, or the mechanism itself raising -- means an
assumption behind this whole section has stopped holding:
a permission changed, a mount went read-only, something is
holding files open.  The next campaign will meet it too.
So a failure is printed when it happens *and* carried to
the end of the run and printed again, because a single line
an hour deep in a log is lost, and this is the class of
thing that must not be.  It is deliberately not made fatal:
the databases and the run log are complete and correct, and
a housekeeping fault is a thing to go and look at, not a
reason to call a good campaign failed.

**What the mechanism had to grow.**  The standalone tool
judges a whole tree at once; (b) judges one directory at a
time.  That is a separation, not a second rule: the
per-directory decision -- resolve the link, refuse what
must be refused, apply the policy, measure -- becomes a
function that the whole-tree planner calls in a loop and
the in-flight prune calls once.  A campaign pruned in
flight and one swept afterwards are then governed by the
same code, which is the same argument that made layer (a)
call the tool's planner instead of reimplementing it.

The one refusal that does not survive being narrowed to a
single directory is the nested-scratch rule, which by its
nature compares one tree against others.  The
per-directory function therefore receives the set to
compare against as an argument: the whole-tree planner
passes every run its walk found, and the in-flight caller
passes the flight's other units.  A unit that has not
started yet contributes nothing to that set, which is
right -- it has no scratch to lose, and it will create its
own when it runs, inside a parent that by then is gone.

### 6.3 makeinput callable build API

This subsection designs the makeinput counterpart of the
imago callable API of 6.1.  It is the makeinput-side twin
of D11/C63: it turns `makeinput.py` from a command-line-
and-cwd-bound script into a script that *also* exposes a
callable "build a run directory" function, with the CLI
reduced to a thin wrapper over it.  It exists to resolve
the one piece 6.1.3 deferred: `run_structure`'s
*structure-and-options* mode promises to "drive makeinput
to build `run_dir`, then call `run_prepared`," but there
is no in-process makeinput entry point to drive.  6.3
supplies it.  The work was folded into C68 as item (a) of
the kaleidoscope prong; this design rung was missing, so
it is captured here before the code lands.

#### 6.3.1 Why a callable build API

The default kaleidoscope wingbeat (`ImagoWingbeat`, 6.2.2)
calls `imago.run_structure(structure, options, run_dir)`
on any unit whose run directory is not already prepared.
`run_structure` in turn must build that directory from a
structure plus a set of makeinput options.  Today the only
way to run makeinput is to invoke its `main()`, which
parses `sys.argv`, operates on the current working
directory, and -- on a missing `$IMAGO_RC` -- calls
`sys.exit`.  None of those is acceptable inside a long-
lived kaleidoscope Parsl worker driving thousands of
builds: a worker has no per-build `argv`, must place
inputs in an arbitrary `run_dir` rather than its own cwd,
and must never be terminated by a `sys.exit` raised deep
in a build (the same hazard 6.1.2 designs against on the
run side).  6.3 therefore mirrors 6.1's split: a callable
core that takes its inputs as arguments and reports
contract faults by raising, wrapped by a thin CLI that is
the only layer touching `argv` or exiting the process.

The boundary on error handling matches 6.1.2 exactly.
*Build-level* faults that are normal outcomes of real
input (a malformed `imago.skl`, an element with no basis
in the database) keep their existing diagnostic *meaning*.
*Contract* faults -- the environment is not configured
(`$IMAGO_RC`/`$IMAGO_DATA` unset), the named structure
file does not exist, the target `run_dir` cannot be
created, or a build is unsupported (an unimplemented
grouping op, a `-pot` override naming an absent database
entry) -- raise a `MakeinputError`.  `MakeinputError` is
the makeinput analog of `ImagoError`: a programmer- or
environment-level fault that no per-unit retry can fix, so
it propagates out of the worker's wingbeat where the
flight records the unit `failed` and continues
(Principle 10).

The one behavior that must change regardless of a fault's
category is **process exit**.  The historical script
signals several faults with `sys.exit`, but `sys.exit`
raises `SystemExit`, which derives from `BaseException`,
*not* `Exception` -- so it slips past the dispatch core's
per-future `except Exception` (6.2.3) and would abort the
flight (the in-process executor) or kill the worker (a
Parsl executor) rather than failing one unit.  Therefore
**no `sys.exit` may remain on the build path**: each
becomes a raised exception -- a `MakeinputError` when no
retry can help, a natural `Exception` otherwise (which the
dispatcher likewise catches) -- and the thin CLI wrapper
is the only layer that exits the process, printing the
message and exiting non-zero.  6.3.5 records the specific
conversions and the audit that the rule is complete.

#### 6.3.2 The build entry point

The API offers a single entry point, because makeinput has
only one granularity of input (a structure plus options);
there is no prepared-directory analog to short-circuit.

- **`build_run_dir(structure, options, run_dir, *,
  settings=None) -> str`** -- given a structure and a set
  of makeinput options, stage the structure into `run_dir`
  and run the full makeinput workflow there, producing the
  staged Imago inputs (`imago.dat`, `structure.dat`,
  `scfV.dat`, the kp files, the `inputs/` tree) that
  `imago.run_prepared` then consumes.  Returns the
  `run_dir` it built (absolute), so a caller can chain
  directly into `run_prepared`.  `structure` is, at this
  design stage, a path to an `imago.skl` -- the same
  commitment 6.1.3 makes for `run_structure`; whether it
  may also be an in-memory `StructureControl` is deferred
  to the ASE-free factory of D12/C64 and is *not* fixed
  here (6.3.7).  When `settings` is omitted it is built
  from `options` via the resolution path of 6.3.3; a
  caller that already holds a reconciled settings object
  (the CLI does) may pass it to avoid rebuilding.

The **CLI wrapper** (`main()`) becomes the outermost layer
and the only one that touches `sys.argv` or exits.  Its
three responsibilities mirror 6.1.3's CLI split:

1. Parse `sys.argv` into makeinput options (the existing
   argparse surface and `reconcile` logic, unchanged in
   meaning).
2. Build the run directory.  A bare `makeinput ...`
   operates on the current working directory, which holds
   `imago.skl` -- today's only behavior -- so the CLI
   calls `build_run_dir` with `run_dir = os.getcwd()` and
   `structure = "imago.skl"`.  No new CLI surface is
   required by this design; the CLI keeps doing exactly
   what it does today, now through the API.
3. On a raised `MakeinputError`, print the message and
   exit non-zero, preserving today's diagnostics.

#### 6.3.3 ScriptSettings split (mirrors C63)

makeinput's `ScriptSettings.__init__` today performs four
steps in the constructor: load the rc defaults, parse
`sys.argv` (`parse_command_line`), `reconcile` the parsed
namespace against the defaults, and `record_clp` (append
the literal `sys.argv` to a `command` file).  Only the
middle two carry meaning the API needs; the first and last
are CLI couplings.  The refactor splits construction the
same way C63 split imago's:

- The constructor loads the rc defaults only and leaves the
  job-type/edge/basis fields unset.
- **`from_command_line()`** -- the CLI path: parse `argv`
  into an `args` namespace, then `reconcile(args)`.
- **`from_options(options)`** -- the API path: turn the
  `options` mapping into the same kind of `args` namespace
  the argparse parser would have produced (every key absent
  from `options` takes its argparse default), then
  `reconcile(args)`.  The keys of `options` are exactly the
  argparse `dest` names (`job`, `edge`, `basis`, `scfkp`,
  `pscfkp`, `reduce`, `target`, `block`, `xanes`, `potdb`,
  `basisdb`, ...), so a client and the CLI describe a run
  identically -- one through a dict, one through flags.

`reconcile` is unchanged: it already takes an `args`
namespace and contains all the option-resolution logic, so
both paths share it verbatim.  This is the single change
that lets settings stop being constructed from `argv`
unconditionally.

#### 6.3.4 cwd discipline and structure staging

makeinput is thoroughly current-working-directory bound:
`initialize_cell` reads the skeleton from the relative name
`imago.skl`, and every file and directory it writes
(`inputs/`, `.inputTemp/`, `imago.dat`, `structure.dat`,
the kp files) is a cwd-relative name.  The API keeps that
internal convention -- rewriting hundreds of relative
paths to be `run_dir`-relative would be invasive and
error-prone -- and instead adopts the **same cwd
discipline 6.1.4 designs for the run core**: `build_run_dir`
treats the cwd as a resource to acquire and release.  It

1. resolves `run_dir` to an absolute path and creates it;
2. stages the structure into it as `run_dir/imago.skl`
   (a copy when `structure` is some other path; a no-op
   when it already *is* `run_dir/imago.skl`);
3. `os.chdir(run_dir)`, runs the workflow, and
4. **restores the original cwd in a `finally`, including on
   failure.**

Without step 4 a single failed build would strand a
flight worker in a stale directory and corrupt every
subsequent build's relative-path resolution -- the same
reentrancy hazard 6.1.4 calls the most important
correctness difference between the one-shot CLI and the
reentrant API.  Because the lock-free makeinput build and
the locked imago run each acquire and release the cwd
around their own scope, the two compose cleanly when
`run_structure` calls them in sequence.

The workflow itself -- the body of today's `main()`:
`setup_environment` -> `initialize_cell` -> `assign_group`
(species, then types) -> optional XANES/EMU passes ->
`print_imago` -> `print_summary` -- is factored out of
`main()` into a callable `build_inputs(settings, sc)` that
both `build_run_dir` and the CLI invoke, so the build
sequence has exactly one definition and the CLI and API
cannot drift apart.  The progress `print`s that `main()`
interleaves are retained (they are harmless and useful in
both modes); only `argv`/exit handling moves to the CLI.

#### 6.3.5 Call provenance and worker-safe errors

Two CLI couplings are retired from the API path:

- **`record_clp`** appends the literal `sys.argv` to a
  `command` file.  In API mode there is no meaningful
  `argv`, so -- exactly as 6.1.3 resolves for imago's
  `recordCLP` -- this becomes CLI-only.  The build records
  the equivalent provenance (the resolved options and the
  `run_dir`) or skips the `command` file entirely; the
  precise choice is an implementation detail with no
  bearing on the produced inputs (6.3.7).
- **Every `sys.exit` on the build path** becomes a raised
  exception, because `SystemExit` is not caught by the
  dispatcher's `except Exception` (6.3.1, 6.2.3).  The three
  in `makeinput.py` are converted to `MakeinputError`: the
  `_load_rc` missing-`$IMAGO_RC` check, the unsupported
  reduce-grouping op, and the `-pot` override naming an
  absent database entry.  The CLI wrapper catches
  `MakeinputError` and exits non-zero, so the only `sys.exit`
  that remains is its own `sys.exit(main())`.  The in-process
  modules the build reaches (`structure_control`,
  `initial_potential_db`, `element_data`) were audited and
  contain no `sys.exit`; the subprocess execs it spawns
  (`makeKPoints`, `contract`) are exempt, because a child
  process's exit cannot kill the parent worker.

#### 6.3.6 Relationship to run_structure (closes 6.1.3)

With 6.3 in place, `imago.run_structure` is finally
implementable as 6.1.3 always intended: stage nothing
itself, delegate the build to makeinput, then run the
prepared directory.

```
function run_structure(structure, options, run_dir,
                       settings=None):
    import makeinput
    makeinput.build_run_dir(structure, options, run_dir)
    return run_prepared(run_dir, settings=settings)
```

The import is local, so `imago.py` keeps importing without
makeinput's environment loaded -- the same lazy-import
courtesy `ImagoWingbeat` already extends to `imago` (6.2.2).
This is the seam that lets a kaleidoscope flight hand a
bare `imago.skl` plus options to the default wingbeat and
have the run directory built and run in one worker call,
which is exactly the dependency the C48.3 potential-DB
producer is waiting on.

Note (6.2.10): the `options` reaching `build_run_dir` here are
**makeinput-only**.  The default wingbeat partitions a unit's
options upstream and forwards the imago-side keys to
`run_prepared` separately, so `run_structure`'s combined-options
form above is the convenience path for a *direct* caller that
already holds makeinput-only options -- the wingbeat does not
use it to split.

#### 6.3.7 Open details (for PSEUDOCODE / implementation)

Deferred to the PSEUDOCODE pass (§14) or to implementation;
none changes the contracts above.

- **`options` dict shape.**  The mapping is keyed by the
  argparse `dest` names, but the exact normalization of
  multi-valued options (`reduce`/`target`/`block`/`xanes`,
  which the CLI accepts as repeatable token lists) into
  dict values must be pinned in pseudocode so a client and
  the CLI produce identical settings.
- **`structure` type.**  Whether `structure` may be an
  in-memory `StructureControl` in addition to an skl path
  depends on the ASE-free factory of D12/C64; 6.3 commits
  only to the skl-path form, matching 6.1.3.
- **Sequential loen for bispectrum.**  A bispectrum
  fingerprint comes from the sequential loen flow (5.10)
  that `makegroups.py` orchestrates as ordinary dispatched
  units -- makeinput no longer re-invokes itself, so there
  is no nested-subprocess concern to handle here.
- **Call-provenance recording.**  What replaces
  `record_clp` in API mode (record the resolved options, or
  skip the `command` file) is an implementation detail.

---

## 7. Historical Guidance Dataspace

### 7.1 Overview and Motivation

This section pins down the schema, data structures, and
algorithms for the historical-guidance dataspace introduced
in VISION Goal 5 and architected in ARCHITECTURE section 10.
The dataspace records, for each converged calculation imago
has run, a feature vector describing the system's chemistry
and structure, the electronic-structure character that
resulted (band gap, spin polarization), and the convergence
settings that worked (initially: k-point density).  A small
two-stage k-nearest-neighbor predictor learns from this
dataspace: given a new system's feature vector, it predicts
the converged operating point and an uncertainty.  New
calculations then run a verification grid around the
prediction whose width tracks the predictor's uncertainty.

**The motivating workflow.**  Today, converging a new
system means deciding on a set of candidate k-point
densities (say 5-7 values), running them all (with all
other knobs fixed at sensible guesses), inspecting the
resulting energy-vs-density curve, and picking the
cheapest density at which the energy has stopped moving.
If the curve has not converged at the top of the range,
the user extends the range and re-runs.  This is correct,
but wasteful: most systems within a chemical family
converge at similar densities, and a researcher who has
worked with the family for years carries that knowledge
in their head.  When a new student takes over -- or when
an automated pipeline (like the C48.3 initial-potential-
database build) tries to converge many systems unattended
-- the embodied knowledge is lost and the wasteful full
scan returns.

**Why a dataspace and a predictor, not a categorical
lookup.**  The first instinct is a categorical database
keyed on the system's elements, stoichiometry, or some
discrete classification (insulator / semiconductor /
metal).  We considered and rejected that shape (DESIGN
7.10).  The driver of converged k-density is *electronic-
structure character* -- gap width, spin polarization,
Fermi-surface complexity -- which is a continuous
quantity, not a category, and which depends on chemistry
in ways smooth enough that a regression / nearest-
neighbor predictor can learn it.  Binning gap into
discrete classes loses signal at the boundaries (a 0.1 eV
narrow-gap semiconductor binned as "metal" predicts wrong
densities); k-NN over a continuous feature space does not.

**The predict-then-verify workflow.**  Given a new
structure:

  1.  Compute the structure's feature vector: composition
      (atom-fraction weighted across 13 element groups),
      lattice family (one-hot, crystalline only), and the
      4-way system_type (crystalline / amorphous /
      nanostructure / molecular) declared by the user.
  2.  Query the predictor for the predicted converged
      k-density and an uncertainty measure (variance of
      the k nearest neighbors at each stage).
  3.  Build a verification grid around the predicted
      point.  Width scales inversely with the predictor's
      confidence: high confidence -> tight 3-point grid;
      low confidence -> wider 5-7 point grid; an
      under-trained predictor -> wide-grid fallback
      (7.9).
  4.  Dispatch the grid through kaleidoscope (DESIGN 6.2)
      using the flight-builder helper (DESIGN 6.2.8).
  5.  Harvest the converged grid point and the measured
      electronic-structure quantities back into the
      dataspace through staging + curator promotion
      (7.8).

**Why a separate artifact from the initial-potential
database (DESIGN 5).**  Both honor Principle 11
(experience as a curated artifact), and both share the
library/producer/consumer discipline.  But they store
different *kinds* of experience -- DESIGN 5 stores
numerical potential coefficients per element; DESIGN 7
stores convergence-settings advice plus the electronic-
structure character that produced it -- with different
lifetimes (DESIGN 5 grows entry-by-entry under
deliberate curation; DESIGN 7 accumulates from every
successful flight) and different consumers (DESIGN 5
feeds `makeinput.py`; DESIGN 7 feeds kaleidoscope's
flight builder).  The two artifacts share only the
curation discipline, not their contents.  Considered and
rejected: cross-referencing them via a `pot_label`
parameter (closed by decision, 7.10).

**Why a separate artifact from kaleidoscope itself.**
Kaleidoscope (DESIGN 6.2) is the dispatch layer that
runs flights; it is domain-agnostic.  The guidance
dataspace is domain-aware (it understands element
groups, lattice families, electronic-structure
characters).  Putting the dataspace inside kaleidoscope
would violate Principle 9 (kaleidoscope stays dumb) and
would couple two artifacts with very different rates of
change.  The clean separation: kaleidoscope dispatches;
the dataspace + predictor advise; the client glues them
together via the flight-builder helper (DESIGN 6.2.8).

**Why this accelerates the initial-potential-database
build (Goal 3).**  The C48.3 producer is itself a
kaleidoscope client.  It must converge SCF calculations
on many reference solids.  Without the dataspace, every
reference solid requires a from-scratch convergence
study, multiplying the cost of populating the potential
DB.  With a seeded dataspace, every reference solid in a
chemistry family the predictor has trained on inherits a
predicted operating point and needs only a small
verification grid.  The accelerator compounds: every
reference solid the producer converges contributes back
into the dataspace, sharpening the predictor for the
next.

### 7.2 TOML Schema (version 1)

The dataspace is a directory of TOML files, one per
converged calculation, partitioned by `system_type`
(ARCH 10.1).  The top-level marker file `SCHEMA_VERSION`
records the current schema version; readers refuse files
whose `schema_version` field disagrees with the marker.

**Per-entry top-level keys (required):**

  Field           Type    Description
  --------------------------------------------------------
  schema_version  int     Currently 1.  Must equal the
                          top-level marker file's
                          contents.
  entry_id        string  Unique within entries/.
                          Conventionally the slug used
                          in the filename (7.5 emitter
                          contract), e.g.
                          `"crystalline-a1b2c3"`.
  generated_at    string  ISO-8601 UTC timestamp of the
                          flight that produced this
                          entry.
  source          string  Either `"flight"` (the entry
                          came from an automated harvest)
                          or `"manual"` (a curator wrote
                          it by hand).

**Signature block, under `[entry.signature]` (required):**

This is the predictor's feature input.

  Field                Type    Description
  --------------------------------------------------------
  system_type          string  One of `"crystalline"`,
                               `"amorphous"`,
                               `"nanostructure"`,
                               `"molecular"`.  Hard
                               partition: the predictor
                               uses a separate sub-model
                               per system_type.
  composition_vector   inline  13-key inline table.  Each
                       table   key is one of the element-
                               group names listed below;
                               each value is the
                               atom-fraction of that
                               group in the system, in
                               [0.0, 1.0].  The 13
                               values sum to 1.0
                               (rule 4).
  lattice_family       string  REQUIRED iff
                               `system_type ==
                               "crystalline"`.  One of
                               `"cubic"`, `"hex"`,
                               `"tet"`, `"ortho"`,
                               `"mono"`, `"tri"`.  Forbidden
                               otherwise.

The 13 element-group names that key
`composition_vector` (lower-case, underscore-separated):

```
alkali           Li, Na, K, Rb, Cs, Fr
alkali_earth     Be, Mg, Ca, Sr, Ba, Ra
halide           F, Cl, Br, I, At
chalcogen        O, S, Se, Te, Po
pnictogen        N, P, As, Sb, Bi
group_iv         C, Si, Ge, Sn, Pb
group_iii        B, Al, Ga, In, Tl
transition_metal Sc..Zn, Y..Cd, Hf..Hg  (lumped 3d/4d/5d)
lanthanide       La..Lu
actinide         Ac..Lr
metalloid        Si and B already in group_iv/iii;
                 metalloid covers the diagonal -- Ge, As,
                 Sb, Te (the canonical metalloids).
                 Si and B are NOT double-counted here;
                 see 7.4 for the assignment rules.
noble_gas        He, Ne, Ar, Kr, Xe, Rn
hydrogen         H  (its own bucket per discussion)
```

The exact element-to-group assignment table lives in
`share/historicalGuidanceDB/elemental_groups.toml` (a checked-
in data file, not code -- Principle 11).  7.4 describes
its layout and how the library consumes it; 7.10 records
the open ambiguity around metalloids and how it is
resolved.

**Measured-quantities block, under `[entry.measured]`
(required):**

These are the values harvested from the converged
calculation -- both the target the predictor learns to
produce (`kpoint_density`) and the intermediate
electronic-structure quantities that drive it (`gap_ev`,
`spin_polarization`).

  Field              Type    Description
  --------------------------------------------------------
  gap_ev             real    Band gap in electron volts.
                             >= 0.  0.0 indicates a metal
                             (no gap).
  gap_kind           string  One of `"direct"`,
                             `"indirect"`, `"none"`.
                             `"none"` if and only if
                             `gap_ev == 0.0` (a metal).
  spin_polarization  real    Dimensionless fractional
                             polarization at the Fermi
                             level (relevant for
                             metals); 0.0 for closed-
                             shell non-magnetic systems.
  total_magnetization  real  Total magnetic moment in
                             Bohr magnetons per formula
                             unit.  Signed real (negative
                             for the opposite spin
                             convention).  0.0 for
                             non-magnetic systems.
  kpoint_density     real    The converged k-point
                             density that produced this
                             entry.  Units: k-points per
                             Bohr^-3 of reciprocal-cell
                             volume (matching the
                             DESIGN 3 convention).
                             This is the predictor's
                             target.

**Context block, under `[entry.context]` (required):**

  Field                       Type    Description
  --------------------------------------------------------
  basis                       string  `"mb"`, `"fb"`, or
                                      `"eb"`.
  functional                  string  The DFT functional
                                      under which the
                                      calculation was
                                      converged (e.g.
                                      `"lda"`, `"gga-
                                      pbe"`).  Free-form
                                      string; the
                                      predictor groups by
                                      this value into
                                      sub-models.
  kpoint_integration          string  The Brillouin-zone integration method
                                      (e.g. `"tetrahedral"`,
                                      `"gaussian-0.1"`).  Part of the
                                      predictor sub-model key with basis and
                                      functional, because the gap and the
                                      converged k-density depend on it.
  scf_threshold               real    The SCF threshold
                                      used (e.g.
                                      `1.0e-6`).
  cell_atom_count             int     Number of atoms in
                                      the unit cell that
                                      was converged.
  cell_volume_per_formula_unit  real  Cell volume divided
                                      by formula units
                                      per cell, in Bohr^3.

**Verification block, under `[entry.verification]`
(required for `source = "flight"`, optional for
`source = "manual"`):**

Records the verification grid that produced the converged
kpoint_density.  Distinct from the parameter-blocks shape
the earlier draft used (one verification per parameter);
v1 has exactly one verified target (kpoint_density), so
the verification block sits at the entry level.

  Field                        Type    Description
  --------------------------------------------------------
  grid_values                  array   The full list of
                                       k-density values
                                       swept by the
                                       verification grid,
                                       sorted ascending.
  grid_energies                array   The total energy
                                       (Hartree) at each
                                       grid point, parallel
                                       to grid_values (same
                                       length, same order).
                                       Recorded so the
                                       curator's auto-promote
                                       rule (7.8) can judge
                                       the converged region's
                                       flatness from the
                                       staging file alone,
                                       without re-reading the
                                       flight workspace.
  converged_at                 real    The value at which
                                       the convergence
                                       metric was first
                                       satisfied.  Must
                                       equal
                                       `measured.kpoint_density`.
  converged_mesh               ints    The resolved axial
                                       counts [n_a,n_b,n_c] of
                                       the converged rung (3.7
                                       / 3.12.4), kept beside
                                       converged_at so the
                                       exact calculation is
                                       auditable where a
                                       density round-trips
                                       only up to rounding.
                                       Optional: absent on
                                       curator-authored or
                                       pre-mesh entries.
  gap_spread                   real    How far `measured.gap_ev`
                                       still moves with the mesh
                                       AT the converged rung: the
                                       largest relative change
                                       between that rung's gap
                                       and its neighbours two
                                       ladder positions away, as
                                       a fraction.  Small means
                                       the gap has settled; large
                                       means the recorded gap is
                                       an accident of where the
                                       ENERGY converged.  See
                                       below.  Optional: absent
                                       when the ladder is too
                                       short to measure it, or on
                                       curator-authored entries.
  metric                       string  Currently
                                       `"total_energy"`.
                                       Reserved:
                                       `"forces"`,
                                       `"density_change"`.
  metric_threshold             real    The threshold the
                                       metric had to cross
                                       to count as
                                       converged.
  predictor_confidence         real    The confidence
                                       score in [0.0, 1.0]
                                       the predictor
                                       returned at the
                                       time this flight
                                       was launched, or
                                       0.0 if launched
                                       without a
                                       prediction (e.g.,
                                       seed flight).
                                       Records the
                                       strength of the
                                       prior that
                                       produced this
                                       verification.
  predictor_neighbor_ids       array   List of entry_id
                                       strings of the
                                       k nearest neighbors
                                       the predictor used,
                                       empty if no
                                       prediction was
                                       made.

**`gap_spread`: how far the recorded gap can be trusted.**

`measured.gap_ev` is read off whichever rung converged the
ENERGY.  Nothing has ever argued that a mesh chosen to flatten
the energy also flattens the gap, and measurement says it often
does not.  Since `gap_ev` is a predictor key -- stage 1's
regression target and confidence source, stage 2's distance
metric (7.6) -- and since nothing downstream re-converges it the
way a consuming SCF re-converges a potential, an entry that
cannot say how settled its gap was is asking to be trusted
further than it has earned.  `gap_spread` is that statement.

**It is a measurement, not a verdict, and deliberately so.**  A
stored boolean would freeze a tolerance we have not chosen, and
entries written under one tolerance would then silently disagree
with entries written under another -- the exact failure
`metric_threshold` exists to prevent for the energy.  A raw
fraction leaves the choice to the consumer and keeps every entry
comparable.  Nothing in the producer or the harvest acts on this
value today: it is recorded so that the decision of what to do
about an unsettled gap can be taken later, from data, rather
than guessed at now.

**Measured against neighbours TWO ladder positions away, not
one.**  This is forced by the data.  A k-point ladder carries a
strong parity sawtooth in the gap: on diamond silicon, adjacent
rungs disagree by 19% (`[11,11,11]` reads 0.9572 eV against
`[12,12,12]`'s 0.8046) even where the gap is settled to about 1%
within a single parity family.  Odd and even meshes sample the
zone differently near the band edges and converge to the same
limit at different rates.  Comparing a rung to its immediate
neighbours would therefore report every ladder as unsettled and
discriminate nothing.

**Relative, not absolute**, for a reason the same data supplies.
Near the top of its ladder si_ia-3's gap moves by 0.010-0.014 eV
per two rungs, which is *smaller* in absolute terms than diamond
silicon's mid-ladder movement -- yet si_ia-3's gap is collapsing
toward zero while silicon's has settled.  As fractions the two
separate cleanly, about 20% against about 1%.

Absent when the ladder holds no rung two positions either side
of the converged one, which a short climb or a metal
short-circuit can leave.  Absent is "not measured", never
"settled".

**Provenance block, under `[entry.provenance]`
(required):**

  Field             Type    Description
  --------------------------------------------------------
  flight_id         string  The kaleidoscope flight
                            identifier that produced
                            this entry.  For `source =
                            "manual"`, the curator
                            records a free-form tag.
  source_structure  string  The structure that the
                            flight converged.  Free-
                            form: a COD id, a Materials
                            Project id, or a relative
                            path under `share/skl/`.
  imago_commit      string  Git SHA of imago at the time
                            of the flight run.
  curator           string  For `source = "manual"`:
                            the curator's name or
                            handle.  For `source =
                            "flight"`: the name of
                            the harvest script.

**The top-level `SCHEMA_VERSION` marker file format.**
The marker is a single line containing a bare decimal
integer followed by a newline (e.g., `1\n`).  No TOML,
no key, no surrounding whitespace.  Simplest possible
form so the reader does not need a TOML parser just to
decide whether to refuse a file.  Day-1 contents:

```
1
```

**Validation rules** (enforced at load time):

1. `schema_version` must equal 1, and must agree with the
   top-level `SCHEMA_VERSION` marker file (parsed as a
   bare decimal integer per the format above).
2. `entry_id` must be unique across all entry files in
   the entries directory.  Collisions are a hard error
   with both filenames listed.
3. `system_type` must be one of the four valid values
   (`"crystalline"`, `"amorphous"`, `"nanostructure"`,
   `"molecular"`).  An entry's file must live under the
   matching `entries/<system_type>/` subdirectory; a
   mismatch is a hard error.
4. `composition_vector` must have exactly the 13 keys
   listed above.  Each value must be in [0.0, 1.0].
   The sum of all 13 values must be `1.0 +/- 1e-6`.
5. If `system_type == "crystalline"`,  `lattice_family`
   must be present and must equal one of `"cubic"`,
   `"hex"`, `"tet"`, `"ortho"`, `"mono"`, `"tri"`.  If
   `system_type != "crystalline"`, `lattice_family`
   must be absent (or empty string).
6. `gap_ev` must be `>= 0`.  `gap_kind` must be one of
   `"direct"`, `"indirect"`, `"none"`.  `gap_kind ==
   "none"` if and only if `gap_ev == 0.0` (a metal).
7. `kpoint_density` must be `> 0`.
8. `basis` must equal one of `"mb"`, `"fb"`, `"eb"`.
   `functional` and `kpoint_integration` must be non-empty.
9. `cell_atom_count` must be `> 0`;
   `cell_volume_per_formula_unit` must be `> 0`.
10. If `[entry.verification]` is present (required for
    `source = "flight"`): `converged_at` must equal
    `measured.kpoint_density`; `grid_values` must be
    sorted ascending and contain `converged_at`;
    `grid_energies`, when present, must have the same
    length as `grid_values` (the two arrays are parallel
    and share an order);
    `metric` must appear in the metric registry
    (initially `{"total_energy"}`);
    `predictor_confidence` must be in [0.0, 1.0];
    `predictor_neighbor_ids` must be a (possibly empty)
    array of strings that refer to existing entry_ids
    (referential integrity is not enforced at load --
    a neighbor entry may have been promoted out --
    but the field is recorded for forensics).
11. `source` must equal `"flight"` or `"manual"`.  For
    `source = "flight"`, the provenance fields
    `flight_id`, `source_structure`, and `imago_commit`
    must all be non-empty; `[entry.verification]` is
    required.  For `source = "manual"`, the curator's
    `flight_id` may be free-form;
    `[entry.verification]` is optional.
12. Every required field listed in the field tables
    above must be present.  A missing field is a hard
    error whose message names the file path, the
    offending block, and the missing field name.  This
    rule mirrors DESIGN 5.2 rule 3: the schema is
    checked before the dataclass is constructed so
    omissions surface as validation failures with full
    context, not as bare TypeError backtraces from the
    constructor.

### 7.3 Sketch (gold, single entry)

This entry is the harvest from a converged TiO2-rutile
calculation, captured at the time the dataspace was being
seeded.  It lives at `share/historicalGuidanceDB/entries/
crystalline/crystalline-a1b2c3.toml`.

```toml
schema_version = 1
entry_id       = "crystalline-a1b2c3"
generated_at   = "2026-05-28T10:30:00Z"
source         = "flight"

[entry.signature]
system_type    = "crystalline"
lattice_family = "tet"

[entry.signature.composition_vector]
alkali           = 0.0000000000000000e+00
alkali_earth     = 0.0000000000000000e+00
halide           = 0.0000000000000000e+00
chalcogen        = 6.6666666666666663e-01
pnictogen        = 0.0000000000000000e+00
group_iv         = 0.0000000000000000e+00
group_iii        = 0.0000000000000000e+00
transition_metal = 3.3333333333333331e-01
lanthanide       = 0.0000000000000000e+00
actinide         = 0.0000000000000000e+00
metalloid        = 0.0000000000000000e+00
noble_gas        = 0.0000000000000000e+00
hydrogen         = 0.0000000000000000e+00

[entry.measured]
gap_ev              = 3.0500000000000000e+00
gap_kind            = "indirect"
spin_polarization   = 0.0000000000000000e+00
total_magnetization = 0.0000000000000000e+00
kpoint_density      = 5.0000000000000000e+01

[entry.context]
basis                        = "fb"
functional                   = "gga-pbe"
kpoint_integration           = "gaussian-0.1"
scf_threshold                = 1.0000000000000000e-06
cell_atom_count              = 6
cell_volume_per_formula_unit = 4.6253846153846157e+02

[entry.verification]
grid_values = [
    2.5000000000000000e+01,
    3.5000000000000000e+01,
    5.0000000000000000e+01,
    7.5000000000000000e+01,
    1.0000000000000000e+02,
]
grid_energies = [
    -1.9512340000000000e+03,
    -1.9512378000000000e+03,
    -1.9512389000000000e+03,
    -1.9512389400000000e+03,
    -1.9512389500000000e+03,
]
converged_at           = 5.0000000000000000e+01
metric                 = "total_energy"
metric_threshold       = 1.0000000000000000e-04
predictor_confidence   = 0.0000000000000000e+00
predictor_neighbor_ids = []

[entry.provenance]
flight_id        = "guidance_seed_2026_05_28"
source_structure = "COD-1530819"
imago_commit     = "6e17c33"
curator          = "guidance_harvest.py"
```

The sketch uses 16-significant-digit float formatting per
the emitter contract in 7.5.  Note that the float values
above are shown in their idealized decimal form for
readability; the bytes the emitter actually writes are the
exact binary64 `%.16e` expansion, so a value like
`gap_ev = 3.05` appears on disk as
`3.0499999999999998e+00`, and `scf_threshold = 1e-6` as
`9.9999999999999995e-07`.  Exactly-representable values
(`0.0`, `5.0000000000000000e+01`, the 2/3 and 1/3
composition weights) are byte-identical either way.
Reading this entry:

- The composition vector says TiO2 = 2/3 chalcogen (O,
  oxygen) + 1/3 transition metal (Ti).  Exactly two
  groups carry non-zero weight; the other eleven are
  0.0.
- `lattice_family = "tet"` is rutile's tetragonal
  Bravais class.
- `gap_ev = 3.05` (eV) and `gap_kind = "indirect"`
  identify TiO2-rutile as a wide-gap indirect
  semiconductor.
- `kpoint_density = 50.0` is the converged density.
- `predictor_confidence = 0.0` and
  `predictor_neighbor_ids = []` record that this
  flight was launched without any prior to lean on --
  a seed run, no prediction was made.  A later
  flight that *did* consult the predictor would
  record a non-zero confidence and a non-empty
  neighbor list (the IDs of the k nearest neighbors
  the prediction interpolated from).

### 7.4 In-Memory Representation

**Purpose of `guidance_db.py`.**  This is the file-format
**and predictor** library: a small, passive helper module
that knows exactly two things -- how to read/validate/write
the per-entry TOML files under
`share/historicalGuidanceDB/entries/`, and how to run the
k-NN predictor over the in-memory dataspace those entries
form.  It contains no orchestration, no kaleidoscope
dispatch, and no harvest logic.  Its only runtime
dependencies are `tomllib` (Python stdlib) and the
existing `structure_control.py` (to compute composition
vectors and lattice families from a `StructureControl`).

It is imported by the flight-builder helper (consumer,
DESIGN 6.2.8), by `guidance_harvest.py` (producer), and
by `guidance_promote.py` (curator helper).  The
library / producer / consumer split keeps read-only
callers from pulling in harvest or curator code they do
not use, and isolates any future schema bump to a single
file (per ARCHITECTURE 10.6).

The module's docstring must capture this purpose
explicitly, per the project's documentation policy.

**Element-group classification table.**  The composition-
vector computation needs an element-to-group lookup.  Per
Principle 11, that table is a checked-in data file rather
than code:

```
share/historicalGuidanceDB/elemental_groups.toml
```

Format:

```toml
schema_version = 1

[groups]
alkali           = ["Li", "Na", "K", "Rb", "Cs", "Fr"]
alkali_earth     = ["Be", "Mg", "Ca", "Sr", "Ba", "Ra"]
halide           = ["F", "Cl", "Br", "I", "At"]
chalcogen        = ["O", "S", "Se", "Te", "Po"]
pnictogen        = ["N", "P", "As", "Sb", "Bi"]
group_iv         = ["C", "Si", "Ge", "Sn", "Pb"]
group_iii        = ["B", "Al", "Ga", "In", "Tl"]
transition_metal = ["Sc", "Ti", "V", "Cr", "Mn", "Fe",
                    "Co", "Ni", "Cu", "Zn",
                    "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
                    "Rh", "Pd", "Ag", "Cd",
                    "Hf", "Ta", "W", "Re", "Os", "Ir",
                    "Pt", "Au", "Hg"]
lanthanide       = ["La", "Ce", "Pr", "Nd", "Pm", "Sm",
                    "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
                    "Tm", "Yb", "Lu"]
actinide         = ["Ac", "Th", "Pa", "U", "Np", "Pu",
                    "Am", "Cm", "Bk", "Cf", "Es", "Fm",
                    "Md", "No", "Lr"]
metalloid        = []   # see DESIGN 7.10
noble_gas        = ["He", "Ne", "Ar", "Kr", "Xe", "Rn"]
hydrogen         = ["H"]
```

The library loads this table at first use and caches it
process-wide.  Every element symbol that the library is
asked to classify must appear in exactly one group; an
unclassified element is a hard error (so a typo in a
structure file fails loudly).  The `metalloid` group is
present in the schema but starts empty pending the
ambiguity resolution recorded in 7.10 (Si, B, Ge already
sit in group_iv / group_iii; whether Ge / As / Sb / Te
should move to metalloid is an open call).

**Canonical orderings.**  Two module-level constants in
`guidance_db.py` pin the index order of the
predictor's feature vectors so every consumer (load(),
save_entry(), compute_signature(), predict()) agrees on
which slot means what.  Schema rule 4 (composition vector
sums to 1.0) and the on-disk TOML representations in 7.2
both use these orderings.

```python
CANONICAL_GROUP_ORDER = (
    "alkali",
    "alkali_earth",
    "halide",
    "chalcogen",
    "pnictogen",
    "group_iv",
    "group_iii",
    "transition_metal",
    "lanthanide",
    "actinide",
    "metalloid",
    "noble_gas",
    "hydrogen",
)   # 13 element groups, the composition-vector slot order

CANONICAL_LATTICE_ORDER = (
    "cubic",
    "hex",
    "tet",
    "ortho",
    "mono",
    "tri",
)   # 6 Bravais classes, the lattice_onehot slot order
```

**Public surface (dataclasses):**

```python
@dataclass(frozen=True)
class Signature:
    """Predictor feature input for one entry."""
    system_type:        str               # "crystalline"
                                          #   / "amorphous"
                                          #   / "nanostructure"
                                          #   / "molecular"
    composition_vector: tuple[float, ...] # 13 floats,
                                          #   ordered by
                                          #   CANONICAL_GROUP_ORDER
    lattice_family:     str               # "" for non-
                                          #   crystalline;
                                          #   one of
                                          #   CANONICAL_LATTICE_ORDER
                                          #   for crystalline
    lattice_onehot:     tuple[float, ...] # 6 floats: the
                                          #   one-hot encoding
                                          #   of lattice_family
                                          #   in CANONICAL_LATTICE_ORDER.
                                          #   All zeros for
                                          #   non-crystalline.
                                          #   Derived field --
                                          #   compute_signature()
                                          #   sets it from
                                          #   lattice_family
                                          #   so the predictor
                                          #   (7.6) can use it
                                          #   directly without
                                          #   re-encoding.

@dataclass(frozen=True)
class Measured:
    """Quantities harvested from the converged calc."""
    gap_ev:              float
    gap_kind:            str             # "direct" | "indirect"
                                         #   | "none"
    spin_polarization:   float
    total_magnetization: float
    kpoint_density:      float

@dataclass(frozen=True)
class Context:
    """Calculation context recorded with each entry."""
    basis:                        str     # "mb" | "fb" | "eb"
    functional:                   str     # e.g. "gga-pbe"
    kpoint_integration:           str     # e.g. "gaussian-0.1"
    scf_threshold:                float
    cell_atom_count:              int
    cell_volume_per_formula_unit: float   # Bohr^3

@dataclass(frozen=True)
class Verification:
    """The grid that validated this entry's k-density."""
    grid_values:            tuple[float, ...]
    grid_energies:          tuple[float, ...] | None  # total
                                          #   energy (Hartree)
                                          #   per grid point,
                                          #   parallel to
                                          #   grid_values; None
                                          #   for a manual
                                          #   entry with no
                                          #   recorded sweep
    converged_at:           float
    converged_mesh:         tuple[int, ...] | None
    #                                      # resolved axial
    #                                      #   counts of the
    #                                      #   converged rung
    #                                      #   (3.12.4); None on
    #                                      #   pre-mesh or
    #                                      #   curator entries
    metric:                 str            # "total_energy"
    metric_threshold:       float
    predictor_confidence:   float          # [0.0, 1.0]
    predictor_neighbor_ids: tuple[str, ...]

@dataclass(frozen=True)
class Provenance:
    """Where this entry came from."""
    flight_id:        str
    source_structure: str
    imago_commit:     str
    curator:          str

@dataclass(frozen=True)
class GuidanceEntry:
    """One datapoint in the dataspace."""
    entry_id:     str
    generated_at: str             # ISO-8601 UTC
    source:       str             # "flight" | "manual"
    signature:    Signature
    measured:     Measured
    context:      Context
    verification: Verification | None   # None permitted
                                        #   only for
                                        #   source = "manual"
    provenance:   Provenance

@dataclass
class Dataspace:
    """The whole dataspace, loaded into memory.

    The predictor (7.6) operates on this object.  Entries
    are partitioned by system_type so the per-system_type
    sub-models can scan only their relevant subset; an
    in-memory dict keyed by system_type makes that O(1).
    """
    schema_version:           int
    entries_by_system_type:   dict[str, list[GuidanceEntry]]
    group_table:              dict[str, str]   # symbol -> group
                                               #   (cached from
                                               #    elemental_groups.toml)

@dataclass(frozen=True)
class PredictionResult:
    """What predict() returns to the flight-builder
    helper (DESIGN 6.2.8).
    """
    predicted_kpoint_density: float
    confidence:               float        # [0.0, 1.0]
    is_under_trained:         bool         # True when the
                                           #   dataspace is
                                           #   too thin for
                                           #   the predictor
                                           #   to trust its
                                           #   own answer
                                           #   (7.6 / 7.9)
    neighbor_entry_ids:       tuple[str, ...]
    predicted_gap:            float | None # None for
                                           #   non-crystalline
    predicted_magnetization:  float | None # intensive moment
                                           #   (muB/atom); None
                                           #   for non-crystalline
```

**Public surface (top-level functions):**

```python
def load(root: Path) -> Dataspace:
    """Read every entry TOML under root/entries/<system_type>/
    and the elemental_groups.toml table.  Validate per 7.2
    rules 1-12, return the loaded Dataspace.  Raises
    GuidanceDataspaceError on any validation failure
    with the filename and the failed rule cited.
    """

def save_entry(entry: GuidanceEntry, root: Path) -> Path:
    """Emit `entry` as TOML into
    root/staging/<system_type>/ using the deterministic
    hand-formatter (7.5).  Returns the written path.
    Raises if a file with the same `entry_id` already
    exists.
    """

def compute_signature(
    structure:   StructureControl,
    system_type: str,
    group_table: dict[str, str],
) -> Signature:
    """Compute the Signature for a given StructureControl.
    Atom-fraction across the 13 element groups using
    group_table; the lattice_family for crystalline is
    read off the StructureControl's Bravais lattice
    detection.  Raises GuidanceDataspaceError if any
    element symbol is missing from group_table (rule
    enforced at compute time so the failure point names
    the structure, not the dataspace load).
    """

def predict(
    dataspace:           Dataspace,
    query:               Signature,
    basis:               str,
    functional:          str,
    kpoint_integration:  str,
) -> PredictionResult:
    """Run the predictor (7.6) for a given query signature
    within the (basis, functional, kpoint_integration)
    sub-model.  The three settings select the sub-model
    (7.6 step 2) so a prediction never interpolates across
    incompatible settings -- in particular a tetrahedral-
    integration density is never mixed with a Gaussian-
    smeared one.  Always returns a PredictionResult; the
    `is_under_trained` flag plus the `confidence` score
    tell the caller how seriously to take the prediction.
    Never returns None: the caller (DESIGN 6.2.8) must
    decide whether to verify-around the prediction, fall
    back to the wide-grid default (7.9), or refuse to
    proceed.
    """
```

### 7.5 Hand-Formatted TOML Emitter

The emitter is hand-written -- not delegated to a third-
party TOML writer -- for the same reasons as DESIGN 5.5:
bit-deterministic output (so version-control diffs are
meaningful), tight control over float formatting (so
numerical comparisons are unambiguous), and freedom from a
third-party dependency.

**Emitter contract:**

- Floats are written with the format string `"%.16e"`,
  yielding 16 significant digits in scientific notation.
- Integers are written as bare decimals.
- Strings are written in TOML basic-string form
  (double-quoted).  Quotes and backslashes inside strings
  are escaped per TOML spec.
- The key order within a block is fixed (per the field-
  list order in 7.2), so the same in-memory entry always
  produces byte-identical TOML output.
- The top-level block sequence is fixed:
  `[entry.signature]` (with its sub-block
  `[entry.signature.composition_vector]` written as
  multi-line inline-table-like form -- one key per line,
  the 13 group keys in canonical order),
  `[entry.measured]`, `[entry.context]`,
  `[entry.verification]` (when present),
  `[entry.provenance]`.
- Arrays of floats are written one element per line, with
  4 leading spaces of indent and a trailing comma after
  every element (including the last), to make per-element
  diffs minimal.  Arrays of strings or integers stay
  inline.
- A blank line separates top-level blocks; no blank lines
  inside a block.

**Slug derivation for entry filenames** (and `entry_id`):

```
slug = system_type + "-" + short_sha
short_sha = first 6 hex digits of SHA-256 over the bytes
            (flight_id || source_structure || generated_at)
```

Two virtues of putting `system_type` in the slug rather
than chemistry: (1) it reflects the on-disk partition
(`entries/<system_type>/<system_type>-<sha>.toml`), so a
human glancing at a single file path can tell what kind
of system it describes; (2) it removes the variable-
length elements_part of the previous design, so every
slug is exactly the same length (about 20 chars).

The `short_sha` is the collision guard discussed in
ARCH 10.8: two flights harvesting an entry at the same
instant produce different hashes (because either
`flight_id` or `source_structure` will differ), so
their files do not collide.  If by extreme coincidence
they do, `save_entry` raises a hard error (rule 2) and
the harvest script retries with a fresh `generated_at`.

### 7.6 Predictor Algorithm

The predictor answers: "given a target system's feature
vector and the (basis, functional) it will be run under,
predict the converged k-density and tell the caller how
confident the prediction is."

The predictor is **k-nearest-neighbor regression with
inverse-distance weighting**, run in **two stages for
crystalline** systems and as a simple per-bucket
canonical for non-crystalline.  The two-stage split for
crystalline exploits the transferability argument from
7.1: chemistry maps to electronic character (stage 1);
electronic character maps to k-density (stage 2).  Each
stage is a separate k-NN, with its own neighbor set,
distance metric, and variance.

**Step 1 -- partition by system_type.**  The Dataspace
(7.4) is partitioned by system_type at load time.  The
predictor first switches on the query's system_type:

- `crystalline`: run the two-stage regression below.
- `amorphous` / `nanostructure` / `molecular`: return the
  canonical entry for that system_type (typically a
  Gamma-floor density driven by the cell-volume
  convention; see 7.9 for the exact day-1 canonical
  values).  Chemistry plays little role here -- the
  density convention dominates -- so the predictor's
  job collapses to a constant.

The rest of this section describes the crystalline path.

**Step 2 -- sub-model selection by (basis, functional,
kpoint_integration).**
The predictor maintains a sub-model per (basis,
functional, kpoint_integration) triple: the k-NN draws only
on entries whose context matches.  Justification: changing
the basis, the functional, or the Brillouin-zone integration
method can shift the converged k-density and the measured
gap meaningfully, and we do not want
interpolation across them to wash out that signal -- a
density converged under analytic tetrahedral integration is
not interchangeable with one converged under Gaussian
smearing.

If the queried (basis, functional, kpoint_integration)
sub-model has fewer than `k_min = 3` entries, the predictor:

- Falls back to the most-populous sub-model under the
  same functional family (e.g. `(mb, gga-pbe, *)` ->
  `(fb, gga-pbe, *)` if mb is sparse).  This degraded
  fallback ignores basis and integration; it is the
  best-effort path when an exact sub-model is too thin.
- If no functional-family fallback exists, falls back
  to the system_type's overall pool (ignoring context).
- If the overall pool also has fewer than `k_min`
  entries, returns
  `PredictionResult(is_under_trained = True, ...)`.

**Step 3 -- Stage 1, chemistry -> electronic character.**
For each entry `E` in the sub-model, define the stage-1
distance to the query `Q`:

```
d1(Q, E) = sqrt(
    w_comp * || Q.composition_vector
              - E.composition_vector ||^2
  + w_latt * || Q.lattice_onehot
              - E.lattice_onehot ||^2 / 2.0
)
```

Both `Q.composition_vector` and `E.composition_vector`
are 13-vectors summing to 1.0 (per schema rule 4), so the
squared Euclidean distance is a well-defined chemistry
similarity in [0.0, 2.0].  Both `Q.lattice_onehot` and
`E.lattice_onehot` are 6-vectors in canonical Bravais
order (cubic, hex, tet, ortho, mono, tri), with exactly
one entry equal to 1.0 and the others 0.0 -- derived from
the entry's `lattice_family` string at
`compute_signature` time (7.4).  The squared Euclidean
distance between two one-hot vectors is 0.0 if they match
and 2.0 if they differ; the `/2.0` normalizes that term
to the same [0.0, 1.0] dynamic range as the
composition-distance contribution, so the `w_latt`
default of 0.25 reads as "lattice contributes up to 25%
of the composition-term weight on a full mismatch."

**Why one-hot for lattice_family rather than a string-
equality test.**  Operationally equivalent for the
single-feature case, but materially more extensible:
adding a second categorical feature later (space-group
class, defect-host indicator, ...) is just a longer
concatenated vector and one more weight, and swapping
the predictor for a fancier model (learned distance,
random forest) treats the one-hot directly as numpy
features without bespoke preprocessing.

Default weights: `w_comp = 1.0`, `w_latt = 0.25`.  These
make composition the dominant signal and let lattice
family separate polytypes (rutile-TiO2 from anatase-TiO2)
without dominating.  Both are tunable; calibration after
the seed flight may shift them.

Find the `k = 5` nearest neighbors by `d1`.  Stage-1
predictions are inverse-distance-weighted means:

```
weights:      w_i = 1.0 / (d1(Q, E_i) + epsilon)
              normalized so sum(w_i) = 1.0
predicted_gap        = sum(w_i * E_i.measured.gap_ev)
predicted_magnetization =
    sum(w_i * |E_i.measured.total_magnetization|
              / E_i.context.cell_atom_count)
```

with `epsilon = 1e-6` to avoid division by zero on an
exact match.

**Why the intensive magnetization (`|M| / N_atoms`) and not
`spin_polarization`.**  The second character feature is a
proxy for "is this a magnetic metal whose Fermi surface needs
a denser mesh."  imago surfaces the magnetic *moment* (the
iteration file's column 6, always written), never a
spin-polarization fraction, so `measured.spin_polarization` is
structurally 0.0 for every harvested entry -- keying on it
would be keying on a dead feature.  The moment itself is
*extensive* (a primitive cell and its N-fold supercell of the
same magnet report N-fold-different totals while needing the
same k-density), so it is divided by the cell atom count to
make it intensive, and taken in magnitude because the up/down
spin labeling is an arbitrary SCF choice.  v1 limitation: an
antiferromagnet has `M = 0` total yet is locally
spin-polarized, so this feature reads it as non-magnetic
(7.10).

Stage-1 **confidence** is derived from the weighted
variance of the neighbors' gap values:

```
gap_variance = sum(w_i * (E_i.measured.gap_ev
                           - predicted_gap)^2)
gap_spread   = sqrt(gap_variance)
confidence_1 = exp(-gap_spread / sigma_gap_ref)
```

with `sigma_gap_ref = 1.0` eV (the gap range over which
the predictor is "comfortably confident").  This gives
`confidence_1` in (0.0, 1.0]: 1.0 when the neighbors
agree perfectly on gap, decaying as they disagree.

**Step 4 -- Stage 2, electronic character -> k-density.**
For each entry `E` in the same sub-model, define a
stage-2 distance over the predicted electronic character:

```
d2(Q, E) = sqrt(
    w_gap  * (predicted_gap - E.measured.gap_ev)^2
              / sigma_gap^2
  + w_spin * (predicted_magnetization
              - |E.measured.total_magnetization|
                / E.context.cell_atom_count)^2
              / sigma_spin^2
)
```

with `sigma_gap = 1.0` eV and `sigma_spin = 0.5` (Bohr
magnetons per atom).  Note the asymmetry: stage 2's
distance uses the **predicted** character (from stage 1),
not the query's chemistry, because the goal is "find
calculations whose gap-and-magnetization look like what
this query is likely to produce."  Default weights:
`w_gap = 1.0`, `w_spin = 0.5`.

Find the `k = 5` nearest neighbors by `d2`.  Predicted
k-density is the inverse-distance-weighted mean:

```
predicted_kpoint_density =
    sum(w_i * E_i.measured.kpoint_density)
```

Stage-2 **confidence** is derived analogously:

```
kpd_variance = sum(w_i * (E_i.measured.kpoint_density
                           - predicted_kpoint_density)^2)
kpd_spread   = sqrt(kpd_variance)
confidence_2 = exp(-kpd_spread / sigma_kpd_ref)
```

with `sigma_kpd_ref = 50.0` (k-points/Bohr^-3 over
which the predictor is "comfortably confident" about
density).

**Combined confidence.**  The two stages compound:

```
confidence = confidence_1 * confidence_2
```

both in (0.0, 1.0], product also in (0.0, 1.0].

**Returning the result.**  The predictor returns a
`PredictionResult` (7.4) carrying
`predicted_kpoint_density`, `confidence`,
`is_under_trained = False`, the union of the stage-1
and stage-2 neighbors' `entry_id`s (deduplicated), and
the intermediate `predicted_gap` and
`predicted_magnetization` for forensics.

**Tuning knobs and their defaults.**  All of these are
named module-level constants in `guidance_db.py`, so
calibration after the seed flight is a one-file
change:

```
k_min         = 3      # below this many entries, refuse
                      #   the sub-model
k_neighbors   = 5      # neighbors used at each stage
epsilon       = 1e-6   # numerical floor on distance
w_comp        = 1.0    # composition weight in d1
w_latt        = 0.25   # lattice-family weight in d1
w_gap         = 1.0    # gap weight in d2
w_spin        = 0.5    # spin weight in d2
sigma_gap     = 1.0    # gap normalization (eV) in d2
sigma_spin    = 0.5    # spin normalization in d2
sigma_gap_ref = 1.0    # gap-spread -> confidence_1
sigma_kpd_ref = 50.0   # kpd-spread -> confidence_2
```

**Why this shape and not something more sophisticated.**
k-NN with inverse-distance weighting is the simplest
predictor that handles sparse, non-stationary, low-
dimensional data well; produces a natural variance-based
confidence; degrades gracefully (with fewer than `k`
neighbors, predict-then-verify just widens the grid via
low confidence rather than failing); and stays auditable
(every prediction is "the weighted average of these 5
listed neighbors").  Linear regression would assume
smoothness we cannot defend across polytypes; neural
networks would obscure the audit trail and require more
data than the seed flight produces.  Calibration after
the seed flight will tell us whether k-NN's accuracy
floor is acceptable or whether a more sophisticated
model is warranted.

### 7.7 Predict-then-Verify Flight Construction

This subsection covers the bridge from a single prediction
to a concrete kaleidoscope flight.  The flight-builder
helper that lives in `src/scripts/kaleidoscope/` (DESIGN
6.2.8) is what calls into this design rung; the algorithm
below is what it executes.

**The verify search itself is specified in 3.12.**  Choosing
the rungs and dispatching them is no longer a single fixed
density grid: it is the adaptive mesh climb (3.12), which
searches in mesh space and stops when the energy is flat.  A
confident prediction instead dispatches a small *fixed mesh*
grid around the predicted rung in one parallel round (3.12.5);
a cold or under-trained prediction climbs from a low start.
What stays in this subsection is the surrounding flow: computing
the query signature, calling the predictor, selecting the sub-
model, and attaching the per-structure `PredictionRecord` the
harvest recovers later (7.8).  Where the algorithm below
constructs a fixed density grid, read that as the density-grid
special case the climb generalizes -- the predicted density is
still what seeds the search, per 3.12.4.

**Inputs:**

- `target`: a `StructureControl` for the system to converge.
- `options`: a dict of `makeinput` options (everything
  *but* the swept knob).
- `system_type`: one of the four valid values, declared by
  the caller.
- `basis`, `functional`, `kpoint_integration`: the
  (basis, functional, kpoint_integration) sub-model under
  which the flight will run; the triple selects the
  predictor's sub-model (7.6 step 2) and is the same triple
  the per-structure record carries (step 5) and the harvest
  reads back (7.8 step 3f).
- `dataspace`: the loaded `Dataspace` (7.4).
- `verify`: optional bool, default True; False triggers
  trust mode (6.2.1) -- a length-1 grid at the predicted
  point.

**Outputs:**

- A `Flight` of `CalcUnit`s (DESIGN 6.2.1) ready to
  dispatch.
- A `PredictionRecord` (7.4-derived) the harvest hook
  recovers later.

**Algorithm:**

```
1.  query_sig = compute_signature(target, system_type,
                                  dataspace.group_table)

2.  prediction = predict(dataspace, query_sig, basis,
                         functional, kpoint_integration)

3.  if not verify:
        # trust mode (DESIGN 6.2.1).
        grid_values = [prediction.predicted_kpoint_density]
        policy      = "trust_no_verify"

    elif prediction.is_under_trained:
        # the dataspace is too thin for the predictor to
        # trust its own answer -- fall back to the wide-
        # grid default (7.9).
        grid_values = default_wide_kpoint_density_grid()
        policy      = "wide_grid_no_prior"

    else:
        grid_values = build_verification_grid(
                          prediction.predicted_kpoint_density,
                          prediction.confidence,
                      )
        policy = "verify_around_prediction"

4.  Round each grid value to an integer k-density and
    dedupe (so the 6.2.4 tag parses back to exactly the
    swept value; the 7.9 wide-grid defaults are already
    integers, and build_verification_grid's logspace
    floats round here -- a degenerate grid where rounding
    merged two close points collapses):

        kpd_grid = sorted(set(round(v) for v in grid_values))
        units = [
            build_calc_unit(target, options, kpd = kpd_int)
            for kpd_int in kpd_grid
        ]

    Two-name convention worth pinning here:
    - `kpd` is the makeinput options-dict key (matches
      makeinput.py's argparse dest for `-kpd`).  Used
      anywhere the value is passed to makeinput.
    - `kpt-density` is the display name used inside
      the calc-tag tree (per 6.2.4's tag convention)
      and inside Flight.sweep.varied_axes.  Used
      anywhere the axis appears as a directory level
      or as a sweep-axis name humans inspect.
    The flight-builder helper (6.2.8) is the
    translation point between the two names.

5.  prediction_record = PredictionRecord(
        policy                   = policy,
        predicted_kpoint_density =
            prediction.predicted_kpoint_density
                              if prediction is not None else None,
        confidence               = prediction.confidence
                              if prediction is not None else 0.0,
        is_under_trained         = (prediction.is_under_trained
                              if prediction is not None
                              else True),
        neighbor_entry_ids       = prediction.neighbor_entry_ids
                              if prediction is not None else [],
        predicted_gap            = prediction.predicted_gap
                              if prediction is not None else None,
        predicted_magnetization  =
            prediction.predicted_magnetization
                              if prediction is not None else None,
        system_type              = system_type,
        feature_vector           = query_sig,
        basis                    = basis,
        functional               = functional,
        kpoint_integration       = kpoint_integration,
    )

    The last three fields -- the (basis, functional,
    kpoint_integration) sub-model the prediction was made
    under -- are recorded ONLY here, on the per-structure
    record; they are deliberately NOT also copied into the
    flight-level `sweep.fixed_axes` (6.2.9).  This is the
    single home for the sub-model, which both lets a producer
    fold many structures into one combined flight even when
    they do NOT share a sub-model -- each structure's harvest
    reads its own sub-model back from its own record (7.8
    step 3f) -- and avoids the confusion of the same fact
    living in two places.

6.  flight = Flight(
        units = units,
        sweep = SweepRecord(
            varied_axes = ("kpt-density",),
            fixed_axes  = {},
        ),
        ...
    )
    attach_prediction_record(flight, prediction_record)
    The `sweep` field (DESIGN 6.2.1) makes
    serialize_flight emit [flight.sweep], so harvest
    (7.8 step 3a) recovers the varied axis without parsing
    run-dir paths.  `fixed_axes` is left empty: the sub-model
    that used to live here now rides on the per-structure
    record (above), so a single-structure and a combined
    multi-structure flight share one shape and the sub-model
    is never duplicated (6.2.9).  (`fixed_axes` remains a
    general SweepRecord field for any future axis a flight
    genuinely holds constant across every unit; it simply has
    no v1 occupant.)

7.  return flight, prediction_record
```

**The verification-grid widening function.**  Now driven
by the predictor's `confidence` rather than a chemistry-
similarity score.  Same shape as before; the semantic
input is different:

```
def build_verification_grid(center: float,
                            confidence: float) -> list[float]:
    """Return a list of k-density values to sweep around
    `center`, with grid width scaling inversely with
    predictor confidence.
    """
    # Width is the multiplicative span: 1.5 means values
    #   from center/1.5 to center*1.5.
    width    = 1.2 + 1.5 * (1.0 - confidence)
    # Number of points scales similarly.
    n_points = round(3 + 4 * (1.0 - confidence))

    # Logarithmically spaced for symmetric coverage.
    lo, hi = center / width, center * width
    return logspace(lo, hi, n_points)
```

Behavior at the extremes:

- `confidence = 1.0` (perfect match, neighbors agree):
  width = 1.2, 3 points, span [center/1.2, center*1.2].
  Tight verification that the prior still holds.
- `confidence = 0.7` (good match, modest spread):
  width = 1.65, ~4 points.
- `confidence = 0.3` (neighbors disagree significantly):
  width = 2.25, 6 points, span [center/2.25, center*2.25].
  Wide enough that if the prediction was a poor guide,
  the true converged point is still likely in range.

The exact constants (`1.2`, `1.5`, `3`, `4`) are starting
heuristics; calibration after the seed flight may
adjust them.  They are tunable knobs in
`src/scripts/kaleidoscope/`, kept in the one function so
calibration is a one-file change.

**The prediction record** is persisted alongside the
flight (in `flight.toml` as `[flight.predictions.<id>]`,
keyed by structure id so one flight can carry many; 6.2.9)
so the harvest hook (7.8) can recover the predicting
neighbors and the confidence score that drove the grid
choice.  Without it, the harvested
`predictor_confidence` and `predictor_neighbor_ids`
fields would be unrecoverable.

**Dispatching the climb: from a mesh to a run.**  The climb
searches in mesh space (3.12), so the unit it dispatches names an
explicit mesh, not a density.  makeinput already accepts one: the
`scfkp = [a, b, c]` option writes a style-code-1 k-point file
(axial counts, shift, and point operations), and imago resolves its
own symmetry shift and irreducible-wedge reduction from those counts
exactly as it would for a density-selected mesh (2.4).  The
requested axial counts therefore ARE the resolved full mesh -- there
is no density-to-mesh rounding to invert -- which is the whole
reason the climb can search in mesh space (3.12.1).

Two-name convention, mirroring the density path's `kpd` /
`kpt-density` above:

- `scfkp` is the makeinput options-dict key (its argparse dest),
  carrying the three counts wherever the value reaches makeinput.
- `kpt-mesh` is the display name in the calc-tag tree (6.2.4) and in
  `Flight.sweep.varied_axes`.  A unit's swept value renders as the
  single tag component `kpt-mesh-<a>-<b>-<c>` (e.g. `kpt-mesh-4-4-4`)
  and reads back to the count triple by the inverse; the
  hyphen-separated counts stay slug-safe because axial counts are
  always positive integers.

Reading a rung back needs two facts from each completed run's
`result.toml` (6.1.2): `total_energy` is the rung's energy (the
basis the flatness test normalizes to eV per atom, 7.8 step 3c) and
`kpoint_mesh` is the resolved axial counts.  Because an explicit
mesh is honoured exactly, `kpoint_mesh` must equal the requested
mesh; the adapter asserts this, so a makeinput or imago change that
silently altered the mesh is caught rather than mis-recorded.  (It
is the same resolved-mesh line whose companion `RESOLVED_KP_CLASSES`
emit validates the producer's axis-class port, 2.7.)

**The climb dispatcher.**  The climb's dispatcher (3.12.5) is a
thin producer-side adapter over the ordinary dispatch layer,
exposing two calls the climb loop drives (4e.5).  `send(mesh_lists)`
takes `{material: [mesh, ...]}`, builds one CalcUnit per mesh
(below), appends them to the one Flight that spans the whole climb,
runs the driver-side prepare pass on just those new units, and
launches them without waiting (`send_off`, 6.2.3).  `next_rung()`
blocks until the next rung lands (`collect_next`, 6.2.3), reads its
completed unit's `(mesh, total_energy)` into a rung, and returns
`(material, rung)` -- or `(material, FAILED)` for a unit that did
not complete.  Because one Flight spans the climb and its unit list
accretes as rungs are decided, `flight.toml` records every rung
asked for; the dispatch core stays domain-ignorant (Principle 12),
running whatever units it is handed.  A mesh already run is a cache
hit (6.2.5), so re-launching it costs nothing and the climb never
has to track what it has already run.

**Splitting the builder.**  The density-era flight builder (6.2.8;
the algorithm above) does prediction AND grid-laying in one call.
The climb separates the two, because it seeds from the prediction
but lays its own rungs:

- `predict_kpoint_density` runs steps 1-2 and 5 above -- the query
  signature, the predictor call, and the per-structure
  `PredictionRecord` -- and returns the predicted density,
  confidence, under-trained flag, and record.  The producer uses the
  density to seed the climb (3.12.4) and the confidence to pick the
  dispatch mode and the persistence (3.12.6); the record travels to
  the harvest as a plain dict the producer stamps per solid (7.8).
- `build_mesh_unit` builds one explicit-mesh CalcUnit -- the `scfkp`
  option, the `kpt-mesh` tag, and the same cache identity (6.2.1)
  the density units used -- and the dispatcher calls it once per
  mesh as the climb sends it.

The single-call grid builder was the density-era special case the
climb generalizes: its confidence-widened grid becomes the confident
mode's small fixed grid (3.12.5).  The producer now drives the climb,
so that builder is retired; `predict_kpoint_density` and
`build_mesh_unit` are the two pieces that remain.

**A rung that fails to run.**  A unit that does not complete has no
`result.toml`, hence no energy, so the adapter returns only the rungs
that completed.  A still-active material that asked to run a mesh and
did not get it back cannot climb further, so the producer stops it as
NON_CONVERGED with a run-failure reason -- distinct from a ceiling
stop (3.12.3) -- and the round loop moves on.  This keeps the loop
from re-dispatching a failing mesh forever.  It is deliberately not a
retry: recovering a flaky run is the runner's and the custodian's
job, above the domain-ignorant dispatch layer (Principle 12), not the
climb's.

### 7.8 Harvest Pipeline (Staging and Promotion)

The harvest hook turns a finished flight into a staged
guidance entry rich enough to feed the predictor.  It
runs after the verification grid has finished (or as a
separate post-step the user invokes).

This is the *guidance* harvest specifically, and it is
distinct from the producer's own potential harvest (5.7).
A guidance entry's whole content is the claim "for a
structure like this, this k-density is converged," and the
predictor's training target *is* that converged k-density.
Convergence is a statement about a neighborhood -- the
total energy has stopped moving as the mesh is refined --
so it is only established by a grid (the two-sided rule of
step 3c needs at least three points).  This is why the
harvest needs a flight that declared a sweep with a varied
k-density axis.

A single one-off calculation -- a curator who does not want
to sweep variants -- is *not* blocked by this.  Such a run
is built as a length-1 sweep (trust mode, or a pinned
`kpoint_spec` override; 6.2.1 / 5.7); the producer's own
harvest (5.7) still extracts its converged potential for
the initial-potential database, and the guidance harvest
simply does not auto-stage a guidance entry from it (step
3a) -- one point is weaker evidence than a grid.  When the
curator already knows a good k-density and wants it in the
guidance dataspace from a one-off, the manual seed path
(7.9, `source = "manual"`) records it directly, with a
human vouching for the convergence claim the automation
could not establish.

**Inputs:**

- A finished flight's workspace directory.
- The flight's `[flight.predictions.<id>]` tables recovered
  from `flight.toml` (7.7), one per structure.

**The three-source rule (Model 1).**  Every entry field is
filled from exactly one of three on-disk inputs, so the
information flow stays simple and homogeneous:

- **`flight.toml`** -- the *plan*: the unit list, the
  `[flight.sweep]` block (the varied axis the swept value is
  read along), and the `[flight.predictions.<id>]` tables
  (one per structure, each carrying that structure's
  prediction *and* the (basis, functional,
  kpoint_integration) sub-model it ran under; 6.2.9).
  Each grid point's swept k-density is read out of its calc
  tag (`kpt-density-<int>`) via the sweep's ordered
  `varied_axes`; the makeinput `options` are not persisted
  in `flight.toml`, so the calc tag is the on-disk source of
  the swept value.
- **each converged run's `result.toml`** -- the *per-run
  facts* the imago callable API (DESIGN 6.1) exposes:
  - `gap_ev`, `gap_kind` (read off the eigenvalue spectrum;
    TODO C76 for the imago-side wire-up).
  - `total_magnetization` (always written; closed-shell runs
    report 0.0).  `spin_polarization` is NOT surfaced -- the
    iteration file carries the magnetic *moment*, not a
    polarization -- so the entry records `spin_polarization =
    0.0` and the predictor keys its spin character on
    `total_magnetization` (DESIGN 7.6).
  - the SCF total energy (used to pick the converged grid
    point).
  - `scf_threshold` -- the SCF criterion the run converged
    to, recorded in the entry's context.  It is a distinct
    criterion from the grid-flatness `metric_threshold`
    (below): one governs a single SCF run, the other judges
    energy versus k-density.
- **the structure `.skl`** -- the *structural facts*: the
  harvest loads it anyway for `compute_signature`, and the
  same load yields `cell_atom_count` (`num_atoms`) and
  `cell_volume_per_formula_unit` (the cell volume in Bohr^3,
  formula-unit count Z = 1 in v1; curator-facing metadata
  the predictor never reads).

When a per-run quantity is absent (a closed-shell run, an
older imago version), the harvest records 0.0 for the
measured value and falls back to `"unknown"` for an absent
`imago_commit` (non-empty, so the schema's rule-11 check
still passes and the curator can spot it on review).

*Where `imago_commit` comes from.*  It is a recorded fact, not a
measured one, and it reaches `result.toml` because the wingbeat
echoes it there out of the unit's `record` mapping (6.2.2 /
6.2.4).  Reading it beside the measured quantities is what keeps
this harvest on three sources rather than four -- it never opens
the dispatch core's `status.toml`.  The build is deliberately not
part of the cache key (6.2.5) and not part of the claim key
(below); recording it and partitioning on it are different acts,
and only the first is wanted.  Until TODO C84 has the binary
report its own build, the value is what the *producer believed*
it launched, which is worth strictly more than `"unknown"` and
strictly less than the engine's own word; C84 replaces it in
place, in this same field.

**Algorithm** (`guidance_harvest.py`):

```
1.  Load the flight report (DESIGN 6.2.6) and the
    [flight.predictions.<id>] tables from flight.toml
    (one per structure).

2.  Keep only the convergence-sweep runs
    (`kind == "convergence"`, 6.2.9) -- other kinds (e.g.
    "fingerprint" loen runs) share a structure id but
    belong to a different harvester -- then group those
    CalcUnits by id (one group per structure).

3.  For each structure group (let `prediction` be its
    entry in the predictions mapping).  If the structure
    has no `prediction` record, SKIP it: the record is now
    the sole source of both `system_type` (step e) and the
    (basis, functional, kpoint_integration) sub-model (step
    f), so a record-less structure cannot be staged.  A
    flight built by the helper (6.2.8) always carries a
    record per structure; a record-less convergence sweep is
    a hand-built flight outside the predict-then-verify path,
    and guidance entries are earned only along that path
    (7.9 covers the by-hand seed route instead).  Otherwise:
      a. Sort the grid by k-density ascending.  A
         single-point grid (trust mode, or a single-point
         curator override) harvests deliverables but is
         NOT staged as a guidance entry -- one converged
         calc is weaker evidence than a grid (6.2.1) --
         and is skipped here before the convergence test.
      b. For each CalcUnit's run, parse result.toml for
         total_energy, gap_ev, gap_kind,
         total_magnetization, scf_threshold, and
         kpoint_mesh (the resolved axial counts, 6.1.2 /
         PSEUDOCODE 4c.6); read the swept k-density out of
         the CalcUnit's calc tag.

         **DROP a point whose SCF did not converge.**  A run
         that finished is not the same as a run that
         converged, and the two questions are answered in
         different places: the flight entry's status says
         whether the job completed, while `converged` in
         result.toml says whether the SCF reached its own
         fixed point.  A non-converged run finishes normally
         and writes a total energy, so nothing distinguishes
         it downstream -- its energy is simply wherever the
         iteration happened to stop.  Feeding that into a
         flatness test is worse than useless: the test asks
         whether the energy has stopped moving with the mesh,
         and this energy was not finished moving for reasons
         that have nothing to do with the mesh.  It can read
         flat by accident and it can break a genuine plateau.
         Only an EXPLICIT `converged = false` drops a point.
         A result.toml carrying no such field cannot be
         judged, and is kept -- the same side this design
         takes on a missing `gap_ev` (3.12.3), and for the
         same reason: a missing reading must never silently
         discard evidence a real run earned.  Dropping is
         reported, never silent.
      c. Collapse duplicate-mesh rungs, then pick the
         converged grid point.  A requested k-density does
         not map one-to-one onto the mesh imago integrates:
         the density-to-mesh map (3.7) is a step function,
         so a finer density can resolve to a mesh already
         seen.  Two rungs with an equal `kpoint_mesh` are
         the same calculation run twice; their energy delta
         is exactly zero, and the two-sided test below would
         read that manufactured zero as convergence
         (DESIGN 3.11).  So first reduce the density-sorted
         grid to one rung per distinct `kpoint_mesh`,
         keeping the lowest-density member -- the cheapest
         request that reaches that mesh, and the right
         convergence density to record.  Because the mesh is
         monotone in density (3.7), equal meshes occupy a
         contiguous density range, so this simply merges
         neighbours.  Equal meshes give equal energies, so
         nothing is lost; a total-energy disagreement
         between two runs of the same `kpoint_mesh` is not
         averaged away but surfaced as an error, since it
         means the runs were not in fact identical.
         Then pick the converged grid point on the collapsed
         grid: the smallest k-density at which the PER-ATOM
         energy change to both neighbours is below the
         k-point threshold -- |E_i - E_{i+1}| /
         cell_atom_count < metric_threshold AND
         |E_i - E_{i-1}| / cell_atom_count < metric_threshold
         for i in (1, len(collapsed_grid)-1), the deltas
         taken in eV.  Requiring both consecutive-pair
         deltas to be small mitigates a single-grid-point
         numerical fluke, and collapsing first guarantees
         those deltas compare genuinely distinct meshes
         rather than a mesh with itself.  Collapsing can
         shrink the grid below the three points the interior
         test needs; when it does, there is no interior
         point and step (d) applies.  `metric_threshold` is
         the k-point convergence threshold in eV/atom,
         resolved from the solid's
         `kpoint_convergence_threshold` (its own value, else
         the manifest `[harvest]` block, else the built-in
         1 meV/atom = 1e-3 eV/atom default; 5.7) --
         separate from `scf_threshold`, which governs one SCF
         run, not the flatness of energy versus k-density.
      d. If no point satisfies the criterion (energy
         still moving at the top of the grid), log a
         warning, tag the flight with
         `prediction_mismatch = true`, and SKIP this
         structure.  Non-converged sweeps do not earn
         an entry.  The user must widen the grid and
         re-run.
      d'. If the structure reads as a metal, log it and
         SKIP: a metal stages no guidance entry.  See
         below.  TWO independent readings can say so, and
         EITHER is sufficient.  The first is the chosen
         run's own gap (`gap_ev` at or below
         `metal_gap_threshold`, the same test the climb's
         metal short-circuit applies, 3.12.3).  The second
         is the caller's LADDER reading: whether any rung
         it computed came back gapless, which this step
         receives as an argument rather than re-deriving.
         Each caller supplies it from the multi-rung
         evidence it has -- the producer from the climb's
         verdict (3.12.3), this standalone harvest from the
         gaps of every point in its own grid.
         **Why one reading is not enough.**  A metal on a
         discrete mesh shows a small artificial gap whose
         size depends on where the mesh points fall (1.6),
         so a SINGLE rung's reading is close to a coin
         toss: fcc Al reads zero at several meshes and
         0.124 eV at another.  The chosen rung alone is
         therefore weak evidence, and the ladder taken
         whole is strong.  Taking either as sufficient
         means neither reading can be outvoted by a worse
         one, and nothing currently caught stops being
         caught.  The threshold
         rides on this structure's `prediction` record,
         the same channel `kpoint_convergence_threshold`
         uses and for the same reason: it is a manifest
         knob (`[kpoint_climb]`, 5.7) and this harvest is
         a standalone tool pointed at a finished
         workspace, so it never sees a manifest.  The
         producer stamps the resolved value onto the
         record when it builds the flight.  An ABSENT
         `gap_ev` is not metallic: a missing reading must
         not silently suppress an entry a real insulator
         earned (3.12.3 takes the same side).
      e. Compute the structure's signature: system_type
         from this structure's `prediction` record
         (which carries it from 7.7); composition_vector
         and lattice_family via compute_signature().
      f. Build a GuidanceEntry (per the three-source rule):
            - signature: from step (e)
            - measured: gap/magnetization from the chosen
              run's result.toml; kpoint_density = the
              chosen calc tag's k-density
            - context: basis/functional/kpoint_integration
              from this structure's `prediction` record
              (which carries the sub-model it ran under;
              7.7 / 6.2.9), NOT from the flight-level
              fixed_axes, so a combined mixed-sub-model
              flight harvests each structure correctly;
              scf_threshold from result.toml; cell_atom_count
              and cell_volume_per_formula_unit from the
              loaded structure (step e)
            - verification: grid_values,
              grid_energies (the parallel total-energy
              array over the COLLAPSED distinct-mesh grid
              of step c, so the auto-promote rule judges
              flatness on genuinely distinct calculations
              and never re-encounters the duplicate-mesh
              zero the harvest already removed),
              converged_at = chosen k-density,
              metric, and metric_threshold -- the
              resolved kpoint_convergence_threshold,
              read from this structure's `prediction`
              record.  That record is the channel: the
              harvest is a standalone tool pointed at a
              finished workspace (7.8) and never sees the
              manifest, so the producer stamps the
              resolved value onto the record when it
              builds the flight, exactly as it does
              system_type, the sub-model, and the
              `metal_gap_threshold` step 3d' applies;
              predictor_confidence and
              predictor_neighbor_ids from this
              structure's `prediction` record.
            - provenance: flight_id, source_structure,
              imago_commit, curator = "guidance_harvest.py".
      g. Write the entry to
         share/historicalGuidanceDB/staging/<system_type>/
         via save_entry().

4.  Print a one-line summary per structure
    (converged / skipped / staged path).
```

**Two thresholds, not one.**  The grid-flatness
`metric_threshold` is deliberately separate from the SCF's
`scf_threshold`.  A single SCF converges its own iteration
residual to `scf_threshold`; a total energy of order tens of
hartree, though, cannot flatten to that same figure versus
k-density -- 1e-6 hartree is ~1e-8 relative, far below real
k-point sampling noise -- so binding the two lets only
accidentally-flat structures harvest while genuinely slow
ones (their energy still moving at the 1e-4 hartree level)
are reported non-converged even though every SCF converged
cleanly.  The threshold is *per atom* because the database
mixes cells of very different size (2 to 16+ atoms), and a
fixed absolute total-energy tolerance would demand
ever-tighter relative convergence as the cell grows.

The 1 meV/atom default is the textbook bar, arrived at from
below rather than chosen for its familiarity: 3.12.3 sets it as
a floor above the scatter of the ladders the test actually
judges.  A purely local flatness test cannot both reject a
low-density false plateau and converge a genuinely slow
structure with a single threshold: a small high-symmetry cell
can sit on a plateau flat to a few hundredths of a meV/atom yet
a few tenths of a meV/atom above its true asymptote, while a
larger cell only truly flattens at a few tenths of a meV/atom.
Since a seed potential is a starting point that every
downstream run re-converges, the choice is to accept a
modestly loose plateau rather than complicate the detector --
a threshold near 1 meV/atom converges the well-behaved and
the slow structures alike, while a genuinely pathological
sweep (an energy swinging by tens of meV/atom across the
grid) is still correctly left unharvested for the curator.
Sharpening the detector to tell a false plateau from a true
asymptote is a deferred refinement, not a v1 concern.

**Metals stage no guidance entry (a held decision).**  A guidance
entry's entire content is the claim "for a structure like this,
this k-density is converged."  A metal cannot make that claim: its
energy does not converge in k-points at any mesh worth paying for,
which is the whole difficulty, and the climb acknowledges this by
short-circuiting at the first gapless rung and settling there as a
deliberately rough potential (3.12.3).  The rung it settles on is
a stopping point, not a converged density, and an entry recording
it as one would be read by the predictor as evidence -- worse,
disproportionate evidence, since a metal is often the only member
of its lattice family in a young collection and would then
dominate every prediction for that family through the distance
weighting of 7.6.  So the harvest skips it.  The producer's own
potential harvest (5.7) is unaffected: the rough potential is
still extracted, because a rough starting potential is exactly
what that database is for.

This is deliberately the small version of the fix.  The richer
one -- record the metal's gap and its lack of a density claim, and
let the predictor's first stage learn to expect a metal -- is the
natural upgrade once metals are more than an occasional member of
a seed collection, and is held until then (7.10).

**Promotion** (`guidance_promote.py`):

A curator helper.  Four modes of operation:

- *Interactive review* (default).  Lists all files
  under `staging/<system_type>/`; for each, prints the
  signature, measured quantities, verification grid,
  and provenance, and asks the curator to PROMOTE,
  SKIP, or DELETE.  Promoted files move to
  `entries/<system_type>/`.  Skipped files stay in
  staging for later review.  Deleted files are
  removed.  On a claim the collection already holds,
  REPLACE joins that list -- the one way an entry ever
  leaves `entries/`, and only by hand (below).
- *Auto-promote rule* (`--auto-promote`).  Promotes
  every staging file that satisfies an objective
  acceptance test:
  - The converged k-density landed in the middle 60%
    of the verification grid (not at either endpoint).
    A converged-at-endpoint result is suspicious: the
    grid may not have been wide enough.
  - The total-energy spread over the top three grid
    points -- their maximum minus their minimum, read
    from the entry's `grid_energies` array (7.2), which
    is why harvest records it -- is below
    `metric_threshold * 10` (the converged region is
    convincingly flat, not just one delta below
    threshold).  The spread is taken per atom and in
    eV, the basis `metric_threshold` is expressed in;
    because `grid_energies` are stored as raw
    total-cell hartree, the test divides by the cell's
    atom count and converts, the same normalization
    `pick_converged` applies (7.8 step c).  A spread --
    a like-for-like linear quantity -- is used rather
    than a variance, whose energy-squared units would
    not match the linear threshold.  A staging entry
    that lacks `grid_energies` (a hand-written manual
    entry) is never auto-promoted on this criterion; it
    falls to interactive review.
  - `gap_ev` and `gap_kind` are consistent
    (`gap_kind == "none"` iff `gap_ev == 0.0`).
  Files failing the rule stay in staging for the
  curator's review.  In practice this auto-promotes
  ~80% of seed-flight entries with the curator
  reviewing only the ~20% outliers.
- *Batch promote all* (`--all`).  Promotes every
  staging file without checking the rule.  Intended
  for one-off cases where the curator has manually
  reviewed the staging directory and decided the lot
  is good.
- *Dry run* (`--dry-run`).  Lists what would happen
  without moving files.

Promotion is a `mv` operation -- the file's contents
do not change.  This keeps provenance intact across
the staging boundary.

**One record per structure per settings.**  This is the whole
uniqueness rule the collection enforces, and promotion is where it
is applied, because promotion is the only stage that sees both the
incoming entry and the promoted corpus.  The harvest cannot: it
writes one entry per converged non-metal solid (metals stage none,
above) and has no view of what a curator accepted months ago.  Nor
does any existing guard
fire -- a re-run mints a fresh `entry_id` by construction,
since the slug hashes the flight id, structure, and
timestamp (7.5), so the collision refusal in `save_entry`
and `move_to_entries` never sees a collision at all.

Without this check a solid run ten times contributes ten
entries the predictor counts as ten independent
observations.  That is worse than untidy.  With
`neighbor_count = 5` (7.6), five copies of one calculation
fill the entire neighbor set; the weighted variance
collapses to zero and the confidence comes out near 1.0.
A single observation is then delivered as certainty, and
that certainty drives the flight builder into its
narrowest, least skeptical search (3.12.4) -- including
past the crystalline opening floor, which the confident
mode skips on the assumption that a real prediction sits
above it.

*What counts as the same claim.*  Two entries are re-runs
of one another when they agree on

    (system_type, basis, functional, kpoint_integration,
     basename(provenance.source_structure))

Each part earns its place.  The three settings fields are
already the predictor's sub-model partition (7.6): a
`gaussian` and a `gaussian-0.1` run of one solid are
different physics and must both survive.  The structure's
*basename* is used rather than its full path, because the
path records only where the structure cache happened to
sit and that location has moved (ARCHITECTURE 8.1); the
basename is `<reference_id>-<cell>.skl` (5.7), which is
the identity actually wanted, since a `full` and a `prim`
run of one COD entry are genuinely different structures.
The engine build is deliberately *not* in the key.  It is a fact
the report prints so a person can weigh it, not a dimension the
collection is partitioned along -- the same stance the run-reuse
cache takes (6.2.5): recorded, never compared.

*The check is an existence test, not a comparison.*  Either the
collection already holds a record for that claim or it does not.
That is the whole question, and it has two answers:

- **Occupied.**  The promoted entry stands untouched, the staged
  one is retired to `superseded/<system_type>/`, and both are
  printed side by side with their converged meshes, their dates,
  and the builds behind them.  Promotion therefore only ever
  *adds* to `entries/` on its own initiative; an entry a curator
  reviewed stays byte-identical for as long as it lives there.
- **Free.**  The ordinary path; the acceptance rule above decides
  it.

The `converged_mesh` (7.2) appears in that report and is never
tested by a branch.  An earlier draft compared the two meshes and
split the occupied case in three -- agreement retired the
newcomer, disagreement refused to act in any mode, and an absent
mesh (a manual entry, 7.9, carries no verification block and so no
mesh) counted as disagreement.  The split bought nothing.  The
tool's action was the same either way, because it has no verb for
retracting a reviewed entry; all the comparison decided was
whether the newcomer was archived or left in `staging/` to be
re-examined and re-reported by every later pass.  Archiving it
always is simpler and kinder to `staging/`, and the report still
shows the person both meshes -- which is where a disagreement was
going to be resolved in any case (VISION 16).

*The curator gets a verb.*  Interactive review therefore offers
**REPLACE** on an occupied claim, alongside PROMOTE / SKIP /
DELETE: retire the promoted entry to `superseded/` and promote the
newcomer in its place.  This is the only route by which anything
leaves `entries/`, and it exists solely behind a per-record prompt
a person answers by hand.  The standing objection -- a tool that
can retract a reviewed entry can do so by accident -- is an
argument against an *automatic* retraction and does not reach an
explicit one.  Without the verb, the report names a situation the
curator can act on only by moving files themselves, which is a
worse place to leave them than a prompt.

*Every mode applies the existence test, `--all` included.*
Refusing to store one claim twice is a correctness guard, not a
quality judgment, and `--all` waives only the latter -- it means
"I have reviewed these," not "store them however many times they
appear."  What the unattended modes do not offer is REPLACE:
`--all` and `--auto-promote` never retract.

*Within one staging batch nothing special happens*, which is the
point.  Two staged files can share a claim before either is
promoted -- the ordinary shape of a re-run harvested twice.
Promotion adds each entry to its in-memory index of the promoted
corpus as it promotes it, so the second file finds the claim
occupied and takes the branch above.  There is no separate
batch-resolution pass and no tie-break on `generated_at`; an
earlier draft carried both, and both were machinery for saving one
rule from having to apply twice in a row.  One consequence is
worth stating because it reverses an older contract: promotion no
longer judges each staged file in isolation, and `staging/` *is* a
uniqueness namespace.

*One documented property is amended.*  The acceptance rule
is specified above as reading the staged file alone, which
is why harvest records `grid_energies` at all.  The existence
test needs the promoted corpus as well.  The intent of
the original rule survives: what it avoided was depending
on the *flight workspace*, which is large, remote, and
reclaimable (6.2.12).  `entries/` is small, local, and
already the thing being written into.

**Why staging exists.**  An automated harvest is not
the same as scientific endorsement.  Bugs in the
harvest script, a verification grid that converged at
an unphysical artifact (e.g., near a numerical
instability), or a structure that was wrongly
classified by the curator could all produce entries
that should not propagate.  Staging gives the curator
a checkpoint to catch these before they influence
future predictions.  The friction is the point; the
`--auto-promote` rule lets the friction scale to a
500-entry seed flight without overwhelming the
curator.

### 7.9 Bootstrap and Day-1 Behavior

The dataspace starts empty.  Several things must work
gracefully in that state and as the dataspace fills.

**Empty-dataspace prediction.**  `predict(dataspace,
query, basis, functional, kpoint_integration)` over an
empty Dataspace returns
`PredictionResult(is_under_trained = True, ...)`.  The
flight-builder helper (DESIGN 6.2.8) then falls back to
the wide-grid default per 7.7 step 3.  The under-trained
path is unified with the no-sub-model and the sparse-
sub-model paths -- 7.6's step 2 fallback decides under
which conditions `is_under_trained` is set.

**Canonical entries for non-crystalline system_types.**
For amorphous, nanostructure, and molecular system_types,
the dataspace ships with one or two **canonical entries
seeded by hand** at day-1 (`source = "manual"`) so the
predictor can return something useful immediately.
Justification: k-density for these system_types is set
by the cell-volume convention to ~Gamma-only regardless
of chemistry, so a single canonical entry per non-
crystalline system_type captures essentially all the
information the predictor needs.  The canonical entries
are committed to git as part of the day-1 deliverable:

```
share/historicalGuidanceDB/entries/amorphous/
  amorphous-canonical.toml      # kpd = 25.0
share/historicalGuidanceDB/entries/nanostructure/
  nanostructure-canonical.toml  # kpd = 25.0
share/historicalGuidanceDB/entries/molecular/
  molecular-canonical.toml      # kpd = 1.0 (Gamma-only)
```

These manual entries carry empty composition_vector
(all zeros) and empty lattice_family -- the predictor
treats them as "the canonical answer for this
system_type."  Day-1 the predictor's non-crystalline
path simply returns the matching canonical entry
verbatim; future evidence may refine them.

**The wide-grid default** (crystalline only, used when
no usable predictor exists):

```
default_wide_kpoint_density_grid() = [
    25.0, 50.0, 100.0, 150.0, 200.0,
    250.0, 300.0, 400.0,
]
```

Eight points spanning a factor of 16, chosen to bracket
the k-density range commonly seen across published OLCAO
results.  This is deliberately broader than any "verify
around a prediction" grid: with no usable predictor, the
flight has to find the converged point unaided.  The
list lives in the flight-builder helper, not in the
dataspace itself (an empty crystalline subtree means no
dataspace content to consult).

**On the chicken-and-egg with Principle 11.**  Hardcoding
the wide-grid default is in mild tension with Principle
11 ("scripts must never silently encode 'experience' as
hardcoded constants").  The justification is that this
list is the *seed* the dataspace starts from, not a knob
meant to encode operating experience: by definition it
cannot itself live in the dataspace when the dataspace
is what is empty.  The list is documented here, kept
short and inspectable, and is *not* updated by
harvested entries.  Once the dataspace's crystalline
subtree holds even a few entries the predictor can use,
the wide-grid path becomes the rare under-trained
fallback rather than the dominant code path -- and the
moment a curator finds the bracket inadequate they can
edit the list explicitly with full audit trail.

**Non-convergence at the top of the grid -- failure mode
and recovery.**  Both the wide-grid default and the
verify-around-prediction grid can fail to converge within
their swept bounds.  Two distinct shapes:

- *Wide-grid default fails to converge.*  The candidate
  system requires a k-density above 400.0.  Diagnostic:
  the energy is still moving between the top two grid
  points.  The harvest hook (7.8) logs a warning, tags
  the flight with `prediction_mismatch = true`, and
  SKIPs the structure -- no entry is staged.  The user
  re-runs with a manually-extended grid (e.g., adds
  `kpoint_density = 600.0` to the flight builder's
  options), the second flight converges, and the
  staged entry then carries the higher value as the
  canonical k-density for that signature.  No automatic
  retry is built in for v1: silent re-dispatch with
  widened bounds would couple kaleidoscope to a notion
  of "still exploring" that Principles 9 and 12
  deliberately avoid.
- *Verify-around-prediction fails to converge.*  The
  prediction was wrong enough that the converged point
  fell outside the widened grid (7.7's
  `build_verification_grid`).  Same diagnostic: energy
  still moving at one end.  Same recovery: harvest hook
  SKIPs and tags `prediction_mismatch = true`.  Many
  such mismatches against the same neighbor set (the
  predictor keeps pointing at neighbors whose value is
  wrong for the queries that follow) signal that one or
  more of those neighbors is misleading and may warrant
  manual attention (delete-and-re-seed) on review.

The principle: **kaleidoscope dispatches once per
CalcUnit and reports.**  Multi-attempt convergence-finding
loops live in Python on the client side, not inside
kaleidoscope.  A researcher (or future tier-3 custom
wingbeat) wraps the dispatch in a re-run loop when needed;
the core stays single-shot.

**The seed flight (TODO C75).**  The first useful
entries for the crystalline subtree come from a
deliberate stratified seed run: ~150-250 calculations
spanning element-group pairs and common stoichiometry
patterns (binary AB, A2B, ABO3 perovskite, etc.) so the
predictor has broad chemistry coverage from day-2 on.
The seed flight uses the wide-grid default per
structure, runs through kaleidoscope, and feeds into the
`--auto-promote` rule of `guidance_promote.py` (7.8)
which lets the curator review only the ~20% outliers
rather than all 250.  Once the seed lands, the C48.3
producer (the first major consumer) sees confident
predictions for most reference solids it converges and
contributes its own entries on top.

**Manual seeding.**  Curators may write entries by hand
(`source = "manual"`) at any point: a researcher who has
converged TiO2-rutile in past work may seed it directly.
Manual entries are subject to the same schema validation
as harvested entries.  They are written directly to
`entries/<system_type>/`, skipping the staging step.

**Schema-version migration.**  When the schema bumps to
v2, `guidance_migrate.py` reads every v1 entry, applies
the v1 -> v2 transformation, and writes the v2 form back
in place.  Old entries are not discarded.  Parallel to
the DESIGN 5.7 build discipline.

### 7.10 Open Design Questions

- **Metals contribute nothing to the collection.**  The
  harvest skips them outright (7.8), which is correct as far
  as it goes -- a metal has no converged k-density to claim --
  but it also means that running one teaches the predictor
  nothing, and that a query for a metal draws its neighbors
  entirely from insulators whose densities are not the right
  answer for it.  The richer treatment records the metal's
  gap and its *absence* of a density claim, so the
  predictor's first stage learns to recognise the case and
  the second stage declines to predict rather than
  predicting badly.  That needs a v2 measured field and a
  predictor branch, and is deferred until metals are more
  than an occasional member of a seed collection.
- **Antiferromagnets are invisible to the spin feature.**
  The predictor's second character feature is the intensive
  magnetization `|M| / N_atoms` (7.6), built from the cell's
  net moment because that is what imago surfaces (the
  iteration file's column 6).  An antiferromagnet has zero
  net moment yet ordered local moments and can need careful
  BZ sampling, so this feature reads it as non-magnetic.
  Capturing it would require a sum of *absolute local*
  moments (`sum |m_i|`), which imago does not currently
  report per atom; a v2 schema could add an
  `abs_local_moment` measured field and switch the spin
  feature to it.  The verification grid is the day-1 safety
  net (an AFM whose density is mispredicted is still caught
  by the sweep), so this is a sharpness limitation, not a
  correctness one.
- **Metalloid assignment in `elemental_groups.toml`.**  The
  canonical metalloids (Si, B, Ge, As, Sb, Te) are also
  members of group_iv (Si, Ge), group_iii (B), pnictogen
  (As, Sb), and chalcogen (Te) by their column.  Each
  element must live in exactly one bucket (rule 4
  enforces the composition vector sums to 1.0 without
  double-counting).  Day-1 we keep them in their
  column-based groups and leave `metalloid` empty,
  documenting the choice in `elemental_groups.toml`'s
  comments.  Real seed-flight data will tell us
  whether metalloid-as-a-group meaningfully separates
  borderline-band semiconductors from their column
  neighbors; if yes, a curator moves them on a v2
  schema bump.
- **k-NN tuning knobs after the seed flight.**  All of
  `k_neighbors`, the distance weights (`w_comp`,
  `w_latt`, `w_gap`, `w_spin`), and the confidence
  normalizations (`sigma_gap_ref`, `sigma_kpd_ref`) are
  named constants in `guidance_db.py` (7.6).  Their
  defaults are educated guesses.  Calibration after the
  seed flight should pick values that minimize the
  predict-then-verify miss rate (how often does the
  verification grid land its converged point at an
  endpoint rather than the middle?).  The calibration
  procedure itself is open: either a one-shot post-seed
  analysis the curator runs, or a recurring auto-tune
  that watches harvest patterns over time.
- **Polytype confusion within the predictor.**
  Composition + lattice_family does not fully separate
  polymorphs that share a Bravais class (e.g.
  alpha-quartz SiO2 vs beta-quartz SiO2, both
  hexagonal).  The k-NN distance metric weights
  composition heavily; close polymorphs may smear into
  one another's predictions.  The verification grid is
  the safety net, but a recurring problem here would
  motivate adding a space-group integer as a side
  feature or extending lattice_family to a finer
  taxonomy.  Deferred until seed flight reveals
  whether this matters in practice.
- **Spin-polarization interpretation across magnetic
  orderings.**  `total_magnetization` is recorded per
  formula unit; an antiferromagnet sums to ~0 even
  though it has substantial local moments.  The
  predictor's stage-2 distance uses total
  magnetization, so AFM systems look like non-magnetic
  systems to it.  A future schema bump may add
  `local_moment_per_atom` (the max local moment
  Imago reports on any site) as a separate field.
- **Cell-size guidance for defect supercells.**
  Adding cell-size guidance to the schema is a
  meaningful future extension.  Most likely shape: a
  separate `defect_cell_size` measured quantity
  recorded when the calculation was a defect-in-host
  supercell (host_lattice and defect_species side
  features identifying the family).  Deferred until
  k-density predictor has proven out.
- **Multi-metric verification.**  Day-1's `metric` is
  `total_energy`.  Future flights may need forces
  or density-change.  The schema's `metric` field is
  registry-keyed (rule 10) so adding new metrics is a
  registry addition.  Open: how to record a *vector*
  of metric thresholds when multiple are required.
- **Functional / basis as sub-model dimensions vs k-NN
  features.**  The predictor conditions on
  (basis, functional) by running separate sub-models
  per pair.  Alternative: treat them as additional
  k-NN features.  The split approach is cleaner (no
  spurious cross-functional interpolation) but
  proliferates sub-models.  Calibration data may
  motivate revisiting.
- **Decay or staleness.**  Should entries from old
  imago commits carry a confidence penalty?  Argument
  for: imago itself evolves; settings that worked at
  commit X may not be optimal at commit Y.  Argument
  against: convergence settings are physical
  properties of the system, not of the code.  Deferred
  until real divergence is observed.

**Closed by decision (2026-05-28, user):** the historical-
guidance dataspace does *not* cross-reference into the
initial-potential database.  Once considered: a future
`pot_label` parameter that would let a guidance entry
say "for this family, use initial-potential DB entry
'default_solid'."  Decided against: the two artifacts
serve different audiences with different update cadences,
and entangling their schemas would couple their lifetimes.
Each artifact stands alone and shares only the curation
discipline (Principle 11), not its contents.

**Closed by decision (2026-05-29, user):** chemistry is
not used as a *signature axis* in v1.  Considered: a
categorical signature shape with `(system_type,
gap_type)` as a discrete partition and chemistry-Jaccard
as a soft refinement.  Replaced by the continuous
feature-space + k-NN design above because (a) gap is a
continuous variable, not a category, and binning it
costs prediction accuracy at the boundaries; (b) the
chemistry-to-electronic-character map is smooth enough
for k-NN regression to learn it, and (c) the resulting
predictor naturally produces a variance-based confidence
score that drives the verification-grid width.
Categorical lookup remains a fall-back option if the
regression approach fails to deliver, but is not the
primary day-1 design.

---

## References

H. J. Monkhorst, J. D. Pack, "Special points for
Brillouin-zone integrations," Phys. Rev. B 13, 5188
(1976). DOI: 10.1103/PhysRevB.13.5188
- The uniform reciprocal-space mesh, its offset (shift),
  and symmetry folding to the irreducible zone that
  section 3.7-3.10 specifies. The mesh `imago` builds is
  a Monkhorst-Pack mesh: a regular grid with equal base
  weights, Gamma-centered or half-shifted, reduced by the
  point group.

P. E. Bloechl, O. Jepsen, O. K. Andersen, "Improved
tetrahedron method for Brillouin-zone integrations,"
Phys. Rev. B 49, 16223 (1994). Key equations:
- Eqs. 14-16: analytic DOS formulas (total per
  tetrahedron)
- Eqs. 18-21: corner weights. Cumulative form
  `cornerIntgWt_LAT` from `bloechlCornerWeights` for
  integrated properties; energy derivatives
  `cornerDOSWt_LAT` from `bloechlCornerDOSWt` for
  energy-resolved DOS/PDOS
- Eq. 22: the curvature correction to the corner
  integration weights, `bloechlCornerCorrection` (1.3.1)
- Eqs. 23-24: a comparison between the true Fermi
  surface and the interpolated polyhedral one. Nothing
  in Imago computes them; they are listed so a reader
  meeting the range "22-24" elsewhere knows the
  correction is eq. 22 alone

A. K. Rappe, C. J. Casewit, K. S. Colwell, W. A.
Goddard III, W. M. Skiff, "UFF, a Full Periodic Table
Force Field for Molecular Mechanics and Molecular
Dynamics Simulations," J. Am. Chem. Soc. 1992, 114,
10024-10035.  DOI: 10.1021/ja00051a040
- Table 1: per-element parameters (r_i, Zstar_i,
  chi_i) used for bond stretching (section 4) and
  as inputs to the angle K heuristic (section 4.8.4)
- Eq. 2: natural bond length with electronegativity
  correction
- Eq. 3: bond stretching force constant formula
- Eq. 8: angle bending potential (cosine Fourier
  expansion -- not used directly; section 4.8.2
  explains why the harmonic approximation is preferred)
- Eq. 13: full UFF angle bending force constant K_IJK
  (not adopted; the geometric-mean heuristic in
  section 4.8.4 is used instead)

W. D. Cornell, P. Cieplak, C. I. Bayly, I. R. Gould,
K. M. Merz Jr., D. M. Ferguson, D. C. Spellmeyer,
T. Fox, J. W. Caldwell, P. A. Kollman, "A Second
Generation Force Field for the Simulation of Proteins,
Nucleic Acids, and Organic Molecules," J. Am. Chem.
Soc. 1995, 117, 5179-5197.  DOI: 10.1021/ja00124a002
- Referenced in section 4.8.4 for calibration context:
  typical harmonic angle force constants for organic
  molecules (C-C-C ~ 40, H-C-H ~ 35 kcal/mol/rad^2)

W. L. Jorgensen, D. S. Maxwell, J. Tirado-Rives,
"Development and Testing of the OPLS All-Atom Force
Field on Conformational Energetics and Properties of
Organic Liquids," J. Am. Chem. Soc. 1996, 118,
11225-11236.  DOI: 10.1021/ja9621760
- Referenced in section 4.8.4 for calibration context:
  independent confirmation that organic angle force
  constants fall in the 30-100 kcal/mol/rad^2 range

A. P. Thompson, H. M. Aktulga, R. Berger, et al.,
"LAMMPS - a flexible simulation tool for particle-based
materials modeling at the atomic, meso, and continuum
scales," Comp. Phys. Comm. 2022, 271, 108171.
DOI: 10.1016/j.cpc.2021.108171
- `angle_style harmonic`: E = K (theta - theta_0)^2
  convention used throughout section 4.8


## 8. Resource & Cost Guidance Dataspace

### 8.1 Overview and Motivation

This section pins the schema, data structures, and algorithms
for the resource-and-cost dataspace introduced in VISION Goal
6 and architected in ARCHITECTURE section 11.  Where the
historical-guidance dataspace (section 7) records what
operating point is *accurate*, this one records what a run
*costs*: for every imago run it stores the problem-size
signature, the parallel execution configuration, the
build/toolchain the binary was compiled with, and the
measured resources (peak memory, disk, walltime).  A
physics-informed regressor learns the cost surface, and the
near-term consumer turns a prediction into a SLURM resource
request that neither overflows memory nor exceeds the
walltime limit.

It is a deliberate sibling of section 7, not an extension
(ARCHITECTURE 11.1).  The two share the
library / producer / consumer discipline, the
staging-then-promote curation, schema versioning, and the
registry-validated-key discipline, but they are independent
artifacts.  The reason they stay separate is portability: a
converged k-density transfers across machines, whereas a
walltime is meaningful only on the machine that produced it.
This dataspace is therefore partitioned by a **hardware
fingerprint** (8.5), and its atomic unit is a single
**execution observation** -- one run under one configuration,
never collapsed to a per-system summary -- so the same
artifact serves provisioning now and configuration
optimization, build comparison, and scaling studies later
with no schema change.

### 8.2 TOML Schema (version 1)

Each observation is one TOML file under
`entries/<hardware_fingerprint>/` (promoted) or
`staging/<hardware_fingerprint>/` (harvested, awaiting a
curator).  It carries top-level keys plus four content blocks
(`[observation.signature]`, `[observation.execution]`,
`[observation.build]`, `[observation.resources]`) and a
`[observation.provenance]` block.

**Top-level keys (required):**

  Field                 Type    Description
  --------------------------------------------------------
  schema_version        int     Schema version integer. Must equal 1. Mirrors
                                the bare-integer SCHEMA_VERSION marker at the
                                dataspace root.
  observation_id        string  Unique slug identifying this observation;
                                equals the file stem and is unique across the
                                whole entries + staging tree.
  generated_at          string  ISO-8601 UTC timestamp of when the observation
                                was harvested or hand-entered.
  source                string  One of `flight` (harvested from a flight run)
                                or `manual` (hand-seeded, e.g. a bootstrap
                                point).
  outcome               string  One of `completed`, `oom`, `timeout`,
                                `failed`. Governs whether the resources block
                                is a measurement or a censored bound (rules
                                below).
  hardware_fingerprint  string  The partition key (8.5). Must be registered in
                                hardware_registry.toml and must equal the
                                entries/<fingerprint>/ directory the file
                                lives under.

**Size signature, under `[observation.signature]` (required):**

The cost-driving dimensions of the problem.  All are known
before the run, derived from the makeinput inputs and the
structure.  `secular_dimension` is the dominant scaling
variable the predictor regresses on.

  Field                    Type  Description
  --------------------------------------------------------
  atom_count               int   Number of atoms in the simulation cell. > 0.
  electron_count           int   Total electron count (all-electron: core +
                                 valence). > 0.
  valence_electron_count   int   SCF-active valence electron count. >= 0 and
                                 <= electron_count. The orthogonalized secular
                                 problem is built on the valence space.
  basis_function_count     int   Number of valence LCAO basis functions,
                                 before the spinor multiplier. > 0.
  wavefunction_components  int   1 for a 1-component (non-relativistic,
                                 Schrodinger) treatment, 4 for a 4-component
                                 (fully relativistic, Dirac) treatment. The
                                 4-spinor structure multiplies the secular
                                 dimension.
  secular_dimension        int   Dimension of the eigenproblem actually solved
                                 -- approximately basis_function_count times
                                 wavefunction_components after core
                                 orthogonalization. The dominant cost driver;
                                 > 0.
  kpoint_count             int   Number of k-points actually computed (IBZ or
                                 full mesh). > 0.
  spin_channels            int   1 for a spin-restricted run, 2 for
                                 spin-polarized.

**Execution configuration, under `[observation.execution]`
(required):**

How the run was launched.  This block is an extensible,
registry-validated key-value table: every key present must
appear in the checked-in `EXECUTION_KNOB_REGISTRY` (8.4), and
a new knob (a GPU count, a NUMA policy) is added by extending
the registry and bumping the schema version -- never by
silently introducing an unrecognized key.

  Field                 Type    Description
  --------------------------------------------------------
  node_count            int     Compute nodes the run used. > 0.
  cores_per_node        int     Physical cores per node engaged. > 0.
  total_cores           int     Total cores across all nodes. > 0. Recorded
                                explicitly, not derived, so a partially-packed
                                node is faithful.
  mpi_ranks             int     Number of MPI processes. > 0.
  omp_threads_per_rank  int     OpenMP threads per MPI rank. > 0; 1 when
                                OpenMP is unused.
  binding               string  Process/thread affinity policy.
                                Registry-validated: e.g. `none`, `core`,
                                `socket`.

**Build configuration, under `[observation.build]`
(required):**

The *coarse layer* of the two-layer build record (ARCHITECTURE
11.3): normalized, bucketed knobs that act as predictor
features.  Like the execution block this is registry-validated
against `BUILD_KNOB_REGISTRY` (8.4) and extensible; the
*fidelity layer* (the full compile string) lives in
provenance.  Values are bucketed on purpose -- an optimization
*level*, not a flag string; a *major* version, not a patch --
so a build is a comparable feature, not a fragmenting one.

  Field               Type    Description
  --------------------------------------------------------
  compiler_family     string  Fortran compiler family: e.g. `gfortran`,
                              `ifort`, `ifx`.
  compiler_version    string  Major (optionally minor) compiler version,
                              bucketed -- e.g. `13` or `13.2`, never a full
                              patch string.
  optimization_level  string  Optimization bucket: e.g. `O0`, `O2`, `O3`,
                              `Ofast`.
  arch_simd           string  Coarse instruction-set tag: e.g. `generic`,
                              `avx2`, `avx512`.
  blas_impl           string  BLAS/LAPACK implementation: e.g. `openblas`,
                              `mkl`, `reference`.
  blas_threading      string  `threaded` or `sequential`.
  scalapack           string  ScaLAPACK presence plus major version, or
                              `none`.
  hdf5                string  HDF5 variant plus major version: `parallel`
                              (MPI-IO) or `serial`.
  mpi_family          string  MPI implementation plus major version: e.g.
                              `openmpi-4`, `intelmpi-2021`, or `none`.

**Measured resources, under `[observation.resources]`
(required):**

What the run actually used.  The metric set is extensible via
`RESOURCE_METRIC_REGISTRY` (8.4).  When `outcome` is not
`completed`, the relevant metric is a *censored bound*, not a
point measurement (see the outcome rules below).

  Field              Type   Description
  --------------------------------------------------------
  peak_memory_bytes  int    Peak resident memory high-water mark. For
                            outcome=completed, the measured peak; for
                            outcome=oom, the memory limit the run hit (a lower
                            bound on the true need -- censored).
  disk_bytes         int    Disk footprint high-water mark (output plus
                            scratch). > 0 for a completed run.
  walltime_seconds   real   Wallclock runtime. For completed, the measured
                            time; for outcome=timeout, the walltime limit (a
                            censored upper bound on the need).
  cpu_seconds        real   Optional aggregate CPU time across ranks and
                            threads. May be absent.
  phase_timings      table  Optional sub-table of per-phase wallclock seconds:
                            setup, scf, eigensolve, postproc. May be absent or
                            partial.

**Provenance, under `[observation.provenance]` (required):**

Where the observation came from, plus the build fidelity
layer.  `compile_string` is always recorded so any flag is
recoverable post-hoc even when it is not a coarse knob.

  Field             Type    Description
  --------------------------------------------------------
  flight_id         string  The flight that produced this run. Non-empty for
                            source=flight.
  source_structure  string  The structure identifier or path the run computed.
  imago_commit      string  Git SHA of imago at run time. Non-empty for
                            source=flight.
  hostname          string  Host or cluster the run executed on (diagnostic;
                            the fingerprint is the canonical machine key).
  compile_string    string  The FULL, verbatim compiler invocation and flag
                            string -- the build fidelity layer (8.4). Always
                            recorded so any flag is recoverable post-hoc.
  library_detail    string  Verbatim detail of the linked libraries (exact
                            HDF5 / ScaLAPACK / BLAS / MPI versions and build
                            options). Free-form provenance.
  curator           string  Who or what produced the entry: e.g.
                            `resource_harvest.py`.

**Validation rules** (enforced at load time; every failure
names the file, block, and field at fault, as in section 7
and DESIGN 5.2):

1.  `schema_version` equals 1, in both the marker file and
    the entry.
2.  `observation_id` is unique across the whole entries +
    staging tree, and equals the file stem.
3.  `hardware_fingerprint` is registered in
    `hardware_registry.toml` AND equals the
    `entries/<fingerprint>/` directory the file lives under.
4.  `source` is one of `flight`, `manual`.
5.  `outcome` is one of `completed`, `oom`, `timeout`,
    `failed`.
6.  Size signature: `atom_count`, `electron_count`,
    `basis_function_count`, `secular_dimension`,
    `kpoint_count` are all > 0; `0 <= valence_electron_count
    <= electron_count`; `wavefunction_components` is 1 or 4;
    `spin_channels` is 1 or 2.
7.  Execution: every key is in `EXECUTION_KNOB_REGISTRY`;
    `node_count`, `cores_per_node`, `total_cores`,
    `mpi_ranks`, `omp_threads_per_rank` are all > 0;
    `binding` is a registered value.
8.  Build: every key is in `BUILD_KNOB_REGISTRY`;
    `compiler_family` and `optimization_level` are non-empty
    and registered.
9.  Resources: every key is in `RESOURCE_METRIC_REGISTRY`.
    For `outcome = completed`, `peak_memory_bytes`,
    `disk_bytes`, and `walltime_seconds` are present and > 0.
    For `outcome = oom`, `peak_memory_bytes` is present and
    interpreted as a lower bound; for `outcome = timeout`,
    `walltime_seconds` is present and interpreted as an upper
    bound (8.7).
10. Provenance: for `source = flight`, `flight_id`,
    `source_structure`, and `imago_commit` are non-empty;
    `compile_string` is present for every source.
11. Registry coupling: an unknown key in any registry-backed
    block is a hard error -- extensibility goes through the
    registry, never through silent key drift.
12. The schema is checked BEFORE the dataclass is built, so
    an omission surfaces as a clear validation failure rather
    than a constructor error.

### 8.3 Sketch (gold, single observation)

A completed run of a 12-atom cell, 1-component, on a 24-core
Haswell node with 4 MPI ranks x 6 OpenMP threads:

```toml
schema_version       = 1
observation_id       = "intel-haswell-24c-128gb-a1b2c3"
generated_at         = "2026-05-29T18:00:00Z"
source               = "flight"
outcome              = "completed"
hardware_fingerprint = "intel-haswell-24c-128gb"

[observation.signature]
atom_count              = 12
electron_count          = 312
valence_electron_count  = 96
basis_function_count    = 348
wavefunction_components = 1
secular_dimension       = 348
kpoint_count            = 84
spin_channels           = 1

[observation.execution]
node_count           = 1
cores_per_node       = 24
total_cores          = 24
mpi_ranks            = 4
omp_threads_per_rank = 6
binding              = "socket"

[observation.build]
compiler_family    = "ifort"
compiler_version   = "2021.5"
optimization_level = "O3"
arch_simd          = "avx2"
blas_impl          = "mkl"
blas_threading     = "threaded"
scalapack          = "mkl-2021"
hdf5               = "parallel-1.14"
mpi_family         = "intelmpi-2021"

[observation.resources]
peak_memory_bytes = 18253611008
disk_bytes        = 2147483648
walltime_seconds  = 4123.7
cpu_seconds       = 98969.0

[observation.resources.phase_timings]
setup      = 88.4
scf        = 3402.1
eigensolve = 2911.6
postproc   = 121.0

[observation.provenance]
flight_id        = "resource_seed_2026_05_29"
source_structure = "COD-1011098"
imago_commit     = "73eb567"
hostname         = "node042.cluster.umkc.edu"
compile_string   = "ifort -O3 -xCORE-AVX2 -qopenmp ..."
library_detail   = "HDF5 1.14.3 parallel (OpenMPI 4.1.5); MKL 2021.5"
curator          = "resource_harvest.py"
```

The emitter is the same hand-formatted, deterministic
discipline as section 7.5 (fixed block sequence, fixed key
order, `%.16e` for real values, byte-identical output for a
given in-memory observation); it is not restated here.

### 8.4 In-Memory Representation

The dataclasses mirror the schema block-for-block.  The
constants and the three registries are named in one place so a
post-seed recalibration or a new knob is a one-file change.

```
SCHEMA_VERSION          = 1
VALID_SOURCES           = ("flight", "manual")
VALID_OUTCOMES          = ("completed", "oom", "timeout",
                           "failed")
WAVEFUNCTION_COMPONENTS = (1, 4)     # Schrodinger | Dirac
SPIN_CHANNELS           = (1, 2)

# Extensible, checked-in registries.  A key not listed here
# is rejected at load (rule 11); a new knob/metric is added
# by extending the registry and bumping SCHEMA_VERSION.
EXECUTION_KNOB_REGISTRY = ("node_count", "cores_per_node",
    "total_cores", "mpi_ranks", "omp_threads_per_rank",
    "binding")
VALID_BINDINGS          = ("none", "core", "socket")
BUILD_KNOB_REGISTRY     = ("compiler_family",
    "compiler_version", "optimization_level", "arch_simd",
    "blas_impl", "blas_threading", "scalapack", "hdf5",
    "mpi_family")
RESOURCE_METRIC_REGISTRY = ("peak_memory_bytes", "disk_bytes",
    "walltime_seconds", "cpu_seconds", "phase_timings")
```

```
dataclass SizeSignature:
    atom_count              : int
    electron_count          : int
    valence_electron_count  : int
    basis_function_count    : int
    wavefunction_components : int    # 1 (Schrodinger) | 4
    secular_dimension       : int    # dominant cost driver
    kpoint_count            : int
    spin_channels           : int    # 1 | 2

dataclass ExecutionConfig:
    knobs : dict          # registry-validated key -> value;
                          #   node_count, mpi_ranks, binding...

dataclass BuildConfig:
    knobs : dict          # registry-validated coarse knobs;
                          #   the verbatim string is in
                          #   Provenance.compile_string

dataclass MeasuredResources:
    metrics : dict        # registry-validated metric -> value;
                          #   phase_timings is a nested dict
    censored : bool       # True when outcome != completed:
                          #   a bound, not a point measurement

dataclass Provenance:
    flight_id        : str
    source_structure : str
    imago_commit     : str
    hostname         : str
    compile_string   : str   # build fidelity layer (verbatim)
    library_detail   : str
    curator          : str

dataclass Observation:
    observation_id       : str
    generated_at         : str
    source               : str        # flight | manual
    outcome              : str        # completed | oom | ...
    hardware_fingerprint : str        # partition key
    signature            : SizeSignature
    execution            : ExecutionConfig
    build                : BuildConfig
    resources            : MeasuredResources
    provenance           : Provenance

dataclass ResourceDataspace:
    schema_version            : int
    observations_by_fingerprint : dict   # fp -> list[Obs]
    hardware_registry         : dict      # fp -> attributes
```

`ExecutionConfig` and `BuildConfig` hold open `dict`s rather
than fixed fields precisely so the registries -- not the
dataclass definition -- are the single source of truth for
which knobs exist.  Promoting a studied compiler flag to a
first-class feature (ARCHITECTURE 11.3) is then a registry
edit, not a dataclass change.

### 8.5 Hardware Fingerprint

The fingerprint is the coarse partition within which cost is
comparable.  The v1 recipe is a normalized slug

```
<cpu_vendor>-<cpu_microarch>-<cores_per_node>c-<mem_per_node_gb>gb
```

e.g. `intel-haswell-24c-128gb`.  The CPU string is normalized
to vendor + microarchitecture family (stepping, base clock,
and exact model number are dropped) so routine BIOS or
microcode churn does not fragment the data -- the granularity
tension flagged in ARCHITECTURE 11.8.  `hardware_registry.toml`
maps each fingerprint to its full probed attributes (exact CPU
model, socket count, memory, interconnect) for diagnostics;
the observation files carry only the fingerprint, never the
repeated attributes, mirroring how section 7 keeps the element
group table out of individual entries.

When a fingerprint is under-populated (below the predictor's
minimum sample count), the predictor falls back to the
nearest related fingerprint by probed attributes, or to a
conservative cold-start request (8.8); it never silently
predicts from one machine for another.

### 8.6 Predictor Algorithm

Within a fixed `(hardware_fingerprint, build-bucket)`, cost is
a smooth, physics-grounded function of size and parallel
configuration.  The model is therefore a **physics-informed
regression** rather than the pure k-NN of section 7.6.  Peak
memory scales roughly as the square of `secular_dimension` and
the eigensolve as its cube; the predictor fits a power law

```
log(resource) = log(A) + p * log(secular_dimension)
                + (parallel and spin correction terms)
```

by least squares per `(fingerprint, build-bucket)` group,
recovering the exponent `p` from the data (expected near 2 for
memory, near 3 for walltime) rather than assuming it.  The
parallel correction captures the speedup from `mpi_ranks` and
`omp_threads_per_rank` and the memory split across ranks.  A
k-NN fallback (over `secular_dimension`, `kpoint_count`,
`spin_channels`, scaled by the parallel config) is used when a
group is too thin to fit a stable exponent.  The exact
functional form and the thin-group threshold are tuning knobs
calibrated after the seed flight (8.8); they are deliberately
*not* required for the artifact to begin accumulating data.

For the near-term consumer -- **provisioning** -- the flight
layer queries the predictor with a proposed parallel config
and the new run's size signature, receives predicted memory /
disk / walltime, applies a safety margin, and emits the SLURM
request.  Because every observation also stores its full
parallel and build configuration, the same fitted surface
later answers *which* configuration or build is cheapest
(configuration optimization, build comparison) with no schema
change.

### 8.7 Capture and Harvest

Each observation is assembled at harvest from the four sources
of ARCHITECTURE 11.4: the dispatch-time size signature and
execution config (recorded by the wingbeat into the run
directory), the CMake-emitted `build_info.toml` (both build
layers), SLURM `sacct` accounting (`MaxRSS`, disk high-water,
`Elapsed`), and the optional imago self-report of per-phase
timings.  `resource_harvest.py` walks a finished flight,
builds one `Observation` per run directory, and writes it to
`staging/<fingerprint>/`; a curator promotes with the same
discipline as section 7.8.

**Censored (non-completed) runs are retained, not discarded.**
A run killed for OOM is positive evidence that its config is
insufficient at that size: it is staged with `outcome = oom`,
`peak_memory_bytes` set to the memory limit it hit, and
`MeasuredResources.censored = True`.  A `timeout` run is staged
with `walltime_seconds` set to the limit.  The regressor (8.6)
treats a censored memory observation as a lower bound and a
censored walltime as an upper bound rather than a point; how
that censoring enters the least-squares fit is an open
question (8.9).  A `failed` run (a Fortran abort unrelated to
resources) carries no usable cost signal and is staged only
for diagnostics, never promoted.

### 8.8 Bootstrap and Day-1 Behavior

A fresh fingerprint has no observations, so the predictor
cannot yet predict for it.  Day-1 behavior on an empty or
under-populated fingerprint: the provisioner falls back to a
conservative resource request (a generous memory and walltime
ceiling, optionally scaled from a related fingerprint by
probed attributes), runs the job, and harvests the result --
which seeds the fingerprint.  A small `manual` seed (a handful
of hand-entered observations spanning the size range on the
local machine) accelerates this, exactly as the section-7 seed
flight (C75) bootstraps convergence guidance.  The artifact
then improves monotonically: each completed flight appends
observations, and the fitted exponents tighten as evidence
accumulates.

### 8.9 Open Design Questions

- **Exact regression form and censored-data handling.**  The
  power-law-in-`secular_dimension` model (8.6) and how OOM /
  timeout bounds enter the fit (a censored / Tobit-style
  regression, or simply weighting them as bounds) are open and
  calibrated after the seed.
- **`secular_dimension` provenance.**  Whether it is recorded
  directly from imago (authoritative) or derived from
  `basis_function_count x wavefunction_components` minus the
  core-orthogonalization reduction (portable but approximate)
  is open; the schema records it directly, with the primitives
  kept for cross-checking.
- **Aggregate vs per-rank memory.**  `peak_memory_bytes` as a
  job aggregate vs per-node vs per-rank changes how the
  parallel correction is modeled; ARCHITECTURE 11.8 flags the
  reconciliation of `sacct` / `time` / self-report sources.
- **Build effects on numerics, not just cost.**  Per
  ARCHITECTURE 11.8, build choices can perturb low-order
  digits of the physics result; whether the build block is
  ever referenced from the section-7 (convergence) side --
  against the no-cross-reference boundary -- is to be settled
  here, not assumed.

## 9. Parallel Decomposition

### 9.1 Overview and Motivation

This section pins down the data-distribution algorithms for
parallelizing a single imago calculation across MPI ranks --
the intra-problem axis of VISION Goal 7, architected in
ARCHITECTURE 6.5-6.8. Two decisions that earlier drafts left
open are now made and recorded here: the distributed
eigensolver is ELPA (ARCHITECTURE 6.6, ruled 2026-08-18 with
the toolchain verified), and the order of implementation
follows the measured cost structure of ARCHITECTURE 6.8
rather than the sibling branch's inherited ordering -- the
three-centre electronic-potential term distribution first
(9.5), the eigensolver boundary second (9.6), the grid
balancer (9.2) LAST, since the loops it serves are 1-3 % of
a medium run. The block-cyclic scheme (9.3) and the
redundant-pair work assignment (9.4) are kept: the solve
requires the first, and the second becomes the pair-level
refinement of the integral stage when rank counts exceed the
term count. Remaining open choices are collected in 9.8.

Two measured inputs anchor this section (dev/PERFORMANCE.md
"Baseline", "Coarse time map", "PA1 cost distributions",
2026-08-18): (a) the three-centre stage and the secular
solve are 72-86 % of every run above toy size; (b) within
the three-centre stage, per-term costs vary only mildly
(max/mean 1.5-1.6, monotone in the alpha exponent), per-row
pair costs likewise (max/mean about 2), and the atom-pair
loop is only 16 % of the stage on a complex multi-k-point
cell -- the balance is the per-term core-orthogonalization
and write, which a term distribution carries with it for
free and a pair distribution would leave serial.

### 9.2 One-Dimensional Grid Load Balance

The real-space site loops in electrostatics and exchange-
correlation are independent per site, so they parallelize by
handing each rank a contiguous range of site indices and
reducing the partial results. Given a quantity `toBalance`
to divide among `mpiSize` ranks:

```
jobsPer   = toBalance / mpiSize       (integer divide)
remainder = mod(toBalance, mpiSize)
```

Each rank receives `jobsPer` sites, and the highest
`remainder` ranks each take one additional site so that no
work is dropped when the division is uneven. The rank then
loops over its `[initialIdx, finalIdx]` range and the
partials are combined with `MPI_REDUCE` under `MPI_SUM`.
This is the `loadBalMPI` algorithm from the sibling branch
and the lowest-risk parallelism imago can adopt. An earlier
draft called it the recommended first increment; the
2026-08-18 measurement demoted it to LAST (ARCHITECTURE 6.8:
the loops it serves are 1-3 % of a medium run). The
algorithm stands unchanged for when its turn comes, and it
doubles as the deal used for the term distribution of 9.5.

### 9.3 Block-Cyclic Matrix Distribution

The interaction-integral and Hamiltonian matrices are
distributed across a two-dimensional process grid in a
block-cyclic pattern -- the layout that distributed dense
linear algebra (ScaLAPACK, ELPA) requires for load balance.
The matrix is tiled into equal blocks, and the blocks are
dealt out to ranks cyclically in both dimensions so that, as
an elimination front sweeps the matrix, every rank stays
busy rather than idling once its corner is consumed. A naive
contiguous split would leave all but one rank idle near the
end of a factorization, which is why the cyclic deal is not
optional.

The process grid is chosen as close to square as the rank
count allows -- a perfect square when possible, otherwise
the most balanced integer factorization of the rank count --
because square grids minimize communication volume in the
factorization. Each rank allocates only its local portion of
the matrix and maintains a descriptor that maps local
(row, col) indices back to global matrix indices and forward
again. The `MatrixDescriptor` type and the most-square grid
helper from the sibling branch's `mpi.f90` are the concrete
starting point.

### 9.4 Work Assignment: Redundant Atom Pairs

A subtlety distinguishes *computing* the matrix from
*distributing* it. Each matrix block draws contributions
from atom-pair orbital interactions, and an atom pair's sub-
matrix generally will not align with block boundaries. Three
strategies were weighed in the sibling branch:

1. Each rank computes only the elements it owns. Rejected:
   an atom pair straddling a block boundary forces different
   ranks to compute different orbital-orbital interactions
   of the *same* pair, demanding intricate partial-
   computation logic.
2. Distribute atom pairs once, then communicate the stray
   elements each rank computed but does not own. Rejected:
   each pair is computed once, but the communication
   bookkeeping is again costly and error-prone.
3. **Adopted:** distribute atom pairs so each rank computes
   every element its own blocks need, accepting that a few
   atom pairs are computed by more than one rank. Each rank
   keeps the elements that fall in its blocks and discards
   the rest.

Strategy 3 trades a little redundant arithmetic for *no
communication and simple logic* during assembly -- the right
trade when integral evaluation is cheap relative to
interconnect cost. This decision is inherited from upolcao's
design and is the recommended starting point; it is revisited
only if profiling shows the redundant computation dominates.
Its place in the sequence: it serves the PAIR-level
refinement of the integral stage (9.5) and the assembly of
distributed H and S for the solve (9.6), not the first
implementation, which distributes whole terms.

### 9.5 Three-Centre Term Distribution

The three-centre electronic-potential integral stage is
distributed BY POTENTIAL TERM. This is the decided first
parallel increment (ARCHITECTURE 6.5/6.8), and the decision
rests on the PA1 measurement rather than symmetry with the
other decompositions: per-term costs are mild and
predictable (max/mean 1.5-1.6, monotone in the alpha
exponent), and 27-84 % of the stage -- the larger figure on
multi-k-point cells -- is per-term work OUTSIDE the atom-pair
loop (the core-orthogonalization and dataset write), which a
term distribution carries to the ranks for free and a pair
distribution would leave serial.

**The serial stage as it exists (seam inventory).** The
stage is `elecPotGaussOverlap(packedVVDims, did, aid)` in
`integrals.F90`, called ONCE from the SCF setup path of
`imago.F90` (after `makeAlphaPotDist`). Its inputs, who
makes them, and when:

- `potDim` (module `O_Potential`), `potTypes` (`O_PotTypes`),
  `potSites` (`O_PotSites`): loaded unconditionally by the
  input-parsing pass before any integral work. A "term" is
  one (potential type, alpha) pair; the term index
  `currentIterCount` runs 1..`potDim` in the doubly nested
  loop over potential sites carrying `firstPotType == 1` and
  their types' alphas. The term's identity is posted in the
  module variables `currPotAlpha`, `currPotNumber`,
  `currPotElement`, `currAlphaNumber`, `currMultiplicity`
  before the per-term worker runs.
- `did(numKPoints, potDim)` = `atomPotOverlap_did` and
  `aid(potDim)` = `atomPotTermOL_aid`: HDF5 dataset and
  completion-attribute handles, allocated and created by
  `initSCFIntegralHDF5` (`hdf5SCFIntg.F90`) inside the
  k-point group during HDF5 setup -- one dataset per (k-point,
  term), one attribute per term.
- The per-term worker `gaussOverlapEP` assembles the full
  valence-valence (plus core-valence, core-core) matrices for
  its term over all atom pairs and lattice cells, then calls
  `ortho(4, packedVVDims, did(:,term), aid(term))`, which
  core-orthogonalizes and writes the term's datasets and
  finally sets the term's completion attribute.
- Restart already exists at term granularity: the term loop
  reads `aid(term)` first and skips completed terms. The
  parallel design inherits this unchanged -- the completion
  attribute is the unit of both restart and distribution.

**The distribution** (deal corrected 2026-08-20: drafting the
pseudocode showed the earlier sentence contradicted itself --
`loadBalMPI` deals CONTIGUOUS blocks, so applying it to a
cost-sorted list would give the first rank every costly term
and the last rank every cheap one, the opposite of the
largest-first balance it claimed). The deal is a SNAKE over
the cost order: sort the `potDim` terms most-diffuse-first
(ascending alpha exponent, the cost proxy the measurement
licensed) and deal them to ranks in boustrophedon rounds
(0,1,..,N-1 then N-1,..,1,0 and repeat), which equalizes rank
loads under any cost monotone in the exponent. ASSIGNMENT
order and ITERATION order are then deliberately different:
each rank walks the terms it was dealt in ORIGINAL term-index
order (type-major, most-diffuse-first within a type), because
that order is what the `anyElecPotInteraction` inheritance
below assumes. Each rank runs the EXISTING term loop body for
its terms: same worker, same `ortho` compute, same datasets.
No partial matrices ever cross ranks. Completion attributes
make the result identical to a serial run's file, term by
term.

**HDF5 discipline (first version; mechanism decided
2026-08-20).** The datasets are disjoint per term, so ranks
never write the same object, but the serial HDF5 library
requires that a FILE not be open for writing by several
processes at once -- and enforces it with its own file lock
(the PA2 acceptance run demonstrated the lock working on this
filesystem: worker opens fail with errno 11 while a writer
holds the file). That lock IS the mutex. Per term, a rank
computes and packs the full k-point set of matrices into a
buffer with the file CLOSED, then opens the file (retrying,
with HDF5's error stack silenced, until the lock is granted),
opens the term's datasets and attribute BY NAME (the same
deterministic names `initSCFIntegralHDF5` created and
`accessSCFIntegralHDF5` already reconstructs for restart),
writes, sets the completion attribute, and closes. Only the
write time is ever serialized, no MPI token protocol exists,
and a one-rank run degenerates to open-write-close with no
contention. For the discipline to work the file must be open
by NOBODY between writes, which sets the stage's file
lifecycle: root (which alone ran the setup writes, see the
run shape below) CLOSES the file before the term loop, every
rank writes its terms under the lock, and root REOPENS it
through the existing restart access path after the stage-end
barrier -- the barrier is required, since an early root
reopen would hold the lock against every rank still writing.
Restart keeps its meaning: root reads the `potDim` completion
attributes before closing and broadcasts the done-mask, and
the snake deal runs over the undone terms only. How large the
serialized write time is has NOT been measured: PA1 stamped
whole `gaussOverlapEP` calls and the pair-loop rows, so what
it bounds is the combined ortho-plus-write remainder of the
stage (27 % on the real/Gamma glass, 84 % on the complex
multi-k cell), and the write-only fraction inside that
remainder is unknown. PA3 therefore stamps the writes and the
lock waits per rank and validates this discipline against
those numbers rather than assuming it; if the write dominates
the remainder on multi-k-point cells, the collective form of
9.7 moves up the queue. The collective parallel-HDF5 form is
otherwise 9.7's calibration, not a prerequisite.

**Run shape at this increment.** Everything that touches the
HDF5 file outside the term stage stays on root. The
replicated, file-free setup (input parsing, lattice,
k-points, basis, alpha distance tables) runs on every rank;
`initHDF5_SCF` and the cheap integral stages that write
through its held handles (overlap, kinetic energy, mass
velocity, nuclear potential -- together a few percent of a
run) are root-only; the term stage is the one distributed
region; and after its barrier the workers have nothing left
-- they skip to the end of `subroutine Imago` and park at the
PA2 certificate barrier while root alone runs the remaining
setup and the whole SCF iteration. Their cores idle there:
that is the honest, accepted cost of distributing one stage
at a time, it is bounded by the stage's share of the run
(ARCHITECTURE 6.8), and PA4 -- the solve, the largest share
at scale -- is what retires it. The acceptance metric for
PA3 is therefore the STAGE stamp falling with rank count,
not the whole-run wall time.

**One cross-term coupling to break.** The serial loop shares
`anyElecPotInteraction` -- the (pair, cell) negligibility
bitmask -- across the alphas of one type: each alpha inherits
the bits its predecessors cleared and prunes its work with
them. The mask is an optimization at the NEGLIGIBILITY FLOOR
(a bit is cleared when a pair/cell contributed nothing above
the threshold -- which is below-threshold, not exactly zero),
so each rank simply rebuilds its own mask for the types it
touches; a rank that holds a type's tighter alphas without
its diffuse ones loses some pruning and no correctness. The
inheritance also survives ownership GAPS, which is why the
iteration order above matters: a bit cleared by a more
diffuse alpha is validly cleared for every tighter alpha of
that type, whether or not the rank owns the alphas in between
-- so a rank walking its owned terms in original order keeps
nearly all of the serial pruning, resetting the mask only
when it crosses into a new type. This is the only coupling
between terms; everything else a term reads is constant input
data.

One measured consequence (2026-08-20, the first multi-rank
acceptance run): because a rank's mask lacks the clearings
that unowned alphas of its types would have made, it COMPUTES
sub-threshold contributions the serial walk skipped -- so a
multi-rank run's term matrices agree with serial to the
negligibility floor (relative 1e-13 observed; slightly MORE
arithmetic than serial, never less), not to the bit. Serial
and the one-rank parallel run keep the full inheritance and
remain bit-identical. The acceptance criterion follows the
physics, and bit-exactness is meaningful only WITHIN one
build: the serial and MPI builds link different OpenBLAS
builds whose results differ in the last bit, an effect
independent of this design (measured 2026-08-20: about 1e-15
absolute, plus the eigenvector gauge it triggers). So: at one
rank, h5diff must be exactly clean against the SAME build's
serial-shaped run (the bare singleton and `mpirun -np 1` must
be bit-identical to each other, and the serial binary must be
bit-identical to the pre-change serial binary); above one
rank, the comparison basis is the same build's one-rank file,
and the term data must agree within a tolerance safely above
the floor and below anything physical (1e-10 relative), with
the derived outputs (iteration trace, energies)
digit-identical against the recorded baselines and the
eigenVECTOR datasets excluded from the file comparison --
floor-level jitter in H rotates degenerate eigenvectors
arbitrarily (the standard gauge freedom) while every
invariant built from them is checked through the
digit-identical outputs.

**Pair-level refinement (later).** When ranks exceed
`potDim`, or a single term must go faster, the atom-pair
loop within a term splits by the packed pair index under 9.4
(each rank computes the pairs its blocks need, redundantly
at block edges) -- or, cheaper and first, OpenMP threads
inside the rank cover the pair loop while ranks stay at term
granularity. The PA1 row measurement (max/mean about 2, flat
in the row index) says a contiguous pair-range split is
adequate when that day comes.

### 9.6 The Eigensolver Boundary

The secular solve sits behind ONE call site with swappable
backends (ARCHITECTURE 6.6): serial LAPACK (today's
`solveZHEGV` / `solveDSYGV`) for one rank, ELPA (decided
2026-08-18) for the distributed case. This subsection names
the seams the boundary must respect.

**The serial solve as it exists (seam inventory).**
`secularEqnSCF(spinDirection, numStates)` in
`secularEqn.F90`, called per spin from the SCF iteration
loop of `imago.F90`:

- It loops k-points INTERNALLY (`do i = 1, numKPoints`),
  and per k-point: reads the packed nuclear-potential matrix
  (`readPackedMatrix` from `atomNPOverlap_did(i)`), then
  accumulates kinetic energy, optionally mass-velocity, and
  the `potDim` potential terms weighted by
  `potCoeffs(j, spinDirection)` (`readPackedMatrixAccum`
  from `atomPotOverlap_did(i,j)`) -- so H is assembled FROM
  THE HDF5 FILE each iteration, with the iteration's
  potential coefficients; the integral datasets of 9.5 are
  its direct input.
- The packed H is unpacked into the module array
  `valeVale(valeDim, valeDim, spin)` (complex) or
  `valeValeGamma` (real), S likewise into `valeValeOL`;
  `solveZHEGV(valeDim, numStates, H, S, eigvals)` (or
  `solveDSYGV`) destroys its inputs and returns the lowest
  `numStates` eigenpairs; eigenvectors and eigenvalues are
  written to the HDF5 eigen datasets
  (`eigenVectors_did(i, spin)`, per-k-point completion
  attributes -- the existing per-k-point restart unit).

**The boundary.** One routine -- working name
`solveSecular(H_desc, S_desc, numStates, eigvals, eigvecs)`
-- owns everything between "packed H and S exist" and
"eigenpairs exist where the writer expects them", including
any redistribution. Backends:

- **Serial (one rank per problem):** exactly today's path;
  OpenBLAS threads inside the rank are the intra-node lever.
- **ELPA (one problem across ranks):** the block-cyclic
  layout of 9.3; the generalized problem is reduced via
  ELPA's Cholesky path and solved by its two-stage
  tridiagonalization; ScaLAPACK supplies the layout
  utilities beneath it. The redistribution of H and S INTO
  the block-cyclic layout (from rank-0 assembly at first;
  from distributed assembly once 9.5's pair refinement
  exists) and of the eigenvectors back OUT lives inside the
  boundary, so its cost is measured with the solve's.

**The outer, communication-free level.** k-points and spins
are independent solves consuming the same integral file.
When `numKPoints * spin >= ranks`, whole solves are dealt to
ranks (each running the serial backend) before any single
solve is distributed -- for the multi-k complex runs this is
the entire win at small rank counts. The per-k-point
completion attributes already make that deal restartable.
The choice between "many serial solves" and "one distributed
solve" is by problem: it is the per-kernel device-placement
principle (VISION 14) applied one level up.

### 9.7 Parallel HDF5 Alignment

The distributed matrices are written to and read from HDF5
collectively. For compression and write efficiency the on-
disk chunk size should align with the block-cyclic block
size, so each rank's write touches whole chunks rather than
splitting them. The exact collective-write pattern against
compressed chunks needs measurement -- one chunk per block
per rank, versus larger chunks filled by several ranks'
collective contributions -- and is flagged as an open
calibration in 9.8. It is NOT a prerequisite for 9.5, whose
first version writes disjoint per-term datasets in turn
through serial HDF5.

### 9.8 Open Design Questions

- ~~**Distributed eigensolver backend.**~~ **DECIDED
  2026-08-18: ELPA** (ARCHITECTURE 6.6; toolchain built and
  the API handshake verified on the cluster). ScaLAPACK
  remains beneath it for layout utilities, not as the
  solver.
- **Device-placement expression in Fortran.** How the per-
  kernel CPU/GPU boundary (VISION Principle 14) is expressed
  -- OpenACC, OpenMP target, CUDA Fortran, or a library
  boundary such as ELPA -- is open and will likely differ
  per kernel.
- **Parallel-HDF5 chunk/block strategy.** The collective-
  write pattern against compressed chunks (9.7) is settled
  by measurement, not assumed.
- **Implementation order** -- SETTLED 2026-08-18 by the
  measurement (ARCHITECTURE 6.8): term distribution (9.5),
  then the eigensolver boundary (9.6), then valence charge
  density on the solve's distribution, then grid work (9.2)
  last. The validation gate between stages is unchanged:
  each stage must reproduce the serial benchmark results
  (energies to print precision, HDF5 content by h5diff)
  before the next begins, and replicate-and-broadcast forms
  are never counted as progress.

---

## 10. Runtime Citation Banner

### 10.1 Why the Program Prints Its Own Citation

VISION principle 15 holds that Imago must be built so its credit
survives being passed on. The licensing machinery of
ARCHITECTURE 1.1 does part of that: `NOTICE` travels with any
derivative work, `CITATION.cff` is machine-readable, and a DOI
gives a reference list something to point at. All of it is
addressed at someone who is redistributing or packaging the
code.

None of it reaches the person writing the paper. That person ran
a calculation, has results in front of them, and is composing a
methods section. They will cite what they can see. Printing the
citation into the output they are already reading is the only
mechanism in this project that reaches attribution at the moment
it is actually decided, which is why LAMMPS, VASP, and Quantum
ESPRESSO all do it. Measured against every other step in the
attribution chain, this is the one with the highest expected
effect.

It also serves a second purpose that has nothing to do with
credit. A log that names the code, its version, and the methods
a run exercised is a provenance record. Six months later it
answers "what produced this number" without recourse to memory.

### 10.2 Two Blocks, Not One

The banner is split in two, and the split follows from when
information becomes available.

The **identity block** prints at startup: the butterfly, the
wordmark, the version, and the citation for Imago itself. All of
it is known before any work begins.

The **methods block** prints at the end of the run: the
references for the specific methods that run actually exercised.
This cannot go at the top, because at startup the program does
not yet know whether the tetrahedron integration will be used or
whether UFF parameters will be consulted. Attempting to print it
early would mean either listing everything Imago could
conceivably do -- which trains the reader to skip it -- or
guessing from the job code, which would be wrong whenever a
branch is not taken.

Splitting also matches how the reader works. The identity block
answers "what is this" on opening the log. The methods block
answers "what do I cite" when the run is finished and the
results are in hand.

### 10.3 The Identity Block

The block is the logo rendered as text: the monarch above, the
`imago` wordmark below, then the version and citation lines.

Three properties are fixed.

**Width is 51 columns.** This is not arbitrary. `O_TimeStamps`
declares its operation labels as `character(len=51)` and prints
its rules at that width, so the banner sits directly above
output already committed to that column. A wider banner would
have to widen every `opLabel` with it.

**The art is literal text, not generated.** It was produced by
sampling the project logo and then hand-kerned -- the letter
spacing between `i` and `m` differs from the spacing between the
remaining letters, because uniform tracking over-separates a
one-column letter. No script reproduces that adjustment. The
generator was a starting point; the checked-in text is the
source of truth, and it is edited directly if it is edited at
all.

**It lives in `src/data/banner.txt`**, installed to `share`
alongside `elements.dat` and the rest, and located at run time
through the `IMAGO_DATA` environment variable. This is the
mechanism the engine already uses -- `elementData.f90` and
`potential.f90` both resolve their data files exactly this way
-- so the banner introduces no new runtime dependency and no
new failure mode. If `IMAGO_DATA` is wrong, `elementData` fails
first and the run never reaches the banner, which means the
file's availability is already guaranteed by the same condition
that makes the run work at all.

Keeping the art in a data file rather than compiling it in has
a second benefit. `imago.py` is what a user actually invokes,
and its output is more visible than the engine's log. A shared
file lets the driver print the same banner without duplicating
the artwork in a second language.

**A parser must not "fix" the whitespace.** Trailing spaces are
insignificant, but leading spaces carry the kerning and the
centring. A reader that strips or normalizes leading whitespace
destroys the alignment, and a well-meaning editor configured to
trim whitespace on save will do the same to the file itself.

### 10.4 The Citation Text

The identity block closes with the citation, whose fields are
the same ones `CITATION.cff` carries: title, author, version,
DOI, and repository URL.

They live in `src/data/banner.txt` with the artwork, below it
and separated by a blank line, rather than being compiled into
the Fortran.

An earlier draft of this section argued the opposite -- that
the version and DOI change only at release, which is already a
rebuild, so a data file would add a failure mode for no gain.
That reasoning rested on a false premise. Reading an installed
data file is not a new mechanism here and carries no new failure
mode, for the reasons given in 10.3, so the cost side of that
argument was close to zero and the conclusion did not follow.

Keeping the citation beside the artwork gives one file a human
edits when the DOI of TODO A12 arrives, instead of a Fortran
literal that must be found and recompiled. It also avoids
duplicating the fields in a second language when `imago.py`
prints the same block.

`CITATION.cff` remains the authoritative record; `banner.txt`
restates it for display, and the two must be updated together.
That duplication is real but it is now between two text files a
maintainer edits in one sitting, rather than between a text file
and compiled source.

Until the DOI of TODO A12 exists, the citation names the
repository and states that a DOI is pending. It must not invent
one.

### 10.5 Methods Actually Exercised

The `## References` section of this document already carries the
citations that Imago's methods rest on: Monkhorst and Pack for
the reciprocal-space mesh, Bloechl and co-workers for the
tetrahedron integration, and Rappe and co-workers for the UFF
parameters. These are the references a methods section needs,
and today a user has to know to go looking for them.

The methods block closes that gap. It is a registry pairing each
reference with a predicate answering "did this run use it," and
at the end of the run it prints only those whose predicate is
true. A run using a Monkhorst-Pack mesh and Gaussian broadening
prints one reference; the same run with tetrahedron integration
prints two.

The registry the engine carries is narrower than the References
list, and must be. UFF appears nowhere in the Fortran: those
parameters belong to the force-field path of section 4.8, which
is Python, as do the Cornell, Jorgensen, and LAMMPS references
beside them. An engine-side entry for any of them could never
be selected. If that path should announce its own citations, it
announces them from `make_reactions.py`.

Keeping the registry beside the References section matters:
adding a method to Imago and adding its citation become one
task, and a reference that no predicate can ever select is
visibly dead. The predicates read state the engine already
holds, such as `kPointIntgCode` distinguishing Gaussian from
tetrahedron integration, so nothing new has to be tracked to
support this.

### 10.6 Where the Banner Goes

**The log (unit 20), and nowhere else.** The identity block
prints immediately after the log is opened in
`parseCommandLine` and before the first operation timestamp;
the methods block prints at the end of the run. The slot at the
head is forced -- the unit does not exist earlier.

**The log is closed once, at the end, after the methods block.**
It is not today. `cleanUpSCF` closes unit 20 when no post-SCF
stage follows, `cleanUpPSCF` closes it unconditionally, and
`loen` closes it before signalling completion, so by the time
the run ends the log is already shut. Writing to a closed unit
does not fail: Fortran reconnects it to `fort.20` and truncates,
which destroys the whole log and leaves behind only the
citation that destroyed it. This was found by running the code,
not by reading it.

The correction is to close the log where it is opened -- once,
at the outermost level, after everything that writes to it has
finished. A file opened in `parseCommandLine` and closed in
three separate cleanup routines is fragile whatever is appended
to the run, and the alternative of having the methods block
reopen the file would leave that fragility in place while
hiding its next symptom.

The block still asks whether the unit is open before writing. If
it is not, it says so on standard output and prints nothing to
the log. That is deliberately not a repair: a later change that
re-introduces an early close should cost the citations and
announce itself, rather than silently costing the run's entire
output the way it did the first time.

**The `fort.2` completion signal moves with it, and for the same
reason.** Section 6.1.2 says `fort.2` certifies that the binary
ran without an abortive error, and `imago.py` treats its
presence as the sole success gate. The engine does not keep that
promise. The file is created at the end of `cleanUpSCF`,
`cleanUpPSCF`, and `loen` -- the same three inner routines -- and
every post-SCF job runs both stages in one invocation, so
`cleanUpSCF` certifies success and the entire post-SCF stage
then runs afterward. A `stop` during that stage leaves the
certificate already on disk, and since `STOP` and
`STOP 'message'` both exit zero under gfortran the return code
does not catch it either. The driver reports success for a run
that died.

Both are the same defect: an event that belongs to the whole run
emitted from a routine that only knows about one stage of it. So
the tail of the run is ordered once, at the outermost level,
where each step means what it says:

1. the methods block -- the last write to the log
2. `close (20)` -- the log is complete and on disk
3. `fort.2` -- every one of the above happened

This also survives the stage combinations the engine is intended
to grow into. One signal at the end certifies whatever set of
stages was requested; a signal per stage cannot, which is what
today's failure demonstrates.

Within that slot the order is also forced. `initVerboseness`
(ARCHITECTURE 12.2) must run before the identity block, because
the block is gated on the mask that call sets. Reversed, the
banner would test an uninitialized mask. The two are adjacent
and the dependency is easy to miss, which is why it is stated
here rather than left for the pseudocode to infer.

**Not the tabular outputs.** The DOS, bond order, and optical
spectra files are consumed by plotting scripts that parse them
positionally. A banner in those files would be a breaking change
to a data format for no gain, since nobody reads them by eye.

**Not the HDF5 files either.** Attaching a citation string to a
root group was considered, on the argument that provenance
should travel with the data rather than with the log. It is not
adopted. A run writes three separate HDF5 files -- `hdf5SCF`,
`hdf5PSCF`, and `hdf5Field` -- so the question is never simply
"the HDF5 file", and answering it would mean deciding which of
them count as a primary result. More to the point, the same
reasoning that keeps the banner out of the tabular outputs
applies with equal force here: these are data files, read by
programs, and citation guidance in the human-readable output is
what the guidance is for. Cluttering the data to restate it buys
nothing a reader will ever see.

### 10.7 Suppression

The identity block is governed by the `banner` category of
ARCHITECTURE 12, and is included in the `normal` default. A
flight that does not want it sets `IMAGO_VERBOSENESS` without
`banner`.

Suppression is expressed through the environment rather than a
command-line flag because the engine has no flags. Its arguments
are positional -- `parseCommandLine` reads bare `getarg` values
in fixed order and `imago.py` builds the invocation as a
positional string -- so introducing an option would mean
changing the argument contract in Fortran and in Python
together, and in every other caller, for what is a cosmetic
toggle. The environment variable avoids that entirely, and a
flight sets it once rather than per unit.

The cost being avoided is real. The seed run of section 7
dispatched 87 units; at roughly thirty lines each the identity
block alone would contribute some 2,600 lines of decoration to a
single flight, and flights will grow.

The methods block is not suppressed. It is a handful of lines,
it is the part a reader is actually meant to copy into a paper,
and a flight that has turned the artwork off has no reason to
discard the citations the run earned.

### 10.8 Open Design Questions

- **Whether the build should generate `banner.txt`'s citation
  lines from `CITATION.cff`.** Doing so would remove the
  duplication accepted in 10.4 and make the `.cff` file the
  single source, at the cost of a generation step and a
  dependency on parsing YAML at configure time. Cheaper than the
  compiled-in variant it replaces, but still not free. Worth
  revisiting once the DOI exists and the fields stop changing.
- **Whether `imago.py` should print the banner too.** The driver
  is what a user invokes and its output is more visible than the
  engine's log. Reading the same `banner.txt` would cost little,
  but printing the block twice in one run would be worse than
  printing it once, so the two would have to agree on which of
  them owns it.
- **How a method with no predicate is caught.** A reference
  whose predicate can never be true is dead weight, and a method
  added without a predicate is silently uncited -- the failure
  that matters. Neither is visible without someone checking, and
  it is not obvious what would check it.
- **Whether `atomSCF`, `gaussFit`, and `contract` print it.**
  They are separate executables with their own log units and
  their own `open(20,...)` calls. They are also rarely run
  directly by a user composing a paper, so the benefit is
  smaller and the duplication real.


## 11. Partial Optical Properties: Decomposition Scheme

This section specifies which decompositions the partial
optical properties offer, how they are numbered, and why the
set stops where it does. It is the first part of the optical
properties code to be brought under the document chain.

### 11.1 What is being decomposed

A total optical calculation forms, for each transition
between an occupied state i and an unoccupied state j, the
momentum matrix element between them, and accumulates its
squared magnitude into a broadened spectrum.

A partial calculation splits that matrix element by where
its two ends sit. Each basis function is assigned to a
group; the matrix element is then resolved into
contributions M(a,b) from initial-state group a and
final-state group b, and the transition probability is
distributed over the resulting pair matrix so that the sum
over (a,b) reproduces the total.

The decisive structural fact is that **the quantity carries
two group indices, not one.** The partial DOS resolves a
single-index quantity, so its storage grows linearly in the
number of groups. The partial optical properties resolve a
pair, so storage grows as the square, and the square is
multiplied by the number of transition pairs, the number of
Cartesian components, the number of k-points and the number
of spins. Everything below follows from that.

### 11.2 The two axes, and the cells worth offering

A decomposition is fixed by two independent choices.

**Grouping** -- what a partial belongs to:

- **type**: every atom of one atomic type contributes to a
  single partial.
- **atom**: every atomic site gets its own partial.

**Resolution** -- how finely a group's basis functions are
split:

- **total**: all basis functions of the group in one
  partial.
- **nl**: one partial per radial function, so per QN_nl
  pair, summed over the m components of each shell.
- **nlm**: one partial per basis function, resolving the
  individual Cartesian components of each shell.

That is a two by three grid. **The nlm column is not
offered**, for two reasons that are independent of each
other, either of which would be sufficient.

The first is cost. At nlm resolution the partial count is
the full valence dimension, so the pair matrix is
`valeDim` by `valeDim` for every transition pair at every
k-point. For a small five atom cell with an extended basis
that is already terabytes. It is not a decomposition that
can be requested and waited for; it is one that cannot be
run.

The second is correctness, and it would bite even if the
storage were free. Individual m components mix under a point
group operation via the representation matrices D^l(R), so
an nlm-resolved quantity cannot be unfolded from an
irreducible wedge by relabeling atoms. Section 2.3 gives the
proof. Note that this applies to the type-grouped nlm cell
as well: summing over every atom of a type does not rescue
it, because the mixing is within a shell rather than between
atoms.

The four remaining cells are all offered:

```
                total            nl
  type      (type, total)   (type, nl)
  atom      (atom, total)   (atom, nl)
```

### 11.3 Numbering, and the principle behind the order

```
  code  grouping   resolution   partial count
  ------------------------------------------------
  0     --         --           no decomposition
  1     type       total        numAtomTypes
  2     type       nl           sum over types of
                                  their radial fns
  3     atom       total        numAtomSites
  4     atom       nl           sum over sites of
                                  their radial fns
```

Code 0 is the default and means an ordinary total optical
calculation with no decomposition performed.

**Grouping is the major key and resolution the minor one.**
The type-grouped cells come first and the atom-grouped cells
after, so the code number is monotone in how finely the
decomposition resolves position. This is not a cosmetic
choice. It makes the IBZ correction a threshold rather than
a memorized set: codes 1 and 2 need nothing, codes 3 and 4
need the atom permutation, and the boundary is a single
comparison. A numbering that interleaved the two groupings
would leave the reader with no way to tell which codes need
the correction except by looking them up, and a reader who
guessed would guess wrong half the time.

**Any cell added later must respect the ordering.** A new
type-grouped cell belongs among the low numbers and a new
atom-grouped one at the end. A change that breaks the
monotonicity costs more than the renumbering it saves,
because it silently invalidates every place that tests the
threshold.

**The numbering is deliberately independent of the partial
DOS detail codes.** The two quantities offer different sets
-- the partial DOS offers nlm and does not offer
(type, total), for the cost reasons above that apply to a
pair matrix and not to a single-index quantity -- so a
shared numbering would have to leave gaps in both. They are
separate schemes and each is stated in full where it is
defined.

### 11.4 Cost, and who is responsible for it

Storage for the pair matrix is

    partials^2 * 3 * transitionPairs * kPoints * spins

in double precision. The `partials` column of the table
above is what to substitute.

The four cells span a wide range. The atom grouped cells
grow with the cell: (atom, total) as the square of the site
count, and (atom, nl) as the square of the site count times
the radial functions per site, which makes it the most
expensive of the four by a wide margin. A few tens of atoms
at nl resolution reaches tens of gigabytes.

The type-grouped cells are bounded by the number of atomic
types instead, and **how much cheaper that is depends
entirely on what the types mean in the system at hand.** For
an ordinary crystal the type count is small and fixed, so a
type-grouped decomposition costs essentially nothing however
large the cell grows. A point defect supercell behaves the
same way and for a reason worth knowing: its types are
assigned from the *pre-defect* symmetry precisely to keep
the count down, so growing the supercell adds sites without
adding types. An amorphous cell does not behave this way at
all. There a type is a bin of locally similar environments
rather than a symmetry orbit, so the type count is set by
how finely the environments were binned and can grow with
the cell. A type-grouped decomposition of an amorphous cell
is therefore not automatically the cheap option, and in the
limit of one type per environment it approaches the cost of
the atom-grouped cell it was chosen instead of. See section
2.3 for what types do and do not mean.

**Storage is not the only cost, and on the cells that
matter it is not the binding one.** The Fortran calculation
produces the pair matrix; `processPOPTC.py` then turns each
pair into its derived spectra, and it does that by invoking
`makePDOS.py` and `imagoKKc` once per pair. Each call to
`makePDOS.py` re-reads the whole raw spectrum file, and that
file is itself quadratic in the partial count. So the
post-processing is quartic in a quantity the calculation is
only quadratic in, and it overtakes the calculation quickly.

Measured on a five atom cell with an extended basis: the
(type, nl) cell has 24 partials and spent 56 minutes in
post-processing, while the (atom, nl) cell has 36 and needs
roughly three hours. The Fortran stage of the latter took
two minutes in 300 MB.

That asymmetry is worth stating plainly, because it inverts
the obvious reading of the storage formula above. A request
can fit in memory comfortably and still be impractical, and
the user who chose it will watch a job that appears to hang
long after the calculation itself has finished. Reducing it
is a question about the post-processor rather than about the
decomposition, so it is not settled here.

**Imago does not currently estimate either cost at input
parse time.** A request that cannot fit is discovered when
the allocation fails rather than when it is made, and a
request that cannot finish in reasonable time is discovered
by waiting for it. Choosing an affordable decomposition is
presently the user's responsibility. A general facility for
projecting the resource cost of a requested calculation
before running it belongs with the resource and cost
guidance of section 8 rather than being built once here for
one quantity.

### 11.5 Relation to IBZ correctness

Which cells need an unfolding correction on a
symmetry-reduced k-point mesh follows from the grid
coordinates alone, and is stated with its proof in section
2.5. In summary: the type-grouped cells need nothing,
because every operation carries each atom onto an atom of
the same type and a type-level sum therefore maps onto
itself; the atom-grouped cells need the atom permutation
applied to both indices of the pair matrix.

Because the nlm column is not offered, **the partial optical
properties do not depend on the deferred D^l(R)
representation matrices at all.** Every offered cell is
correctable with the atom permutation alone. This is not
true of the partial DOS, whose nlm mode still waits on them.

That the type-grouped cells need no correction is a
statement about arithmetic, and it holds whatever the types
happen to mean physically, because `buildAtomPerm` verifies
the closure at startup rather than inferring it (section
2.3). The exception is style code 0, where no symmetry maps
are built and so nothing is verified.

### 11.6 What a type-grouped partial means

Choosing a grouping is not only a cost decision, because
what a type-grouped partial *is* varies with the system.

In a crystal it is what it appears to be: a sum over
symmetry-equivalent atoms, which are physically
indistinguishable, so nothing is lost by summing them.

In an amorphous cell a type is a bin of locally similar
environments carrying no symmetry content, so a type-grouped
partial is an average over a population of genuinely
different sites. That may be exactly what is wanted -- it is
the spectroscopic analogue of asking what a *kind* of
environment contributes -- but it is an average over
inequivalent things, not a redundant sum over equivalent
ones, and the spread within the bin is not recoverable from
the output.

In a point defect supercell the trap is sharper. Types are
assigned from the pre-defect symmetry, so the atoms
neighbouring the defect usually carry the same type as
chemically identical atoms far from it. A type-grouped
decomposition then averages the perturbed neighbourhood into
the unperturbed bulk, which dilutes the defect signature by
roughly the ratio of the supercell size to the neighbourhood
size -- and the larger the supercell, the more thoroughly it
is diluted. **A defect study almost always wants an
atom-grouped cell**, and wants it for a reason that has
nothing to do with symmetry: the decomposition must be able
to separate atoms that the type assignment deliberately does
not.

---

## 12. Optical Properties: Brillouin-Zone Integration

Section 11 specifies how the optical properties are
DECOMPOSED. This section specifies how they are INTEGRATED
over the Brillouin zone, which until now had no design-level
home at all: the optical properties arrived from OLCAO with
Gaussian broadening wired in and no alternative.

The section covers the total spectra and the partial ones
together, because they share one integration and differ only
in what is carried through it. Everything here applies to
both unless it says otherwise.

### 12.1 What is being integrated, and why it is not a DOS

For each pair of an occupied state i and an unoccupied state
j, the calculation forms the momentum matrix element between
them and accumulates its squared magnitude at the transition
energy. The quantity wanted is

    eps2(E) ~ (1/E^2) sum_{i,j} Int_BZ dk |M_ij(k)|^2
              delta( e_j(k) - e_i(k) - E )

This is a JOINT density of states weighted by the matrix
element. The distinction from section 1 is the whole of the
design problem and is worth stating before anything else:
the LAT machinery of section 1 finds the surface where ONE
band takes the value E, while this needs the surface where
the DIFFERENCE of two bands takes the value E.

That difference is what makes the extension tractable. Write

    epsDiff_ij(k) = e_j(k) - e_i(k)

and the integral above is exactly the section 1.4 problem
with `epsDiff_ij` in place of the band energy and
`|M_ij|^2` in place of the Mulliken projection. No new
analytic geometry is required. The existing corner-weight
routines apply unchanged, because they are functions of four
corner VALUES and do not care what quantity produced them.

### 12.2 Two pathways, and how one is chosen

Gaussian broadening is retained as a full alternative rather
than replaced. This follows the DOS precedent exactly, and
the reason is not only backward compatibility: the two
schemes fail differently, so a disagreement between them is
diagnostic, and losing one would lose that.

**The switch already exists.** `kPointIntgCode` in O_KPoints
is 0 for the histogram/Gaussian method and 1 for LAT. It is
already exposed to the user as `-scfkpint`, `-pscfkpint` and
`-kpint`, already read from the k-point file, and already
consumed by `computeDOS`, `computeBond`, `valeCharge` and
the SCF occupation path. What is new here is only that the
optical properties begin to honour it; today they ignore it
entirely and are Gaussian unconditionally. **No new input
option is introduced**, which keeps this clear of the
option-contract and cache-key consequences that a new switch
would carry.

**Which dispatch shape.** The DOS path uses both shapes and
they are not interchangeable. `computeTDOS_LAT` is a
separate routine that the caller selects with an explicit
`if (kPointIntgCode == 1)`; `computeDOS` instead branches
internally. The difference is whether the two methods share
the surrounding loop structure. They do not here: section
1.3 records that the LAT loop INVERTS relative to Gaussian
-- outer over bands and tetrahedra rather than over k-points
-- and the same inversion applies to the transition pairs.
**So the optical path follows the `computeTDOS_LAT` shape:
a separate accumulation routine, selected by the caller.**
An internal branch would be a branch around the whole body,
which is a separate routine wearing a disguise.

**Where the selection sits, and the renaming it forces.**
`computeTDOS_LAT` is selected in `subroutine dos` in
`imago.F90`, the top level of the program. The optical
counterpart is not currently reachable from the equivalent
place: `subroutine optc` calls `printOptcResults`, and that
routine builds the conversion factors, the energy grid and
`kPointFactor`, allocates the spectrum arrays, calls the
accumulators, calls the printers, and deallocates. Printing
is the last third of it, so the name already describes a
minority of what it does, and hanging a pathway choice on it
would make that worse.

The optical path is therefore restructured to match the DOS
path rather than merely imitating it:

- `optcCond` and `optcCondPOPTC` are promoted to module
  scope in O_OptcSpectra. This is the house pattern already --
  `transitionProb` lives at module scope in
  O_OptcTransitions -- and they are local today by habit
  rather than by design.
- `printOptcResults` splits into `computeOptcSpectra`, which
  performs the setup and selects the pathway, and
  `printOptcSpectra`, which writes the files.
- `subroutine optc` in `imago.F90` calls the two in turn,
  which is exactly where and how `subroutine dos` makes the
  same choice.

The accumulators are renamed with it, since the same
principle applies: `getOptcCond` does not get anything, it
accumulates a broadened spectrum into an array it is handed.
They become `accumulateOptcCond` and
`accumulateOptcCondPOPTC`, with `_LAT` suffixed counterparts
following the convention of section 1.5 -- the suffix marks
the integration method and the unsuffixed name is the
Gaussian one, as with `electronPopulation` and
`electronPopulation_LAT`. That the unsuffixed name means
Gaussian is a convention worth knowing, since nothing in the
name says so.

**What the two pathways share, and where they diverge.** An
earlier draft of this section claimed that only the
accumulation differs and that `computePairs` and
`computePOPTCPairs` were untouched. Section 12.4 shows that
is not so: the array those routines produce cannot be
consumed by a tetrahedron loop at all. The division is:

- **Shared:** the momentum matrix elements themselves, which
  are the expensive part, and the decomposition index of
  section 11. Both pathways need the same physics from the
  same eigenvectors.
- **Divergent:** how the resulting transition strengths are
  indexed, filtered, ranged and occupied. Five things differ
  and only one of them is the energy sort. The others are
  that the storage slot is a running counter over accepted
  pairs rather than a band pair; that pairs failing the
  transition-energy cutoff are dropped, per k-point, when a
  tetrahedron may need a pair that fails at one corner and
  passes at three; that the band range is per k-point where
  the tetrahedron path needs the union; and that the
  occupation factors come from `electronPopulation` rather
  than `electronPopulation_LAT`.

**So the shared physics is extracted rather than
duplicated.** The construction of `conjWaveMomSum` -- the
sum over basis functions of the conjugated wave function
against the momentum matrix -- becomes its own routine,
called by both pathways' producers. Each producer then
carries only its own bookkeeping. Duplicating the momentum
construction instead would put the one calculation whose
errors are hardest to see in two places; threading five
conditionals through a single producer instead would leave a
routine that does two different jobs and can honestly be
named for neither.

### 12.3 The Gaussian pathway, as it stands

Recorded here because it has never been written down, and
because the LAT pathway is specified as a departure from it.

`getOptcCond` walks spin, then IBZ k-points, then the
transitions at that k-point, then energy points. Each
transition deposits `transitionProb` into every energy bin,
scaled by a Gaussian in the distance from the transition
energy to the bin, and by `kPointFactor(i)`, which carries
`kPointWeight(i)` and the normalization. `getOptcCondPOPTC`
is the same loop over the pair matrix instead of the total.

Two properties of this arrangement matter later. The k-point
weight is the ONLY place Brillouin-zone geometry enters, so
the star multiplicity is applied as a scalar; and the
broadening width `sigmaOPTC` is a user input that controls
both the numerical convergence and the appearance of the
result at once. Section 12.5 returns to that conflation.

### 12.4 The LAT pathway

For each transition pair (i, j) and each tetrahedron T,
gather the four corner values of `epsDiff_ij`, sort them,
and obtain the four corner densities from the existing
`bloechlCornerDOSWt` in O_MathSubs. The contribution is

    dEps2(E) = (V_T / V_BZ) sum_{c=1..4}
               cornerDOSWt_LAT(c) * |M_ij(k_c)|^2

**The normalization, derived rather than asserted**, since
the two pathways must land on one scale and the factors do
not correspond one-to-one. The Gaussian path multiplies
each transition by

    kPointFactor(i) = kPointWeight(i) * 0.5
                      / (sigma * sqrt(pi)) / hartree / spin

and the LAT path replaces that with

    tetraVol * sum(kPointWeight) * 0.5 / hartree / spin

multiplying the corner density weight. Term by term:

- `kPointWeight(i)`, summed over the irreducible points,
  becomes `tetraVol` summed over tetrahedra. The first sums
  to 2 by Imago's convention and the second to 1, so
  `sum(kPointWeight)` restores the scale. This is the factor
  section 1.3 requires of every LAT accumulation.
- `1 / (sigma * sqrt(pi))` is the normalized Gaussian
  standing in for a delta function. It is DROPPED, because
  the corner density weight is that delta function evaluated
  exactly and already carries units of inverse energy.
- `0.5`, `1 / hartree` and `1 / spin` are properties of the
  quantity rather than of the integration, so all three
  survive unchanged. The `0.5` in particular is not a
  geometric factor: the transition probabilities already
  account for two electrons per state in the
  spin-unpolarized case, and it prevents the k-point weights
  from multiplying by two a second time. Dropping it on the
  assumption that it belonged to the k-point sum would halve
  every LAT spectrum.

Note that `integratePDOS_LAT` uses
`cornerDOSWt * tetraVol * kpWtSum / hartree` with no `0.5`
and no `1 / spin`, and is NOT a template to copy here. Its
spin division is folded into the projections upstream
(section 1.4), and the partial DOS has no `0.5` to carry
because its Gaussian counterpart has none either.

Since `kPointFactor` is built by the caller and passed in,
the caller builds the LAT factor instead when
`kPointIntgCode` selects that pathway; the Gaussian factor
is meaningless there, containing a `sigma` the pathway does
not use.

**The stored transition probabilities cannot be reused, and
this is the structural precondition for everything else.**
`transitionProb` is indexed `(component, pair, kpoint,
spin)`, where `pair` is a position in a list that
`computePairs` sorts by transition ENERGY before storing.
The band identity (i, j) is discarded by that sort. Worse,
the list is built from `firstOccupiedState` and its
companions, which are dimensioned `(numKPoints, spin)`, so
both the number of pairs and the range of bands entering
them vary from one k-point to the next.

A tetrahedron needs the SAME band pair at all four of its
corners. A sorted pair index cannot supply it: position p is
a different (i, j) at each corner, and at a metal's Fermi
surface a pair present at one corner may not be enumerated
at another at all. So the LAT pathway requires the matrix
elements stored under a band-pair index, `(component, i, j,
kIBZ, spin)`, rather than the energy-sorted index the
Gaussian pathway uses. The transition energy is then
`e_j - e_i` recomputed from the eigenvalues rather than
looked up, which is what the corner values need in any case.

This is the analogue of the two-pass requirement that
section 1.4 imposes on the PDOS, and it is sharper: the PDOS
needed its projections merely to survive to a second pass,
while this needs them re-indexed.

**The index ORDER is a performance decision and must be
stated, not left to look arbitrary.** Fortran stores the
leftmost index fastest, so a slice is contiguous only when
the colons sit on the left. The existing arrays already
respect this -- `transitionProb(dim3, pair, kpoint, spin)`
is read as `(:, j, i, h)` -- and the obvious imitation for
the banded store, `(dim3, i, j, kIBZ)`, would break it. The
reason is the access pattern rather than the shape: with the
band pair fixed, the tetrahedron loop fetches four different
k-points per tetrahedron, so putting `kIBZ` last strides by
`dim3 x nOcc x nUnocc` on every corner and misses cache on
each one.

The correct order is `(dim3, kIBZ, i, j)`. With the band
pair fixed, the whole block `(:, :, i, j)` is contiguous:
three components by the k-point count, a few tens of
kilobytes for a typical mesh. It loads once and stays
resident while every tetrahedron for that band pair is
processed, which turns the innermost fetch from a strided
miss into a cache hit.

The partial array needs the same care and does not get the
same result. It is far too large to hold one band pair's
slice in cache, so the goal there is narrower: the pair
matrix over partials, for one component at one corner, must
be contiguous, and the loop over partials must run the
leftmost index innermost. The destination remains scattered
whatever is done, because the permutation of section 12.6 is
the entire point of the operation.

A later reader will be tempted to reorder these into
something that looks tidier. The orders are chosen against
the loops that consume them; changing one means changing the
other.

Note that the occupation factors are not a new problem. The
Gaussian path already multiplies each transition by the
initial state's occupancy and the final state's vacancy, so
the f_i(1 - f_j) structure exists; a band-pair index simply
lets it be evaluated per corner. That is most of what
item (c) below asks for.

Beyond that, four things must be got right, and each is a
silent wrong answer rather than a failure.

**(a) The sort permutation applies to the matrix element
too.** `bloechlCornerDOSWt` returns its four weights in
SORTED corner order. The matrix elements must be carried
through the same permutation before they are paired with
those weights. Pairing a sorted weight with an unsorted
matrix element is the single easiest mistake here, and it
produces a plausible spectrum rather than an obviously
broken one. Section 1.4 states the same requirement for the
Mulliken projections and is the model to follow.

**(b) The corners are FULL-mesh points; the matrix elements
are stored at IBZ points.** This is section 1.4's
"fundamental constraint" again. The resolution is the same:
store per IBZ k-point and map at corner assembly through
`fullKPToIBZKPMap`, applying the operation in
`fullKPToIBZOpMap`. For the partial properties the
decomposition index must be permuted at the same moment --
which is precisely what `partialPerm` already does (section
2.5, PSEUDOCODE 7a). See 12.6.

**(c) Occupied and unoccupied are properties of a corner,
not of a tetrahedron.** The Gaussian path classifies bands
per k-point, through `firstOccupiedState` and its
companions, and never has to reconcile two k-points. A
tetrahedron spans four, and near a Fermi surface a band may
be occupied at some corners and empty at others. Once the
matrix elements are stored per band pair the enumeration
problem goes away -- every (i, j) in range is present at
every corner -- and what remains is that the occupation
factors must be evaluated corner by corner rather than
carried in from the pair's own k-point.

**For a gapped system the distinction is inert**, since
every corner classifies the same way, and that is the case
to implement and validate first. For a metal two things
change together: the factors vary across the tetrahedron,
and the occupations themselves should come from
`electronPopulation_LAT` rather than from the Gaussian
path's `electronPopulation`, so that one scheme determines
both the geometry and the filling. That pairing is the same
argument section 1.6(a) makes for the Fermi level, and for
the same reason -- mixing the two leaves an error the
calculation partly absorbs. It is a real extension and
should be specified separately rather than assumed to fall
out of this one.

**`electronPopulation_LAT` is not available on the optical
path, and this must be arranged rather than assumed.**
`computeElectronPopulation_LAT` is called from `subroutine
bond` in `imago.F90`, and `subroutine optc` is a sibling
that never calls it. A run that asks only for optical
properties never enters `bond` at all, so the array is
simply absent. The optical path must therefore call
`computeElectronPopulation_LAT` itself when
`kPointIntgCode` selects the tetrahedron pathway, at the
same point in its sequence that `bond` does -- after the
eigenvalues have been shifted to put the Fermi level at
zero, so that the energy argument is zero for the same
reason it is there. Reading an unallocated array is the
loud failure here; reading a stale one from an earlier
`bond` call in a combined run would be the quiet one.

**(d) Degenerate corners.** When `epsDiff_ij` is nearly
constant over a tetrahedron the sorted values coincide and
the analytic denominators vanish. Section 1.3 already
requires guards for this in the band case, and they carry
over -- but the frequency is different and the difference is
worth flagging. Parallel bands make `epsDiff` flat by
construction, and parallel bands are exactly what produce
the sharp critical-point structure that optical spectra are
computed to show. The degenerate branch is therefore not a
rare edge case here but the physically interesting one, and
its guards deserve direct testing rather than inheritance.

### 12.5 What the broadening parameter means under each path

A reader who knows section 1.6(e) will expect the LAT path
to ignore the broadening input, as the SCF occupation path
ignores `thermalSigma`. **That expectation is wrong here,
and the difference is worth being explicit about.**

Under Gaussian broadening `sigmaOPTC` does two jobs at once.
It is the numerical device that turns a sum over discrete
k-points into a continuous spectrum, and it is also the
physical broadening that a measured spectrum genuinely has,
from finite lifetimes and instrumental resolution. The two
are indistinguishable in the output, which is why a
converged Gaussian result is a surface over (sigma, N_k)
rather than a number.

LAT removes the first job and leaves the second untouched.
The tetrahedron integration produces the spectrum without
any smearing parameter, but almost nobody wants to LOOK at
that spectrum: an unbroadened joint density of states is a
spiky object that no measurement resembles. So the LAT path
does not ignore `sigmaOPTC`; it changes what the parameter
means, from a numerical requirement to a physical model
applied afterwards, and one the user may legitimately set to
zero to see the raw result.

Two consequences follow. A LAT run and a Gaussian run at the
same `sigmaOPTC` are NOT the same calculation broadened two
ways, and should not be compared as though the parameter had
one meaning. And the convergence question changes shape: the
LAT answer converges in mesh density alone, which is the
whole point of adopting it.

### 12.6 Where this leaves the IBZ correction

The partial optical properties currently correct for IBZ
reduction by averaging the pair matrix over the star of each
IBZ k-point (section 2.5, PSEUDOCODE 7a). Under LAT that
correction does not move -- it DISAPPEARS, and is replaced
by the corner assembly of 12.4(b).

The reason is structural. The star average exists because
the Gaussian path visits only IBZ points and must
redistribute each one's contribution over the members of its
star. LAT visits full-mesh corners directly, so there is
nothing to redistribute; the permutation is applied once per
corner as the matrix element is fetched, which is the same
arithmetic arriving at the same answer by a shorter route.

**This is the sequencing argument, and it is the reason this
section is being written before TODO O3.** O3 must rotate
the Cartesian components of the momentum operator, and
PSEUDOCODE 7a already records that doing so requires lifting
the star average up to the complex matrix element. Under LAT
that rotation belongs at corner assembly instead. Both
changes therefore land in the same place, and settling the
integration question first means writing that code once. It
does NOT mean LAT solves O3: the components mix under a
point operation however the zone is integrated. LAT decides
only where the rotation is applied.

### 12.7 Cost

The Gaussian accumulation is O(numKPoints x pairs x
energyPoints) over IBZ points. The LAT accumulation is
O(numTetrahedra x pairs x energyPoints), and
`numTetrahedra` is `24 x numFullMeshKP` (section 1.2 cuts
every box four ways, so it is 24 and not 6). So the inner
work grows by twenty-four times the IBZ reduction factor,
which is the same full-mesh scaling section 1.6(b) notes for
the SCF occupation path, and it is a substantial increase
rather than a rounding error.

That factor is large enough to be worth attacking rather
than accepting, and there is an obvious way in. The
accumulation currently touches each corner once per
containing tetrahedron and multiplies by the whole partial
matrix each time, so a k-point's matrix is walked 96 times.
Accumulating the corner WEIGHTS per (k-point, energy) first
and multiplying by the partial matrix once afterwards does
the same arithmetic with the multiplications divided by 96,
at the cost of holding a weight array. That restructure is
recorded as its own task rather than folded in here, because
it changes the loop structure of a routine this section
specifies and is not needed for correctness.

Set against it: LAT is expected to reach a converged answer
at a lower mesh density, which is the reason to accept the
per-mesh cost. Whether the trade is favourable for the
optical properties specifically is a measurement, not a
prediction, and the comparison should be made at equal
converged accuracy rather than at equal mesh.

Note also that the transition pair count already dominates
both expressions and is quadratic in the band count, so
neither path is cheap and the decomposition cost of section
11.4 sits on top of whichever is chosen.

### 12.8 Open questions

  1. **Metals**, per 12.4(c). The corner-occupation
     weighting is sketched, not specified. A gapped system
     should be implemented and validated first, with metals
     as a separate increment.
  2. **Whether the Kramers-Kronig path needs anything.**
     `imagoKKc` consumes eps2 and produces eps1 and the
     derived spectra. It should be indifferent to how eps2
     was integrated, but this has not been checked, and
     section 12.5 changes the spectrum's character near
     sharp features, which is where a quadrature is most
     easily embarrassed. See TODO O6 for what that routine's
     integration is already known to do.
  3. **The per-axis columns.** These are unverified on a
     reduced mesh today (TODO O3) and this section does not
     repair them. It only relocates where the repair goes.
  4. **Validation target.** The natural first check is a
     gapped cubic system where the two pathways must agree
     in the converged limit, run against the same mesh
     ladder. A disagreement that persists with mesh density
     is the informative outcome, and per 12.5 the comparison
     must hold the meaning of `sigmaOPTC` fixed rather than
     its value.

---

## 13. Optical Cartesian Components Under Symmetry

### 13.1 The defect, and the one thing that survives it

The momentum operator is a vector. Under a point group
operation its Cartesian components mix:

  P_i(Rk) = sum_j R_ij P_j(k)

The IBZ unfolding of section 2.4 handles the SITE index --
`atomPerm` says which atom becomes which -- and nothing
handles the COMPONENT index. Every member of a k-point's
star is therefore credited with the representative's
orientation, unrotated.

Exactly one combination survives, and saying which is worth
the space, because it is why this went unnoticed. Rotating a
vector does not change the sum of the squares of its
components:

  sum_c |P'^c|^2 = sum_c |P^c|^2

The isotropic column is that sum divided by three, so it is
correct on any mesh for any crystal. The three individual
columns are a redistribution error.

**This axis is independent of the atom decomposition, and
the output makes that easy to miss.** `printSpectrumPOPTC`
writes the same header for the undecomposed unit and for
every atom-pair unit:

```
COL_LABELS 4
TOTAL x y z
```

"TOTAL" there means direction-AVERAGED, not undecomposed. A
run with `detailCodePOPTC = 3` on the five-atom cell writes
26 units, and every one of them -- the undecomposed unit and
all 25 atom pairs -- carries its own isotropic column and
its own x, y, z columns. So this section is about those
three columns wherever they appear, and says nothing about
the atom decomposition, which is sections 11 and 1.7's
concern.

**Every isotropic column in the file is correct, including
the partials', and that needs its own argument.** The
undecomposed case is the sum of squares above. A partial
stores a cross-term rather than a modulus, but the same
orthogonality carries it:

  sum_c M'_c(o,n) conjg(S'_c)
      = sum_{c,d,e} R_cd R_ce M_d(o,n) conjg(S_e)
      = sum_d M_d(o,n) conjg(S_d)

using sum_c R_cd R_ce = delta_de. So the defect costs the
x, y and z columns of every unit and none of the isotropic
ones.

Measured on cubic KNbO3, 4x4x4 shifted mesh, peak of eps2:

```
                          x        y        z    isotropic
unreduced, 64 k-points  48.470   48.470   48.470   48.470
reduced,    4 k-points  37.607   72.701   85.237   48.470
```

The unreduced run needs no unfolding and comes out isotropic
on its own to every printed digit. That is both the proof
that the momentum matrix elements themselves are sound and
the reference the reduced run fails to reproduce.

### 13.2 What can be recovered, and what cannot

One trustworthy number is available per energy point. How
far it goes is decided by how many independent directional
values the crystal's symmetry allows.

```
  cubic                    1 value    recoverable
  tetragonal, hexagonal,   2 values   not recoverable
    trigonal
  orthorhombic             3 values   not recoverable
  monoclinic, triclinic    3 values, plus off-diagonal
                                      terms not computed at
                                      all today
```

A cubic crystal is fully recoverable, and not by any clever
averaging: symmetry forces the three directions equal, so
there is one unknown and the reliable sum supplies it.
Averaging the three columns over the point group amounts
there to printing the isotropic value three times.

**Orthorhombic is the trap, and the reason group averaging
must never be offered as a partial remedy here.** Its
operations carry each axis onto plus or minus ITSELF and
never onto another axis, so averaging over the point group
changes nothing at all, leaves the three columns exactly as
wrong as they were, and gives no sign that it did nothing.
An apparent remedy that is a silent no-op is worse than
none.

**The workaround that works today** is to not reduce the
mesh: with every k-point diagonalized there is no unfolding
and nothing to rotate. Worth knowing before judging urgency.
Symmetry reduction saves most in cubic crystals, where
directional spectra carry no information beyond the
isotropic one, and least in low-symmetry crystals, where
they are physically interesting -- a monoclinic cell reduces
by about four. The workaround is cheapest exactly where it
is needed.

### 13.3 Why this is a storage change

The obvious repair -- build the Cartesian rotations and
apply them where each star member is credited -- cannot
work, and the reason decides the whole design.

`computePairs` forms the squared modulus the instant the
matrix element exists (`optc.F90:1399`), and the tetrahedron
producer does the same:

```fortran
valeValeXMom(k) = sum(valeVale(:,i,1) &
      & * conjWaveMomSum(:,finalStateIndex,k))
transitionProbTemp(k,transPairCount) = ( &
      & real(valeValeXMom(k),double)**2 &
      & + aimag(valeValeXMom(k))**2) &
      & * initStateFactor*finStateFactor
```

A squared modulus cannot be rotated. Under R,

  |M'^c|^2 = sum_{d,e} R_cd R_ce M^d conjg(M^e)

which needs the OFF-DIAGONAL products M^d conjg(M^e). The
code keeps only the three diagonal entries of a rank-two
tensor, so by the time anything could rotate, what it needs
has been discarded.

**Decision: keep the complex matrix element.** Store M^c as
three complex numbers rather than three squared moduli,
rotate at the deposit, and square there. Correct for every
crystal system and both integration pathways.

Storing the full Hermitian tensor M^d conjg(M^e) instead
also works and is strictly worse: nine complex numbers where
three suffice, since the tensor is reconstructible from the
vector at the point of use.

**Occupancy factors stay separate.** They are real scalars
per (initial band, final band, k-point) and are folded into
the stored product today. Folding them into the stored M as
a square root would be wrong: the final-state vacancy
`1 - occupancy` can go slightly negative through round-off
near a Fermi surface, and the square root of a negative
number is a crash rather than a rounding error. Store M
unscaled and apply the real factor at the deposit, where the
band and k-point indices are already in hand.

### 13.4 The Cartesian rotation matrices

The rotations exist only in FRACTIONAL form:
`convAbcPointOps` as read from the k-point file, and
`abcRealPointOps(3,3,numPointOps)` conjugated into the
loaded direct basis by `computeRealPointOps` (section 2.7).
There is no Cartesian form anywhere in `src/`.

  R_cart = L R_abc L^-1

with L the real lattice vectors as COLUMNS, which is
`O_Lattice`'s `realVectors`. Note this is the transpose of
the row layout the k-point file's `CONV_LATTICE` block uses;
section 1.2's module data comment records the same trap.
This is the conjugation `computeRealPointOps` already
performs, into one more basis, so the new array belongs
beside it in `O_KPoints` as `xyzRealPointOps(3,3,numPointOps)`,
built in `initializeKPoints` on the same branches.

Section 2.6 used to name `xyzPointOps` and `xyzFracTrans`
among the arrays the SYBD branch skips allocating. Neither
identifier has ever existed in the source; that prose is
corrected alongside this section, and its point survives --
`xyzRealPointOps` is built on the same branches as its
fractional sibling and is therefore absent on the SYBD path
and for k-point style code 0.

### 13.5 Where the rotation is applied

The two pathways reach IBZ correctness by different routes
(section 12.6), so the rotation lands in a different place
on each.

**Tetrahedron pathway.** It already visits full-mesh corners
directly and applies that corner's operation as it fetches,
so the rotation joins the fetch. The corner supplies both
the IBZ k-point and the operation index; the stored complex
vector is rotated by `xyzRealPointOps(:,:,opR)` and squared
there.

**Gaussian pathway, totals.** The harder one, and the
difficulty is structural rather than incidental: it never
visits star members at all. It accumulates over IBZ points
with `kPointFactor` carrying the star multiplicity through
`kPointWeight`. There is no loop for a per-member rotation
to live in.

Two ways to give it one:

  1. Add a star loop, mirroring the tetrahedron pathway.
     Correct, and it multiplies the accumulation by the
     reduction factor, 4 to 48, for a quantity whose
     isotropic part was already right.
  2. Precompute the star-summed rotation. What is wanted is
     the sum over the star of |M'^c|^2, which is
     sum_{d,e} (sum_R R_cd R_ce) M^d conjg(M^e). The inner
     sum depends only on which operations make up that
     star -- not on the band, the energy, or the matrix
     element -- so it is a fixed object per IBZ k-point,
     computed once. The accumulation applies it to the
     tensor formed from the stored vector, and the star
     loop disappears.

**Option 2 is chosen.** It leaves the Gaussian pathway's
cost where it is, and it writes down explicitly what that
pathway has always done implicitly -- sum over the star --
instead of hiding it inside a multiplicity factor.

**The object carries FOUR component indices, not three.**
The expression above is written for a diagonal entry, where
both factors of the squared modulus carry the same component
`c`. Direction code 2 (section 13.7) asks for the off
diagonal entries of the symmetric tensor as well, and there
the two factors carry DIFFERENT components:

```
  sum over star of  M'^c1 conjg(M'^c2)
     = sum_{d,e} ( sum_R R_c1,d R_c2,e ) M^d conjg(M^e)
```

So the precomputed object is `R_c1,d R_c2,e` summed over the
star -- 81 numbers per IBZ k-point rather than 27. Setting
`c2 = c1` recovers the diagonal form exactly, so codes 0 and
1 read their entries out of the same array and nothing about
them changes.

This was missed when the section was first written, and the
consequence is worth naming because it is not visible in any
result: a three index object gives correct spectra for
direction codes 0 and 1, which is every case that had been
run. It fails only for code 2, and it fails by having no
entry to read rather than by producing a wrong number, so
the omission surfaces as a missing capability rather than as
a defect in an answer.

### 13.6 The decomposed case

The partial store holds, for each ordered pair of partials,
the real part of a product between that pair's matrix
element and the total over all pairs:

  transProbPOPTCBanded(o,n,c) =
      Re[ M_c(o,n) conjg( S_c ) ],
      S_c = sum_{o',n'} M_c(o',n')

Summing over (o,n) returns the undecomposed transition
probability, which is the sum rule section 11 depends on.

Both factors carry the component index, so both must be
rotated, and the rotated quantity is Re[ M'_c(o,n)
conjg(S'_c) ]. The store therefore keeps the complex
M_c(o,n) and, separately, the complex total S_c. S_c is one
vector per (initial band, final band, k-point, spin) rather
than one per pair, so it costs nothing beside the pair
matrix.

**The star average of section 7a must move inside the
component loop.** As written, `computePOPTCPairs` walks the
star with the component index held FIXED outside it
(`optc.F90:2596`), which is correct only while an operation
cannot mix components. Once it can, the component index has
to sit inside the star walk and the scratch slab has to
carry it. This is the one place in the change where an
existing loop nest is reordered rather than extended, and it
is the first place to look if the sum rule breaks.

### 13.7 Direction resolution, its control, and its cost

**The cost ladder is blunter than it looks.** Per
transition, per partial pair:

```
  what is wanted            stored        against today
  ----------------------------------------------------
  isotropic only            1 real        one third
  diagonal x, y, z          3 complex     twice
  full symmetric tensor     3 complex     twice
  full Hermitian tensor     3 complex     twice
```

Anything directional costs the SAME to store, because it all
comes from the same three complex numbers at the deposit.
So the storage decision is binary -- isotropic or
directional -- and choosing between the diagonal and the
full tensor costs accumulation time and output width only.

The other end is the useful surprise: an isotropic-only mode
is a third of today's storage rather than a compromise, and
the isotropic column is the one most users read.

**Three levels, controlled from `OPTC_INPUT_DATA`.** This is
an optical output choice rather than an integration one, so
it sits beside `detailCodePOPTC` rather than in the k-point
file where sections 1.2 and 1.7 put their settings. It forms
a THIRD axis on section 11's grouping-by-resolution grid
rather than crossing into it: crossing would turn five
decomposition codes into fifteen.

```
  0   isotropic only.  One column per unit.  Correct on any
        mesh today, and the cheapest thing the code can do.
  1   diagonal.  The TOTAL x y z columns written today,
        made correct on a reduced mesh.
  2   full symmetric tensor.  Six components, xx yy zz xy
        xz yz, which Imago cannot express at all today and
        which monoclinic and triclinic crystals require.
```

Level 2 is nearly free once level 1 exists, since it needs
no extra storage -- only a wider accumulation and a wider
output record.

**Totals and partials are set independently.** The
undecomposed store carries no `sumNumPartials` factor and is
small; `transProbPOPTCBanded` is sized as

  sumNumPartials^2 x 3 x numKPoints x nOcc x nUnocc x spin

so it grows as the SQUARE of the atom count and is already
the binding array -- which is why TODO O8 exists. Doubling
it matters most exactly where it is already the limit.
Tying the two settings together would make anyone who wants
directional totals on a large system pay for directional
partials they may not want, so they are separate fields.

The natural combination on a large low-symmetry system is
therefore level 2 on the totals and level 0 on the partials:
the full dielectric tensor where it is affordable, and the
direction-averaged decomposition where it is not.

**What level 2 opens.** A dielectric response is a
three-by-three tensor whose off-diagonal parts are
physically real for monoclinic and triclinic crystals and
are not computed at all today. Once the complex vector is
stored they are products of quantities already in hand.
Writing them changes the output record width and therefore
`imagoKKc` and `processPOPTC.py`, which is why level 2 is a
level rather than the default.

**What is NOT added.** No diagonalization on either
pathway, and the Gaussian pathway's star-summed rotation is
built once per IBZ k-point rather than once per transition.
Section 11.4's cost projection and TODO O8 need the storage
factor above.

**Level 0 changes what SOME of the derived spectra mean, and
the split is not where one would first guess.** What matters
is not whether a quantity is linear, but whether `imagoKKc`
forms its isotropic value by AVERAGING the three
per-direction values or by computing it from the isotropic
epsilon directly. Level 0 has no per-direction values to
average, so only the first group can move.

```
  from the isotropic values directly -- unaffected
    epsilon-1, epsilon-1i, the energy loss function
  averaged over the three directions -- level 0 differs
    refractive index, extinction coefficient,
    reflectivity, absorption
```

Measured on cubic KNbO3, level 0 against level 1:

```
  epsilon-1              1e-5   agrees
  energy loss function   1e-6   agrees
  refractive index       5e-2   genuinely different
```

Epsilon-1 agrees because Kramers-Kronig is linear, so
transforming the average equals averaging the transforms;
the residual is interpolation on the fine grid. The energy
loss function agrees because it is built from the total
epsilon-1 and epsilon-2 rather than averaged -- which is
also the correct construction, since the average of
Im(-1/eps) is not Im(-1/eps_avg).

The refractive index differs because level 1 reports the
average of n(eps_x), n(eps_y) and n(eps_z) while level 0
reports n(eps_avg), and n is a square root. Neither is
wrong; they answer different questions. But a user moving
between the two levels will see it shift by several percent,
and an unexplained five percent step is indistinguishable
from a bug, which is why it is written down here.

**Level 2 is a strict superset of level 1.** It emits
everything level 1 emits, plus epsilon-2 and epsilon-1 for
the three off-diagonal entries. Nothing that level 1
produces is withdrawn or changed.

```
  epsilon-2, epsilon-1     xx yy zz xy xz yz
  ELF, n, k, R, alpha      total, x, y, z
```

The asymmetry in that table is the whole of the design
question, and it deserves stating rather than hiding.
Kramers-Kronig acts on each tensor component independently,
because causality holds element by element, so epsilon-1
generalizes to all six without argument. The other five do
not. They are functions of a SCALAR complex dielectric
function: the energy loss function is properly
`Im(-eps^-1)`, a matrix inverse rather than an element-wise
one, and for an anisotropic medium the refractive index
follows from the Fresnel equation, which admits two
eigenmodes per propagation direction and no "index along
xy" at all.

**So why emit the per-axis five at level 2?** Because
level 1 already emits them, for every crystal, including the
monoclinic and triclinic ones where they rest on x, y and z
being principal axes -- which is exactly the assumption a
non-zero off-diagonal entry denies. Level 2 does not make
that worse. It makes it VISIBLE: the off-diagonal columns
sitting beside the diagonal ones are the measurement of how
far the assumption is from holding, and a reader who
compares `eps_xy` against `eps_xx` learns something no
level 1 run can tell them.

The alternative considered was to withhold the five at
level 2 and offer the isotropic column alone. It was
rejected because it would make level 2 give LESS than
level 1 for the same run, which reads as a regression rather
than as a caution, and because a user who wants the
diagonal values would simply run level 1 and get them with
no warning attached at all.

**Diagonalization is deferred, not refused on the merits.**
The physically complete answer is to diagonalize
`eps(omega)` at each frequency and report along its
principal axes. That is real work with real subtlety: the
tensor is complex symmetric, so the transformation is
complex-orthogonal rather than unitary; the principal axes
rotate with frequency in monoclinic and triclinic cells; and
eigenvalue branches have to be tracked to stay continuous.
It is worth doing when a study needs principal-axis optics,
and it is a separate task rather than a condition on level 2.

### 13.8 Seam inventory

Every quantity the new code consumes or produces: where it
comes from, who builds it, and when.

```
quantity            supplied by           lifetime
------------------------------------------------------------
abcRealPointOps     computeRealPointOps,  built in
                      from convAbcPointOps  initializeKPoints
                      (section 2.7)         for style codes
                                            0, 1, 2; NOT on
                                            the SYBD branch
realVectors         O_Lattice, read from  resident for the
                      the structure file    whole run; its
                                            COLUMNS are the
                                            lattice vectors
xyzRealPointOps     NEW. Built beside     same branches and
                      abcRealPointOps as    lifetime as
                      L R_abc L^-1          abcRealPointOps
fullKPToIBZOpMap    the IBZ fold in       before optc runs;
                      initializeKPointMesh  absent for style
                                            code 0
transitionProb,     computePairs and      become COMPLEX;
  transProbBanded     computeTransProb-     allocation sites
                      Banded                and lifetimes
                                            otherwise unchanged
transProbPOPTC-     computeTransProb-     becomes COMPLEX and
  Banded              POPTCBanded           gains a companion
                                            holding the
                                            per-pair total S_c
starRotationSum     NEW. One 3x3x3 per    built once after
                      IBZ k-point from      the IBZ fold,
                      that point's star     released with the
                                            k-point set
cartesianCodeOPTC,  NEW. readOptcControl, read once per run,
  cartesianCode-      from fort.5's         before any store
  POPTC               OPTC_INPUT_DATA;      is sized -- both
                      O_Input               stores are
                                            dimensioned from
                                            them
```

Two consequences bind the code. The permutation tables are
absent for k-point style code 0, and `xyzRealPointOps` will
be too since it is built on the same branches, so every
consumer guards on allocation and writes to fort.20 when it
skipped -- exactly as section 1.7 requires of the
symmetrization. And the SYBD branch builds no point-op
machinery at all (section 2.6), so a band structure run must
not reach any of this.

### 13.9 Relation to the partial DOS restriction

Section 1.4 refuses `detailCodePDOS == 3` on the tetrahedron
pathway because per-lm projections "require D^l(R) rotation
matrices". The three Cartesian components of a vector ARE
the l = 1 representation, so the optical columns are the
l = 1 instance of that same gap: the unfolding handles the
site index and nothing handles the representation index.

For l = 1 the matrix is R itself, which makes this the cheap
end of one piece of work rather than a separate problem. A
later section closing the PDOS restriction should generalize
`xyzRealPointOps` to the D^l(R) it needs rather than
introduce a second mechanism beside it.

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
parallelepipeds. Each parallelepiped has 8 corners and is
decomposed into exactly 6 tetrahedra that tile without overlap.
The standard decomposition (Bloechl 1994) shares the main
diagonal M1-M8:

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

The mesh is periodic (the BZ is a torus), so indices wrap with
modular arithmetic: `mod(a, nA) + 1` etc. The total count is:

  numTetrahedra = 6 * nA * nB * nC

All tetrahedra span an equal fraction of the BZ:

  tetraVol = 1 / numTetrahedra

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
expressions are in Bloechl 1994, equations 14-16. The Bloechl
correction terms (eqs. 22-24) substantially improve accuracy at
lower k-point densities and should be included.

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
For mode 0 (per atom-type, per l-shell), R preserves
species, so no channel permutation is needed.

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
                          invariant under R
  1     per-atom total    invAtomPerm(R, atomIdx)
  2     per-atom, per-l   invAtomPerm remaps atom;
                          l-shell offset unchanged
                          (same species ⇒ same
                          orbital structure)
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
  ────────────────────────────────────────────────────

*Mode 0 sums projections over all atoms of the same
type. Since point-group operations permute atoms only
within the same species, the type-level sum is
automatically invariant -- no correction needed.

**Mode 3 resolves individual Cartesian Gaussian
components (px, py, pz separately). Under rotation these
mix via D^l(R), so atom permutation alone does not
suffice. Correct unfolding requires the full
representation matrices for each l-shell. This is
deferred -- mode 3 is rarely used in practice.

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
  order, PDOS) will not be correct unless the user has
  taken extreme care to provide a symmetric k-point mesh.
  For style codes 1 and 2, Imago builds the full mesh
  and all symmetry maps internally, so the atom permutation
  fix works for both Gaussian and LAT integration paths.
  makeinput.py no longer produces style code 0 files;
  mesh mode (`-kp`) now writes style code 1.

**SYBD path bypasses atomPerm.**  Symmetric band structure
(`-sybd`, `-pscfsybd`) replaces the loaded k-point set with
a 1-D path between user-specified high-symmetry vertices.
On that path every k-point is its own end product -- band-
structure output is per-k-point eigenvalues, and the planned
partial decomposition (future work) is a direct per-atom
projection at the very k-point being plotted, with no star
to unfold.  There are simply no shell-summed quantities for
atomPerm to reconstruct, so the table is not needed.

This matters in practice because the SYBD branch of
`initializeKPoints` calls `makePathKPoints` and skips all of
the point-ops setup that the style-code 0/1/2 branches do
(`numPointOps` assignment, `xyzPointOps`/`xyzFracTrans`
allocation, `computeRealPointOps`).  `abcRealPointOps` and
`abcRealFracTrans` therefore stay unallocated.
Consequently, the calls to `buildAtomPerm` and
`buildInvAtomPerm` in `setupSCF` (SCF path) and `intgPSCF`
(PSCF path) are guarded with `if (doSYBD_SCF /= 1)` and
`if (doSYBD_PSCF /= 1)` respectively.

The downstream consumers of `atomPerm` and `invAtomPerm`
are `computeBond` (effective charge / bond order star
distribution) and the LAT PDOS channel-permutation step
in `dos.F90`.  Both are themselves gated by their own
`doBond_*` / `doDOS_*` flags, so a pure `-sybd` (or
`-pscfsybd`) run with no decomposition flag never reaches
them.  Combining `-sybd` with `-bond` or `-dos` is not
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
`convAbcFracTrans`, reads the `CONV_LATTICE` block as
`M_conv`, takes `M_loaded = realVectors`, and forms a
single change-of-basis matrix

```
C = M_loaded^{-1} * M_conv = invRealVectors^T * M_conv
```

once, before the per-operation loop.
`computeRealPointOps` then transforms each operation
as

```
R_real_abc = C * R_conv_abc * C^{-1}
t_real_abc = C * t_conv_abc
```

producing `abcRealPointOps` and `abcRealFracTrans` in
the basis of whichever cell `O_Lattice` holds -- which
is exactly the basis `buildAtomPerm` uses for atom
positions.  In `full` mode the `CELL_MODE` flag tells
the routine that `M_loaded == M_conv`, so `C = I` and
the loop collapses to a copy from the on-disk arrays
into the runtime arrays without invoking the
conjugation kernel.  In `prim` mode the flag selects
the full conjugation path.  Both identities hold for
any non-singular `M_loaded` and `M_conv`; no cell-
shape, centering, or symmetry assumption is built
into the math.

`computeRecipPointOps` does the matching reciprocal-
space transform with the same change-of-basis matrix
`C`, so kpoint folding via `abcRecipPointOps` and atom
permutation via `abcRealPointOps` continue to descend
from a single source of operations.  No full/prim
branching exists anywhere outside the two transform
routines, where the `CELL_MODE` flag selects identity
versus full conjugation.

**Generality.**  This design works for every cell type
imago describes.  Old uolcao had a latent version of
the same basis-mismatch bug for non-cubic systems:
`makeKPoints` used a Cartesian-assumed conjugation but
consumed conv-abc input that was only equivalent to
Cartesian for cubic-aligned conventional axes.  The
pre-`buildAtomPerm` pipeline never propagated the bug
because the SCF binary received pre-folded explicit
k-points and never saw a rotation matrix, so the
wrongness silently cancelled.  Imago carries
operations all the way to the SCF binary and uses them
for both reciprocal-space folding and real-space atom
permutation, so the consumer-side conjugation is
necessary for correctness across cell types.

**Diagnostic history.** The issue surfaced when
diamond/prim (SG 227 Fd-3m reduced to its primitive
rhombohedral cell) stopped at `buildAtomPerm: no atom
match found` after the new IBZ-correctness machinery
landed.  The Fd-3m mirror perpendicular to cubic x,
applied directly to the primitive-cell atom positions,
sent atom 2 from `(0.25, 0.25, 0.25)` to
`(0.25, -0.25, -0.25)`, wrapped to
`(0.25, 0.75, 0.75)` -- a vacancy.  The non-symmorphic
translation `(0.25, 0.25, 0.25)` happens to round-trip
through the transform with the same numerical value
(the body-diagonal-of-conv-cube coincidence shared
with the primitive-cell body diagonal), but rotations
like the cubic-x mirror require the actual conjugation
to produce the right primitive-basis matrix.  An
earlier iteration of this section moved the
conjugation onto the producer side via a Cartesian-xyz
intermediate; the current design relocates it to the
consumer with `C = M_loaded^{-1} * M_conv` and keeps
the on-disk values tied to the spaceDB entries.  Under
either form the same operation maps atom 2 to itself
(a self-image fixed point of the mirror), and
diamond/prim runs cleanly through `buildAtomPerm` and
the rest of the SCF.

---

## 3. Density-Based K-Point Input Pipeline

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
   the space group database (`_extract_point_ops`) and
   writes kpoint files (`kp-scf.dat` and/or `kp-pscf.dat`)
   directly, bypassing the `makeKPoints` executable.
2. Each file uses `kPointStyleCode=2` and contains:
   - `KPOINT_STYLE_CODE` = 2
   - `KPOINT_INTG_CODE` = 0 (default; histogram)
   - `MIN_KP_LINE_DENSITY` = the user's volume density
     (label is historical; the value is a volume density)
   - `KP_SHIFT_A_B_C` = the shift (auto or user-specified)
   - `NUM_POINT_OPS` = number of point group operations
   - `POINT_OPS` = the 3x3 rotation matrices (abc coords)
3. At runtime, `imago` reads this file in `readKPoints`
   (including the point group operations), then in
   `initializeKPoints` calls:
   - `computeAxialKPoints` — density to axial counts
   - `computeRecipPointOps` — convert abc point ops to
     reciprocal-space abc operations
   - `initializeKPointMesh(1)` — build the uniform mesh
     and fold it to the IBZ using the reciprocal-space
     point group operations
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
declares.

**Top-level keys (required):**

  Field            Type    Description
  --------------------------------------------------------
  schema_version   int     Currently 2.  The reader
                           rejects any other value.
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
  coefficients     array   Length = num_gaussians.
                           Gaussian expansion
                           coefficients (col 1 of the
                           legacy coeff1 file).
  alphas           array   Length = num_gaussians.
                           Authoritative explicit alpha
                           values per basis function
                           (col 2 of coeff1).  May
                           diverge from the geometric
                           series implied by alpha_min
                           and alpha_max if a future
                           entry uses a non-geometric
                           layout.

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

Additional fields on the record are matcher-specific and
not validated by the schema.  As examples for the two
matchers Phase 2 ships with:

- `"bispectrum"` records carry `values` (array of reals,
  length `2 * twoj2 + 1`).
- `"reduce"` records carry a `shell_code` inline table
  encoding per-level distance, neighbor count, and
  neighbor element/species multiset.

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
  convergence_threshold  float   SCF convergence
                                 threshold of the
                                 reference run.
  scf_iterations         int     Iteration count of the
                                 run that produced this
                                 potential.  Feeds the
                                 validation harness
                                 (5.8).

**Validation rules** (enforced at load time):

1. `schema_version` must equal 2.
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
   hard error.  No implicit fallback to `"isolated"`:
   the curator must declare the selection explicitly.
8. Within any one `[[potential]]` entry's fingerprint
   array, the pair `(method, sub_spec)` must be unique.
   Two fingerprint records with the same `method` and
   the same `sub_spec` keys-and-values are a hard
   error.  Records with the same `method` but
   differing `sub_spec` are explicitly allowed (e.g.,
   one bispectrum fingerprint at `twoj1=8, twoj2=8`
   alongside another at `twoj1=6, twoj2=4`).
9. `method` must be one of the matcher names
   registered in `makeinput.py` at load time (8.9).
   Unknown methods are a hard error rather than a
   silent skip, so that a typo in the manifest fails
   loudly rather than quietly omitting the fingerprint
   from the lookup.

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
source                = "Imago"
commit                = "fedcba2"
generated_at          = "2026-05-08T14:30:00Z"
reference_id          = "COD-1011098"
atom_site             = 1
kpoint_spec           = "12 12 12 0 0 0"
convergence_threshold = 1.0e-6
scf_iterations        = 28

[[potential.fingerprint]]
method   = "bispectrum"
sub_spec = { twoj1 = 8, twoj2 = 8 }
values   = [
   1.2345678901234567e-01,
   ...  (9 entries; length = 2 * twoj2 + 1)  ...
]

[[potential.fingerprint]]
method   = "bispectrum"
sub_spec = { twoj1 = 6, twoj2 = 4 }
values   = [
   ...  (5 entries)  ...
]
```

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
handling or SCF-runner code they don't use,
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
diff of a regenerated file.  It does **not** imply
file-level byte-identity across pipeline runs: the
regeneration pipeline (5.7) refreshes provenance
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
-pot LABEL              Manual override.  Apply LABEL
                        uniformly across the structure
                        (where no other scheme picks
                        an entry).  Optional.  See
                        precedence in 5.6.3.

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
4. **Coverage check** (only when an environment-based
   scheme is in play).  Confirm that at least one entry
   in the loaded database carries a fingerprint record
   matching the requested `(method, sub_spec)` for the
   active scheme.  If not, abort with a message naming
   the element and the missing `(method, sub_spec)`.
   This fails fast before the expensive nested-makeinput
   bootstrap (5.10) starts.

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
- **Environment-based** (`-reduce scope=N`, or
  `-bispec scope=N`).  Determine the atom set:
    - No `scope=`: all atoms of the active element.
    - `scope=N`: only atoms inside the spatial region
      named `N`.
    - `scope=~N`: only atoms outside the named region.
  Compute the per-atom fingerprint vectors over that
  set (in-Python for `-reduce`; via the nested loen
  bootstrap of 5.10 for `-bispec`).  Bucket atoms whose
  vectors agree within tolerance into species (the
  matcher's `distance` and `default_similarity_floor`,
  ARCHITECTURE 8.9).  Atoms outside the scope keep
  whatever assignment earlier flags produced.

After the species pass, every atom has a final
`atom_species_id[atom]`.

#### 5.6.5 Manifest-entry pick per species

For each `(element, species)` pair appearing in the
final assignment, pick exactly one `PotentialEntry`
from that element's database.  Precedence, top to
bottom:

1. **`-pot LABEL` (manual override).**  If given, every
   species in every element uses the entry with
   `label == LABEL`.  `KeyError` for any element that
   lacks the label is fatal — `-pot` is a deliberate
   override and silent fallback would mask user intent.
2. **Fingerprint match** (only for atoms assigned by an
   environment-based scheme).  For each species' atoms,
   ask the matcher to summarize them into one
   representative fingerprint via its `representative`
   method (8.9).  The matcher chooses semantics
   appropriate to its descriptor space: `BispecMatcher`
   returns the element-wise mean of the member vectors;
   `ReduceMatcher` returns the first member's shell-code
   (intra-species reduce fingerprints agree within the
   matcher's tolerance by construction, so any member
   is equally good).  Future matchers may use a medoid
   or another scheme; the protocol does not pin the
   choice.  Among the element's database entries that
   carry a fingerprint record matching
   `(method, sub_spec)`, pick the entry whose recorded
   fingerprint minimizes
   `distance(representative, entry_fingerprint)`.  If
   the best distance exceeds the matcher's
   `default_similarity_floor`, fall through to step 3
   with a warning naming the species and the best-but-
   rejected entry.  Each matcher carries a heuristic
   default for its similarity floor; the concrete
   numbers (e.g., 0.05 for `ReduceMatcher`, 0.10 for
   `BispecMatcher`) are starting values intended to be
   tuned during the Phase-2 validation pass (TODO C61)
   against the benchmark systems' actual fingerprint-
   distance distributions, and users may override them
   per scheme on the CLI when a particular system
   warrants a tighter or looser tolerance.
3. **Default tag.**  For any species not yet assigned a
   potential — atoms grouped by a position-based flag,
   atoms outside any environment-scheme scope, or
   atoms whose fingerprint match was rejected at the
   similarity floor — use `default_entry(db)` for that
   element.  Guaranteed to succeed by rule 7.

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

### 5.7 Regeneration Pipeline Algorithm

**Purpose of `build_initial_potentials.py`.**  This
is the script that **produces** the augmented
per-element database files
(`share/atomicPDB/<elem>/s_gaussian_pot.toml`).  It
takes a curated set of reference solids as input,
runs Imago SCF on each (or loads cached results),
harvests converged potentials at named atom sites,
and writes the results into per-element database
files via `initial_potential_db.save()`.  It is the
**producer** half of the library / producer /
consumer split documented in 5.4.

**The curation manifest -- what it is, what it
contains, why.**  The pipeline's primary input is
a single human-readable file: the **curation
manifest**.  It declares *which solids, which atom
sites, under which labels, with which SCF
settings*.  Its job is threefold:

1. **Declare the curation set.**  The manifest
   *is* the curation strategy made explicit.
   Adding a reference solid means adding a
   manifest entry; removing one means removing
   an entry; reviewing the curation set means
   reading the manifest.
2. **Tell the pipeline what to harvest.**  For
   each reference solid, which atom sites'
   converged potentials go into the database,
   and under what labels.
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
curation choice captured in one version-controlled
file alongside the structure files it points at,
so the regeneration becomes a deterministic
function of (manifest, structure files, Imago
build).

**What it contains** (per the schema sketched
below): per-solid fields the SCF run needs
(`structure_path`, `kpoint_spec`,
`convergence_threshold`) plus a stable
`reference_id`; per-entry harvest declarations
(`atom_site`, expected `element`, `label`,
`description`).

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
  regeneratable build artifact, not an archival
  hand-curated table.
- **Provenance metadata (free).**  Timestamps,
  commits, and similar fields refresh on every run
  and are exempt from any reproducibility
  guarantee.

**Inputs:**

- A curation manifest (TOML; schema v1 specified below).  Default
  location is `share/atomicBDB/manifest.toml`; alternate via the
  `--manifest` flag.
- The existing `share/atomicPDB/` tree (for `pot1` and `coeff1`
  reads when refreshing each element's `"isolated"` baseline
  entry).
- An Imago build location, for running reference SCF calculations.
- Network access to the Crystallography Open Database (COD).  Only
  consulted on cache miss for `[[reference_solid]]` entries that
  declare a `cod_id`; cache hits never touch the network.

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

[[reference_solid]]
reference_id          = "au_fcc"
# Exactly one of (cod_id, cod_revision) or structure_path is set
# per [[reference_solid]] (validation rule 4 below).
cod_id                = 9008463
cod_revision          = "2023-04-12"
# structure_path      = "au_fcc.skel"     # alternative form
kpoint_spec           = { density = 60.0, shift = [0.0, 0.0, 0.0] }
convergence_threshold = 1.0e-6

  [[reference_solid.entry]]
  element     = "Au"
  atom_site   = 1
  label       = "default_solid"
  default     = true
  description = "Au in fcc bulk (Fm-3m)."

    [[reference_solid.entry.fingerprint]]
    method   = "bispectrum"
    sub_spec = { twoj1 = 8, twoj2 = 8 }

    [[reference_solid.entry.fingerprint]]
    method   = "bispectrum"
    sub_spec = { twoj1 = 6, twoj2 = 4 }

    [[reference_solid.entry.fingerprint]]
    method   = "reduce"
    sub_spec = { level = 2, thick = 0.5, cutoff = 5.0,
                 tolerance = 0.05 }

[[reference_solid]]
# ... another solid ...
```

**Per-solid fields.**

- `reference_id` (string): stable, human-readable identifier for
  the reference solid.  Used as the cache directory name; must
  match `[A-Za-z0-9_-]+` and be unique across the manifest.
- `cod_id` (positive integer, optional iff `structure_path` is
  set): Crystallography Open Database entry ID.  The pipeline
  fetches the structure at regeneration time on cache miss.
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
- `convergence_threshold` (real): SCF convergence threshold for
  the reference run.  Recorded in provenance.

**Per-entry fields (`[[reference_solid.entry]]`).**

- `element` (string): element symbol the entry contributes to.
  Cross-checked against the species at `atom_site` after Imago
  loads the structure.
- `atom_site` (positive integer): 1-based site index into the
  structure, matching Imago's site-indexing convention everywhere
  else in the codebase.
- `label` (string): the label this entry is written under in the
  element's `s_gaussian_pot.toml`.  `(element, label)` must be
  unique across the entire manifest (rule 6).
- `default` (bool): whether this entry should be tagged
  `default = true` in the element's `s_gaussian_pot.toml`.
  Exactly one entry per element across the entire manifest is
  marked `true` (rule 7), giving each per-element file its
  required-single default entry per 5.2 rule 7.
- `description` (string): one-sentence prose explanation of the
  chemical environment; copied verbatim into the database
  entry's `description` field (5.2).

**Per-entry fingerprint declarations
(`[[reference_solid.entry.fingerprint]]`).**  Each
declaration tells the producer to compute and harvest one
fingerprint record alongside the numerical potential.  The
record is written into the database entry's
`[[potential.fingerprint]]` array per 5.2.

- `method` (string): matcher name.  Must be a method known
  to the matcher registry in `makeinput.py` (8.9).
- `sub_spec` (inline table): method-specific parameters.
  Two declarations on the same entry with the same `method`
  and the same `sub_spec` keys-and-values are a hard error
  (rule 8); same `method` with differing `sub_spec` is
  explicitly allowed.

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
   `kpoint_spec`, `convergence_threshold`, and *exactly one* of
   `{(cod_id, cod_revision), structure_path}` (see rule 4 for
   details).
3. Every `[[reference_solid.entry]]` carries `element`,
   `atom_site`, `label`, `default`, `description`.
4. Exactly one of `cod_id` or `structure_path` is set on each
   `[[reference_solid]]`.  If `structure_path`, it resolves to
   an existing file under the manifest's directory.  If `cod_id`,
   it parses as a positive integer *and* `cod_revision` is
   present as a non-empty string.
5. `reference_id` is unique across the manifest.
6. `(element, label)` is unique across the manifest.  Two solids
   cannot both produce, e.g., `("Au", "default_solid")` — the
   database holds one entry per `(element, label)`, so silent
   overwrite would mask a curation mistake.
7. For every element appearing on any
   `[[reference_solid.entry]]`, exactly one such entry across
   the entire manifest carries `default = true`.  Zero or
   multiple defaults for the same element is a hard error.
   This rule mirrors per-element-database rule 7: the
   manifest is the single source of truth for the default
   tag.
8. Within any one `[[reference_solid.entry]]`'s fingerprint
   declarations, `(method, sub_spec)` is unique.  Same
   `(method, sub_spec)` declared twice on the same entry is
   a hard error.  Same `method` with differing `sub_spec`
   coexists on the same entry by design.
9. Every `method` in a `[[reference_solid.entry.fingerprint]]`
   declaration must be a name known to the matcher registry
   in `makeinput.py` (8.9).  Unknown methods are a hard
   error rather than a silent skip.

**Cache layout and contract.**

SCF runs are expensive (minutes to hours).  The pipeline caches
results per reference solid so that adding or editing
`[[reference_solid.entry]]` declarations does not re-trigger SCF
— only the cheap harvest step (`extract_potential` in PSEUDOCODE
11.4) re-runs.  The cache lives under `share/atomicBDB/cache/scf/`,
one directory per solid keyed by `reference_id`:

```
share/atomicBDB/cache/scf/
  au_fcc/
    inputs.toml        Snapshot of the per-solid SCF inputs:
                       kpoint_spec, convergence_threshold,
                       imago_commit (plus cod_id + cod_revision,
                       or structure_path, for diagnostic
                       reporting on a miss).
    structure.<ext>    Byte-for-byte copy of the structure file
                       used.  For COD entries this is the file
                       fetched from COD; for structure_path
                       entries it is a copy of the on-disk file.
                       Compared byte-for-byte on every cache
                       check.
    imago.out          Cached SCF output (and any peer files
                       extract_potential needs).
    loen-<m>-<s>.out   One file per declared fingerprint:
                       `<m>` is the matcher name (e.g.,
                       `bispec`); `<s>` is a deterministic
                       slug derived from the `sub_spec` keys
                       and values (e.g., `twoj1_8-twoj2_8`).
                       Only present for Fortran-side
                       matchers; Python-side matchers
                       (`reduce`) compute in-process and
                       leave no cache file.
```

`is_cached(ref, imago_commit)` (PSEUDOCODE 11.4) is defined as:

1. If `cache/scf/<reference_id>/inputs.toml` does not exist,
   miss.
2. Read `inputs.toml`; compare each scalar field against `ref`
   and `imago_commit`.  If any field differs, miss; log which
   field changed (e.g., "`convergence_threshold` 1e-6 → 1e-7 —
   re-running SCF").
3. Materialize the current structure file for `ref` (fetch from
   COD if `cod_id`, read from disk if `structure_path`).  Compare
   it byte-for-byte against
   `cache/scf/<reference_id>/structure.<ext>`.  If different,
   miss; log "structure file changed for `<reference_id>`."
4. Otherwise, hit; return cached `imago.out`.

**Direct comparison over hashing.**  At our scale (~100 reference
solids, regeneration run by hand) a content hash buys nothing —
storage is trivial, lookup is one disk read — and costs
debuggability.  A hash miss tells the curator "different hash";
a direct-comparison miss can name the exact field that changed.
Storing a literal copy of the structure file also makes the
cache directory fully self-describing: open it and see exactly
what produced the cached result.

Explicitly excluded from the cache comparison:

- The `entries` list.  Adding a new `[[reference_solid.entry]]`
  must reuse the cached SCF result; only the harvest step
  re-runs.
- `reference_id` itself.  Used only as the cache directory name;
  the comparison logic never reads it from `inputs.toml`.
  Renaming a reference solid in the manifest will *appear* as a
  miss because the cache directory is new — the curator can
  rename the cache directory by hand if they want to preserve
  the hit.
- `structure_path` as a string.  We compare file *contents*, not
  paths, so a curator can rename `au_fcc.skel` →
  `gold_fcc.skel` on disk without re-running SCF.

**COD-fetch contract.**

- On cache miss for a `cod_id` reference solid, the pipeline
  fetches the structure at the pinned `cod_revision` from the
  Crystallography Open Database and writes it to
  `cache/scf/<reference_id>/structure.<ext>` before running SCF.
- On cache hit the network is never touched.  The cached
  `structure.<ext>` is the source of truth for the cached SCF
  result, by construction.
- Fetch failures (network down, COD outage, revision missing)
  are strict: the pipeline errors out, names the failing fetch,
  and refuses to fall back to any other revision.  Silent
  fallback would produce an SCF result inconsistent with the
  pinned manifest — exactly the failure mode the cache contract
  is designed to prevent.

**Procedure:**

1. Load and validate the manifest (nine rules above).
2. For every element with a directory in `share/atomicPDB/`,
   refresh (or create) the `"isolated"` entry of that element's
   `s_gaussian_pot.toml` directly from the current `pot1` and
   `coeff1` files.  Guarantees the baseline is always present
   (rule 6 of 5.2) and tracks any changes in atomSCF output.
3. For each `[[reference_solid]]` in the manifest:
   a. If `is_cached(ref, imago_commit)` succeeds (and `--force`
      was not passed), load
      `cache/scf/<reference_id>/imago.out`.  Otherwise:
         - For `cod_id` entries, fetch the structure at
           `cod_revision` from COD; for `structure_path`
           entries, read it from disk.  Write the bytes to
           `cache/scf/<reference_id>/structure.<ext>`.
         - Run Imago SCF on that structure using `kpoint_spec`
           and `convergence_threshold` from the manifest.
         - Write `cache/scf/<reference_id>/inputs.toml` and
           `cache/scf/<reference_id>/imago.out`.
   b. Record SCF iteration count and convergence metrics in the
      run log.
   c. For each `[[reference_solid.entry]]`:
      i.   Extract converged potential coefficients and alphas
           for the named `atom_site`.
      ii.  For each `[[reference_solid.entry.fingerprint]]`
           declaration:
              - Compute the fingerprint at `atom_site` for the
                requested `(method, sub_spec)`.  Python-side
                matchers (e.g., `reduce`) compute in-process
                from the reference structure.  Fortran-side
                matchers (e.g., `bispectrum`) run `imago.py
                -loen -scf no` on the reference structure,
                cached as
                `cache/scf/<reference_id>/loen-<m>-<s>.out`,
                and parse the row for `atom_site` from
                `fort.21`.
              - Build a `FingerprintRecord` (5.4) and attach
                it to the entry-in-progress.
      iii. Construct a `PotentialEntry` (5.2) with the
           manifest-supplied `label`, `default`, and
           `description`, the run-supplied numerical fields,
           the run-supplied provenance, and the
           `FingerprintRecord` list assembled in step ii.
      iv.  Insert the entry into the in-memory `ElementDatabase`
           for its element.  If an entry with the same label
           already exists, replace it.
4. Save each affected `ElementDatabase` to disk via
   `initial_potential_db.save()` (5.5).
5. Write `share/curation/run_log.toml` capturing the manifest
   snapshot, per-run iteration counts, and the Imago commit.
   The validation harness (5.8) reads this log.

**Flags:**

- `--force`: re-run every manifest entry from scratch, bypassing
  the cache check.  Cache entries are still written afresh
  afterwards, so the next ordinary run benefits from the warm
  cache.
- `--manifest PATH`: alternate manifest location (default:
  `share/atomicBDB/manifest.toml`).
- `--element ELEM`: restrict regeneration to a single element's
  `s_gaussian_pot.toml`.  Reference solids whose entries
  contribute only to other elements are skipped at the harvest
  step, but their SCF cache is still warmed so a follow-up run
  without `--element` benefits.

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

### 5.10 Nested-makeinput Bootstrap for Fortran-side Matchers

Matchers split into two families by where the
descriptor computation lives:

- **Python-side** (`reduce`).  Computes from
  `StructureControl` in-process during the species
  pass of 5.6.4.  No external program invocation; no
  intermediate files.
- **Fortran-side** (`bispectrum`).  Computes inside
  `imago.f90`'s `O_LocalEnv::bispec` path, which needs
  a populated `imago.dat` to read structure, k-points,
  and the `(twoj1, twoj2)` parameters.  The Imago build
  produces `fort.21` machine-readable output: one row
  per potential site, `2 * twoj2 + 1` real
  bispectrum-component values plus a sum.

The producer (5.7) and the consumer (5.6) both need to
run Fortran-side matchers but for different
structures:

- The **producer** runs them on each curated reference
  structure during harvest, harvesting one row from
  `fort.21` per declared
  `[[reference_solid.entry.fingerprint]]`.  This happens
  at database-build time, in `build_initial_potentials.py`.
- The **consumer** runs them on the *user's* structure
  during the species pass, harvesting one row from
  `fort.21` per atom of the active element.  This
  happens at input-generation time, inside
  `makeinput.py`.

The producer already has a full Imago-run orchestration
in place (the SCF caching loop of 5.7).  Adding a
follow-on `imago.py -loen -scf no` invocation is a
small extension of that loop.

The consumer's orchestration is the new piece this
section specifies.  `makeinput.py` is itself the script
that produces `imago.dat` files; to run loen on the
user's structure during the species pass it must, in
effect, generate a throwaway `imago.dat` and then run
loen against it -- a nested self-invocation.

#### 5.10.1 When the bootstrap runs

A bootstrap is triggered when both conditions hold:

1. The active environment-based scheme is a
   Fortran-side matcher (`needs_loen_run = true`
   on the matcher's protocol object, 8.9).
2. The per-element preflight (5.6.3 step 4) has
   already confirmed that every element in the
   structure has at least one manifest entry carrying
   a fingerprint matching the active `(method,
   sub_spec)`.

Python-side matchers (`reduce`) skip the bootstrap
entirely; they compute fingerprints in-process and
proceed to 5.6.4's bucketing.

#### 5.10.2 Bootstrap procedure

The bootstrap runs from inside `makeinput.py`'s species
pass, before any environment-based grouping decision is
made:

1. **Set up the scratch directory.**  Use the existing
   `.inputTemp/` scratch convention; create a
   `.inputTemp/loen_bootstrap/` subdirectory so the
   bootstrap inputs are isolated from the main run's
   eventual outputs.
2. **Build a minimal `imago.dat` in the scratch
   directory.**  Re-invoke `makeinput.py` as a
   subprocess with a stripped-down argument list:
      - The user's structure file (same skeleton, same
        lattice, same atom list).
      - No `-pot LABEL` flag.  The nested call's
        per-element preflight (5.6.3) loads each
        element's database file and falls through to
        the default-tagged entry on the no-scheme
        branch -- exactly the path we want.  Each
        element gets its curated default potential
        without any per-element `-pot` machinery on
        the CLI.
      - No grouping flags (`-target`, `-block`,
        `-reduce`, `-bispec`, `-xanes` all suppressed).
        Every atom in the nested call is its own
        species and type.  The grouping doesn't matter
        for loen, which iterates over potential sites
        regardless.  With no environment-based flag
        active, the nested call's selection algorithm
        never reaches step 2 of 5.6.5; the
        default-tag fallback in step 3 is the sole
        manifest-pick branch exercised.
      - A trivial k-point spec (Γ-only is sufficient).
        `loen` reads the k-point block but does not
        use it for fingerprint output.
      - The active matcher's loen-side parameters
        (`BispecMatcher.to_loen_input(sub_spec)`, etc.)
        flow into the existing
        `LOEN_INPUT_DATA` block that `makeinput.py`
        already emits (currently with hardcoded
        defaults; see 5.10.5 for the parameter
        contract and C58 for the wiring).

The **recursion guard** -- both the no-environment-
flag construction and the explicit
`--no-loen-bootstrap` flag of 5.10.3 -- is the sole
mechanism preventing the nested call from triggering
its own bootstrap.  No special `-pot`-driven bypass
is needed.
3. **Run loen.**  Invoke
   `imago.py -loen -scf no` against the scratch
   `imago.dat`, writing `fort.21` into the scratch
   directory.  Errors abort the parent makeinput run
   with a message naming the scratch directory so the
   user can inspect what was attempted.
4. **Parse `fort.21`.**  Read the per-site row.  For
   each potential site, build a fingerprint vector of
   length `2 * twoj2 + 1` (the matcher's
   `parse_loen_output` method, 8.9).
5. **Hand back to 5.6.4.**  The species pass now has a
   per-atom fingerprint vector for the active matcher
   and proceeds to bucketing.

The nested-call boundary is the key architectural
property: the bootstrap is a self-contained subprocess
whose only output is `fort.21`.  The parent run is
unaffected by any side-effects of the nested call
beyond the scratch directory's contents.

#### 5.10.3 Recursion guard

The nested `makeinput.py` invocation must not itself
trigger a bootstrap.  This is enforced two ways:

- The nested call uses `-pot LABEL` rather than any
  environment-based scheme, so no matcher is selected
  and no bootstrap fires by the trigger condition of
  5.10.1.
- As belt-and-suspenders, the nested call carries an
  internal `--no-loen-bootstrap` flag (or equivalent
  environment variable) that suppresses the bootstrap
  unconditionally.  If a future code path somehow
  enables a matcher during the nested call, the flag
  prevents infinite recursion.

#### 5.10.4 Caching deferred

`fort.21` recomputes from scratch on every makeinput
run.  Per the parked memory note's
"caching-deferred" decision, the rerun cost is judged
acceptable for the first cut.  If user feedback
indicates otherwise, a content-keyed cache (structure
bytes + `(method, sub_spec)` + Imago commit) under
`.inputTemp/loen_cache/` can be added without
algorithmic change to the bootstrap itself.

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
                                       `2 * twoj2 + 1`.
  max_neigh      `sub_spec.get(        Integer cap on
                  "max_neigh", 20)`    the per-site
                                       neighbor list
                                       length.
  cutoff         `sub_spec.get(        Real radial
                  "cutoff", 5.0)`      cutoff in Bohr
                                       on the
                                       neighbor list.
  angleSqueeze   `sub_spec.get(        Real angular
                  "angle_squeeze",     compression
                  0.85)`               factor (see
                                       `loen.f90`
                                       notes on
                                       `angleSqueeze`).

The required `sub_spec` keys are `twoj1` and `twoj2`.
The remaining three are optional, with the defaults
shown -- matching the current hardcoded values
`makeinput.py` emits today.  Two fingerprints whose
`sub_spec` differs in *any* of these five values
produce different bispectrum vectors and must coexist
as separate fingerprint records per DESIGN 5.2 rule 8.
When comparing `sub_spec`s for rule-8 uniqueness, the
canonical form is the post-default-resolution dict
(omitting a key is equivalent to specifying the
matcher's documented default for that key).

When element-aware bispectrum lands (TODO D10 / C62),
`BispecMatcher.to_loen_input` gains a `bispecByElement`
flag from `sub_spec.get("by_element", False)`.

Centralizing the mapping in the matcher keeps producer
and consumer aligned by construction.

---

## 6. High-Throughput Calculation Campaigns

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
layer reaches through; the kaleidoscope runner (D13),
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
   campaign must learn this as *data* (a result it can
   record and skip), not as an exception that aborts
   the whole batch.
2. **Where is the converged potential?**  An absolute
   path to the converged `scfV` output file (the
   `<edge>_initPot-<basis>.dat` that today's
   `manage_output` writes from `fort.8`).  The producer
   reads the Gaussian coefficients out of that file
   directly; the alphas it already knows, because they
   are an *input* (the min/max/number it fed into
   makeinput).  "Converged `scfV` matches input `scfV`"
   (5.7) is exactly this: coefficients come from the
   output, alphas from the input, taken together.
3. **Under what conditions did it run?**  The SCF
   settings actually used -- basis, k-point spec,
   convergence threshold, Imago build commit -- so the
   producer can fill the provenance fields of 5.2 and
   so kaleidoscope can form its run-reuse cache key
   (`kpoint_spec` + `convergence_threshold` +
   `imago_commit` + structure bytes, ARCHITECTURE 9.6).
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
  outputs           dict[str, str]: logical name ->
                      absolute path for each output file
                      produced (e.g. "scfV", "energy",
                      "iteration", plus property-specific
                      keys like "tdos", "bond"); the
                      producer reads outputs["scfV"]
  job               echo of identity: edge, job_name,
                      basis_scf, basis_pscf
  runtime_seconds   float: wall-clock time of the run,
                      for kaleidoscope's status.toml
  message           human-readable summary or error text
```

`RunStatus` is an enum with four members, chosen so the
campaign layer can branch on outcome without parsing
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
  raised) so the campaign can record-and-continue.
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

The boundary on error handling is deliberate and
important for Principle 10.  *Run-level* failures
(non-convergence, a Fortran abort, a missing input
file) are reported as a returned `ImagoResult` with the
appropriate status -- they are normal outcomes of
running real calculations and must not abort a
campaign.  *Contract* failures (the environment is not
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
  This is the mode kaleidoscope's default runner uses,
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
leave a campaign worker stranded in a stale temp
directory and corrupt every subsequent run's relative
path resolution.  This is the single most important
correctness difference between the CLI's one-shot
assumption and the API's reentrant requirement.

#### 6.1.5 Lock-file and checkpoint behavior preserved

Both existing robustness mechanisms carry over
unchanged in *meaning*; the design only clarifies how
they behave under concurrent, in-process use.

**The lock file is already per-run-directory, so
campaign concurrency is safe by construction.**  Today
the lock (`imagoLock`) lives in the `temp` directory,
which `get_temp_dir` derives by mirroring the run
directory's path under `$IMAGO_TEMP`.  Two different run
directories therefore mirror to two different temp
directories and two different lock files.  A
kaleidoscope campaign running thousands of independent
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
is never a process-global or campaign-global lock.

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

### 6.2 kaleidoscope campaign runner

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
below it (the runner, 6.2.2) or above it (the client,
6.2.6).  The two other load-bearing principles are 8
(the runner seam keeps dispatch independent of the
execution adapter) and 10 (complete-and-report: one
failed unit never aborts the campaign).

#### 6.2.1 The unit of work and the campaign

Kaleidoscope's data model is two plain, domain-agnostic
records.

```
CalcUnit
  id            stable per-structure key (6.2.4); the
                  curation reference_id for the producer,
                  a COD id for an acquisition campaign
  calc          optional calc-variant tag (6.2.4); None
                  when the structure hosts one calc
  structure     path to an imago.skl (or a structure
                  handle the chosen runner understands)
  options       makeinput options for this unit
  runner        which runner executes it (6.2.2);
                  defaults to the campaign default
  key_fields    client-declared cache identity (6.2.5):
                  scalar fields + names of key files

Campaign
  root          workspace root directory (6.2.4)
  units         list[CalcUnit]
  default_runner   runner used when a unit names none
  parsl_config  the Parsl Config (deployment, 6.2.3)
  on_outcome    optional per-unit callback (6.2.6)
```

A client builds a `Campaign` in process -- kaleidoscope
is a library first (Principle 9), not a CLI -- and hands
it to the dispatch entry point.  Kaleidoscope serializes
the campaign to `campaign.toml` in the workspace root so
a campaign is inspectable and a resume has an
authoritative record of *what was asked for*, separate
from `status.toml`'s record of *what happened* (6.2.4).
Whether `campaign.toml` may also be hand-authored as the
primary surface, rather than always generated from the
in-process `Campaign`, is left open (6.2.7).

The producer (C48) is the worked example throughout: it
builds one `CalcUnit` per curated reference solid, with
`id = reference_id`, `structure` the curated skl,
`options` the makeinput flags from the manifest, the
default (Imago) runner, and `key_fields` declaring
`kpoint_spec` + `convergence_threshold` + `imago_commit`
as scalars and the structure file as a key file (6.2.5).

#### 6.2.2 The pluggable runner seam

A *runner* is the seam (Principle 8) that isolates
kaleidoscope's dispatch core from how a unit actually
executes.  It is a small protocol:

```
Runner.run(unit, run_dir) -> RunOutcome
```

The runner receives a unit and the prepared run
directory, executes the calculation however it likes,
and returns a **domain-agnostic** `RunOutcome`:

```
RunOutcome
  ok        bool: did the unit complete (not "succeed
              scientifically" -- see detail)
  detail    short opaque string the runner chooses and
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
while converged / non-converged are runner-supplied
`detail` strings.  This is what lets kaleidoscope
surface convergence in `status.toml` *and* stay ignorant
of what convergence means.

- The **default runner** (`ImagoRunner`) calls the 6.1
  API: `run_structure(unit.structure, unit.options,
  run_dir)` (or `run_prepared` when inputs are already
  staged).  It maps the returned `ImagoResult` into a
  `RunOutcome`: `ok = status in {CONVERGED,
  NOT_CONVERGED, SKIPPED}` (the binary *ran*),
  `detail = status.name.lower()`.  It also **persists
  the full `ImagoResult` into the run directory** as
  `result.toml` (6.2.6), so the Imago-native detail
  survives for the client to reload without
  kaleidoscope ever parsing it.
- An **ASE runner** wraps `ImagoCalculator` (D12) for
  units that need ASE-MD or ASE-relaxation semantics; it
  too ultimately calls the 6.1 API underneath.
- A single campaign may **blend runners** per unit, so
  plain SCFs and adapter-wrapped calculations dispatch
  under one campaign (ARCHITECTURE 9.4).  New adapters
  slot in by implementing the protocol; the dispatch
  core never changes (Principle 8).

#### 6.2.3 Parsl dispatch and complete-and-report

Each unit becomes one Parsl app: a `python_app` that
runs `unit.runner.run(unit, run_dir)` on a worker.
Kaleidoscope's `parsl_config` (a Parsl `Config`, supplied
by the client/deployment) maps those apps onto SLURM via
a `HighThroughputExecutor` and a SLURM provider, so the
same dispatch code serves a laptop, an interactive node,
and a batch allocation -- only the `Config` changes.

The two workload shapes ARCHITECTURE 9.4 calls out are
both expressed in this one model:

- **Embarrassingly parallel sweeps** (thousands of
  independent SCFs): each unit is an independent app
  future; Parsl schedules them across the executor's
  workers.
- **Tightly iterative inner loops** (adaptive
  convergence, future AIMD): the *iteration* lives
  inside the unit's runner (it calls the 6.1 API in a
  loop, or drives ASE's optimizer), so kaleidoscope still
  dispatches one unit; it does not need to model the
  inner loop as a DAG.  If a future client genuinely
  needs cross-unit data flow, Parsl's own futures compose
  -- but that is not required by D13.

**Complete-and-report (Principle 10)** is the dispatch
core's contract.  Kaleidoscope gathers all futures and
**catches exceptions per future** rather than letting
one propagate: a unit whose app raised, or whose Parsl
task was lost cluster-side, is recorded with the
appropriate status (6.2.4) and the campaign continues.
No single unit failure aborts the batch.  When all
futures have resolved, kaleidoscope returns a
`CampaignReport` (6.2.6); deciding whether the aggregate
is scientifically acceptable is the client's job, never
kaleidoscope's.

#### 6.2.4 Workspace layout (resolves ARCHITECTURE 9.8)

This pins the strawman of ARCHITECTURE 9.6 into a
committed scheme.

```
<root>/
  campaign.toml          generated from the Campaign
                           (6.2.1): what to run.
  structures/<id>/        acquired/curated inputs.
  runs/<id>[/<calc>]/      one working dir per calc:
      <staged makeinput inputs + run outputs>
      cache_key.toml      identity snapshot (6.2.5).
      result.toml         runner-persisted native result
                            (6.2.6); Imago for ImagoRunner.
      status.toml         lifecycle + outcome (below).
  results/                client-aggregated outputs.
  logs/
```

**Stable-id convention.**  `<id>` is the client-supplied
stable per-structure key.  Kaleidoscope requires it to
be filesystem-safe and unique within the campaign:
lowercased, restricted to `[a-z0-9_-]`, with any other
character rejected at `Campaign` build time (not
silently rewritten -- a surprising rewrite would break
the cache hit-test, 6.2.5).  The producer uses the
curation `reference_id`; an acquisition campaign uses the
COD id.  Uniqueness collisions abort the build with the
two offending units named.

**`<calc>` tag format.**  The optional second level
exists only when one structure hosts more than one
calculation (different bases, a property run vs. its SCF).
A unit with `calc = None` runs directly in `runs/<id>/`
with no second level.  When present, `<calc>` obeys the
same `[a-z0-9_-]` rule and must be unique among the calcs
sharing an `id`.  If the client supplies no tag but an
`id` ends up hosting multiple units, kaleidoscope derives
a default tag from the runner's job identity (for the
Imago runner, `"<job_name>-<basis_scf>"`, e.g.
`"scf-mb"`), and errors only if that derived tag still
collides.

**`status.toml` schema.**  One file per run directory,
rewritten as the unit moves through its lifecycle:

```
id               = "<id>"
calc             = "<calc>"     # omitted when None
status           = "queued" | "running" | "done"
                   | "failed" | "lost"
detail           = "<runner string>"  # e.g. "converged"
runner           = "imago" | "ase" | ...
submitted_at     = <iso8601>
started_at       = <iso8601>    # omitted until running
finished_at      = <iso8601>    # omitted until terminal
runtime_seconds  = <float>      # omitted until terminal
message          = "<text>"
```

The five `status` values are kaleidoscope-owned and
generic.  `queued` / `running` are lifecycle;
`done` / `failed` are terminal runner outcomes
(`done` iff `RunOutcome.ok`); `lost` is the
kaleidoscope-only category for a Parsl-side
disappearance (worker died, allocation expired) where no
`RunOutcome` ever came back.  Convergence does **not**
appear as a status -- it rides in `detail`, per 6.2.2.

#### 6.2.5 The run-reuse cache

The cache is the general kaleidoscope mechanism of
ARCHITECTURE 9.6, split into mechanism (kaleidoscope) and
policy (client) so generality does not cost correctness.

**Mechanism (kaleidoscope).**  Before launching a unit,
kaleidoscope resolves its `run_dir = runs/<id>[/<calc>]/`
and performs the hit-test:

1. If the directory exists, holds a `cache_key.toml`
   that matches the unit's *current* key (below), and
   its `status.toml` reads `status = "done"`, the unit
   is a **hit**: skip the launch, and report the existing
   outcome straight from `status.toml` / `result.toml`.
2. Otherwise (no directory, key mismatch, or a
   non-`done` status) it is a **miss**: write a fresh
   `cache_key.toml`, set `status = "queued"`, dispatch,
   and update `status.toml` through the lifecycle.

Resuming a campaign is therefore *nothing more than
re-running it*: the hit-test over every unit naturally
skips the completed ones and re-dispatches the rest.

**The key has two parts**, mirroring the producer's
existing `is_cached_v2` (DESIGN 5.7) and generalizing it:

- **Scalar fields** -- written verbatim into
  `cache_key.toml` as TOML and compared field-by-field
  (the producer's `kpoint_spec`, `convergence_threshold`,
  `imago_commit`).
- **Key files** -- declared by name in `key_fields`;
  compared by **byte-comparison against the copy already
  staged in the run directory** (the producer's structure
  file).  This deliberately keeps DESIGN 5.7's
  "byte-compared structure file copies, no hashing, for
  debuggability" property: a developer can diff the files
  to see *why* a cache missed, which a hash would hide.

**Policy (client).**  The client supplies the key fields
in `CalcUnit.key_fields`; only it knows which inputs
define identity for its calculations.  Kaleidoscope never
guesses -- a too-broad key risks false hits and wrong
science; a too-narrow key risks needless re-runs.  This
mechanism subsumes the producer's bespoke
`share/atomicBDB/cache/scf/<reference_id>/`; C69 folds
that producer cache into this one.

The boundary with 6.1's checkpointing is clean and worth
restating: `imago.py` resumes *within* a run directory
(skip completed integrals, skip an already-done basis
SCF); kaleidoscope decides whether to *launch* the run
directory at all.  The two never overlap.

#### 6.2.6 Harvest handoff and the campaign report

Kaleidoscope returns a `CampaignReport`: one entry per
unit, each carrying `id`, `calc`, `status`, `detail`,
`run_dir`, `runtime_seconds`, and `message` -- exactly
the generic `status.toml` fields, nothing domain-specific.
An optional per-unit `on_outcome` callback (6.2.1) fires
as each unit reaches a terminal state, so a client can
stream-process rather than wait for the whole batch.

**Harvest stays on the client side.**  Kaleidoscope does
not read domain data out of run directories (Principle
9).  The handoff is the run directory itself: the runner
persisted its native result there (`result.toml`), so the
client walks the report and, for each unit it deems
acceptable, opens `run_dir` and reads what it needs.  For
the producer that means: keep units whose `detail ==
"converged"`, reload the 6.1.2 `ImagoResult` from
`result.toml`, read the converged `scfV` from
`result.outputs["scfV"]`, and pair its coefficients with
the input alphas (5.7 / ARCHITECTURE 9.7).  Non-converged
or failed units are simply skipped -- recorded in the
report, never harvested.

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
  kaleidoscope re-dispatch on the next campaign run via
  its non-`done` status) versus left for the client to
  re-launch.  The cache mechanism already makes a plain
  re-run safe; the open question is only whether to
  retry *eagerly*.
- **`campaign.toml` as an authoring surface.**  Whether
  it may be hand-written as the primary input rather than
  always generated from the in-process `Campaign`.
- **Concurrency limits for tightly-iterative units.**
  Whether such units need a distinct executor or a
  resource cap so a few long inner loops do not starve a
  parallel sweep sharing the same allocation.
- **`result.toml` for non-Imago runners.**  The Imago
  runner persists an `ImagoResult`; what a future
  non-Imago runner persists (and how a mixed-runner
  client reads it back) is that runner's contract, set
  when the runner is added, not here.

---

## References

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
- Eqs. 22-24: correction terms for improved accuracy

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

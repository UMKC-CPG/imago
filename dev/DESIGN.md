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

By contrast, a `1 1 1` mesh is *not* Gamma: it is one k-point
with the usual (auto or `-kpshift`) shift applied -- a single
shifted, mean-value sample -- and runs on the general complex
executable `imago`. The no-option default is likewise a single
shifted point, so prior behavior is unchanged; only the
explicit `0` sentinel is new. A density of 1.0 stays a genuine
density request, distinct from the `0` Gamma sentinel.

A zero mixed with positive mesh counts (e.g. `0 0 1`) is
rejected as a fatal typo for the all-zero sentinel.

The rule lives in `makeinput.reconcile` (per-group resolution
of `kp_gamma` and the display `kp_note`) and `_make_kp` (the
canonical Gamma write). It supersedes an earlier behavior in
which a `1 1 1` mesh was labeled `(Gamma)` while a nonzero
auto-shift was still written -- a mismatch that routed the job
to the complex executable despite the label.

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
  preferred  bool          Optional, default false.  Marks
                           the one record per matcher family
                           that the consumer's file-dictated
                           (crystalline) match uses (5.6.5
                           step 2).  Exactly one record per
                           present family carries `preferred
                           = true` (rule 10); the preferred
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
  list length).  The neighbor multiset is element-only --
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
   registered in `matchers.py` at load time (8.9).
   Unknown methods are a hard error rather than a
   silent skip, so that a typo in the manifest fails
   loudly rather than quietly omitting the fingerprint
   from the lookup.
10. For each `method` that appears in the file's
    fingerprint records, exactly one record carries
    `preferred = true` -- the record the consumer's
    file-dictated match uses (5.6.5 step 2).  Zero (a
    family present but none preferred) or two-or-more
    preferred for one `method` is a hard error.  The
    database-wide constraint that all elements' preferred
    records of a family share one `sub_spec` (manifest
    rule 11) is cross-file and cannot be checked from a
    single per-element file; the loader trusts the
    producer to have written a consistent database, in
    keeping with VISION Principle 5 (the database is
    regenerated from the manifest, never hand-edited).

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

**Dedup on insert, weight by multiplicity.**  Harvesting
is therefore an *insert-or-merge*.  When a harvested
environment duplicates one already stored, the producer
does not append a new entry; it increments a
**multiplicity** weight on the existing one (how many
atoms, across how many models, have collapsed into it).
The database grows with environment *diversity*, which is
sub-linear in atoms and plateaus as the space fills in.
The dedup tolerance is the producer-side mirror of the
consumer-side similarity floor (5.6.5, TODO C61): the
consumer asks "is this query close enough to a stored
environment to *import* its potential?", and the producer
asks "is this new environment close enough to a stored one
that storing it *adds nothing*?" -- one tolerance, two
uses.

**Dedup subsumes the symmetry gate.**  No special case is
needed for crystals.  A symmetry-assigned type's atoms are
identical, so they dedup to a single entry automatically
(multiplicity equal to the type's size) -- exactly "one
representative per type."  A disordered type's atoms vary,
so they dedup to however many genuinely distinct
environments exist, generally far fewer than the atom
count.  The crystalline case is just the degenerate
collapse of the same rule.

**Dedup must be conservative -- the union of per-method
distinctness.**  An environment is a duplicate only when
it is close under *every* method at once; if it is novel
under *any* method it must be kept.  This matters because,
within a reduce-assigned type, all atoms share one reduce
fingerprint but carry *different* bispectrum witnesses:
deduping by the reduce metric alone would collapse the
type and discard that bispectrum coverage irrecoverably.
Merging only on an all-methods duplicate keeps every
method's index complete, and the multiplicity counts only
true all-methods duplicates.

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
exactly the deduplicated, weighted, coverage-oriented
corpus the dedup model produces:

- At inference the network pays no per-query search cost;
  the lookup is amortized into its weights.  So the search
  concern of 5.2.3 does not apply to it.
- Its *training* is harmed by raw near-duplicates: they
  teach nothing new and *bias* the model toward
  over-represented environments (a network trained on raw
  per-atom data would be dominated by bulk-like sites and
  under-learn the rare defect environments that matter
  most).  A deduplicated set with multiplicity available
  as a *sample weight* -- rather than baked in as
  duplication -- is the correct corpus.
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

**Schema implications (deferred to implementation).**  The
native/witness model adds the `type_assignment` provenance
field (already listed in 5.2's provenance table), and the
dedup model adds a per-entry `multiplicity` integer
(default 1, still to be threaded through 5.2, 5.4, 5.5,
and 5.7).  The per-atom-pairing safeguard costs nothing
until normalization, when the origin label is carried onto
split records.  As an implementation status: the present
producer -- the C60 harvest -- emits neither
`type_assignment` nor `multiplicity` and writes one entry
per declared site, so the loader does not yet enforce
them; the schema above specifies the target.  The
insert-or-merge harvest, the conservative dedup rule, and
the validation that enforces both fields are scheduled
work (TODO C88), not part of the current bispectrum half.

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
   the whole run by `-nofingerprint`, which sends every species
   straight to step 3).  Match is a single best-effort lookup,
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

### 5.7 Regeneration Pipeline Algorithm

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
`scf_threshold`) plus a stable
`reference_id`; per-entry harvest declarations
(`atom_site`, expected `element`, `description`,
and an optional `label` override -- absent labels
are derived at harvest per 5.2.1).

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
the authoring tool share one definition: the
producer imports it to *run* a manifest,
`expand_manifest.py` imports it to *write* one, and
the library itself depends only on the lower
libraries it validates against
(`initial_potential_db`, `guidance_db`).  A curator
does not hand-write a complete manifest from
nothing: `cod_fish.py` discovers and pins structures
and prints a complete sketch -- a `schema_version`
header plus one `[[reference_solid]]` stub per
structure, each with a `reference_id` auto-derived
from the CIF metadata
(`<formula>_<symbol>_<number>_<year>`);
`expand_manifest.py` reads that sketch and fills in
the shared method defaults and the per-structure
harvest curation, either interactively or by
stamping the defaults and emitting a fill-in
template.  Each sketch stub also carries two
discovery hints cod_fish read from the CIF -- the
composition (``elements``) and a
``source_description`` -- which the interactive flow
offers as the element and description defaults, so
the curator invents neither; they are not schema
fields (the producer ignores them and the finished
manifest omits them).  `cod_fish.py` stays a pure
discovery tool and never writes a manifest.  The writer emits
human-readable TOML -- shortest round-trippable
floats, inline `sub_spec` tables in their authored
order, `label` only when present and `preferred`
only when true -- and its output round-trips through
`load_manifest_v2`.

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

- A curation manifest (TOML; schema v2 specified below).  Default
  location is `share/atomicBDB/manifest.toml`; alternate via the
  `--manifest` flag.
- The historical guidance dataspace (`share/historicalGuidanceDB/`,
  DESIGN 7).  Consulted per reference solid to predict the
  converged k-point density and to size the verification grid
  around that prediction.  When the dataspace cannot predict for a
  solid (too few neighbors), the producer falls back to the wide
  default grid defined in the flight builder (DESIGN 6.2.8).
- The existing `share/atomicPDB/` tree (for `pot1` and `coeff1`
  reads when refreshing each element's `"isolated"` baseline
  entry).
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

[[reference_solid]]
reference_id          = "au_fcc"
system_type           = "crystalline"
# Exactly one of (cod_id, cod_revision) or structure_path is set
# per [[reference_solid]] (validation rule 4 below).
cod_id                = 9008463
cod_revision          = "2023-04-12"
# structure_path      = "au_fcc.skel"     # alternative form
kpoint_spec           = { density = 60.0, shift = [0.0, 0.0, 0.0] }
scf_threshold         = 1.0e-6
basis                 = "fb"
functional            = "wigner"
kpoint_integration    = "linear-tetrahedral"

  [[reference_solid.entry]]
  element     = "Au"
  atom_site   = 1
  label       = "default_solid"
  default     = true
  description = "Au in fcc bulk (Fm-3m)."

    # preferred = true is the record the consumer uses for a
    # file-dictated match; one preferred sub_spec per family,
    # uniform database-wide (see the convention note below).
    [[reference_solid.entry.fingerprint]]
    method    = "bispectrum"
    sub_spec  = { twoj1 = 8, twoj2 = 8 }
    preferred = true

    # A non-preferred alternate sub_spec is allowed (stored for
    # validation comparison); it just may not also be preferred.
    [[reference_solid.entry.fingerprint]]
    method   = "bispectrum"
    sub_spec = { twoj1 = 6, twoj2 = 4 }

    [[reference_solid.entry.fingerprint]]
    method    = "reduce"
    sub_spec  = { level = 2, thick = 0.5, cutoff = 5.0,
                  tolerance = 0.05 }
    preferred = true

[[reference_solid]]
# ... another solid ...
```

**Per-solid fields.**

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

**Per-entry fields (`[[reference_solid.entry]]`).**

- `element` (string): element symbol the entry contributes to.
  Cross-checked against the species at `atom_site` after Imago
  loads the structure.
- `atom_site` (positive integer): 1-based site index into the
  structure, matching Imago's site-indexing convention everywhere
  else in the codebase.
- `label` (string, *optional*): the label this entry is written
  under in the element's `s_gaussian_pot.toml`.  When omitted, the
  producer derives it at harvest per 5.2.1
  (`<reference_id>-<element><species>-t<type>-a<site>`); when
  present it overrides that derived default.  An explicit
  `(element, label)` must be unique across the entire manifest;
  a derived label's `(reference_id, element, atom_site)` must be
  unique (rule 6).
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
  to the matcher registry in `matchers.py` (8.9).
- `sub_spec` (inline table): method-specific parameters.
  Two declarations on the same entry with the same `method`
  and the same `sub_spec` keys-and-values are a hard error
  (rule 8); same `method` with differing `sub_spec` is
  explicitly allowed -- the database stores as many sub_specs
  per family as the curator wants.
- `preferred` (bool, *optional*, default false): marks the one
  record per matcher family the consumer uses for a
  file-dictated (crystalline) match (5.6.5 step 2).  Exactly
  one declaration per `method` is preferred for each element
  that carries any fingerprint of that family, and the
  preferred `sub_spec` for a family is *uniform across the
  whole database* -- declared once and stamped onto the
  matching record in every element file.  Storing additional
  non-preferred sub_specs of the same family is always allowed;
  what the validator forbids is a second preferred record, or a
  preferred `sub_spec` that diverges from the database-wide one
  for that family.

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
   `system_type`, `basis`, `functional`, `kpoint_integration`,
   `kpoint_spec`, `scf_threshold`, and *exactly one* of
   `{(cod_id, cod_revision), structure_path}` (see rule 4 for
   details).  `system_type` must be one of the four allowed
   values `{"crystalline", "amorphous", "nanostructure",
   "molecular"}`; any other value is a hard error, since the
   guidance predictor (DESIGN 7) switches its sub-model on it and
   the produced entry records it for forensics.  `basis`,
   `functional`, and `kpoint_integration` are likewise required:
   they select the predictor sub-model (DESIGN 7.6) and are
   recorded on every produced entry's context, so nothing the
   producer emits depends on an implicit default (VISION
   Principle 5).
3. Every `[[reference_solid.entry]]` carries `element`,
   `atom_site`, `default`, `description`.  `label` is *optional*
   (5.2.1): when present it is an explicit curator override; when
   absent the producer derives it at harvest from the run's site
   identity, so the species and type numbers (unknown until the
   grouping pass runs) need not be authored ahead of time.
4. Exactly one of `cod_id` or `structure_path` is set on each
   `[[reference_solid]]`.  If `structure_path`, it resolves to
   an existing file under the manifest's directory.  If `cod_id`,
   it parses as a positive integer *and* `cod_revision` is
   present as a non-empty string.
5. `reference_id` is unique across the manifest *and* label-safe
   (lowercase letters, digits, `-`, `_`): it is embedded verbatim
   in every derived entry label and typed into `-pot`, so a
   non-conforming id is a hard error (5.2.1).
6. No two entries may produce the same database entry.  For an
   entry with an explicit `label`, `(element, label)` is unique
   across the manifest — two solids cannot both produce, e.g.,
   `("Au", "default_solid")`.  For an entry with a derived label,
   `(reference_id, element, atom_site)` is unique, since the
   derived label is built from exactly those three; two such
   entries would collide on the identical label.  Either way a
   silent overwrite that would mask a curation mistake is
   refused.
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
10. For every `(element, method)` that appears with any
    fingerprint declaration across the manifest, *exactly one*
    of those declarations carries `preferred = true`.  Zero
    preferred (a family present but none marked) or two-or-more
    preferred for the same `(element, method)` is a hard error.
    This mirrors rule 7 for the entry-level `default` tag:
    the manifest is the single source of truth for which
    record the consumer picks.
11. The preferred `sub_spec` for a family is *uniform across
    the whole manifest*: all `preferred = true` declarations
    sharing a `method` must carry the identical `sub_spec`,
    regardless of element.  A preferred record whose `sub_spec`
    diverges from the family's database-wide preferred one is a
    hard error naming both.  (Non-preferred declarations may
    use any `sub_spec`; only the preferred ones are constrained
    to agree, which is what guarantees a single loen run per
    structure -- DESIGN 5.6.5 step 2.)

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

Earlier drafts of this section gave the producer its own
per-solid SCF cache under `share/atomicBDB/cache/scf/`, with an
`is_cached(ref, imago_commit)` check that compared a snapshot of
the SCF inputs plus a byte-for-byte copy of the structure file.
That machinery is **removed**.  Now that the producer delegates
SCF to kaleidoscope, the avoid-recompute responsibility moves to
kaleidoscope's run-reuse cache (DESIGN 6.2.5), which keys each run
on the structure file's contents together with the makeinput
options and the Imago build identity.  The properties the old
producer cache provided are preserved by that one mechanism:

- *Edits to harvest declarations stay cheap.*  Adding or changing
  a `[[reference_solid.entry]]` changes neither the structure nor
  the flight options, so the run is a cache hit and only the
  harvest step (`extract_potential`, PSEUDOCODE 11.4) re-runs.
- *Content, not path.*  The cache compares structure file
  *contents*, so renaming a `structure_path` file on disk does
  not force a re-run.
- *Build- and threshold-sensitive.*  Changing `scf_threshold`,
  the k-point grid, or the Imago build identity invalidates the
  cached run, as it must.

The reason the old design preferred direct comparison over a
content hash still holds — at ~100 reference solids run by hand,
naming the exact field that changed beats a bare "different hash"
— and DESIGN 6.2.5 carries that reasoning for the kaleidoscope
cache.

**Procedure.**

The pipeline runs in three phases — *build*, *dispatch*,
*harvest* — so that every reference solid's verification runs are
launched as one flat parallel batch (the predict-then-verify
shape of DESIGN 6.2.1) rather than one solid at a time.

1. Load and validate the manifest (the rules above).
2. For every element with a directory in `share/atomicPDB/`,
   refresh (or create) the `"isolated"` entry of that element's
   `s_gaussian_pot.toml` directly from the current `pot1` and
   `coeff1` files.  Guarantees the baseline is always present
   (rule 6 of 5.2) and tracks any changes in atomSCF output.
3. **Build.**  For each `[[reference_solid]]` in the manifest:
   a. `materialize_structure(ref)` → a local structure file
      (read from disk for `structure_path`, or fetched once from
      COD at the pinned `cod_revision` for `cod_id`).
   b. Build the two inputs the next step needs.  `options =
      make_producer_options(ref)` translates the manifest physics
      into each tool's coded settings (`functional` → `xccode`,
      `kpoint_integration` → `scfkpint`, `basis` → `scf_basis`,
      `scf_threshold` → `converg`, `kpoint_spec.shift` →
      `kpshift`; 6.2.10), while `submodel = {basis, functional,
      kpoint_integration}` keeps the human names the predictor and
      record speak.  Then `build_kpoint_convergence(structure,
      options, dataspace, system_type, submodel, …)` (DESIGN
      6.2.8 / 7) consults the guidance dataspace and returns the
      **verification grid**: a set of k-point densities bracketing
      the predicted converged density, widened by the prediction's
      uncertainty (a length-1 grid in trust mode; the wide default
      grid when the dataspace cannot predict).  A manifest
      `kpoint_spec.density` is passed as the `center` argument -- a
      curator override that pins or centres the grid and bypasses
      the predictor.
   c. Emit one `CalcUnit` per grid point — each carrying the
      materialized structure, the coded run settings (the SCF
      convergence limit `converg` among them), and the k-density
      for that point — plus one structure-only
      `imago.py -loen -scf no` unit per declared Fortran-side
      fingerprint (the bispectrum fingerprint depends on geometry
      alone, not on SCF convergence, so it need not wait for the
      converged grid point).  Tag the `<calc>` level by k-density
      per DESIGN 6.2.4 and collect every solid's units into a
      single `Flight`.
4. **Dispatch.**  Hand the whole `Flight` to kaleidoscope, which
   runs and tracks the batch through the wingbeat seam and its
   run-reuse cache (DESIGN 6.2.5).  The producer runs no SCF or
   loen calculation itself.
5. **Harvest.**  For each `[[reference_solid]]`:
   a. Pick the converged grid point from the solid's units with
      the two-sided delta-below-threshold rule (DESIGN 7.8 / the
      C72 harvest).  If no grid point converges (for example,
      non-convergence at the top of the range), skip the solid's
      entries and flag it in the run log rather than harvesting a
      non-converged potential.
   b. Record that run's SCF iteration count and convergence
      metrics in the run log.
   c. For each `[[reference_solid.entry]]`:
      i.   `extract_potential` for the named `atom_site` from the
           converged run.  The converged `scfV` output lists every
           potential type (a `NUM_TYPES` header + per-type Gaussian
           blocks under a `TOTAL__OR__SPIN_UP` channel); the harvest
           selects the site's type block -- its type number read
           from `datSkl.map` (ARCH 9.7) -- and takes each term's
           coefficient and alpha (columns 1 and 2) together.
      ii.  For each `[[reference_solid.entry.fingerprint]]`
           declaration, compute the fingerprint at `atom_site`
           for the requested `(method, sub_spec)`.  Python-side
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
           `FingerprintRecord` (5.4) and attach it to the
           entry-in-progress.
      iii. Construct a `PotentialEntry` (5.2) with the
           manifest-supplied `label`, `default`, and
           `description`, the run-supplied numerical fields, the
           run-supplied provenance (which records the solid's
           `system_type` for forensics, per rule 2), and the
           `FingerprintRecord` list assembled in step ii.
      iv.  Insert the entry into the in-memory `ElementDatabase`
           for its element.  If an entry with the same label
           already exists, replace it.
   d. **Guidance contribution.**  Harvest the same converged grid
      point into the historical guidance dataspace's staging area
      (`share/historicalGuidanceDB/staging/<system_type>/`, via
      the C72 `harvest_flight` hook), so every reference solid the
      producer converges becomes training data that sharpens the
      predictor for the next solid.  Trust-mode (length-1 grid)
      runs harvest the potential but do *not* auto-stage a
      guidance entry — a single point is weaker evidence than a
      converged grid (DESIGN 7).
6. Save each affected `ElementDatabase` to disk via
   `initial_potential_db.save()` (5.5).
7. Write `share/curation/run_log.toml` capturing the manifest
   snapshot, per-run iteration counts, the converged k-density
   chosen for each solid, and the Imago commit.  The validation
   harness (5.8) reads this log.

**Flags:**

- `--force`: forward a cache-bypass to kaleidoscope so every
  dispatched unit re-runs from scratch instead of reusing the
  run-reuse cache (DESIGN 6.2.5).  Fresh results are still written
  into that cache afterwards, so the next ordinary run is a
  warm-cache hit.
- `--manifest PATH`: alternate manifest location (default:
  `share/atomicBDB/manifest.toml`).
- `--element ELEM`: restrict regeneration to a single element's
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
   a provisional `imago.dat`.  The potential each atom
   receives is irrelevant -- the bispectrum is geometric
   -- so the per-element default-tagged entry (5.6.5
   step 3) is fine and no `-pot` is needed.  The
   `LOEN_INPUT_DATA` block of this `imago.dat` carries the
   matcher's `(twoj1, twoj2, ...)` parameters via
   `BispecMatcher.to_loen_input` (5.10.5).
2. **Run loen.**  Invoke `imago.py -loen -scf no` against
   that `imago.dat`, producing `fort.21` (5.10.3).
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

There is no self-invocation and nothing to guard against
recursion: each step is a separate, ordinary process the
orchestrator runs in order, and the only state passed
between the two makeinput runs is the skeleton file
itself.

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
   producer can fill the provenance fields of 5.2 and
   so kaleidoscope can form its run-reuse cache key
   (`kpoint_spec` + `scf_threshold` +
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
  measured          MeasuredQuantities | None: scalar
                      electronic-structure quantities
                      harvested from the converged SCF
                      output (see below).  None for
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
`key_fields` (scalar `scf_threshold` and
`imago_commit`; the structure file as a key file); they
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
  API.  It first **partitions the unit's options** into the
  makeinput-side and imago-side settings (6.2.10), then builds
  the run directory with `makeinput.build_run_dir` and runs it
  with `imago.run_prepared` (or `run_prepared` alone when the
  inputs are already staged).  It maps the returned
  `ImagoResult` into a `WingbeatOutcome`: `ok = status in
  {CONVERGED, NOT_CONVERGED, SKIPPED}` (the binary *ran*),
  `detail = status.name.lower()`.  It also **persists
  the full `ImagoResult` into the run directory** as
  `result.toml` (6.2.6), so the Imago-native detail
  survives for the client to reload without
  kaleidoscope ever parsing it.
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
  -- but that is not required by D13.

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
```

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

Resuming a flight is therefore *nothing more than
re-running it*: the hit-test over every unit naturally
skips the completed ones and re-dispatches the rest.

**The key has two parts**, mirroring the producer's
existing `is_cached_v2` (DESIGN 5.7) and generalizing it:

- **Scalar fields** -- written verbatim into
  `cache_key.toml` as TOML and compared field-by-field
  (the producer's `kpoint_spec`, `scf_threshold`,
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
                               #   converg, kpshift, imago_commit,
                               #   ... held fixed across the grid and
                               #   copied verbatim into every
                               #   CalcUnit.  Carries NO physics-name
                               #   keys -- the wingbeat forwards it
                               #   to the tools as-is (DESIGN 6.2.10).
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

*The routing rule (three buckets).*  The wingbeat splits a
unit's options into:

- **imago run options** -- the fixed, known imago key set
  `{job, edge, scf_basis, pscf_basis, serialxyz, valgrind}`
  (6.1), exported as `imago.OPTION_KEYS` so the wingbeat does not
  hard-code it -> handed to `from_options`.
- **bookkeeping / cache-only** -- the keys in `CACHE_ONLY_KEYS`
  (today just `imago_commit`, the build identity that busts the
  run-reuse cache across imago versions, 6.2.5), exported from
  `kaleidoscope.wingbeats` because this is kaleidoscope's own
  bookkeeping, not an imago tool input.  Such a key is *dropped
  before forwarding* and reaches neither tool.  This is safe
  because the cache identity is captured
  separately in `unit.key_fields` at build time
  (`standard_key_fields` copies the scalars out of `options`), so
  removing the key from the forwarded options does not change the
  cache key.
- **makeinput build options** -- everything else, handed to
  `build_run_dir`.  makeinput **keeps its strict unknown-key
  check**, which now serves as the typo backstop: a key that is
  neither an imago key nor bookkeeping nor a real makeinput dest
  still raises, so the safety that strictness buys is preserved.

| Producer emits | Bucket | Tool dest / note |
| --- | --- | --- |
| `scf_basis` | imago | basis selection, `fb -> 2` |
| `kpd` | makeinput | k-point density (already correct) |
| `xccode` | makeinput | XC functional code (wigner = 100) |
| `scfkpint` | makeinput | k-point integration (LAT = 1) |
| `kpshift` | makeinput | gamma-centred mesh offset |
| `converg` | makeinput | SCF threshold (new dest; below) |
| `imago_commit` | dropped | cache identity only (6.2.5) |

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
scalar (6.2.5), so `standard_key_fields`' `_KEY_SCALAR_NAMES`
becomes `("converg", "imago_commit")` -- the cache keys on the
dest the run actually used.

*Related robustness fix (surfaced by the failed smoke run).*
When every unit fails this way, no `result.toml` is written, and
the harvest's `pick_converged_unit` raised `FileNotFoundError`
trying to open it.  A failed or result-less unit must be treated
as **non-converged** (logged and skipped, 5.7), never an
uncaught crash: the harvest reads `status.toml` and skips any
unit whose status is not a completed run before it reaches for
that unit's `result.toml`.

*Follow-on code (for TODO).*  (a) add makeinput `-converg`;
(b) rewrite `make_producer_options` to emit the dest-keyed, coded
vocabulary above; (c) export `imago.OPTION_KEYS` plus a
`CACHE_ONLY_KEYS` set from `kaleidoscope.wingbeats`; (d) move the
partition into `ImagoWingbeat.run` and retire the single-shared-
options call to `run_structure`; (e) update `_KEY_SCALAR_NAMES`;
(f) harden `pick_converged_unit` against a missing `result.toml`.

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
- `memory_per_node` / `memory_per_worker` -- guards so a run
  neither overflows memory nor under-requests it.

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
all.  The everyday path is a single command (captured in the
`command` log the scripts already keep); for a reproducible
record, the client may also write the resolved configuration
as a human-readable file in the run directory, beside the
manifest.

The command-line default is `slurm-per-job`, because the whole
point of this work is that the producer and the database seed
*should* reach the scheduler rather than run on a login node
(ARCHITECTURE 9.7): on a cluster, dispatching one scheduler
job per unit is the right thing to do with no flags at all.
A run launched where no site settings file is present is
therefore a configuration error, surfaced up front, rather
than a quiet fall-back to a serial local run -- the cluster
behaviour is the default, and a local run is the deliberate
opt-out.  That opt-out is `--dispatch local`, which needs no
site facts and builds no `Config` at all; the test suite, a
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
    to, reused as the entry's `metric_threshold` (the v1
    convention).
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
      b. For each CalcUnit's converged run, parse
         result.toml for total_energy, gap_ev, gap_kind,
         total_magnetization, and scf_threshold; read the
         swept k-density out of the CalcUnit's calc tag.
      c. Pick the converged grid point: the smallest
         k-density at which |E_i - E_{i+1}| <
         metric_threshold AND |E_i - E_{i-1}| <
         metric_threshold for i in (1, len(grid)-1).
         Stricter than the original "single delta below"
         rule: requires both consecutive-pair deltas to
         be small, mitigating a single-grid-point
         numerical fluke.
      d. If no point satisfies the criterion (energy
         still moving at the top of the grid), log a
         warning, tag the flight with
         `prediction_mismatch = true`, and SKIP this
         structure.  Non-converged sweeps do not earn
         an entry.  The user must widen the grid and
         re-run.
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
              array gathered in step b, so the
              auto-promote rule can judge flatness from
              the staged file alone),
              converged_at = chosen k-density,
              metric/metric_threshold,
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

**Promotion** (`guidance_promote.py`):

A curator helper.  Four modes of operation:

- *Interactive review* (default).  Lists all files
  under `staging/<system_type>/`; for each, prints the
  signature, measured quantities, verification grid,
  and provenance, and asks the curator to PROMOTE,
  SKIP, or DELETE.  Promoted files move to
  `entries/<system_type>/`.  Skipped files stay in
  staging for later review.  Deleted files are
  removed.
- *Auto-promote rule* (`--auto-promote`).  Promotes
  every staging file that satisfies an objective
  acceptance test:
  - The converged k-density landed in the middle 60%
    of the verification grid (not at either endpoint).
    A converged-at-endpoint result is suspicious: the
    grid may not have been wide enough.
  - The total-energy variance over the top three grid
    points -- read from the entry's `grid_energies`
    array (7.2), which is why harvest records it -- is
    below `metric_threshold * 10` (the converged region
    is convincingly flat, not just one delta below
    threshold).  A staging entry that lacks
    `grid_energies` (a hand-written manual entry) is
    never auto-promoted on this criterion; it falls to
    interactive review.
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
the DESIGN 5.7 regeneration discipline.

### 7.10 Open Design Questions

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
ARCHITECTURE 6.5-6.7. It transcribes the block-cyclic scheme
designed in the sibling upolcao branch onto imago's terms,
records the work-assignment decision that scheme rests on,
and specifies the grid load balancer and parallel-I/O
alignment that are already proven. Backend choices that
remain open -- which distributed eigensolver, how device
placement is expressed -- are collected in 9.6 rather than
committed here, consistent with VISION's intent to keep
every parallel axis open.

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
and the lowest-risk parallelism imago can adopt; it is the
recommended first increment.

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

### 9.5 Parallel HDF5 Alignment

The distributed matrices are written to and read from HDF5
collectively. For compression and write efficiency the on-
disk chunk size should align with the block-cyclic block
size, so each rank's write touches whole chunks rather than
splitting them. The exact collective-write pattern against
compressed chunks needs measurement -- one chunk per block
per rank, versus larger chunks filled by several ranks'
collective contributions -- and is flagged as an open
calibration in 9.6.

### 9.6 Open Design Questions

- **Distributed eigensolver backend.** ScaLAPACK `PZHEGVX`
  (the interface upolcao declared) versus ELPA (a faster
  two-stage solver with a single API across CPU and GPU)
  versus a GPU vendor solver. ELPA is the leading candidate
  because it reuses the block-cyclic layout of 9.3 and
  unifies the CPU/GPU axis behind one call, but the choice
  is deferred to a benchmark on representative secular
  dimensions. See ARCHITECTURE 6.6.
- **Device-placement expression in Fortran.** How the per-
  kernel CPU/GPU boundary (VISION Principle 14) is expressed
  -- OpenACC, OpenMP target, CUDA Fortran, or a library
  boundary such as ELPA -- is open and will likely differ
  per kernel.
- **Parallel-HDF5 chunk/block strategy.** The collective-
  write pattern against compressed chunks (9.5) is settled
  by measurement, not assumed.
- **Replicate-and-broadcast retirement.** The order in which
  the interim replicate-and-broadcast paths (ARCHITECTURE
  6.5) are replaced by genuine distribution -- grid work
  first, then integral assembly, then the solve -- and the
  validation gate between each stage, is sequenced here as
  implementation begins.

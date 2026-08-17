# Pseudocode

> **Document hierarchy:** VISION -> ARCHITECTURE -> DESIGN
> -> **PSEUDOCODE** -> Code. For the design rationale behind
> these algorithms, see `DESIGN.md`.

---

## 1. Generate Tetrahedra (DESIGN 1.2)

Every box is cut once per long diagonal, and each resulting
tetrahedron carries an equal share of the weight.  DESIGN 1.2
gives the reason: six tetrahedra sharing ONE diagonal is not
a decomposition the crystal point group carries onto itself,
and the average of all four cuts is.

`numTetraDiagonals` says how many to use.  It comes from the
k-point file, defaults to 4, and may be 1 for the cheaper
single-cut decomposition, which is safe for a totals-only run
and for the SCF occupation path (DESIGN 1.2).  The routine
below is written once and takes the count as data; it must
not branch on 1 versus 4.

The share of the weight is not written anywhere.  `tetraVol`
is `1 / numTetrahedra` as before, so it follows the count
automatically and every consumer stays unchanged.  Nothing
downstream needs to know how many diagonals were used, and
nothing downstream should ask.

```
# A box corner is a triple of bits. The diagonal from a
#   corner to its opposite is walked one coordinate at a
#   time; the four corners visited form a tetrahedron, and
#   the six orders of the three coordinates give the six
#   tetrahedra of that diagonal. Starting from (0,0,0) this
#   reproduces the M1..M8 table of DESIGN 1.2 exactly, which
#   is the check to run first when transcribing this.
# The four diagonals, named by the corner each starts from.
#   The opposite ends are (1,1,1), (0,1,1), (1,0,1) and
#   (1,1,0), so these four cover every diagonal once. Taking
#   the first `numTetraDiagonals` of them means the single-cut
#   setting reproduces the historical decomposition exactly,
#   which is what makes a one-diagonal run comparable with
#   anything computed before this section existed.
DIAGONAL_STARTS = [(0,0,0), (1,0,0), (0,1,0), (0,0,1)]

function generateTetrahedra(nA, nB, nC, numTetraDiagonals):
    numTetrahedra = 6 * numTetraDiagonals * nA * nB * nC
    allocate tetrahedra(4, numTetrahedra)
    t = 0
    for a = 1 to nA:
        for b = 1 to nB:
            for c = 1 to nC:
                for start in the first numTetraDiagonals
                        entries of DIAGONAL_STARTS:
                    opposite = (1,1,1) - start
                    for (first, second) in the six ordered
                            pairs of distinct axes:
                        second_corner = start with
                              coordinate `first` changed
                        third_corner = second_corner with
                              coordinate `second` changed

                        t = t + 1
                        tetrahedra(:, t) = [
                            corner(a, b, c, start),
                            corner(a, b, c, second_corner),
                            corner(a, b, c, third_corner),
                            corner(a, b, c, opposite)]

function corner(a, b, c, offset):
    # One box corner, with periodic wrapping.
    return idx(a + offset(1),
               b + offset(2),
               c + offset(3))

function idx(a, b, c):
    # Periodic wrapping, 1-based indexing
    return getIndexFromIndices(
        mod(a-1, nA) + 1,
        mod(b-1, nB) + 1,
        mod(c-1, nC) + 1)
```

**Check on the construction, not on the physics.**  Three
properties hold for any mesh and each fails loudly if the
transcription is wrong, so test them before testing any
spectrum: the number of DISTINCT tetrahedra is
`6 * numTetraDiagonals * nA * nB * nC`; no tetrahedron has a
repeated mesh point; and every mesh point belongs to the same
number of tetrahedra (24 per diagonal, so 96 at the default
on a 4x4x4 mesh).  The reference implementation of these
checks is `jobs/knbo3/o9_pdos/tetra_symmetry_remedies.py`,
which also verifies the generic diagonal walk against the six
explicit assignments of DESIGN 1.2.

**And one check on the setting itself.**  With
`numTetraDiagonals = 1` the routine must reproduce the
historical decomposition tetrahedron for tetrahedron, since
it is then walking the same M1-M8 diagonal by a general rule
instead of an explicit table.  That equality is the cheapest
possible test of the rewrite and should be the first one run,
because it separates a mistake in the general construction
from a change of physics: if it fails, the walk is wrong, and
no result from the four-diagonal setting means anything until
it passes.

---

## 2. LAT TDOS (DESIGN 1.3)

The TDOS at each energy grid point is the sum of per-
corner DOS weights from `bloechlCornerDOSWt` (section
2a), summed over all bands and tetrahedra. The per-corner
weights are also used individually by the PDOS (section
8.3).

```
function computeTDOS_LAT(eigenValues, tetrahedra,
        numTetrahedra, tetraVol,
        energyGrid, numEnergyPoints,
        numStates, numSpins,
        fullKPToIBZKPMap):
    allocate tdos(numEnergyPoints, numSpins) = 0.0

    for spin = 1 to numSpins:
        for n = 1 to numStates:
            for T = 1 to numTetrahedra:
                # Map full-mesh corners to IBZ
                # eigenvalues.
                for c = 1 to 4:
                    kFull = tetrahedra(c, T)
                    kIBZ = fullKPToIBZKPMap(kFull)
                    eps(c) = eigenValues(
                        n, kIBZ, spin)

                # Sort eigenvalues ascending.
                sortedEps = sort(eps)

                for iE = 1 to numEnergyPoints:
                    E = energyGrid(iE)
                    if E < sortedEps(1) or
                            E >= sortedEps(4):
                        cycle

                    # Per-corner DOS weights. The
                    # TDOS uses only their sum.
                    cornerDOSWt_LAT(1:4) =
                        bloechlCornerDOSWt(
                            E, sortedEps)

                    tdos(iE, spin) +=
                        sum(cornerDOSWt_LAT)
                        * tetraVol / spin
                        / hartree

    # Diagnostic: integrated area should equal
    # the number of spin states in the energy
    # range. Use deltaDOS * hartree because
    # deltaDOS is in Hartree but the TDOS is
    # in states/eV.
    integratedArea = trapezoid(tdos)
        * deltaDOS * hartree

    return tdos
```

---

## 2a. Bloechl Corner DOS Weights (DESIGN 1.3)

`bloechlCornerDOSWt` computes the per-corner DOS
density weights `cornerDOSWt_LAT(1:4)` for one
tetrahedron at energy E. These are the energy
derivatives of the cumulative corner integration
weights `cornerIntgWt_LAT` (section 3a):

  cornerDOSWt_LAT(c) = d/dE [ cornerIntgWt_LAT(c) ]

Their sum equals the total per-tetrahedron DOS
(the `bloechlDOS` value from the original TDOS
implementation). This identity provides a built-in
self-consistency check.

The derivation follows from the product rule applied
to the cumulative weight expressions in section 3a.
For each case, we reuse the same intermediate
variables (t_j, s_j, a, b, c, d, v_I, v_II, v_III)
and also compute the total DOS `gTotal`, then
decompose it across the four corners.

```
function bloechlCornerDOSWt(E, eps):
    # eps = [e1, e2, e3, e4] sorted ascending.
    # Returns cornerDOSWt_LAT(1:4): the per-corner
    # DOS density weights at energy E.
    #
    # cornerDOSWt_LAT(c) is the spectral density
    # (units: 1/energy) attributed to sorted
    # corner c. Their sum equals the total DOS
    # per unit BZ volume for this tetrahedron.

    e1 = eps(1);  e2 = eps(2)
    e3 = eps(3);  e4 = eps(4)
    tol = 1.0e-12

    # Case 0: outside eigenvalue range.
    if E < e1 or E >= e4:
        return [0, 0, 0, 0]

    # ------------------------------------------
    # Case 1: e1 <= E < e2
    # ------------------------------------------
    # From section 3a, the cumulative weights are:
    #   w(j) = f * t_j / 4     for j = 2,3,4
    #   w(1) = f - w(2) - w(3) - w(4)
    # where f = t2*t3*t4, t_j = (E-e1)/(e_j-e1).
    #
    # Applying the product rule:
    #   d(f*t_j)/dE = df/dE * t_j + f * dt_j/dE
    #               = gTotal * t_j + f/(e_j - e1)
    # where gTotal = df/dE = 3*(E-e1)^2 / denom
    #   is the total DOS for this case.
    #
    if E < e2:
        denom = (e2-e1) * (e3-e1) * (e4-e1)
        if abs(denom) < tol:
            return [0, 0, 0, 0]
        t2 = (E - e1) / (e2 - e1)
        t3 = (E - e1) / (e3 - e1)
        t4 = (E - e1) / (e4 - e1)
        f = t2 * t3 * t4
        gTotal = 3.0 * (E - e1)**2 / denom

        g(2) = (gTotal*t2 + f/(e2-e1)) / 4
        g(3) = (gTotal*t3 + f/(e3-e1)) / 4
        g(4) = (gTotal*t4 + f/(e4-e1)) / 4
        g(1) = gTotal - g(2) - g(3) - g(4)
        return g

    # ------------------------------------------
    # Case 2: e2 <= E < e3 (middle range)
    # ------------------------------------------
    # From section 3a, the cumulative weights use
    # three sub-tetrahedra volumes v_I, v_II,
    # v_III with intersection parameters a, b,
    # c, d. Their derivatives are:
    #
    #   da/dE = 1/e31,  db/dE = 1/e41
    #   dc/dE = 1/e32,  dd/dE = 1/e42
    #
    #   dv_I/dE   = da*b + a*db
    #             = b/e31 + a/e41
    #   dv_II/dE  = da*d*(1-b) + a*dd*(1-b)
    #             + a*d*(-db)
    #             = d*(1-b)/e31 + a*(1-b)/e42
    #             - a*d/e41
    #   dv_III/dE = (-da)*c*d + (1-a)*dc*d
    #            + (1-a)*c*dd
    #             = -c*d/e31 + (1-a)*d/e32
    #             + (1-a)*c/e42
    #
    # Each cumulative weight w(j) is a linear
    # combination of v_I, v_II, v_III with
    # coefficients that also depend on a, b, c, d.
    # Applying the product rule to each term
    # and collecting:
    #
    if E < e3:
        e31 = e3-e1;  e41 = e4-e1
        e32 = e3-e2;  e42 = e4-e2
        if e31*e41 < tol or e32*e42 < tol:
            return [0, 0, 0, 0]

        a = (E-e1) / e31
        b = (E-e1) / e41
        c_var = (E-e2) / e32
        d_var = (E-e2) / e42

        # Sub-tetrahedra volumes.
        v_I   = a * b
        v_II  = a * d_var * (1 - b)
        v_III = (1 - a) * c_var * d_var

        # Volume derivatives.
        dv_I   = b/e31 + a/e41
        dv_II  = d_var*(1-b)/e31
                 + a*(1-b)/e42
                 - a*d_var/e41
        dv_III = -c_var*d_var/e31
                 + (1-a)*d_var/e32
                 + (1-a)*c_var/e42

        # Parameter derivatives.
        da = 1/e31;  db = 1/e41
        dc = 1/e32;  dd = 1/e42

        # Corner 1: w(1) = [v_I*(3-a-b)
        #   + v_II*(2-a-b) + v_III*(1-a)] / 4
        g(1) = (dv_I*(3-a-b)
                + v_I*(-da - db)
                + dv_II*(2-a-b)
                + v_II*(-da - db)
                + dv_III*(1-a)
                + v_III*(-da)) / 4

        # Corner 2: w(2) = [v_I + v_II*(2-d)
        #   + v_III*(3-c-d)] / 4
        g(2) = (dv_I
                + dv_II*(2-d_var)
                + v_II*(-dd)
                + dv_III*(3-c_var-d_var)
                + v_III*(-dc - dd)) / 4

        # Corner 3: w(3) = [v_I*a + v_II*a
        #   + v_III*(a+c)] / 4
        g(3) = (dv_I*a + v_I*da
                + dv_II*a + v_II*da
                + dv_III*(a+c_var)
                + v_III*(da + dc)) / 4

        # Corner 4: w(4) = [v_I*b
        #   + v_II*(b+d) + v_III*d] / 4
        g(4) = (dv_I*b + v_I*db
                + dv_II*(b+d_var)
                + v_II*(db + dd)
                + dv_III*d_var
                + v_III*dd) / 4

        return g

    # ------------------------------------------
    # Case 3: e3 <= E < e4
    # ------------------------------------------
    # From section 3a, the cumulative weights use
    # the unoccupied sub-tet fraction f_un and
    # parameters s_j = (e4-E)/(e4-e_j). Their
    # derivatives are:
    #   ds_j/dE = -1/(e4-e_j)
    #   df_un/dE = -gTotal  (where gTotal is the
    #     total DOS for this case)
    #
    # Applying the product rule to
    #   w(j) = 1/4 - f_un*s_j/4 for j=1,2,3:
    #   dw(j)/dE = -(df_un*s_j + f_un*ds_j)/4
    #            = (gTotal*s_j + f_un/(e4-e_j))/4
    #
    denom = (e4-e1) * (e4-e2) * (e4-e3)
    if abs(denom) < tol:
        return [0, 0, 0, 0]
    s1 = (e4 - E) / (e4 - e1)
    s2 = (e4 - E) / (e4 - e2)
    s3 = (e4 - E) / (e4 - e3)
    f_un = s1 * s2 * s3
    gTotal = 3.0 * (e4 - E)**2 / denom

    g(1) = (gTotal*s1 + f_un/(e4-e1)) / 4
    g(2) = (gTotal*s2 + f_un/(e4-e2)) / 4
    g(3) = (gTotal*s3 + f_un/(e4-e3)) / 4
    g(4) = gTotal - g(1) - g(2) - g(3)
    return g
```

---

## 3. electronPopulation_LAT (DESIGN 1.5)

```
function computeElectronPopulation_LAT(
        eigenValues, tetrahedra,
        numTetrahedra, tetraVol,
        eFermi, numStates,
        numKPoints, numSpins):
    # Computes electronPopulation_LAT(n, k, spin):
    # the LAT analog of electronPopulation for
    # integrated properties (effective charge,
    # bond order). Each entry gives the fractional
    # electron occupation of state (n, k) as
    # determined by tetrahedron integration.
    allocate electronPopulation_LAT(
        numStates, numKPoints, numSpins) = 0.0

    for spin = 1 to numSpins:
        for n = 1 to numStates:
            for T = 1 to numTetrahedra:
                corners(1:4) = tetrahedra(1:4, T)
                eps_raw(1:4) =
                    eigenValues(n, corners, spin)

                # Sort and track permutation
                sigma = argsort(eps_raw)
                eps(1:4) = eps_raw(sigma)

                # Corner integration weights for
                # occupied-state integration
                # (Bloechl eqs. 18-21, evaluated
                # at the Fermi energy)
                cornerIntgWt_LAT(1:4) =
                    bloechlCornerWeights(
                        eFermi, eps)

                # Bloechl's curvature correction
                # (DESIGN 1.3.1), added rather than
                # folded in, so the weights above stay
                # the paper's uncorrected expressions.
                # It sums to zero across the four
                # corners, so this loop's contribution
                # to the electron COUNT is unchanged
                # and only the distribution over
                # k-points moves.
                cornerCorrWt(1:4) =
                    bloechlCornerCorrection(
                        eFermi, eps)

                for i = 1 to 4:
                    ki = corners(sigma(i))
                    electronPopulation_LAT(
                        n, ki, spin) +=
                        (cornerIntgWt_LAT(i)
                         + cornerCorrWt(i))
                            * tetraVol

    return electronPopulation_LAT
```

This one routine serves BOTH consumers of the corrected
weights: the bond/effective-charge path calls it after the
SCF (DESIGN 1.5), and `populateLAT` calls it at each
converged trial Fermi level inside the SCF (DESIGN 1.6).
Adding the correction here therefore reaches both, and there
is no second insertion point to keep in step.  The DOS path
is NOT a consumer -- it uses `bloechlCornerDOSWt`, and eq. 22
corrects `w` rather than `dw/dE` (DESIGN 1.3.1).

---

## 3a. Bloechl Corner Integration Weights (DESIGN 1.5)

### Motivation

Section 2 (LAT TDOS) computes a single number g(E)
for each tetrahedron: the total DOS contribution at
energy E. Section 3 (electronPopulation_LAT) calls
`bloechlCornerWeights(E, eps)` to decompose the
tetrahedron's occupation into four separate corner
weights. This section derives those weights from first
principles and presents the pseudocode.

**Why corners need separate weights.** Each corner of
a tetrahedron corresponds to a different k-point with
different eigenvector projections (Mulliken populations,
orbital character). The total DOS does not need this
decomposition because it depends only on eigenvalues.
Partial properties (effective charge, bond order, PDOS)
require knowing how much of each corner's projection
to include. The corner weights provide exactly this
decomposition.

### Definitions

Consider a tetrahedron with four corners having sorted
eigenvalues e1 <= e2 <= e3 <= e4. Within the
tetrahedron, the eigenvalue is linearly interpolated
via barycentric coordinates:

  epsilon(r) = lambda_1 * e1 + lambda_2 * e2
             + lambda_3 * e3 + lambda_4 * e4

where lambda_i >= 0 and sum(lambda_i) = 1.

The **corner integration weight** w_j(E) is defined as
the integral of the j-th barycentric coordinate over
the occupied region {epsilon <= E}, normalized by the
tetrahedron volume V_T:

  w_j(E) = (1/V_T) * integral_{epsilon<=E} lambda_j dV

Since sum(lambda_j) = 1, the four weights sum to the
total occupied fraction:

  f(E) = sum_j w_j(E) = (1/V_T) * Vol({epsilon <= E})

### Key property: vertex averaging

The integral of any linear function L(r) over a
tetrahedron equals the volume times the average of L
at the four vertices:

  integral_T L dV = V_T * [L(v1)+L(v2)+L(v3)+L(v4)]/4

Since barycentric coordinates are linear functions,
this means: if we decompose the occupied region into
sub-tetrahedra, each corner weight contribution from
a sub-tetrahedron S is:

  w_j(S) = (V_S / V_T)
         * [sum of lambda_j at S's 4 vertices] / 4

### Case 0: trivial bounds

  E < e1:   w_j = 0 for all j  (empty region)
  E >= e4:  w_j = 1/4 for all j  (full tetrahedron)

### Case 1: e1 <= E < e2

Only corner 1 has eigenvalue below E. The occupied
region is a small tetrahedron with apex at corner 1,
cut by the iso-energy surface epsilon = E. The surface
intersects the three edges from corner 1 at:

  edge 1->j at parameter t_j = (E-e1)/(ej-e1)
                                for j = 2, 3, 4

The sub-tetrahedron has four vertices with barycentric
coordinates (lambda_1, lambda_2, lambda_3, lambda_4):

  corner 1:          (1,     0,   0,   0)
  edge 1->2 at t_2:  (1-t_2, t_2, 0,   0)
  edge 1->3 at t_3:  (1-t_3, 0,   t_3, 0)
  edge 1->4 at t_4:  (1-t_4, 0,   0,   t_4)

Volume ratio: f = t_2 * t_3 * t_4

Applying vertex averaging (summing lambda_j across
the 4 vertices, multiplying by f/4):

  w_2 = f * t_2 / 4
  w_3 = f * t_3 / 4
  w_4 = f * t_4 / 4
  w_1 = f - w_2 - w_3 - w_4

Verification: sum(w_j) = f.

### Case 4: e3 <= E < e4 (complement of Case 1)

The *unoccupied* region is a small tetrahedron near
corner 4. Define:

  s_j = (e4 - E) / (e4 - ej)    for j = 1, 2, 3

Unoccupied sub-tetrahedron vertices:

  corner 4:          (0,   0,   0,   1)
  edge 4->1 at s_1:  (s_1, 0,   0,   1-s_1)
  edge 4->2 at s_2:  (0,   s_2, 0,   1-s_2)
  edge 4->3 at s_3:  (0,   0,   s_3, 1-s_3)

Unoccupied fraction: f_unocc = s_1 * s_2 * s_3
Occupied fraction:   f = 1 - f_unocc

The occupied weights are the whole-tetrahedron weights
(1/4 each) minus the unoccupied contributions:

  w_1 = 1/4 - f_unocc * s_1 / 4
  w_2 = 1/4 - f_unocc * s_2 / 4
  w_3 = 1/4 - f_unocc * s_3 / 4
  w_4 = f - w_1 - w_2 - w_3

Verification: sum(w_j) = f.

### Case 2: e2 <= E < e3 (middle range)

Corners 1 and 2 lie below E; corners 3 and 4 lie
above. The iso-energy surface cuts four edges:

  edge 1->3 at  a = (E-e1)/(e3-e1)     point A
  edge 1->4 at  b = (E-e1)/(e4-e1)     point B
  edge 2->3 at  c = (E-e2)/(e3-e2)     point C
  edge 2->4 at  d = (E-e2)/(e4-e2)     point D

The occupied region is a pentahedron with vertices
{corner 1, corner 2, A, B, C, D}. We decompose it
into three sub-tetrahedra:

  T_I   = (corner 1, corner 2, A, B)
  T_II  = (corner 2, A, B, D)
  T_III = (corner 2, A, C, D)

The volume ratios follow from the determinant of
the 4x4 barycentric coordinate matrix for each
sub-tetrahedron:

  v_I   = a * b
  v_II  = a * d * (1 - b)
  v_III = (1 - a) * c * d

Occupied fraction: f = v_I + v_II + v_III

The barycentric coordinates at each vertex:

  T_I:
    corner 1   (1,     0,     0,   0)
    corner 2   (0,     1,     0,   0)
    A          (1-a,   0,     a,   0)
    B          (1-b,   0,     0,   b)

  T_II:
    corner 2   (0,     1,     0,   0)
    A          (1-a,   0,     a,   0)
    B          (1-b,   0,     0,   b)
    D          (0,     1-d,   0,   d)

  T_III:
    corner 2   (0,     1,     0,   0)
    A          (1-a,   0,     a,   0)
    C          (0,     1-c,   c,   0)
    D          (0,     1-d,   0,   d)

Summing lambda_j over the four vertices of each
sub-tetrahedron, multiplying by v_k/4, and summing
over sub-tetrahedra gives the corner weights:

  w_1 = [v_I*(3-a-b) + v_II*(2-a-b)
         + v_III*(1-a)] / 4

  w_2 = [v_I + v_II*(2-d)
         + v_III*(3-c-d)] / 4

  w_3 = [v_I*a + v_II*a
         + v_III*(a+c)] / 4

  w_4 = [v_I*b + v_II*(b+d)
         + v_III*d] / 4

Verification: for each sub-tetrahedron, the sum of
all four lambda_j at any vertex is 1, so the sum over
all four vertices is 4. Therefore:
  sum(w_j) = (4*v_I + 4*v_II + 4*v_III) / 4
           = v_I + v_II + v_III = f.

### Continuity between cases

The formulas are continuous at the case boundaries:

- At E = e2: Case 2 reduces to Case 1 because
  c = d = 0, so v_II = v_III = 0 and
  f = v_I = a*b = (e2-e1)^2 / [(e3-e1)(e4-e1)].
  All four corner weights match.

- At E = e3: Case 2 reduces to Case 4 because
  a = 1, so v_III = 0 and the occupied fraction
  equals 1 - (e4-e3)^2 / [(e4-e1)(e4-e2)].
  All four corner weights match.

### Derivative consistency with TDOS

The energy derivative of sum(w_j) must equal the
total per-tetrahedron DOS. This relationship is now
built into `bloechlCornerDOSWt` (section 2a): the
four `cornerDOSWt_LAT` values are defined so that
`sum(cornerDOSWt_LAT) = gTotal`. This was verified
numerically for the middle range: with e1=0, e2=1,
e3=3, e4=5, both the derivative of
f = v_I + v_II + v_III and the TDOS formula give
g(2) = 51/120.

### Pseudocode

```
function bloechlCornerWeights(E, eps):
    # eps = [e1, e2, e3, e4] sorted ascending.
    # Returns w(1:4): the integrated corner
    # weights at energy E for one tetrahedron.
    #
    # w(i) is the fraction of the tetrahedron's
    # occupation attributed to sorted corner i.
    # sum(w) = f(E), the occupied volume fraction.

    e1 = eps(1);  e2 = eps(2)
    e3 = eps(3);  e4 = eps(4)
    tol = 1.0e-12

    # Case 0: trivial bounds
    if E < e1:
        return [0, 0, 0, 0]
    if E >= e4:
        return [0.25, 0.25, 0.25, 0.25]

    # Case 1: e1 <= E < e2
    if E < e2:
        denom = (e2-e1) * (e3-e1) * (e4-e1)
        if abs(denom) < tol:
            return [0, 0, 0, 0]
        t2 = (E - e1) / (e2 - e1)
        t3 = (E - e1) / (e3 - e1)
        t4 = (E - e1) / (e4 - e1)
        f = t2 * t3 * t4
        w(2) = f * t2 / 4
        w(3) = f * t3 / 4
        w(4) = f * t4 / 4
        w(1) = f - w(2) - w(3) - w(4)
        return w

    # Case 2: e2 <= E < e3
    if E < e3:
        e31 = e3-e1;  e41 = e4-e1
        e32 = e3-e2;  e42 = e4-e2
        if e31*e41 < tol or e32*e42 < tol:
            return [0, 0, 0, 0]
        a = (E-e1) / e31
        b = (E-e1) / e41
        c = (E-e2) / e32
        d = (E-e2) / e42

        v_I   = a * b
        v_II  = a * d * (1 - b)
        v_III = (1 - a) * c * d

        w(1) = (v_I*(3-a-b) + v_II*(2-a-b)
                + v_III*(1-a)) / 4
        w(2) = (v_I + v_II*(2-d)
                + v_III*(3-c-d)) / 4
        w(3) = (v_I*a + v_II*a
                + v_III*(a+c)) / 4
        w(4) = (v_I*b + v_II*(b+d)
                + v_III*d) / 4
        return w

    # Case 3: e3 <= E < e4
    denom = (e4-e1) * (e4-e2) * (e4-e3)
    if abs(denom) < tol:
        return [0.25, 0.25, 0.25, 0.25]
    s1 = (e4 - E) / (e4 - e1)
    s2 = (e4 - E) / (e4 - e2)
    s3 = (e4 - E) / (e4 - e3)
    f_un = s1 * s2 * s3
    w(1) = 0.25 - f_un * s1 / 4
    w(2) = 0.25 - f_un * s2 / 4
    w(3) = 0.25 - f_un * s3 / 4
    w(4) = (1 - f_un) - w(1) - w(2) - w(3)
    return w
```

### The curvature correction (DESIGN 1.3.1)

Linear interpolation of eps(k) inside a tetrahedron misplaces
the iso-energy surface.  Bloechl's correction compensates it
by shifting weight between the four corners:

!! dw_i = (1/40) * D_T(E_F) * sum_{j=1..4} (eps_j - eps_i)
!!      = (1/10) * D_T(E_F) * (epsBar - eps_i)

with `D_T(E_F)` this tetrahedron's DOS at the Fermi level and
`epsBar` the mean of its four corner eigenvalues.  The second
form is the one implemented -- it is the same quantity, and it
makes the zero-sum property self-evident rather than something
a reader has to derive before trusting the routine near a
converged SCF.

`D_T(E_F)` needs no new formula: it is the sum of the corner
DOS weights already specified in section 2a, evaluated at the
same energy, before any `tetraVol` factor.  The routine forms
it internally so that it stays a pure function of one
tetrahedron's corners, exactly like its two siblings, and
inherits their degenerate-corner guards by calling through
them.

```
function bloechlCornerCorrection(E, eps):
    # eps(1:4) sorted ascending, as for the two
    # routines above.  Returns the four dw_i.

    # A tetrahedron with no iso-energy surface through
    # it has D_T = 0, so the correction vanishes.  The
    # two bounds are handled explicitly rather than
    # left to arithmetic, because they are the common
    # case (every fully occupied or fully empty
    # tetrahedron in the mesh) and because it is the
    # property that keeps insulators untouched.
    if E < eps(1) or E >= eps(4):
        return [0, 0, 0, 0]

    cornerDOSWt(1:4) = bloechlCornerDOSWt(E, eps)
    D_T = sum(cornerDOSWt(1:4))
    epsBar = sum(eps(1:4)) / 4

    for i = 1 to 4:
        dw(i) = D_T * (epsBar - eps(i)) / 10
    return dw
```

**The self-check this buys.**  `sum(dw) == 0` to rounding, for
ANY energy and any four corner energies, because
`sum_i (epsBar - eps_i)` is identically zero.  That is a
complete test of the routine requiring no reference values and
no other quantity -- which is the reason the correction is a
separate routine rather than folded into
`bloechlCornerWeights`, where it could only be checked against
hand-computed expectations.

**Equation 22 is the whole correction.**  The reference list in
DESIGN once cited "eqs. 22-24", which invited the belief that
two further terms were owed here.  They are not: 23 and 24
compare the true Fermi surface against the interpolated
polyhedral one, which the paper does as an assessment rather
than as a formula anything computes.  This routine is complete
as written.

---

## 3a. LAT in the SCF Occupation Path (DESIGN 1.6)

Section 3 builds `electronPopulation_LAT` at a Fermi level
someone else determined, which suits a property computed
after the SCF has finished.  Inside the SCF the Fermi level
moves every iteration and must be determined by the SAME
integration scheme that supplies the weights, or the two
disagree about how many electrons are present (DESIGN 1.6a).
This section adds that determination, its call site, and the
one substitution it needs in the charge accumulation.

Three facts about the existing Gaussian path are load-bearing
here and are mirrored rather than reinvented:

- ONE Fermi level serves both spin channels, constrained by
  the total `numElectrons`.  The Gaussian path merges every
  (state, spin, kpoint) triplet into one sorted list and
  fills it in ascending energy, so the magnetic moment is an
  OUTCOME of which channel sorts lower, never an input.  LAT
  keeps that; two per-channel levels would be a physics
  change wearing an integration-scheme costume.
- `populateStandard` runs unconditionally and produces
  `occupiedEnergyIndex`, which is what brackets the
  smearing search.  The LAT search seeds its bracket the same
  way, so the initial guess has one source.
- The excited-state (XANES) correction is NOT automatic.
  `populateSmearing` raises `numElectrons` by one, populates,
  then calls `correctCorePopulation`, which removes the core
  electron's occupancy AND restores the count.  A LAT search
  carries the same obligation or it counts the core hole's
  electron.

  **v1 REFUSES this combination rather than guessing it.**  The
  Gaussian correction addresses the flat, sorted occupation
  array through `indexEnergyEigenValues`, and its band
  arithmetic does not carry over to the `(band, kpoint, spin)`
  array without a derivation this work has not done -- the
  orbital-state counts are scaled by `spin` at initialization,
  so the mapping is not the obvious one.  A wrong correction
  misplaces exactly one electron, which is small enough to
  read as a convergence problem and never be questioned.  So
  `populateLAT` stops with a message naming Gaussian
  integration as the supported path for an excited run.  The
  ground-state path is unaffected, and that is the path the
  metals work needs.

```
function populateStates:
    # DESIGN 1.6. Unchanged through populateStandard, which
    # every path still needs for its bracket and for the
    # degeneracy-averaged starting occupations.
    if excitedQN_n /= 0 and coreStructInit == 0:
        call initCoreStateStructures
    call populateStandard

    # LAT and thermal smearing are ALTERNATIVES, not layers
    # (DESIGN 1.6e): tetrahedron integration determines
    # occupations geometrically and needs no broadening, so a
    # run that set both would be asking for two answers.  LAT
    # wins and thermalSigma is ignored for the SCF occupation;
    # say so rather than applying one and discarding the other.
    if kPointIntgCode == 1:
        call populateLAT
    else if thermalSigma /= 0:
        call populateSmearing

    write occupiedEnergy
```

```
function populateLAT:
    # Find the Fermi level from the TETRAHEDRON integral and
    # fill electronPopulation_LAT at it (DESIGN 1.6a).
    #
    # Bracket exactly as populateSmearing does: populateStandard
    # has already put occupiedEnergy between the highest
    # occupied and lowest unoccupied sorted eigenvalues, and
    # fermiSearchLimit widens that to a range no real Fermi
    # level escapes.
    minEnergy = sortedEnergyEigenValues(occupiedEnergyIndex)
                - fermiSearchLimit
    maxEnergy = sortedEnergyEigenValues(occupiedEnergyIndex+1)
                + fermiSearchLimit

    # Seed from the previous SCF iteration when there is one:
    # the level moves little between iterations, so Newton
    # usually converges in two or three passes.
    E = (previousFermi if previousFermi is set
         else 0.5 * (minEnergy + maxEnergy))

    # An excited run does not reach here (refused above).  When it
    # is implemented, the target is a LOCAL variable rather than a
    # mutation of the module's numElectrons: the Gaussian form
    # increments and restores across two routines, which is correct
    # but reads as unbalanced inside either one, and a search loop
    # is the worst place to leave module state temporarily wrong.
    targetCount = numElectrons

    for i = 1 to MAX_FERMI_ITER:
        (N, dNdE) = latElectronCount(E)

        if abs(N - targetCount) < smallThresh:
            break

        # Maintain the bracket from every evaluation, so the
        # safeguard below always has a valid one.  N(E) is
        # monotone non-decreasing, which is what makes this
        # sound.
        if N < targetCount: minEnergy = E
        else:               maxEnergy = E

        # Newton where the derivative is usable, bisection
        # where it is not.  dN/dE is the density of states, so
        # it VANISHES inside a gap -- an insulator's first
        # guess lands there and a bare Newton step would
        # diverge.  Fall back whenever the step leaves the
        # bracket or the slope is negligible (DESIGN 1.6a).
        if dNdE > slopeThresh:
            E_next = E - (N - targetCount) / dNdE
        else:
            E_next = 0.5 * (minEnergy + maxEnergy)
        if E_next <= minEnergy or E_next >= maxEnergy:
            E_next = 0.5 * (minEnergy + maxEnergy)
        E = E_next

    occupiedEnergy = E
    previousFermi  = E
    call computeElectronPopulation_LAT(..., eFermi = E)  # 3
```

```
function latElectronCount(E):
    # The integrated electron count at trial energy E and its
    # derivative, in ONE pass (DESIGN 1.6a).  The derivative is
    # free: cornerDOSWt_LAT is by construction d/dE of
    # cornerIntgWt_LAT (section 2), so the same corner sort
    # yields both.
    #
    # Cost note (DESIGN 1.6b): the loop is over the FULL mesh's
    # tetrahedra, not the IBZ.  Most (T, n) pairs are wholly
    # below or wholly above E and contribute a constant with no
    # corner work, so only the straddling pairs are re-evaluated
    # as E moves.
    N = 0 ; dNdE = 0
    for spin = 1 to numSpins:
        for n = 1 to numStates:
            for T = 1 to numTetrahedra:
                eps(1:4) = sorted eigenValues(n, corners(T), spin)
                if eps(4) <= E:            # wholly occupied
                    N = N + tetraVol * spinFactor
                    continue
                if eps(1) > E:  continue   # wholly empty
                w  = cornerIntgWt_LAT(eps, E)   # section 3
                dw = cornerDOSWt_LAT(eps, E)    # section 2
                N    = N    + tetraVol * spinFactor * sum(w)
                dNdE = dNdE + tetraVol * spinFactor * sum(dw)
    return (N, dNdE)
```

**The normalization must be checked, not assumed (DESIGN
1.6d).**  `electronPopulation` carries the k-point weight
folded in and divides by `spin`, so each state holds one
electron when spin-polarized and two when not.  The LAT sum
carries its own Brillouin-zone fraction instead: a fully
occupied band gives `sum over T of tetraVol = 1` per channel,
so `spinFactor = 2 / numSpins` reproduces the same
convention.  That equality is stated here as the thing to
verify, not as a derivation to trust: the calibration test is
an INSULATOR, where `latElectronCount(E)` for any E inside
the gap must return exactly `numElectrons`.  A wrong factor
is a constant error on the valence charge that an SCF partly
absorbs, so it will not announce itself.

```
function correctCorePopulation_LAT:
    # NOT IMPLEMENTED IN v1 -- populateLAT refuses an excited run
    # instead (above).  Kept as the specification for when it is
    # written, because the shape below is the part that is already
    # settled; what is NOT settled is the band range, since
    # numOrbitalStates is scaled by `spin` at initialization and the
    # Gaussian routine's arithmetic over the flat sorted array does
    # not transfer to this array's band index unexamined.
    #
    # The LAT-shaped core-hole correction (DESIGN 1.6).  Same
    # physics as correctCorePopulation, but simpler: that
    # routine walks a flat array through indexEnergyEigenValues
    # because its occupations are stored sorted, while
    # electronPopulation_LAT is already (n, k, spin) and the
    # excited band index addresses it directly.
    #
    # It also restores numElectrons, exactly as the Gaussian
    # form does -- the increment in the search above and this
    # decrement are ONE pair split across two routines.
    for each core band n in the excited (QN_n, QN_l) orbital:
        for k, spin:
            electronPopulation_LAT(n, k, spin) -=
                electronPopulation_LAT(n, k, spin)
                * numSpins / 2 / numOrbitalStates(...)
    numElectrons = numElectrons - 1
```

```
function valeCharge (LAT branch only):
    # DESIGN 1.6c.  The substitution point is the UNPACK, not
    # the accumulation below it.  The Gaussian path reads a
    # flat electronPopulation walked in (kpoint, spin, state)
    # order -- deliberately NOT the eigenvalue sort order --
    # and unpacks it into structuredElectronPopulation.
    # electronPopulation_LAT is already in that target shape,
    # so the LAT path REPLACES the unpack rather than adding a
    # reordering step.
    #
    # Routing the LAT array through the flat-index loop is the
    # one mistake available here: it would scramble the
    # occupations silently, since both arrays hold plausible
    # numbers of the right size.
    # The copy also CONVERTS CONVENTION, which an earlier draft of
    # this section omitted.  electronPopulation carries the
    # kPointWeight convention, whose weights sum to 2.0 so a
    # non-polarized run holds two electrons per state;
    # electronPopulation_LAT holds pure BZ volume fractions summing
    # to 1.0 per occupied band per spin.  computeBond already applies
    # exactly this factor at its own point of use, so the array keeps
    # ONE convention and every consumer converts (DESIGN 1.6d).
    if kPointIntgCode == 1:
        structuredElectronPopulation =
            electronPopulation_LAT * 2 / numSpins
    else:
        energyLevelCounter = 0
        for i = 1 to numKPoints:
          for j = 1 to numSpins:
            for k = 1 to numStates:
              energyLevelCounter += 1
              structuredElectronPopulation(k,i,j) =
                  electronPopulation(energyLevelCounter)
    # ... accumulation into potRho unchanged ...
```

**Lifecycle.**  `computeElectronPopulation_LAT` runs once per
SCF iteration, not once per run, because the weights depend on
a Fermi level that moves.  Its array has the same lifecycle as
`electronPopulation`: allocated on first use, freed by
`cleanUpPopulation`.  No permutation table is involved --
`potRho` is indexed by potential TYPE, so every component is
already an orbit sum and invariant under the IBZ reduction
(DESIGN 1.6), and LAT changes only a scalar occupation.

---

## 4. Build Atom Permutation Table (DESIGN 2.4, 2.7)

The atom permutation table records, for each point group
operation R and each atom A, which atom B = R(A) the
operation maps A to.  This is the single piece of
infrastructure needed for correct IBZ unfolding of all
shell-summed quantities (Q*, bond order, PDOS modes 0-2).

The algorithm works in fractional (abc) coordinates of
the loaded real lattice (whichever cell ended up in
O_Lattice -- full conventional or primitive reduction).
The rotation matrices and per-operation translations
used here are abcRealPointOps and abcRealFracTrans --
the loaded-cell-abc forms produced by computeRealPointOps
(section 4b below) from the conv-abc operations that
arrive on disk.  Atom Cartesian positions are converted
to fractional using invRealVectors (= recipVectors /
2*pi) of the same loaded lattice, dotting its COLUMNS
(reciprocal vectors), not its rows.

```
function buildAtomPerm(numPointOps, abcRealPointOps,
                       abcRealFracTrans, numAtomSites,
                       atomSites, invRealVectors):
    # Returns atomPerm(numPointOps, numAtomSites)
    #   where atomPerm(R, A) = B means operation R
    #   maps atom A to atom B.  Both R and atom
    #   positions are in the loaded real lattice abc
    #   basis after computeRealPointOps has run.

    allocate atomPerm(numPointOps, numAtomSites)

    # Convert all atom positions from Cartesian (xyz)
    # to fractional (abc) coordinates of the loaded
    # real lattice.
    allocate abcAtomPos(3, numAtomSites)
    # invRealVectors holds reciprocal vectors as COLUMNS,
    #   so the fractional component along axis i is the
    #   dot with column i, not row i.
    for A = 1 to numAtomSites:
        for i = 1 to 3:
            abcAtomPos(i, A) =
                sum(invRealVectors(:,i)
                    * atomSites(A)%cartPos(:))

    # For each operation and atom, apply {R|t} in the
    # loaded-cell abc basis and find the matching atom.
    for R = 1 to numPointOps:
        for A = 1 to numAtomSites:

            # Apply the rotation + translation in the
            # loaded-cell abc basis.
            for i = 1 to 3:
                rotPos(i) =
                    sum(abcRealPointOps(i,:,R)
                        * abcAtomPos(:,A))
                    + abcRealFracTrans(i, R)

            # Wrap the rotated position into [0,1).
            for i = 1 to 3:
                rotPos(i) = modulo(rotPos(i), 1.0)

            # Search for the atom at the rotated
            # position, restricted to atoms of the same
            # type. This restriction is what MAKES the
            # type-level sums invariant rather than a
            # consequence of them: a type need not be a
            # symmetry orbit (see DESIGN 2.3), so the
            # requirement that R stay within a type is
            # imposed here and the stop below enforces
            # it.
            atomPerm(R, A) = -1  # sentinel
            for B = 1 to numAtomSites:
                if atomType(B) != atomType(A):
                    cycle

                # Compute the difference, wrapped
                # into [-0.5, 0.5) on each axis.
                for i = 1 to 3:
                    diff(i) = rotPos(i)
                             - abcAtomPos(i, B)
                    diff(i) = diff(i)
                             - nint(diff(i))

                if all(|diff(:)| < threshold):
                    atomPerm(R, A) = B
                    exit  # found the match

            # Safety check: every atom must have a
            # match. If not, the point group or the
            # atom positions are inconsistent.
            if atomPerm(R, A) == -1:
                error("No match for atom", A,
                      "under operation", R)

    deallocate abcAtomPos
    return atomPerm
```

---

## 4a. Build Inverse Atom Permutation (DESIGN 1.4, 2.4)

The inverse atom permutation invAtomPerm(R, B) gives
the atom A such that atomPerm(R, A) = B, i.e.,
A = R^{-1}(B). It is used during LAT PDOS tetrahedron
corner assembly to map channel indices from full-mesh
k-points back to their IBZ representatives (see
section 8). Built in O_AtomicSites alongside atomPerm.

```
function buildInvAtomPerm(numPointOps,
                          numAtomSites,
                          atomPerm):
    allocate invAtomPerm(numPointOps,
                         numAtomSites)

    for R = 1 to numPointOps:
        for A = 1 to numAtomSites:
            B = atomPerm(R, A)
            invAtomPerm(R, B) = A

    return invAtomPerm
```

---

## 4b. Conv-abc On-Disk Operations and Lattice
       Conjugation (DESIGN 2.7)

Symmetry operations cross two boundaries on the way from
the space-group database to `buildAtomPerm`: a producer-
side write step in `makeinput.py` that emits each
operation in its native conventional-cell-abc fractional
form (the spaceDB convention) into the kp file, and a
consumer-side step in imago (Fortran) that conjugates
those conv-abc operations into the basis of the lattice
currently loaded in O_Lattice (full conventional cell or
primitive reduction, depending on the skeleton's `full`
/ `prim` flag).  The on-disk format also carries two
small metadata blocks -- `CONV_LATTICE` (the
conventional-cell matrix in Bohr) and `CELL_MODE`
(`full` or `prim`) -- which give the consumer the
inputs it needs to form the change-of-basis matrix and
to choose between the full conjugation path and a
`full`-mode identity shortcut.  Both boundaries live
behind `realVectors`-type lattice matrices and require
no special-casing by cell type, centering, or full-vs-
prim mode beyond the identity shortcut.  See DESIGN 2.7
for the motivation and full background.

### 4b.1 Writer Additions (makeinput.py kp-file writer)

The previous design used a producer-side similarity
helper (`_to_cartesian_ops`) that converted spaceDB
operations into a Cartesian xyz intermediate before
writing.  Under the new design that helper is removed:
each spaceDB operation is written into the kp file
exactly as it appears in `share/spaceDB/<sg>` -- three
rotation lines plus one translation line per operation
-- with no producer-side math applied to the matrix
entries or fractional translations.

Two small additions accompany the existing block of
operation lines: `CONV_LATTICE` and `CELL_MODE`.

```
function writeKPointSymmetryBlock(point_ops, frac_trans,
                                  conv_lattice,
                                  cell_mode):
    # point_ops, frac_trans:  conv-abc fractional
    #     entries lifted straight from
    #     share/spaceDB/<sg>
    # conv_lattice:           sc.full_cell_real_lattice,
    #     the conventional-cell snapshot captured at the
    #     top of apply_space_group() before any
    #     primitive reduction may overwrite the
    #     in-memory lattice
    # cell_mode:              'full' or 'prim' from the
    #     skeleton's lattice-mode flag

    write 'POINT_OPS', numPointOps
    for each (R_conv_abc, t_conv_abc) in
            (point_ops, frac_trans):
        write 3 rows of R_conv_abc (one per line)
        write the 3-component t_conv_abc on a 4th line

    # New: emit the conventional lattice in Bohr so the
    # consumer can form M_loaded^{-1} * M_conv without
    # carrying implicit cell-choice knowledge.
    write 'CONV_LATTICE'
    write 3 rows of conv_lattice (Bohr)

    # New: emit the cell-mode flag so the consumer can
    # take the identity shortcut when the loaded cell
    # equals the conventional cell.
    write 'CELL_MODE'
    write cell_mode    # 'full' or 'prim'
```

The `POINT_OPS` block is byte-identical to the spaceDB
entries for every cell type (cubic, hex, monoclinic,
triclinic, ...); only `CONV_LATTICE` and `CELL_MODE`
are new on-disk content.

### 4b.2 Reader Additions (imago readKPoints)

For style codes 1 and 2, the imago Fortran reader parses
the existing operations block into the renamed arrays
`convAbcPointOps(3, 3, numPointOps)` and
`convAbcFracTrans(3, numPointOps)` -- same on-disk
layout as before, only the destination array names
change to reflect the basis the entries live in.  Two
new parse steps follow the operations block:

```
function readKPointSymmetryBlock(file):
    # Existing: POINT_OPS plus per-operation
    # translation.
    read 'POINT_OPS', numPointOps
    for i = 1 to numPointOps:
        read 3 rows into convAbcPointOps(:,:,i)
        read 3-vector into convAbcFracTrans(:,i)

    # New: conventional-cell matrix (Bohr), 3 rows.
    read 'CONV_LATTICE'
    read 3 rows into convLattice(:,:)

    # New: cell-mode flag, single string token.
    read 'CELL_MODE'
    read string into cellMode    # 'full' or 'prim'

    return (numPointOps, convAbcPointOps,
            convAbcFracTrans, convLattice, cellMode)
```

Style code 0 still synthesizes identity-only operations
in memory and does not require the new blocks; `cellMode`
defaults to `full` and `convLattice` defaults to
`realVectors` for that path so the consumer-side
shortcut applies trivially.

### 4b.3 Consumer-Side Lattice Conjugation
       (imago Fortran)

Run inside `initializeKPoints` once per kp-file load,
right after `readKPoints` deposits the on-disk
operations into `convAbcPointOps` / `convAbcFracTrans`
and the metadata into `convLattice` / `cellMode`.  Two
siblings:

```
function computeRealPointOps(numPointOps,
                             convAbcPointOps,
                             convAbcFracTrans,
                             cellMode, convLattice,
                             realVectors):
    # Conjugate conv-abc operations into the basis of
    # the loaded real lattice for use by buildAtomPerm.
    # Let L (= realVectors) and Lc be the loaded and
    # conventional lattices with vectors as COLUMNS.
    # convLattice stores Lc's vectors as ROWS, so
    # Lc = transpose(convLattice).  The change of basis
    # carrying loaded fractional to conv fractional is
    #   T          = Lc^{-1} * L         (r_conv = T r_loaded)
    # Direct-space coordinates are covariant, so:
    #   R_loaded   = T^{-1} * R_conv * T
    #   t_loaded   = T^{-1} * t_conv
    allocate abcRealPointOps(3, 3, numPointOps)
    allocate abcRealFracTrans(3, numPointOps)

    if cellMode == 'full':
        # Identity shortcut: L == Lc, so T = I and the
        # conjugation collapses to a copy.  Also serves
        # style 0, whose convLattice is a column-layout
        # copy of realVectors the T path must not touch.
        for i = 1 to numPointOps:
            abcRealPointOps(:,:,i) =
                convAbcPointOps(:,:,i)
            abcRealFracTrans(:,i) =
                convAbcFracTrans(:,i)
    else:
        # Prim path: form T once, then conjugate per op.
        T     = inverse_3x3(transpose(convLattice)) * realVectors
        T_inv = inverse_3x3(T)
        for i = 1 to numPointOps:
            abcRealPointOps(:,:,i) =
                T_inv * convAbcPointOps(:,:,i) * T
            abcRealFracTrans(:,i) =
                T_inv * convAbcFracTrans(:,i)

    return (abcRealPointOps, abcRealFracTrans)
```

```
function computeRecipPointOps(numPointOps,
                              abcRealPointOps):
    # Build the reciprocal-space operations for IBZ
    # folding.  k-point (reciprocal) coordinates are
    # contravariant -- the dual of the covariant direct
    # coordinates -- so the reciprocal representation of
    # an operation is the INVERSE TRANSPOSE of its
    # direct-space representation:
    #   R_recip = (R_real)^{-T}
    # This consumes the abcRealPointOps already built by
    # computeRealPointOps (which therefore must run
    # first) and needs no lattice matrices of its own.
    allocate abcRecipPointOps(3, 3, numPointOps)

    for i = 1 to numPointOps:
        abcRecipPointOps(:,:,i) =
            transpose(inverse_3x3(abcRealPointOps(:,:,i)))

    return abcRecipPointOps
```

`computeRealPointOps` runs in every style-code branch
(style 0 sets up a trivial identity op; styles 1 and 2
read real symmetry from the kp file).
`computeRecipPointOps` runs only for styles 1 and 2
(IBZ folding) and must be called AFTER
`computeRealPointOps`.  The `cellMode` flag selects the
identity shortcut versus the full conjugation path in
the real-space routine -- no other `full`-vs-`prim`
branching exists outside these two routines.

---

## 4c. K-Point Mesh Selection and Reduction
        (DESIGN 3.7-3.11)

This section turns a requested k-point volume density into
a resolved Monkhorst-Pack mesh, then either keeps that mesh
whole or folds it to the irreducible zone.  The mesh is the
single primary product; the two downstream uses are
selected by one flag:

- **Full mesh** (`reduce = false`): the pure Monkhorst-Pack
  grid, kept intact for Linear Analytic Tetrahedral
  integration (sections 1-3).  Tetrahedra tile the full
  grid, so no point may be folded away.
- **Reduced mesh** (`reduce = true`): the same grid folded
  to the IBZ by the point group (4c.5), for the symmetry-
  reduced eigenvalue problem (sections 5-7).

In density mode `makeinput.py` writes the symmetry block and
density into the kp file (section 4b), and the map below runs
inside `imago` (`kpoints.f90`) once the reciprocal lattice is
known.  The `applySymmetry` argument already carried by
`initializeKPointMesh` is exactly the `reduce` flag.

**Placement note (for the follow-on shift task).** The
automatic shift (4c.3) depends on the parity of the
resolved counts, and in density mode the counts are
resolved here, inside `imago`, not in `makeinput.py`.  The
automatic shift must therefore be chosen here as well, not
assigned by crystal-system name upstream.  The removal of
the upstream by-name assignment across the makeinput ->
kp-file -> `imago` chain is specified separately; this
section specifies only the algorithm that replaces it.

### 4b.4 Reading the on-disk operations (shared, symmetry.py)

Both the kp-file writer (4b.1) and the producer's axis-class
sourcing (4c.7) need a space group's point operations off
disk, in the conventional-abc fractional form
`share/spaceDB/<sg>` stores.  ONE reader serves both, in the
focused `symmetry.py` module (ARCHITECTURE 2, 7), so a single
parser of the spaceDB operation file cannot drift and the
classes the producer seeds from are derived from the very
operations the kp writer emits.

The spaceDB file layout is a description line, a root
space-group number / sub-number line, a `numSpaceOps
numShifts` line, then per operation a blank line, three rows
of a 3x3 rotation matrix, and one fractional-translation row.
The first `numSpaceOps / numShifts` operations are the pure
point group (no centering translation); those are the ones
read.

```
function read_conv_abc_point_ops(space_group_db,
                                 space_group_name):
    # space_group_db is the spaceDB directory and
    #   space_group_name the file in it -- both carried by a
    #   loaded StructureControl (space_group_db /
    #   space_group_name) and by makeinput's settings (space_db
    #   / space_group_name), which adapt the two names.
    open join(space_group_db, space_group_name)
    skip the description line
    skip the root space-group number line
    (num_space_ops, num_shifts) = read two ints
    num_point_ops = num_space_ops / num_shifts    # integer div

    point_ops  = []           # each a 3x3 rotation, rows
    frac_trans = []           # each a length-3 translation
    for i = 1 to num_point_ops:
        skip the blank line before the operation
        rows  = read 3 rows of 3 floats each
        trans = read 1 row of 3 floats
        append rows  to point_ops
        append trans to frac_trans
    return (point_ops, frac_trans)
```

The kp writer consumes both returns -- it emits the rotations
and carries the fractional translations for `buildAtomPerm`
(4b.1); the producer's `build_climb_config` (11.4) consumes
only the rotations, feeding them to `axisClassesForCell`
(4c.7).  `StructureControl` exposes a thin `point_ops`
accessor that delegates here and returns just the rotations,
for the producer's "ask the loaded cell" path; `makeinput`'s
former `_extract_point_ops` becomes a thin call to this
reader, so the parse lives in exactly one place.

### 4c.1 Axis classes (DESIGN 3.8)

Two reciprocal axes are *coupled* when some operation
connects them off the diagonal; the transitive closure of
coupling gives the axis classes that must share a count.
The reciprocal-space operations are used, since the mesh
transforms under them.

```
function computeAxisClasses(recipPointOps, numPointOps):
    # Union-find over the three axes.
    parent = [1, 2, 3]

    for m = 1 to numPointOps:
        R = recipPointOps(:,:,m)
        for i = 1 to 3:
            for j = 1 to 3:
                if i != j and R(i,j) != 0:
                    union(parent, i, j)

    # Class label of each axis is its union-find root.
    return [find(parent, i) for i in 1..3]
```

### 4c.2 Axial count selection (DESIGN 3.7)

Choose the most isotropic integer mesh, per class, that
meets the density floor.  Axes in one class share a count
by construction, so the result is symmetry-compatible.

```
function selectAxialCounts(density, recipMag,
                           recipCellVolume, classes):
    # recipMag(i) = |b_i|.  density is the volume density D.
    # A non-positive density is the Gamma sentinel (3.6),
    # resolved to a single point upstream but guarded here.
    if density <= 0:
        return [1, 1, 1]

    # Continuous isotropic counts at a common spacing h:
    #   h = (prod|b_i| / (recipCellVolume * D))^(1/3).
    h = (recipMag(1) * recipMag(2) * recipMag(3)
         / (recipCellVolume * density)) ^ (1/3)
    n = countsAtSpacing(h, recipMag, classes)

    # Raise WHOLE classes until the full-mesh product meets
    # the floor.  Never raise a single axis inside a multi-
    # axis class -- that would break symmetry compatibility.
    floor = density * recipCellVolume
    while (n(1) * n(2) * n(3)) < floor:
        # Among the classes, pick the one whose increment
        # leaves the three axis spacings |b_i|/n_i closest
        # together (best isotropy).
        bestClass = argmin over classes c of
            spacingSpread(n with class c incremented by 1,
                          recipMag)
        for i in members(bestClass): n(i) = n(i) + 1

    return n

function countsAtSpacing(h, recipMag, classes):
    # The isotropic integer mesh at k-spacing h: one count per
    # |b_i| / h on each axis, shared per class and rounded to a
    # positive integer.  Shared by the density path
    # (selectAxialCounts, above) and the crystalline floor
    # (crystallineFloorMesh, below), so both distribute and round
    # a spacing identically -- the two only differ in the spacing
    # they pass and whether a product-floor bump follows.
    x = [recipMag(i) / h for i in 1..3]
    # Force one shared real count per class.  Coupled axes
    # already have equal |b_i| (hence equal x); the mean
    # guards against round-off before rounding.
    for each distinct class label c:
        members = [i for i in 1..3 if classes(i) == c]
        shared  = mean(x(i) for i in members)
        for i in members: x(i) = shared
    # Nearest positive integer, per class (already equal
    # within a class, so the class stays uniform).
    return [max(1, round(x(i))) for i in 1..3]

function crystallineFloorMesh(recipMag, classes, floorAxisCount):
    # The crystalline climb's opening-floor rung (DESIGN 3.12.4).
    # Set the k-spacing so the DENSEST reciprocal axis (the largest
    # |b_i|) gets exactly floorAxisCount points; every other axis,
    # coarser in reciprocal space, is sampled to that SAME spacing
    # and so gets fewer -- never more -- floored at one point.  A
    # cubic cell floors at [4,4,4]; an anisotropic one at [4,4,2],
    # [4,3,2], [4,1,1], ..., never exceeding the cap on any axis.
    # No product-floor bump follows: the cap is the whole point.
    h = max(recipMag) / floorAxisCount
    return countsAtSpacing(h, recipMag, classes)

function spacingSpread(n, recipMag):
    s = [recipMag(i) / n(i) for i in 1..3]
    return max(s) - min(s)
```

### 4c.3 Shift selection (DESIGN 3.9)

Default to the Gamma-centered mesh, which is invariant
under every point group.  Prefer a half-shift on even,
multi-point axes only where the whole shift vector passes
the invariance test; among the invariant candidates take
the one with the most half-components, which folds deepest.

```
function selectShift(recipPointOps, numPointOps, counts):
    D    = diag(counts)
    Dinv = diag(1 / counts(1), 1 / counts(2), 1 / counts(3))

    # Per-axis offset options.  A half-shift is a candidate
    # only on an axis with an even count > 1 (a single-point
    # axis takes no shift, per 3.6).
    for i = 1 to 3:
        if counts(i) > 1 and counts(i) is even:
            options(i) = {0, 1/2}
        else:
            options(i) = {0}

    best       = [0, 0, 0]     # Gamma-centered: always valid
    bestHalves = 0
    for each s in options(1) x options(2) x options(3):
        if isShiftInvariant(s, D, Dinv,
                            recipPointOps, numPointOps):
            halves = count of components of s equal to 1/2
            if halves > bestHalves:
                best = s; bestHalves = halves
    return best

function isShiftInvariant(s, D, Dinv, recipPointOps,
                          numPointOps):
    # (D M D^{-1} - I) s must be integral for every op M.
    for m = 1 to numPointOps:
        v = (D * recipPointOps(:,:,m) * Dinv - I) * s
        for k = 1 to 3:
            if abs(v(k) - nint(v(k))) > 1e-9:
                return false
    return true
```

### 4c.4 Full Monkhorst-Pack mesh (DESIGN 3.2, 3.9)

The uniform grid in fractional abc coordinates with equal base
weight, measured from the origin: `k = (m + shift)/counts`,
`m = 0..counts-1` per axis (DESIGN 3.9).  With this convention
`shift = 0` contains Gamma and `shift = 1/2` excludes it for
every count, odd or even, and a single point with zero shift
IS Gamma; the classic Monkhorst-Pack grid is the `shift = 1/2`
member.  There is deliberately no `-1/2` offset: points are
compared modulo a reciprocal lattice vector everywhere (4c.5),
so the offset would change no physics for even counts, but it
would swap which shift contains Gamma for odd counts and put a
lone point at the zone corner.  This is the primary product of
the section; 4c.5 optionally folds it.

```
function generateFullMesh(counts, shift):
    delta(i) = 1 / counts(i)
    numFull  = counts(1) * counts(2) * counts(3)
    weightSum = 2                     # 2 electrons/state (3.2)
    baseWeight = weightSum / numFull

    p = 0
    for i = 1 to counts(1):
        for j = 1 to counts(2):
            for k = 1 to counts(3):
                p = p + 1
                mesh(:,p) = ([i, j, k] - 1 + shift) * delta
    return mesh, baseWeight, numFull
```

### 4c.5 Reduction with reciprocal-lattice periodicity
         (DESIGN 3.10)

Fold the full mesh to the IBZ.  A merge is accepted when a
rotated point coincides with a partner MODULO a reciprocal
lattice vector -- i.e. their fractional difference is
integral -- so wrapped coincidences are never missed.  The
merge is exact (a true symmetry gives epsilon(Mk) =
epsilon(k)), so this only shrinks the IBZ, never biases it.
The op-map store (section 5) is retained.

```
function reduceToIBZ(mesh, numFull, baseWeight,
                     recipPointOps, numPointOps):
    identityOpIndex = 1               # guaranteed first (5)
    kpThresh = 1e-5

    tracker = [1, 2, ..., numFull]    # each pt its own rep
    numIBZ  = 0
    for a = 1 to numFull:
        if tracker(a) != a: continue          # already folded
        numIBZ = numIBZ + 1
        tracker(a) = -numIBZ
        ibzPoint(:,numIBZ) = mesh(:,a)
        ibzWeight(numIBZ)  = baseWeight
        fullKPToIBZOpMap(a) = identityOpIndex

        for m = 1 to numPointOps:
            img = recipPointOps(:,:,m) * mesh(:,a)
            for b = a+1 to numFull:
                if tracker(b) != b: continue
                if isPeriodicMatch(img, mesh(:,b), kpThresh):
                    tracker(b) = -numIBZ
                    ibzWeight(numIBZ) += baseWeight
                    fullKPToIBZOpMap(b) = m

    for a = 1 to numFull:
        fullKPToIBZKPMap(a) = -tracker(a)

    # Postconditions (checkable invariants, DESIGN 3.10):
    #   sum of star sizes           == numFull
    #   each star size divides       numPointOps
    #   sum(ibzWeight(1..numIBZ))   == weightSum
    return ibzPoint, ibzWeight, numIBZ,
           fullKPToIBZKPMap, fullKPToIBZOpMap

function isPeriodicMatch(u, v, kpThresh):
    # Coincidence modulo a reciprocal lattice vector: each
    # component of the difference is integral.  nint = round
    # to nearest integer.  Basis-independent, no interval
    # convention needed, correct for non-orthogonal cells.
    for k = 1 to 3:
        d = u(k) - v(k)
        if abs(d - nint(d)) > kpThresh:
            return false
    return true
```

Note on the op map under periodicity: `fullKPToIBZOpMap(b)`
records the operation carrying the IBZ representative to
`b`, now possibly modulo a reciprocal lattice vector.  That
is physically exact -- `k` and `k + G` are the same point --
and the downstream atom permutation (section 6, 7) uses only
the operation index, so nothing else changes.

### 4c.6 Driver: density to resolved mesh

Ties the pieces together and records the resolved mesh so
the convergence ladder (DESIGN 3.11 / 7.8) and the resource
dataspace (DESIGN 8.2, `kpoint_count`) can read what was
actually integrated.  The emission of that record into
`result.toml` is specified with the 7.8 guard, not here.

```
function buildMeshFromDensity(density, shiftRequest,
                              recipMag, recipCellVolume,
                              recipPointOps, numPointOps,
                              reduce):
    classes = computeAxisClasses(recipPointOps, numPointOps)
    counts  = selectAxialCounts(density, recipMag,
                                recipCellVolume, classes)

    if shiftRequest is AUTO:
        shift = selectShift(recipPointOps, numPointOps,
                            counts)
    else:
        shift = shiftRequest          # explicit user override

    mesh, baseWeight, numFull =
        generateFullMesh(counts, shift)

    record resolvedCounts = counts       # 3.11 duplicate guard
    record kpointCount    = numFull      # 8.2 size signature

    if reduce:
        return reduceToIBZ(mesh, numFull, baseWeight,
                           recipPointOps, numPointOps)
    else:
        # Full MP mesh for tetrahedron integration; identity
        # maps make every point its own IBZ representative.
        return mesh, baseWeight, numFull,
               identityMaps(numFull)
```

### 4c.7 Producer-side axis-class sourcing (DESIGN 2.7)

The adaptive climb (4e) runs in the producer, in Python, and
needs a cell's axis `classes` before it dispatches anything --
`initial_meshes` (4e.4) seeds a grid or a descent from the
predicted density, and both call `selectAxialCounts` (4c.2),
which needs `classes`.  imago computes `abcRecipPointOps` and
`computeAxisClasses` at runtime, but the producer cannot wait
for a run to learn how a cell's axes couple.  So the producer
sources the classes itself, from the same two ingredients imago
uses: the space-group operations on disk and the cell the
structure loads in.

This is a Python mirror of the runtime chain
`computeRealPointOps -> computeRecipPointOps ->
computeAxisClasses` (DESIGN 2.7 / 4c.1).  The raw operations
come from `share/spaceDB/<sg>` in conventional-abc fractional
form, through the shared `read_conv_abc_point_ops` reader
(4b.4) -- the SAME parser the kp writer uses, so the classes
the producer seeds from are derived from the operations imago
will run under.  The loaded and conventional lattices come
from the loaded `StructureControl` (`real_lattice`,
`full_cell_real_lattice`), which also already holds the
reciprocal lattice the climb needs for `recipMag` and
`recipCellVolume`.

```
function axisClassesForCell(convAbcPointOps,
                            loadedLattice, convLattice, cellMode):
    # loadedLattice, convLattice: 3x3, lattice vectors as ROWS
    #   (the StructureControl layout), in a common length unit.
    # convAbcPointOps: the rotation matrices R acting as r' = R*r
    #   on conventional-abc fractional coordinates -- iterated
    #   directly, so the count is not a separate argument (this is
    #   a producer-only Python routine, unlike the shared
    #   computeAxisClasses of 4c.1).

    # Change of basis loaded-fractional -> conventional-fractional
    #   (DESIGN 2.7).  Vectors-as-columns are the transpose of the
    #   row-major storage: Lc = transpose(convLattice),
    #   L = transpose(loadedLattice).  In `full` mode L == Lc so
    #   T = I and the conjugation collapses to a copy.
    if cellMode == "full":
        T    = identity(3)
        Tinv = identity(3)
    else:                                   # "prim"
        Lc   = transpose(convLattice)
        L    = transpose(loadedLattice)
        T    = matmul(inverse(Lc), L)       # r_conv = T * r_loaded
        Tinv = inverse(T)

    recipOps = []
    for each Rconv in convAbcPointOps:
        # Direct-space op in the loaded basis: ordinary similarity
        #   (covariant fractional coords), DESIGN 2.7.
        Rloaded = matmul(Tinv, matmul(Rconv, T))
        # Reciprocal-space twin: inverse transpose (contravariant),
        #   DESIGN 2.7.  This is what the mesh folds under.
        Rrecip  = transpose(inverse(Rloaded))
        recipOps.append(Rrecip)

    # 4c.1 union-find over the reciprocal ops' off-diagonal
    #   couplings gives the axis classes.
    return computeAxisClasses(recipOps, len(recipOps))
```

The producer rounds each `Rloaded`/`Rrecip` to the nearest
integer before the union-find: a valid lattice automorphism is
integer with determinant +/-1 (DESIGN 2.7), so the conjugation's
floating-point residue is round-off to be cleaned, and a
non-integer result after rounding-tolerance check is a corrupt
cell or space group and is raised, not silently classed.

**Validation against the runtime.**  The producer's classes must
match the ones imago resolves for the same cell, or the climb
would seed and step a differently-coupled ladder than imago
integrates.  To make the agreement checkable rather than
asserted, imago's mesh emit (4d.5) also prints its resolved
`abcRecipPointOps`-derived class vector as `RESOLVED_KP_CLASSES`;
a producer self-test runs `axisClassesForCell` on a spread of
seed cells (cubic, hexagonal `full`, and a `prim` reduction) and
asserts it reproduces that emitted vector.  The emit is a
validation hook, not a runtime dependency -- the climb never
reads it back, since it needs the classes before the first run.

---

## 4d. Wiring the Mesh Map Through the makeinput ->
        kp-file -> imago Chain (DESIGN 3.2, 3.4, 3.9)

Section 4c gives the mesh map as algorithms; this section
places them across the three stages that produce a k-point
mesh, and specifies that the shift is no longer chosen by
crystal-system name.  The controlling fact: the shift depends
on the parity of the resolved counts (4c.3), and in density
mode the counts are resolved inside imago (4c.2).  So the
automatic shift is resolved in imago, after the counts;
makeinput records only the request.

### 4d.1 makeinput: record the shift request, select nothing

makeinput writes the shift REQUEST into the kp file and does
not choose a shift.  An unset shift is the AUTO sentinel
`-1 -1 -1`; an explicit `-kpshift` is passed through verbatim.
No code in makeinput maps a crystal system to a shift.

```
function shift_request_for_kpfile(settings):
    # DESIGN 3.2 / 3.9.  The shift written to the kp file is the
    # request, not a resolved value.  A whole-Gamma group is
    # handled separately (3.6): its canonical 1x1x1 / "0 0 0"
    # form is written by the Gamma path, never the sentinel.
    return settings.kp_shift        # AUTO sentinel or explicit
```

The AUTO sentinel is distinguishable from any real shift,
which lies in [0, 1); `-1 -1 -1` never denotes an offset.

### 4d.2 imago readKPoints: parse the request, flag AUTO

```
function readKPointShift(file):
    read 'KP_SHIFT_A_B_C'
    read kPointShift(1..3)                    # three floats
    # The sentinel is a request to select later (4c.3), not a
    # usable offset.  Flag it; resolve it in initializeKPoints
    # once the counts are known.
    isAutoShift = all(kPointShift(i) == -1 for i in 1..3)
```

### 4d.3 imago initializeKPoints: counts, then shift, then mesh

This realizes the driver of 4c.6 in imago's module-state
style, wired to the two mesh style codes.  The essential
change from the earlier flow is order: the reciprocal
operations are formed FIRST, because the axis classes and the
shift selection both need them, and the count selection needs
the classes.

```
# Style code 2 (density mode).  DESIGN 3.2 step 3.
elif kPointStyleCode == 2:
    computeRealPointOps()
    computeRecipPointOps()
    classes = computeAxisClasses(abcRecipPointOps,          # 4c.1
                                 numPointOps)
    numAxialKPoints = selectAxialCounts(                   # 4c.2
        minKPointDensity, recipMag, recipCellVolume, classes)
    resolveShift(classes)
    initializeKPointMesh(applySymmetry)     # 4c.4 + 4c.5

# Style code 1 (explicit mesh).  numAxialKPoints came from the
# file; the shift may still be an AUTO request.
elif kPointStyleCode == 1:
    computeRealPointOps()
    computeRecipPointOps()
    classes = computeAxisClasses(abcRecipPointOps, numPointOps)
    resolveShift(classes)
    initializeKPointMesh(applySymmetry)


function resolveShift(classes):
    # DESIGN 3.9 / 4c.3.  Resolve an AUTO request from the
    # counts; else honor the explicit shift, zeroing any single-
    # point axis (3.6) and warning if it is not invariant -- an
    # override is the user's prerogative, but a non-invariant
    # shift reduces poorly.
    if isAutoShift:
        kPointShift = selectShift(abcRecipPointOps,
                                  numPointOps, numAxialKPoints)
    else:
        for i in 1..3:
            if numAxialKPoints(i) == 1:
                kPointShift(i) = 0
        if not isShiftInvariant(kPointShift,
                diag(numAxialKPoints), diag(1 / numAxialKPoints),
                abcRecipPointOps, numPointOps):
            warn("explicit k-point shift is not invariant under"
                 + " the cell's point group; symmetry reduction"
                 + " may be incomplete (DESIGN 3.9)")
```

`initializeKPointMesh` is the 4c.4 + 4c.5 routine: it builds
the full Monkhorst-Pack mesh, then folds it to the IBZ with
the reciprocal-lattice-periodic match when `applySymmetry`
is 1.  The full mesh is retained (`numFullMeshKP`,
`fullKPToIBZKPMap`) so tetrahedron integration (LAT) tiles the
whole mesh even when the eigenvalue problem uses the IBZ --
this is how one built mesh serves both the reduced and the
full (LAT) consumers.

### 4d.4 Parallel consumer: the makeKPoints executable

The standalone `makeKPoints` program (`makekpoints.F90`, its
own `foldMesh`) pre-generates explicit k-point lists on the
non-density path, and folds with the same raw comparison that
4c.5 replaces.  It is a second consumer of the reduction and
needs the same reciprocal-lattice-periodic match, so it does
not under-reduce non-orthogonal cells while imago's own path
does.  Its fold adopts `isPeriodicMatch` (4c.5) identically;
the count and shift selection (4c.1-4c.3) do not apply there,
since that path receives explicit counts and shift.

### 4d.5 Emit the resolved mesh to the run output

Once `initializeKPoints` has resolved the mesh --
`numAxialKPoints` set, and `numKPoints` after the fold --
imago writes two labeled records to the main run output (the
settled SCF output on Fortran unit 20, `gs_scf-fb.out`), so
imago.py can recover the mesh without re-deriving it
(DESIGN 6.1.2; the parse is 12.5).  The resolved mesh is a
per-run fact, constant across SCF cycles, so it does NOT
belong in the per-cycle iteration file where the energy and
gap live; a labeled record in the settled output is the robust
home (DESIGN 6.1.6 prefers settled files to stdout scraping).

Each record follows imago's own label/value convention: the
all-caps tag on its own line, the value on the next line (as
`KPOINT_STYLE_CODE` and the other kp tags are written and
echoed).

```
# In initializeKPoints.  builtAxialMesh is a local logical,
# false on entry, set true ONLY inside the style-code 1 and 2
# branches, which are the branches that build a uniform mesh
# AND compute the axis classes.  The file's style code alone
# is not the test: the SYBD and MTOP branches take precedence
# over it and build their own k-point sets without touching
# axisClass, so a record keyed on the style code would report
# a class vector nobody computed.
if builtAxialMesh:
    write unit 20: "RESOLVED_KP_MESH"
    write unit 20: numAxialKPoints(1), numAxialKPoints(2),
                   numAxialKPoints(3)
    write unit 20: "RESOLVED_KP_COUNT"
    write unit 20: numKPoints
    # Validation hook for the producer's axis-class port (4c.7):
    #   emit the resolved axis-class label of each reciprocal
    #   axis, so a producer self-test can assert its own
    #   axisClassesForCell reproduces what imago resolved.
    write unit 20: "RESOLVED_KP_CLASSES"
    write unit 20: axisClass(1), axisClass(2), axisClass(3)
```

`RESOLVED_KP_MESH` is the full uniform mesh's axial counts;
`RESOLVED_KP_COUNT` is the number of k-points actually
computed -- the IBZ size when symmetry was applied, the full-
mesh size otherwise.  `RESOLVED_KP_CLASSES` is the axis-class
label vector (`axisClass`, filled by `computeAxisClasses`,
4c.1) -- three integers whose equality pattern says which axes
the point group couples; it exists to validate the producer's
4c.7 mirror, not for the run itself.  Only the mesh style-code
branches (1 and 2) emit these.  An explicit-list run (style 0)
builds no axial mesh and emits none; a band-structure (SYBD)
or polarization (MTOP) run replaces the loaded set with its
own path or full string mesh (DESIGN 2.6) and emits none
either.  In every such case imago.py records them as absent
(6.1.2).

---

## 4e. Adaptive Mesh Climb (DESIGN 3.12)

The convergence search, in mesh space.  It reuses the mesh
selection of 4c.2 as its rung primitive and replaces the fixed
verify grid the producer built for predict-then-verify (11.4 /
7.7).  The producer predicts and records in density, but
searches by climbing meshes: each rung is a distinct symmetry-
compatible mesh, and the climb stops when the energy is flat.
Two dispatch modes, gated by the prediction's confidence, share
one rung rule and one stop test.

Everything below is producer-side (Python).  The mesh primitives
`selectAxialCounts` and `spacingSpread` are the 4c.2 functions,
re-expressed here on an explicit counts vector rather than the
imago module state, since the producer reasons about meshes for
many materials at once.  `recipMag`, `recipCellVolume`, and the
axis `classes` for a material are computed once from its cell
(4c.1) and carried in its `config`.

### 4e.1 Rung mechanics

```
function climbOneRung(counts, classes, recipMag):
    # The next distinct mesh up the climb: one step of the 4c.2
    # floor loop -- increment the axis class that most evens the
    # three inter-point spacings (DESIGN 3.12.2).
    bestClass = None
    bestSpread = +infinity
    for c in distinct(classes):
        trial = bump(counts, classes, c, +1)
        s = spacingSpread(trial, recipMag)          # 4c.2
        if s < bestSpread:
            bestSpread = s
            bestClass = c
    return bump(counts, classes, bestClass, +1)

function climbNRungs(counts, n, classes, recipMag):
    # Advance `n` rungs up the ladder: the rung rule applied n times
    # (DESIGN 3.12.2, the stride).  n = 1 is climbOneRung.  A stride
    # crosses n ladder positions for the cost of one calculation --
    # the bracket phase (4e.3) uses geometric strides so an unknown
    # convergence distance is bracketed in log steps.
    mesh = counts
    repeat n times:
        mesh = climbOneRung(mesh, classes, recipMag)
    return mesh

function descendOneRung(counts, classes, recipMag):
    # A lower mesh that climbOneRung steps back up to `counts`
    # (DESIGN 3.12.4).  A mesh can have more than one such lower
    # neighbour (e.g. hexagonal [4,4,2] is reached from both
    # [4,4,1] and [3,3,2]); the loop returns the first in class
    # order, so climbOneRung(descendOneRung(counts)) == counts
    # always holds, but descendOneRung need not undo a particular
    # climb path.  Returns `counts` unchanged when it is already
    # minimal for the cell.
    for c in distinct(classes):
        if any axis of class c in `counts` equals 1:
            continue                       # cannot go below 1
        trial = bump(counts, classes, c, -1)
        if climbOneRung(trial, classes, recipMag) == counts:
            return trial
    return counts

function bump(counts, classes, c, step):
    # Add `step` to every axis whose class is c (keeps a class
    # equal, so the mesh stays symmetry-compatible, DESIGN 3.8).
    return [n + step if classes[i] == c else n
            for i, n in enumerate(counts)]
```

### 4e.2 The stop test and the ceiling

```
function pick_converged_climb(rungs, cell_atom_count, threshold,
                              flat_needed):
    # rungs: the distinct meshes run so far, sorted ascending, each
    # a {mesh, energy}.  Return the index of the converged rung, or
    # None to keep climbing.  A rung converges when its per-atom
    # energy is within `threshold` of BOTH neighbours (the 7.8 step
    # 3c two-sided test), and that flatness PERSISTS over
    # `flat_needed` consecutive interior rungs (DESIGN 3.12.3):
    # flat_needed is 1 for a confident search, 2 for cold/bootstrap
    # (Q1).  With flat_needed = 1 this is exactly pick_converged.
    e = [per_atom_ev(r.energy, cell_atom_count) for r in rungs]

    def two_sided_flat(j):
        return (abs(e[j] - e[j - 1]) < threshold and
                abs(e[j] - e[j + 1]) < threshold)

    for i in range(1, len(e) - 1):
        top = i + flat_needed - 1          # last interior to check
        if top > len(e) - 2:               # not enough rungs above
            break                          #   to confirm i yet
        if all(two_sided_flat(j) for j in range(i, top + 1)):
            return i
    return None


function stride_is_flat(lo_rung, hi_rung, cell_atom_count, threshold):
    # The bracket flatness test (DESIGN 3.12.3): a stride is flat when
    # its two endpoints' per-atom energies are within `threshold`.  The
    # caller passes the LOOSER bracket threshold here, not the strict
    # convergence threshold (config.stride_threshold, 4e.3): a stride
    # that has nearly settled brackets one geometric step sooner,
    # shaving the far larger endpoint the strict threshold would demand.
    # Because a stride adds many k-points, a small change across it is
    # strong evidence the energy has settled -- but only evidence: the
    # refine phase VERIFIES with the full two-sided test at the strict
    # threshold, so a stride that reads loosely flat but has not truly
    # converged (an oscillating near-metal, or a nearly-settled stride
    # above the real convergence) is caught there, not trusted here.
    lo = per_atom_ev(lo_rung.energy, cell_atom_count)
    hi = per_atom_ev(hi_rung.energy, cell_atom_count)
    return abs(hi - lo) < threshold


function is_gapless(rung, gap_threshold):
    # The metal test (DESIGN 3.12.3): true when a rung's computed
    # band gap is essentially zero -- at or below `gap_threshold`
    # (config.metal_gap_threshold, 4e.3), an eV value low enough that
    # no real insulator crosses it yet high enough to catch a true
    # metal's near-zero reading.  A metal has no gap: its energy
    # oscillates as the mesh crosses the Fermi surface and never
    # settles, so chasing k-point convergence on it is futile.  The
    # gap is read straight from the rung's result -- a DIRECT metal
    # signal, unlike the retired proxy that inferred metallicity from
    # a finer mesh raising the energy and so missed the common
    # small-amplitude oscillator.  This judges ONE rung, not a stride.
    #
    # A thin wrapper over the scalar rule (is_gapless_value, 15.7):
    # this side knows how to find a gap on a RUNG, that side holds
    # what the gap means -- including that an unknown gap (None) is
    # NOT metallic, so a missing reading never stops a climb that was
    # converging.  The guidance harvest calls the same core on a
    # result dict, so the climb's metal short-circuit and the
    # harvest's metal skip cannot drift apart (DESIGN 7.8).
    return is_gapless_value(rung.gap, gap_threshold)


function at_ceiling(mesh, max_count):
    # The fixed per-axis backstop (DESIGN 3.12.3): the LARGEST axial
    # count reaching max_count, not the product or the stride.  A cost
    # ceiling from the resource dataspace (16) layers on later; the
    # climb stops at whichever bites first.
    return max(mesh) >= max_count
```

### 4e.3 One material's next action

Three search shapes (DESIGN 3.12.3 / 3.12.5) share the two-sided
stop test (4e.2), the metal test (4e.2) and the rung rule (4e.1);
they differ only in which rungs they compute.  `climb_next`
dispatches on the mode, threading a per-material search `state`
for the stateful bracket-refine shape; the grid and the unit-step
climb ignore it.

The metal test sits in `climb_next` itself, ABOVE the dispatch,
which is what makes it shared rather than replicated (DESIGN
3.12.3).  Recognising a metal is a classification, and a
classification cannot belong to a search shape: the shapes exist
to disagree about which rungs are worth computing, not about what
a computed rung means.  Placing it here also means there is
exactly one copy of the rule, so no shape can grow its own
variant of it.
Every shape returns one of: `RUN(mesh)` -- run one more mesh;
`CONVERGED(rung)` -- done, an insulator settled at that rung;
`METAL(rung)` -- stop, a metal recognised by its vanishing gap
(4e.2) and settled at that rung as a deliberately rough starting
potential; `CEILING` -- stop, non-converged (a hard insulator that
hit the count ceiling still steep).  `CONVERGED` and `METAL` both
RECORD a result rung -- the metal's is rough by intent, not a
k-converged energy -- and differ only in the reason tag; `CEILING`
alone is the non-converged verdict the producer (4e.5) retires.

```
function climb_next(rungs, state, config):
    # rungs: the material's computed {mesh, energy, gap}, sorted
    #   ascending; gap is the rung's band gap (eV), read for the
    #   metal test (4e.2), alongside the energy the stop test reads.
    # state: the per-material search state (bracket-refine only).
    # Returns (action, state').
    #
    # The metal test, ahead of every shape's own logic and ahead of
    #   any convergence work (DESIGN 3.12.3): if any computed rung
    #   reads gapless, stop and settle ON THAT RUNG -- a rough metal
    #   potential, not a k-converged energy.  A metal reads gap ~ 0
    #   from the floor up, so this usually fires on the opening rung;
    #   a near-metal that shows a small gap coarse and closes it
    #   finer fires on the rung where the gap first vanishes.
    #
    # rungs is sorted ascending by mesh, so scanning it in order and
    #   taking the FIRST match settles on the COARSEST gapless rung.
    #   That is the rule DESIGN 3.12.3 states, and it is not the same
    #   as the densest rung on the ladder: a confident opening grid
    #   resolves several rungs at once and a refine fill lands below
    #   rungs already computed, so the gapless rung need not be the
    #   last one.  Settling on the densest instead could settle on a
    #   rung that read a GAP, and the harvest re-reads that single
    #   rung (7.8 / 15.7) -- it would then see an insulator where
    #   this test saw a metal.  Taking the gapless rung makes the two
    #   agree by construction rather than by coincidence.
    #
    # A NEGATIVE metal_gap_threshold can never fire, since no band
    #   gap is negative -- that is how a curator asks for a known
    #   metal's every rung to be computed (DESIGN 3.12.3 / 3.12.6).
    #   No special case is needed for it here; it falls out of the
    #   comparison.
    metal_rung = first r in rungs with
        is_gapless(r, config.metal_gap_threshold)
    if metal_rung is not None:
        return METAL(metal_rung), state
    if config.mode == BRACKET_REFINE:
        return bracketRefineNext(rungs, state, config)
    # UNIT_STEP, or a GRID opening that did not converge and now
    #   continues as a unit-step climb (4e.4): stateless in rungs.
    return climbAction(rungs, config), state
```

```
function climbAction(rungs, config):
    # The unit-step climb (DESIGN 3.12.3): walk one rung at a time,
    # testing the whole accumulated ladder.  RUN one more, CONVERGED
    # at the flat interior rung, or CEILING.  It returns no METAL of
    # its own -- climb_next tests the gap before dispatching here, so
    # a metal never reaches this function.
    idx = pick_converged_climb(
        rungs, config.cell_atom_count, config.threshold,
        config.flat_needed)
    if idx is not None:
        return CONVERGED(rungs[idx])
    if at_ceiling(rungs[-1].mesh, config.max_count):
        return CEILING
    return RUN(climbOneRung(rungs[-1].mesh, config.classes,
                            config.recipMag))
```

The **bracket-refine** search (DESIGN 3.12.3) is a two-phase state
machine.  Its per-material `state`:

```
# newBracketRefineState(seed) = {
#   phase:     BRACKET,     # BRACKET | REFINE
#   stride:    1,           # current geometric stride (BRACKET)
#   endpoints: [seed],      # bracket endpoint meshes so far,
#                           #   ascending; endpoints[0] is the seed,
#                           #   endpoints[-1] the newest
#   lo: None, hi: None,     # the interval REFINE fills, as meshes
#   from_cap:  False }      # True iff this bracket runs up to the
#                           #   ceiling, so an empty refine is a
#                           #   CEILING stop, not a false bracket
```

```
function bracketRefineNext(rungs, state, config):
    # No metal test here: climb_next above applies it to every shape
    #   before dispatching, so by the time this runs no rung on the
    #   ladder is gapless.  The recursive resume at the foot of this
    #   function may therefore skip it too -- the rungs it re-judges
    #   are the same ones climb_next already read.
    if state.phase == BRACKET:
        top = state.endpoints[-1]
        if length(state.endpoints) == 1:
            # Only the seed is computed; launch the first stride
            #   (stride 1).  No flatness to test yet.
            return strideUp(state, 1, config)
        # Two or more endpoints.  The bracket test uses the LOOSER
        #   stride_threshold (4e.2); the refine below keeps the strict
        #   convergence threshold.  (No near-metal check here: a metal
        #   is caught by climb_next's gap test before this shape is
        #   dispatched at all, so before any stride is judged.)
        prev = state.endpoints[-2]
        if stride_is_flat(rung_at(rungs, prev), rung_at(rungs, top),
                          config.cell_atom_count,
                          config.stride_threshold):
            # First flat stride.  The converged rung lies at or just
            #   above the bottom of the flat stride, prev.  Fill
            #   flat_needed + 1 rungs above prev, so the persistence
            #   test has flat_needed interior candidates each with a
            #   computed neighbour on both sides (the two-sided test
            #   excludes a block's endpoints).  The + 1 matters because
            #   prev itself need not be settled: its own lower
            #   neighbour may still be moving, so the first rung the
            #   test can confirm is often prev + 1, and confirming
            #   flat_needed rungs from there needs a computed neighbour
            #   up through prev + flat_needed + 1.  Filling fewer would
            #   expose too few interior candidates, which a two-
            #   consecutive-flat search could never confirm.  If prev
            #   is the seed (the very first stride was flat) there is
            #   no lower endpoint either, so lo extends one rung below
            #   prev.
            hi = climbNRungs(prev, config.flat_needed + 1,
                             config.classes, config.recipMag)
            if length(state.endpoints) >= 3:
                lo = state.endpoints[-3]
            else:
                lo = descendOneRung(prev, config.classes,
                                    config.recipMag)
            return enterRefine(rungs, state, lo, hi, False, config)
        # Not flat: grow the stride geometrically and step up, unless
        #   the next endpoint would pass the ceiling -- then refine
        #   from the top endpoint up to the ceiling (DESIGN 3.12.3),
        #   so a convergence a stride jumped over just below the cap
        #   is still found.
        nextStride = min(2 * state.stride, config.max_stride)
        nxt = climbNRungs(top, nextStride, config.classes,
                          config.recipMag)
        if at_ceiling(nxt, config.max_count):
            return enterRefine(rungs, state, top,
                               ceilingMesh(top, config), True, config)
        return strideUp(state, nextStride, config)

    else:  # REFINE -- fill [lo, hi] lowest-first, TESTING as we go
        # Judge after each fill so a convergence low in the bracket
        #   stops the fill before the wide rungs above it are computed
        #   (DESIGN 3.12.3).  Test only the CONSECUTIVE block anchored
        #   at lo -- the run of computed rungs each one step up from
        #   the last -- not the whole [lo, hi] range: mid-fill the
        #   range still has gaps (a sparse bracket endpoint sitting
        #   above an unfilled rung), and the two-sided test would
        #   compare non-neighbours across such a gap and could read a
        #   false convergence.  lo is always computed by the time the
        #   first fill lands (nextFillMesh fills it first), so the
        #   anchor is safe.  Because the fill climbs from the bottom
        #   and this returns the SMALLEST passing rung, the mesh it
        #   converges on is exactly the one a full fill would have --
        #   only the rungs above it go uncomputed.
        block = consecutive_block(rungs, rung_at(rungs, state.lo),
                                  config.classes, config.recipMag)
        idx = pick_converged_climb(block, config.cell_atom_count,
                                   config.threshold, config.flat_needed)
        if idx is not None:
            return CONVERGED(block[idx]), state
        # Not verified yet: fill the next-lowest gap if any remains.
        gap = nextFillMesh(rungs, state.lo, state.hi,
                           config.classes, config.recipMag)
        if gap is not None:
            return RUN(gap), state              # keep filling
        # Interval fully filled and still nothing verified.
        if state.from_cap:
            return CEILING, state               # steep even at the cap
        # A coincidentally flat stride (an oscillating energy): no
        #   rung verified.  Resume striding from hi (DESIGN 3.12.3).
        return bracketRefineNext(
            rungs, newBracketRefineState(state.hi), config)
```

```
function strideUp(state, stride, config):
    # Launch the next bracket endpoint `stride` rungs above the top,
    # recording it on the endpoint list (DESIGN 3.12.2 / 3.12.3).
    nxt = climbNRungs(state.endpoints[-1], stride, config.classes,
                      config.recipMag)
    state' = state with endpoints = state.endpoints + [nxt],
             stride = stride
    return RUN(nxt), state'

function enterRefine(rungs, state, lo, hi, from_cap, config):
    # Switch to REFINE on [lo, hi] and launch the first fill mesh; if
    # the interval is already fully computed (its endpoints and
    # nothing between), judge it at once (DESIGN 3.12.3).
    state' = state with phase = REFINE, lo = lo, hi = hi,
             from_cap = from_cap
    gap = nextFillMesh(rungs, lo, hi, config.classes, config.recipMag)
    if gap is None:
        return bracketRefineNext(rungs, state', config)
    return RUN(gap), state'

function nextFillMesh(rungs, lo, hi, classes, recipMag):
    # The lowest ladder position in [lo, hi] not yet computed, or None
    # when the interval is fully filled (DESIGN 3.12.3, the fill).
    mesh = lo
    loop:
        if mesh not among [r.mesh for r in rungs]:
            return mesh
        if mesh == hi:
            return None
        mesh = climbOneRung(mesh, classes, recipMag)

function ceilingMesh(fromMesh, config):
    # The first mesh at or above the per-axis ceiling, climbing from
    # `fromMesh` -- the upper bound of an up-to-ceiling refine (3.12.3).
    mesh = fromMesh
    while not at_ceiling(mesh, config.max_count):
        mesh = climbOneRung(mesh, config.classes, config.recipMag)
    return mesh

function rung_at(rungs, mesh):
    # The computed {mesh, energy} at `mesh` (a bracket endpoint is
    # always already computed when its stride is tested).
    return the r in rungs with r.mesh == mesh
```

`meshSize(mesh)` is the full-mesh point count (product of the axial
counts), the same monotone key the ladder is sorted by (4e.6), so a
point-count range selects exactly the consecutive block the fill
produced.

### 4e.4 Seeding, and the three search modes

```
function initial_meshes(density, config):
    # The mesh(es) to run in the first round (DESIGN 3.12.4/3.12.5).
    # The seed DENSITY -- the prediction, or the wide-grid floor when
    # under-trained (7.9) -- picks a mesh via 4c.2; the mode decides
    # whether the first round is a parallel grid or a single rung.
    # (The producer resolves the confidence policy into `config`
    # first, so seeding takes only the density, not the prediction.)
    seed = selectAxialCounts(density, config.recipMag,
                             config.recipCellVolume, config.classes)

    if config.mode == PARALLEL_GRID:        # confident (Q3)
        # The seed plus grid_width rungs on each side, all laid
        # down in the first round and judged as one grid.
        meshes = [seed]
        lo = seed
        hi = seed
        repeat config.grid_width times:
            lo = descendOneRung(lo, config.classes, config.recipMag)
            hi = climbOneRung(hi, config.classes, config.recipMag)
            meshes = [lo] + meshes + [hi]
        return distinct(meshes)

    else:                     # BRACKET_REFINE or UNIT_STEP: a climb
        # Begin below the prediction so the seed acquires a lower
        # neighbour, then climb upward.  start_offset grows as
        # confidence falls (a weaker prediction starts lower).  Both
        # climbs open on this single rung; they differ only in the
        # action (4e.3) -- bracket-refine strides, unit-step walks.
        start = seed
        repeat config.start_offset times:
            start = descendOneRung(start, config.classes,
                                   config.recipMag)
        # Crystalline floor (DESIGN 3.12.4): no crystalline material
        # converges on a mesh coarser than the floor rung, so a cold
        # bootstrap that would open below it is lifted up to it -- out
        # of the unreliable coarse regime the near-metal bail (4e.3)
        # must otherwise guard against.  config.opening_floor is the
        # floor mesh for a crystalline solid (built once in the config,
        # 11.4) or None for a non-crystalline one, which seeds at or
        # near Gamma by convention (7.9) and must not be floored up.  A
        # warm seed already above the floor is untouched.
        if config.opening_floor is not None
           and meshSize(config.opening_floor) > meshSize(start):
            start = config.opening_floor
        return [start]


function newSearchState(config, opening_meshes):
    # The per-material search state converge_by_climb (4e.5) threads.
    # Only the bracket-refine climb is stateful (4e.3); the grid and
    # the unit-step climb carry an empty state they never read.  A
    # climb's seed is its single opening mesh; a grid seeds no state.
    if config.mode == BRACKET_REFINE:
        return newBracketRefineState(opening_meshes[0])
    return EMPTY_STATE
```

`config.mode`, `config.flat_needed`, `config.grid_width`,
`config.start_offset`, and -- for the bracket-refine climb --
`config.max_stride` are set from `prediction.confidence` and the
curator's chosen climb shape by a small policy the implementation
tunes (DESIGN 3.12.6): high confidence selects `PARALLEL_GRID`,
`flat_needed = 1`, a narrow grid; low or under-trained confidence
selects a climb -- `BRACKET_REFINE` by default, or `UNIT_STEP`
when the curator pins the fine shape -- with `flat_needed = 2`
and a lower start.  An under-trained prediction (7.6) starts the
climb from the wide-grid floor rather than a predicted seed (7.9).
`max_stride` caps the geometric stride, which keeps the bracket a
refine has to fill small (4e.3 / DESIGN 3.12.3).

`config.stride_threshold` -- the LOOSER threshold the bracket test
uses (4e.2 / 4e.3) -- is assembled where the producer builds the
config (11.4 / build_climb_config): it is the solid's strict
`kpoint_convergence_threshold` multiplied by the database-wide
`stride_flatness_multiple` knob (>= 1, default a small multiple to
be fixed by experiment, DESIGN 3.12.6).  The strict threshold
stays on `config.threshold` for the refine; only the bracket
stride test reads the looser one.  `config.metal_gap_threshold` --
the eV band gap at or below which a rung reads as gapless (4e.2 /
4e.3) -- is the database-wide `metal_gap_threshold` knob taken
directly: an absolute gap in eV, not a multiple of the convergence
threshold (DESIGN 3.12.3 / 3.12.6).  `config.opening_floor` -- the
crystalline
climb's floor rung (4e.4 / DESIGN 3.12.4) -- is built once where the
producer builds the config (11.4 / build_climb_config): for a
crystalline solid it is `crystallineFloorMesh`, whose densest axis
holds `crystalline_floor_axis_count` points and whose other axes
scale down by their `|b_i|` ratios, so it is `[4,4,4]` on a cubic
cell and never exceeds the cap on any axis of an anisotropic one;
for a non-crystalline solid it is None.  Because the climb never
opens below that rung, the gap test reads a trustworthy gap at
every rung and needs no coarse-mesh guard of its own (4e.3).

The numeric knobs that policy reads -- the `confidence_high`
threshold, the two `flat_needed` counts, `grid_width`, the two
`start_offset` values, `max_stride`, the climb-shape choice, the
`stride_flatness_multiple`, the `metal_gap_threshold`, the
`crystalline_floor_axis_count`, and the per-axis `max_count`
ceiling -- are config, not constants (Principle 11).  They are sourced from the
manifest `[harvest.kpoint_climb]` sub-table (DESIGN 5.7 / 3.12.6),
each knob falling back to a documented provisional default when the
sub-table or that knob is omitted.  The producer resolves them once
per run by merging the sub-table over the provisional defaults into
a `PolicyThresholds` bundle plus `max_count`; the manifest reader
validates the sub-table's keys against the known knob names
(`KPOINT_CLIMB_KEYS`), so a mistyped knob fails loudly at load
rather than silently taking a default.  The one knob with a
restricted VALUE, `climb_shape`, is checked too: the merge rejects
any value that is not one of the known climb shapes (`BRACKET_REFINE`
or `UNIT_STEP`), so a typo like `"unit-step"` fails loudly rather
than falling through to a default shape.  The
`stride_flatness_multiple` is likewise checked to be `>= 1`, since
a value below one would invert its meaning (a stricter-than-
convergence bracket); and `crystalline_floor_axis_count` is
checked `>= 1` as a per-axis count.  `metal_gap_threshold` is
NOT range-checked: any real value is meaningful.  It is an
absolute band gap in eV, and a negative one is the documented way
to disable the metal test for a diagnostic ladder, since no band
gap is negative (DESIGN 3.12.3 / 3.12.6).  A `> 0` check here
would reject exactly the setting the design tells a curator to
use.  The provisional default
values themselves are still to be fixed by the seed experiment
(3.12.6).

### 4e.5 Concurrent orchestration across materials

```
function converge_by_climb(materials, configs, seed_densities,
                           dispatcher, on_non_converged):
    # Drive every material through the climb to a verdict -- one of
    # CONVERGED (the energy went flat), METAL (a rung read gapless,
    # 3.12.3) or NOT_CONVERGED (a ceiling, or a rung that failed to
    # run).  Returns (outcomes, rungs, verdicts): outcomes[m] is the
    # settled Rung or the NON_CONVERGED sentinel, verdicts[m] is the
    # REASON, kept because CONVERGED and METAL both produce a settled
    # rung and are NOT interchangeable downstream (DESIGN 5.7).
    # Discarding the reason here is what forced every later stage to
    # re-derive the classification from whatever evidence it happened
    # to hold -- and the harvest holds only ONE rung, whose apparent
    # gap on a discrete mesh is close to a coin toss (DESIGN 1.6).
    # Serial within a material, concurrent across, and NO
    # material waits on another: a chain climbs on the instant its
    # own rung lands (DESIGN 3.12.5).  The injected `dispatcher`
    # owns the in-flight set so this loop tracks only its per-
    # material ladders (Principle 12); it exposes two calls (4e.7):
    #   dispatcher.send(mesh_lists) -- launch one calc per (material,
    #     mesh) WITHOUT waiting (send_off, 13.5).
    #   dispatcher.next_rung() -- block until the next rung lands and
    #     return (material, result), where result is a {mesh, energy,
    #     gap}
    #     Rung or the FAILED marker for a rung that did not complete
    #     (7.7).
    # seed_densities[m] is m's round-0 seed density (the prediction,
    # or the wide-grid floor when under-trained, 3.12.4 / 7.9); the
    # confidence policy is already resolved into configs[m].
    # on_non_converged(m) tags a prediction mismatch (7.8 3d); it is
    # injected so this loop stays free of the workspace.

    rungs    = { m: [] for m in materials }
    search   = {}                  # per-material search state (4e.4)
    outcomes = {}
    verdicts = {}                  # m -> CONVERGED | METAL |
                                   #      NOT_CONVERGED
    active   = set(materials)
    in_air   = {}                  # rungs still in flight, per m
    opening  = set(materials)      # still in the opening (grid) phase

    # retire m with a settled rung (or the NON_CONVERGED sentinel)
    #   AND its reason, then drop it from the active set.  One place
    #   writes both, so an outcome can never be recorded without the
    #   reason that produced it.  A non-converged stop also tags the
    #   prediction mismatch (7.8 3d).
    function retire(m, outcome, verdict):
        outcomes[m] = outcome
        verdicts[m] = verdict
        if verdict is NOT_CONVERGED:
            on_non_converged(m)
        active.discard(m)

    # judge m's ladder (climb_next, 4e.3) and either retire it or
    #   launch its single next rung.  climb_next threads the per-
    #   material search state, so the bracket-refine phase persists
    #   across landings; the grid and unit-step climbs ignore it.
    function judge(m):
        (action, search[m]) = climb_next(rungs[m], search[m],
                                          configs[m])
        if action is CONVERGED or action is METAL:
            # Both record a settled rung -- a METAL's is a rough metal
            #   potential (DESIGN 3.12.3), a CONVERGED's a k-converged
            #   insulator energy -- and both leave active with no
            #   mismatch tag.  They are told apart by the verdict, and
            #   ONLY by it: nothing about the rung itself says which
            #   kind of stop produced it.
            retire(m, action.rung,
                   METAL if action is METAL else CONVERGED)
        elif action is CEILING:
            retire(m, NON_CONVERGED, NOT_CONVERGED)  # 7.8 3d; a hard
            #   insulator still steep at the count ceiling, dropped.
        else:                                        # RUN(mesh)
            dispatcher.send({ m: [action.mesh] })
            in_air[m] += 1

    # Seed every material's opening rung or grid at once
    #   (initial_meshes, 4e.4): one rung for a climb, a small grid
    #   for the confident mode (3.12.6).
    first = { m: initial_meshes(seed_densities[m], configs[m])
              for m in materials }
    for m in materials:
        # The search state (4e.4) is stateful only for the bracket-
        #   refine climb, seeded from its opening rung.
        search[m] = newSearchState(configs[m], first[m])
        in_air[m] = length(first[m])
    dispatcher.send(first)

    # Collect rungs as they land, in landing order, until nothing is
    #   in flight.  Each landing advances exactly the one material it
    #   belongs to; the others are untouched, so no chain is paced by
    #   another.
    while sum(in_air) > 0:
        (m, result) = dispatcher.next_rung()
        in_air[m] -= 1
        if result is not FAILED:
            rungs[m] = merge_distinct(rungs[m], [result])

        if m in opening:
            # The confident mode's opening grid is judged as a group,
            #   so wait until the WHOLE grid has resolved, then judge
            #   on whatever landed.  A material whose entire opening
            #   failed has no rung to stand on (run failure, 7.7).  A
            #   climb's opening is a single rung, so it is judged at
            #   once.
            if in_air[m] > 0:
                continue
            opening.discard(m)
            if rungs[m] is empty:
                retire(m, NON_CONVERGED,
                       NOT_CONVERGED)                # run failure
            else:
                judge(m)
        else:
            # A continuation is exactly one rung.  If it failed to
            #   run the climb cannot advance, so stop the material
            #   rather than re-dispatch it forever (7.7); otherwise
            #   judge the extended ladder.
            if result is FAILED:
                retire(m, NON_CONVERGED,
                       NOT_CONVERGED)                # run failure
            else:
                judge(m)

    return outcomes, rungs, verdicts
```

No chain ever waits on another: a material advances the instant
its own rung lands and drops out of `active` the moment it
converges or hits its ceiling, so a late, expensive chain never
holds back the ones already done (DESIGN 3.12.5).  A material has
at most its opening grid in flight and, after that, exactly one
rung, so once a material is judged nothing further lands for it.
The dispatch core (13) stays domain-ignorant: it launches the
meshes handed to it and reports energies as they finish; the
choice of the next mesh lives entirely here.

### 4e.6 Recording the converged rung

```
function record_converged(rung, rungs, config):
    # Build the DENSITY / MESH / GRID harvest inputs for a converged
    # CLIMB material; build_entry (15.7) adds the gap, magnetization,
    # sub-model, and provenance around them.  `rung` is the converged
    # rung, `rungs` its ascending distinct-mesh ladder.  The dataspace
    # key is a DENSITY; the mesh is stored exact alongside it (DESIGN
    # 3.12.4 / Q4).  The density a mesh represents is its full-mesh
    # volume density, product(mesh) / recipCellVolume -- self-
    # consistent with 4c.2, so a future prediction of this density
    # reproduces this mesh in this cell.  (m and prediction are the
    # caller's to thread into build_entry; record_converged needs
    # only the rungs and the cell's reciprocal volume.)
    converged_density = product(rung.mesh) / config.recipCellVolume
    # The stored flatness trace is the CONSECUTIVE block of the ladder
    # around the converged rung -- the rungs the two-sided test
    # actually compared.  For the unit-step climb and the grid that is
    # the whole of `rungs`; for the bracket-refine climb it is the
    # filled bracket, dropping the sparse bracket endpoints below it
    # (search scaffolding, DESIGN 3.12.3).  Recording only the
    # consecutive block is what lets auto_promote_ok (15.7 / 7.8)
    # re-judge on adjacent meshes -- a sparse endpoint left in could
    # make the two-sided test read a false early convergence.
    trace = consecutive_block(rungs, rung, config.classes,
                              config.recipMag)
    return {
        converged_kpoint_density = converged_density,
        converged_mesh           = rung.mesh,       # 7.2 (Q4)
        # Ascending because `trace` is, and product(mesh) rises along
        # the ladder.
        grid_values   = [product(r.mesh) / config.recipCellVolume
                         for r in trace],
        grid_energies = [r.energy for r in trace],
    }


function consecutive_block(rungs, rung, classes, recipMag):
    # The maximal run of `rungs` whose meshes are consecutive ladder
    # positions (each one climbOneRung from the last) and that contains
    # `rung` (DESIGN 3.12.3).  Walks down then up from `rung` while the
    # immediate neighbour is present; climbOneRung(descendOneRung(m)) ==
    # m makes descendOneRung the exact one-below step (4e.1).
    by_mesh = { r.mesh : r for r in rungs }
    block   = [rung]
    m = rung.mesh                               # walk down
    loop:
        below = descendOneRung(m, classes, recipMag)
        if below == m or below not in by_mesh:
            break
        prepend by_mesh[below] to block
        m = below
    m = rung.mesh                               # walk up
    loop:
        above = climbOneRung(m, classes, recipMag)
        if above not in by_mesh:
            break
        append by_mesh[above] to block
        m = above
    return block
```

`build_entry` (15.7) gains the `converged_mesh` field on the
verification block; it reads the exact mesh from the chosen rung's
`result.toml` -- the same source in both harvests, since both hand
`build_entry` that rung's result -- and the emitter (15.4) writes
it when present.  `record_converged` therefore need not surface a
mesh into the entry; it supplies only the density and the flatness
ladder, and build_entry recovers the exact mesh from the result.
Everything else in the guidance schema and predictor is unchanged
(DESIGN 3.12.6).

### 4e.7 Dispatching a mesh (DESIGN 7.7)

The climb searches in mesh space, so a rung is dispatched as an
EXPLICIT mesh, not a density.  Three small pieces bridge a mesh to
a run and back: the calc-tag encoding, the predict-only builder that
seeds the climb, and the climb dispatcher `converge_by_climb` (4e.5)
drives.  All three are specified in DESIGN 7.7.

```
function encodeMeshValue(mesh):        # DESIGN 6.2.4 / 7.7
    # A mesh's calc-tag value token: the three axial counts joined
    # by hyphens.  Counts are positive integers, so the token stays
    # slug-safe ([a-z0-9-]).
    return str(mesh[0]) + "-" + str(mesh[1]) + "-" + str(mesh[2])

function decodeMeshValue(token):
    # Invert encodeMeshValue: "4-4-4" -> [4, 4, 4].
    return [int(part) for part in split(token, "-")]


function build_mesh_unit(structure, options, mesh, id, record):
    # One explicit-mesh convergence unit (DESIGN 7.7 / 6.2.1).
    # `scfkp` is the makeinput key for an explicit axial-count mesh
    # (a style-code-1 k-point file); `kpt-mesh` is its calc-tag axis.
    # The cache identity is the same one the density units used
    # (6.2.1), so a mesh re-run in a later round is a cache hit and
    # costs nothing.
    #
    # `record` is where the build identity travels (DESIGN 6.2.4 /
    # 6.2.10): a fact ABOUT the run, not an input TO it, so it must
    # not ride in `options` -- every key there is a real tool input,
    # which is what keeps makeinput's strict unknown-key check a
    # pure typo backstop.  The driver stamps it into status.toml at
    # launch (13.5) and the wingbeat echoes it into result.toml
    # (13.2).  Without it set HERE the whole recorded-not-compared
    # path is empty: the reuse plan would name no build and a
    # guidance entry's provenance would read "unknown" forever.
    unit_options = copy(options)
    unit_options["scfkp"] = mesh                 # [a, b, c]
    calc = buildCalcTag({ "kpt-mesh": encodeMeshValue(mesh) })
    return CalcUnit(id=id, calc=calc, structure=structure,
                    options=unit_options, wingbeat="imago",
                    record=record,
                    key_fields=standardKeyFields(structure, options))


function predict_kpoint_density(structure, dataspace, system_type,
                                submodel, center,
                                harvest_thresholds):
    # The prediction HALF of the former grid builder (7.7 steps 1-2
    # and 5): signature -> predict -> PredictionRecord.  It lays NO
    # grid -- the climb seeds from the density and picks its mode and
    # persistence from the confidence (3.12.4 / 3.12.6).  A curator-
    # pinned `center` bypasses the predictor exactly as before
    # (5.7 / 6.2.9).  Returns the seed density, the confidence and
    # under-trained flag the policy reads, and the record the harvest
    # recovers (7.8).
    sig = computeSignature(load(structure), system_type,
                           dataspace.group_table)
    if center is not None:
        result = None                            # curator override
        density = center
        confidence = 1.0
        under_trained = false
        policy = "curator_override"
    else:
        result = predict(dataspace, sig, submodel.basis,
                         submodel.functional,
                         submodel.kpoint_integration)
        density = result.predicted_kpoint_density
        confidence = result.confidence
        under_trained = result.is_under_trained
        policy = "predict_then_climb"
    # Named `prediction`, not `record`: `record` is the unit's
    #   build-identity bookkeeping in the neighbouring builders
    #   (4e.7 / DESIGN 6.2.4), and the two are unrelated.
    #   `harvest_thresholds` carries the two resolved manifest knobs
    #   the harvest cannot look up for itself -- the grid-flatness
    #   tolerance and the metal gap cut (15.6) -- so they are stamped
    #   here, on the one per-structure record the harvest recovers.
    prediction = buildPredictionRecord(policy, density, confidence,
                                       under_trained, result, sig,
                                       system_type, submodel,
                                       harvest_thresholds)
    return density, confidence, under_trained, prediction


function make_climb_dispatcher(structures, options_by_material,
                               workspace, parsl_config, executor,
                               force, record, tidy_run = False,
                               scratch_root = "",
                               prune_problems = None):
    # Build the dispatcher converge_by_climb (4e.5) drives, closing
    # over each material's structure and options, the workspace root,
    # the resolved Config (13.7), and the ONE shared executor every
    # send runs under (make_executor, 13.5): the pre-flight loen batch
    # and every climb rung run beneath the SAME executor, so the whole
    # run rides one warm pool (DESIGN 6.2.11) and lands in one tree.
    # `force` bypasses the run-reuse cache exactly as the pre-flight
    # dispatch does.  The material key doubles as the unit id
    # (materials ARE the reference ids the producer already uses).
    # `record` is the producer's per-run bookkeeping, in practice
    # `{imago_commit: <sha>}`, handed to every unit this dispatcher
    # builds; ONE value for the whole build, since it describes the
    # producer's run and not any one material (DESIGN 6.2.4).
    #
    # ONE flight spans the whole climb -- its root is the workspace --
    #   and its unit list ACCRETES as rungs are decided.  Each send
    #   re-serializes the flight, so flight.toml records every rung
    #   asked for (13.1), and any mesh already run is a cache hit
    #   (6.2.5).  The dispatcher owns the in-flight set and a small
    #   map from each launched unit to its (material, mesh), so an
    #   energy routes back without re-decoding the calc tag.
    flight = Flight(units = [], root = workspace,
                    parsl_config = parsl_config,
                    sweep = Sweep(varied = ("kpt-mesh",), fixed = {}))
    outstanding = []          # (unit, future) -- the in-flight rungs
    origin = {}               # unit identity -> (material, mesh)

    # Layer (b) (DESIGN 6.2.12).  With --tidy-run the flight carries
    #   the producer's prune callback (11.4), so each rung's scratch
    #   goes as that rung lands rather than after the climb ends.
    #   This is the flight that matters for it: the climb dispatches
    #   many rungs per solid, and is the phase long enough to fill
    #   scratch before it finishes.  The callback reads flight.units
    #   at call time, which is why the accreting list above is safe
    #   to hand it now.
    if tidy_run:
        flight.on_outcome = make_prune_callback(
            flight, workspace, scratch_root, prune_problems)

    # send: build one unit per requested mesh, remember its origin,
    #   append it to the flight's growing list, prepare ONLY the new
    #   units (6.2.5), and launch just them (send_off re-serializes
    #   the whole flight, 13.5).
    function send(mesh_lists):
        new_units = []
        for m in mesh_lists:
            for mesh in mesh_lists[m]:
                unit = build_mesh_unit(structures[m],
                                       options_by_material[m], mesh,
                                       id = m, record = record)
                origin[identity(unit)] = (m, mesh)
                append(new_units, unit)
                append(flight.units, unit)
        prepareUnits(flight, new_units)          # driver-side (6.2.5)
        launched = send_off(flight, new_units,   # phase 1, 13.5
                            executor, force)
        outstanding.extend(launched)

    # next_rung: block until the next rung lands (collect_next, 13.5),
    #   translate it to (material, Rung-or-FAILED), and return.  A
    #   unit that did not complete is a run failure the climb stops on
    #   (7.7); a landed mesh must equal the one requested (7.7), or
    #   makeinput/imago silently changed it -- fail loudly.
    function next_rung():
        (unit, entry, remaining) = collect_next(flight, outstanding)
        outstanding = remaining
        (m, mesh) = origin[identity(unit)]
        # TWO different questions, answered in two places.  The
        #   flight entry says whether the JOB completed; res.converged
        #   says whether the SCF reached its own fixed point.  A
        #   NOT_CONVERGED run does both -- exits cleanly AND writes a
        #   total energy -- so it passes the first test and its energy
        #   is indistinguishable from a real one once it is a number
        #   on a ladder (DESIGN 5.7).
        if entry.status != "done":
            return (m, FAILED)
        res = readResult(unit)                   # result.toml (6.1.2)
        assert res.kpoint_mesh == mesh           # honoured exactly
        # An unconverged energy is wherever the iteration happened to
        #   stop, so a flatness test over it asks the wrong question:
        #   it wants energy that has stopped moving with the MESH.
        #   Such a rung is treated as one that did not run, which
        #   stops the material.  Dropping it and climbing on was
        #   rejected -- the next mesh is chosen FROM the ladder, so a
        #   ladder that does not grow re-requests the same mesh
        #   forever.  Only an EXPLICIT false drops it: a result.toml
        #   with no such field cannot be judged and is kept, the same
        #   side taken on a missing gap (4e.2).
        if res.converged is false:
            return (m, FAILED)
        return (m, Rung(mesh, res.total_energy))

    return dispatcher(send = send, next_rung = next_rung)
```

The producer builds each rung's mesh with `build_mesh_unit`
(4e.7) as the climb chooses it, so the convergence search needs
no single up-front flight of grid points: `predict_kpoint_density`
(4e.7) supplies the seed, and the dispatcher above launches and
collects the rungs one at a time.

---

## 5. Save fullKPToIBZOpMap (DESIGN 2.4)

This augments the existing IBZ folding loop in
`initializeKPointMesh`.  The current code saves
`fullKPToIBZKPMap(k_full) = k_IBZ` (which IBZ point
does this full-mesh point fold onto).  We additionally
save `fullKPToIBZOpMap(k_full) = R` (which point group
operation maps the IBZ representative to k_full, in
the forward direction: R(k_IBZ) = k_full).

The change is a single integer store at the point where
a match is found, plus identity for the IBZ
representative itself.

```
# Inside initializeKPointMesh, after allocating
# fullKPToIBZOpMap(numFullMeshKP):

# When a new IBZ representative i is found:
    fullKPToIBZOpMap(i) = identityOpIndex

# When mesh point j matches IBZ point i under
# operation m (the existing isMatch==1 branch):
    fullKPToIBZOpMap(j) = m
```

**Identity operation index.**  The space group
database guarantees that the identity is always the
first point group operation (verified across all 759
space group files in share/spaceDB).  Therefore
`identityOpIndex = 1` -- no runtime search is needed.

---

## 6. Corrected Effective Charge (DESIGN 2.4)

The current code accumulates Q* directly into
`atomCharge(A)` using only the IBZ k-point's Mulliken
projection.  The fix loops over the star of each IBZ
k-point and distributes the projection into the
permuted atom index.

The star of IBZ k-point k_IBZ is the set of full-mesh
k-points that fold onto it.  This set is not stored
explicitly -- it is traversed by scanning
fullKPToIBZKPMap.

The key change is in the innermost accumulation.
The outer loop structure (spin h, kpoint i, band j,
basis function l, atom k) remains the same.  The
Mulliken projection `oneValeRealAccum` for atom k at
IBZ k-point i is computed exactly as before.  The
difference is where it is accumulated.

```
# Current code (incorrect with IBZ):
#   atomCharge(k, h) += oneValeRealAccum
#                       * statePopulation

# Corrected code:
#
# After the band loop (j) completes for kpoint i,
# we have accumulated per-atom projections for
# this kpoint:
#   ibzAtomProj(k) = sum over bands of
#       oneValeRealAccum(k,j) * statePopulation(j)
#
# Count the star size (number of full-mesh points
# that fold to this IBZ kpoint):
#   starSize = count(fullKPToIBZKPMap(:) == i)
#
# Distribute across the star:
#   for each full-mesh kpoint f where
#           fullKPToIBZKPMap(f) == i:
#       R = fullKPToIBZOpMap(f)
#       for A = 1 to numAtomSites:
#           atomCharge(atomPerm(R, A), h) +=
#               ibzAtomProj(A) / starSize

# The normalizer is starSize, not numFullMeshKP.
# statePopulation already encodes the full BZ-
# integration weight for the star of kpoint i
# (kPointWeight for Gaussian, or the summed
# tetrahedron corner weights for LAT -- both
# proportional to starSize).  ibzAtomProj is
# therefore the total contribution from kpoint i.
# Dividing by starSize distributes this total
# equally among the star members.  The sum across
# the star recovers ibzAtomProj, preserving the
# total charge.
```

**Alternative (equivalent, avoids scanning):**  The
same result can be achieved within the existing
kpoint loop by changing only the accumulation target.
Instead of accumulating at index k, accumulate at
atomPerm(R, k) for each operation R in the star.
But the star is implicit -- it requires collecting
all full-mesh k-points for IBZ index i.

The cleanest implementation collects the per-atom
projection for one IBZ k-point, then distributes it
in a separate inner loop over the star.

---

## 7. Corrected Bond Order (DESIGN 2.4)

The same star-distribution pattern applies to bond
order.  The Mulliken overlap between atoms A and B
at IBZ k-point i is computed as before.  The
correction distributes it into the permuted atom
pair.

```
# Current code (incorrect with IBZ):
#   bondOrder(A_bonded, k) +=
#       oneValeRealAccum * statePopulation

# Corrected code:
#
# After computing bondOrderRaw(A, B) for IBZ
# kpoint i (accumulated over bands j):
#
#   starSize = count(fullKPToIBZKPMap(:) == i)
#
#   for each full-mesh kpoint f where
#           fullKPToIBZKPMap(f) == i:
#       R = fullKPToIBZOpMap(f)
#       for each bonded pair (A, B):
#           A_rot = atomPerm(R, A)
#           B_rot = atomPerm(R, B)
#           bondOrder(A_rot, B_rot) +=
#               bondOrderRaw(A, B) / starSize
```

**Integration with existing loop structure.**  The
current `computeBond` accumulates bond order inside
the band loop (j) interleaved with charge
accumulation.  The IBZ correction requires a second
pass over the star after all bands are processed for
a given IBZ k-point.  This suggests restructuring:

  1. For each IBZ kpoint i:
     a. Read eigenvectors and overlap (unchanged)
     b. Loop over bands j, accumulate raw per-atom
        charge ibzAtomProj(A) and raw per-pair bond
        order ibzBondRaw(A, B) at IBZ indices
     c. Distribute ibzAtomProj and ibzBondRaw across
        the star of i using atomPerm

Step (c) is the only new code.  Steps (a) and (b)
are the existing computation, with the accumulation
target changed from the final arrays to temporary
per-IBZ-kpoint buffers.

**Weight convention.**  The raw projection is
weighted by statePopulation (from either
electronPopulation_LAT or electronPopulation).
Both already encode the full star-weighted BZ-
integration weight for each IBZ kpoint:

- Gaussian: statePopulation includes
  kPointWeight(i) = 2 * starSize(i) /
  numFullMeshKP, so ibzAtomProj and ibzBondRaw
  are proportional to starSize.
- LAT: statePopulation includes
  electronPopulation_LAT(j,i,h), which sums
  tetrahedron corner weights from all full-mesh
  points in the star, again proportional to
  starSize.

The star distribution divides by starSize to
extract the per-full-mesh-point contribution, then
deposits it at the permuted atom (or atom pair).
The sum across starSize members recovers the
original ibzAtomProj (or ibzBondRaw), so the
overall charge and bond order totals are preserved.

---

## 7a. POPTC IBZ Unfolding (DESIGN 2.5)

The partial optical properties decompose the momentum matrix
element of a single transition between a *pair* of groups,
so what each k-point carries is a matrix over group indices
rather than a vector over them.  Section 7 already
distributes a two-index quantity across the star; this
section applies the same distribution to a quantity that
additionally carries a Cartesian component and a
transition-pair index.

**This section is the GAUSSIAN pathway's correction.** It
exists because that pathway visits only irreducible
k-points and must spread each one's contribution over the
members of its star.  Section 19 adds a tetrahedron pathway
that visits full-mesh corners directly and permutes once per
corner as it fetches, so this block does not run there at
all -- doing both would count the symmetry twice (section
19.5).  Guard it on `kPointIntgCode == 0` when that pathway
lands.  Both pathways must be correct under reduction; they
reach it by different routes.

### Which detail codes this touches

```
code  grouping       resolution   action here
-------------------------------------------------------
1     type           total        nothing
2     type           nl           nothing
3     atom           total        star average
4     atom           nl           star average
```

DESIGN 11.3 numbers the codes with grouping as the major
key, so the test is a threshold rather than a set of cases:
**codes 1 and 2 need nothing, codes 3 and 4 need the star
average, and the boundary is `detailCodePOPTC >= 3`.**

The type-grouped codes are already correct on a reduced mesh
because every operation Imago reduces by carries each atom
onto an atom of the same type -- which `buildAtomPerm`
enforces at startup rather than infers from a type being a
symmetry orbit (DESIGN 2.3) -- so a type-level sum maps onto
itself.  The resolution axis does not enter for either
grouping: both offered resolutions sum over complete shells,
which is the condition equation (4) of DESIGN 2.3 needs.

**The distinction matters because types are often not
orbits.**  In an amorphous cell a type is a bin of locally
similar environments with no symmetry content; in a point
defect supercell types come from the *pre-defect* symmetry.
Both land safely here, but by verification rather than by
assumption: the amorphous cell carries only the identity to
reduce by, and the defect supercell is typed *coarser* than
the true orbits, which is the safe direction because a union
of orbits stays closed.  Typing *finer* than the orbits is
the unsafe direction, and that is the case `buildAtomPerm`
aborts on.  DESIGN 2.3 has the argument; DESIGN 11.6 has
what it means for reading the output, which is a separate
question from whether the arithmetic is sound.

### The permuted partial

For code 3 the partial index *is* the atom site index, and
`atomPerm` could index the pair matrix directly.  For code 4
it is not: section 18.2 lays a partial out as a slot within
a segment, so a partial carries a site and a QN_nl slot
together.  The correction therefore permutes PARTIALS, not
atoms, and the atom permutation is what the partial
permutation is built from:

```
# The image of each partial under operation R.  Built once
#   per run, not per k-point -- it depends only on the
#   permutation table and the layout of section 18.2.
for R = 1, numPointOps:
   for site = 1, numAtomSites:
      siteRot = atomPerm(R, site)
      for slot = 1, slotsPerSegment(site):
         partialPerm(R, segmentBase(site) + slot) =
               segmentBase(siteRot) + slot
```

**The slot survives the permutation unchanged**, and that is
the whole reason this works.  `buildAtomPerm` only ever maps
an atom onto an atom of the same type, and atoms of one type
share a basis, so the permuted site's slots stand in
one-to-one correspondence with the original's and carry the
same QN_nl meaning.  This is the same argument `computeBond`
already relies on for orbital-resolved charge.

For code 3 there is one slot per segment and `partialPerm`
reduces to `atomPerm` exactly, so the block below is written
once for both codes rather than special-cased.

### The correction

```
# In computePOPTCPairs, after the transition double loop has
#   filled transitionProbTemp(:,:,c,1..transPairCount) for
#   IBZ k-point i, and BEFORE the mergeSort copy into
#   transitionProbPOPTC.  Atom-grouped codes only.

if (detailCodePOPTC < 3)        skip
if (.not. allocated(atomPerm))  skip, after the style code 0
                                warning described below

# The star of this IBZ k-point, counted exactly as in
#   section 7.
starSize = count(fullKPToIBZKPMap(:) == i)

allocate pairSlabSym(sumNumPartials, sumNumPartials)

for each transition pair p = 1, transPairCount:
   for each Cartesian component c = initComponent,
                                   finComponent:

      pairSlabSym(:,:) = 0

      for each full-mesh kpoint f with
              fullKPToIBZKPMap(f) == i:
         R = fullKPToIBZOpMap(f)
         for a = 1, sumNumPartials:
            aRot = partialPerm(R,a)
            for b = 1, sumNumPartials:
               bRot = partialPerm(R,b)
               pairSlabSym(aRot,bRot) =
                     pairSlabSym(aRot,bRot)
                     + transitionProbTemp(a,b,c,p)
                       / starSize

      transitionProbTemp(:,:,c,p) = pairSlabSym(:,:)

deallocate pairSlabSym
```

The deposit-forward shape (permute the source index and add
into the destination) is written to match section 7 rather
than the gather form `M(a,b) <- M(invAtomPerm(R,a),
invAtomPerm(R,b))` that DESIGN 2.5 states.  The two are the
same map read in opposite directions; one table suffices,
and section 7's is the one already in the codebase.

### Why the permutation alone is exact, and for what

The momentum operator is a vector, so an operation does two
things at once: it relabels the atoms, and it mixes the
Cartesian components, P_c(Rk) = sum_d R_cd P_d(k).  The
permutation above handles the first and not the second.
What makes it exact anyway -- for the isotropic column, and
only for that column -- is that the mixing cancels when the
three components are summed.

Write M^c_ab for the decomposed matrix element of one
transition at the IBZ point and M^c_tot for its sum over all
(a,b).  The stored quantity is

    T^c_ab = Re(M^c_ab) Re(M^c_tot)
           + Im(M^c_ab) Im(M^c_tot)

whose sum over (a,b) is |M^c_tot|^2, which is the "sum
squared to sum of squares" construction in the source.  At a
star member Rk both factors pick up the mixing, so

    T^c_ab(Rk) = sum_{d,e} R_cd R_ce
                 [ Re(M^d_a'b') Re(M^e_tot)
                 + Im(M^d_a'b') Im(M^e_tot) ]

with a' = invAtomPerm(R,a) and b' = invAtomPerm(R,b).
Summing over c and using the orthogonality of R,
sum_c R_cd R_ce = delta_de, kills every cross term:

    sum_c T^c_ab(Rk) = sum_d T^d_a'b'(k)

So the component-summed pair matrix transforms by pure index
permutation, exactly as bond order does.  Permuting each
component slab separately -- which is what the block above
does -- yields precisely this sum, because a relabeling of
(a,b) commutes with summing over c.

The per-component slabs it leaves behind are not the correct
per-component slabs; they are wrong in the same way and to
the same degree that the per-axis columns of the TOTAL
spectra are already wrong on a reduced mesh, which is the
open question recorded as TODO O3.  This change neither
repairs that nor worsens it.

Two consequences worth stating plainly:

- The isotropic column that `printSpectrumPOPTC` writes (the
  three components summed and divided by three) becomes
  correct per atom pair on a symmetry-reduced mesh.
- The per-axis POPTC columns remain unverified, exactly as
  before, until O3 is settled.

### Why the sum rule cannot check this

`transitionProb`, the total that `getOptcCond` broadens, is
formed as the sum of the pair matrix over (a,b).  A
permutation of (a,b) does not change that sum, and neither
does averaging permuted copies of it.  The total spectra are
therefore bit-for-bit unchanged by this correction, and the
identity that the partials sum to the total holds both
before and after.  It is exactly blind to whether the
unfolding is present, absent, or wrong.

Verify instead by running one structure whose
symmetry-equivalent atoms are inequivalently oriented, once
on a full mesh (`applySymmetry` off) and once IBZ-reduced,
and requiring the two isotropic per-atom-pair spectra to
agree.

Compare the broadened spectra rather than the probability of
any single transition pair.  Within a degenerate multiplet
the eigenvectors are fixed only up to a rotation among the
degenerate bands, so the decomposition of one band taken on
its own is basis dependent, and the diagonalizer at a full
mesh point need not choose the same basis as at the IBZ
representative.  The sum over the multiplet is invariant,
and since its members land at the same energy under the same
broadening, the spectra agree even where individual pairs do
not.  This is a property of every per-band Mulliken
decomposition in the codebase, not something introduced
here, but it is the failure that a per-pair comparison would
report as a bug.

### Placement, and what O3 will do to it

The star sum is a fixed linear map on the pair matrix.  It
does not depend on energy, and Gaussian broadening is
linear, so it is applied once per transition per k-point,
before broadening, and the weighted accumulation over IBZ
points in `getOptcCondPOPTC` is left exactly as it stands.
Placing the star loop inside the energy loop instead would
multiply the innermost work -- a sumNumPartials by
sumNumPartials by three slab per transition per energy point
-- by the reduction factor of 4 to 48 for the same answer.

The cost as written is transPairCount * 3 * starSize *
sumNumPartials^2 additions per IBZ k-point, against a
transition loop that is already valeDim * sumNumPartials^2
per pair.  The scratch slab is one sumNumPartials by
sumNumPartials double array.  Code 4 makes that partial
count the larger of the two atom-grouped cells, so this is
the cell where the star loop's cost is felt; DESIGN 11.4
gives the sizes.

When O3 is taken up, this block moves.  Rotating the
Cartesian components correctly cannot be done on T at all:
the R_cd R_ce cross terms above do not factor through T^d,
so the star average has to be formed from the complex
M^c_ab, one star member at a time, with T built afterwards
from the rotated copy.  That means lifting the star loop up
into the transition double loop and splitting the component
loop, since all three components of M must exist before any
one of them can be rotated.  It is a restructuring rather
than an insertion, and it is deliberately not attempted
here: this block is self-contained and switchable, which is
what makes the full-mesh comparison above interpretable.

### Guards

Detail codes 1 and 2 skip the block entirely, because a
type-grouped sum is already invariant.  Nothing is deferred:
every offered cell is either correct as computed or made
correct here, since DESIGN 11.2 withdraws the one resolution
that would have needed D^l(R).

Style code 0 supplies an explicit k-point list, from which
Imago cannot build the symmetry maps, so neither this
unfolding nor the `buildAtomPerm` closure check runs.  The
guard is `allocated(atomPerm)`, and the case is covered by
the standing warning of DESIGN 2.6.  A hand-supplied,
already-reduced k-point list is the one configuration where
a wrong atom-level decomposition can pass silently, so the
warning is the only protection available and must not be
dropped.

An unreduced mesh reaches the block with starSize = 1 and
the identity operation, so it costs one extra copy of the
slab and changes nothing.

---

## 8. LAT PDOS (DESIGN 1.4)

The LAT PDOS requires Mulliken projections at all four
corners of each tetrahedron simultaneously. Since
eigenvectors exist only at IBZ k-points, a two-pass
design is required: first compute and store projections
at IBZ k-points, then integrate over tetrahedra with
on-the-fly IBZ unfolding of the channel index.

### 8.1 Channel Permutation Table

For efficiency the channel permutation is precomputed
as a lookup table channelPermTable(R, alpha) so the
inner loop avoids repeated decode/encode. Mode 0
needs no permutation. Mode 1 uses invAtomPerm
directly. Mode 2 remaps the atom index while
preserving the l-shell offset within the atom.

```
function buildChannelPermTable(
        detailCodePDOS, numPointOps,
        cumulDOSTotal, cumulNumDOS,
        numAtomSites, invAtomPerm):

    allocate channelPermTable(numPointOps,
                              cumulDOSTotal)

    if detailCodePDOS == 0:
        # Per-type, per-l: identity (type-level
        # sums are invariant under R).
        for R = 1 to numPointOps:
            for alpha = 1 to cumulDOSTotal:
                channelPermTable(R, alpha) = alpha
        return channelPermTable

    if detailCodePDOS == 1:
        # Per-atom total: channel = atom index.
        for R = 1 to numPointOps:
            for A = 1 to numAtomSites:
                channelPermTable(R, A) =
                    invAtomPerm(R, A)
        return channelPermTable

    if detailCodePDOS == 2:
        # Per-atom, per-l: remap atom index,
        # preserve l-shell offset.
        for R = 1 to numPointOps:
            for A = 1 to numAtomSites:
                permA = invAtomPerm(R, A)
                baseOld = cumulNumDOS(A)
                baseNew = cumulNumDOS(permA)
                nOrbitals = cumulNumDOS(A+1)
                          - cumulNumDOS(A)
                for off = 1 to nOrbitals:
                    channelPermTable(R,
                        baseOld + off) =
                        baseNew + off
        return channelPermTable
```

### 8.2 Pass 1: Compute Projections

Stream through IBZ k-points, read eigenvectors and
overlap from HDF5, compute Mulliken projections, and
store into projArray(channel, band, kIBZ). The
Mulliken computation is identical to the existing code
in computeDOS (waveFnSqrd, oneValeRealAccum).

```
function computeProjections(inSCF, h,
        numKPoints, numStates, numAtomSites,
        numAtomStates, pdosIndex, valeDim,
        cumulDOSTotal, spin):
    allocate projArray(cumulDOSTotal,
                       numStates, numKPoints)
    projArray = 0.0

    for i = 1 to numKPoints:
        # Read eigenvectors + overlap for this
        # IBZ kpoint and spin orientation.
        readData(h, i, numStates, overlapCode=1)

        for j = 1 to numStates:
            valeDimIndex = 0
            for k = 1 to numAtomSites:
                for l = 1 to numAtomStates(k):
                    valeDimIndex += 1

                    # Compute Mulliken projection
                    # (existing waveFnSqrd *
                    # valeValeOL dot product).
                    oneValeRealAccum =
                        mullikenProjection(
                            valeDimIndex, j)

                    # Accumulate into the channel
                    # determined by pdosIndex.
                    ch = pdosIndex(valeDimIndex)
                    projArray(ch, j, i) +=
                        oneValeRealAccum
                        / real(spin)

    return projArray
```

### 8.3 Pass 2: Tetrahedron Integration

Loop over bands and tetrahedra. For each tetrahedron,
sort corner eigenvalues with tracked permutation,
compute `bloechlCornerDOSWt` at each energy point,
and accumulate weighted projections into pdosComplete.
The channel permutation table handles IBZ unfolding of
the projection index.

Note: this uses `bloechlCornerDOSWt` (section 2a),
which returns per-corner DOS density weights (units:
1/energy). The cumulative corner weights from
`bloechlCornerWeights` (section 3a) are NOT used here
-- those are for integrated properties only (section
3, `electronPopulation_LAT`).

```
function integratePDOS_LAT(projArray,
        channelPermTable,
        eigenValues, tetrahedra,
        numTetrahedra, tetraVol,
        fullKPToIBZKPMap, fullKPToIBZOpMap,
        energyScale, numEnergyPoints,
        numStates, cumulDOSTotal, spin):
    allocate pdosComplete(cumulDOSTotal,
                          numEnergyPoints)
    pdosComplete = 0.0

    for n = 1 to numStates:
        for T = 1 to numTetrahedra:
            # Look up corner info from the full
            # mesh, mapping to IBZ eigenvalues.
            for c = 1 to 4:
                kFull(c) = tetrahedra(c, T)
                kIBZ(c) =
                    fullKPToIBZKPMap(kFull(c))
                opIdx(c) =
                    fullKPToIBZOpMap(kFull(c))
                eps(c) =
                    eigenValues(n, kIBZ(c), h)

            # Sort eigenvalues ascending, tracking
            # permutation: sigma(i) = original
            # corner index in sorted position i.
            sigma = argsort(eps)
            sortedEps = eps(sigma)

            for iE = 1 to numEnergyPoints:
                E = energyScale(iE)

                # Skip if outside eigenvalue range.
                if E < sortedEps(1) or
                        E >= sortedEps(4):
                    cycle

                # Per-corner DOS density weights
                # for the sorted eigenvalues.
                cornerDOSWt_LAT(1:4) =
                    bloechlCornerDOSWt(
                        E, sortedEps)

                # Accumulate weighted projections
                # into pdosComplete. Each sorted
                # corner c maps back to original
                # corner sigma(c), whose IBZ kpoint
                # and operation index determine the
                # projection lookup.
                for c = 1 to 4:
                    orig = sigma(c)
                    R = opIdx(orig)
                    kIc = kIBZ(orig)

                    for alpha = 1 to cumulDOSTotal:
                        permA =
                            channelPermTable(
                                R, alpha)
                        pdosComplete(alpha, iE) +=
                            cornerDOSWt_LAT(c)
                            * tetraVol * kpWtSum
                            / hartree
                            * projArray(
                                permA, n, kIc)

    return pdosComplete
```

`kpWtSum` is `sum(kPointWeight)`, and it is not optional.
DESIGN 1.3 requires it of every LAT accumulation, for the
reason given there: `tetraVol` sums to 1 while
`kPointWeight` sums to 2, so without it the LAT result sits
at half the Gaussian one and the two paths cannot be
compared. It is easy to leave out, because unlike the
`/hartree` conversion it corrects no unit -- both paths are
dimensionally correct without it and simply disagree by a
factor of two.

### 8.4 Normalization

For the LAT path, the corner DOS weights from
`bloechlCornerDOSWt` provide exact BZ integration
(no broadening artifacts). The electronFactor ratio
(currentPopulation / totalElectronsComputed) should
be ≈ 1.0. Compute and log it as a diagnostic but do
not apply it to pdosComplete. A ratio significantly
different from 1.0 signals an integration bug.

The "Spin States Calculated" diagnostic (integral of
totalSystemDos over the energy grid) must use
`deltaDOS * hartree` in the trapezoidal rule because
deltaDOS is stored in Hartree while the DOS is in
states/eV. This applies to both the Gaussian and LAT
paths.

---

## 9. UFF Bond Parameter Computation (DESIGN 4.2)

Given two atomic numbers, compute the UFF equilibrium
bond length and harmonic force constant.  The per-element
parameters (covalent radius r_i, effective charge Zstar_i,
GMP electronegativity chi_i) are read once from
`bond_parameters.dat` and stored in arrays indexed by
atomic number Z.

The prefactor 332.06 = 664.12 / 2 converts from the UFF
spring constant convention E = (1/2) k (r-r0)^2 to the
LAMMPS `bond_style harmonic` convention E = K (r-r0)^2.

```
# -------------------------------------------------
# Data structures (populated once by init_bond_data
# from bond_parameters.dat):
#
#   num_uff_elements : int
#       Number of elements in the table
#       (= maximum Z covered).
#   uff_r(Z)     : covalent radius (Angstroms)
#   uff_Zstar(Z) : effective charge
#   uff_chi(Z)   : GMP electronegativity (eV)
#
# These arrays are indexed by atomic number Z
# (1-based: uff_r(1) = hydrogen, etc.).
#
# The reader uses the Z column on each data line
# as the array index (not the sequential row
# number).  This makes the file order-independent.
# -------------------------------------------------

UFF_K_PREFACTOR = 332.06

function get_bond_params(z1, z2):
    # Compute UFF equilibrium bond length and
    # LAMMPS harmonic force constant for the
    # element pair (z1, z2).
    #
    # Inputs:
    #   z1, z2 : atomic numbers (order irrelevant;
    #            the formula is symmetric)
    #
    # Returns:
    #   K_ij : force constant (kcal/mol/A^2)
    #   r_ij : equilibrium bond length (Angstroms)
    #
    # Requires uff_r, uff_Zstar, uff_chi arrays
    # to be initialized from bond_parameters.dat.

    # --- Validate element coverage ---
    if z1 < 1 or z1 > num_uff_elements:
        error("Element Z =", z1,
              "not in bond_parameters.dat")
    if z2 < 1 or z2 > num_uff_elements:
        error("Element Z =", z2,
              "not in bond_parameters.dat")

    # --- Look up per-element parameters ---
    r1    = uff_r(z1)
    r2    = uff_r(z2)
    Zs1   = uff_Zstar(z1)
    Zs2   = uff_Zstar(z2)
    chi1  = uff_chi(z1)
    chi2  = uff_chi(z2)

    # --- Electronegativity correction ---
    # r_EN shortens the bond between elements of
    # unequal electronegativity.  For homonuclear
    # bonds (chi1 == chi2), r_EN = 0 and the bond
    # length is simply r1 + r2.
    denom_EN = chi1 * r1 + chi2 * r2
    if denom_EN > 0:
        r_EN = r1 * r2
               * (sqrt(chi1) - sqrt(chi2))**2
               / denom_EN
    else:
        r_EN = 0.0

    # --- Equilibrium bond length ---
    r_ij = r1 + r2 - r_EN

    # --- Force constant ---
    # Guard against zero or near-zero bond length
    # (should not occur for physical elements, but
    # protects against corrupt data).
    if r_ij <= 0:
        error("Non-positive bond length for Z =",
              z1, z2, "; check bond_parameters.dat")

    K_ij = UFF_K_PREFACTOR * Zs1 * Zs2
           / r_ij**3

    return K_ij, r_ij
```

**Usage in create_lammps_files and normalize_types.**
Both output paths contain a linear scan over
`hooke_bond_coeffs` to match element pairs to force
constants.  Both are replaced by a direct call to
`get_bond_params`:

```
# Current code (linear scan, in both paths):
#   for hb = 1 to num_hooke_bonds:
#       if atom1_z == hbc[hb][1]
#              and atom2_z == hbc[hb][2]:
#           k  = hbc[hb][3]
#           r0 = hbc[hb][4]
#
# New code (direct computation, both paths):
    K_ij, r_ij = bond_data.get_bond_params(
                     atom1_z, atom2_z)
    K_ij = K_ij * self.bond_parameter_scale
```

The `bond_parameter_scale` multiplier (default 1.0,
defined in `condenserc.py`, overridable in condense.in)
is applied after the UFF computation in **both**
output paths.  It scales only K_ij, not r_ij.

---

## 10. Angle Clustering and Force Constants (DESIGN 4.8)

Replace the `angles.dat` database lookup with a two-phase
procedure: (a) cluster observed bond angles by element
triplet to discover angle types, then (b) compute a force
constant for each type from UFF bond stiffnesses.

### 10a. Cluster Observed Angles by Triplet (DESIGN 4.8.3)

For each molecule, the bond analysis produces a list of
angles with atom indices and observed angle values.  The
`create_lammps_files` method already iterates over these
and extracts the element triplet (Z1, Zv, Z2).  The new
clustering step replaces the `angles.dat` lookup.

```
# -------------------------------------------------
# Data structures:
#
#   angle_cluster_tolerance : float
#       Maximum deviation (degrees) between an
#       observed angle and a cluster's running mean
#       for the angle to join that cluster.
#       Default: 5.0.  Read from condense.in.
#
#   spread_cap : float
#       Maximum allowed total span (max - min) of
#       theta values within any one cluster.  Set
#       to 2.0 * angle_cluster_tolerance.  Prevents
#       a long chain of closely-spaced observations
#       from silently sweeping values from opposite
#       ends of a wide distribution into a single
#       cluster.  The same cap is applied in 10e
#       for cross-source clustering.
#
#   Input: a list of angle observations, each
#       being (Z1, Zv, Z2, theta_obs, base_tag)
#       where Z1 <= Z2 (canonicalized).  The
#       base_tag is the producer's tag prefix
#       (element names, species ids, molecule
#       ids) for this specific atom triple.
#
#   Output:
#       angle_types : list of
#           (Z1, Zv, Z2, theta_0, obs_count,
#            base_tag)
#       angle_type_map : maps each observation
#           index to its angle_type index
#
#       obs_count is the number of observations
#       merged into the cluster.  The
#       representative base_tag is taken from
#       the first observation in the cluster.
#       The slot ordering (obs_count at slot 5,
#       base_tag at slot 6) matches 10e's
#       local_records and final_types tuples,
#       so 10c/10d/10e/10f use consistent slot
#       indices throughout.
# -------------------------------------------------

function cluster_angles(observations, tolerance):
    # Group observations by element triplet.
    # Each entry carries the observed theta, the
    # original observation index (for the
    # angle_type_map), and the producer's
    # representative base_tag for that observation.
    groups = {}
    for each (idx, obs) in enumerate(observations):
        key = (obs.Z1, obs.Zv, obs.Z2)
        groups[key].append(
            (obs.theta, idx, obs.base_tag))

    angle_types = []
    angle_type_map = array(len(observations))
    spread_cap = 2.0 * tolerance

    for each key in groups:
        # Sort angles within this triplet group.
        entries = groups[key]
        sort entries by theta ascending

        # Greedy clustering: walk the sorted list
        # and merge the candidate into the current
        # cluster while BOTH of these hold:
        #   (a) |theta - running_mean| <= tolerance
        #   (b) resulting (max - min) <= spread_cap
        # If either fails, finalize the current
        # cluster as a type and start a new cluster
        # at the candidate.  cluster_rep_base_tag is
        # captured from the first observation in
        # the cluster and propagated to the emitted
        # angle_type record so 10c/10d can build
        # tag tails and 10e can carry it across
        # sources.
        cluster_rep_base_tag = entries[0].base_tag
        cluster_sum   = entries[0].theta
        cluster_count = 1
        cluster_min   = entries[0].theta
        cluster_max   = entries[0].theta
        cluster_members = [entries[0].idx]

        for i = 1 to len(entries) - 1:
            cluster_mean = cluster_sum / cluster_count
            candidate_theta = entries[i].theta
            new_max = max(cluster_max, candidate_theta)
            new_min = min(cluster_min, candidate_theta)

            within_tol =
                |candidate_theta - cluster_mean|
                    <= tolerance
            within_cap =
                (new_max - new_min) <= spread_cap

            if within_tol and within_cap:
                # Merge into current cluster.
                cluster_sum += candidate_theta
                cluster_count += 1
                cluster_min = new_min
                cluster_max = new_max
                cluster_members.append(entries[i].idx)
            else:
                # Finalize current cluster as a type.
                # cluster_count is the number of raw
                # observations that fed this local
                # cluster (slot 5 = obs_count); the
                # first member's base_tag is the
                # representative prefix (slot 6).
                theta_0 = cluster_sum / cluster_count
                type_id = len(angle_types) + 1
                angle_types.append(
                    (key.Z1, key.Zv, key.Z2, theta_0,
                     cluster_count,
                     cluster_rep_base_tag))
                for m in cluster_members:
                    angle_type_map[m] = type_id

                # Start a new cluster at entry i.
                cluster_rep_base_tag =
                    entries[i].base_tag
                cluster_sum = candidate_theta
                cluster_count = 1
                cluster_min = candidate_theta
                cluster_max = candidate_theta
                cluster_members = [entries[i].idx]

        # Finalize the last cluster.  Slot
        # ordering matches the finalize above:
        # slot 5 = obs_count, slot 6 = base_tag.
        theta_0 = cluster_sum / cluster_count
        type_id = len(angle_types) + 1
        angle_types.append(
            (key.Z1, key.Zv, key.Z2, theta_0,
             cluster_count,
             cluster_rep_base_tag))
        for m in cluster_members:
            angle_type_map[m] = type_id

    return angle_types, angle_type_map
```

### 10b. Angle Force Constant (DESIGN 4.8.4)

Compute the harmonic angular spring constant for a
given angle type from the UFF bond stiffnesses of its
two arms.  This reuses `get_bond_params()` (section 9).

```
# -------------------------------------------------
# Data structures:
#
#   angle_stiffness_coeff : float
#       Dimensionless calibration constant that
#       converts the geometric mean of bond
#       stiffnesses into an angular stiffness.
#       Default: 0.15.  Read from condense.in.
#
#   angle_parameter_scale : float
#       Global multiplier on all angle force
#       constants.  Default: 1.0.
#       Read from condense.in.
# -------------------------------------------------

function get_angle_k(z1, zv, z2,
                     angle_stiffness_coeff,
                     angle_parameter_scale):
    # Compute the bond force constants for the
    # two arms of the angle: (z1, zv) and (zv, z2).
    K_arm1, _ = get_bond_params(z1, zv)
    K_arm2, _ = get_bond_params(zv, z2)

    # Geometric mean of arm stiffnesses, scaled
    # by the calibration constant and the global
    # user scale factor.
    K_angle = angle_stiffness_coeff
              * sqrt(K_arm1 * K_arm2)
              * angle_parameter_scale

    return K_angle
```

### 10c. Integration into create_lammps_files (DESIGN 4.8.8)

The existing angle loop in `create_lammps_files` extracts
(Z1, Zv, Z2) triplets and observed angles, then searches
`angles.dat` for a match.  The replacement runs the same
collect-cluster-emit structure as 10d, scoped to the
single lammps.dat file produced by `create_lammps_files`.
Both producers (10c and 10d) invoke the identical
`cluster_angles` helper from 10a, so local clustering
semantics are byte-identical.  Any residual theta_0
differences between 10c and 10d outputs are resolved by
10e during cross-source clustering inside
`normalize_types`.

`create_lammps_files` writes one further file this section
does not cover: the SLURM submission file that runs the
condensation.  That is a job-class question rather than a
force-field one, so it lives with the other submission-file
generator, in 13.7 (`condense_write_submission`).

```
# Phase 1: Collect all angle observations.  Replaces
# the angles.dat lookup loop.  The base_tag is built
# exactly as the current code builds it -- element
# names, species ids, and molecule ids -- with no
# rest-angle or type-id suffix appended yet (Phase 3
# adds those).
observations = []
for each atom with bond angles:
    for each angle_idx of atom:
        # Extract end atoms a1, a2 and vertex atom.
        z1 = element_z(a1)
        zv = element_z(atom)
        z2 = element_z(a2)
        if z1 > z2:
            swap z1, z2
            swap a1, a2
        theta_obs = bond_angles_ext[atom][angle_idx]
        base_tag = tag_string_for(a1, atom, a2)
        observations.append(
            (z1, zv, z2, theta_obs, base_tag,
             a1, atom, a2))

# Phase 2: Cluster locally using the shared helper
# from 10a.  Identical call signature and tolerance
# value as 10d.
angle_types, angle_type_map =
    cluster_angles(observations,
                   self.angle_cluster_tolerance)

# Phase 3: Build local type records with the
# cluster-mean theta_0 carried in the tag tail.
# These tables are local to this lammps.dat;
# normalize_types may merge them with types from
# reaction templates via 10e and rewrite both the
# tags and the per-angle ids in 10f.
num_local_angle_types = len(angle_types)
local_angle_tags   =
    [None] * (num_local_angle_types + 1)
local_angle_coeffs =
    [None] * (num_local_angle_types + 1)

for t = 1 to num_local_angle_types:
    atype = angle_types[t - 1]
    K = get_angle_k(
        atype.Z1, atype.Zv, atype.Z2,
        self.angle_stiffness_coeff,
        self.angle_parameter_scale)
    local_angle_coeffs[t] =
        [None, K, atype.theta_0]
    local_angle_tags[t] = (
        f"{atype.base_tag} "
        f"{atype.theta_0:.4f} {t}")

# Phase 4: Record per-atom angle entries with the
# local type ids.  normalize_types walks these and
# remaps the ids in 10f.  angle_bonded_atoms and
# ordered_angle_type follow the existing flat
# per-angle layout that the LAMMPS writer expects.
for i, obs in enumerate(observations):
    local_type_id = angle_type_map[i]
    angle_bonded_atoms.append(
        [None, obs.a1, obs.atom, obs.a2])
    ordered_angle_type.append(local_type_id)

# Export to normalize_types:
#   source tag = "lammps.dat"
#   local_angle_tags, local_angle_coeffs
#   angle_bonded_atoms, ordered_angle_type
#   per-local-type obs_count (slot 5 of each
#       entry in angle_types)
```

### 10d. Integration into make_reactions.py (DESIGN 4.8.8)

Mirrors 10c for the template-emission side of the
pipeline, but with two deliberate differences that follow
from DESIGN 4.8.8 item 3 and DESIGN 4.8.10:

1. **Local clustering tolerance is fixed at 0.**  The
   `cluster_angles` call here uses tolerance = 0.0 so
   only bit-identical observations (after the 0.5-degree
   rounding that `_read_angle_data` already applies to
   entries from `bondAnalysis.ba`) are collapsed into the
   same local type.  This preserves the template
   reusability property (DESIGN 4.8.10): any downstream
   `condense.py` run can apply any
   `angle_cluster_tolerance` value to the records
   emitted here, because no non-identical observations
   have been fused at the producer.  The obs_count weighting in 10e is
   associative under this pre-merge (DESIGN 4.8.8 item
   4a), so collapsing identical duplicates does not
   change the final cross-source theta_0.

2. **No K_angle computation.**  Reaction template files
   carry only connectivity, per-atom angle entries, and
   the tag tail "{theta_0_local} {t}" -- no K value is
   ever written to a template.  `normalize_types()`
   recomputes K authoritatively from the triplet in 10e
   / DESIGN 4.8.8 item 4b, which does not depend on any
   producer-side K.  `make_reactions.py` therefore does
   not call `get_angle_k` and does not need `BondData`.

The existing Python port's angle construction loop
(around line 2518) iterates over angles in a reaction
template and searches `hooke_angle_coeffs` for a matching
row to build the tag tail.  The replacement runs the
collect-cluster-emit structure below, scoped to one
reaction template at a time.

```
# Phase 1: Collect all angle observations for this
# reaction template.  Replaces the hooke_angle_coeffs
# scan near line 2518.  The base_tag is built exactly
# as today -- element names, species ids, and
# molecule ids, with no rest-angle or type-id suffix
# appended yet (Phase 3 adds those).
observations = []
for each vertex atom v in the template:
    for each (a1, a2) angle arm pair through v:
        z1 = element_z(a1)
        zv = element_z(v)
        z2 = element_z(a2)
        if z1 > z2:
            swap z1, z2
            swap a1, a2
        theta_obs = bond_angle(a1, v, a2)
        base_tag = tag_string_for(a1, v, a2)
        observations.append(
            (z1, zv, z2, theta_obs, base_tag,
             a1, v, a2))

# Phase 2: Cluster locally using the shared helper
# from 10a with tolerance = 0.0 (identity-only merge).
# This collapses bit-identical theta_obs values for the
# same (Z1, Zv, Z2) triplet into a single local record
# with obs_count > 1, keeping the template file compact
# without performing any interpretive merging.
# Non-identical observations -- even those differing by
# just 0.5 degrees -- remain as separate local types, so
# `normalize_types` in 10e sees the full raw resolution
# and can apply any angle_cluster_tolerance value the
# downstream condense.py simulation chooses.  DESIGN
# 4.8.10 explains why the tolerance is not a tunable
# parameter on the make_reactions.py side.
angle_types, angle_type_map =
    cluster_angles(observations, 0.0)

# Phase 3: Build local type tags with the cluster-mean
# theta_0 carried in the tag tail.  No K_angle is
# computed or stored and no local_angle_coeffs table is
# built here -- reaction templates do not carry angle
# coefficients, and normalize_types recomputes K
# authoritatively in 10e / DESIGN 4.8.8 item 4b from
# the triplet alone.  This is the only place where the
# template producer's output shape differs from 10c's:
# the lammps.dat producer builds local_angle_coeffs as
# intermediate storage for the Angle Coeffs section
# (consumed by the LAMMPS writer just below it), but
# the template producer has no analogous consumer.
num_local_angle_types = len(angle_types)
local_angle_tags   =
    [None] * (num_local_angle_types + 1)

for t = 1 to num_local_angle_types:
    atype = angle_types[t - 1]
    local_angle_tags[t] = (
        f"{atype.base_tag} "
        f"{atype.theta_0:.4f} {t}")

# Phase 4: Record per-atom angle entries with the
# local type ids.  normalize_types walks these
# and remaps the ids in 10f.
for i, obs in enumerate(observations):
    local_type_id = angle_type_map[i]
    angle_bonded[obs.v].append(
        [None, obs.a1, obs.a2])
    angle_tag_id[obs.v].append(local_type_id)

# Export to normalize_types:
#   source tag = "template:{name}"
#   local_angle_tags (no local_angle_coeffs -- see
#       Phase 3 rationale)
#   angle_bonded, angle_tag_id
#   per-local-type obs_count (slot 5 of each entry in
#       angle_types -- at tolerance=0 this counts only
#       bit-identical observations, typically small for
#       a single template but >1 wherever the template
#       has geometric duplicates such as a benzene
#       ring's six identical C-C-C angles)
```

### 10e. Cross-Source Angle Clustering (DESIGN 4.8.8 item 4a)

The first phase of angle handling inside
`normalize_types`.  Takes the per-source local cluster
centers emitted by 10c (lammps.dat) and 10d (each
reaction template) and merges any whose theta_0 values
represent the same physical angle.  This is what makes
bond/react type IDs consistent across sources.  The
algorithm is greedy merge with the same
`2 * tolerance` spread cap that 10a applies locally --
so local and cross-source clustering are semantically
consistent -- and adds observation-count weighting on
top, so a cluster anchored by many observations pulls
the final mean more strongly than a sparse one.

```
# -------------------------------------------------
# Data structures:
#
#   local_records : list, one entry per
#       (source, local_type_id) pair:
#           (z1, zv, z2,
#            theta_0_local,
#            obs_count,          # raw observations
#                                # feeding this
#                                # local cluster
#            base_tag,           # representative
#                                # tag prefix
#            source,             # "lammps.dat" or
#                                # "template:<name>"
#            local_type_id)
#
#   tolerance : float
#       angle_cluster_tolerance (default 5.0).
#
#   spread_cap : float
#       Max allowed total span (max-min) of
#       theta_0_local values within one final
#       cluster.  Default: 2.0 * tolerance.
#       Prevents greedy chaining from sweeping a
#       broad distribution into a single cluster.
#
#   Output:
#       final_types : list of
#           (z1, zv, z2,
#            theta_0_final,
#            obs_count_total,
#            representative_base_tag)
#       remap : dict
#           (source, local_type_id)
#               -> final_type_id
# -------------------------------------------------

function cross_source_cluster(local_records,
                              tolerance):
    # Group local records by canonical triplet.
    groups = {}
    for each rec in local_records:
        key = (rec.z1, rec.zv, rec.z2)
        groups[key].append(rec)

    final_types = []
    remap = {}
    spread_cap = 2.0 * tolerance

    for each key in groups:
        entries = groups[key]
        sort entries by theta_0_local ascending

        # Greedy merge, weighted by obs_count.
        cluster_w_sum = (entries[0].theta_0_local
                         * entries[0].obs_count)
        cluster_w     = entries[0].obs_count
        cluster_min   = entries[0].theta_0_local
        cluster_max   = entries[0].theta_0_local
        members       = [entries[0]]

        for i = 1 to len(entries) - 1:
            running_mean = cluster_w_sum / cluster_w
            candidate_theta = entries[i].theta_0_local
            new_max = max(cluster_max, candidate_theta)
            new_min = min(cluster_min, candidate_theta)

            within_tol =
                |candidate_theta - running_mean|
                    <= tolerance
            within_cap =
                (new_max - new_min) <= spread_cap

            if within_tol and within_cap:
                cluster_w_sum +=
                    candidate_theta * entries[i].obs_count
                cluster_w += entries[i].obs_count
                cluster_min = new_min
                cluster_max = new_max
                members.append(entries[i])
            else:
                finalize(key, members,
                         cluster_w_sum,
                         cluster_w,
                         final_types, remap)
                # Start a new cluster at entry i.
                cluster_w_sum = (candidate_theta
                    * entries[i].obs_count)
                cluster_w = entries[i].obs_count
                cluster_min = candidate_theta
                cluster_max = candidate_theta
                members = [entries[i]]

        finalize(key, members,
                 cluster_w_sum, cluster_w,
                 final_types, remap)

    return final_types, remap

function finalize(key, members,
                  cluster_w_sum, cluster_w,
                  final_types, remap):
    theta_0_final   = cluster_w_sum / cluster_w
    obs_count_total = cluster_w
    final_id = len(final_types) + 1
    # Take the first member's base_tag as the
    # representative prefix for the final type.
    # base_tag carries species_id / molecule_id
    # metadata that Z alone cannot reconstruct.
    representative_base_tag = members[0].base_tag
    final_types.append(
        (key.z1, key.zv, key.z2,
         theta_0_final, obs_count_total,
         representative_base_tag))
    for m in members:
        remap[(m.source, m.local_type_id)] =
            final_id
```

Decision notes embedded in the algorithm above:
- **Weighting.**  The running mean is `obs_count`-
  weighted, so a local cluster built from 200
  observations anchors the final theta_0 more
  strongly than one from 3.  This matches the
  physical intuition that the larger sample is a
  better estimator.
- **Spread cap.**  Greedy merge alone can chain
  across a wide distribution (e.g., observations at
  105, 107, 109, 111, 113 with tolerance 2.5 all
  collapse into one cluster spanning 8 degrees).
  The `spread_cap = 2 * tolerance` rule forces a
  split once the total span would exceed that
  bound, producing tighter clusters at distribution
  boundaries.
- **Canonical base_tag.**  Representative tag is
  taken from the first-encountered member rather
  than reconstructed from Z values, because the
  tag prefix carries species_id and molecule_id
  fields that Z alone does not encode.

### 10f. Tag Rewrite and Type-ID Remap (DESIGN 4.8.8 item 4c)

The second phase of angle handling inside
`normalize_types`, executed once 10e has produced
`(final_types, remap)`.  Every angle reference in every
source file is rewritten: the per-angle type id is
remapped to the global id, and the tag tail is replaced
with the final canonical theta_0 so any downstream tool
that inspects the tag sees a consistent value.  The
rewrite is deterministic given the cluster map, so
repeated runs on identical inputs produce byte-identical
output.

```
function rewrite_angles(sources, remap,
                        final_types,
                        angle_stiffness_coeff,
                        angle_parameter_scale):
    # Phase A: rewrite per-angle type ids in every
    # source.  lammps.dat carries an Angles section
    # with explicit type ids; each reaction
    # template carries a per-atom angle_tag_id
    # array.
    for each src in sources:
        if src is lammps.dat:
            for each angle entry in src.Angles:
                old_id = entry.type_id
                entry.type_id =
                    remap[(src.source_tag, old_id)]
        else:  # reaction template
            for each vertex atom v in src:
                for i in range(
                        len(src.angle_tag_id[v])):
                    old_id = src.angle_tag_id[v][i]
                    src.angle_tag_id[v][i] =
                        remap[(src.source_tag,
                               old_id)]

    # Phase B: build the unified global
    # unique_angle_tags table from final_types.
    # Each entry is
    #   "{representative_base_tag} {theta_0} {t}"
    # carrying the final canonical theta_0 and the
    # global type id.
    unique_angle_tags =
        [None] * (len(final_types) + 1)
    for t = 1 to len(final_types):
        ft = final_types[t - 1]
        unique_angle_tags[t] = (
            f"{ft.representative_base_tag} "
            f"{ft.theta_0_final:.4f} {t}")

    # Phase C: build the unified global
    # unique_angle_coeffs table via get_angle_k.
    # K_angle depends only on the triplet, so
    # recomputation here yields the same value
    # that any producer's local 10c/10d phase
    # computed -- cross-source merging does not
    # alter K_angle, only theta_0.
    unique_angle_coeffs =
        [None] * (len(final_types) + 1)
    for t = 1 to len(final_types):
        ft = final_types[t - 1]
        K = get_angle_k(
            ft.z1, ft.zv, ft.z2,
            angle_stiffness_coeff,
            angle_parameter_scale)
        unique_angle_coeffs[t] =
            [None, K, ft.theta_0_final]

    # Phase D: emit the cluster-map diagnostic.
    # For each final cluster, write:
    #   - global id
    #   - canonical theta_0_final
    #   - (z1, zv, z2)
    #   - every contributing
    #       (source, local_type_id,
    #        theta_0_local, obs_count) tuple
    # See DESIGN 4.8.8 item 4d.  This file is the
    # primary debuggability payback for routing
    # all clustering through normalize_types.
    write_cluster_map(final_types, remap,
                      local_records)

    return unique_angle_tags, unique_angle_coeffs
```

`normalize_types()`'s angle handling is therefore:
1. Gather `local_records` from every source.
2. `cross_source_cluster(local_records, tolerance)`
   -> `(final_types, remap)`.  (10e)
3. `rewrite_angles(sources, remap, final_types,
   angle_stiffness_coeff, angle_parameter_scale)`
   -> `(unique_angle_tags, unique_angle_coeffs)`.
   (10f)

No other changes are required inside
`normalize_types()`.  Bond handling (section 9) and
other type tables are unchanged.

---

## 11. Initial SCF Potential Database (DESIGN 5)

Five algorithms support the augmented initial-potential
database: the TOML reader (11.1), the hand-formatted
emitter (11.2), the runtime lookup invoked from
`makeinput.py` (11.3), the build pipeline
(11.4), and the validation harness (11.5). All live in
Python under `src/scripts/`. The Fortran side does not
change.

### 11.1 TOML Reader (DESIGN 5.2, 5.2.5, 5.4)

Parses a per-element `s_gaussian_pot.toml` file and
applies the validation rules from DESIGN 5.2
(schema v2).  Returns an `ElementDatabase`; raises a
clear error on any rule violation, naming the file
path, label, and field at fault.

Rule 9 (method must be a registered matcher) requires
knowledge of the active matcher registry, which lives
in `makeinput.py` (ARCHITECTURE 8.9).  To keep
`initial_potential_db.py` free of any import from
`makeinput.py`, `load()` accepts an optional
`known_methods` parameter -- a set of matcher-name
strings.  Callers that have a registry (`makeinput.py`,
`build_initial_potentials.py`) pass it in; callers
that do not (isolated unit tests) pass `None` and
rule 9 is skipped.  This decouples the library from
the registry without weakening the rule for real
runs.

Rule 1 (the schema version) is enforced by the version
gate of DESIGN 5.2.5 rather than by a bare equality
test, and it runs BEFORE the general required-field
sweep.  The ordering is the whole point: a file older
than this build is very likely to be missing a field
that sweep requires, and the gate exists so that such a
file is told it is out of date instead of being told a
symptom.  Only `schema_version`'s own presence is
checked ahead of the gate, since the gate cannot speak
without it.

```
CURRENT_SCHEMA_VERSION = 2

# The version table (DESIGN 5.2.5).  One entry per bump,
# keyed by the version that bump PRODUCES.  Each entry
# lists the fields that bump newly requires, mapping
# each to a derivation -- a function rewriting older
# parsed data into the newer shape -- or to
# NOT_DERIVABLE when no honest fill exists.  A field is
# NOT_DERIVABLE whenever an older file simply does not
# record what the new field asserts; a plausible default
# is not a derivation (DESIGN 5.2.5).
#
# Empty today: version 2 is current, and the v1 -> v2
# bump predates any database worth carrying forward, so
# no migration was ever written for it.  A v1 file
# therefore takes the "no migration path" refusal below,
# which is the reject-and-regenerate behaviour DESIGN
# 5.2 already specified -- now reported in version terms
# with a recovery named.
SCHEMA_MIGRATIONS = {}


function apply_schema_migrations(raw, path):
    found = raw["schema_version"]

    # Outcome 1: current.  Nothing to do.
    if found == CURRENT_SCHEMA_VERSION:
        return raw

    # Outcome 4: the file is ahead of the code.  The
    # recovery is the opposite of outcome 3's -- the
    # database is right and this build is behind -- so
    # say so rather than inviting a regeneration that
    # would DESTROY correct data.
    if found > CURRENT_SCHEMA_VERSION:
        fail(path, "schema_version " + str(found)
             + " was written by a newer Imago; this "
             + "build reads and writes version "
             + str(CURRENT_SCHEMA_VERSION)
             + ".  Update Imago to read this file; do "
             + "not regenerate it.")

    # Outcomes 2 and 3: replay every bump between the
    # file's version and ours, in ascending order, so a
    # file may cross several versions in one load.
    for version in range(found + 1,
                         CURRENT_SCHEMA_VERSION + 1):
        if version not in SCHEMA_MIGRATIONS:
            fail(path, "schema_version " + str(found)
                 + "; this build writes version "
                 + str(CURRENT_SCHEMA_VERSION)
                 + ", and no migration to version "
                 + str(version) + " exists.  "
                 + "Regenerate this file with "
                 + "build_initial_potentials.py.")

        for (field, rule) in
                SCHEMA_MIGRATIONS[version].derivations:
            # Outcome 3: refuse rather than guess.  An
            # invented value would leave the file
            # well-formed and WRONG, which no later
            # check could catch (DESIGN 5.2.5).
            if rule is NOT_DERIVABLE:
                fail(path, "schema_version " + str(found)
                     + "; this build writes version "
                     + str(CURRENT_SCHEMA_VERSION)
                     + ".  Field '" + field + "', "
                     + "required from version "
                     + str(version) + ", cannot be "
                     + "inferred from a version-"
                     + str(found) + " file.  "
                     + "Regenerate this file with "
                     + "build_initial_potentials.py.")
            raw = rule(raw)

    # Outcome 2 completed.  Stamp the migrated data with
    # the current version so the next save (11.2), which
    # emits db.schema_version, writes the file forward.
    # The producer refreshes the isolated baseline on
    # every run (11.4), so any file it touches is
    # rewritten current and the migration is paid once.
    raw["schema_version"] = CURRENT_SCHEMA_VERSION
    return raw


function load(path, known_methods = None):
    raw = tomllib.load(path)

    # Rule 1, via the version gate.  Its own presence
    # check comes first and separately: the gate speaks
    # in versions, so it needs the version to speak at
    # all.
    require("schema_version" in raw, path,
        "missing top-level field: schema_version")
    raw = apply_schema_migrations(raw, path)

    # Rule 3 (top-level half): every required
    # top-level key must be present.  Check before
    # any value-level rule so the error message
    # names the missing field rather than failing
    # later with a value mismatch.  Safe to run after
    # the gate, and only safe there: `raw` is now known
    # to be at the current version, so a field missing
    # HERE is a corrupt or hand-edited file, not an old
    # one, and the bare message is the right report.
    for f in ("element_symbol",
              "nuclear_z", "nuclear_alpha",
              "covalent_radius"):
        require(f in raw, path,
            "missing top-level field: " + f)
    # Rule 2: element symbol must match parent dir.
    expected_elem = basename(dirname(path))
    require(lower(raw["element_symbol"])
            == lower(expected_elem),
        path, "element_symbol does not match dir")

    db = ElementDatabase(
        schema_version  = raw["schema_version"],
        element_symbol  = raw["element_symbol"],
        # Z is coerced to a real: nominally integral, but
        # Imago consumes it as a real number.
        nuclear_z       = float(raw["nuclear_z"]),
        nuclear_alpha   = raw["nuclear_alpha"],
        covalent_radius = raw["covalent_radius"],
        potentials      = [])

    seen_labels    = set()
    default_count  = 0
    # File-wide half of rule 10: the canonical sub_spec each
    # method's preferred records must agree on.  Filled by the
    # first preferred record of a method, checked against by
    # every later one.
    preferred_subspec = {}          # method -> canonical sub_spec

    for entry_dict in raw.get("potential", []):
        # Per-entry required fields (rule 3,
        # per-entry half).  Check "label" first so
        # that subsequent error messages can name
        # the entry; the "label" check itself can
        # only cite the file path and the
        # [[potential]] index.
        require("label" in entry_dict, path,
            "[[potential]] missing field: label")
        lbl = entry_dict["label"]
        for f in ("default", "description",
                  "num_gaussians", "alpha_min",
                  "alpha_max", "coefficients",
                  "alphas", "provenance"):
            require(f in entry_dict, path, lbl,
                "missing field: " + f)

        # Length consistency (rule 4): coefficients and
        # alphas are per-coefficient arrays of length
        # num_gaussians.
        n = entry_dict["num_gaussians"]
        require(len(entry_dict["coefficients"]) == n
                and len(entry_dict["alphas"]) == n,
            path, lbl,
            "coefficients/alphas length"
            + " != num_gaussians")

        # Label uniqueness (rule 5)
        require(lbl not in seen_labels, path,
            "duplicate label: " + lbl)
        seen_labels.add(lbl)

        # Default tag counting (rule 7).  We'll
        # check the total after the loop so we can
        # report "zero" and "multiple" with the
        # same message structure.
        if entry_dict["default"]:
            default_count += 1

        # Provenance fields
        require_provenance(entry_dict["provenance"],
            path, lbl)

        # Fingerprint sub-blocks (rules 8, 9, 10).
        # Each [[potential.fingerprint]] record
        # contributes one FingerprintRecord.  The
        # method/sub_spec pair is unique per entry
        # (rule 8); the method must be known if a
        # registry was supplied (rule 9); and the
        # preferred flag is counted PER ENTRY, per
        # method (rule 10, below the loop).
        fingerprints = []
        seen_method_subspec = set()
        methods_on_entry = set()      # rule 10
        preferred_on_entry = {}       # method -> count
        entry_preferred_subspec = {}  # method -> canon sub_spec
        for fp_dict in entry_dict.get(
                "fingerprint", []):
            require("method" in fp_dict, path, lbl,
                "fingerprint missing field: method")
            require("sub_spec" in fp_dict, path, lbl,
                "fingerprint missing field:"
                + " sub_spec")
            method = fp_dict["method"]
            sub_spec = fp_dict["sub_spec"]
            canon = canonicalize_sub_spec(sub_spec)

            # Rule 8: per-entry (method, sub_spec)
            # uniqueness.
            key = (method, canon)
            require(key not in seen_method_subspec,
                path, lbl,
                "duplicate fingerprint"
                + " (method=" + method
                + ", sub_spec=" + str(canon)
                + ")")
            seen_method_subspec.add(key)

            # Rule 9: method must be registered.
            # Skipped when known_methods is None
            # (test contexts without a registry).
            if known_methods is not None:
                require(method in known_methods,
                    path, lbl,
                    "fingerprint method '" + method
                    + "' not in matcher registry")

            # Rule 10 bookkeeping.  `preferred` is
            # optional and defaults to false: a record
            # nobody flagged is an alternate sub_spec
            # riding along, never the canonical one.
            # Both halves of the rule are checked below
            # the loop, per entry, in that order.
            preferred = fp_dict.get("preferred", False)
            methods_on_entry.add(method)
            if preferred:
                preferred_on_entry[method] = (
                    preferred_on_entry.get(method, 0) + 1)
                entry_preferred_subspec[method] = canon

            # Payload = all keys other than method,
            # sub_spec, and preferred.  Matchers validate
            # their own payload shape at lookup time.
            payload = {k: v for (k, v) in fp_dict
                       if k not in ("method", "sub_spec",
                                    "preferred")}
            fingerprints.append(FingerprintRecord(
                method    = method,
                sub_spec  = sub_spec,
                preferred = preferred,
                payload   = payload))

        # Rule 10, first half: PER ENTRY.  For each method
        # present on THIS entry, exactly one of its records is
        # preferred.  The flag marks this entry's canonical
        # record for that family: the consumer reads any one to
        # learn which sub_spec to query at, and the dedup
        # (5.2.3) asks EACH entry for its own canonical
        # bispectrum -- so every harvested entry flags its own.
        #
        # An entry with no fingerprints (the "isolated"
        # baseline) has no method present and is vacuously
        # exempt: the loop below does not run.
        #
        # This is checked BEFORE the file-wide half, and the
        # order is load-bearing.  Rule 8 already forbids two
        # records sharing a (method, sub_spec) on one entry, so
        # two FLAGGED records of one method necessarily differ
        # in sub_spec.  Were the file-wide agreement checked
        # first, it would fire on them and blame a disagreement
        # "across the file" for what is ambiguity within a
        # single entry.  Each failure gets its own accurate
        # message this way.
        for method in sorted(methods_on_entry):
            count = preferred_on_entry.get(method, 0)
            require(count == 1, path, lbl,
                "method '" + method + "' is present on this"
                + " entry with " + str(count) + " record(s)"
                + " flagged preferred = true; exactly one is"
                + " required so the entry has an unambiguous"
                + " canonical record for that family")

        # Rule 10, second half: FILE-WIDE.  Every preferred
        # record of a method must name the SAME sub_spec.  That
        # agreement is what the flag MEANS -- it names the
        # settings the consumer computes its query with (DESIGN
        # 5.6.5 step 2) -- and two entries flagging different
        # settings would leave no canonical answer.  The first
        # entry to flag a method fixes its sub_spec; every later
        # entry is compared against that.
        for method in sorted(entry_preferred_subspec):
            canon = entry_preferred_subspec[method]
            if method in preferred_subspec:
                require(preferred_subspec[method] == canon,
                    path, lbl,
                    "preferred '" + method + "' records"
                    + " disagree on sub_spec across the file;"
                    + " the flag names the canonical settings,"
                    + " so every preferred record of a method"
                    + " must share one sub_spec")
            else:
                preferred_subspec[method] = canon

        db.potentials.append(PotentialEntry(
            label         = lbl,
            default       = entry_dict["default"],
            description   = entry_dict["description"],
            num_gaussians = n,
            alpha_min     = entry_dict["alpha_min"],
            alpha_max     = entry_dict["alpha_max"],
            coefficients  = entry_dict["coefficients"],
            alphas        = entry_dict["alphas"],
            provenance    = entry_dict["provenance"],
            fingerprints  = fingerprints))

    # File-level rule 6
    require("isolated" in seen_labels, path,
        "missing required 'isolated' baseline entry")

    # File-level rule 7: exactly one default tag.
    require(default_count == 1, path,
        "expected exactly one entry with"
        + " default = true; found "
        + str(default_count))

    return db


function canonicalize_sub_spec(sub_spec):
    # Sub-spec equality for rule 8 is a deep dict
    # comparison after canonicalization: keys sorted,
    # numeric types normalized (int versus float
    # treated as equal when their value is equal),
    # nested dicts canonicalized recursively.  The
    # returned object must be hashable so it can
    # live in a set: a frozenset of (key, value)
    # pairs with values canonicalized to tuples for
    # nested structures.
    return freeze_dict(sub_spec)


function require_provenance(prov, path, lbl):
    for f in ("source", "commit", "generated_at"):
        require(f in prov, path, lbl,
            "provenance missing: " + f)
    require(prov["source"] in ("atomSCF", "Imago"),
        path, lbl,
        "provenance.source must be"
        + " 'atomSCF' or 'Imago'")
    if prov["source"] == "Imago":
        for f in ("reference_id", "atom_site",
                  "kpoint_spec",
                  "scf_threshold",
                  "scf_iterations"):
            require(f in prov, path, lbl,
                "Imago provenance missing: " + f)


function lookup(db, label):
    for entry in db.potentials:
        if entry.label == label:
            return entry
    raise KeyError(label)


function baseline(db):
    # The "isolated" entry, guaranteed by rule 6.
    # Used by the validation harness (11.5).
    return lookup(db, "isolated")


function default_entry(db):
    # The entry with default == true, guaranteed
    # unique by rule 7.  Used by the consumer
    # (11.3) on the no-scheme fallback branch.
    for entry in db.potentials:
        if entry.default:
            return entry
    error("internal: load() must enforce rule 7")


function find_fingerprint(entry, method, sub_spec):
    # Return the FingerprintRecord on `entry` whose
    # (method, sub_spec) matches; sub_spec
    # comparison uses the same canonicalization as
    # rule 8.  Raises KeyError on miss; callers
    # decide whether the miss is fatal (consumer
    # preflight, 11.3) or expected (producer
    # checking whether to overwrite, 11.4).
    target = canonicalize_sub_spec(sub_spec)
    for fp in entry.fingerprints:
        if (fp.method == method
                and canonicalize_sub_spec(fp.sub_spec)
                    == target):
            return fp
    raise KeyError((method, sub_spec))
```

---

### 11.2 TOML Emitter (DESIGN 5.5)

Deterministic hand-formatted writer: given an
`ElementDatabase`, produces byte-identical file
contents.  This is the bit-level guarantee --
formatting never introduces spurious diff churn.
Determinism is achieved by fixed key ordering, fixed
indentation, fixed float format (`%.16e`), and
per-block `=` alignment.  At the pipeline level
(11.4), file-level byte-identity across runs is not
promised: provenance timestamps refresh, and SCF /
fit numerical drift can perturb the numbers
themselves.  Any real diff between two builds
isolates such changes from formatting noise.

```
function save(db, path):
    out = []

    # Top-level block.  schema_version always writes
    # as 2 -- the emitter is paired with the schema
    # version it understands.  Use db.schema_version
    # rather than a literal so v1->v2 migration is a
    # single-place change.
    top_keys = ["schema_version", "element_symbol",
                "nuclear_z", "nuclear_alpha",
                "covalent_radius"]
    out.extend(format_block(db.__dict__, top_keys))
    out.append("")

    for entry in db.potentials:
        out.append("[[potential]]")

        # `default` slots between `label` and
        # `description` (matches DESIGN 5.3 sketch).
        body_keys = ["label", "default", "description",
                     "num_gaussians", "alpha_min",
                     "alpha_max"]
        # Width spans body keys plus the array keys
        # so the array openers align with the rest.
        align_keys = body_keys + ["coefficients",
                                  "alphas"]
        width = max(len(k) for k in align_keys)

        for k in body_keys:
            out.append(format_kv(
                k, entry.__dict__[k], width))
        out.append(format_array_open(
            "coefficients", width))
        for x in entry.coefficients:
            out.append("   " + fmt_float(x) + ",")
        out.append("]")
        out.append(format_array_open(
            "alphas", width))
        for x in entry.alphas:
            out.append("   " + fmt_float(x) + ",")
        out.append("]")
        out.append("")

        # Provenance block
        out.append("[potential.provenance]")
        prov_keys = ordered_provenance_keys(
            entry.provenance)
        out.extend(format_block(
            entry.provenance, prov_keys))
        out.append("")

        # Fingerprint sub-blocks (v2).  Emitted in
        # insertion order, which matches the reader's
        # parse order; this keeps round-trip
        # load(save(db)) byte-deterministic without
        # imposing a method-and-sub_spec ordering
        # the producer would otherwise have to
        # canonicalize.  Empty fingerprint list
        # produces no blocks and no extra blank
        # lines.
        for fp in entry.fingerprints:
            emit_fingerprint_block(out, fp)

    # Trim trailing blanks; ensure exactly one newline
    while out and out[-1] == "":
        out.pop()
    write_file(path, "\n".join(out) + "\n")


function emit_fingerprint_block(out, fp):
    out.append("[[potential.fingerprint]]")

    # method and sub_spec come first, in that order,
    # then `preferred` when the record carries it, then
    # the payload keys in the payload dict's iteration
    # order.  Width alignment spans method/sub_spec plus
    # the payload's scalar and multi-line keys, so all
    # `=` signs align within the block.
    #
    # `preferred` is emitted ONLY when true, and joins the
    # alignment key set only then: a non-preferred record
    # omits the line entirely (the flag defaults to false
    # on read, 11.1) and must not pad the block to a width
    # it never uses.  It must be emitted, though -- rule 10
    # is a statement about the FILE, so a flag the emitter
    # dropped would make every saved database fail its own
    # loader on the next read.
    fixed_keys   = ["method", "sub_spec"]
    if fp.preferred:
        fixed_keys.append("preferred")
    payload_keys = list(fp.payload.keys())
    align_keys   = fixed_keys + payload_keys
    width        = max(len(k) for k in align_keys)

    out.append(format_kv("method", fp.method, width))
    out.append(format_kv(
        "sub_spec",
        format_inline_table(fp.sub_spec),
        width))
    if fp.preferred:
        out.append(format_kv("preferred", "true", width))

    for k in payload_keys:
        v = fp.payload[k]
        if is_float_list(v):
            # Multi-line float array (matches the
            # coefficients/alphas layout).
            out.append(format_array_open(k, width))
            for x in v:
                out.append("   " + fmt_float(x) + ",")
            out.append("]")
        elif isinstance(v, dict):
            # Inline-table payload field (e.g.,
            # reduce shell_code).  Single line.
            out.append(format_kv(
                k, format_inline_table(v), width))
        else:
            out.append(format_kv(k, v, width))

    out.append("")


function format_inline_table(d):
    # Deterministic inline-table emission: keys in
    # alphabetical order, one space inside braces,
    # ` = ` between key and value, `, ` between
    # pairs.  Nested values use format_scalar
    # recursively; nested tables emit as inline
    # tables themselves.
    parts = []
    for k in sorted(d.keys()):
        v = d[k]
        if isinstance(v, dict):
            rhs = format_inline_table(v)
        else:
            rhs = format_scalar(v)
        parts.append(k + " = " + rhs)
    return "{ " + ", ".join(parts) + " }"


function format_block(d, keys):
    width = max(len(k) for k in keys)
    return [format_kv(k, d[k], width) for k in keys]


function format_kv(key, value, width):
    # `value` here is either a TOML-ready string
    # (from format_inline_table) or a Python scalar.
    if isinstance(value, str) and starts_with_brace(
            value):
        rhs = value          # already a TOML literal
    else:
        rhs = format_scalar(value)
    return pad_right(key, width) + " = " + rhs


function format_array_open(key, width):
    return pad_right(key, width) + " = ["


function format_scalar(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, int):   return str(v)
    if isinstance(v, float): return fmt_float(v)
    if isinstance(v, str):   return toml_quote(v)
    error("unsupported scalar type: " + type(v))


function fmt_float(x):
    return sprintf("%.16e", x)


function ordered_provenance_keys(prov):
    base = ["source", "commit", "generated_at"]
    extras = ["reference_id", "atom_site",
              "kpoint_spec", "scf_threshold",
              "scf_iterations"]
    if prov["source"] == "Imago":
        return base + extras
    return base


function is_float_list(v):
    return (isinstance(v, list)
            and len(v) > 0
            and all(isinstance(x, float) for x in v))
```

---

### 11.3 makeinput.py Lookup (DESIGN 5.6, 5.10)

The consumer side of the augmented database splits into
seven concerns, one per sub-section below.  The driver
(11.3.g) chains them together; the matcher protocol
(11.3.a) is the dispatch surface that lets the species
pass and the entry pick stay agnostic of which
descriptor family is in play.

Fingerprint matching runs by default (DESIGN 5.6.5, the
C93 decoupling): the entry pick applies to *every* species,
no matter how it was grouped.  A reduced subset of this flow
is still common, though -- sub-section 11.3.0 pins it: the
path active when fingerprint matching is turned off
(`-nofingerprint`) or no element's database carries a usable
(preferred) fingerprint, so the pick collapses to the `-pot`
override and the default tag.  It is also the path the first
consumer milestone (C47) implemented before C54 onward layered
in the matcher machinery; reading it first makes the seven
full-flow sub-sections easier to place.

The mapping from DESIGN to PSEUDOCODE here is one-to-
one:

  PSEUDOCODE   DESIGN
  -----------------------------------------------------
  11.3.a       8.9 + 5.10.5 (matcher protocol and
               LOEN parameter contract)
  11.3.b       5.6.3        (per-element preflight,
                             including the rule-4
                             coverage note)
  11.3.c       5.6.4        (species pass and scope
                             resolution)
  11.3.d       5.6.5        (manifest-entry pick per
                             species)
  11.3.e       5.6.6        (type pass and electronic-
                             state perturbation)
  11.3.f       5.10         (makegroups.py bispectrum
                             grouping)
  11.3.g       5.6.7        (driver and on-the-wire
                             emit)

---

#### 11.3.0 Reduced Flow (no fingerprint pick; DESIGN 5.6)

The seven sub-sections above specify the full Phase-2
selection flow.  A reduced path skips the fingerprint pick
entirely and collapses to the `-pot` override plus the
default tag.  It is taken when EITHER the user passed
`-nofingerprint` (an explicit opt-out, DESIGN 5.6.1) OR no
element in the structure carries a usable (preferred)
fingerprint, so there is nothing to match against.  This
sub-section names which branches are live in that path so the
reduced consumer has a self-contained spec, separate from the
matcher machinery that C54 and onward layer in.

In the reduced path:

  - The preflight (11.3.b) still loads each element's
    database and still marks missing ones for the legacy
    path.  `noteCoverage` may still fire (info-level) to
    report an element with no fingerprints, unless
    `-nofingerprint` suppressed it (DESIGN 5.6.3 step 4).
  - Grouping (11.3.c) is unaffected -- crystallographic,
    position-based (`-target`, `-block`), or an explicit
    environment scheme still partition the atoms; the reduced
    path is about the PICK, not the grouping.
  - The entry pick (11.3.d) skips precedence 2: with the
    fingerprint match disabled or unmatchable, only
    precedence 1 (`-pot LABEL`) and precedence 3
    (`default_entry`) survive.

An element whose augmented file is absent is handled the
same way as in the full flow: the preflight marks
`databases[elem] = None` and the driver (11.3.g, steps 6
and 8) emits it via the legacy `pot1`/`coeff1` path, with
no library entry consulted.  There is no schema-v1 case
to consider here -- the reader's version gate (11.1)
either brings an older file up to the current version or
refuses it outright, and the producer (11.4) writes every
on-disk file at the current version, so a loaded database
is always current and always carries a `default` tag.

The reduced entry pick is the fingerprint-disabled
restriction of 11.3.d:

```
function pickEntryReduced(db, pot_override):
    # Precedence 1: manual -pot override.  A KeyError
    # here is fatal, exactly as in 11.3.d -- a
    # deliberate user choice must never silently fall
    # back to a different potential.
    if pot_override is not None:
        try:
            return lookup(db, pot_override)
        except KeyError:
            error("-pot " + pot_override + " not found"
                + " in this element's database; a"
                + " manual override must match a label")

    # Precedence 3: the default-tagged entry, guaranteed
    # to exist and be unique by validation rule 7.  With
    # the fingerprint pick off (or nothing to match), there
    # is no precedence-2 step to attempt first.
    return default_entry(db)
```

**Carry-forward to the full flow.**  The full 11.3.d simply
inserts precedence 2 between the override and the default;
this reduced pick is its `-nofingerprint`-or-no-match
restriction.  No branch written for the reduced flow is
rewritten -- the full flow only adds the middle precedence.

---

#### 11.3.a Matcher Protocol (DESIGN 5.6, 5.10.5; ARCH 8.9)

Each matcher knows one descriptor family.  The species
pass calls `compute_query` and `distance` to bucket
atoms (11.3.c); the entry pick calls `representative`
and `distance` against per-entry fingerprints (11.3.d);
the producer (11.4) and the `makegroups` bispectrum flow
(11.3.f) both call `to_loen_input` and `parse_loen_output`
on matchers whose `needs_loen_run` is true.  The protocol
isolates Imago's Fortran side from manifest-schema
growth: a new descriptor family is a new class plus a
new `MATCHERS` registry entry.

The entry pick also needs a query for a `needs_loen_run`
matcher in its file-dictated branch (11.3.d): when the
database's preferred record for a crystalline/pre-assigned
species is a bispectrum one, the pick obtains the per-atom
descriptors through `loen_descriptors(structure, matcher,
sub_spec)` -- a thin reuse-or-run wrapper over the same
loen seam (`to_loen_input` -> run `imago.py -loen -scf no`
-> `parse_loen_output`).  It reuses a loen run already
performed for grouping (the makegroups `fort.21`, 11.3.f)
when one exists for this `(method, sub_spec)`, and
otherwise triggers one fast loen run for the pick.  Python-
side matchers (reduce) never need it -- `compute_query`
serves the pick in-process.

```
class Matcher:
    # Protocol surface; subclasses fill in everything.
    name                     = ""        # set below
    needs_loen_run           = False     # set below
    default_similarity_floor = 0.0       # set below
    active_sub_spec          = None      # bound by
                                         # argparse to
                                         # the user's
                                         # CLI sub_spec

    function compute_query(structure, sub_spec):
        # Return one fingerprint vector per atom of
        # the WHOLE structure, in site-index order.
        # The list length always equals the structure's
        # atom count; element filtering happens at
        # the call site (11.3.c) where atoms_in_scope
        # already encodes element and spatial scope.
        # Python-side matchers (reduce) compute
        # in-process.  Loen-side matchers (bispectrum) do
        # NOT implement this: their vectors come from the
        # makegroups sequential loen flow (11.3.f), read
        # off fort.21 by parse_loen_output, which produces
        # one row per potential site of the whole
        # structure.
        abstract

    function distance(vec_a, vec_b):
        # Symmetric, non-negative scalar distance in
        # this matcher's descriptor space.  Used for
        # both species bucketing (11.3.c) and the
        # manifest-entry similarity test (11.3.d).
        abstract

    function representative(members):
        # Reduce a list of member-atom fingerprint
        # vectors into one fingerprint that represents
        # the whole species (11.3.d step 2).  Each
        # subclass chooses semantics appropriate to
        # its descriptor space; the protocol pins
        # only the shape (members in, one vector out).
        abstract

    function extract_query_vector(payload):
        # Read the query vector out of a
        # FingerprintRecord payload (5.4).  Each
        # matcher knows its own field name -- DESIGN
        # 5.2 documents the per-matcher payload
        # shape (bispectrum uses `values`, reduce
        # uses `shell_code`).  Returned in whatever
        # form this matcher's distance() expects.
        abstract

    function build_payload(query_vector):
        # Inverse of extract_query_vector.  Wraps a
        # freshly-harvested vector in the payload
        # dict the matcher's records use, so the
        # producer (11.4) can attach it to a
        # FingerprintRecord and the consumer (11.3.d)
        # can read it back symmetrically.
        abstract

    function to_loen_input(sub_spec):
        # Translate the user's sub_spec into the
        # LOEN_INPUT_DATA parameter dict that
        # O_Input::readLoEnControl consumes.  Only
        # meaningful when needs_loen_run is true.
        abstract

    function parse_loen_output(path, sub_spec):
        # Read fort.21 written by a loen run; return
        # per-site fingerprint vectors in site-index
        # order, one per potential site of the whole
        # structure.  Only meaningful when
        # needs_loen_run is true.
        abstract


class ReduceMatcher extends Matcher:
    name                     = "reduce"
    needs_loen_run           = False
    default_similarity_floor = 0.05      # matches the
                                         # tolerance of
                                         # the existing
                                         # group_reduce
                                         # path

    function compute_query(structure, sub_spec):
        # Shells over a PERIODIC NEIGHBOUR LIST
        # (DESIGN 5.11), for EVERY atom in the
        # structure (not just one element).  Returned
        # in site-index order so 11.3.c can index by
        # full-structure atom index.  The species-pass
        # filter on method.element handles per-element
        # selection at the call site; computing all
        # atoms keeps the matcher contract uniform
        # across Python-side and loen-side families
        # (loen naturally writes one row per site of
        # the whole structure).
        fingerprints = [None] * (structure.num_atoms + 1)
        for atom in 1 .. structure.num_atoms:
            fingerprints[atom] = shellCode(
                structure, atom, sub_spec)
        return fingerprints


function shellCode(structure, atom, sub_spec):
    # One atom's concentric shells (DESIGN 5.11).

    # THE NEIGHBOUR LIST.  Every periodic IMAGE within
    # the cutoff is a neighbour, counted once per image
    # -- including images of the central atom itself,
    # which are ordinary neighbours in space (in an fcc
    # lattice the whole second shell of a site is images
    # of that site).  Only the atom at distance zero is
    # excluded: an atom is not its own neighbour.
    #
    # Counting IMAGES rather than cell atoms is what
    # makes the descriptor transferable, which is the
    # property 5.2 relies on when it matches a stored
    # shell code against another structure.  The
    # extended coordinates and the image -> central-atom
    # map are already built by createMinDistMatrix, so
    # no new geometry is computed here.
    neighbours = []          # (distance, central atom)
    for image in 1 .. structure.num_atoms_ext:
        distance = norm(structure.direct_xyz[atom]
                        - structure.ext_direct_xyz_list[image])
        if 0 < distance <= sub_spec.cutoff:
            neighbours.append((distance,
                structure.ext_to_central_item_map[image]))

    # THE WALK.  Build shells outward, seeding each level
    # at the closest neighbour not yet assigned and
    # sweeping the [seed, seed + thick] band into it.
    assigned       = [0] * len(neighbours)
    level_distance = [None] * (sub_spec.level + 1)
    levels         = [None] * (sub_spec.level + 1)

    for level in 1 .. sub_spec.level:
        seed = index of the unassigned neighbour with the
               smallest distance, or None if all assigned

        # EXHAUSTION.  A cutoff too small to reach the
        # requested levels leaves one with nothing to
        # seed it.  REFUSE, naming the level and the
        # cutoff: an empty shell is a value the walk did
        # not find, and inventing one would store a
        # descriptor no structure produced (5.11).
        if seed is None:
            error("reduce level " + level + " has no "
                  "neighbour to seed it within cutoff "
                  + sub_spec.cutoff)

        level_distance[level] = neighbours[seed].distance
        for i in 0 .. len(neighbours) - 1:
            d = neighbours[i].distance
            if (d >= level_distance[level] and
                    d <= level_distance[level]
                         + sub_spec.thick and
                    d <= sub_spec.cutoff):
                assigned[i] = level

    # Each shell records its seed distance and the
    # neighbours in it.  Within one structure the
    # multiset carries (element, species) -- species
    # distinguishes atoms here; the multiset STORED in
    # the database carries element symbols only, since
    # species numbering is local to a structure and
    # would not transfer (5.2).
    for level in 1 .. sub_spec.level:
        members = [(structure.atom_element_id[c],
                    structure.atom_species_id[c])
                   for (d, c) in neighbours
                   where assigned[that neighbour] == level]
        names   = [structure.atom_element_name[c] ...]
        levels[level] = ReduceShellLevel(
            level_distance[level], members, names)

    return ReduceShellCode(
        structure.atom_element_id[atom],
        structure.atom_element_name[atom],
        sub_spec.tolerance, levels)

    function distance(a, b):
        # Hamming-like comparison of two shell-code
        # vectors: count of slots that differ beyond
        # sub_spec["tolerance"].  Matches the
        # comparison group_reduce performs internally
        # so refactoring it behind this surface
        # preserves behavior exactly.
        return shell_code_distance(a, b,
            sub_spec.get("tolerance", 0.05))

    function representative(members):
        # All intra-species members of a reduce bucket
        # agree within tolerance by construction.
        # Returning the first member is correct;
        # which member never affects downstream
        # distance comparisons in the matcher's space.
        return members[0]

    function extract_query_vector(payload):
        # Reduce records carry the shell-code in the
        # `shell_code` field (DESIGN 5.2 / 5.4).
        return payload["shell_code"]

    function build_payload(shell_code):
        # Serialize the in-memory shell code into the
        # element-only, cross-structure form stored on
        # disk (DESIGN 5.2): the central atom's element
        # symbol plus one entry per level holding the
        # shell distance and the list of neighbor element
        # symbols.  The structure-local integer ids and
        # the neighbor species are dropped here -- species
        # numbering does not transfer across structures,
        # so only the transferable element symbols are
        # kept.  All symbols are lowercased.
        # shell_code.levels is 1-indexed with a None
        # placeholder at slot 0, so iterate from slot 1.
        return {"shell_code": {
            "element": lower(shell_code.element_name),
            "levels": [
                {"distance":  level.distance,
                 "neighbors": lower_each(
                     level.member_names)}
                for level in shell_code.levels[1:]]}}

    function to_loen_input(sub_spec):
        error("ReduceMatcher is Python-side; no LOEN"
            + " input is built")

    function parse_loen_output(path, sub_spec):
        error("ReduceMatcher is Python-side; no"
            + " fort.21 is parsed")


class BispecMatcher extends Matcher:
    name                     = "bispectrum"
    needs_loen_run           = True
    default_similarity_floor = 0.10      # heuristic;
                                         # overridable
                                         # per scheme
                                         # on the CLI

    # No compute_query: bispectrum vectors are not
    # produced in-process.  The makegroups sequential
    # loen flow (11.3.f) runs Imago and reads the
    # resulting fort.21 through parse_loen_output below;
    # the orchestrator, not this matcher, drives that
    # sequence.  (Element-aware mode is deferred to TODO
    # C62 / D10, where to_loen_input gains by_element.)

    function distance(a, b):
        # Euclidean distance between bispectrum
        # vectors of length twoj2 + 1.  Symmetric,
        # cheap to compute, and consistent with the
        # element-wise mean used by `representative`
        # below.
        return l2_norm(vector_subtract(a, b))

    function representative(members):
        # Element-wise arithmetic mean of the member
        # vectors.  All members share the same length
        # (twoj2 + 1) because they come from one loen
        # run under one sub_spec, so the mean is
        # well-defined slot by slot.
        n = len(members)
        return [sum(m[i] for m in members) / n
                for i in range(len(members[0]))]

    function to_loen_input(sub_spec):
        # Parameter contract per DESIGN 5.10.5.
        # Required keys: twoj1, twoj2.  Optional
        # keys carry the defaults that match the
        # currently-hardcoded LOEN_INPUT_DATA block
        # makeinput.py emits today.
        require("twoj1" in sub_spec,
            "BispecMatcher requires sub_spec[twoj1]")
        require("twoj2" in sub_spec,
            "BispecMatcher requires sub_spec[twoj2]")
        return {
            "loenCode"     : 1,
            "twoj1"        : sub_spec["twoj1"],
            "twoj2"        : sub_spec["twoj2"],
            # cutoff (Bohr) must enclose every atom's
            #   first shell or that atom gets an all-zero
            #   descriptor; max_neigh caps the per-site
            #   list and must fit that reach (DESIGN 5.10.5).
            "max_neigh"    : sub_spec.get(
                "max_neigh", 50),
            "cutoff"       : sub_spec.get(
                "cutoff", 9.0),
            "angleSqueeze" : sub_spec.get(
                "angle_squeeze", 0.85)}

    function parse_loen_output(path, sub_spec):
        # fort.21 (DESIGN 5.10.3): a HEADER line, then one
        # row per potential site of the whole structure in
        # site-index order.  Each row leads with identity
        # columns -- site#, element, species,
        # type_in_species, type_flat -- then twoj2 + 1
        # real bispectrum values, then a trailing sum the
        # matcher ignores.  Skip the header; from each data
        # row return the identity fields plus the
        # twoj2 + 1 components.  (The orchestrator,
        # 11.3.f, uses the identity fields to map a row to
        # its atom/type without a separate datSkl.map.)
        n_slots   = 2 * sub_spec["twoj2"] + 1
        data_rows = drop_header(read_text_rows(path))
        return [parse_identity_and_components(r, n_slots)
                for r in data_rows]

    function extract_query_vector(payload):
        # Bispectrum records carry the vector in the
        # `values` field (DESIGN 5.2 / 5.3 sketch).
        return payload["values"]

    function build_payload(query_vector):
        # Stored under `values` (DESIGN 5.2 / 5.4); copy
        # into a fresh list so the record never aliases a
        # vector the caller may still reuse.
        return {"values": list(query_vector)}


# Module-level registry (ARCHITECTURE 8.9).
# initial_potential_db.load() consults MATCHERS.keys()
# when enforcing per-element-database rule 9 ("method
# must be a known matcher").  Adding a new descriptor
# family is a new class plus a new entry here; no
# other code path needs to change.
MATCHERS = {
    "reduce"     : ReduceMatcher,
    "bispectrum" : BispecMatcher,
}
```

---

#### 11.3.b Per-Element Preflight (DESIGN 5.6.3)

Runs once before the species pass starts.  Loads the
augmented database for every element in the structure,
marks elements without a database for the legacy
fallback path, and -- when fingerprint matching is enabled
(the default; `-nofingerprint` turns it off) -- notes
(info-level, never fatal) any element whose database
carries no fingerprint records at all.  Such an element
still groups normally; its species simply fall through to
the default-tagged entry at the per-species pick (11.3.d
step 3).  We do not abort: the fingerprint pick is a bonus
layered on grouping that the default entry always backstops
(DESIGN 5.6.3 step 4).  (The bispectrum *grouping* path
reports its own loen-coverage condition in makegroups,
11.3.f; that is about grouping, not this pick coverage.)

```
function perElementPreflight(structure,
        fingerprinting_enabled):
    # `fingerprinting_enabled` is False under -nofingerprint.
    # Returns a mapping elem -> ElementDatabase | None; None
    # marks the legacy-fallback path that 11.3.g consumes.
    databases = {}
    elements = unique_element_symbols(structure)

    for elem in elements:
        path = ("share/atomicPDB/" + lower(elem)
                + "/s_gaussian_pot.toml")

        if not file_exists(path):
            info("augmented database not yet"
                + " populated for " + elem
                + "; using legacy pot1/coeff1 for"
                + " this element")
            databases[elem] = None    # legacy marker
            continue

        # Passing MATCHERS.keys() enables rule 9 (any
        # fingerprint method must be a registered
        # matcher).  Without it the loader skips that
        # rule -- desirable for unit tests but not
        # for a real consumer run.
        databases[elem] = load(path,
            known_methods = MATCHERS.keys())

        # Coverage note (5.6.3 step 4).  Independent of the
        # grouping scheme now (C93): the pick runs for every
        # species, so the only question is whether this
        # element has any fingerprint to match.  Suppressed
        # under -nofingerprint, where the default is the
        # deliberate choice.  Never fatal.
        if fingerprinting_enabled:
            noteCoverage(databases[elem], elem, path)

    return databases


function noteCoverage(db, elem, path):
    # Emit an info note when NO entry in `db` carries ANY
    # FingerprintRecord (of any registered matcher), so its
    # species can only take the default potential.  This does
    # not abort: the per-species pick (11.3.d step 3) falls
    # through to the default-tagged entry.
    for entry in db.potentials:
        if len(entry.fingerprints) > 0:
            return       # coverage exists; nothing to note
    info("element " + elem + " has no fingerprint records"
        + " in " + path + "; its species will use the"
        + " default potential.  To enable a fingerprint"
        + " match, add a fingerprint declaration to the"
        + " curation manifest and re-run the producer"
        + " (DESIGN 5.7).")
```

---

#### 11.3.c Species Pass (DESIGN 5.6.4)

Walks `settings.methods` in CLI order, exactly the
dispatch surface today's `assign_group` uses.  Position-
based flags carry a `name=NAME` keyword that registers
the spatial region for later `scope=NAME` references;
environment-based flags resolve their scope, request
per-atom fingerprints from the active matcher, and
bucket the in-scope atoms by descriptor distance.
Atoms outside the active scope keep whatever species ID
earlier flags produced.

With the C93 decoupling the entry pick (11.3.d) runs for
*every* species regardless of how it was grouped, so the
species pass no longer tracks an `env_species_ids` set: it
returns only the per-atom species assignment and the named
spatial regions.  The per-atom grouping descriptors stay
with the driver (11.3.g), which passes them to the pick so
the explicit-scheme regime can reuse them for its
representative.

```
function speciesPass(structure, settings, databases,
        atom_fingerprints):
    # `atom_fingerprints` is None when the user selected
    # no environment-based grouping scheme.  When non-None
    # it holds one vector per atom of the whole structure
    # (in site-index order) from the grouping matcher; the
    # bucketing step indexes into it by full-structure
    # index, and the driver later hands the same vectors
    # to the pick (11.3.d) for descriptor reuse.
    n_atoms          = len(structure.atoms)
    atom_species_id  = [1] * n_atoms
    named_regions    = {}     # name -> atom-index set

    for method in settings.methods:
        if method.op == "spatial":
            # Position-based (-target, -block).
            # Existing geometric grouping; the new
            # `name=` keyword is consumed here.
            in_region = compute_spatial_membership(
                structure, method)
            assign_new_species(atom_species_id,
                in_region)
            if method.name is not None:
                named_regions[method.name] = \
                    in_region

        elif method.op == "environment":
            # Environment-based (-reduce, -bispec).
            # Resolve scope, bucket the in-scope
            # atoms by matcher distance.  Mutual
            # exclusion (DESIGN 5.6.2) means at most
            # one environment method appears in this
            # loop.  The grouped species flow on to
            # the pick like any other -- no env-only
            # gate (C93).
            scope = resolve_scope(method.scope,
                named_regions, n_atoms)
            atoms_in_scope = \
                atoms_of_element_in_scope(structure,
                    method.element, scope)
            assign_species_by_bucketing(
                atom_species_id,
                atoms_in_scope,
                atom_fingerprints,
                method.matcher)

        elif method.op == "electronic":
            # -xanes etc. defer to the type pass
            # (11.3.e); no species-level effect.
            continue

    return (atom_species_id, named_regions)


function resolve_scope(scope_spec, named_regions,
        n_atoms):
    # `scope_spec`: None, "NAME", or "~NAME".
    #   None    -> every atom.
    #   "NAME"  -> atoms in named_regions[NAME].
    #   "~NAME" -> complement.
    if scope_spec is None:
        return set(range(n_atoms))
    if scope_spec.startswith("~"):
        name = scope_spec[1:]
        require(name in named_regions,
            "scope=~" + name + " references unknown"
            + " spatial region; declare it with"
            + " name=" + name + " on an earlier"
            + " -target or -block flag")
        return (set(range(n_atoms))
                - named_regions[name])
    require(scope_spec in named_regions,
        "scope=" + scope_spec + " references unknown"
        + " spatial region; declare it with name="
        + scope_spec + " on an earlier -target or"
        + " -block flag")
    return named_regions[scope_spec]


function assign_species_by_bucketing(atom_species_id,
        atoms_in_scope, atom_fingerprints, matcher):
    # Greedy single-pass bucketing.  Each in-scope
    # atom either joins an existing bucket whose
    # representative is within the matcher's
    # similarity floor, or starts a new bucket.
    # Matches the behavior of today's group_reduce.
    buckets = []     # list of Bucket(representative,
                     #                atom_indices)
    for atom_i in atoms_in_scope:
        vec = atom_fingerprints[atom_i]
        merged = False
        for b in buckets:
            if (matcher.distance(vec,
                    b.representative)
                <= matcher.default_similarity_floor):
                b.atom_indices.append(atom_i)
                # Refresh the representative so it
                # tracks the running set; cheap
                # enough at our atom counts.
                b.representative = \
                    matcher.representative([
                        atom_fingerprints[j]
                        for j in b.atom_indices])
                merged = True
                break
        if not merged:
            buckets.append(Bucket(
                representative = vec,
                atom_indices   = [atom_i]))

    # Assign a fresh species ID to each bucket.  Atoms
    # outside `atoms_in_scope` keep whatever species ID
    # earlier flags produced.  Mutates atom_species_id in
    # place; the caller needs no return value (the pick no
    # longer gates on which species were env-grouped, C93).
    next_id = max(atom_species_id) + 1
    for b in buckets:
        for i in b.atom_indices:
            atom_species_id[i] = next_id
        next_id += 1
```

---

#### 11.3.d Manifest-Entry Pick per Species (DESIGN 5.6.5)

For each `(element, species)` pair, chooses exactly
one `PotentialEntry` from the element's database.
Precedence: `-pot LABEL` manual override; then a single
best-effort fingerprint match (disabled by
`-nofingerprint`); then the default-tag fallback, the one
point that always succeeds (rule 7).  The match (DESIGN
5.6.5 step 2, the C93 model) fixes exactly one descriptor
family and one `sub_spec` -- the user's when they grouped
with `-reduce`, otherwise the database's `preferred` record
-- computes one query, and accepts a miss.  It never
searches across families or sub_specs, and it never aborts:
a database that lacks the chosen `(method, sub_spec)` is a
silent fall-through to the default, the database never
overruling the user.

```
function pickManifestEntry(species_atoms, element, db,
        pot_override, fingerprinting_enabled,
        building_loen_input,
        user_scheme, structure, grouping_descriptors):
    # species_atoms          -- atom indices in this
    #                           (element, species) bucket
    # db                     -- the element's ElementDatabase
    # pot_override           -- the -pot LABEL value, or None
    # fingerprinting_enabled -- False under -nofingerprint
    # building_loen_input    -- True when makeinput is building
    #                           the loen input itself (a
    #                           -loeninput run, DESIGN 5.10.2);
    #                           skips the match so the build
    #                           cannot invoke itself (below)
    # user_scheme            -- the active matcher object the
    #                           user grouped with (carries
    #                           active_sub_spec), or None.  In
    #                           makeinput this is the reduce
    #                           scheme; a -bispec run is
    #                           grouped upstream by makegroups
    #                           (11.3.f) and arrives as
    #                           pre-assigned species -> the
    #                           file-dictated branch below.
    # structure              -- for computing a query on demand
    # grouping_descriptors   -- per-atom vectors from the
    #                           grouping pass (reused in the
    #                           user-scheme regime), or None

    # Precedence 1: manual override.  KeyError is fatal;
    # -pot is a deliberate choice and a silent fallback
    # would mask the intent.
    if pot_override is not None:
        try:
            return lookup(db, pot_override)
        except KeyError:
            error("-pot " + pot_override + " not found in "
                + "share/atomicPDB/" + lower(element)
                + "/s_gaussian_pot.toml; manual override"
                + " must match an existing label")

    # -nofingerprint (and the reduced flow, 11.3.0): skip
    # the match entirely.
    if not fingerprinting_enabled:
        return default_entry(db)

    # A loen-descriptor build skips the match too (DESIGN
    # 5.6.5 step 2 / 5.10.2).  In the file-dictated regime the
    # match may run a loen descriptor computation, whose own
    # first step is a -loeninput build exactly like this one --
    # so matching here would invoke the build within itself,
    # without end.  Take the default entry: the bispectrum is
    # geometric, so the potential is irrelevant to the
    # descriptor this build feeds.
    if building_loen_input:
        return default_entry(db)

    # Precedence 2: a single best-effort fingerprint match.
    # Fix exactly one (matcher, sub_spec) and one query,
    # chosen by regime (DESIGN 5.6.5 step 2).
    if user_scheme is not None:
        # The user grouped with an environment scheme: honor
        # it.  Match that family at the user's sub_spec,
        # reusing the per-atom descriptors grouping computed.
        matcher  = user_scheme
        sub_spec = user_scheme.active_sub_spec
        per_atom = grouping_descriptors
    else:
        # File-dictated species (crystalline / pre-assigned):
        # the database decides via its preferred records --
        # bispectrum if it has a preferred bispectrum record,
        # else reduce.  One family only, no cascade.
        pref = find_preferred(db, "bispectrum")
        if pref is None:
            pref = find_preferred(db, "reduce")
        if pref is None:
            return default_entry(db)   # nothing to match
        matcher  = MATCHERS[pref.method]
        sub_spec = pref.sub_spec
        # Compute the one query this family needs: reduce
        # in-process, bispectrum via the loen seam (reuse a
        # grouping/loen run if one exists, else a fast loen
        # run for the pick; 11.3.a, 11.3.f).
        if matcher.needs_loen_run:
            per_atom = loen_descriptors(structure, matcher,
                sub_spec)
        else:
            per_atom = matcher.compute_query(structure,
                sub_spec)

    rep = matcher.representative(
        [per_atom[i] for i in species_atoms])

    # Shared match: nearest entry fingerprint at (method,
    # sub_spec), accepted only within the similarity floor.
    best_entry    = None
    best_distance = +infinity
    for entry in db.potentials:
        try:
            fp = find_fingerprint(entry, matcher.name,
                sub_spec)
        except KeyError:
            continue      # this entry has no record at the
                          # chosen (method, sub_spec); skip
        # extract_query_vector reads the matcher-specific
        # payload field (bispec `values`, reduce `shell_code`)
        # so this stays descriptor-agnostic.
        d = matcher.distance(rep,
            matcher.extract_query_vector(fp.payload))
        if d < best_distance:
            best_distance = d
            best_entry    = entry

    if (best_entry is not None
            and best_distance
                <= matcher.default_similarity_floor):
        return best_entry
    if best_entry is not None:
        # A near miss: records exist but none is close enough.
        warn("species in " + element + " best fingerprint"
            + " match " + best_entry.label + " at distance "
            + str(best_distance) + " (> floor "
            + str(matcher.default_similarity_floor)
            + "); using the default tag")
    # best_entry is None -> no comparable record at all (e.g.
    # the user ran at a sub_spec the database lacks): a silent
    # best-effort miss, the database never overruling the user.

    # Precedence 3: default tag.  Guaranteed by rule 7.
    return default_entry(db)


function find_preferred(db, method):
    # Return A FingerprintRecord flagged preferred = true for
    # `method` in this element's database, or None if the family
    # is absent from the file entirely.
    #
    # Several entries flag their own canonical record (rule 10 is
    # per entry), so there are as many hits as there are
    # harvested entries.  Any of them will do, and the first is
    # the cheapest: the caller reads only `method` and `sub_spec`
    # off the result -- never the payload -- to learn which
    # settings to compute its query with, and rule 10's file-wide
    # half guarantees every preferred record of a method names
    # the SAME sub_spec.  A None result therefore means "this
    # family is absent," never "present but unpreferred."
    for entry in db.potentials:
        for fp in entry.fingerprints:
            if fp.method == method and fp.preferred:
                return fp
    return None
```

---

#### 11.3.e Type Pass (DESIGN 5.6.6)

Types are subdivisions of a species made on electronic
grounds, not geometric grounds.  Every species starts
with one inherited type; electronic-state flags split
off new types whose potentials come from existing
machinery (the core-hole potential for `-xanes` today).
From Imago's perspective the type pass produces the
flattened per-type list emitted in 11.3.g.

```
function typePass(structure, atom_species_id,
        species_potentials, settings):
    # Start with one type per species: every atom in
    # species S becomes (S, 1) and inherits the
    # species' chosen PotentialEntry.
    n_atoms          = len(structure.atoms)
    atom_type_id     = [(s, 1)
                        for s in atom_species_id]
    type_potential   = {}    # (species, type)
                             #   -> PotentialEntry
    for s, entry in species_potentials.items():
        type_potential[(s, 1)] = entry

    # Apply electronic-state flags in CLI order.
    # Today only -xanes splits types; future flags
    # layer in the same way without touching the
    # matcher protocol or the species pass.
    for method in settings.methods:
        if method.op != "electronic":
            continue
        if method.name == "xanes":
            apply_xanes_type_split(structure,
                atom_species_id, atom_type_id,
                type_potential, method)
        else:
            error("unknown electronic-state flag: "
                + method.name)

    return atom_type_id, type_potential


function apply_xanes_type_split(structure,
        atom_species_id, atom_type_id,
        type_potential, method):
    # The XANES core-hole atom and its in-sphere
    # neighbors are split off into new types within
    # their parent species.  This function handles
    # only the species-to-type bookkeeping; the
    # potentials themselves come from the existing
    # XANES core-hole machinery, unchanged.
    affected = compute_xanes_affected_atoms(structure,
        method)
    species_to_next_type = {}    # species
                                 #   -> next type id

    for atom_i in affected.core_hole:
        s = atom_species_id[atom_i]
        t = next_type_id(species_to_next_type, s)
        atom_type_id[atom_i] = (s, t)
        type_potential[(s, t)] = \
            build_xanes_core_hole_potential(
                structure, atom_i, method)

    for atom_i in affected.in_sphere_neighbors:
        s = atom_species_id[atom_i]
        t = next_type_id(species_to_next_type, s)
        atom_type_id[atom_i] = (s, t)
        type_potential[(s, t)] = \
            build_xanes_neighbor_potential(
                structure, atom_i, method)


function next_type_id(species_to_next_type, s):
    # Type 1 is always the inherited type from
    # 11.3.d.  Each subsequent split takes the next
    # unused id within the parent species.
    current = species_to_next_type.get(s, 1)
    species_to_next_type[s] = current + 1
    return current + 1
```

---

#### 11.3.f makegroups.py: bispectrum grouping (DESIGN 5.10)

A Fortran-side descriptor can only come from a completed
Imago run, so bispectrum grouping is a *sequence* run from
*outside* makeinput by `makegroups.py` -- never by
makeinput re-invoking itself.  `makegroups` is dual-mode:
an importable `group_by_bispectrum` the producer
(`build_initial_potentials.py`) calls, plus a `__main__`
CLI for manual use.  It runs the loen flow and rewrites the
skeleton with explicit per-element species tags; makeinput
then reads those tags like any other explicit assignment.
Grouping applies only to non-crystalline (P1) systems; a
symmetry-bearing skeleton is refused up front, since
rewriting its types in P1 would drop the space group that
k-point folding depends on (DESIGN 5.10.1).

```
function group_by_bispectrum(skeleton_path, sub_spec,
        similarity_floor):
    matcher = MATCHERS["bispectrum"]()

    # 1. P1 guard.  Grouping rewrites every atom's type
    #    from its fingerprint, which only makes sense for a
    #    non-crystalline cell.  Read the skeleton and refuse
    #    unless its space group resolves to number 1 (the
    #    `space` line reads `1_a`) and its supercell is
    #    1 1 1.  Regrouping a real crystal would force it
    #    into P1 and discard the space group the
    #    Brillouin-zone k-point folding relies on, so a
    #    symmetry-bearing skeleton is rejected here rather
    #    than silently corrupted; a crystal is harvested on
    #    the witness path instead (DESIGN 5.10.1, 5.10.4).
    structure = read_skeleton(skeleton_path)
    require(structure.space_group_num == 1
            and structure.supercell == [1, 1, 1],
        skeleton_path + " is not P1 (space group "
        + structure.space_group + ", supercell "
        + str(structure.supercell) + "); bispectrum "
        + "grouping rewrites types in P1 and would drop "
        + "the space group k-point folding needs.  Group "
        + "only non-crystalline (P1, `1_a`) skeletons; "
        + "harvest a crystal on the witness path instead.")

    # 2. First makeinput: a provisional imago.dat with no
    #    grouping.  The LOEN_INPUT_DATA block carries the
    #    sub_spec via matcher.to_loen_input.  This is a
    #    -loeninput build, so the fingerprint match is skipped
    #    (building_loen_input, pickManifestEntry above) and each
    #    atom takes the default entry -- the potential is
    #    irrelevant (bispectrum is geometric), and the skip is
    #    what keeps this build from invoking itself.
    run_makeinput(skeleton_path,
        loen_params = matcher.to_loen_input(sub_spec))

    # 3. Run loen.  -scf no skips the SCF; loen needs only
    #    the structure and the LOEN block.  Produces a
    #    self-describing fort.21 (DESIGN 5.10.3).
    run_imago(flags = ["-loen", "-scf", "no"])

    # 4. Read fort.21.  Each row carries its own identity
    #    (site#, element, species, type_in_species,
    #    type_flat) and the bispectrum vector, so the
    #    row -> atom mapping is read off the file -- no
    #    separate datSkl.map lookup (DESIGN 5.10.3).
    rows = matcher.parse_loen_output("fort.21", sub_spec)

    # 5. Bucket atoms by fingerprint distance within the
    #    floor, per element, refreshing each bucket's
    #    representative as it grows (the same bucketing as
    #    11.3.c, but run here in the orchestrator rather
    #    than inside makeinput).
    species_of = bucketByFingerprint(rows, matcher,
        similarity_floor)

    # 6. Rewrite the skeleton with explicit per-element
    #    species tags -- Si1,Si2,...,O1,O2,... restarting
    #    at 1 for each element (DESIGN 5.10.4).  A
    #    round-trip test guards the numbering.
    write_skeleton_with_species(skeleton_path, species_of)
    return species_of
```

The producer then runs makeinput on the rewritten skeleton
(now explicitly typed) and proceeds to SCF and harvest.
There is no recursion to guard against: each step is an
ordinary process the orchestrator runs in order.

---

#### 11.3.g Driver (DESIGN 5.6.7)

Top-level orchestrator.  Chains preflight, species
pass, manifest-entry pick, type pass, and emit (no
bootstrap step -- bispectrum atoms arrive pre-grouped
from makegroups, 11.3.f).  The driver is
matcher-agnostic; all descriptor-family knowledge lives
in the matcher classes (11.3.a) and the `MATCHERS`
registry.

```
function emitInitialPotentials(structure, settings,
        imago_input):
    # 1. Identify the in-makeinput environment scheme
    #    the user grouped with (reduce), if any.  DESIGN
    #    5.6.2 mutual exclusion is enforced at argparse
    #    time, so at most one is active -- often none, when
    #    species are crystallographic, spatial, or
    #    pre-assigned (bispectrum, grouped upstream by
    #    makegroups, 11.3.f).  `user_scheme` is None in
    #    that file-dictated case.
    user_scheme = first_environment_matcher(
        settings.methods)

    # 2. Whether the fingerprint pick runs at all.  On by
    #    default (the C93 decoupling); off under
    #    -nofingerprint.
    fingerprinting_enabled = not settings.nofingerprint

    # 3. Per-element preflight.  Loads each element's
    #    database, emits the (info, never fatal) coverage
    #    note when fingerprinting is enabled, and marks
    #    elements without a database for the legacy fallback.
    databases = perElementPreflight(structure,
        fingerprinting_enabled)

    # 4. Compute per-atom grouping descriptors once when the
    #    user grouped with an environment scheme, so the
    #    species-pass bucketing and the pick's representative
    #    see the same vectors (the pick reuses them).  Only
    #    Python-side matchers (reduce) group here;
    #    bispectrum grouping is done ahead of makeinput by
    #    makegroups (11.3.f), so those atoms arrive typed and
    #    take the file-dictated branch in the pick.
    grouping_descriptors = None
    if user_scheme is not None:
        grouping_descriptors = \
            user_scheme.compute_query(structure,
                user_scheme.active_sub_spec)

    # 5. Species pass.  Position-based and environment-based
    #    flags compose in CLI order; output is a per-atom
    #    species ID array and the dict of named regions used
    #    in scope resolution.  No env-only set -- the pick
    #    runs for every species (C93).
    atom_species_id, named_regions = speciesPass(
        structure, settings, databases,
        grouping_descriptors)

    # 6. Manifest-entry pick per species.  Atoms in
    #    elements with no database file
    #    (databases[elem] is None) bypass the
    #    manifest machinery and emit via the legacy
    #    pot1/coeff1 path in step 8.
    species_potentials = {}
    for species in unique_species(atom_species_id):
        atoms   = atoms_of_species(species,
            atom_species_id)
        element = structure.atoms[atoms[0]].element
        if databases[element] is None:
            species_potentials[species] = \
                LegacyEntry(element)    # marker
                                        # consumed in
                                        # step 8
            continue
        species_potentials[species] = \
            pickManifestEntry(
                species_atoms        = atoms,
                element              = element,
                db                   = databases[element],
                pot_override         = settings
                                        .pot_override,
                fingerprinting_enabled =
                    fingerprinting_enabled,
                building_loen_input  =
                    (settings.loeninput is not None),
                user_scheme          = user_scheme,
                structure            = structure,
                grouping_descriptors =
                    grouping_descriptors)

    # 7. Type pass.  Inherit from the species pass
    #    and apply electronic-state flags (XANES
    #    today; future flags layer in unchanged).
    atom_type_id, type_potential = typePass(
        structure, atom_species_id,
        species_potentials, settings)

    # 8. Emit per-Imago-type blocks in today's on-
    #    the-wire format.  Imago is unaware of the
    #    manifest or the matcher; it sees only the
    #    resolved per-type numbers.
    for (species, type_id), entry in \
            sorted_type_iter(type_potential):
        element = element_of_species(species,
            atom_species_id, structure)
        if isinstance(entry, LegacyEntry):
            emitLegacyElementBlock(imago_input,
                entry.element)
        else:
            db = databases[element]
            emitElementBlock(imago_input, db, entry)


function emitElementBlock(out, db, entry):
    # Order matches today's on-the-wire Imago input
    # format.  Field names below are TOML-side; the
    # Imago input writer keeps its existing tags.
    write_pot_block(out,
        nuclear_z       = db.nuclear_z,
        nuclear_alpha   = db.nuclear_alpha,
        covalent_radius = db.covalent_radius,
        num_gaussians   = entry.num_gaussians,
        alpha_min       = entry.alpha_min,
        alpha_max       = entry.alpha_max)
    write_coeff_block(out,
        coefficients = entry.coefficients,
        alphas       = entry.alphas)


function first_environment_matcher(methods):
    # DESIGN 5.6.2: at most one environment-based
    # scheme per run; argparse rejects multiples.
    # The matcher object on each method record was
    # built at parse time with its active_sub_spec
    # already bound from the CLI.
    for m in methods:
        if m.op == "environment":
            return m.matcher
    return None
```

---

### 11.4 Build Pipeline (DESIGN 5.7, schema v2)

Incremental build of the affected
`s_gaussian_pot.toml` files (see DESIGN 5.7 for the
layered reproducibility contract: bit-level emitter,
precision-level numerics, free metadata).  The
producer is a **kaleidoscope client**: it runs no SCF
itself.  It works in three phases -- *build*,
*converge*, *harvest*.  The build phase loads each
existing database and refreshes its "isolated" entry
from current `pot1`/`coeff1` files (so atomSCF changes
propagate), then, per reference solid, materializes the
structure, asks the guidance predictor for a **seed
density** (`predict_kpoint_density`, 4e.7 -- prediction
only, no grid) and builds the solid's `ClimbConfig` from
that prediction's confidence and the run's climb policy.
It also emits one structure-only `imago.py -loen -scf no`
unit per Fortran-side fingerprint; these are geometry-only,
so they dispatch once in a small pre-flight batch before
the climb.  The converge phase drives every solid through
the adaptive mesh climb (`converge_by_climb`, 4e.5): serial
within a solid, concurrent across, each solid climbing on
the instant its own rung lands so no solid waits on another,
each rung's mesh (`build_mesh_unit`, 4e.7) launched as the
climb chooses it, until each solid's energy flattens (its
converged rung) or hits the `max_count` ceiling
(NON_CONVERGED).  The
harvest phase, for each converged solid, locates the
converged rung's `kpt-mesh-<a>-<b>-<c>` run, then for each
distinct environment the run discovered takes one
representative, extracts its potential, computes the
database-wide `[characterization]` fingerprints (plus any
rare per-entry override), and insert-or-skips the result
into the per-element database (DESIGN 5.2.3).  It then
contributes the same converged rung back to the guidance
dataspace **in memory** -- `record_converged` (4e.6) feeds
the shared `build_entry` (15.7) -- with no re-read of the
workspace.

**The injected-step convention.**  Five of the pipeline's steps
are taken as *parameters* rather than called by name:
`prepare_fn`, `dispatch_fn`, `extract_fn`, `identity_fn`, and
`fingerprint_fn`.  Each defaults to the real thing -- the
driver-side prepare pass, kaleidoscope's dispatch, the `scfV`
potential reader, the `datSkl.map` site-identity reader, and the
fingerprint harvest -- so a production run names none of them and
reads exactly as the phases above describe.

They exist because every one of those five steps needs a live
Imago run underneath it, and the producer's *orchestration* is
what the tests need to exercise: whether the manifest resolves,
whether one flight is built, whether a converged rung is picked,
whether the insert-or-skip is idempotent.  Injecting the five
lets that orchestration be tested end to end with the toolchain
mocked, and keeps the seam explicit at the signature rather than
hidden behind patched module globals.  The convention is named
once here so a reader meets it once rather than five times; the
pseudocode below calls each step by its plain name.

The v2 manifest reader (`load_manifest_v2` below)
enforces the validation rules of DESIGN 5.7.
That reader, the relaxed structure-only reader
(`load_structure_sources`), and the writer
(`format_manifest`) live in the shared schema
library `curation_manifest.py`, imported by both the
producer and the `expand_manifest.py` authoring tool
that writes manifests from `cod_fish.py` sketches.
Producer-side fingerprint
harvest splits Python-side matchers (in-process,
descriptor-agnostic) from Fortran-side matchers, which
read the `fort.21` of the `-loen` unit kaleidoscope
already dispatched (no separate loen cache; the
kaleidoscope run-reuse cache of DESIGN 6.2.5 subsumes
it).

```
function prepare_units(flight, workspace, units = None):
    # The driver-side prepare pass (DESIGN 6.2.5).  The cache
    # keys on makeinput's OUTPUT, so makeinput must run before
    # the hit-test can read it -- and it runs HERE, in the
    # driver, which is what lets a hit be decided from local
    # files and never reach the scheduler.
    #
    # `units` selects which of the flight's units to prepare.
    # It defaults to every unit (a one-shot flight); the climb
    # passes only the rungs it has newly decided (4e.7), so an
    # accreting flight re-prepares nothing it already staged.
    targets = flight.units if units is None else units
    for unit in targets:
        # Each unit is staged into its own area, deliberately
        # SEPARATE from its run directory: the hit-test is about
        # to byte-compare against the PRIOR run's staged copy,
        # so building over that copy first would destroy the
        # very reference the comparison needs (the "must not
        # clobber" rule).
        staging = join(workspace, "prepare", unit.id, *unit.calc)
        makeinput_options, _ = partition_options(unit.options)
        build_run_dir(unit.structure, makeinput_options, staging)
        unit.prepared_dir = staging

        # Re-point each KeyFile's source at the copy just built.
        # standard_key_fields left it provisional because only
        # this pass knows where the unit was staged.
        #
        # The join is the unit's directory plus the declared
        # PATH, and it is the same join the hit-test makes on the
        # run directory (13.4).  That symmetry is the point: the
        # declared path already carries the `inputs/` level, so
        # this pass adds nothing to it and has no layout of its
        # own to keep in step (DESIGN 6.2.5).
        for key_file in unit.key_fields.files:
            key_file.source = join(staging, key_file.path)
```

```
function buildInitialPotentials(manifest_path,
        force, single_element, dispatch_shape,
        partition, nodes, walltime, profile,
        save_config, clean_after, tidy_run):
    manifest     = load_manifest_v2(manifest_path)
    # Fill each sparse solid from the shared blocks before the
    #   pipeline reads it (11.4): run settings from [defaults],
    #   the k-point flatness tolerance from [harvest] or its
    #   built-in default (DESIGN 5.7).  After this pass every
    #   setting is populated -- harvest (Phase 3) reads the
    #   resolved kpoint_convergence_threshold per solid.
    apply_manifest_defaults(manifest)
    # Resolve the climb policy ONCE for the run (4e.4 / DESIGN
    #   3.12.6): the manifest's optional [harvest.kpoint_climb] knobs
    #   merged over the provisional defaults into the confidence
    #   `thresholds` bundle and the per-axis `max_count` ceiling that
    #   every solid's ClimbConfig reads.  A mistyped knob already
    #   failed loudly at load (load_manifest_v2 validates the keys).
    thresholds, max_count = climb_policy_from_manifest(
        manifest.harvest.kpoint_climb)
    dataspace    = guidance_db.load(
        "share/historicalGuidanceDB/")
    imago_commit = git_sha("HEAD")
    timestamp    = iso8601_now_utc()
    workspace    = curation_workspace_root()

    # The producer's per-run bookkeeping, hung on every unit it
    #   dispatches (DESIGN 6.2.4).  It is what the reuse plan prints
    #   behind a reused result (13.5) and what a guidance entry's
    #   provenance records (15.7), and it is never compared: a
    #   rebuilt engine does not invalidate a stored starting
    #   potential (DESIGN 6.2.5).  One value for the whole build --
    #   it describes this run of the producer, not any one solid --
    #   and deliberately NOT a member of any solid's `options`,
    #   because it reaches neither makeinput nor imago.
    record = {"imago_commit": imago_commit}

    # ===== Phase 1: build =============================
    # Step 1a: refresh "isolated" entries.  atomSCF
    # changes propagate every run.  The isolated
    # entry's `default` flag is computed at build
    # time per is_isolated_default_for() below: true
    # iff the manifest contributes no other entry
    # for this element.
    elements = list_dirs("share/atomicPDB/")
    if single_element is not None:
        elements = [single_element]

    databases = {}      # elem -> ElementDatabase
    for elem in elements:
        path = element_path(elem)
        if file_exists(path):
            db = load(path,
                known_methods = MATCHERS.keys())
        else:
            pot1 = read_pot1(elem)
            db = ElementDatabase(
                schema_version  = 2,
                element_symbol  = elem,
                nuclear_z       = pot1.nuclear_z,
                nuclear_alpha   = pot1.nuclear_alpha,
                covalent_radius = pot1.covalent_radius,
                potentials      = [])
        # INCREMENTAL (DESIGN 5.7): the existing database was
        # loaded above (or started empty), so every environment
        # harvested by earlier runs is preserved -- there is no
        # reset.  Refresh only the "isolated" baseline from the
        # current pot1/coeff1 (so atomSCF changes propagate):
        # drop the old isolated entry by label and append the
        # fresh one, leaving all harvested entries untouched.
        # The harvest phase (below) then inserts-or-skips this
        # run's solids on top; re-running an unchanged manifest
        # moves nothing, because every duplicate is skipped
        # (5.2.3), so no count can inflate.
        iso = build_isolated_entry(elem, imago_commit,
            timestamp, manifest)
        db.potentials = [p for p in db.potentials
                         if p.label != "isolated"]
        db.potentials.append(iso)
        databases[elem] = db

    # Step 1b: per reference solid, materialize the structure,
    # PREDICT its seed density (no grid), and gather everything the
    # climb needs for it.  The CONVERGENCE units are NOT built here
    # -- the climb builds each round's meshes as it runs (Phase 2).
    # Only the geometry-only fingerprint units are built now, all
    # collected into one pre-flight batch (they are mesh-independent,
    # so they belong to no climb round -- DESIGN 5.7 / Q3).
    struct_of      = {}   # reference_id -> local struct path
    options_of     = {}   # reference_id -> tool-facing options
    configs        = {}   # reference_id -> ClimbConfig
    seed_densities = {}   # reference_id -> round-0 seed density
    predictions    = {}   # reference_id -> PredictionRecord dict
    loen_units     = []   # every solid's fingerprint units, one batch
    for ref in manifest.reference_solids:
        struct = materialize_structure(ref)
        struct_of[ref.reference_id] = struct

        # Fixed RUN SETTINGS in each tool's coded vocabulary
        # (DESIGN 6.2.10): make_producer_options maps the manifest
        # physics -- functional -> xccode, kpoint_integration ->
        # scfkpint (a "gaussian-0.1" width -> thermsmear /
        # THERMAL_SMEARING_SIGMA), basis -> scf_basis, scf_threshold
        # -> converg, shift -> kpshift.  EVERY key here is a real
        # tool input; the build identity is not one and travels on
        # unit.record instead (DESIGN 6.2.10, and see `record` just
        # below).  Kept for the whole climb: every round's mesh unit
        # copies them (build_mesh_unit, 4e.7).
        options = make_producer_options(ref)
        options_of[ref.reference_id] = options

        # The predictor and the PredictionRecord speak the human
        # physics names, not the codes, so the sub-model travels in
        # its OWN dict -- never mixed into the tool-facing options
        # (DESIGN 6.2.8 / 6.2.10), which would both duplicate the
        # basis and make makeinput reject "functional" /
        # "kpoint_integration".
        submodel = {"basis":              ref.basis,
                    "functional":         ref.functional,
                    "kpoint_integration": ref.kpoint_integration}

        # PREDICTION ONLY (4e.7): signature -> predict -> record.
        # It lays NO grid; the climb seeds from the density and
        # picks its mode and persistence from the confidence
        # (3.12.4 / 3.12.6).  `center` carries a curator-pinned
        # density (bypasses the predictor); otherwise None and
        # predict runs, returning is_under_trained when the
        # dataspace has no useful prior (7.9).
        # The two resolved knobs the harvest cannot look up travel
        #   with the prediction (15.6): this solid's per-atom
        #   flatness tolerance, and the database-wide metal gap cut
        #   the climb reads as config.metal_gap_threshold (4e.3), so
        #   the climb's short-circuit and the harvest's metal skip
        #   apply the SAME resolved number.  Named `prediction`
        #   here, since `record` above is the unit bookkeeping.
        harvest_thresholds = {
            "kpoint_convergence_threshold":
                ref.kpoint_convergence_threshold,
            "metal_gap_threshold": thresholds.metal_gap_threshold}
        density, confidence, under_trained, prediction = \
            predict_kpoint_density(
                struct, dataspace, ref.system_type,
                submodel, center = ref.kpoint_spec.density,
                harvest_thresholds = harvest_thresholds)
        seed_densities[ref.reference_id] = density

        # Everything the climb needs for THIS solid, gathered once
        # (build_climb_config below): the reciprocal-cell geometry
        # the rung mechanics read, the confidence-derived mode /
        # persistence / grid from the resolved policy, and the
        # energy / ceiling knobs.
        configs[ref.reference_id] = build_climb_config(
            ref, struct, confidence, under_trained,
            thresholds, max_count)

        # Store it as a plain dict (metadata must be TOML-
        # serializable).  Both resolved thresholds are already
        # fields of the record itself (15.6) rather than stamped on
        # afterwards, so there is exactly one list of what a
        # prediction record carries and the harvest's reads cannot
        # outrun it (DESIGN 7.8 / 5.7).
        predictions[ref.reference_id] = as_dict(prediction)

        # Geometry-only fingerprint units: one structure-only
        # `-loen -scf no` unit per Fortran-side declaration, tagged
        # kind = "fingerprint" (DESIGN 6.2.9).  The bispectrum
        # fingerprint depends on geometry alone, so these need not
        # wait for a converged mesh; they dispatch in the pre-flight
        # below, and their run dirs persist for the harvest.
        loen_units.extend(build_loen_units(
            ref, struct, options, record,
            manifest.characterization))

    # Fail fast, before anything is dispatched (DESIGN 5.10.6).
    # Every Fortran-side declaration the Phase 3 harvest could
    # read must already have a unit above.  This can only fire
    # when something upstream is already broken -- which is
    # precisely why it is worth asserting here, where it costs a
    # set comparison, rather than discovering it in Phase 3 after
    # the whole convergence sweep has been paid for.
    assertLoenCoverage(loen_units, manifest.reference_solids,
                       manifest.characterization)

    # ===== Phase 1b: loen pre-flight (DESIGN 5.7 / Q3) ==
    # The fingerprint units are geometry-only and mesh-independent,
    # so they belong to no climb round: dispatch them ONCE as a
    # small flat batch, and their run dirs persist for the Phase 3
    # fingerprint harvest.  makeinput runs HERE in the driver
    # (DESIGN 6.2.5): the cache keys on structure.dat, so it must
    # exist before the hit-test; prepare_units stages each unit
    # into its own `prepare/` area, separate from its run dir (the
    # "must not clobber" rule).  The CONVERGENCE units are prepared
    # as the climb dispatcher sends each rung instead (4e.7), not
    # here.
    #
    # Resolve the dispatch Config once for the whole run (13.7):
    # local -> None (the driver runs in process); a cluster shape ->
    # a real Config.  Both the pre-flight and every climb round
    # dispatch under it; `force` bypasses the run-reuse cache and
    # fresh results still repopulate it (6.2.5).
    parsl_config, choices = resolve_dispatch(
        dispatch_shape, partition, nodes, walltime, profile)
    if save_config and choices is not None:
        write_resolved_dispatch(workspace, choices, profile)

    # One executor for the whole run's dispatch (DESIGN 6.2.11's
    # pooled shape, 13.5): the pre-flight and every climb round
    # share one warm pool.  Build it ONCE and close it ONCE, in a
    # finally so a mid-climb error still releases the SLURM
    # allocation.  Phase 3 harvest reads run dirs only and needs no
    # dispatch, so the pool is freed the moment convergence ends.
    executor = make_executor(parsl_config)            # 13.5

    # Layer (b), prune-as-you-go (DESIGN 6.2.12).  When --tidy-run
    # is asked for, each flight is handed a callback that kaleido-
    # scope fires as each unit reaches a terminal state, and that
    # callback prunes that unit's scratch there and then, while the
    # rest of the campaign is still running.  Nothing is added to
    # kaleidoscope for this: on_outcome already reports a landing
    # and already carries the run directory.
    #
    # A campaign that is not asked to tidy passes None and behaves
    # exactly as before.  $IMAGO_TEMP unset means the check that
    # scratch lies where it should cannot run, so the prune is
    # declined out loud rather than performed unchecked -- the same
    # rule layer (a) applies after the harvest.
    #
    # Prune failures are collected across BOTH flights into one
    # list and reported again at the very end (below).  A failure
    # cannot be allowed to stop the campaign, but it must not be
    # able to scroll away either: it says a filesystem assumption
    # of 6.2.12 has broken, and that is the user's to act on.
    scratch_root = env("IMAGO_TEMP")
    prune_enabled = tidy_run and scratch_root != ""
    prune_problems = []
    if tidy_run and not prune_enabled:
        report("--tidy-run skipped: $IMAGO_TEMP is unset")

    try:
        # The flight is built whether or not there is anything to
        # dispatch: Phase 3's fingerprint harvest reads it below.
        loen_flight = Flight(
            units        = loen_units,
            root         = workspace,
            parsl_config = parsl_config,
            sweep        = SweepRecord(
                varied_axes = (), fixed_axes = {}))
        if prune_enabled:
            loen_flight.on_outcome = make_prune_callback(
                loen_flight, workspace, scratch_root,
                prune_problems)
        if loen_units is not empty:
            prepare_units(loen_flight)    # driver-side makeinput
            dispatch(loen_flight, executor = executor,
                     force = force)

        # ===== Phase 2: converge ======================
        # The climb dispatcher (4e.7) closes over each solid's
        # structure, its coded options, and the shared executor,
        # and owns one flight whose unit list accretes as rungs are
        # decided.  Its send(mesh_lists) builds one mesh unit per
        # requested mesh (build_mesh_unit), prepares the new units
        # driver-side, and launches them under that one executor
        # WITHOUT waiting (send_off); its next_rung() blocks until
        # the next rung lands (collect_next) and reads its
        # result.toml back into a Rung -- returning the FAILED
        # marker for any mesh whose run did not complete (7.7).
        # converge_by_climb (4e.5) drives every solid to a verdict,
        # serial within a solid and concurrent across, judging each
        # solid the moment its own rung lands.  A mesh re-run later
        # is a cache hit (6.2.5).  on_non_converged tags the solid's
        # workspace with a prediction mismatch (7.8 3d); it is
        # injected so the climb loop stays free of the workspace.
        #
        # The climb's flight gets the same prune callback.  This is
        # where prune-as-you-go earns its keep: the pre-flight is a
        # handful of geometry-only runs, while the climb dispatches
        # many rungs per solid and is the phase long enough to
        # exhaust scratch before it ends.
        dispatcher = make_climb_dispatcher(
            struct_of, options_of, workspace,
            parsl_config = parsl_config, executor = executor,
            force = force, record = record,
            tidy_run = prune_enabled,
            scratch_root = scratch_root,
            prune_problems = prune_problems)
        materials = [ref.reference_id
                     for ref in manifest.reference_solids]
        outcomes, rungs, verdicts = converge_by_climb(
            materials, configs, seed_densities, dispatcher,
            on_non_converged =
                lambda m: tag_prediction_mismatch(workspace, m))
    finally:
        executor.close()      # release the pool (also on error)

    # ===== Phase 3: harvest ===========================
    log = []
    for ref in manifest.reference_solids:
        struct  = struct_of[ref.reference_id]
        outcome = outcomes[ref.reference_id]
        if outcome is NON_CONVERGED:
            # A ceiling stop or a run failure (4e.5): flag it,
            # harvest neither a potential nor a guidance entry.
            log.append(make_nonconverged_log_entry(ref))
            continue

        # record_converged (4e.6) turns this solid's converged rung
        # and its ladder into the guidance density, exact mesh, and
        # flatness trace IN MEMORY.  Both the run log and the guidance
        # entry read it, so it is computed once, here.
        harvest = record_converged(outcome, rungs[ref.reference_id],
                                   configs[ref.reference_id])

        # The converged rung names the mesh whose run carries the
        # converged potential (Q4).  Rebuild the SAME mesh unit the
        # climb dispatched (build_mesh_unit, 4e.7) to point at that
        # run dir -- a cache hit (6.2.5), so it costs nothing -- and
        # read its result.toml for the SCF iteration count and the
        # measured character the guidance entry needs.
        converged = build_mesh_unit(
            struct, options_of[ref.reference_id], outcome.mesh,
            id = ref.reference_id, record = record)
        converged_result = read_result_toml(workspace, converged)
        # The run log records the settled mesh AND its k-density
        # (both from record_converged), the SCF iteration count, and
        # the climb's VERDICT verbatim.  The row's `converged` flag
        # is derived from the verdict rather than asserted: it means
        # k-point converged, so it is FALSE for a metal even though
        # the row names a mesh and a potential WAS harvested (DESIGN
        # 5.7).  Asserting it true here is what made a metal
        # indistinguishable from an insulator on disk.
        log.append(make_run_log_entry(
            ref, harvest, converged_result,
            verdicts[ref.reference_id]))

        # Harvest one representative per DISTINCT ENVIRONMENT the
        # converged run discovered (DESIGN 5.7).  The grouping pass
        # gave every atom a species; the assigning method's
        # partition (symmetry for an ordered solid, bispectrum for a
        # disordered one) defines the environments.  When a curator
        # customization addresses an environment it pins the
        # representative atom_site and supplies label / default /
        # description; otherwise the representative is the order-
        # independent one (DESIGN 5.6.5) and label / description are
        # derived / auto-composed from ref.source_description.
        #
        # INTERIM (C88): environment auto-discovery is built
        # (discover_environments below partitions the run by
        # (element, species, type) and yields one representative per
        # group, with customizations layered on); the cross-run
        # insert_or_skip dedup is not, so it still degrades to
        # append-or-replace-by-label.
        for env in discover_environments(converged, ref):
            elem = env.element
            if elem not in databases:
                continue   # filtered out by --element
            # extract_potential reads the representative site's TYPE
            # block from the converged scfV (NUM_TYPES + per-type;
            # the type number comes from datSkl.map, 9.7).
            coeffs, alphas = extract_potential(
                converged, env.atom_site)
            # Every environment computes the database-wide
            # [characterization] preferred fingerprints; a rare
            # per-entry override adds extra non-preferred sub_specs.
            # Fortran-side matchers read the pre-flight loen run
            # (loen_flight, dispatched in Phase 1b); Python-side
            # matchers compute in-process (DESIGN 5.7 step 5c ii).
            fingerprints = harvestFingerprints(loen_flight,
                ref, env, struct,
                manifest.characterization)
            new = PotentialEntry(
                label         = env.label,
                default       = env.default,
                description   = env.description,
                num_gaussians = len(coeffs),
                alpha_min     = min(alphas),
                alpha_max     = max(alphas),
                # The entry stores this representative's harvested
                # potential verbatim (DESIGN 5.2.3): no mean,
                # spread, or counts -- the leaner skip-on-match
                # model keeps one representative per environment and
                # nothing more.  insert_or_skip appends it when the
                # environment is new and drops it on a match.
                coefficients    = coeffs,
                alphas          = alphas,
                # Provenance records ref.system_type for forensics
                # (DESIGN 5.7 rule 2) and the converged run's SCF
                # iteration count (from its result.toml).
                provenance    = make_imago_provenance(
                    imago_commit, timestamp,
                    ref, env.atom_site,
                    converged_result.iterations),
                fingerprints  = fingerprints)
            insert_or_skip(databases[elem], new, ref)

        # Step 3d: in-memory guidance contribution (DESIGN 5.7).
        # The harvest facts record_converged produced above (the
        # density and the flatness ladder), with the converged run's
        # result.toml (gap / magnetization / SCF threshold / exact
        # mesh / commit), feed the SHARED build_entry (15.7) -- the
        # same builder the standalone density harvest uses, so the
        # two paths cannot diverge -- and save_entry stages it.
        # build_entry reads the exact converged mesh from the
        # result.toml (15.7), the same source in both paths, so the
        # climb hands only the density and ladder here.  A climb
        # that stopped on its TWO-SIDED test carries at least the
        # three distinct rungs that test required (4e.3), so its
        # ladder is long enough for the curator to re-judge.  The
        # metal short-circuit is the exception -- it stops at the
        # first gapless rung, so its ladder can be a single point
        # and its density is no convergence claim -- and build_entry
        # returns None for it (15.7).  Not every converged solid
        # contributes an entry; every converged NON-METAL does.
        #
        # ladder_is_metal hands build_entry the climb's VERDICT
        # rather than making it re-derive the classification from the
        # one settled rung it can see (DESIGN 7.8 d').  The climb read
        # every rung; the builder reads one, and one rung's apparent
        # gap on a discrete mesh is close to a coin toss (DESIGN 1.6).
        entry = build_entry(
            workspace, struct, predictions[ref.reference_id],
            dataspace, load_structure(struct),
            ref.kpoint_convergence_threshold,
            harvest.grid_values, harvest.grid_energies,
            harvest.converged_kpoint_density,
            converged_result,
            ladder_is_metal =
                (verdicts[ref.reference_id] is METAL))
        if entry is not None:
            save_entry(entry, "share/historicalGuidanceDB/")
        else:
            log(ref.reference_id + ": metal -- potential harvested,"
                + " no guidance entry")

    # ===== Write outputs ==============================
    # All affected element files via the deterministic
    # emitter (5.5).
    for elem, db in databases.items():
        save(db, element_path(elem))

    # Run log capturing the manifest snapshot, per-run
    # iteration counts, the converged mesh and its k-density
    # per solid, and the Imago commit.  The validation
    # harness (11.5) reads this log.
    write_run_log(
        "share/curation/run_log.toml",
        manifest_snapshot = manifest,
        imago_commit      = imago_commit,
        per_run_log       = log)

    # Layer (a), the post-harvest sweep (DESIGN 6.2.12).  Placed
    # deliberately AFTER the databases and the run log are on disk:
    # everything harvested is then safe, so a reclamation that went
    # wrong could not cost us the results.  It calls the standalone
    # tool's own planner rather than reimplementing the walk, and
    # supplies only the workspace -- the default policy is the same
    # "spent once `done` with a result.toml" rule this producer
    # would apply anyway, so the two cannot diverge.
    #
    # It remains useful alongside --tidy-run, because the in-flight
    # prune deliberately leaves some units behind: one still nested
    # inside another (REFUSAL 5), and one whose scratch appeared
    # after it landed.  This is the second pass that takes them.
    if clean_after:
        reclaim_run_scratch(workspace)

    # Say again, at the end and all together, anything that went
    # wrong while tidying in flight.  Each was already printed when
    # it happened, but a campaign is hours long and that line is
    # far up the log by now.  This is the summary the user is
    # actually looking at when the run finishes, so a broken
    # filesystem assumption gets named where it will be seen -- and
    # named as something to look into, not as a run that failed,
    # because every database and every log above is intact.
    if prune_problems is not empty:
        report("")
        report("tidy-run: " + count(prune_problems) + " scratch "
               "tree(s) could NOT be removed.  The results above "
               "are complete and unaffected, but reclamation is "
               "not working -- check permissions and the mount, "
               "and sweep with tidy_scratch.py when fixed:")
        for (label, message) in prune_problems:
            report("    " + label + ": " + message)


function make_prune_callback(flight, workspace, scratch_root,
                             problems):
    # Layer (b) (DESIGN 6.2.12): build the on_outcome callback a
    # flight is handed when --tidy-run is asked for.  Kaleidoscope
    # fires it once per unit, in LANDING order, as that unit
    # reaches a terminal state.
    #
    # It closes over the FLIGHT, not over a fixed unit list,
    # because the climb's list ACCRETES as rungs are decided
    # (4e.7): reading flight.units at call time therefore sees
    # every rung decided so far, which is exactly the set that
    # could own scratch right now.
    #
    # `problems` is the shared list every prune failure is appended
    # to, so the end of the run can say plainly that tidying went
    # wrong.  A REFUSAL is not a problem and never lands there: a
    # unit skipped because it is unfinished, nested, or too recent
    # is the mechanism working.  A FAILURE is a removal that was
    # attempted and did not happen, or the mechanism itself
    # throwing -- and that means an assumption of 6.2.12 no longer
    # holds, which the user needs to know about.
    function on_landed(entry):
        # The comparison set for REFUSAL 5 is this flight's OTHER
        # units.  A unit not yet run owns no scratch and so
        # contributes nothing, which is right: it has nothing to
        # lose, and it will create its own tree when it runs.
        others = []
        for unit in flight.units:
            run_dir = unit_run_dir(flight, unit)
            if run_dir == entry.wingbeat_dir:
                continue
            (target, _) = scratch_target(run_dir, scratch_root)
            if target != None:
                others.append(target)

        # Nothing raised here may escape.  This runs inside
        # kaleidoscope's collect (13.5), which does NOT guard the
        # hook, so an exception would propagate out through the
        # climb and abandon the whole campaign -- losing a run of
        # cluster time to a piece of housekeeping.  Tidying may
        # fail; a campaign may not fail because tidying did.
        #
        # Contained is not the same as ignored.  A failure is
        # printed the moment it happens AND recorded in `problems`,
        # because a mid-campaign line is easily lost in hours of
        # scrollback and this is exactly the kind of thing that
        # must not be lost: a scratch tree that will not remove
        # means a permission, mount, or lock assumption has broken,
        # and the next campaign will hit it too.
        label = relpath(entry.wingbeat_dir,
                        join(workspace, "wingbeats"))
        try:
            record = reclaim_one_dir(
                entry.wingbeat_dir, scratch_root,
                other_targets = others, label = label)
        except Exception as exc:
            # The mechanism itself threw -- a stronger signal than
            # a removal that merely failed, since 6.2.12's refusals
            # are meant to leave no uncaught path.
            problems.append((label, "prune raised: " + str(exc)))
            report("tidy-run: FAILED to prune " + label + ": " +
                   str(exc))
            return

        if record.ok and not record.removed:
            problems.append((label, record.failure))
            report("tidy-run: FAILED to prune " + label + ": " +
                   record.failure)
        else:
            # A removal, or an ordinary refusal.  Both are the
            # mechanism working as designed, so both are
            # NARRATION under the 5.7 reporting rule: silent
            # unless --verbose asked to see them.  The two
            # FAILURE branches above print either way, because
            # those say an assumption has stopped holding.
            if record.removed:
                narrate("tidy-run: pruned " + label
                        + ", freeing " + human_bytes(record.bytes))
            else:
                narrate("tidy-run: kept " + label + ": "
                        + record.reason)
    return on_landed


function build_climb_config(ref, structure, confidence,
                            under_trained, thresholds, max_count):
    # Assemble one solid's ClimbConfig (DESIGN 3.12 / 5.7) from
    # three sources: the loaded cell's reciprocal geometry (the rung
    # mechanics search there), the confidence-derived policy, and
    # the energy / ceiling knobs.  Called once per solid in Phase 1b,
    # BEFORE any run -- the climb needs the axis classes to seed and
    # step, so they are recomputed here, never read back from a run.
    cell = load_structure(structure)          # StructureControl

    # Reciprocal-cell geometry (the density <-> mesh map, 4c.2): the
    # three reciprocal-axis magnitudes |b_i| and the reciprocal cell
    # volume, from the reciprocal lattice the mesh writer uses
    # (structure_control.make_inv_or_recip_lattice, which fills the
    # cell's recip_lattice and recip_cell_volume).  UNLIKE the real
    # lattice (whose ROWS are the cell vectors), the reciprocal is
    # stored as an inverse, so its reciprocal vectors b_i are its
    # COLUMNS: the magnitude of reciprocal axis i is the norm DOWN
    # column i.  (A row-wise read transposes the vectors and silently
    # mis-scales every non-cubic cell -- a cubic cell hides it, since
    # its row and column norms match.)
    make_inv_or_recip_lattice(cell)           # fills recip_lattice
    recip = cell.recip_lattice
    recip_mag = [norm(column i of recip) for i in 0..2]
    recip_cell_volume = cell.recip_cell_volume

    # Axis classes (4c.7 / DESIGN 2.7): which reciprocal axes must
    # share a k-point count.  Recomputed from the SAME two
    # ingredients imago uses at runtime -- the conventional-abc
    # space-group rotations and the loaded / conventional cells --
    # so a run's emitted RESOLVED_KP_CLASSES (4d.5) can later
    # self-check this value (the 6c self-test).  The rotations come
    # from the shared reader (4b.4) via the cell's `point_ops`
    # accessor -- the SAME parser the kp writer uses, so the climb
    # seeds from the operations imago will run under.  cell_mode is
    # "prim" for a primitive reduction, else "full".
    classes = axis_classes_for_cell(
        cell.point_ops,
        cell.loaded_lattice, cell.conventional_lattice,
        cell.cell_mode)

    # Confidence -> the shape of the search (4e.4): dispatch mode,
    # flatness persistence, grid width, climb start offset.  The
    # knob thresholds were resolved once for the run from the
    # manifest (climb_policy_from_manifest, above).
    policy = resolve_climb_policy(confidence, under_trained,
                                  thresholds)

    # The crystalline floor rung (DESIGN 3.12.4 / 4e.4): the mesh
    # whose densest axis holds crystalline_floor_axis_count points
    # and whose other axes scale DOWN from there by their |b_i|
    # ratios (never exceeding the cap), so [4,4,4] on a cubic cell
    # and [4,4,2] / [4,3,2] / [4,1,1] on an anisotropic one.  None
    # for a non-crystalline solid, which seeds at or near Gamma (7.9)
    # and must not be floored up off it.  Built from the same recip
    # geometry and classes the seeding will read, so the floor is a
    # rung the climb can actually open on.
    if ref.system_type == "crystalline":
        opening_floor = crystallineFloorMesh(
            recip_mag, classes,
            thresholds.crystalline_floor_axis_count)
    else:
        opening_floor = None

    return ClimbConfig(
        classes           = classes,
        recip_mag         = recip_mag,
        recip_cell_volume = recip_cell_volume,
        mode              = policy.mode,
        flat_needed       = policy.flat_needed,
        grid_width        = policy.grid_width,
        start_offset      = policy.start_offset,
        opening_floor     = opening_floor,
        cell_atom_count   = cell.num_atoms,
        threshold         = ref.kpoint_convergence_threshold,
        max_count         = max_count)


function discover_environments(converged, ref):
    # DESIGN 5.7 / 5.2.3: yield one record per DISTINCT
    # ENVIRONMENT in the converged run.  The run's site-identity
    # map (datSkl.map) partitions every atom by
    # (element, species, type); atoms sharing that key are
    # equivalent under the assigning method (symmetry for a
    # crystalline reference) and carry the same potential, so
    # one representative speaks for the whole group.  The
    # representative is order-independent -- the lowest skeleton
    # index in the group (DESIGN 5.6.5).
    #
    # A customization annotates the environment that contains
    # its pinned atom_site, supplying the curator's label /
    # default / description / fingerprint overrides and the
    # representative site to harvest; an auto-discovered
    # environment derives its label, takes default=false, and
    # auto-composes its description from ref.source_description.
    # A site-less customization cannot yet be matched and is
    # skipped (INTERIM C88; matching by label is later work).
    site_identity = read_site_identity_map(converged)
    sites_by_env = {}
    for skeleton_atom, key in site_identity.items():
        sites_by_env.setdefault(key, []).append(skeleton_atom)

    # Map each pinned customization to its environment; two on
    # one environment is ambiguous and is an error.
    custom_by_env = {}
    for spec in ref.entries:
        if spec.atom_site is None:
            continue
        require(spec.atom_site in site_identity, ...)
        key = site_identity[spec.atom_site]
        require(key not in custom_by_env, ...)
        custom_by_env[key] = spec

    envs = []
    for key in sorted(sites_by_env):
        (site_element, species, type_number) = key
        spec = custom_by_env.get(key)
        atom_site = (spec.atom_site if spec and spec.atom_site
                     else min(sites_by_env[key]))
        element = resolve_element(spec, site_element)   # cross-check
        label = (spec.label if spec and spec.label
                 else assemble_entry_label(ref.reference_id,
                     element, species, type_number, atom_site))
        description = (spec.description if spec and spec.description
                       else compose_auto_description(
                           ref.source_description, ref.reference_id,
                           element, species, atom_site))
        envs.append(Environment(
            element     = element,
            atom_site   = atom_site,
            species     = species,
            type_number = type_number,
            label       = label,
            default     = spec.default if spec else false,
            description = description,
            overrides   = spec.fingerprints if spec else []))
    return envs


function insert_or_skip(db, new, ref):
    # DESIGN 5.2.3: the per-element database stores DISTINCT
    # ENVIRONMENTS, not atoms.  An explicit customization
    # label replaces by label first; otherwise the dedup
    # LOCATES a match on the BISPECTRUM descriptor at the
    # preferred sub_spec -- the transferable one every entry
    # carries (native or witness, 5.2.2); symmetry has no
    # meaning across structures and reduce is the disposable
    # witness, so neither gates the dedup.  A match is
    # SKIPPED (the stored representative stands); a new
    # environment is appended.  Skipping is what makes the
    # build idempotent: re-running an unchanged manifest finds
    # every environment already present and moves nothing.
    if new.label is not None:
        existing = find_by_label(db, new.label)
        if existing is not None:
            replace(db, existing, new)
            return
    # INTERIM (C88): the fingerprint dedup is not built yet.
    # The present build appends (the C60 witness path);
    # duplicates are not yet detected.
    dup = find_bispectrum_duplicate(db, new)        # C88
    if dup is None:
        db.potentials.append(new)
    # else: a duplicate environment is already stored, so skip
    # `new` entirely -- the first representative's potential
    # is kept (DESIGN 5.2.3).  Reconciling duplicates into a
    # statistical mean (with the spread, counts, and the
    # alpha-set-equality assertion that averaging needs) is
    # the deferred upgrade (TODO C103).


function load_manifest_v2(path):
    # Strict-refusal validator implementing manifest
    # rules 1-11 of DESIGN 5.7.  Every failure names
    # the failing rule and the offending entry; no
    # warning-and-continue path exists.
    raw = tomllib.load(path)

    # Rule 1: schema_version must equal 2.
    require(raw.get("schema_version") == 2, path,
        "manifest rule 1: schema_version must equal"
        + " 2 (found "
        + str(raw.get("schema_version")) + ")")

    # Rule 10: the [characterization] block declares the
    # database-wide preferred recipe -- at most one
    # fingerprint per method, each a known matcher (rule 9).
    # That single declaration per method IS the family's
    # preferred record; per-entry declarations may not be
    # preferred (checked in the entry loop below).
    char_methods = set()
    for fp in raw.get("characterization", {}).get(
            "fingerprint", []):
        require("method" in fp and "sub_spec" in fp, path,
            "manifest rule 10: [characterization]"
            + " fingerprint needs method and sub_spec")
        require(fp["method"] in MATCHERS, path,
            "manifest rule 9: unknown matcher method '"
            + fp["method"] + "' in [characterization]")
        require(fp["method"] not in char_methods, path,
            "manifest rule 10: method '" + fp["method"]
            + "' declared twice in [characterization]"
            + " (the preferred recipe has one home)")
        char_methods.add(fp["method"])

    # Rule 2: the recipe is required.  A manifest must declare a
    # [characterization] block with at least one fingerprint, so
    # the build cannot silently produce a database with no
    # preferred descriptors for the consumer to match against
    # (VISION Principles 5 and 11).  The relaxed
    # load_structure_sources reader does NOT apply this -- it
    # only materializes structures and never harvests.
    require(len(char_methods) > 0, path,
        "manifest rule 2: a [characterization] block"
        + " declaring at least one fingerprint is required")

    # The optional [defaults] block: shared run settings a solid
    #   inherits when it omits them (DESIGN 5.7).  Only the five
    #   run-setting keys are meaningful.  Solids stay sparse here
    #   (a missing setting is None) and are resolved after load by
    #   apply_manifest_defaults (below) -- the writer (11.6) needs
    #   them sparse to emit the compact [defaults]-plus-overrides
    #   form.
    raw_defaults = { k : raw["defaults"][k]
                     for k in RUN_SETTING_KEYS
                     if k in raw.get("defaults", {}) }
    # The optional [harvest] block: shared HARVEST settings, read
    #   back after the runs finish rather than fed into any run
    #   (DESIGN 5.7).  kpoint_convergence_threshold carries a built-in
    #   producer default, and the [harvest.kpoint_climb] sub-table's
    #   knobs each do too (3.12.6), so the block -- and every key --
    #   may be omitted entirely.  The sub-table's keys are validated
    #   against KPOINT_CLIMB_KEYS so a mistyped knob fails loudly at
    #   load rather than silently taking its default.
    raw_harvest = raw.get("harvest", {})
    climb = raw_harvest.get("kpoint_climb", {})
    for knob in climb:
        require(knob in KPOINT_CLIMB_KEYS, path,
            "unknown [harvest.kpoint_climb] knob: " + knob)

    solids               = raw.get("reference_solid",
                                   [])
    seen_ref_ids         = set()
    seen_element_label   = set()
    default_per_element  = {}    # elem -> count of
                                 #   default=true
                                 #   customizations

    for ref in solids:
        # Rule 2: reference_id and system_type are always named
        # per solid.  The five run settings may be named here OR
        # inherited from [defaults], but each must be RESOLVABLE
        # one way or the other -- together with system_type they
        # select the guidance predictor's sub-model (DESIGN 7.6)
        # and land on every produced entry, so nothing emitted
        # rides on an implicit default (VISION Principles 5, 11).
        for f in ("reference_id", "system_type"):
            require(f in ref, path,
                "manifest rule 2: [[reference_solid"
                + "]] missing field: " + f)
        for f in RUN_SETTING_KEYS:
            if f in EXEMPT_RUN_SETTING_KEYS:
                continue
            require(f in ref or f in raw_defaults, path,
                "manifest rule 2: [[reference_solid "
                + ref.get("reference_id", "?")
                + "]] run setting " + f + " not resolvable"
                + " (absent here and from [defaults])")
        # `cell` is EXEMPT: it carries a built-in default
        # (DEFAULT_CELL) and is recorded nowhere -- it selects no
        # predictor sub-model and no harvested value depends on
        # it -- so a manifest that never names a cell leaves no
        # provenance gap, only accepts the conventional cell.
        # The exemption ends the moment `cell` is recorded on an
        # entry: it becomes emitted knowledge and rejoins the
        # rule (DESIGN 5.7).
        # The harvest setting kpoint_convergence_threshold is
        # EXEMPT from this resolvability rule: it carries a
        # built-in default (2e-3 eV/atom; DESIGN 5.7 / 7.8 and the
        # noise-floor rule of 3.12.3), so a solid naming neither it
        # nor a [harvest] block is accepted --
        # apply_manifest_defaults supplies the default.

        rid = ref["reference_id"]

        # Rule 2 (domain): a named cell must be one of the two
        # valid values, wherever it is named.  A typo would
        # otherwise reach structure_control, which refuses any
        # token that is neither `full` nor `prim` -- but only
        # after the fetch, so catching it at load is both earlier
        # and clearer.
        for source in (ref, raw_defaults):
            if "cell" in source:
                require(source["cell"] in VALID_CELLS, path,
                    "manifest rule 2: cell must be one of "
                    + " / ".join(VALID_CELLS) + " (found "
                    + str(source["cell"]) + ")")

        # Rule 2 (domain): system_type must be one of
        # the four valid values; the guidance predictor
        # (DESIGN 7) switches its sub-model on it.
        require(ref["system_type"] in (
                "crystalline", "amorphous",
                "nanostructure", "molecular"),
            path,
            "manifest rule 2: [[reference_solid "
            + rid + "]] system_type must be one of"
            + " crystalline / amorphous /"
            + " nanostructure / molecular (found "
            + str(ref["system_type"]) + ")")

        # Rule 4: exactly one of cod_id /
        # structure_path; cod_revision required and
        # non-empty whenever cod_id is set.
        has_cod = "cod_id" in ref
        has_pth = "structure_path" in ref
        require(has_cod != has_pth, path,
            "manifest rule 4: [[reference_solid "
            + rid + "]] must set exactly one of"
            + " cod_id or structure_path")
        if has_cod:
            require("cod_revision" in ref
                    and len(ref["cod_revision"]) > 0,
                path,
                "manifest rule 4: cod_revision"
                + " required (non-empty) when"
                + " cod_id is set (" + rid + ")")
        else:
            sp = path_join(dirname(path),
                ref["structure_path"])
            require(file_exists(sp), path,
                "manifest rule 4: structure_path"
                + " resolves to a missing file: "
                + sp)

        # Rule 5: reference_id uniqueness.
        require(rid not in seen_ref_ids, path,
            "manifest rule 5: duplicate"
            + " reference_id: " + rid)
        seen_ref_ids.add(rid)

        # Per-entry checks (rules 3, 6, 7, 8, 9, 10).
        # Entries are OPTIONAL customizations (rule 3); a
        # solid may carry none, and every customization
        # field is optional.  atom_site, when given, is the
        # representative-atom handle; element is cross-
        # checked at harvest; label / default / description
        # override the derived / auto-composed values.
        for entry in ref.get("entry", []):
            elem  = entry.get("element")
            label = entry.get("label")

            # Rule 6: only an EXPLICIT customization label needs
            # the cross-manifest (element, label) check;
            # derived labels are unique by construction per
            # environment (DESIGN 5.7).
            if label is not None and elem is not None:
                key = (elem, label)
                require(key not in seen_element_label,
                    path,
                    "manifest rule 6: duplicate"
                    + " (element, label): " + str(key))
                seen_element_label.add(key)

            # Rule 7 (load half): no element may carry two
            # default=true customizations.  The "exactly one
            # harvested element" half is checked at harvest,
            # once the run reveals which elements exist.
            if entry.get("default", False) and (
                    elem is not None):
                default_per_element[elem] = (
                    default_per_element.get(elem, 0) + 1)

            # Per-entry fingerprints are RARE overrides
            # (extra NON-preferred sub_specs); the preferred
            # recipe lives in [characterization] (rule 10).
            seen_fp = set()
            for fp in entry.get("fingerprint", []):
                require("method" in fp
                        and "sub_spec" in fp, path,
                    "manifest rule 8: fingerprint"
                    + " declaration must carry both"
                    + " method and sub_spec (" + rid + ")")
                # Rule 10: a per-entry declaration may not
                # be preferred -- the preferred record is
                # the [characterization] one.
                require(not fp.get("preferred", False),
                    path,
                    "manifest rule 10: per-entry"
                    + " fingerprint may not be preferred ("
                    + rid + "); set it in"
                    + " [characterization]")
                # Rule 9: method must be a known matcher.
                require(fp["method"] in MATCHERS, path,
                    "manifest rule 9: unknown matcher"
                    + " method '" + fp["method"] + "' ("
                    + rid + ")")
                # Rule 8: (method, sub_spec) unique within
                # this entry.  Same canonicalization the
                # per-element-database reader uses (11.1).
                canon = canonicalize_sub_spec(
                    fp["sub_spec"])
                k2 = (fp["method"], canon)
                require(k2 not in seen_fp, path,
                    "manifest rule 8: duplicate"
                    + " (method, sub_spec) in " + rid)
                seen_fp.add(k2)

    # Rule 7 (load half): no element may carry more than
    # one default=true customization.  Zero is allowed here --
    # element may take its default at harvest, or the
    # isolated baseline is the default when the manifest
    # adds nothing for it.  The "exactly one per HARVESTED
    # element" check runs at harvest time (DESIGN 5.7).
    for elem, count in default_per_element.items():
        require(count <= 1, path,
            "manifest rule 7: element " + elem
            + " has " + str(count)
            + " default=true customizations (at most one)")

    # The manifest object carries the two shared blocks alongside
    #   the sparse solids -- manifest.defaults (the five run
    #   settings) and manifest.harvest (the harvest block) -- so
    #   apply_manifest_defaults can resolve each solid after load.
    return parse_manifest_object(raw, raw_defaults, raw_harvest)


# The six run settings that may live in [defaults] and be
#   inherited per solid (DESIGN 5.7).  system_type is NOT among
#   them -- it is structure metadata, always named per solid.
RUN_SETTING_KEYS = ("basis", "functional",
    "kpoint_integration", "kpoint_spec", "scf_threshold",
    "cell")

# The cell a reference run computes in (DESIGN 5.7): the
#   conventional cell of the structure's space group, or its
#   primitive reduction.  A cost setting, not a physics one --
#   the harvested potential and every fingerprint are
#   cell-invariant -- so it is NOT a predictor sub-model
#   selector.  The default is the PRIMITIVE reduction, which
#   costs about half as much per converged calculation and
#   agrees with the conventional cell to 0.002 meV/atom.
VALID_CELLS  = ("full", "prim")
DEFAULT_CELL = "prim"

# The run settings exempt from rule 2's resolvability requirement,
#   because they are recorded nowhere and so leave no provenance
#   gap when omitted (DESIGN 5.7).  A key leaves this set the
#   moment the producer starts emitting it.
EXEMPT_RUN_SETTING_KEYS = ("cell",)

# The producer's built-in k-point flatness tolerance, used when a
#   solid names neither its own kpoint_convergence_threshold nor a
#   [harvest] block (DESIGN 5.7 / 7.8).  Per atom, in eV.
#
# The value is a FLOOR set by the ladder, not a taste (DESIGN
#   3.12.3).  It must sit ABOVE the rung-to-rung scatter of the
#   energies it judges, because a bar beneath the noise can be
#   cleared only by two coincidences in a row.  The former 5.0e-4
#   sat below all measured scatter and converged two of thirteen
#   seed solids.
#
# The scatter that once argued for 2e-3 -- 0.0008 to 0.0047
#   eV/atom -- was measured on METALS, whose energies oscillate as
#   the mesh crosses the Fermi surface.  Those ladders no longer
#   reach this test at all: the gap test stops a metal on every
#   search shape before any convergence work (4e.3).  What is left
#   to judge is insulators, which settle rather than oscillate,
#   and for the six ordinary si_fd-3m seeds 1e-3 and 2e-3 pick the
#   IDENTICAL mesh, [10,10,10].
#
# Tightening buys gap quality: gap_ev is read off whichever rung
#   the climb stopped on and is a predictor key nothing downstream
#   re-converges, so a looser bar records a coarser-mesh gap
#   (DESIGN 7.6, and the defect carried as D22).
DEFAULT_KPOINT_CONVERGENCE_THRESHOLD = 1.0e-3    # 1 meV/atom

# The adaptive-climb tuning knobs that may live in the optional
#   [harvest.kpoint_climb] sub-table (DESIGN 5.7 / 3.12.6).  All but
#   max_count name the confidence-to-policy PolicyThresholds (4e.4):
#   the confidence split, the grid width and two start offsets, the
#   two flat-rung persistence counts, the stride cap and climb shape,
#   the stride-flatness multiple and the metal_gap_threshold, and
#   the crystalline_floor_axis_count that sets the opening floor
#   (3.12.4).  max_count is the per-axis ceiling (4e.2).  Each carries
#   a provisional built-in default (mesh_climb), so the sub-table --
#   and any knob -- may be omitted.  Database-wide: no per-solid
#   override.
KPOINT_CLIMB_KEYS = ("confidence_high", "grid_width",
    "start_offset_moderate", "start_offset_cold",
    "flat_needed_confident", "flat_needed_cold",
    "max_stride", "climb_shape", "stride_flatness_multiple",
    "metal_gap_threshold", "crystalline_floor_axis_count",
    "max_count")


function apply_manifest_defaults(manifest):
    # Producer-side resolve pass (build_initial_potentials), run
    #   once after load so the pipeline reads fully-populated
    #   solids while the library keeps them sparse for the compact
    #   writer (11.6).  One shared step fills BOTH shared blocks.
    manifest.reference_solids = [
        resolve_settings(solid, manifest.defaults,
                         manifest.harvest)
        for solid in manifest.reference_solids]


function resolve_settings(solid, defaults, harvest):
    # DESIGN 5.7: return a copy of `solid` with its settings
    #   filled from the shared blocks.  Precedence is the solid's
    #   own value first, then the shared block; the harvest setting
    #   alone falls back further, to a built-in producer default.
    #
    # Run settings have NO built-in fallback: the loader already
    #   proved each one resolvable from the solid or [defaults]
    #   (rule 2), so pick() never returns None for them.
    function pick(key):
        own = getattr(solid, key)
        return own if own is not None else defaults.get(key)

    # The one harvest setting: solid's own -> [harvest] -> the
    #   built-in default.  This is the arm C112 adds; the run-
    #   setting arm above already ships as resolve_run_settings.
    threshold = solid.kpoint_convergence_threshold
    if threshold is None:
        threshold = harvest.get(
            "kpoint_convergence_threshold",
            DEFAULT_KPOINT_CONVERGENCE_THRESHOLD)

    return copy_of(solid,
        basis              = pick("basis"),
        functional         = pick("functional"),
        kpoint_integration = pick("kpoint_integration"),
        kpoint_spec        = pick("kpoint_spec"),
        scf_threshold      = pick("scf_threshold"),
        kpoint_convergence_threshold = threshold)


function structure_cache_dir(pdb_root):
    # Where materialized reference structures are cached
    # (ARCHITECTURE 8.1): under the producer's `curation`
    # tree beside its workspace and run log, NOT under a
    # database.  Everything there is reconstructible by
    # re-running the producer, so the whole campaign
    # footprint clears in one gesture.  Derived from
    # pdb_root's parent so a relocated data root carries
    # the cache with it.
    return dirname(pdb_root) + "/curation/structures/"


function materialize_structure(ref):
    # Option A (DESIGN 5.7): guarantee that the
    # reference solid's structure exists as a local
    # file and return its path.  This is the producer's
    # ONLY network access and is deliberately decoupled
    # from any run cache -- it carries no SCF state and
    # makes no hit/miss decision.  Recompute avoidance
    # belongs to kaleidoscope's run-reuse cache (DESIGN
    # 6.2.5), which keys on this file's contents.
    if ref.structure_path is not None:
        # Disk read; the loader already resolved the
        # path under the manifest directory (rule 4).
        # No network.  The curator's skeleton carries its
        # OWN full/prim token and that token stands --
        # ref.cell governs only what the producer writes,
        # and the producer never rewrites a curator's file
        # (DESIGN 5.7).
        return ref.structure_path

    # cod_id ref: fetch the pinned revision once to a
    # plain local location.  Strict on failure (network
    # down / COD outage / pinned revision missing) --
    # never falls back to another revision, because a
    # silent fallback would desync the build from the
    # pinned manifest (DESIGN 5.7).
    cache = "share/curation/structures/"
    cif = cache + ref.reference_id + cod_extension(ref)
    # The CACHED SKELETON's name carries every manifest
    # setting that changes what the conversion writes --
    # today just `cell` (DESIGN 5.7).  Without the
    # qualifier a later run under a different cell is
    # handed the earlier run's file: no error, a
    # well-formed skeleton, and the wrong answer reported
    # as a success.  The CIF needs no qualifier, since
    # cod_id + cod_revision already pin its bytes.
    skl = cache + ref.reference_id + "-" + ref.cell + ".skl"

    if not file_exists(skl):
        if not file_exists(cif):
            fetch_cod_structure(
                cod_id       = ref.cod_id,
                cod_revision = ref.cod_revision,
                dest         = cif)
        # Convert with the space group PRESERVED, writing
        # ref.cell as the skeleton's full/prim token.  A
        # CIF whose space group cannot be resolved to a
        # spaceDB setting is a hard error -- no silent P1
        # fallback for a crystal (DESIGN 5.7).
        cif_to_skeleton(cif, skl, cell = ref.cell,
                        title = ref.reference_id)
    return skl


function load_structure_sources(path):
    # The relaxed read behind --materialize-only (DESIGN
    # 5.7), and the read expand_manifest (11.6) uses to
    # complete a sketch.  It enforces only what materializing
    # a structure needs -- schema version (rule 1), a
    # label-safe unique reference_id (rule 5), and exactly one
    # structure source (rule 4) -- and captures system_type
    # when present.  Run and harvest fields are left at
    # placeholders; this view is never dispatched, only
    # materialized, so the full rule set does not apply.
    raw = toml_load(path)
    require(raw.schema_version == 2)               # rule 1
    seen = empty_set
    sources = []
    for ref in raw.reference_solid:
        require("reference_id" in ref)
        rid = ref.reference_id
        require(is_label_safe(rid))                # rule 5
        require(rid not in seen); seen.add(rid)
        require(has(ref,"cod_id") != has(ref,"structure_path"))
        # cod_id (positive int + non-empty cod_revision) or a
        # structure_path; the path's existence is checked at
        # materialize time, not here.
        sources.append(ReferenceSolid(
            reference_id   = rid,
            system_type    = ref.get("system_type", ""),
            cod_id         = ref.get("cod_id"),
            cod_revision   = ref.get("cod_revision"),
            structure_path = ref.get("structure_path"),
            <run + harvest fields left at placeholders>))
    return sources


function materialize_only(manifest_path, pdb_root,
        cache_dir = None):
    # The --materialize-only pre-flight (DESIGN 5.7): fetch
    # and convert every reference structure, then STOP -- no
    # SCF dispatch.  Lets a curator get a freshly pinned set
    # materializing cleanly before the run and harvest fields
    # are filled in.  A per-solid failure does NOT abort the
    # batch: it is recorded and the next solid is tried.
    sources = load_structure_sources(manifest_path)
    if cache_dir is None:
        cache_dir = structure_cache_dir(pdb_root)
    report = []
    for ref in sources:
        try:
            skl = materialize_structure(ref)   # fetch + convert
            report.append(ok_row(ref, skl))
        except materialize_error as e:
            report.append(fail_row(ref, e))
    return report   # the CLI reports it: failures always, the
                    #   rest only when asked (print_materialize_
                    #   report, below)


# ---- Reporting (DESIGN 5.7) ----------------------------------
# Outcomes and problems always; narration only on request.  The
#   setting is module-level, established once from the parsed
#   flags before any work starts, because how the process talks
#   to its user is not a property of how a flight converges --
#   threading it down into the prune hook would put a reporting
#   concern inside four functions about physics and dispatch.

module_state narrate_progress = false

function set_verbosity(verbose):
    # Called by main() as its FIRST action, before the manifest
    #   is read or a structure fetched, so that every later
    #   report -- including one from a failure during start-up --
    #   is already governed by what the user asked for.
    narrate_progress = verbose

function narrate(message):
    # Per-item progress.  Silent unless --verbose was given.
    if narrate_progress:
        print(message)


function print_materialize_report(report):
    # A clean pre-flight is SILENT.  Eight structures that all
    #   arrived tell the curator nothing to act on, and a line
    #   printed on every successful run is a line that stops
    #   being read -- which is how a real failure would come to
    #   hide among them.
    failures = [row for row in report if not row.ok]

    for row in report:
        if row.ok:
            narrate("  [ok  ] " + row.reference_id + ": "
                    + row.source)
            narrate("          -> " + row.skl_path)
        else:
            # Always, and regardless of --verbose: this is the
            #   fetch error the curator needs in order to act.
            print("  [FAIL] " + row.reference_id + ": "
                  + row.source)
            print("          " + row.message)

    # The tally earns its line only when it says something: how
    #   much of the set survived a partial failure, or the
    #   confirmation someone asked for with --verbose.
    if failures or narrate_progress:
        print("materialize: " + count(ok rows) + "/"
              + count(report) + " reference structures fetched "
              + "and converted")


function main_submit_mode(argv, args, data_root):
    # `argv` is the flag vector main() parsed, threaded through
    # rather than re-read from the process (13.7 explains why).
    # --submit: run the producer as its OWN batch job (DESIGN
    # 6.2.11; the sbatch generator is 13.7).  Materialize every
    # structure HERE, on the login node, because it is the one
    # step that needs the network and a compute node may have
    # none.  Only then submit the batch job that re-runs this
    # producer, which finds every structure already cached.
    #
    # A materialize failure STOPS the submission: a batch job
    # that cannot read its structures would burn a queue slot to
    # fail, and the curator needs the fetch error, not a job id.
    report = materialize_only(args.manifest, args.pdb_root)
    print_materialize_report(report)
    if not all(row.ok for row in report):
        return error("materialize failed; not submitting")
    job_id = submit_orchestrator_batch(argv, args, data_root)
    print("submitted orchestrator batch job " + job_id)
    return ok

# --submit and --materialize-only are mutually exclusive: each
#   names a different stopping point (submit the driver's job
#   versus stop after the fetch), so the CLI rejects both at once.


function fingerprintDeclarations(characterization,
                                 overrides):
    # THE declaration-set rule (DESIGN 5.10.6).  Defined once
    # and used by BOTH ends of the producer: the harvest, to
    # know which fingerprints to store, and the build, to know
    # which `-loen -scf no` units to dispatch.  Writing it once
    # is what keeps the set that is harvested identical to the
    # set that was run; the drift this prevents cost a whole
    # flight before the rule existed.
    #
    # An environment harvests the database-wide
    # [characterization] preferred recipe (one sub_spec per
    # method, marked preferred = true), plus any rare per-entry
    # override the customization added (extra NON-preferred
    # sub_specs).  DESIGN 5.7.  Each record already carries its
    # own `preferred` flag, so composing is plain concatenation,
    # recipe first.
    return characterization + overrides


function producerFingerprintDeclarations(ref,
                                         characterization):
    # Every declaration ANY environment of this solid could
    # present at harvest (DESIGN 5.10.6).  The build side needs
    # this superset because the environments themselves are not
    # known yet -- they are discovered from the converged run --
    # so it applies the one rule above to each case that could
    # arise and unions the results:
    #   * an environment with no customization  -> no overrides;
    #   * an environment a customization annotates -> that
    #     entry's overrides.
    # A declaration for an entry no environment turns out to
    # claim (a site-less customization, which environment
    # discovery cannot yet match) merely builds one extra
    # geometry-only run.  That is the right direction to err:
    # a spare descriptor costs one short run, a missing one
    # costs the flight.  Callers dedup by calc tag, so the
    # repeated recipe collapses to one unit per (method,
    # sub_spec).
    decls = fingerprintDeclarations(characterization, [])
    for entry in ref.entries:
        decls = decls + fingerprintDeclarations(
            characterization, entry.fingerprints)
    return decls


function harvestFingerprints(flight, ref, env,
        result_toml, characterization):
    # This environment's declarations, by the one rule above --
    # the same rule the build side used to decide what to
    # dispatch, so what is read here was certainly run.
    decls = fingerprintDeclarations(
        characterization, env.overrides)
    if decls is empty:
        return []

    # The two matcher families split here (DESIGN 5.10).  A
    # Python-side declaration (reduce) computes in process from
    # the run's expanded structure, below; a Fortran-side one
    # (bispectrum) reads the descriptor of the `-loen -scf no`
    # unit the pre-flight already dispatched, via
    # harvestLoenFingerprint.  The split is by the matcher's
    # own needs_loen_run flag, so a new family joins the right
    # branch without editing this function.
    python_decls = [d for d in decls
                    if not MATCHERS[d["method"]]().needs_loen_run]

    # The run's EXPANDED full-cell structure is read only when a
    # Python-side declaration needs it -- a bispectrum-only
    # entry never touches it -- and then only once.  Read it
    # (outputs["structure"], makeinput's imago.fract-mi) and
    # build its minimum-image distance matrix ONCE, sized to
    # the LARGEST cutoff any declaration requests, then reuse
    # the one structure for them all.  Reading the run's own
    # expansion -- not re-expanding the materialized source --
    # reuses the exact geometry and numbering the run computed
    # and avoids duplicating applySpaceGroup.  No subprocess
    # and no on-disk cache: recomputing the shells in process
    # is cheaper than the cache bookkeeping would be.  Sharing
    # one matrix built to the max cutoff is safe because
    # compute_query independently trims neighbors to each
    # declaration's own sub_spec cutoff, so a smaller request
    # ignores the matrix's extra reach (periodic boundary
    # conditions enter only here).  The [characterization]
    # recipe plus any override can contribute several reduce
    # sub_specs, so building once to the max cutoff genuinely
    # applies -- but only over the PYTHON-side declarations,
    # since a loen-side sub_spec's cutoff sizes the engine run,
    # not this matrix.

    # The map step comes first, because BOTH families need it.
    # The expanded skeleton is ordered by the run's sorted
    # (dat) numbering, but atom_site is a skeleton index, so
    # map it to the structure row through datSkl.map; the map
    # yields both the row and that row's element symbol.
    # Resolve it once and reuse the pair for every declaration.
    (dat_index, map_element) = skeleton_to_dat(
        result_toml.outputs["datSkl_map"])[env.atom_site]

    structure = None
    if python_decls is not empty:
        max_cutoff = max(d["sub_spec"]["cutoff"]
                         for d in python_decls)
        structure = read_structure(
            result_toml.outputs["structure"])
        build_min_dist_matrix(structure, max_cutoff)

        # Guard the numbering assumption: the structure row and
        # the map must name the same element, or the expansion
        # and the map have desynced and the fingerprint would
        # describe the wrong atom.  Strict refusal beats a
        # silent mismatch.  (The loen branch guards the same
        # way, against the descriptor's own identity column.)
        if lower(structure.atom_element_name[dat_index])
                != lower(map_element):
            raise ValueError(
                f"site {env.atom_site}: datSkl.map names "
                f"{map_element} but the expanded structure row "
                f"{dat_index} is a different element; "
                f"numbering desync")

    fingerprints = []
    for d in decls:
        method   = d["method"]
        sub_spec = d["sub_spec"]
        matcher  = MATCHERS[method]()
        if matcher.needs_loen_run:
            # Fortran-side: read the descriptor of the unit the
            # pre-flight already dispatched.  No engine run
            # happens here (DESIGN 5.10).
            payload = harvestLoenFingerprint(
                flight, ref, dat_index, map_element,
                matcher, sub_spec)
        else:
            # Python-side: in-process compute against the shared
            # structure; compute_query trims to this
            # declaration's own sub_spec cutoff.
            vectors = matcher.compute_query(structure, sub_spec)
            payload = matcher.build_payload(vectors[dat_index])
        # Both branches wrap their vector via
        # matcher.build_payload, so the per-matcher payload
        # field name (DESIGN 5.2: bispec uses `values`, reduce
        # uses `shell_code`) flows through one accessor and the
        # two families stay symmetric on field naming.
        fingerprints.append(FingerprintRecord(
            method    = method,
            sub_spec  = sub_spec,
            preferred = d["preferred"],
            payload   = payload))
    return fingerprints


function buildLoenUnits(ref, struct_path, options, record,
                        characterization):
    # One structure-only `imago -loen -scf no` unit per distinct
    # Fortran-side declaration (DESIGN 5.10).  A bispectrum
    # descriptor is computed by the engine, so each declaration
    # needs its own dispatched run -- but the bispectrum is
    # geometry-only, so these runs need no converged SCF and
    # belong to no climb round.
    #
    # The declaration set comes from the SAME rule the harvest
    # applies, widened to every environment this solid could
    # present (DESIGN 5.10.6).  That is the whole point: the
    # build set is a superset of the harvest set by
    # construction, not by inspection.
    declarations = producerFingerprintDeclarations(
        ref, characterization)

    units     = []
    seen_tags = {}
    for d in declarations:
        matcher = MATCHERS[d["method"]]()
        if not matcher.needs_loen_run:
            continue      # Python-side: harvested in process.
        # One run serves every site sharing a (method,
        # sub_spec), because the descriptor table holds one row
        # per atom.  The calc tag is the dedup key, and it is
        # the same tag harvestLoenFingerprint rebuilds to find
        # the run directory.
        calc_tag = "loen-" + matcher.name + "-" + \
                   sub_spec_slug(d["sub_spec"])
        if calc_tag in seen_tags:
            continue
        seen_tags.add(calc_tag)

        # Layer the loen overrides on a COPY of the solid's
        # options, so the convergence units are untouched.
        loen_options = copy(options)
        loen_options["job"]       = "loen"
        loen_options["scf_basis"] = "no"
        loen_options["loeninput"] = loen_input_values(
            matcher, d["sub_spec"])
        # `record` carries the build identity here too (DESIGN
        #   6.2.4): a reused loen run is reported in the reuse plan
        #   with the build behind it, exactly as a reused rung is.
        units.append(CalcUnit(
            id        = ref.reference_id,
            structure = struct_path,
            calc      = [calc_tag],
            options   = loen_options,
            record    = record,
            kind      = "fingerprint"))
    return units


function assertLoenCoverage(units, refs, characterization):
    # The pre-dispatch invariant (DESIGN 5.10.6).  Structural
    # agreement is an argument, and an argument can be wrong, so
    # assert it directly: every Fortran-side declaration the
    # harvest could read must already have a unit in the flight.
    #
    # Called once the units are assembled and BEFORE any is
    # sent.  It is a set comparison over calc tags, so it costs
    # nothing -- and it converts build/harvest drift from a
    # failure found after minutes of cluster SCF time into one
    # raised before a single job is submitted.
    dispatched = {(u.id, u.calc[0]) for u in units
                  if u.kind == "fingerprint"}
    for ref in refs:
        for d in producerFingerprintDeclarations(
                ref, characterization):
            matcher = MATCHERS[d["method"]]()
            if not matcher.needs_loen_run:
                continue
            calc_tag = "loen-" + matcher.name + "-" + \
                       sub_spec_slug(d["sub_spec"])
            if (ref.reference_id, calc_tag) not in dispatched:
                # Name the solid AND the sub_spec: the curator
                # needs to know which declaration went unrun,
                # not merely that one did.
                raise ValueError(
                    f"{ref.reference_id}: fingerprint "
                    f"{d['method']} {d['sub_spec']} has no "
                    f"dispatched loen unit; the harvest would "
                    f"read a descriptor that was never "
                    f"computed")


function harvestLoenFingerprint(flight, ref,
        dat_atom, element, matcher, sub_spec):
    # Read the descriptor of the `-loen -scf no` unit that
    # kaleidoscope already dispatched for this
    # (solid, method, sub_spec) back in step 1b.  No loen run
    # happens here and there is no separate loen cache --
    # kaleidoscope's run-reuse cache (DESIGN 6.2.5) already
    # owns recompute avoidance.  The unit's run directory
    # follows the calc-tag convention (DESIGN 6.2.4): id =
    # reference_id, calc = "loen-<method>-<slug>"; the slug
    # encodes the sub_spec so two declarations differing in
    # any key or value land in different run directories by
    # construction.  Rebuilding the tag here from the SAME
    # helper the build side used is what makes the two agree.
    slug     = sub_spec_slug(sub_spec)
    calc_tag = "loen-" + matcher.name + "-" + slug
    run_dir  = unit_run_dir(flight.root,
        ref.reference_id, calc_tag)
    # The descriptor file is located by makegroups' finder
    # rather than named literally: the engine writes
    # `<edge>_loen<basis>.plot`, whose name varies with the
    # run's edge and basis (DESIGN 5.10.3).
    out_path = find_loen_descriptor(run_dir)

    rows = matcher.parse_loen_output(out_path, sub_spec)

    # The descriptor table holds one row per atom in dat order,
    # and `dat_atom` is the row this site maps to -- the CALLER
    # resolved it through datSkl.map, so both families index by
    # the same resolved row rather than each re-deriving it.
    # Index by the row's own site column instead of by position:
    # the table is self-describing (5.10.3), so a table that
    # omits or reorders a site is caught here rather than
    # silently returning a neighbour's fingerprint.
    rows_by_site = {row.site: row for row in rows}
    if dat_atom not in rows_by_site:
        raise ValueError(
            f"loen descriptor {out_path} has no row for dat "
            f"site {dat_atom} (it holds {sorted(rows_by_site)}); "
            f"the descriptor and the run's numbering desynced")
    row = rows_by_site[dat_atom]
    # The same element cross-check the Python-side branch makes,
    # against the descriptor's own identity column: a numbering
    # desync must fail loudly, never store the wrong atom's
    # fingerprint.
    if lower(row.element) != lower(element):
        raise ValueError(
            f"dat site {dat_atom}: datSkl.map names {element} "
            f"but the loen descriptor row is {row.element}; "
            f"numbering desync")
    # The matcher's build_payload accessor wraps the vector in
    # the per-matcher payload shape (DESIGN 5.2: bispec uses
    # `values`, reduce uses `shell_code`) so producer and
    # consumer stay symmetric on field naming.
    return matcher.build_payload(row.vector)


function sub_spec_slug(sub_spec):
    # Deterministic slug for the -loen unit's calc tag
    # (DESIGN 6.2.4).  Keys in alphabetical order,
    # joined as "key_value" segments, hyphen-separated.
    # Floats format as "%.6g" -- long enough to
    # disambiguate the parameters humans actually pick,
    # short enough to fit on a calc-tag line.
    parts = []
    for k in sorted(sub_spec.keys()):
        v = sub_spec[k]
        if isinstance(v, float):
            parts.append(k + "_"
                + sprintf("%.6g", v))
        else:
            parts.append(k + "_" + str(v))
    return "-".join(parts)


function build_isolated_entry(elem, commit, ts,
        manifest):
    # The isolated entry is rebuilt from current
    # atomSCF output every run.  Its `default` flag
    # is true iff no customization designates a default
    # environment for `elem`, so the isolated
    # baseline is the FALLBACK default and the per-
    # element database always has exactly one
    # default-tagged entry (rule 7 of 5.2).  As a
    # It stores the single atomSCF potential verbatim,
    # like any other environment (DESIGN 5.2.3).
    pot1   = read_pot1(elem)
    coeff1 = read_coeff1(elem)
    return PotentialEntry(
        label         = "isolated",
        default       = is_isolated_default_for(
                            elem, manifest),
        description   = ("Single isolated " + elem
                       + " atom (from atomSCF)."),
        num_gaussians = pot1.num_gaussians,
        alpha_min     = pot1.alpha_min,
        alpha_max     = pot1.alpha_max,
        coefficients    = coeff1.coefficients,
        alphas          = coeff1.alphas,
        provenance    = {
            "source"       : "atomSCF",
            "commit"       : commit,
            "generated_at" : ts},
        fingerprints  = [])


function is_isolated_default_for(elem, manifest):
    # True iff no customization designates a default
    # environment for `elem`.  The default tag goes
    # to the customized environment when one
    # exists; otherwise the isolated baseline is the
    # fallback default, so the per-element file always
    # carries exactly one default (rule 7 of 5.2) with
    # no missing-default error.  This also covers the
    # common case of an element with no manifest
    # contribution, whose only entry is the isolated
    # baseline.  `default` is an optional customization
    # field, false when absent.
    for ref in manifest.reference_solids:
        for entry in ref.entries:
            if (entry.element == elem
                    and entry.default):
                return False
    return True


function element_path(elem):
    return ("share/atomicPDB/" + lower(elem)
            + "/s_gaussian_pot.toml")
```

---

### 11.5 Validation Harness (DESIGN 5.8)

Compares iteration counts under "isolated" vs
"default_solid" across a benchmark set, computes the
mean reduction, and gates PASS/FAIL on the >=20%
threshold from VISION Principle 7.

```
function benchInitialPotential(benchmark_path):
    manifest = load_benchmark(benchmark_path)

    # Held-out sanity check
    curated_ids = curation_reference_ids()
    held_out = [t for t in manifest.tests
                if t.reference_id not in curated_ids]
    require(len(held_out) >= 1,
        "benchmark has no held-out systems"
        + " (would only measure training set)")

    results = []
    for test in manifest.tests:
        iter_iso = run_imago(
            test, label="isolated").iterations
        iter_def = run_imago(
            test, label="default_solid").iterations
        pct = ((iter_iso - iter_def)
               / iter_iso * 100.0)
        results.append({
            "test_id"       : test.id,
            "reference_id"  : test.reference_id,
            "iter_isolated" : iter_iso,
            "iter_default"  : iter_def,
            "pct_reduction" : pct,
            "held_out"      : (test in held_out)})

    mean_pct = mean(r["pct_reduction"]
        for r in results)
    held_out_mean_pct = mean(r["pct_reduction"]
        for r in results if r["held_out"])

    verdict = "PASS" if mean_pct >= 20.0 else "FAIL"

    write_report(
        "share/curation/bench_report.md",
        per_system        = results,
        mean_pct          = mean_pct,
        held_out_mean_pct = held_out_mean_pct,
        verdict           = verdict)

    exit(0 if verdict == "PASS" else 1)
```

### 11.6 Manifest authoring (DESIGN 5.7)

The *write* side of the schema library
(`curation_manifest.py`); `load_manifest_v2` (11.4) is the
*read* side.  `cod_fish` prints a complete manifest (or, with
`--sketch-only`, bare `[[reference_solid]]` stubs), the curator
collects them, and `expand_manifest` completes a sketch into a
manifest the writer serializes.  The writer emits human-readable
TOML -- shortest round-trippable floats, inline `sub_spec`
tables in their authored order, and every optional customization
field only when it is set.

**The round-trip is the writer's contract**: whatever
`load_manifest_v2` reads, `format_manifest` must write back, or a
curator's setting silently disappears on the next authoring pass
(DESIGN 5.7).  This is why the two *shared* blocks are emitted
whenever they are non-empty, and why the reader's four top-level
concerns -- `schema_version`, `[characterization]`, `[defaults]`,
`[harvest]` -- all have a writing counterpart below.

Two things the writer deliberately does NOT emit:

- **A run setting a solid did not name.**  A solid carries a run
  setting only to *override* `[defaults]`; an unset one is
  `None` and is left out, so the file stays sparse and the
  inheritance is visible rather than copied onto every solid.
- **The `preferred` flag, anywhere.**  Preference is
  *structural*, recovered from the block a fingerprint
  declaration lands in: `true` for a `[characterization]`
  record, `false` for a per-entry one (DESIGN 5.7, rule 11).
  Writing it as a key would let it contradict its own position.

```
function format_manifest(manifest):
    lines = ["schema_version = " + manifest.schema_version]

    # The database-wide preferred recipe: one sub_spec per method,
    # applied to every harvested environment (rule 11).  Each
    # method is its own sub-table; no preferred key is written.
    if manifest.characterization:
        emit "[characterization]"
        for fp in manifest.characterization:
            emit "[[characterization.fingerprint]]"
            emit method, sub_spec (inline table)

    # The shared RUN settings, inherited by every solid that does
    # not override them.  Emitted in RUN_SETTING_KEYS order (not
    # dict order) so the file is stable across authoring tools.
    if manifest.defaults:
        emit "[defaults]"
        for key in RUN_SETTING_KEYS:
            if key in manifest.defaults:
                emit key = manifest.defaults[key]

    # The shared HARVEST settings -- how finished runs are read
    # BACK, as opposed to how they are run (DESIGN 5.7).  Emitted
    # on the same when-non-empty rule as [defaults]: silent when
    # the manifest leans on the built-in default, faithful when a
    # curator wrote one down.  Omitting this block would drop an
    # authored kpoint_convergence_threshold on every rewrite.
    if manifest.harvest:
        emit "[harvest]"
        # The scalar settings first; the kpoint_climb sub-table is a
        #   dict, emitted as its own [harvest.kpoint_climb] block.
        for key in HARVEST_SETTING_KEYS:
            if key == "kpoint_climb":
                continue
            if key in manifest.harvest:
                emit key = manifest.harvest[key]
        climb = manifest.harvest.get("kpoint_climb", {})
        if climb:
            emit "[harvest.kpoint_climb]"
            for knob in KPOINT_CLIMB_KEYS:
                if knob in climb:
                    emit knob = climb[knob]

    for solid in manifest.reference_solids:
        emit "[[reference_solid]]"
        emit reference_id, system_type
        # Exactly one structure source is set (rule 4).
        if solid.cod_id is not None:
            emit cod_id, cod_revision
        else:
            emit structure_path
        # A table's scalars (and the inline kpoint_spec) precede
        # its sub-tables, as TOML requires.  Each run setting is
        # emitted ONLY as an override of [defaults].
        for key in RUN_SETTING_KEYS:
            if getattr(solid, key) is not None:
                emit key = getattr(solid, key)
        if solid.source_description is not None:
            emit source_description

        for entry in solid.entries:
            emit "[[reference_solid.entry]]"
            # Every customization field is optional (DESIGN 5.2.2):
            # emit each only when set, and `default` only when
            # true, since an absent flag reads as false.
            if entry.element     is not None: emit element
            if entry.atom_site   is not None: emit atom_site
            if entry.default:                 emit default
            if entry.description is not None: emit description
            if entry.label       is not None: emit label

            # A per-entry fingerprint is a RARE override: an extra,
            # NON-preferred sub_spec harvested for this environment
            # alongside the database-wide preferred recipe above.
            for fp in entry.fingerprints:
                emit "[[reference_solid.entry.fingerprint]]"
                emit method, sub_spec (inline table)

    return join(lines, "\n") + "\n"   # floats via shortest repr

function write_manifest(manifest, path):
    write_file(path, format_manifest(manifest))
```

The four run-setting values an authoring tool writes into a
fresh `[defaults]` block live here, in the schema library, so
that `cod_fish` and `expand_manifest` cannot drift apart.  Note
what these are NOT: the loader requires every run setting to
resolve from the solid or from `[defaults]` (11.4, rule 2), so
none of them is a resolve-time fallback.  They are the answer to
"what should a newly authored manifest say?", and a manifest
that says something else is honoured as written.

```
# The full basis.  These are reference-quality potentials that
#   every later calculation starts from, so the producer pays
#   for the larger basis once rather than have every consumer
#   inherit a minimal-basis starting guess.
DEFAULT_BASIS = "fb"

# The Wigner interpolation functional -- Imago's own default,
#   and a predictor sub-model selector, so a database built on
#   anything else cannot inform a run built on this.
DEFAULT_FUNCTIONAL = "wigner"

# Linear tetrahedral Brillouin-zone integration, with the
#   Bloechl correction (DESIGN 5.7, 1.6).  The default because
#   the producer must choose an integration scheme BEFORE it
#   knows whether the solid is a metal -- that is what the
#   k-point ladder discovers, often several rungs up -- and this
#   is the choice that is safe under both answers.  In a metal
#   the tetrahedron method varies the occupied volume
#   continuously as the mesh refines, where unsmeared Gaussian
#   integration moves whole states across the Fermi level and
#   rattles the energy by amounts that do not shrink with the
#   mesh spacing.  In a gapped system it costs nothing at all:
#   every tetrahedron is wholly occupied or wholly empty, the
#   Bloechl weights reduce to a quarter per corner, and the
#   result is the Gaussian answer exactly (measured on
#   si_fd-3m_227_2001 at mesh 6-6-6).
#
# This default is the GROUND STATE's.  A core-level
#   spectroscopy run names "gaussian" and gets it, because the
#   core-hole correction is not written against the tetrahedral
#   occupation array (DESIGN 5.7).
DEFAULT_KPOINT_INTEGRATION = "linear-tetrahedral"

# The SCF self-consistency threshold.  Distinct from the k-point
#   flatness tolerance in 11.4: this one governs a SINGLE run's
#   iteration to its own fixed point, and it IS part of the
#   cache key.
DEFAULT_SCF_THRESHOLD = 1.0e-6

function default_run_settings():
    # kpoint_spec is EMPTY on purpose: the producer predicts a
    #   starting density and verifies it by climbing the ladder,
    #   so a pinned density would only override
    #   predict-then-verify.  cell is emitted even though it is
    #   exempt from rule 2 -- an authoring tool writing a fresh
    #   file should say what it chose rather than leave a reader
    #   to know the built-in.
    return { "basis"              : DEFAULT_BASIS,
             "functional"         : DEFAULT_FUNCTIONAL,
             "kpoint_integration" : DEFAULT_KPOINT_INTEGRATION,
             "kpoint_spec"        : {},
             "scf_threshold"      : DEFAULT_SCF_THRESHOLD,
             "cell"               : DEFAULT_CELL }
```

`cod_fish` writes the common COD case straight through this
emitter: it pairs the shared-library `default_characterization()`
and `default_run_settings()` with one sparse solid per pinned
structure, so `cod_fish pin <ids> > manifest.toml` is the whole
authoring step.  Neither authoring tool populates `harvest`, so
neither emits a `[harvest]` block -- the harvest setting has a
built-in default and the resolved value is recorded on every
guidance entry (DESIGN 5.7), so leaving it unwritten loses
nothing.  A curator who wants the tolerance visible and editable
adds the block by hand, and the writer above preserves it.

`expand_manifest.py` -- the sketch-to-manifest authoring tool --
has two modes.  Both hoist the shared settings into `[defaults]`
ONCE and leave each solid sparse, rather than stamping five
settings onto every solid; the two helpers below are what make
that split explicit.

```
function shared_defaults(basis, functional,
                         kpoint_integration, scf_threshold):
    # The top-level [defaults] run settings, emitted once.
    # kpoint_spec is left EMPTY: the producer predicts a starting
    # density and verifies it by a convergence sweep, so pinning a
    # density here would only override predict-then-verify.
    # system_type is NOT among these -- it is an intrinsic property
    # of a structure, named per solid, never a shared run setting.
    return { "basis"              : basis,
             "functional"         : functional,
             "kpoint_integration" : kpoint_integration,
             "kpoint_spec"        : {},
             "scf_threshold"      : scf_threshold }

function sparse_solid(source, system_type, entries = []):
    # Copy the sketch stub's identity and structure source, resolve
    # its system_type, and leave EVERY run setting unset (None) so
    # it inherits [defaults].  The source_description cod_fish read
    # from the CIF rides along, so the harvest can compose each
    # environment's auto-description from it.
    return ReferenceSolid(
        reference_id       = source.reference_id,
        system_type        = system_type,
        cod_id             = source.cod_id,
        cod_revision       = source.cod_revision,
        structure_path     = source.structure_path,
        source_description = source.source_description,
        entries            = copy(entries))
        # basis / functional / kpoint_integration / kpoint_spec /
        #   scf_threshold all stay None -> inherited.

function build_mechanical(sources, **shared):
    # Carry every sketched structure through as a sparse solid,
    # leaving the customizations empty.  The CLI appends a
    # commented fill-in template.  The result is structurally
    # valid and loadable -- the recipe block is present and every
    # run setting resolves from [defaults] (rule 2) -- it simply
    # harvests every environment with auto-composed descriptions
    # until customizations are added by hand.
    solids = [sparse_solid(s, s.system_type or system_type_default)
              for s in sources]
    return CurationManifest(schema_version = 2,
        characterization = default_characterization(),
        defaults         = shared_defaults(**shared),
        reference_solids = solids)

function build_interactive(sources, ask, **shared):
    # Walk the curator through completing each structure.
    # ask(prompt, default) returns the reply or the default;
    # injecting it (rather than calling input) keeps the flow
    # testable with scripted answers.
    confirm the shared defaults once via ask
    default_elements = empty_set     # elements already given a
                                     #   default entry (rule 7)
    solids = []
    for source in sources:
        announce the structure (a printed header, not a prompt)
        system_type = ask(..., source.system_type or default)
        entries = []
        # The FIRST "add an entry?" defaults to yes; once one is
        # added, "add another" defaults to no, so pressing Enter
        # ends this structure rather than looping forever.
        while ask_yes_no("add an entry?", default = (entries == [])):
            # Default the element to the structure's composition
            # (the next not-yet-entered element, so a one-element
            # solid auto-fills) and the description to the
            # CIF-derived hint -- both from cod_fish.
            element, atom_site, description,
                label  <- ask(...)        # blank label -> None
            # Default to yes only until this element has its one
            # default entry, so accepting the defaults yields
            # exactly one default per element (rule 7).
            is_default = ask_yes_no("  default for this element?",
                default = element not in default_elements)
            if is_default: default_elements.add(element)
            # NO per-entry fingerprints are authored here: the
            # preferred recipe is the database-wide
            # [characterization] block, set once below, and a
            # per-entry declaration is a rare hand-edit override.
            entries.append(ReferenceEntry(element, atom_site,
                is_default, description, label, fingerprints = []))
        solids.append(sparse_solid(source, system_type, entries))
    return CurationManifest(schema_version = 2,
        characterization = default_characterization(),
        defaults         = shared_defaults(**shared),
        reference_solids = solids)
```

## 12. imago.py Callable API (DESIGN 6.1)

The refactor of `imago.py` from a command-line-only
driver into a callable Python API, per DESIGN 6.1.  Five
pieces: the result object and status enum (12.1); the
single-source-of-truth output-name table that both the
output writer and the result collector consult (12.2,
which resolves the DESIGN 6.1.6 open detail); the two
entry points plus the thin CLI wrapper (12.3); the
private run core with its lock lifecycle, cwd-restore
discipline, and returned-status-vs-raised-error boundary
(12.4); and the harvest helpers that read the result
fields off the settled output files (12.5).

The behavior of an actual run is unchanged from today's
`main()`; the structural change is that the orchestration
becomes a function that *returns a value* and reports
failure by *status* rather than calling `sys.exit`, so a
long-lived caller (a kaleidoscope worker, §13) can drive
many runs in one process.

### 12.1 Result object and status (DESIGN 6.1.2)

```
enum RunStatus:
    CONVERGED       # ran; SCF reached its threshold
    NOT_CONVERGED   # ran cleanly; hit the iteration ceiling
    FAILED          # did not complete (abort / missing
                    #   success file / missing input)
    SKIPPED         # nothing to do; checkpoint found the
                    #   requested work already complete

dataclass ImagoResult:
    status            : RunStatus
    success           : bool      # status == CONVERGED
    run_dir           : str       # absolute project home
    temp_dir          : str       # absolute IMAGO_TEMP mirror
    scf_iterations    : int|None  # None when no SCF ran
    converged         : bool      # SCF met threshold
    reused_checkpoint : bool      # work was short-circuited
    total_energy      : float|None  # Hartree, when available
    outputs           : dict      # logical key -> abs path
    job               : JobIdentity  # edge, job_name,
                                     #   basis_scf, basis_pscf
    runtime_seconds   : float
    message           : str
```

`success` is a derived convenience (`status ==
CONVERGED`).  `outputs["scfV"]` is the converged potential
the database producer harvests (DESIGN 6.1.1); it is only
trustworthy when `status == CONVERGED`.

A contract-level failure raises instead of returning:

```
class ImagoError(Exception):
    # Raised for programmer/environment faults that no
    # per-job retry can fix: $IMAGO_RC / $IMAGO_TEMP /
    # $IMAGO_BIN unset; run_dir missing or holding no
    # inputs; the per-run-dir lock already held by another
    # process.  Run-level failures (non-convergence, a
    # Fortran abort, a missing-at-run-time input) are NOT
    # raised -- they come back as a FAILED / NOT_CONVERGED
    # ImagoResult so a flight can record-and-continue
    # (VISION Principle 10).
    pass
```

### 12.2 Output-name table (resolves DESIGN 6.1.6)

The single source of truth for the project-home output
filenames.  Today `manage_output` *moves* `fort.*` files
to these names inline; factoring the names into one table
means the writer and the API's result collector (12.5)
cannot drift apart.  Both consult `project_home_outputs`.

The names reuse the existing `FileNames` tokens (`scfV`,
`enrg`, `iter`, the property tags `dos`/`bond`/...) and
the `edge_`, `basis` tags `manage_output` already builds.

```
function project_home_outputs(settings):
    # Returns {logical_key: filename}.  Filenames are
    # relative to run_dir; the collector makes them
    # absolute and keeps only those that exist (some are
    # conditional on spin or job type).
    edge  = settings.edge
    jn    = settings.job_name
    jid   = settings.job_id

    # The basis tag mirrors manage_output exactly.
    if jid < 200:        basis = "-" + settings.basis_scf
    elif jid < 300:      basis = "-" + settings.basis_pscf
    else:                basis = "-fb"

    out = {}

    # --- SCF-always block (the producer's keys) ---
    # Written whenever an SCF ran: job_id < 200, or a
    # post-SCF job whose basis_scf is not "no".
    if jid < 200 or settings.basis_scf != "no":
        out["scfV"]      = f"{edge}_scfV{basis}.dat"   # fort.8
        out["energy"]    = f"{edge}_enrg{basis}.dat"   # fort.14
        out["iteration"] = f"{edge}_iter{basis}.dat"   # fort.7
        # iterTDOS plot is emitted only if fort.1000 existed.
        out["iterTDOS"]  = f"{edge}_{jn}{basis}.iterTDOS.plot"

    # --- All-tasks block ---
    out["out"] = f"{edge}_{jn}{basis}.out"             # fort.20

    # --- Property-specific block, by job_id % 100 ---
    # Each property contributes the file family its
    # _manage_<prop>_output helper writes to the project
    # home.  Spin-polarized runs add ".up"/".dn" variants
    # of the same keys; the collector keeps whichever
    # exist.  Tag = the FileNames token for the property.
    prop = jid % 100
    out.update(property_outputs(prop, edge, jn, basis))
    return out
```

```
function property_outputs(prop, edge, jn, basis):
    # tag -> (key family).  Compact transcription of the
    # _manage_*_output destinations; ".t"/".p" = total/
    # partial, ".up"/".dn" = spin, suffixes are the
    # quantity tags (".cond", ".eps1", ...).
    #   1  dos   : "dos"   (.t/.p tot+partial, .loci)
    #   2  bond  : "bond"  (.raw, .3c three-center)
    #   3  dimo  : "dimo"  (.t total moment)
    #   4  optc  : "optc"  (.t + .p partial; .cond,
    #               .eps1, .eps2, .elf, .nref, .kext,
    #               .aabs, .Rref, .eps1i families)
    #   5  pacs  : "pacs"  (.plot)
    #   6  nlop  : "optc"  (.chi1, .chi2)
    #   7  sige  : "sige"  (.cond)
    #   8  sybd  : "sybd"  (.plot) + "vdim" (.raw)
    #   9  force : "force" (.dat)
    #  10  field : "field" (.prof profiles, .rho, .xdmf3)
    #  11  mtop  : "mtop"  (.t total)
    # job_id == 311 also adds loen : "loen" (.plot)
    # Build {key: filename} for the matching tag from the
    # FileNames tokens, exactly as the helper would name
    # them.  (Implementation mirrors the helper bodies;
    # the table above is the authoritative key set.)
    ...
```

The producer (the first client) reads only `scfV`,
`energy`, and `iteration`; the property keys exist so
later clients (DOS sweeps, bond-order flights) reach
their outputs through the same contract.

### 12.3 Entry points and CLI wrapper (DESIGN 6.1.3)

Two API entry points and the CLI, all funneling into the
private core of 12.4.

```
function run_prepared(run_dir, settings = None):
    # Prepared-directory mode: run_dir already holds the
    # staged inputs (imago.dat, structure.dat, scfV.dat,
    # kp files).  No makeinput call.
    if settings is None:
        settings = ScriptSettings.from_options({})  # rc
                                          # defaults only
    require_contract(is_dir(run_dir),
        "run_dir does not exist: " + run_dir)
    return _run_core(run_dir, settings)
```

```
function run_structure(structure, options, run_dir,
                       settings = None):
    # Structure-and-options mode: build the run directory
    # with makeinput first, then run it.  `structure` is a
    # path to an imago.skl (a StructureControl handle is
    # deferred to D12/C64; see DESIGN 6.1.6).  The build API
    # is §14; run_prepared is the run entry it chains into.
    import makeinput              # local: imago.py imports
                                  #   without makeinput's env
    if settings is None:
        settings = ScriptSettings.from_options(options)
    makeinput.build_run_dir(structure, options, run_dir)
    return run_prepared(run_dir, settings = settings)
```

```
function cli_main(argv):
    # The thin wrapper: the ONLY layer that touches argv
    # or exits the process.
    # 1. Parse argv into run options (today's argparse
    #    surface + reconcile logic, unchanged in meaning).
    settings = ScriptSettings.from_command_line(argv)
    # 2. Pick the entry mode.  A bare `imago ...` runs the
    #    current working directory as a prepared dir --
    #    today's only behavior.
    try:
        result = run_prepared(getcwd(), settings)
    except ImagoError as e:
        log_runtime(e.message)
        return 1
    # 3. Translate the result into an exit code.
    if result.status in (CONVERGED, SKIPPED):
        return 0
    log_runtime(result.message)
    return 1   # NOT_CONVERGED or FAILED
```

`ScriptSettings` is split so argv is no longer mandatory
(DESIGN 6.1.3): `from_command_line(argv)` keeps today's
behavior (argparse -> `reconcile`), while
`from_options(mapping)` builds the same reconciled
settings from a plain dict with no argv and no
`command`-file side effect.  Both share the existing
`reconcile()`; only the source of the `args` namespace
differs.

### 12.4 The private run core (DESIGN 6.1.4, 6.1.5)

One core performs today's `main()` sequence, but
returns an `ImagoResult` and is reentrant.

```
function _run_core(run_dir, settings):
    start_clock = now()
    original_cwd = getcwd()           # for the finally
    temp = mirror_under_imago_temp(run_dir)  # get_temp_dir
    fn = FileNames()

    # Contract checks raise (not return); they are
    # environment/programmer faults (DESIGN 6.1.2).
    require_contract(env("IMAGO_RC") and env("IMAGO_TEMP")
                     and env("IMAGO_BIN"),
        "Imago environment not configured")

    makedirs(temp, exist_ok = True)
    lock_path = join(temp, fn.imago_lock)

    # Per-run-dir lock.  Because temp mirrors run_dir, two
    # different run dirs take two different locks, so a
    # flight of parallel runs never collides (DESIGN
    # 6.1.5).  An already-held lock is a contract fault
    # in API mode -> raise.
    if exists(lock_path):
        raise ImagoError(
            "lock already held in " + temp
            + " (another run owns this directory)")
    write_lock(lock_path)

    try:
        chdir(temp)                   # cwd is a resource

        # Within-run-dir checkpoint assessment (DESIGN
        # 6.1.5).  Reads the SAME completed-calculation
        # markers the current script/Fortran already use;
        # this surfaces their state, it does not redesign
        # the mechanism.
        ckpt = assess_checkpoint(temp, run_dir, settings)
        if ckpt == COMPLETE:
            # All requested work already done: short-
            # circuit without invoking the binary.
            return _build_result(
                SKIPPED, run_dir, temp, settings,
                reused = True,
                seconds = now() - start_clock,
                message = "checkpoint: already complete")

        # Stage inputs, run the binary + immediate
        # secondary jobs (SYBD post-pass, optical KK),
        # exactly as today.  manage_input + execute mirror
        # the current flow; execute returns whether the
        # fort.2 success file appeared.
        manage_input(settings, fn, run_dir, temp)
        ran_ok = execute_program(build_job_clp(settings),
                                 settings, fn, temp)

        if not ran_ok:
            # Fortran aborted / no success file: a run-
            # level failure, RETURNED not raised.
            return _build_result(
                FAILED, run_dir, temp, settings,
                reused = (ckpt == PARTIAL),
                seconds = now() - start_clock,
                message = "Fortran success file missing")

        # Collect outputs into run_dir (the writer also
        # consults project_home_outputs, 12.2) and build
        # the result by harvesting the settled files.
        manage_output(settings, fn, run_dir)
        return _harvest_result(
            run_dir, temp, settings,
            reused  = (ckpt == PARTIAL),
            seconds = now() - start_clock)

    except ImagoError:
        raise                          # contract fault
    except Exception as e:
        # Unexpected: report as FAILED, do not kill the
        # caller's process.
        return _build_result(
            FAILED, run_dir, temp, settings,
            reused = False,
            seconds = now() - start_clock,
            message = "unexpected error: " + str(e))
    finally:
        # Always release the lock and restore cwd, even on
        # failure -- the single most important reentrancy
        # difference from the one-shot CLI (DESIGN 6.1.4).
        remove_if_exists(lock_path)
        chdir(original_cwd)
```

### 12.5 Result harvesting (resolves DESIGN 6.1.6)

`_harvest_result` reads the result fields off the settled
output files (the robust default chosen in DESIGN 6.1.6,
over scraping stdout).  The convergence verdict and the
total energy both come from a single read of the
iteration file's last line, which closes the DESIGN 6.1.6
open detail with no Fortran change.

**The iteration file's shape matters here.**  It is
`fort.7` with one header line, written only when the file
is first created (`safe_append`'s full-copy branch);
because `safe_append`'s `skip_lines` is 1-based
(`tail -n +N`), reruns append `fort.7` from line 2 on, so
they contribute data rows with no extra header.  Two
consequences: (1) there is exactly one header line, ever;
(2) successive SCF runs in the same run directory append
their cycles back-to-back, so the file may hold several
runs' worth of rows.  The last data row is therefore the
most recent SCF cycle of the most recent run -- exactly
the row to inspect.

```
function _harvest_result(run_dir, temp, settings, reused,
                         seconds):
    names   = project_home_outputs(settings)      # 12.2
    outputs = { key: join(run_dir, fname)
                for key, fname in names.items()
                if exists(join(run_dir, fname)) }

    iters = energy = mag = gap_ev = gap_kind = None
    scf_threshold = None
    conv = False
    if "iteration" in outputs:
        # One read of the last data row yields several fields at
        # once.  The row is a fixed 8-column form (all 1-based):
        #   1 iter        4 convergence     7 gap (hartree)
        #   5 total_E     6 magnetic moment 8 gap-kind code
        # Columns 6-8 are length-gated so a shorter iteration
        # file (a property pass, or a run before the gap columns)
        # still parses cleanly.
        row = last_data_row(outputs["iteration"])
        scf_threshold = read_scf_threshold(
                            join(run_dir, "imago.dat"))
        # Column 4 is the SCF convergence metric; converged iff
        # below the imago.dat criterion.  Column 5 is the last
        # iteration's total energy.  Column 1 is a per-run cycle
        # counter that resets each SCF invocation, so it is THIS
        # run's iteration count even though reruns append rows.
        conv   = (column(row, 4) < scf_threshold)
        energy = column(row, 5)
        iters  = int(column(row, 1))
        if len(row) >= 6:
            mag = column(row, 6)
        if len(row) >= 8:
            gap_ev   = column(row, 7) * HARTREE_TO_EV
            gap_kind = GAP_KIND_BY_CODE[int(column(row, 8))]

    # Resolved mesh: two labeled lines in the SCF output (4d.5).
    # Absent for an explicit-list run (style 0) or an older
    # binary, in which case both stay None (DESIGN 6.1.2), which
    # leaves the k-density guard inert (15.7).
    kpoint_mesh  = None
    kpoint_count = None
    if "out" in outputs:
        kpoint_mesh  = read_labeled_ints(outputs["out"],
                          "RESOLVED_KP_MESH", 3)
        kpoint_count = read_labeled_int(outputs["out"],
                          "RESOLVED_KP_COUNT")

    status = CONVERGED if conv else NOT_CONVERGED
    # No SCF at all (e.g. -scf no post-SCF property run):
    # there is nothing to converge; treat a clean run as
    # CONVERGED so success reflects "ran as asked".
    if iters is None:
        status = CONVERGED

    return ImagoResult(
        status = status, success = (status == CONVERGED),
        run_dir = run_dir, temp_dir = temp,
        scf_iterations = iters, converged = conv,
        reused_checkpoint = reused, total_energy = energy,
        total_magnetization = mag, gap_ev = gap_ev,
        gap_kind = gap_kind, scf_threshold = scf_threshold,
        kpoint_mesh = kpoint_mesh, kpoint_count = kpoint_count,
        outputs = outputs, job = job_identity(settings),
        runtime_seconds = seconds,
        message = status.name)


function read_labeled_ints(path, label, n):
    # imago's label/value convention (4d.5): the tag is alone on
    # a line and its value is on the NEXT line.  Return the first
    # n integers on the line after the first line whose stripped
    # text is exactly `label`, or None if the label is absent.
    # The exact whole-line match is robust to surrounding output
    # and to the label appearing as a substring elsewhere.
    lines = read_lines(path)
    for i in range(len(lines)):
        if strip(lines[i]) == label:
            tokens = split(lines[i + 1])
            return [int(tokens[k]) for k in range(n)]
    return None

function read_labeled_int(path, label):
    result = read_labeled_ints(path, label, 1)
    return None if result is None else result[0]
```

```
function read_scf_threshold(imago_dat_path):
    # The SCF convergence criterion is the value on the
    # line immediately following the "CONVERGENCE_TEST"
    # label in imago.dat (the run's own input), so the
    # verdict uses the same criterion the run was held to.
    lines = read_lines(imago_dat_path)
    i = index_of_line(lines, "CONVERGENCE_TEST")
    return float(lines[i + 1])
```

This resolves the DESIGN 6.1.6 open detail in full.  The
convergence verdict, the total energy, and the per-run
iteration count all come from one read of the iteration
file's last data row (columns 4, 5, and 1 respectively),
compared against the `CONVERGENCE_TEST` criterion in
`imago.dat` -- no new Fortran signal is needed, because
the verdict reuses the convergence metric the SCF already
writes per cycle.  Everything else in §12 is a faithful
restructuring of behavior that `imago.py` already has.

## 13. kaleidoscope Flight Dispatch (DESIGN 6.2)

The Parsl-based dispatcher that drives a *set* of Imago
calculations, per DESIGN 6.2.  It builds on §12: the
default wingbeat calls the §12 callable API and persists
its `ImagoResult`.  The pieces, helpers first then the
driver: the data model and `flight.toml` (13.1); the
wingbeat protocol and the Imago wingbeat (13.2); the
workspace paths, id rules, and `status.toml` (13.3); the
cache hit-test (13.4); the dispatch driver with
complete-and-report (13.5); and the report plus the
client-side harvest handoff (13.6).

The governing rule (VISION Principle 9): kaleidoscope is
domain-agnostic.  It dispatches, tracks, and caches; it
never interprets what a run computed.  The wingbeat's
`detail` string is recorded verbatim and never parsed;
all domain harvest is client-side (13.6).

### 13.1 Data model and flight.toml (DESIGN 6.2.1)

```
dataclass KeyFile:
    path   : str     # where the staged copy sits, RELATIVE to
                     #   the unit's directory.  Usually carries a
                     #   directory part -- the producer declares
                     #   `inputs/structure.dat` -- and the same
                     #   path is joined onto the prepare dir and
                     #   the run dir alike (13.4 / DESIGN 6.2.5)
    source : str     # the current input byte-compared against
                     #   that staged copy

dataclass KeyFields:
    scalars : dict   # verbatim-compared identity fields; for the
                     #   producer just {converg} (DESIGN 6.2.5).
                     #   The engine build is NOT among them: it is
                     #   recorded per run and never compared
    files   : list   # KeyFile entries to byte-compare (path +
                     #   source); naming both keeps the core from
                     #   guessing how inputs map onto staged files

dataclass CalcUnit:
    id          : str          # stable per-structure key
    calc        : tuple[str,...]  # per-axis directory
                               #   components (DESIGN
                               #   6.2.1); () = no second
                               #   level; one element per
                               #   varied sweep axis
    structure   : str          # path to an imago.skl
    prepared_dir: str | None   # per-unit staging dir the
                               #   driver's prepare step fills
                               #   with the built inputs
                               #   (structure.dat, imago.dat...):
                               #   the structure.dat KeyFile's
                               #   source points here for the
                               #   hit-test, and the wingbeat
                               #   commits it into the run dir on
                               #   a miss (DESIGN 6.2.5, Model A).
                               #   None until prepared.
    options     : dict         # makeinput options
    wingbeat    : str | None   # wingbeat name; None -> the
                               #   flight default
    kind        : str          # run role (DESIGN 6.2.9): a
                               #   short label the core stores
                               #   and round-trips but never
                               #   interprets.  Default
                               #   "convergence"; each
                               #   harvester reads only the
                               #   kinds it understands (e.g.
                               #   "fingerprint" for loen runs)
    key_fields  : KeyFields    # client-declared identity
    record      : dict         # free-form facts ABOUT the run
                               #   that are not inputs TO it --
                               #   the engine build identity is
                               #   the standing case.  Copied
                               #   verbatim into status.toml's
                               #   [record] at launch and never
                               #   compared, interpreted, or put
                               #   in the report (DESIGN 6.2.4/
                               #   6.2.5).  Default {}

dataclass Flight:
    root             : str     # workspace root directory
    units            : list    # list[CalcUnit]
    default_wingbeat : str     # wingbeat name for None units
    parsl_config     : object  # a Parsl Config (deployment)
    sweep            : SweepRecord | None  # varied/fixed axes
                               #   when built by the predict-
                               #   then-verify helper (DESIGN
                               #   6.2.8); None otherwise
    on_outcome       : callable | None  # per-unit callback
    metadata         : dict    # opaque per-key tables the
                               #   core round-trips verbatim
                               #   as [flight.<key>] and never
                               #   reads (Principle 9); 6.2.8/
                               #   6.2.9 stash the per-structure
                               #   predictions mapping here

dataclass SweepRecord:
    varied_axes : tuple[str,...]  # axis names, in the order
                               #   they appear at each level
                               #   of CalcUnit.calc
    fixed_axes  : dict         # axis -> value for axes held
                               #   constant across the flight
```

```
function serialize_flight(flight):
    # Write flight.toml: the authoritative record of
    # WHAT was asked for, separate from each run's
    # status.toml record of WHAT HAPPENED (13.3).  A resume
    # (13.5) reads the units back from here.  The optional
    # [flight.sweep] block is emitted only when the
    # flight was built by the predict-then-verify helper
    # (DESIGN 6.2.1/6.2.8); each unit's calc tuple serializes
    # as a TOML array of directory-component strings.
    record = { root = flight.root,
               default_wingbeat = flight.default_wingbeat,
               units = [ as_dict(u) for u in flight.units ] }
    if flight.sweep is not None:
        record["sweep"] = as_dict(flight.sweep)
    # Each metadata[key] becomes a verbatim [flight.<key>]
    # table; the core never reads the contents (Principle 9).
    # The 6.2.8/6.2.9 builder stashes the predictions mapping
    # this way, emitted as [flight.predictions.<id>] sub-tables.
    for key, table in flight.metadata.items():
        record[key] = table
    write_toml(join(flight.root, "flight.toml"), record)
```

### 13.2 Wingbeat protocol and ImagoWingbeat (DESIGN 6.2.2)

The wingbeat is the seam (Principle 8) between dispatch and
execution.  It returns a *domain-agnostic* outcome.

```
dataclass WingbeatOutcome:
    ok              : bool     # did the unit COMPLETE
                               #   (not "succeed
                               #   scientifically")
    detail          : str      # opaque string the wingbeat
                               #   chooses; recorded, never
                               #   interpreted by kaleido-
                               #   scope (e.g. "converged")
    runtime_seconds : float
    message         : str

protocol Wingbeat:
    function run(unit, wingbeat_dir) -> WingbeatOutcome
```

```
class ImagoWingbeat implements Wingbeat:
    function run(unit, wingbeat_dir):
        # Default wingbeat: drive the §12 API.  First make the run
        # dir hold runnable inputs (stage_inputs, below), then run
        # them with the unit's imago-side settings.  Those settings
        # (job / edge / scf_basis) are RUNTIME options NOT baked
        # into the staged imago.dat (DESIGN 6.2.10), so they must
        # be re-applied on EVERY launch.  The job type and the SCF
        # suppression live only in these settings, so if they are
        # dropped imago no longer sees the unit's `-loen -scf no`
        # request and falls back to its DEFAULT job, a ground-state
        # SCF.  (`-loen -scf no` never runs an SCF itself; the
        # unwanted SCF is purely the dropped-settings fallback --
        # the "SCF after loen" the seed run hit.)
        stage_inputs(unit, wingbeat_dir)
        imago_opts = { key : value
                       for key, value in unit.options.items()
                       if key in imago.OPTION_KEYS }
        settings = ScriptSettings.from_options(imago_opts)
        result = imago.run_prepared(
                     wingbeat_dir, settings = settings)

        # Persist the §12.1 ImagoResult for the client to
        # reload (13.6).  kaleidoscope never reads it; it is
        # the wingbeat -> client handoff, kept domain-side.
        #
        # One RECORDED fact rides along with the measured ones: the
        # build identity out of unit.record (DESIGN 6.2.4), written
        # as `imago_commit`.  A guidance entry's provenance reads it
        # here (15.7), which keeps that harvest on the three sources
        # it already has and off the core's status.toml.  The
        # engine's own word wins when it has one -- imago does not
        # report its build yet (TODO C84), and when it does this
        # `or` stops preferring the producer's belief without any
        # other change.
        fields = as_dict(result)
        if "imago_commit" not in fields:
            fields["imago_commit"] = unit.record.get(
                                         "imago_commit")
        write_toml(join(wingbeat_dir, "result.toml"), fields)

        # Map the Imago-native status onto the generic
        # outcome.  "Ran" covers CONVERGED / NOT_CONVERGED
        # / SKIPPED; only a hard FAILED is not-ok.  The
        # status name becomes the opaque detail string
        # (e.g. "converged", "not_converged", "skipped").
        ok = result.status in (CONVERGED, NOT_CONVERGED,
                               SKIPPED)
        return WingbeatOutcome(
            ok = ok,
            detail = lower(result.status.name),
            runtime_seconds = result.runtime_seconds,
            message = result.message)
```

```
function stage_inputs(unit, wingbeat_dir):
    # Ensure the run dir holds runnable inputs.  Three cases:
    #   - the driver's prepare step (11.4) already built them into
    #     unit.prepared_dir -> COMMIT that staged copy into the run
    #     dir (the Model-A producer path; DESIGN 6.2.5);
    #   - the run dir already holds a staged imago.dat -- a re-run
    #     of a dir a prior launch built -> nothing to do;
    #   - neither (a client that did not prepare) -> BUILD from the
    #     unit's structure and makeinput-side options.
    # Two buckets, not three (DESIGN 6.2.10): every key in
    # `options` is a real tool input, so what is not an imago key
    # is a makeinput dest and makeinput's strict unknown-key check
    # is once again a pure typo backstop.  Bookkeeping that reaches
    # neither tool rides on unit.record instead of being carried
    # here and dropped again.
    if unit.prepared_dir is not None:
        commit_prepared_inputs(unit.prepared_dir, wingbeat_dir)
    else if not is_prepared(wingbeat_dir):
        mk_opts = { k : v for k, v in unit.options.items()
                    if k not in imago.OPTION_KEYS }
        makeinput.build_run_dir(unit.structure, mk_opts,
                                wingbeat_dir)


function is_prepared(wingbeat_dir):
    # A run dir is 'prepared' when it already holds the primary
    # imago.dat (directly or under inputs/), so run_prepared can
    # run it as-is without a makeinput build.
    return exists(join(wingbeat_dir, "imago.dat")) \
        or exists(join(wingbeat_dir, "inputs", "imago.dat"))
```

```
function commit_prepared_inputs(prepared_dir, wingbeat_dir):
    # Hand the driver-staged inputs (structure.dat, imago.dat,
    # scfV, kp files -- DESIGN 6.2.5) to the run dir so
    # run_prepared finds them.  The staging area is transient
    # (the prepare pass, 11.4, rebuilds it each producer run),
    # so the commit reads from it and never writes back.
    #
    # TWO steps, and the first is the one that is easy to omit.
    # A commit lands on a run dir that MAY ALREADY HOLD a prior
    # calculation's files, and copying the new ones in is not by
    # itself enough to make the old ones stop being read.  Clear
    # first, then copy (DESIGN 6.2.5, "What a commit owes a
    # surviving run directory").
    clear_superseded_root_copies(prepared_dir, wingbeat_dir)

    # Copy every staged entry across, merging into whatever the
    # run dir already holds.  Directories are merged rather than
    # replaced, so the staged inputs/ refreshes the run dir's
    # inputs/ without discarding anything else parked there.
    make_dirs(wingbeat_dir, exist_ok = True)
    for name in list_entries(prepared_dir):
        source = join(prepared_dir, name)
        target = join(wingbeat_dir, name)
        if is_dir(source):
            copy_tree(source, target, merge = True)
        else:
            copy_file(source, target)


function clear_superseded_root_copies(prepared_dir, wingbeat_dir):
    # Remove the run dir's ROOT copy of every name this commit is
    # about to stage under inputs/.
    #
    # Why this exists.  makeinput writes only inputs/; it is
    # imago.py that copies each file up to the run-dir root, on
    # the first run and only when the root copy is absent, and
    # that thereafter reads the ROOT copy in preference to the
    # staged one (DESIGN 6.2.5).  The root copy is a cache of
    # inputs/ that nothing invalidates.  Without this step a
    # commit refreshes inputs/, the root copies keep the previous
    # calculation's contents, and the engine runs the OLD physics
    # while the key file, the run's summary, and the flight
    # report all describe the new.
    #
    # DELETE, do not overwrite.  Overwriting would require this
    # step to know which staged names get flattened to the root
    # and under what names -- knowledge that already lives in
    # imago.py, and that would drift the moment it lived twice.
    # An absent root copy is unambiguous: imago.py's own copy-up
    # refills it from the staged file, leaving exactly one writer
    # of the root copy and one source for its contents.
    #
    # The removal list is exactly the staged inputs/ names, which
    # is what keeps a prior run's OUTPUTS untouched.  A converged
    # potential (gs_scfV-*.dat) is an output name, absent from
    # inputs/, so it survives by construction rather than by a
    # carve-out -- and it should survive, being a starting point
    # every later SCF re-converges (DESIGN 6.2.5).  fort.* units,
    # the intermediate/ link and the logs are likewise not staged
    # names and are likewise left alone.  This is NOT "wipe the
    # run directory": 6.1's within-directory checkpointing and
    # the stored potential both depend on it surviving a commit.
    staged_inputs = join(prepared_dir, makeinput.INPUTS_DIR)
    if not is_dir(staged_inputs):
        return          # a staging dir with no inputs/ stages
                        #   no name whose root copy could be stale
    for name in list_files(staged_inputs):
        root_copy = join(wingbeat_dir, name)
        if is_file(root_copy):
            remove_file(root_copy)
    # Removing nothing is the normal first-run case (a clean run
    # dir has no root copies yet), so an empty pass is silent and
    # a missing file is never an error.
```

An ASE wingbeat (D12) and future adapters implement the
same protocol; the dispatch core (13.5) never changes
when one is added.  `commit_prepared_inputs` is
ImagoWingbeat's own step; another wingbeat stages its
inputs however its toolchain requires.  So is the clearing
it now does first: the root-copy precedence being undone
there is `imago.py`'s, so a wingbeat driving some other
toolchain has no such copies to clear and needs no
equivalent.

### 13.3 Workspace paths, ids, status.toml (DESIGN 6.2.4)

```
function unit_run_dir(flight, unit):
    base = join(flight.root, "wingbeats", unit.id)
    # The optional <calc> level(s) exist only when a
    # structure hosts more than one calculation.  calc is a
    # tuple of per-axis directory components (one level per
    # varied sweep axis, DESIGN 6.2.1), so splat it onto the
    # path; an empty tuple leaves the unit directly in base.
    return join(base, *unit.calc) if unit.calc else base
```

```
function validate_flight(flight):
    # Enforce the id/<calc> scheme of DESIGN 6.2.4 at build
    # time; abort (raise) on any violation, naming the
    # offenders -- a silent rewrite would break the cache
    # hit-test (13.4).
    seen = {}                       # id -> set of calc tuples
    for unit in flight.units:
        require_slug(unit.id)       # lowercased [a-z0-9_-]
        for component in unit.calc:  # each directory level
            require_slug(component)  #   is its own slug
        tag = unit.calc             # a tuple (possibly empty)
        # Derive a default <calc> when an id ends up hosting
        # multiple units but a unit gave no tag (DESIGN
        # 6.2.4): a one-element tuple holding
        # "<job_name>-<basis_scf>" for the Imago wingbeat.
        if tag == () and id_hosts_multiple(flight,
                                           unit.id):
            tag = (derive_calc_tag(unit),)
            unit.calc = tag
        require(unit.id not in seen
                or tag not in seen[unit.id],
            "duplicate run dir for id="
            + unit.id + " calc=" + str(tag))
        seen.setdefault(unit.id, set()).add(tag)
```

```
function require_slug(s):
    # Filesystem-safe and unique-friendly: lowercase,
    # [a-z0-9_-] only.  Reject anything else rather than
    # rewriting it.
    require(matches(s, "^[a-z0-9_-]+$"),
        "id/calc not a slug: " + s)
```

```
function write_status(wingbeat_dir, **fields):
    # One file per run dir, rewritten through the
    # lifecycle.  status is kaleidoscope-owned and generic;
    # convergence rides in `detail`, never in `status`.
    # Omit started_at/finished_at/runtime_seconds until
    # they exist; omit calc when it is the empty tuple.
    #
    # A [record] table already present in the file SURVIVES a
    # rewrite (DESIGN 6.2.4): the lifecycle rewrites the status
    # fields many times, but the record is stamped once at launch
    # and describes the run, so re-reading and carrying it forward
    # is what keeps it from being erased at the first transition.
    prior = read_status(wingbeat_dir)
    if prior is not None and "record" in prior \
            and "record" not in fields:
        fields["record"] = prior["record"]
    write_toml(join(wingbeat_dir, "status.toml"), fields)

function read_status(wingbeat_dir):
    p = join(wingbeat_dir, "status.toml")
    return read_toml(p) if exists(p) else None
```

The five `status` values are `queued`, `running`, `done`,
`failed`, `lost` -- the first four are the unit lifecycle
(`done` iff `WingbeatOutcome.ok`); `lost` is the
kaleidoscope-only category for a Parsl-side disappearance
where no `WingbeatOutcome` came back (13.5).

### 13.4 Cache hit-test (DESIGN 6.2.5)

Mechanism owned by kaleidoscope; the key *fields* are
supplied by the client on each `CalcUnit`.

```
function is_cache_hit(unit, wingbeat_dir):
    # Hit iff the dir exists, its recorded key still
    # matches the unit's current key, AND its status is
    # "done".  Anything else is a miss -> (re)launch.
    if not is_dir(wingbeat_dir):
        return False
    st = read_status(wingbeat_dir)
    if st is None or st["status"] != "done":
        return False
    return cache_key_matches(unit, wingbeat_dir)
```

```
function cache_key_matches(unit, wingbeat_dir):
    saved = read_toml(join(wingbeat_dir, "cache_key.toml"))
    if saved is None:
        return False
    # Scalar fields: verbatim field-by-field compare.
    if saved["scalars"] != unit.key_fields.scalars:
        return False
    # Each key file is checked TWICE, for two different things
    # (DESIGN 6.2.5).  No hashing anywhere -- a developer can diff
    # the files to see why a cache missed (DESIGN 6.2.5 / 5.7).
    #
    # files_byte_equal READS BOTH FILES, every call.  It must not
    # memoize on a stat signature: size and mtime are exactly what
    # a same-size in-place rewrite leaves unchanged, and mtime
    # resolution is coarse enough that two writes a few
    # microseconds apart routinely share one tick -- measured at
    # 169 of 200 on this filesystem.  A memoized comparison then
    # answers "equal" for files that differ, which is a false HIT:
    # the stored result is returned for a calculation whose inputs
    # have changed.  Python's `filecmp.cmp` memoizes exactly this
    # way, so an implementation built on it must clear that cache
    # before use (TODO D23).
    for key_file in unit.key_fields.files:

        # (1) IDENTITY.  Both sides at the SAME relative path,
        # which for the producer is under inputs/ -- the one
        # surface makeinput writes for every unit, whatever that
        # unit's job reads.  Declaring a name the run directory
        # carries only for some job kinds makes the other kinds
        # permanently uncacheable (TODO D23).
        #
        # BOTH sides are guarded, and a file that cannot be read
        # is a MISS rather than an error (DESIGN 6.2.5).  Either
        # may be absent -- a prepare directory reclaimed as
        # scratch, a structure cache that moved, a run directory
        # left half written by a job that died -- and all of those
        # mean the one thing: this unit's identity cannot be
        # established, so re-run it rather than trust it.  Raising
        # here would let a single unreadable file abort a campaign
        # that had already paid for hours of converged rungs.
        staged = join(wingbeat_dir, key_file.path)
        if not exists(staged) or not exists(key_file.source):
            return False
        if not files_byte_equal(key_file.source, staged):
            return False

        # (2) AGREEMENT.  The flattened root copy is the file
        # imago actually reads, so it must not say something the
        # key does not.  ABSENT is not a fault: it means this
        # unit's job does not read that file, and identity has
        # already been settled on inputs/.  PRESENT and
        # disagreeing means the engine would run inputs the key
        # does not describe -- miss, and re-run (DESIGN 6.2.5).
        #
        # A client declaring a bare name (no directory part)
        # lands root_copy on staged itself, so the compare is
        # self against self and the test is a no-op.  That is the
        # intended behaviour, not a case to special-case away.
        root_copy = join(wingbeat_dir, basename(key_file.path))
        if exists(root_copy):
            if not files_byte_equal(root_copy, staged):
                return False
    return True
```

```
function write_cache_key(wingbeat_dir, unit):
    # The identity snapshot, written on launch (13.5).  Only the
    # key-file PATHS are recorded (for inspection); the byte
    # compare reads the staged files themselves, not this list.
    # So changing what the declared paths are called does NOT
    # invalidate a stored unit -- nothing compares this array.
    write_toml(join(wingbeat_dir, "cache_key.toml"),
        { scalars = unit.key_fields.scalars,
          files   = [kf.path for kf in unit.key_fields.files] })
```

Each `KeyFile` names both halves the compare needs -- the
`source` (the current input) and the `path` (where the staged
copy sits, relative to the unit's directory) -- so the core
stays oblivious to how a client's inputs map onto staged files
(DESIGN 6.2.5).  The field is `path` rather than `name`
precisely because it routinely carries a directory part:
the producer declares `inputs/structure.dat` and
`inputs/kp-scf.dat`, and a reader who takes it for a bare
filename will put the compare back on the run-directory root,
which is the D23 defect.  For the producer, the driver's
prepare step (11.4, Phase 1b) points each KeyFile's `source`
at the copy it builds.

### 13.5 Dispatch driver (DESIGN 6.2.3)

One task per unit; per-future exception capture so a single
failure never aborts the flight (Principle 10).  Resuming a
flight is just re-running it: the hit-test skips the `done`
units and re-dispatches the rest -- unless the re-run switch
(`force`) is set, which bypasses the cache so every unit
re-launches (DESIGN 6.2.5).  Before any of that, the driver
prints its **reuse plan**.  The closing counts -- how many units
it will reuse, how many it will run -- always print, because that
is the decision being announced.  The per-unit lines behind them
(reuse or run, and on a reuse the finish time and recorded build
behind the result) are that count's evidence and print only under
`verbose`, per the reporting rule of DESIGN 5.7: the climb calls
`send_off` once per round (4e.5), so an unconditional line per
unit would refill the screen C131 cleared, on the very path it
cleared.  `force` is the ordinary way a re-run is asked for, and
`preview` prints the plan in full -- lines and counts, `verbose`
or not, since reading them one by one is what a preview is for --
and dispatches nothing.

The driver runs each unit through an *executor* -- the seam
that hides where the work actually lands.  A `LocalExecutor`
runs units in the current process (tests, a laptop, the
materialize pre-flight); a `ParslExecutor` dispatches them as
Parsl tasks onto whatever its `Config` describes (a laptop
thread pool, or a cluster).  When the caller pins no
executor, the driver builds one from the flight with
`make_executor` (below): a `ParslExecutor` if the flight carries
a `parsl_config` (13.7), else a `LocalExecutor`.  Both kinds of
future share one contract -- `result()`, which returns the
outcome or re-raises the worker's exception, and `done()`, true
once the result is ready -- so the gather logic below is
identical across them.  A cache hit yields an already-done future
too (`completed_future`), so hits and misses sit uniformly in the
outstanding set.

Dispatch is *two phases, and both are public* (DESIGN 6.2.3): a
`send_off` that launches a chosen set of units and returns one
future per unit without waiting, and a `collect` that resolves a
single future into its terminal status and report entry.  The
one-shot `dispatch` below is exactly `send_off` then `collect`
every future in unit order -- the convenience form, so every
existing caller is unchanged.  A control loop instead calls the
two directly: the climb (4e.5) sends the rungs it has decided,
then uses `collect_next` to take whichever rung lands first and
send that chain's successor at once, never waiting out a whole
batch.  `collect_next` is domain-ignorant -- it polls `done()`
and knows only futures, never k-points -- so the choice of what
to send next stays in the producer (Principle 12).

A client that runs *many* flights under one config builds the
executor ONCE and pins it to every call.  A Parsl `Config`'s
executor is single-use: constructing it starts one coordinator
and its pool of SLURM workers (13.7), and closing it tears them
down, so a second `dispatch()` handed the same config would try
to restart a spent pool.  The producer's climb is exactly this
-- a pre-flight batch and then a stream of rungs sent and
collected one at a time (DESIGN 3.12.5), all under one config --
so it builds one executor with `make_executor(config)`, passes it
to every `send_off`/`collect`, and closes it once after the last.
The whole run then shares one warm pool that grows and shrinks
with demand (DESIGN 6.2.11's pooled shape) instead of rebuilding
it per rung.

```
function reuse_plan(flight, units, force):
    # What the driver is ABOUT to do, decided from local files and
    # nothing else, and computed without touching a thing (DESIGN
    # 6.2.5).  This is what stands in for the automatic staleness
    # guard the cache key no longer applies: the build behind a
    # reused result is REPORTED, so a curator who has since fixed
    # that build can re-run on purpose, rather than SILENTLY
    # COMPARED, which would discard every stored result on every
    # rebuild.  Read-only, so `preview` and the real send share it.
    plan = []
    for unit in units:
        wingbeat_dir = unit_run_dir(flight, unit)
        if not force and is_cache_hit(unit, wingbeat_dir):   # 13.4
            prior = read_status(wingbeat_dir)
            plan.append((unit, "reuse", {
                "finished_at": prior.get("finished_at"),
                "record":      prior.get("record", {})}))
        else:
            # `force` is why it runs when a hit was available; the
            # plan says so rather than leaving the reader to guess.
            plan.append((unit, "run",
                {"reason": "forced" if force else "no usable result"}))
    return plan


function print_reuse_plan(plan, per_unit=False):
    # The counts are the decision and always print; the per-unit
    # lines are the evidence for them and are held back unless
    # asked for (DESIGN 5.7 / 6.2.5).  `per_unit` is true when the
    # module-level verbosity is on OR the caller is a preview,
    # whose entire purpose is those lines -- so the flag is passed
    # in rather than read here, keeping this a pure printer.
    #
    # On a reuse the line carries the facts a judgment would want:
    # when the result finished, and the build recorded behind it.
    if per_unit:
        for (unit, action, detail) in plan:
            print(unit.id, unit.calc, action, detail)
    print(count_of(plan, "reuse"), "to reuse,",
          count_of(plan, "run"), "to run")


function send_off(flight, units, executor, force):
    # Phase 1, made callable on its own: launch a chosen set of
    # units and return one future per unit WITHOUT waiting on any
    # of them (DESIGN 6.2.3).  `units` is what to launch now --
    # flight.units for a one-shot fan-out, one round of newly
    # decided rungs for the climb (4e.5).  The flight's WHOLE unit
    # list is (re)serialized, so flight.toml records every rung
    # asked for even as the climb's list grows (13.1).
    validate_flight(flight)            # 13.3
    makedirs(flight.root, exist_ok=True)
    serialize_flight(flight)           # 13.1
    # Announce before spending: counts always, the per-unit lines
    # only when the module-level verbosity switch is on (DESIGN
    # 5.7).  The plan is recomputed inside dispatch_unit below
    # rather than threaded through, because a hit-test is a few
    # local file reads and a stale plan would be worse than a
    # repeated one.
    print_reuse_plan(reuse_plan(flight, units, force),
                     per_unit = is_verbose())
    outstanding = []                   # list of (unit, future)
    for unit in units:
        outstanding.append(
            (unit, dispatch_unit(flight, unit, executor, force)))
    return outstanding
```

`is_verbose()` reads a module-level switch the driver's own
reporting helper owns, set once by whichever client is driving
(the producer's `main` sets it alongside its own, DESIGN 5.7).  It
is deliberately NOT a `send_off` argument: verbosity describes how
the process talks to its user, not how a flight dispatches, and
threading it would put a reporting concern into the signature of
every function between `main` and the printer -- the same
reasoning, and the same conclusion, as the producer's side of the
boundary.  Reporting is not domain knowledge, so owning a
verbosity switch costs kaleidoscope none of its ignorance about
k-points (Principle 9).

```
function collect(flight, unit, fut):
    # Phase 2 for ONE unit: resolve its future, write its terminal
    # status, build its report entry, and fire the stream hook
    # (DESIGN 6.2.3, 6.2.6).  Per-future exception capture so one
    # failure never propagates (Principle 10).  Public, so a
    # control loop can collect rungs one at a time.
    wingbeat_dir = unit_run_dir(flight, unit)
    try:
        fut.result()       # re-raises any worker exception
    except ParslTaskLost:
        # No WingbeatOutcome ever came back: cluster-side loss.
        write_status(wingbeat_dir, id=unit.id, calc=unit.calc,
            status="lost", finished_at=now(),
            message="cluster-side loss")
    except Exception as e:
        # App raised on the worker but the failure returned:
        # status.toml may already say running; force failed.
        write_status(wingbeat_dir, id=unit.id, calc=unit.calc,
            status="failed", finished_at=now(),
            message=str(e))
    # In every case status.toml is now terminal; build the report
    # entry from it (single source of truth), then stream it in
    # LANDING order (not unit order) so a control-loop consumer
    # sees each rung the moment it is collected.
    entry = report_entry_from_status(unit, wingbeat_dir)
    if flight.on_outcome is not None:
        flight.on_outcome(entry)
    return entry
```

```
function collect_next(flight, outstanding):
    # Wait for WHICHEVER outstanding rung lands first, collect it,
    # and return it with the shrunken outstanding list (DESIGN
    # 6.2.3).  Domain-ignorant: it polls done() and knows only
    # futures, never k-points, so the climb's decision of what to
    # send next stays in the producer (Principle 12).  A cache hit
    # is an already-done future, so it is returned first, no wait.
    loop forever:
        for (index, (unit, fut)) in enumerate(outstanding):
            if fut.done():
                entry = collect(flight, unit, fut)
                remaining = outstanding with index dropped
                return (unit, entry, remaining)
        sleep(POLL_INTERVAL)   # nothing ready yet; brief pause
```

```
function dispatch(flight, executor=None, force=False,
                  preview=False):
    # The one-shot convenience form: send every unit off, then
    # collect them all in unit order (DESIGN 6.2.3).  Behaviour is
    # identical to the pre-split driver, so every existing caller
    # is unchanged; the climb (4e.5) uses send_off + collect_next
    # directly instead.  We tear the executor down at the end only
    # if we built it here (a caller-supplied executor is the
    # caller's to close).
    #
    # `preview` prints the plan and stops -- no executor is even
    # built, so the decision to spend can be made BEFORE a flight
    # starts rather than watched going past during it (DESIGN
    # 6.2.5).  Here the per-unit lines print whether or not
    # verbosity is on: a preview that showed only the counts would
    # answer nothing the caller could not already guess.  It
    # returns an empty report, not a partial one: no unit ran, so
    # there is nothing to report on.
    if preview:
        print_reuse_plan(reuse_plan(flight, flight.units, force),
                         per_unit = True)
        return FlightReport(entries = [])
    owns_executor = (executor is None)
    if executor is None:
        executor = make_executor(flight.parsl_config)  # below
    try:
        outstanding = send_off(flight, flight.units,
                               executor, force)
        # Collect in send order, which is unit order; on_outcome
        # fires inside collect, once per entry.
        entries = [collect(flight, unit, fut)
                   for (unit, fut) in outstanding]
    finally:
        if owns_executor:
            executor.close()      # tear down Parsl if we built it
    return FlightReport(entries = entries)
```

```
function make_executor(parsl_config):
    # The one seam that turns a resolved config into a live
    # executor, shared so the choice lives in exactly one place.
    # A cluster config -> a ParslExecutor (constructing it loads
    # one coordinator and its SLURM worker pool, 13.7); the local
    # opt-out (config None) -> an in-process LocalExecutor.  A
    # client that dispatches many flights under one config calls
    # this ONCE and pins the result to every dispatch (the
    # producer's climb, 11.4), so the whole run shares one pool; a
    # one-shot flight lets dispatch build and close its own.
    if parsl_config is not None:
        return ParslExecutor(parsl_config)      # 13.7
    return LocalExecutor()
```

```
function dispatch_unit(flight, unit, executor, force):
    wingbeat_dir = unit_run_dir(flight, unit)
    # The re-run switch bypasses the run-reuse cache: with force
    # set, even a still-valid `done` unit re-launches.  force is
    # a driver-level instruction about the cache (DESIGN 6.2.5),
    # independent of which executor runs the unit -- so it rides
    # here, not on a worker.
    if not force and is_cache_hit(unit, wingbeat_dir):   # 13.4
        # Hit: no task submitted.  Return an already-done future
        # (done() true, result() a no-op) so the hit sits in the
        # outstanding set exactly like a miss; collect reads the
        # entry back from the existing status.toml.
        return completed_future()
    # Miss: prepare the dir, snapshot the key, mark queued,
    # and hand the unit to the executor.  The unit's `record`
    # is stamped here, once, alongside the key snapshot: the key
    # holds what is COMPARED, the record holds what is only ever
    # read by a person (DESIGN 6.2.4/6.2.5).  It is written on the
    # miss only, so a later hit leaves it describing the run that
    # produced the stored result, not the flight that reused it.
    makedirs(wingbeat_dir, exist_ok=True)
    write_cache_key(wingbeat_dir, unit)          # 13.4
    write_status(wingbeat_dir, id=unit.id, calc=unit.calc,
        status="queued",
        wingbeat=(unit.wingbeat or flight.default_wingbeat),
        submitted_at=now(),
        record=unit.record)
    return executor.submit_unit(
        unit, wingbeat_dir, flight.default_wingbeat)
```

```
function execute_wingbeat_task(unit, wingbeat_dir, default_wingbeat):
    # The unit of work both executors run: LocalExecutor calls it
    # in process, ParslExecutor wraps it as a Parsl task.  It
    # takes default_wingbeat (not the whole flight, which need not
    # travel to a worker).  Returns the WingbeatOutcome; raising
    # here surfaces to collect as a worker-side failure.
    write_status(wingbeat_dir, id=unit.id, calc=unit.calc,
        status="running", started_at=now())
    wingbeat = resolve_wingbeat(unit, default_wingbeat)  # ->Wingbeat
    outcome = wingbeat.run(unit, wingbeat_dir)         # 13.2
    write_status(wingbeat_dir, id=unit.id, calc=unit.calc,
        status=("done" if outcome.ok else "failed"),
        detail=outcome.detail,
        finished_at=now(),
        runtime_seconds=outcome.runtime_seconds,
        message=outcome.message)
    return outcome
```

### 13.6 Report and client-side harvest (DESIGN 6.2.6)

```
dataclass ReportEntry:
    id, calc, status, detail, wingbeat_dir,
    runtime_seconds, message

dataclass FlightReport:
    entries : list      # list[ReportEntry]

function report_entry_from_status(unit, wingbeat_dir):
    st = read_status(wingbeat_dir)
    return ReportEntry(
        id=unit.id, calc=unit.calc,
        status=st["status"], detail=st.get("detail"),
        wingbeat_dir=wingbeat_dir,
        runtime_seconds=st.get("runtime_seconds"),
        message=st.get("message"))
```

Harvest is *not* kaleidoscope's job (Principle 9).  The
handoff is the run directory: the wingbeat persisted its
native `result.toml` there (13.2), so the client walks the
report and reads what it needs from the dirs it deems
acceptable.  The C48 producer's harvest, which lives in
`build_initial_potentials.py` and is shown here only to
fix the contract:

```
# CLIENT side (build_initial_potentials.py), NOT
# kaleidoscope.  This is the precise C48.3 shape.
function harvest_converged_potentials(report, manifest):
    for entry in report.entries:
        # Keep only scientifically acceptable units; the
        # client owns this judgment, not kaleidoscope.
        if entry.detail != "converged":
            continue                 # skip, recorded in report
        result = read_toml(
            join(entry.wingbeat_dir, "result.toml"))   # §12.1
        scfV_path = result["outputs"]["scfV"]
        # The converged scfV output lists every potential
        # type (NUM_TYPES + per-type blocks under a
        # TOTAL__OR__SPIN_UP channel; the producer runs
        # non-spin, so that channel is the total potential).
        # Select the harvested site's type block and take each
        # term's coefficient and alpha (columns 1-2) together
        # (DESIGN 5.7 / ARCHITECTURE 9.7).
        coeffs, alphas = read_scfV_type_block(
            scfV_path, site_type(entry))
        store_potential_entry(entry, coeffs, alphas)
```

This closes the loop with §12 and with the producer: the
flight runs and tracks the batch and owns the cache;
the client declares the units and the key, then harvests
converged potentials from the run directories the report
points at.

### 13.7 Cluster dispatch configuration (DESIGN 6.2.11)

The generator turns (site facts + per-run choices) into the
`flight.parsl_config` the driver (13.5) loads.  It lives in
kaleidoscope so every client shares one copy.  The
command-line default is `slurm-pooled`, taken from the site's
`default_topology` (not hardcoded on the flag): the producer
and the seed are *meant* to reach the scheduler, so on a
cluster they do so with no flags, and a run with no settings
file present is a configuration error rather than a quiet
local fall-back (DESIGN 6.2.11, decision 2).  `local` is the deliberate
opt-out -- it builds no `Config` (returns None), needs no
settings file, and the driver (13.5) runs it under a
`LocalExecutor` in process.  The test suite, a laptop, and the
materialize pre-flight all request `local` explicitly, so they
neither read a settings file nor touch the scheduler.

Three layers feed it: the per-site resource-control file
(stable facts), the per-run choices (CLI flags), and -- the
deferred third -- the per-unit size, which is one uniform
slice for now (DESIGN 6.2.11, decision 3).  TODO C81 layers
predictive per-unit sizing on top later without disturbing
this code.

```
# Layer 1: the tiered per-site resource-control file.  A tiny
# required core; every other key is optional and falls back to
# the default shown.  Same *rc.py convention as the other
# resource-control files -- PURE DATA: a module returning one
# dict, nothing more.  The two required fields ship as None (a
# REQUIRED comment marks them); load_site_config refuses an
# unfilled one.  cluster_probe.py generates a filled starter.
#
# Most of the file is cluster FACT (queues, account, per-node
# capacity, environment bring-up), true no matter what runs.  The
# rest sizes a job, and there are exactly three job CLASSES: a
# worker (one calculation, the per-worker keys), an orchestrator
# (a driver that prepares units and fans them out), and an md job
# (an external MD program run as many MPI ranks on one node) --
# the last two as the grouped blocks at the end.  The file grows
# by job class, never by builder (ARCHITECTURE 9.4).
function clusterrc.parameters_and_defaults():
    return {
        # --- Required core (enough to dispatch at all) ---
        "partitions"        : None,      # REQUIRED list; [0]=default
        "worker_init"       : None,      # REQUIRED shell bring-up
        "account"           : None,      # set where required
        # --- Performance tuning (optional) ---
        "cores_per_node"    : None,      # None -> 1 worker/node
        "workers_per_node"  : None,      # None -> derive below
        "cores_per_worker"  : 1,         # serial imago today
        "nodes"             : 1,
        "walltime"          : "01:00:00",
        "default_topology"  : "slurm-pooled",
        "max_blocks"        : 1,         # pooled growth cap
        # Two DISTINCT memory concepts, deliberately not merged.
        #   memory_per_node is a node's physical capacity, in
        #   MEGABYTES -- a ceiling held for future packing and
        #   estimation checks, never spent as a request.
        #   memory_per_worker is what ONE calculation needs, in
        #   GIGABYTES -- the figure scheduler_options turns into a
        #   request.  None on either means "let the scheduler
        #   apply its own default" (DESIGN 6.2.11).
        "memory_per_node"   : None,
        "memory_per_worker" : 10,
        # --- Advanced / forward-looking (power users) ---
        "launcher"          : "single",  # MPI/GPU seam later
        "ranks_per_worker"  : 1,         # OpenMP/MPI balance:
        "threads_per_rank"  : 1,         #   product = cores held
        "binding"           : None,      # core/socket/NUMA pin
        "omp_places"        : None,      # finer thread placement
        "omp_proc_bind"     : None,
        "gpus_per_node"     : 0,
        "queue_overrides"   : {},        # per-queue key tweaks
        "profiles"          : {},        # named site profiles
        "extra_scheduler_options" : [],  # raw passthrough
        # --- The orchestrator job class (DESIGN 6.2.11) ---
        # The resources the DRIVER process asks for when it is
        #   wrapped in its own batch job.  ONE default shape shared
        #   by every orchestrator, overridden key by key per run
        #   via --orchestrator-{cores,memory,walltime}
        #   (resolve_orchestrator, below).  Under a fan-out
        #   dispatch the driver only prepares units and submits
        #   them, so it is modest; under --dispatch local it runs
        #   the SCFs in process, and that run raises the shape on
        #   the command line rather than editing this file.
        "orchestrator" : {
            "cores"    : 2,
            "memory"   : "8G",
            "walltime" : "24:00:00",     # outlast the flight
        },
        # --- The md job class (DESIGN 6.2.11) ---
        # An external molecular-dynamics program (LAMMPS today)
        #   run under MPI: many ranks filling ONE node, so it is
        #   sized by rank count rather than by the cores-per-task
        #   an orchestrator asks for.  `ranks` None falls back to
        #   cores_per_node, and where the site recorded neither,
        #   to a single rank -- visibly wrong to whoever opens the
        #   generated file, which is the intent.
        #
        # `init` is this class's OWN bring-up: worker_init starts
        #   imago, which an external program neither needs nor is
        #   served by, and holding the two apart is what keeps
        #   that program's install location out of src/ entirely.
        #   It ships blank and REQUIRED, as worker_init does: the
        #   sizing keys above merely size a job, whereas an absent
        #   bring-up leaves nothing on the path to run at all.
        #   But it is deliberately NOT checked by _require_core --
        #   that check guards every dispatch, and a site that
        #   flies calculations and never condenses must not be
        #   refused a flight over a setting no flight reads.  The
        #   md generator enforces it instead (build_md_sbatch).
        "md" : {
            "ranks"    : None,           # None -> cores_per_node
            "walltime" : "01:00:00",
            "memory"   : None,           # None -> no --mem
            "init"     : None,           # REQUIRED shell bring-up
        },
    }
```

```
function merge_settings(base, overlay):
    # The ONE merge every overlay uses -- profile, queue, and the
    # per-run flags alike (DESIGN 6.2.11, decision 1).  Per key, and
    # one level down: when a setting is itself a BLOCK of settings
    # (`orchestrator` and `md`), the overlay names only the keys it
    # means to change and the rest keep the value the layer beneath
    # gave them.
    #
    # Replacing the whole block instead would silently discard facts
    # the curator never mentioned -- "the driver needs 2G on the
    # debug queue" would also drop its cores and walltime, and they
    # would reappear as plausible-looking fallbacks rather than as
    # an error.  A block holds plain values, never further blocks,
    # so one level is the whole of the descent.  The md block's
    # `init` keeps that intact: a list of shell lines is a value,
    # exactly as top-level worker_init is, not a nested block.
    result = copy(base)
    for key, value in overlay.items():
        if is_block(base.get(key)) and is_block(value):
            result[key] = merge(base[key], value)   # one level
        else:
            result[key] = value
    return result


function load_site_config(profile=None, partition=None):
    # Reading the settings file and overlaying it are ONE
    #   operation (DESIGN 6.2.11, decision 1).  The loader takes
    #   the queue so it can apply the queue overlay itself; there
    #   is deliberately no way to obtain un-overlaid settings, so
    #   no reader can forget to overlay and end up with the
    #   cluster-wide walltime where a queue's cap belongs.
    #
    # Resolve clusterrc.py by precedence: the working directory
    #   first (a per-run override), then $IMAGO_RC (the global
    #   default).  cluster_probe.py uses the same order.
    site = clusterrc.parameters_and_defaults()
    # Overlay 1: a named profile (advanced tier), so a user with
    # several clusters selects one by name.
    if profile is not None:
        site = merge_settings(site, site["profiles"][profile])
    # The required core must be present; everything else has a
    # default.  A gap here is a config error raised up front,
    # never a crash mid-flight (the strict-contract discipline
    # the producer already follows, DESIGN 6.3.1).  It is checked
    # BEFORE the queue overlay because picking the default queue
    # reads `partitions`.
    require_core(site)               # partitions, worker_init
    # Overlay 2: the selected queue.  The queue is a per-run
    # choice, defaulting to the first entry of the (now
    # profile-overlaid) partitions list.
    queue = partition or site["partitions"][0]
    site = apply_queue_overrides(site, queue)
    # An override may set worker_init, so re-check the core it
    # could have emptied.  (It may not set partitions.)
    require_core(site)
    return site
```

The queue overlay itself (DESIGN 6.2.11, decision 1).  It is
called only from the loader above; it is kept a named step so
its two guards can be read -- and tested -- on their own.

```
function apply_queue_overrides(site, partition):
    # Overlay ONLY the selected queue's settings.  A file may
    # carry overrides for every queue on the cluster; the ones
    # this run does not use are simply not applied.
    override = site["queue_overrides"].get(partition)
    if override is None:
        return site

    for key in override:
        # A key naming no known setting is a TYPO, and a silently
        # ignored typo in a resource request is exactly what this
        # file exists to prevent.  Refuse it up front.
        if key not in site:
            raise ConfigError(
                "queue override for " + partition +
                " names unknown setting " + key)
        # partitions / profiles choose WHICH overlay applies, so an
        # overlay that rewrote them would refer to itself.
        if key in ("partitions", "profiles"):
            raise ConfigError(
                "queue override may not set " + key)
        # The guard descends one level, exactly as far as the merge
        # below does, so the merge cannot reach a place the guard
        # cannot see.  A typo INSIDE a block is the quieter fault:
        # it leaves the real key standing at its old value beside
        # the stray one, so the run uses the number the curator
        # meant to change and nothing says otherwise ("rank" for
        # "ranks" runs the job at the site's width).
        if is_block(site[key]) and is_block(override[key]):
            for inner in override[key]:
                if inner not in site[key]:
                    raise ConfigError(
                        "queue override for " + partition +
                        " names unknown setting " + key +
                        "." + inner)

    # Per key, one level down: a queue that names only the driver's
    # memory keeps the site's driver cores and walltime.
    return merge_settings(site, override)
```

```
# Layer 2: per-run choices.  Each CLI flag defaults from the
# site file -- the dispatch shape from default_topology
# (slurm-pooled) -- so a fully configured site needs no flags
# at all.  The flag surface is --dispatch {local, slurm-pooled,
# slurm-per-job}, --partition, --nodes, --walltime.  This runs
# only for a cluster shape; the local opt-out skips it (13.7
# run_flight short-circuits before load_site_config).  These four
# size the WORKER job class; the orchestrator's own shape is
# resolved separately below.
function resolve_choices(site, cli):
    return {
        "dispatch"  : cli.dispatch  or site["default_topology"],
        "partition" : cli.partition or site["partitions"][0],
        "nodes"     : cli.nodes     or site["nodes"],
        "walltime"  : cli.walltime  or site["walltime"],
    }
```

The orchestrator's shape resolves on its own, because it sizes a
different job class (a driver, not a calculation) and the two
must not be conflated (DESIGN 6.2.11).  This is what makes the
single site-default block *bounded*: a second orchestrator with
different needs overrides the shape for its run rather than
earning a second block in the settings file (ARCHITECTURE 9.4).

```
# Layer 2b: the driver's own per-run shape.  The flag surface is
# --orchestrator-cores, --orchestrator-memory,
# --orchestrator-walltime (DESIGN 6.2.11, decision 2).  An
# explicit flag wins over the site's orchestrator block, KEY BY
# KEY: overriding the memory must leave the site's cores and
# walltime standing, which a whole-block replacement would
# silently discard.  A key nobody sets stays absent, and
# build_orchestrator_sbatch decides what an absent key means --
# cores and memory simply go unrequested, while walltime falls
# back once more, to the run's resolved --walltime, so a driver
# job always carries a time limit.
#
# NOTE the worker flags do NOT reach here: --walltime and --nodes
# size the WORKER class.  A curator shortening --walltime is
# speaking about the calculations, not about the process that
# submits them.
function resolve_orchestrator(site, cli):
    shape = copy(site["orchestrator"])      # may be {} at a site
    for key in ("cores", "memory", "walltime"):
        flag = getattr(cli, "orchestrator_" + key, None)
        if flag is not None:
            shape[key] = flag
    return shape
```

```
function build_dispatch_config(site, choices):
    # Local: no Parsl Config at all.  The driver runs the flight
    # in process (13.5), one unit at a time, exactly as before.
    if choices["dispatch"] == "local":
        return None
    if choices["dispatch"] == "slurm-pooled":
        return build_pooled_config(site, choices)
    if choices["dispatch"] == "slurm-per-job":
        return build_per_job_config(site, choices)
    raise ConfigError("unknown dispatch " + choices["dispatch"])
```

Both cluster shapes are the *same* SlurmProvider wiring with
different block geometry, so one helper turns site facts plus
per-run choices into the provider and the two builders only
differ in how blocks map to units.

```
function slurm_provider(site, choices, nodes_per_block,
                        init_blocks, min_blocks, max_blocks,
                        workers_per_block = 1):
    # Guard the deferred parallel seam first: refuse any parallel
    # knob set away from its serial default (see _require_serial_only).
    _require_serial_only(site)
    # worker_init is the site's bring-up script, so a worker can
    # find imago; account/partition/walltime come from the resolved
    # choices; the memory, core, and GPU knobs ride along as raw
    # scheduler directives (scheduler_options).  workers_per_block is
    # how many calculations share a node, so the per-node memory and
    # core requests both scale with it.
    #
    # exclusive is stated, not defaulted.  The provider's own default
    # claims the whole node, which would undo the slice the request
    # just asked for (DESIGN 6.2.11): a one-core block would hold
    # every core on the node, and sibling blocks would each queue for
    # a node of their own rather than sharing one.  scheduler_options
    # names the cores that this then leaves us, so the two belong
    # together -- neither is correct without the other.
    return SlurmProvider(
        partition         = choices["partition"],
        account           = site["account"],
        walltime          = choices["walltime"],
        nodes_per_block   = nodes_per_block,
        init_blocks       = init_blocks,
        min_blocks        = min_blocks,
        max_blocks        = max_blocks,
        worker_init       = join_lines(site["worker_init"]),
        launcher          = make_launcher(site),
        exclusive         = false,
        scheduler_options = scheduler_options(
                                site, workers_per_block))
```

```
function build_pooled_config(site, choices):
    # One (optionally auto-scaled) allocation; many units stream
    # through its workers.  Size the block by the per-run nodes
    # and the site's per-node worker packing.  max_blocks lets
    # the pool grow when work backs up.  The SAME packed-worker
    # count caps the executor and scales the memory request.
    packed_workers = workers_per_node(site)
    provider = slurm_provider(site, choices,
        nodes_per_block = choices["nodes"],
        init_blocks = 1, min_blocks = 1,
        max_blocks  = site["max_blocks"],
        workers_per_block = packed_workers)
    executor = HighThroughputExecutor(
        label                = "imago-pooled",
        provider             = provider,
        cores_per_worker     = site["cores_per_worker"],
        max_workers_per_node = packed_workers)
    return Config(executors = [executor])
```

```
function build_per_job_config(site, choices):
    # One scheduler submission per unit: each unit maps to its
    # own one-node, one-worker block, so calculations queue and
    # run independently.  max_blocks bounds how many run at once.
    provider = slurm_provider(site, choices,
        nodes_per_block = 1,
        init_blocks = 0, min_blocks = 0,
        max_blocks  = site["max_blocks"],
        workers_per_block = 1)      # one calc/node -> one worker's mem
    executor = HighThroughputExecutor(
        label                = "imago-per-job",
        provider             = provider,
        max_workers_per_node = 1)   # exactly one unit per block
    return Config(executors = [executor])
```

```
function workers_per_node(site):
    # Explicit override wins; else derive from the node's cores
    # and the per-worker core count; else fall back to one worker
    # per node (the no-cores_per_node default).
    if site["workers_per_node"] is not None:
        return site["workers_per_node"]
    if site["cores_per_node"] is not None:
        return max(1, site["cores_per_node"] //
                      site["cores_per_worker"])
    return 1
```

```
function _require_serial_only(site):
    # The deferred parallel seam.  Today imago is serial: one
    # calculation on one core through the single-node launcher.  The
    # knobs that describe a parallel calculation -- a non-serial
    # launcher, the ranks_per_worker x threads_per_rank split, and the
    # binding / omp_* placement -- cannot be realized yet, so any of
    # them set away from its serial default is a clear error rather
    # than a silently dropped setting.  A real MPI launcher (TODO C100
    # / C81) replaces this guard and consumes those knobs at launch.
    serial_defaults = {launcher: "single", ranks_per_worker: 1,
        threads_per_rank: 1, binding: None, omp_places: None,
        omp_proc_bind: None}
    for knob, default in serial_defaults:
        if site[knob] != default:
            raise NotImplementedError(knob + " is the deferred "
                "parallel-imago seam; leave it at its serial default")
```

```
function make_launcher(site):
    # Serial today: one process per worker.  The MPI launcher that
    # would honour the ranks_per_worker x threads_per_rank split and
    # the binding / omp_* placement is the deferred seam (already
    # refused by _require_serial_only), so a non-serial launcher
    # raises rather than returning a launcher that does not exist yet.
    if site["launcher"] == "single":
        return SingleNodeLauncher()
    raise NotImplementedError("MPI launcher is the deferred seam")
```

```
function scheduler_options(site, workers_per_block = 1):
    # Assemble the raw scheduler directives the site facts imply --
    # the memory guard and the GPU request, which are allocation-time
    # requests and so belong in #SBATCH -- then append
    # extra_scheduler_options verbatim so a power user is never
    # blocked by the schema.  CPU/NUMA binding is NOT here: pinning is
    # a launch-time concern (srun --cpu-bind, OMP_PLACES), applied by
    # the launcher in the deferred parallel path, not a batch
    # directive.  Returns the directives joined into one string.
    #
    # The memory guard is DERIVED, not copied.  memory_per_worker is
    # what one calculation needs, while the scheduler's --mem is a
    # per-NODE figure, and a node runs workers_per_block calculations
    # at once (one under the per-job shape, the node's packed worker
    # count under the pooled shape).  memory_per_node is never spent
    # here: it is capacity, not a request (DESIGN 6.2.11).
    #
    # The core request is derived the SAME way and for the same
    # reason: the block asks for its own workers' slices and no more
    # (DESIGN 6.2.11, "a block asks for its slice, not for the node").
    # It must be stated -- a block that omits it takes the scheduler's
    # one-core default, which silently starves a packed pool.  One
    # task per node holds the whole worker pool (the dispatch driver
    # submits tasks_per_node = 1), so the node's cores are the task's
    # cores and --cpus-per-task carries the count.
    lines = []
    if site["memory_per_worker"]:
        node_memory_gb = (site["memory_per_worker"]
                          * workers_per_block)
        lines.append("#SBATCH --mem=" + node_memory_gb + "G")
    node_cores = site["cores_per_worker"] * workers_per_block
    lines.append("#SBATCH --cpus-per-task=" + node_cores)
    if site["gpus_per_node"] > 0: lines.append(gres_line(site))
    lines.extend(site["extra_scheduler_options"])
    return join_lines(lines)
```

The producer change-over (DESIGN 6.2.11; TODO C100) needs no
client-specific worker-builder.  Because the generator lives
in kaleidoscope and the driver auto-selects the executor
(13.5), every client -- the producer, the validation harness,
future flights -- uses the SAME two steps: attach the
generated config, then call dispatch.  The producer's old
`curation_executor` (which always returned a `LocalExecutor`
and carried the re-run switch) goes away; the re-run switch
becomes a `dispatch` argument, where the run-reuse cache it
governs already lives.

Those two steps are themselves shared, in kaleidoscope, so no
client copies them.  `resolve_dispatch` turns a dispatch choice
into `(parsl_config, choices)`, and `write_resolved_dispatch`
records the resolved choices beside the run; a client wires
them around its own dispatch call:

```
# SHARED in kaleidoscope -- the run_flight steps every client uses.
function resolve_dispatch(dispatch, partition, nodes, walltime,
                          profile):
    if dispatch == "local":
        # The deliberate opt-out: no settings file is read and no
        #   Config is built, so a laptop or the test suite runs in
        #   process without a clusterrc present at all.
        return (None, None)
    # A cluster shape: the settings file IS required here, and a gap
    #   in its required core is a config error raised up front
    #   (load_site_config), never a quiet local fall-back.  Passing
    #   the queue is what makes the loader overlay it; the remaining
    #   per-run choices then default from the overlaid site (DESIGN
    #   6.2.11, decision 1: defaults -> profile -> queue -> flags).
    site    = load_site_config(profile, partition)
    choices = resolve_choices(site, {dispatch, partition, nodes,
                                     walltime})
    return (build_dispatch_config(site, choices), choices)

function write_resolved_dispatch(run_dir, choices, profile):
    # A small human-readable record for reproducibility: the dispatch
    #   shape, queue, nodes, and walltime actually used, plus the
    #   profile when one was selected.  The stable site facts are NOT
    #   duplicated -- they live in clusterrc.py, and the profile name
    #   pins which overlay fed this run.
    write_lines(run_dir + "/resolved_dispatch.toml",
                ["profile" if profile, dispatch, partition,
                 nodes, walltime])
```

```
# CLIENT side -- the SAME shape for every flight; the client adds
#   only its own build (before) and harvest (after) around this.
function run_flight(flight, cli):
    flight.parsl_config, choices = resolve_dispatch(cli.dispatch,
        cli.partition, cli.nodes, cli.walltime, cli.profile)
    if cli.save_config and choices is not None:
        write_resolved_dispatch(flight.root, choices, cli.profile)
    # No explicit executor: the driver picks Local or Parsl from
    #   flight.parsl_config.  The re-run switch rides as a driver
    #   argument, not bundled into a worker (DESIGN 6.2.5).
    return dispatch(flight, force=cli.force)
```

The run-reuse cache (13.4) is untouched: a worker executes
each unit in its own run directory on the shared filesystem
exactly as the in-process local path does, so a cluster run
and a local run share one cache.

**The driver's own batch job.**  Everything above places each
*unit*.  It leaves open where the *driver* runs -- the
orchestrator process that reads the manifest, prepares every
unit (13.2 / 11.4), decides cache hits from local files, and
submits and awaits the rest.  Because the driver now does real
per-unit work before any SCF (a makeinput build, plus a fast
`imago -loen` when a solid's species assignment needs one, once
per unit *including cache hits*), at scale it would occupy a
login node's terminal for the whole flight.  So the driver may
be wrapped in a scheduler job of its own, sized from the site's
`orchestrator` block -- a job class distinct from the per-worker
sizing, because one process that fans work out and one process
that runs a calculation are different shapes and conflating them
would missize both (DESIGN 6.2.11; ARCHITECTURE 9.4).

The generator lives in kaleidoscope beside the `Config`
builders, so every future orchestrator renders its job the same
way.  Note the `extra_scheduler_options` passthrough: those
entries are already complete `#SBATCH` lines, exactly as
`scheduler_options` forwards them to Parsl, so they are copied
verbatim rather than given a second directive marker.

```
function build_orchestrator_sbatch(site, choices, command,
                                   orchestrator = None):
    # The driver is ONE process, so the header asks for one node
    # with the orchestrator shape.  `orchestrator` is the shape
    # resolve_orchestrator merged from the site block and this
    # run's flags; when a caller passes none, the site block
    # stands alone.  A missing walltime falls back to the run's
    # resolved walltime, so the driver's job always carries a
    # time limit.  `command` is the already-quoted command line
    # the batch job runs.  Test for ABSENCE, not falsiness: a caller
    # that deliberately passes an empty shape means "request nothing
    # but the fallbacks," and must not silently re-inherit the site
    # block.
    orch     = (orchestrator if orchestrator is not None
                else site.get("orchestrator", {}))
    cores    = orch.get("cores", 1)
    memory   = orch.get("memory")                 # None -> no --mem
    walltime = orch.get("walltime") or choices["walltime"]

    # A login shell, so worker_init below runs in a shell whose
    # profile has been read.  Where that bring-up uses `module`,
    # a plain shell would work only when the submitting shell had
    # already been set up, and fail from cron, a workflow driver,
    # or --export=NONE (DESIGN 6.2.11).
    lines = ["#!/bin/bash -l",
             "#SBATCH --job-name=imago-orchestrator"]
    if site["account"]:                           # some sites need none
        lines.append("#SBATCH --account=" + site["account"])
    lines.append("#SBATCH --partition=" + choices["partition"])
    lines.append("#SBATCH --nodes=1")
    lines.append("#SBATCH --cpus-per-task=" + cores)
    if memory:
        lines.append("#SBATCH --mem=" + memory)
    lines.append("#SBATCH --time=" + walltime)
    lines.extend(site["extra_scheduler_options"])   # already #SBATCH

    # The bring-up runs first so the batch job can find imago, then
    # the producer command itself.
    lines.append("")
    lines.extend(site["worker_init"])
    lines.append("")
    lines.append(command)
    return join_lines(lines)
```

`join_lines` finishes with a trailing newline, in this generator
and in the md one below, so the file ends the way a text file is
expected to end and the last command is a complete line.

**The md generator.**  The second submission-file generator sits
beside the first, so a reader who finds one finds the other and
the two stay alike where they can -- both open with a login
shell, for the reason DESIGN 6.2.11 gives.  It differs in four
ways, each following from the job being many MPI ranks of an
external program rather than one driver process: it asks for
`--ntasks` instead of `--cpus-per-task`; it runs the md block's
own bring-up instead of `worker_init`, and refuses to write a
file at all when that bring-up is missing; and it pins one
thread per rank, because ranks sized to fill a node must each
hold a single core.

```
function build_md_sbatch(site, choices, command, md = None):
    # Test for ABSENCE, not falsiness, exactly as the orchestrator
    # generator does: a caller passing an empty shape means
    # "request nothing but the fallbacks" and must not silently
    # re-inherit the site block.  The SIZING keys fall back; the
    # bring-up below does not, so a site with no md block at all
    # is refused rather than handed a job that cannot start.
    shape    = (md if md is not None else site.get("md", {}))
    memory   = shape.get("memory")                # None -> no --mem
    walltime = shape.get("walltime") or choices["walltime"]

    # The bring-up is required HERE rather than in the loader's
    #   required core (DESIGN 6.2.11).  Without it nothing puts the
    #   MD program on the path and the job cannot start -- but a
    #   flight that never condenses must not be refused over it, so
    #   this generator is the one that insists.  Same emptiness test
    #   the loader uses: None and [] are both unfilled.
    bring_up = shape.get("init")
    if is_empty(bring_up):
        raise ConfigError("the cluster settings file records no md "
                          "bring-up (md.init), so the generated job "
                          "would have no way to find the MD program; "
                          "fill it in clusterrc.py")

    # Ranks come from the node, never from a number written into
    #   the source: a hard-coded count is what once asked for 125
    #   tasks while naming a partition whose nodes hold 48.  Where
    #   the site recorded no core count there is nothing to derive
    #   from, so ask for ONE rank and say so in the file itself.
    #   A one-rank MD job is visibly wrong to whoever opens it,
    #   whereas a guessed count would run, and run wrong, without
    #   ever announcing that the site was never configured.
    ranks   = shape.get("ranks") or site["cores_per_node"]
    unsized = (ranks is None)
    if unsized:
        ranks = 1

    lines = ["#!/bin/bash -l",
             "#SBATCH --job-name=lmp"]
    if site["account"]:                           # some sites need none
        lines.append("#SBATCH --account=" + site["account"])
    lines.append("#SBATCH --partition=" + choices["partition"])
    lines.append("#SBATCH --nodes=1")
    lines.append("#SBATCH --ntasks=" + ranks)
    if memory:
        lines.append("#SBATCH --mem=" + memory)
    lines.append("#SBATCH --time=" + walltime)
    lines.extend(site["extra_scheduler_options"])   # already #SBATCH

    # Directives are done; a comment here cannot swallow one.
    if unsized:
        lines.append("")
        lines.append("# One rank only: this site's settings file")
        lines.append("#   records no cores_per_node, so there was")
        lines.append("#   nothing to size this job from.  Set it,")
        lines.append("#   or the md block's ranks, and rerun.")

    # The bring-up runs first so the batch job can find the MD
    #   program.  The thread pin comes AFTER it, so that a module
    #   setting a thread count of its own cannot overwrite it: the
    #   ranks were sized to fill the node, so each must hold one
    #   core, and a threaded BLAS left alone would start a thread
    #   per core in EVERY rank (DESIGN 6.2.11).
    lines.append("")
    lines.extend(bring_up)
    lines.append("")
    lines.append("export OMP_NUM_THREADS=1")
    lines.append("")
    lines.append(command)
    # Deliberately NOT written, both to stay alike with the
    #   orchestrator generator: --output/--error, leaving the
    #   scheduler's own default naming; and a cd to the submit
    #   directory, which SLURM has already done.
    return join_lines(lines)
```

**`condense.py` calls it for the condensation run.**  The script
reads `condenserc.py` for the settings that are its own business
and the site file for the cluster facts that are not (DESIGN
6.2.11).  It has no dispatch flags of its own, so its per-run
choices are simply the site's defaults -- but it reaches them
through the same resolver a flight uses, handed an empty flag
set, rather than reading `partitions[0]` and `walltime` out of
the site itself.  That is the case `resolve_choices` already
serves for a fully configured site that passes no flags, and
going through it means the day `condense.py` grows a
`--partition` of its own, the resolution needs no redesign.

It also inherits the loader's refusal.  Where the site's
required core is unfilled, `load_site_config` raises and no
submission file is written, rather than one being assembled
from guesses that would fail later at the scheduler with
nothing pointing back at the settings (DESIGN 6.2.11).

**The site is read at settings time, not at writing time**
(DESIGN 6.2.11).  The natural place to *use* the site file is
the last step of `create_lammps_files`, and that is the wrong
place to *read* it, for two reasons.  A refusal there arrives
after the bonds, the angles, the clustering and the whole
LAMMPS input have been computed, throwing that work away and
meeting the user again on the rerun; and by then the script has
entered the `lammps/` directory, so a read that searches the
current directory first (13.7 `load_site_config`) can resolve a
different settings file than the run started under.  So
`ScriptSettings` reads the site beside its `condenserc.py`
read, before any work, and hands the result down.  It reads it
*after* parsing the command line, so `--help` and a misspelled
flag still answer as they always did rather than being met by a
site-configuration error.

The submission file lands beside `lammps.in` and `lammps.dat`,
named `slurm`, exactly where the hand-built template used to be
written, so nothing downstream of it moves.

```
# Settings time, in ScriptSettings: read the site beside the
#   script's own rc file, before any work has been done and
#   while the current directory is still the one the user
#   launched from.
function condense_read_settings(cli):
    settings      = read_rc("condenserc.py")
    settings      = reconcile(settings, parse_command_line())
    # Raises where the required core is unfilled, and that
    #   refusal is inherited deliberately (DESIGN 6.2.11): a
    #   file naming the wrong queue and carrying no bring-up
    #   fails later, on another machine, with nothing pointing
    #   back to here.  Raised HERE so it costs a user nothing
    #   but the rerun of a run that has not started.
    settings.site = load_site_config()
    return settings
```

```
function condense_write_submission(site):
    # Called at the close of create_lammps_files, which has already
    #   entered the lammps/ directory, so the file lands beside the
    #   input it submits.  The site arrives already loaded, from
    #   settings time; nothing is read from disk here.
    #
    # No dispatch flags of its own, so every choice falls through
    #   to the site default -- reached through the shared resolver
    #   rather than by reading partitions[0] out of the site.
    choices = resolve_choices(site, flags_with_none_set())
    # The one line naming the MD program; everything above it in
    #   the file came from the site.
    command = 'mpirun -np "$SLURM_NTASKS" lmp -in lammps.in'
    write_file("slurm", build_md_sbatch(site, choices, command))
```

**Materialize on the login node, then submit.**  The one step
that needs the network is the structure fetch, because compute
nodes may have no internet.  So a `--submit` run materializes
every structure on the login node first, then submits the batch
job; inside that job the prepare step (11.4) consumes the
already-fetched skeletons and touches no network.  The batch job
re-invokes the producer with the same arguments *minus*
`--submit`, so it runs the build rather than submitting again.
`--submit` and `--materialize-only` are therefore mutually
exclusive: each names a different stopping point.

```
# In build_initial_potentials.py (the client), after the
# materialize pre-flight has fetched every structure.
function submit_orchestrator_batch(argv, args, data_root):
    # `argv` is the FLAG vector this run parsed -- the same list
    # main() was handed, with no program name in it.  It is passed
    # in rather than read from the process, because a library caller
    # may drive main(argv) with a vector that has nothing to do with
    # the process's own arguments; reading the process would then
    # submit a batch job running whatever launched us.
    # The driver's OWN job is sized from the same overlaid site the
    # units are, so the queue rides into the loader here too: a
    # debug queue that caps walltime caps the driver's job as well.
    site         = load_site_config(args.profile, args.partition)
    choices      = resolve_choices(site, args)
    orchestrator = resolve_orchestrator(site, args)
    # Re-run THIS producer inside the batch job, dropping --submit
    # so the batch invocation builds instead of resubmitting.  The
    # --orchestrator-* flags may ride along harmlessly: they are
    # read only when submitting, and the inner run does not submit.
    # The script names ITSELF by path, so the command is correct
    # however this process happened to be launched.
    inner   = [item for item in argv if item != "--submit"]
    command = shell_quote_all(
                  [interpreter, this_script_path, *inner])

    script_path = join(data_root, "orchestrator.sbatch")
    write_file(script_path,
               build_orchestrator_sbatch(site, choices, command,
                                         orchestrator))
    output = run(["sbatch", script_path], check = true)
    return last_token(output)      # "Submitted batch job <id>"
```

Which shape the driver's job uses is the same per-run
`--dispatch` choice: `local` inside the orchestrator job at seed
scale (the driver runs the SCFs in process, so the orchestrator
block is sized compute-heavy), `slurm-per-job` or `pooled` later
-- a flag, not a rewrite.  One thing rides with it, deferred: at
seed scale the driver prepares every unit serially inside its
job, which is exactly what keeps a cache hit off the scheduler
(13.4).  When that serial prepare becomes the bottleneck,
prepare-and-hit-test can move onto dispatched worker units, at
the cost that a hit then occupies a cheap worker slot rather
than being decided driver-local.

The discovery tool `cluster_probe.py` (a SEPARATE program, not
part of the pure-data clusterrc.py; DESIGN 6.2.11) reads what
the machine can report and writes a *starter* settings file, so
a newcomer does not hand-assemble the performance and advanced
tiers.  It is best-effort and scheduler-specific: every query
is wrapped so a missing tool or an unparseable line drops that
fact rather than aborting, and the result is a draft the user
reviews -- never an authority.

```
# In cluster_probe.py (a top-level script), NOT in clusterrc.py.
# ONLY the scheduler is queried -- never the login node's own
#   hardware (lscpu/numactl), which would describe the wrong machine.
function probe_site():
    facts = {}
    rows = parse_sinfo_rows(run_query("sinfo ... %R %c %m %G"))
    if rows:
        facts["partitions"] = distinct_partitions(rows)   # queue list
        # Each per-node number: fill it in only when every node AGREES;
        #   if the nodes disagree (a heterogeneous cluster) do NOT
        #   guess -- record the distinct values under *_options and
        #   leave the setting itself unset for the user.
        for setting in [cores, memory, gpus]:
            values = sorted(distinct values across rows)
            if len(values) == 1: facts[setting]       = values[0]
            else:                facts[setting+"_options"] = values
    # Accounts the user may charge (sacctmgr) -- a hint, not the answer.
    facts["accounts"] = user_associations()      # guarded; may be absent
    return facts

function render_starter_clusterrc(facts):
    # Self-contained: _starter_schema() is THIS tool's own copy of the
    #   settings (NOT imported from clusterrc -- a test keeps the two
    #   identical), so the tool writes a full starter without reading
    #   any clusterrc.py.  Each setting gets a one-line plain-language
    #   note.  Overlay the discovered values; leave blank + FILL IN
    #   both the required core (worker_init, partitions-if-none) AND
    #   any per-node number the nodes disagreed on, listing the values
    #   seen in a "nodes vary" note.  No login-node facts are written.
    #
    # The copy carries every job-class BLOCK the settings file has,
    #   orchestrator and md alike, in the same order and with the
    #   same defaults -- the drift test compares whole dictionaries,
    #   so a block added to one file and not the other fails it.
    #   "md.init" names a key one level inside a block, which is as
    #   deep as a block ever goes; it is blanked for the same reason
    #   worker_init is, being site convention no query can report
    #   (where an MD program was installed is policy, not fact).
    settings = _starter_schema()
    for key in ["partitions", "cores_per_node",
                "memory_per_node", "gpus_per_node"]:
        if key in facts: settings[key] = facts[key]
    return render_full_dict(settings, notes = _SETTINGS,
                            options = facts, account_hint = facts,
                            blanks = ["partitions", "worker_init",
                                      "md.init"])
```

### 13.8 Scratch reclamation (DESIGN 6.2.12)

The mechanism half of reclamation: walk a root, decide
what is reclaimable by the CLIENT's policy, and remove
only scratch.

Two kinds of root are recognized, and one call handles
exactly one of them (DESIGN 6.2.12).  A **workspace**
holds `wingbeats/`, and its units prove they are finished
in `status.toml`; the default policy is the conservative
one -- spent once `done` with a `result.toml` -- and a
client passes its own when it needs the working files for
longer.  A **job tree** holds ordinary `imago.py` run
directories with no flight above them; they prove they
are finished by holding no `imagoLock` and ending their
`runtime` log with the completion marker.

A policy therefore receives `(run_dir, target)` -- the run
directory and its already-resolved scratch -- because the
job-tree test needs to look inside the scratch for the
lock, and a client policy is no worse off for being handed
the path it is deciding about.

Everything here is defined by what it refuses.  Refusals
1-3, 5 and 6 are checked per unit, each recording a reason
so the report can explain a skip rather than silently
passing over it; refusals 4 and 7 are properties of the
walks themselves.

There are two ways in, and they share their judgment.
`plan_reclamation` walks a whole root and is what the
standalone tool and the producer's post-harvest sweep call;
`reclaim_one_dir` judges a SINGLE run directory and removes
it, and is what the in-flight prune of layer (b) calls as
each unit lands.  Both reach the per-directory decision
through one function, `plan_one_dir`, so a campaign pruned
in flight and one swept afterwards cannot be governed by
different rules (DESIGN 6.2.12).

```
function default_reclaim_policy(run_dir, target):
    # WORKSPACE contract.  A unit is spent when it finished
    # AND left a result.  `done` alone is not enough: a run
    # that completed but wrote no result is the state a
    # curator most wants to look at, so it is preserved
    # (DESIGN 6.2.12).  `target` is unused here; the status
    # file is authority enough.
    status = read_status(run_dir)
    if status is None or status["status"] != "done":
        return (False, "not done")
    if not file_exists(join(run_dir, "result.toml")):
        return (False, "no result.toml")
    return (True, "")


function hand_run_policy(run_dir, target):
    # JOB-TREE contract (REFUSAL 6).  A hand run writes no
    # status.toml, but it is not silent: the CLI is a thin
    # wrapper over the same callable core, so it leaves the
    # same two traces every unit does.  BOTH are required.
    #
    # The three names below -- the lock file, the log file,
    # and the completion marker -- are imago.py's to define
    # (imago.LOCK_FILE, imago.RUNTIME_FILE,
    # imago.COMPLETION_MARKER), imported rather than
    # re-spelled here so a change to how the engine names or
    # marks a run cannot silently desync this recognizer
    # (DESIGN 6.2.12).
    #
    # The lock lives in the SCRATCH and is taken before any
    # work begins, so its presence means the run owns the
    # directory now or died without releasing it -- and it
    # is what makes a reclamation racing a just-started run
    # refuse rather than mis-time.
    if file_exists(join(target, LOCK_FILE)):
        return (False, "imagoLock present: running or died")

    # The runtime log is opened in APPEND mode, so a
    # directory run four times holds four markers and only
    # the LAST non-blank line describes the current state.
    # Read the tail; never grep for the marker anywhere.
    lines = read_nonblank_lines(join(run_dir, RUNTIME_FILE))
    if lines is empty:
        return (False, "no runtime log")
    if lines[-1] != COMPLETION_MARKER:
        return (False, "runtime log ends mid-run")
    return (True, "")


function scratch_target(run_dir, scratch_root):
    # Resolve the `intermediate` symlink imago.py created,
    # and REFUSE anything that does not land under the
    # scratch root.  A cleanup tool that a symlink can
    # redirect is a hazard: imago.py renames a stale link to
    # `intermediateFIXME`, and a hand-edited workspace can
    # point anywhere at all.
    link = join(run_dir, "intermediate")
    if not is_symlink(link):
        return (None, "no intermediate link")
    target = resolve_real_path(link)
    if not exists(target):
        return (None, "already reclaimed")
    if not is_under(target, scratch_root):
        return (None, "target outside scratch root: " + target)
    return (target, "")


function find_run_dirs(root):
    # A run directory is one carrying a status.toml: the <calc>
    # level is optional (DESIGN 6.2.4), so a unit may sit
    # directly under its id or one or more levels below, and
    # keying on the status file finds either without assuming
    # a depth.
    #
    # REFUSAL 4: never descend through a symlink.  `intermediate`
    # IS a symlink into the scratch area, so a walk that followed
    # links would wander out of the workspace and could plan a
    # removal outside it entirely.
    for (directory, subdirs, files) in walk(join(root,
            "wingbeats")):
        subdirs = [d for d in subdirs
                   if not is_symlink(join(directory, d))]
        if "status.toml" in files:
            yield directory


function detect_root_kind(root):
    # Decide what KIND of thing the user pointed at, by
    # looking rather than by being told (DESIGN 6.2.12).
    # One call handles one kind, so this decision is made
    # once, up front, and fixes the contract for the run.
    if is_dir(join(root, "wingbeats")):
        return "workspace"
    if any run directory is found by find_job_run_dirs(root):
        return "job-tree"
    return None            # neither; the caller refuses


function find_job_run_dirs(root):
    # A hand-run directory is one carrying an `intermediate`
    # symlink -- there is no status.toml to key on, and no
    # fixed depth, since job trees are organized however
    # their author liked.
    #
    # REFUSAL 4: never descend through a symlink (as in
    # find_run_dirs; `intermediate` IS one).
    #
    # REFUSAL 7: never descend into a workspace.  A job tree
    # may hold one far below it, and walking in would judge
    # provable units by the presumption-based contract.  The
    # directory is yielded as a SkippedWorkspace so the
    # bytes left behind stay visible in the report, and is
    # then pruned from the walk.
    for (directory, subdirs, files) in walk(root):
        if is_dir(join(directory, "wingbeats")):
            yield SkippedWorkspace(directory)
            subdirs = []                  # refusal 7
            continue
        subdirs = [d for d in subdirs
                   if not is_symlink(join(directory, d))]
        if is_symlink(join(directory, "intermediate")):
            yield RunDir(directory)


function plan_one_dir(run_dir, label, scratch_root, policy,
        other_targets = [], older_than = None):
    # Judge ONE run directory and return its plan record.  The
    # whole per-directory decision -- resolve, refuse, apply the
    # policy, measure -- lives here so that the whole-tree planner
    # below and the in-flight prune of layer (b) reach it by the
    # same path and cannot drift apart (DESIGN 6.2.12).
    #
    # `other_targets` is the comparison set for REFUSAL 5, supplied
    # by the CALLER because containment is the one refusal a single
    # directory cannot judge on its own: the whole-tree planner
    # passes every run its walk found, the in-flight caller passes
    # the flight's other units.
    #
    # `label` is only how the report names this run -- the path
    # relative to the walk's base for a whole-tree plan, the unit's
    # own directory name for a single prune.  Nothing is decided
    # from it.
    (target, why) = scratch_target(run_dir, scratch_root)
    if target is None:
        return Skipped(run_dir, label, why)

    # REFUSAL 5: an outer tree holding another run's scratch is
    # deferred, never removed.  Taking it would delete the inner
    # run's working files as collateral, which would turn the
    # "never touch an unfinished run" refusal into a formality.
    # Once the inner ones are gone a later pass takes the outer, so
    # nothing is lost -- only deferred.
    nested = [t for t in other_targets
              if t != target and is_under(t, target)]
    if nested is not empty:
        return Skipped(run_dir, label,
            "holds " + count(nested) + " nested run's scratch; "
            "reclaim those first")

    (spent, why) = policy(run_dir, target)
    if not spent:
        return Skipped(run_dir, label, why)

    # Age is the NEWEST mtime anywhere in the tree, never the top
    # directory's: that moves only when entries are added or
    # removed, so a job that has spent a week writing into an
    # already-created HDF5 still presents a week-old directory.
    # The size walk visits every file, so one pass yields both.
    (bytes, newest) = tree_stats(target)
    if older_than is not None and days_since(newest) < older_than:
        return Skipped(run_dir, label, "too recent")

    return Reclaimable(run_dir, label, target, bytes)


function reclaim_one_dir(run_dir, scratch_root, policy = None,
        other_targets = [], older_than = None, label = None):
    # Layer (b)'s entry point: judge ONE finished run and remove
    # its scratch when the policy calls it spent.  Returns the plan
    # record, carrying whether the removal actually happened and
    # the message when it did not, so a caller can report a prune
    # the way the standalone tool reports a sweep.
    #
    # Nothing about flights appears here.  The caller decides WHEN
    # to call it -- on a unit landing -- and WHICH policy applies;
    # this is only the mechanism (DESIGN 6.2.12).  The default is
    # the workspace policy, since a unit landing in a flight is
    # what layer (b) exists for.
    if policy is None:
        policy = default_reclaim_policy
    record = plan_one_dir(run_dir, label or basename(run_dir),
                          scratch_root, policy, other_targets,
                          older_than)
    if not record.ok:
        return record
    (removed, freed, failures) = apply_reclamation([record])
    record.removed = (removed == 1)
    # apply_reclamation returns (unit, message) pairs, so the
    # message is the second element of the one failure a
    # single-directory plan can produce.
    record.failure = (failures[0][1] if failures else None)
    return record


function plan_reclamation(root, scratch_root, policy = None,
        ids = None, calc_pattern = None, older_than = None,
        kind = None, match = None):
    # Build the plan; remove NOTHING.  Planning and applying are
    # separate calls rather than one function with an `apply`
    # flag, because that split is what lets the preview, the
    # standalone run, and the producer's --clean-after all share
    # one plan and one code path (DESIGN 6.2.12).
    #
    # The filters are concrete, not a general predicate: they are
    # what a user selects on at the command line.  A workspace
    # selects on `ids` and `calc_pattern`; a job tree, which has
    # neither concept, selects on `match` -- a glob over the run
    # directory's path relative to the root.
    #
    # The root's kind fixes the contract for the whole call, and
    # with it the policy that fits when none is given.
    if kind is None:
        kind = detect_root_kind(root)
    if kind is None:
        return []            # neither root: nothing to reclaim,
                             # and no contract under which to try
    if policy is None:
        policy = (default_reclaim_policy if kind == "workspace"
                  else hand_run_policy)

    if kind == "workspace":
        base = join(root, "wingbeats")
        found = [("run", d) for d in find_run_dirs(root)]
    else:
        base = root
        found = find_job_run_dirs(root)   # ("run"|"workspace", d)

    # Resolve EVERY run directory's scratch before judging any of
    # it.  Scratch mirrors the run path, so a run nested in
    # another has its scratch nested too, and refusal 5 needs the
    # whole set -- including runs this call filtered out, since an
    # excluded inner run is exactly the one an outer removal would
    # take as collateral.
    plan = []
    resolved = []                    # (directory, relative, target)
    for (item_kind, directory) in sorted(found by directory):
        relative = relpath(directory, base)
        # A workspace a job-tree walk declined to enter is not a
        # candidate; it is recorded so the report shows where the
        # untouched bytes went (REFUSAL 7).
        if item_kind == "workspace":
            plan.append(SkippedWorkspace(directory, relative))
            continue
        (target, _) = scratch_target(directory, scratch_root)
        resolved.append((directory, relative, target))

    every_target = [t for (_, _, t) in resolved if t != None]

    for (directory, relative, _) in resolved:
        # A workspace selects on unit id and calc tag; a job tree
        # on the path relative to the root.
        if not selected(relative, kind, ids, calc_pattern, match):
            continue

        # Every refusal from here down is plan_one_dir's, so the
        # sweep and the in-flight prune apply one set of rules.
        # `every_target` is passed as the comparison set for
        # REFUSAL 5 -- every run the WALK found, not the selected
        # subset, so a filtered-out inner run still stops the outer
        # removal, which is precisely the case the refusal exists
        # to stop.
        #
        # The link is resolved twice: once above, to build that
        # comparison set, and once inside plan_one_dir.  Those are
        # two cheap stats, and paying them keeps plan_one_dir
        # usable on its own, with no caller obliged to hand it a
        # pre-resolved target.
        plan.append(plan_one_dir(directory, relative, scratch_root,
                                 policy, every_target, older_than))
    return plan


function selected(relative, kind, ids, calc_pattern, match):
    # The CLI filters, applied to one run directory.  A workspace
    # selects on the concepts its layout supplies -- the stable id
    # (first path component under wingbeats/) and the calc tag
    # (the directory's own name).  A job tree has neither, so it
    # selects on the whole relative path, the only handle its
    # free-form layout offers.
    if kind == "workspace":
        if ids != None and first_component(relative) not in ids:
            return False
        if calc_pattern != None and
                not glob_match(basename(relative), calc_pattern):
            return False
        return True
    return match is None or glob_match(relative, match)


function apply_reclamation(plan):
    # Remove the scratch TREE only.  The run directory is the
    # record of the calculation and is never touched, and the
    # `intermediate` link is left dangling on purpose so the run
    # still shows where its scratch was (DESIGN 6.2.12).
    #
    # A tree that will not remove -- a permission problem, a busy
    # filesystem -- is COLLECTED and reported, not raised: one
    # stuck directory must not abandon the rest of a campaign's
    # reclamation.
    removed = 0;  freed = 0;  failures = []
    for item in plan where item is Reclaimable:
        try:
            remove_tree(item.target)
        except OSError as exc:
            failures.append((item.unit, str(exc)));  continue
        removed += 1;  freed += item.bytes
    return (removed, freed, failures)
```

The makeinput-side twin of §12.  It turns `makeinput.py`
from an argv-and-cwd-bound script into one that also
exposes a callable `build_run_dir`, with the CLI a thin
wrapper, so `imago.run_structure` (§12.3) finally has an
in-process makeinput entry point to drive.  The pieces:
the `ScriptSettings` split (14.1); the build orchestration
with cwd discipline (14.2); the thin CLI wrapper (14.3);
and the `run_structure` body it unblocks (14.4).

The governing rules mirror §12.  Build-level faults
(a malformed skl, an element with no basis) are makeinput's
existing behavior, unchanged.  *Contract* faults (the
environment is unconfigured, the structure file is missing,
the run dir cannot be created) raise a `MakeinputError`
(the analog of `ImagoError`) instead of calling `sys.exit`,
so they cannot kill a long-lived kaleidoscope worker
(DESIGN 6.3.1).  And the cwd is a resource acquired and
released around the build (DESIGN 6.3.4).

### 14.1 ScriptSettings split (DESIGN 6.3.3)

The constructor loads rc defaults only; two builders supply
the `args` namespace that the existing `reconcile()`
consumes, exactly as §12.3 splits imago's settings.

```
function ScriptSettings.from_command_line(argv):
    # The CLI path: today's behavior, unchanged in meaning.
    s = ScriptSettings()            # rc defaults only
    args = s.parse_command_line(argv)   # argparse surface
    s.reconcile(args)
    return s

function ScriptSettings.from_options(options):
    # The API path: the same reconciled settings from a
    # plain dict, with no argv and no command-file side
    # effect (record_clp is CLI-only, 14.3 / DESIGN 6.3.5).
    s = ScriptSettings()            # rc defaults only
    args = build_args_namespace(options)
    s.reconcile(args)
    return s
```

```
function build_args_namespace(options):
    # Turn the options mapping into the SAME args namespace
    # argparse would have produced, so reconcile cannot tell
    # which builder called it.  Keys are the argparse `dest`
    # names (job, edge, basis, scfkp, pscfkp, kp, potdb,
    # basisdb, reduce, target, block, xanes, ...).
    args = empty_namespace()
    for dest in ALL_ARGPARSE_DESTS:
        # Absent keys take the argparse default for that
        # dest, so an empty options dict reproduces a bare
        # `makeinput` invocation.
        args[dest] = options.get(dest, argparse_default(dest))
    return args
```

The one subtlety (resolves the DESIGN 6.3.7 open detail):
the **multi-valued flags** -- `reduce`, `target`, `block`
(argparse `action="append"`) and `xanes` (`nargs=
REMAINDER`) -- are repeatable token lists on the command
line, and `reconcile` already turns each into its parsed
form via `_parse_reduce` / `_parse_target` / `_parse_block`
/ `_parse_xanes`.  `from_options` therefore expects the
client to supply each as the *same list-of-token-lists
shape argparse yields* (e.g. `options["reduce"] = [["0.3",
"...","..."], ...]`), and `build_args_namespace` places it
under `args.reduce` verbatim.  The default for an absent
multi-valued flag is the argparse default (`None` or `[]`),
so reconcile's existing "skip when empty" logic applies
unchanged.  This keeps a dict-described run and a
flag-described run byte-identical after reconcile.

### 14.2 The build orchestration (DESIGN 6.3.2, 6.3.4)

`main()`'s body becomes a callable `build_inputs`, and
`build_run_dir` wraps it with structure staging and the
cwd discipline.

```
function build_inputs(settings, sc):
    # The exact sequence today's main() runs inline, minus
    # argv/exit handling.  One definition shared by the CLI
    # and the API so they cannot drift (DESIGN 6.3.4).
    setup_environment(settings)
    initialize_cell(settings, sc)        # reads imago.skl
                                         #   from the cwd
    assign_group(settings, sc, "species")
    assign_group(settings, sc, "types")
    if settings.xanes == 1:
        assign_xanes_types(settings, sc)
    if settings.emu == 1:
        initialize_emu(settings, sc)
    print_imago(settings, sc)            # writes inputs/...
    print_summary(settings, sc)
```

```
function build_run_dir(structure, options, run_dir,
                       settings = None):
    # Build the staged Imago inputs in run_dir from a
    # structure + makeinput options, then return run_dir so
    # a caller can chain into run_prepared (§12.3).
    if settings is None:
        settings = ScriptSettings.from_options(options)

    # Contract checks raise (DESIGN 6.3.1), never sys.exit.
    require_contract(env("IMAGO_RC") or local_makeinputrc(),
        "makeinput environment not configured")
    require_contract(exists(structure),
        "structure file not found: " + structure)

    run_dir = abspath(run_dir)
    makedirs(run_dir, exist_ok = True)

    # Stage the skeleton as run_dir/imago.skl, because
    # makeinput reads the relative name "imago.skl" from the
    # cwd (initialize_cell).  A no-op when structure already
    # IS run_dir/imago.skl.
    staged_skl = join(run_dir, "imago.skl")
    if abspath(structure) != staged_skl:
        copy_file(structure, staged_skl)

    # cwd discipline: acquire the cwd for the build and
    # restore it on EVERY exit, so a failed build cannot
    # strand a flight worker in run_dir (DESIGN 6.3.4).
    original_cwd = getcwd()
    sc = StructureControl()
    try:
        chdir(run_dir)
        build_inputs(settings, sc)
    finally:
        chdir(original_cwd)
    return run_dir
```

Note `build_run_dir` takes no lock: makeinput is a pure
input-staging step writing only into its own `run_dir`,
and the per-run-dir lock that guards concurrent execution
is taken later by `_run_core` (§12.4 / DESIGN 6.1.5).  When
`run_structure` calls the two in sequence, the lock-free
build and the locked run each acquire and release the cwd
around their own scope, so they compose cleanly.

### 14.3 The CLI wrapper (DESIGN 6.3.2, 6.3.5)

`main()` becomes the only layer that touches argv or exits.

```
function cli_main(argv):
    # 1. Parse argv into settings (today's surface).
    settings = ScriptSettings.from_command_line(argv)
    settings.record_clp(argv)   # append argv to `command`;
                                #   CLI-only (DESIGN 6.3.5)
    # 2. Build the cwd as the run dir, holding imago.skl --
    #    today's only behavior, now through the API.
    try:
        build_run_dir("imago.skl", options = {},
                      run_dir = getcwd(), settings = settings)
    except MakeinputError as e:
        log_runtime(e.message)
        return 1            # preserve today's diagnostics
    return 0
```

`record_clp` moves out of the constructor and is called
only here (DESIGN 6.3.5): in API mode there is no
meaningful argv, so `from_options` records the resolved
options as provenance or skips the `command` file -- an
implementation detail with no bearing on the produced
inputs.  The `_load_rc` `sys.exit` on a missing `$IMAGO_RC`
likewise becomes a raised `MakeinputError` the wrapper
catches.

### 14.4 run_structure, completed (DESIGN 6.3.6)

With 14.2 in place, `imago.run_structure` (already shown in
§12.3) is the seam that joins the two APIs:

```
function run_structure(structure, options, run_dir,
                       settings = None):
    import makeinput              # local: imago.py imports
                                  #   without makeinput's env
    if settings is None:
        settings = ScriptSettings.from_options(options)
    makeinput.build_run_dir(structure, options, run_dir)
    return run_prepared(run_dir, settings = settings)
```

The default wingbeat (§13.2) no longer calls this combined
form: it partitions a unit's options and calls `build_run_dir`
and `run_prepared` itself (DESIGN 6.2.10), so the `options`
reaching `build_run_dir` here are makeinput-only.
`run_structure` remains the one-call convenience path for a
*direct* caller that already holds makeinput-only options --
still the shape the C48.3 producer's seam is built on.

## 15. Historical Guidance Dataspace (DESIGN 7)

The accumulation prong: a dataspace of converged
calculations plus a small two-stage k-NN predictor that
turns a new system's chemistry into a predicted converged
k-density and an uncertainty, so a new flight verifies a
small grid around the prediction instead of scanning a wide
one (DESIGN 7.1).  Four blocks, helpers first then drivers:
the file-format-and-predictor library `guidance_db.py`
(15.1 shapes, 15.2 signatures, 15.3 reader, 15.4 emitter,
15.5 predictor); the flight-builder helper inside
`src/scripts/kaleidoscope/` (15.6); the harvest and curator
producers `guidance_harvest.py` / `guidance_promote.py`
(15.7).  All Python under `src/scripts/`; the Fortran side
changes only to expose gap/spin/dos in `result.toml`
(TODO C76).

The governing discipline (VISION Principle 11): the
dataspace is a curated artifact, not tribal knowledge.  The
library reads and validates; the producers stage and
promote; the consumer (the kaleidoscope helper) predicts.
The element-group classification lives in a checked-in data
file (`elemental_groups.toml`), never hardcoded.

### 15.1 Constants and in-memory shapes (DESIGN 7.4)

```
# Schema + partition constants.
SCHEMA_VERSION              = 1
VALID_SYSTEM_TYPES          = ("crystalline", "amorphous",
                               "nanostructure", "molecular")
NON_CRYSTALLINE_TYPES       = ("amorphous", "nanostructure",
                               "molecular")
VALID_BASES                 = ("mb", "fb", "eb")
VALID_GAP_KINDS             = ("direct", "indirect", "none")
METRIC_REGISTRY             = ("total_energy",)  # 7.2 rule 10

# Canonical slot orderings (DESIGN 7.4).  These pin which
# vector slot means which group / Bravais family, so the
# reader, emitter, compute_signature, and predictor all
# agree.  The composition vector sums to 1.0 (7.2 rule 4).
CANONICAL_GROUP_ORDER       = (   # 13 element groups
    "alkali", "alkali_earth", "halide", "chalcogen",
    "pnictogen", "group_iv", "group_iii", "transition_metal",
    "lanthanide", "actinide", "metalloid", "noble_gas",
    "hydrogen")
CANONICAL_LATTICE_ORDER     = (   # 6 Bravais families
    "cubic", "hex", "tet", "ortho", "mono", "tri")

# Predictor tuning knobs (DESIGN 7.6).  All named here so a
# post-seed-flight calibration is a one-file change.

# How many entries a sub-model needs before it is trusted, and
# how many neighbors each k-NN stage averages over.
min_submodel_entries = 3
neighbor_count       = 5

# Keeps the inverse-distance weight finite when a neighbor sits
# exactly on the query point.
distance_floor       = 1e-6

# Relative weights of the two terms in each stage's distance:
# d1 compares chemistry, d2 the electronic character d1 predicts.
composition_weight    = 1.0
lattice_family_weight = 0.25
gap_weight            = 1.0
magnetization_weight  = 0.5

# The physical scale each d2 term is measured against, so a gap
# difference in eV and a moment difference in Bohr magnetons per
# atom become comparable dimensionless numbers.
gap_distance_scale           = 1.0   # eV
magnetization_distance_scale = 0.5   # Bohr magnetons per atom

# The spread at which a stage's confidence has fallen to 1/e: a
# tight neighborhood means a trustworthy prediction, so the
# spread of the neighbors' values is what confidence measures.
gap_confidence_scale            = 1.0    # eV
kpoint_density_confidence_scale = 50.0
```

The dataclasses mirror DESIGN 7.4 exactly; restated here in
field order so the reader (15.3) and emitter (15.4) have a
single target.  `Verification.grid_energies` is the array
the harvest records so the curator's auto-promote rule
(15.7) reads flatness from a staging file alone.

```
dataclass Signature:        # the predictor's feature input
    system_type        : str            # one of the four
    composition_vector : tuple[float]   # 13, group order
    lattice_family     : str            # "" if non-crystalline
    lattice_onehot     : tuple[float]   # 6, lattice order;
                                        #   all zeros if non-
                                        #   crystalline

dataclass Measured:
    gap_ev              : float
    gap_kind            : str           # direct|indirect|none
    spin_polarization   : float
    total_magnetization : float
    kpoint_density      : float         # predictor target

dataclass Context:
    basis                        : str  # mb|fb|eb
    functional                   : str  # e.g. "gga-pbe"
    kpoint_integration           : str  # e.g. "gaussian-0.1"
    scf_threshold                : float
    cell_atom_count              : int
    cell_volume_per_formula_unit : float   # Bohr^3

dataclass Verification:
    grid_values            : tuple[float]
    grid_energies          : tuple[float] | None  # parallel
                                          #   to grid_values
    converged_at           : float
    converged_mesh         : tuple[int] | None    # resolved axial
                                          #   counts of the
                                          #   converged rung
                                          #   (3.12.4); None on
                                          #   pre-mesh / curator
                                          #   entries
    metric                 : str          # "total_energy"
    metric_threshold       : float
    predictor_confidence   : float        # [0.0, 1.0]
    predictor_neighbor_ids : tuple[str]
    gap_spread             : float | None = None
                                          #   how far measured.gap_ev
                                          #   still moves with the
                                          #   mesh AT the converged
                                          #   rung, as a fraction of
                                          #   it (7.2).  Recorded,
                                          #   never acted on.  Last
                                          #   and defaulted because
                                          #   it is OPTIONAL: a
                                          #   hand-written entry need
                                          #   not carry one.  None
                                          #   means NOT MEASURED,
                                          #   never "settled".

dataclass Provenance:
    flight_id        : str
    source_structure : str
    imago_commit     : str
    curator          : str

dataclass GuidanceEntry:
    entry_id     : str
    generated_at : str               # ISO-8601 UTC
    source       : str               # flight|manual
    signature    : Signature
    measured     : Measured
    context      : Context
    verification : Verification | None   # None only for
                                         #   source=manual
    provenance   : Provenance

dataclass Dataspace:
    schema_version         : int
    entries_by_system_type : dict        # system_type ->
                                         #   list[GuidanceEntry]
    group_table            : dict        # symbol -> group name

dataclass PredictionResult:    # what predict() returns
    predicted_kpoint_density : float
    confidence               : float     # [0.0, 1.0]
    is_under_trained         : bool
    neighbor_entry_ids       : tuple[str]
    predicted_gap            : float | None  # None if non-
    predicted_magnetization  : float | None  #   crystalline
#                              (intensive moment, muB/atom)
```

### 15.2 Element groups and compute_signature (DESIGN 7.4)

`elemental_groups.toml` is the checked-in element-to-group table
(Principle 11).  The loader inverts it into a symbol ->
group dict and refuses an element that lands in two groups
(a data-file typo must fail loudly, not silently win the
last assignment).

```
function load_elemental_groups(path):
    raw = tomllib.load(path)
    require(raw["schema_version"] == SCHEMA_VERSION, path,
        "elemental_groups.toml schema_version != "
        + str(SCHEMA_VERSION))
    table = {}                       # symbol -> group name
    for group in CANONICAL_GROUP_ORDER:
        # Every group key must be present (even metalloid,
        # which ships empty per DESIGN 7.4 / 7.10).
        require(group in raw["groups"], path,
            "elemental_groups.toml missing group: " + group)
        for symbol in raw["groups"][group]:
            require(symbol not in table, path,
                "element " + symbol + " assigned to two"
                + " groups (" + table.get(symbol, "?")
                + ", " + group + ")")
            table[symbol] = group
    return table
```

`compute_signature` turns a `StructureControl` into the
predictor's feature input.  Composition is atom-fraction
weighted across the 13 groups; lattice family is read off
the structure's Bravais detection (crystalline only).  An
element symbol missing from the table is a hard error here
-- at compute time -- so the message names the offending
structure, not the dataspace load (DESIGN 7.4).

```
function compute_signature(structure, system_type,
                           group_table):
    require(system_type in VALID_SYSTEM_TYPES,
        "unknown system_type: " + system_type)

    # Composition vector: count atoms per group, normalize
    # to atom-fraction, lay out in CANONICAL_GROUP_ORDER.
    counts = { g: 0 for g in CANONICAL_GROUP_ORDER }
    total_atoms = 0
    for site in structure.atom_sites:
        symbol = element_symbol_of(site)
        require(symbol in group_table,
            "element " + symbol + " (in structure "
            + structure.name + ") not in elemental_groups.toml")
        counts[group_table[symbol]] += 1
        total_atoms += 1
    require(total_atoms > 0, "structure has no atoms")
    composition = tuple(
        counts[g] / total_atoms for g in CANONICAL_GROUP_ORDER)

    # Lattice family + one-hot (crystalline only).  For non-
    # crystalline the family is "" and the one-hot all zeros
    # (DESIGN 7.4); the predictor's stage-1 distance then
    # never sees a lattice term for those system_types.
    if system_type == "crystalline":
        family = bravais_family_of(structure)   # 15.2 note
        require(family in CANONICAL_LATTICE_ORDER,
            "unrecognized lattice family: " + family)
        onehot = tuple(
            1.0 if f == family else 0.0
            for f in CANONICAL_LATTICE_ORDER)
    else:
        family = ""
        onehot = tuple(0.0 for f in CANONICAL_LATTICE_ORDER)

    return Signature(
        system_type        = system_type,
        composition_vector = composition,
        lattice_family     = family,
        lattice_onehot     = onehot)
```

`bravais_family_of(structure)` maps the structure's detected
crystal system to one of the six families.  It reuses the
StructureControl's existing space-group / Bravais detection
(the same machinery the applySpaceGroup path drives).  The
v1 mapping lumps trigonal into `hex`:

```
CRYSTAL_SYSTEM_TO_FAMILY = {
    "triclinic":    "tri",   "monoclinic":   "mono",
    "orthorhombic": "ortho", "tetragonal":   "tet",
    "trigonal":     "hex",   "hexagonal":    "hex",
    "cubic":        "cubic" }
```

The trigonal -> hex lumping (six families, not seven) is a
v1 simplification flagged in DESIGN 7.10; isolating it in
this one table makes a future split a one-line change.

### 15.3 TOML reader load() (DESIGN 7.2 rules 1-12)

`load(root)` reads the `SCHEMA_VERSION` marker, the
`elemental_groups.toml` table, and every entry under
`entries/<system_type>/`, validating each against the 12
rules and partitioning by system_type.  Like
`initial_potential_db.load` (§11.1), every failure names the
file, block, and field at fault.

```
function load(root):
    # Rule 1 (marker half): the bare-integer marker file.
    marker = strip(read_file(join(root, "SCHEMA_VERSION")))
    require(marker == str(SCHEMA_VERSION),
        join(root, "SCHEMA_VERSION"),
        "marker " + marker + " != " + str(SCHEMA_VERSION))

    group_table = load_elemental_groups(
        join(root, "elemental_groups.toml"))

    entries_by_type = { t: [] for t in VALID_SYSTEM_TYPES }
    seen_ids = {}                       # entry_id -> path
    for system_type in VALID_SYSTEM_TYPES:
        subdir = join(root, "entries", system_type)
        if not exists(subdir):
            continue
        for path in sorted(glob(subdir, "*.toml")):
            entry = load_entry(path, system_type, seen_ids)
            entries_by_type[system_type].append(entry)

    return Dataspace(
        schema_version         = SCHEMA_VERSION,
        entries_by_system_type = entries_by_type,
        group_table            = group_table)
```

`load_entry` is the per-file validator.  It checks the
schema BEFORE building the dataclass (rule 12) so an
omission surfaces as a clear validation failure, not a bare
constructor TypeError -- the same discipline as DESIGN 5.2
rule 3.

```
function load_entry(path, system_type_dir, seen_ids):
    raw = tomllib.load(path)

    # Rule 12 (top-level half): required top-level keys.
    for f in ("schema_version", "entry_id", "generated_at",
              "source"):
        require(f in raw, path, "missing top-level: " + f)

    # Rule 1 (entry half): version agrees with the marker.
    require(raw["schema_version"] == SCHEMA_VERSION, path,
        "schema_version " + str(raw["schema_version"])
        + " != " + str(SCHEMA_VERSION))

    # Rule 11: source domain + provenance/verification
    # coupling.
    source = raw["source"]
    require(source in ("flight", "manual"), path,
        "source must be flight|manual, got " + source)

    # Rule 2: entry_id unique across the whole entries tree.
    eid = raw["entry_id"]
    require(eid not in seen_ids, path,
        "duplicate entry_id " + eid + " (also in "
        + seen_ids.get(eid, "?") + ")")
    seen_ids[eid] = path

    # --- signature block -------------------------------
    require("entry" in raw and "signature" in raw["entry"],
        path, "missing [entry.signature]")
    sig = raw["entry"]["signature"]
    require("system_type" in sig, path,
        "missing signature.system_type")
    st = sig["system_type"]

    # Rule 3: system_type valid AND matches the directory
    # the file lives under.
    require(st in VALID_SYSTEM_TYPES, path,
        "invalid system_type: " + st)
    require(st == system_type_dir, path,
        "system_type " + st + " under entries/"
        + system_type_dir + "/")

    # Rule 4: composition vector has exactly the 13 keys,
    # each in [0,1], summing to 1.0 +/- 1e-6.
    require("composition_vector" in sig, path,
        "missing signature.composition_vector")
    cv = sig["composition_vector"]
    require(set(cv.keys()) == set(CANONICAL_GROUP_ORDER),
        path, "composition_vector keys != the 13 groups")
    composition = tuple(cv[g] for g in CANONICAL_GROUP_ORDER)
    for g, x in zip(CANONICAL_GROUP_ORDER, composition):
        require(0.0 <= x <= 1.0, path,
            "composition_vector[" + g + "] out of [0,1]")
    require(abs(sum(composition) - 1.0) <= 1e-6, path,
        "composition_vector sums to "
        + str(sum(composition)) + " != 1.0")

    # Rule 5: lattice_family present+valid iff crystalline.
    family = sig.get("lattice_family", "")
    if st == "crystalline":
        require(family in CANONICAL_LATTICE_ORDER, path,
            "crystalline entry: lattice_family must be one"
            + " of the six families, got '" + family + "'")
        onehot = tuple(
            1.0 if f == family else 0.0
            for f in CANONICAL_LATTICE_ORDER)
    else:
        require(family == "", path,
            "non-crystalline entry must not set"
            + " lattice_family")
        onehot = tuple(0.0 for f in CANONICAL_LATTICE_ORDER)

    signature = Signature(st, composition, family, onehot)

    # --- measured block --------------------------------
    require("measured" in raw["entry"], path,
        "missing [entry.measured]")
    m = raw["entry"]["measured"]
    for f in ("gap_ev", "gap_kind", "spin_polarization",
              "total_magnetization", "kpoint_density"):
        require(f in m, path, "missing measured." + f)

    # Rule 6: gap_ev >= 0; gap_kind valid; none iff metal.
    require(m["gap_ev"] >= 0.0, path, "gap_ev < 0")
    require(m["gap_kind"] in VALID_GAP_KINDS, path,
        "invalid gap_kind: " + m["gap_kind"])
    is_metal = (m["gap_ev"] == 0.0)
    require((m["gap_kind"] == "none") == is_metal, path,
        "gap_kind=='none' iff gap_ev==0.0 violated")
    # Rule 7: kpoint_density > 0.
    require(m["kpoint_density"] > 0.0, path,
        "kpoint_density must be > 0")

    measured = Measured(
        gap_ev              = m["gap_ev"],
        gap_kind            = m["gap_kind"],
        spin_polarization   = m["spin_polarization"],
        total_magnetization = m["total_magnetization"],
        kpoint_density      = m["kpoint_density"])

    # --- context block ---------------------------------
    require("context" in raw["entry"], path,
        "missing [entry.context]")
    c = raw["entry"]["context"]
    for f in ("basis", "functional", "kpoint_integration",
              "scf_threshold", "cell_atom_count",
              "cell_volume_per_formula_unit"):
        require(f in c, path, "missing context." + f)
    # Rule 8: basis valid; functional + kpoint_integration
    # non-empty.
    require(c["basis"] in VALID_BASES, path,
        "invalid basis: " + c["basis"])
    require(len(c["functional"]) > 0, path,
        "functional must be non-empty")
    require(len(c["kpoint_integration"]) > 0, path,
        "kpoint_integration must be non-empty")
    # Rule 9: cell counts/volumes positive.
    require(c["cell_atom_count"] > 0, path,
        "cell_atom_count must be > 0")
    require(c["cell_volume_per_formula_unit"] > 0.0, path,
        "cell_volume_per_formula_unit must be > 0")
    context = Context(
        c["basis"], c["functional"], c["kpoint_integration"],
        c["scf_threshold"], c["cell_atom_count"],
        c["cell_volume_per_formula_unit"])

    # --- verification block ----------------------------
    # Required for source=flight (rule 11); optional for
    # source=manual.
    verification = None
    if "verification" in raw["entry"]:
        v = raw["entry"]["verification"]
        verification = load_verification(v, measured, path)
    require(verification is not None or source == "manual",
        path, "source=flight requires [entry.verification]")

    # --- provenance block ------------------------------
    require("provenance" in raw["entry"], path,
        "missing [entry.provenance]")
    p = raw["entry"]["provenance"]
    for f in ("flight_id", "source_structure",
              "imago_commit", "curator"):
        require(f in p, path, "missing provenance." + f)
    if source == "flight":
        # Rule 11: flight entries need non-empty source +
        # commit + flight id.
        for f in ("flight_id", "source_structure",
                  "imago_commit"):
            require(len(p[f]) > 0, path,
                "source=flight needs non-empty " + f)
    provenance = Provenance(
        p["flight_id"], p["source_structure"],
        p["imago_commit"], p["curator"])

    return GuidanceEntry(
        entry_id     = eid,
        generated_at = raw["generated_at"],
        source       = source,
        signature    = signature,
        measured     = measured,
        context      = context,
        verification = verification,
        provenance   = provenance)
```

```
function load_verification(v, measured, path):
    # Rule 10: verification internal consistency.
    for f in ("grid_values", "converged_at", "metric",
              "metric_threshold", "predictor_confidence",
              "predictor_neighbor_ids"):
        require(f in v, path, "missing verification." + f)
    grid = v["grid_values"]
    require(grid == sorted(grid), path,
        "grid_values not sorted ascending")
    require(v["converged_at"] in grid, path,
        "converged_at not present in grid_values")
    require(v["converged_at"] == measured["kpoint_density"],
        path, "converged_at != measured.kpoint_density")
    require(v["metric"] in METRIC_REGISTRY, path,
        "unknown metric: " + v["metric"])
    require(0.0 <= v["predictor_confidence"] <= 1.0, path,
        "predictor_confidence out of [0,1]")
    # grid_energies optional; if present, parallel length.
    energies = v.get("grid_energies")
    if energies is not None:
        require(len(energies) == len(grid), path,
            "grid_energies length != grid_values length")
    # converged_mesh optional (absent on pre-mesh / curator
    #   entries); when present it is the three axial counts (7.2).
    mesh = v.get("converged_mesh")
    if mesh is not None:
        require(len(mesh) == 3, path,
            "converged_mesh must be three axial counts")
    # gap_spread optional (absent when the ladder was too short to
    #   measure it).  It is a relative change, so a present one
    #   below zero is corrupt.  NOT bounded above: a gap that halves
    #   between two rungs gives a spread over 1, which is a real
    #   reading rather than an error (7.2).
    spread = v.get("gap_spread")
    if spread is not None:
        require(spread >= 0.0, path, "gap_spread must be >= 0")
    return Verification(
        grid_values            = tuple(grid),
        grid_energies          = (tuple(energies)
                                  if energies is not None
                                  else None),
        converged_at           = v["converged_at"],
        converged_mesh         = (tuple(mesh)
                                  if mesh is not None else None),
        gap_spread             = spread,
        metric                 = v["metric"],
        metric_threshold       = v["metric_threshold"],
        predictor_confidence   = v["predictor_confidence"],
        predictor_neighbor_ids = tuple(
            v["predictor_neighbor_ids"]))
```

### 15.4 Hand-formatted emitter save_entry() (DESIGN 7.5)

Deterministic hand-formatter, same philosophy as §11.2:
fixed block sequence, fixed key order, `%.16e` floats,
float arrays one-per-line with a trailing comma after every
element.  Byte-identical output for a given in-memory entry
so version-control diffs are meaningful.

```
function fmt_float(x):      return format(x, ".16e")
function fmt_string(s):
    # TOML basic string: escape backslash and quote.
    return '"' + s.replace("\\","\\\\").replace('"','\\"')
           + '"'
```

```
function short_sha(flight_id, source_structure,
                   generated_at):
    # DESIGN 7.5 slug guard: first 6 hex of SHA-256 over the
    # three provenance fields concatenated.  Two simultaneous
    # harvests differ in flight_id or source_structure, so
    # their hashes (and files) differ.
    blob = (flight_id + source_structure
            + generated_at).encode("utf-8")
    return sha256_hex(blob)[:6]

function slug_for(entry):
    return (entry.signature.system_type + "-"
            + short_sha(entry.provenance.flight_id,
                        entry.provenance.source_structure,
                        entry.generated_at))
```

```
function save_entry(entry, root):
    # Emit into staging/<system_type>/<slug>.toml.  Refuse a
    # collision (7.2 rule 2 / 7.5): the caller retries with a
    # fresh generated_at on the rare hash clash.
    slug = slug_for(entry)
    subdir = join(root, "staging", entry.signature.system_type)
    make_dirs(subdir)
    path = join(subdir, slug + ".toml")
    require(not exists(path), "save_entry: " + path
        + " already exists (entry_id collision)")
    write_file(path, format_entry(entry, slug))
    return path
```

```
function format_entry(entry, slug):
    out = []
    # Top-level block.  entry_id always equals the slug.
    out.append("schema_version = " + str(SCHEMA_VERSION))
    out.append("entry_id       = " + fmt_string(slug))
    out.append("generated_at   = "
        + fmt_string(entry.generated_at))
    out.append("source         = " + fmt_string(entry.source))
    out.append("")

    # [entry.signature] + the multi-line composition vector.
    out.append("[entry.signature]")
    out.append("system_type    = "
        + fmt_string(entry.signature.system_type))
    if entry.signature.system_type == "crystalline":
        out.append("lattice_family = "
            + fmt_string(entry.signature.lattice_family))
    out.append("")
    out.append("[entry.signature.composition_vector]")
    width = max(len(g) for g in CANONICAL_GROUP_ORDER)
    for g, x in zip(CANONICAL_GROUP_ORDER,
                    entry.signature.composition_vector):
        out.append(pad(g, width) + " = " + fmt_float(x))
    out.append("")

    # [entry.measured].
    out.append("[entry.measured]")
    emit_kv(out, "gap_ev",              entry.measured.gap_ev)
    emit_kv(out, "gap_kind",            entry.measured.gap_kind)
    emit_kv(out, "spin_polarization",
            entry.measured.spin_polarization)
    emit_kv(out, "total_magnetization",
            entry.measured.total_magnetization)
    emit_kv(out, "kpoint_density",
            entry.measured.kpoint_density)
    out.append("")

    # [entry.context].
    out.append("[entry.context]")
    emit_kv(out, "basis",      entry.context.basis)
    emit_kv(out, "functional", entry.context.functional)
    emit_kv(out, "kpoint_integration",
            entry.context.kpoint_integration)
    emit_kv(out, "scf_threshold", entry.context.scf_threshold)
    emit_kv(out, "cell_atom_count",
            entry.context.cell_atom_count)
    emit_kv(out, "cell_volume_per_formula_unit",
            entry.context.cell_volume_per_formula_unit)
    out.append("")

    # [entry.verification] when present.  grid_values and
    # grid_energies are one-float-per-line, trailing comma.
    if entry.verification is not None:
        v = entry.verification
        out.append("[entry.verification]")
        emit_float_array(out, "grid_values", v.grid_values)
        if v.grid_energies is not None:
            emit_float_array(out, "grid_energies",
                             v.grid_energies)
        emit_kv(out, "converged_at",     v.converged_at)
        # converged_mesh: an inline int array (the resolved axial
        #   counts), beside converged_at; omitted when None (7.2).
        if v.converged_mesh is not None:
            out.append("converged_mesh = ["
                + join_csv(str(n) for n in v.converged_mesh)
                + "]")
        # gap_spread: how settled the recorded gap was at that same
        #   rung (7.2).  Omitted when None, so an absent field reads
        #   as NOT MEASURED rather than as a settled gap.
        if v.gap_spread is not None:
            emit_kv(out, "gap_spread", v.gap_spread)
        emit_kv(out, "metric",           v.metric)
        emit_kv(out, "metric_threshold", v.metric_threshold)
        emit_kv(out, "predictor_confidence",
                v.predictor_confidence)
        out.append("predictor_neighbor_ids = ["
            + join_csv(fmt_string(i)
                       for i in v.predictor_neighbor_ids)
            + "]")
        out.append("")

    # [entry.provenance].
    out.append("[entry.provenance]")
    emit_kv(out, "flight_id",      entry.provenance.flight_id)
    emit_kv(out, "source_structure",
            entry.provenance.source_structure)
    emit_kv(out, "imago_commit",     entry.provenance.imago_commit)
    emit_kv(out, "curator",          entry.provenance.curator)

    return "\n".join(out) + "\n"
```

```
function emit_kv(out, key, value):
    # Render one key = value line.  Floats use %.16e; ints
    # bare; strings TOML-quoted.  Block-internal alignment
    # follows the 7.3 sketch's hand-aligned '=' columns.
    if value is a float:    text = fmt_float(value)
    elif value is an int:   text = str(value)
    else:                   text = fmt_string(value)
    out.append(key + " = " + text)

function emit_float_array(out, key, values):
    out.append(key + " = [")
    for x in values:
        out.append("    " + fmt_float(x) + ",")
    out.append("]")
```

(The exact `=`-column alignment per block matches the
DESIGN 7.3 gold sketch; a tiny `pad`/width pass like
§11.2's `format_block` produces it.  Omitted here for
brevity -- the byte-determinism that matters comes from the
fixed key order, the fixed float format, and the
one-element-per-line arrays.)

### 15.5 Predictor predict() (DESIGN 7.6)

`predict` switches on system_type: non-crystalline returns
the canonical entry; crystalline runs the two-stage k-NN.
It always returns a `PredictionResult` (never None); the
`is_under_trained` flag plus `confidence` tell the caller
how seriously to take it (DESIGN 7.4).

```
function predict(dataspace, query, basis, functional,
                 kpoint_integration):
    pool = dataspace.entries_by_system_type.get(
               query.system_type, [])

    if query.system_type in NON_CRYSTALLINE_TYPES:
        return predict_non_crystalline(pool)

    entries, under_trained = select_submodel(
        pool, basis, functional, kpoint_integration)
    if under_trained:
        # No usable sub-model: the caller (15.6) falls back
        # to the wide-grid default (DESIGN 7.9).  The
        # density field is unused in this branch.
        return PredictionResult(
            predicted_kpoint_density = 0.0,
            confidence               = 0.0,
            is_under_trained         = True,
            neighbor_entry_ids       = (),
            predicted_gap            = None,
            predicted_magnetization  = None)

    pgap, pmag, conf1, n1 = stage1(query, entries)
    pkpd, conf2, n2       = stage2(pgap, pmag, entries)
    return PredictionResult(
        predicted_kpoint_density = pkpd,
        confidence               = conf1 * conf2,
        is_under_trained         = False,
        neighbor_entry_ids       = dedup(n1 + n2),
        predicted_gap            = pgap,
        predicted_magnetization  = pmag)
```

```
function predict_non_crystalline(pool):
    # k-density is set by the cell-volume convention, not
    # chemistry, so the single hand-seeded canonical entry
    # (DESIGN 7.9) captures essentially all the signal.
    canon = [e for e in pool if e.source == "manual"]
    if len(canon) == 0:
        return PredictionResult(
            0.0, 0.0, True, (), None, None)  # under-trained
    # Day-1 there is exactly one; if several accumulate, the
    # most recent canonical wins (deterministic by
    # generated_at).
    entry = max(canon, key = lambda e: e.generated_at)
    return PredictionResult(
        predicted_kpoint_density = entry.measured.kpoint_density,
        confidence               = 1.0,
        is_under_trained         = False,
        neighbor_entry_ids       = (entry.entry_id,),
        predicted_gap            = None,
        predicted_magnetization  = None)
```

Sub-model selection is the
(basis, functional, kpoint_integration) ->
functional-family -> overall-pool fallback chain of DESIGN
7.6 step 2.  `is_under_trained` is set only when even the
overall pool is too thin.

```
function select_submodel(pool, basis, functional,
                         kpoint_integration):
    # 1. Exact (basis, functional, kpoint_integration)
    #    sub-model.
    exact = [e for e in pool
             if e.context.basis == basis
             and e.context.functional == functional
             and e.context.kpoint_integration
                 == kpoint_integration]
    if len(exact) >= min_submodel_entries:
        return exact, False

    # 2. Most-populous (basis, functional) sub-model within
    #    the same functional family (DESIGN 7.6: (mb,gga-pbe)
    #    -> (fb,gga-pbe)).
    fam = functional_family(functional)
    family = [e for e in pool
              if functional_family(e.context.functional) == fam]
    best = most_populous_submodel(family)
    if len(best) >= min_submodel_entries:
        return best, False

    # 3. The whole system_type pool, context ignored.
    if len(pool) >= min_submodel_entries:
        return pool, False

    # 4. Too thin everywhere -> under-trained.
    return pool, True

function functional_family(functional):
    # v1: the token before the first hyphen ("gga-pbe" ->
    # "gga", "lda" -> "lda").  Whether functional/basis are
    # sub-models or features is open (DESIGN 7.10); isolating
    # the rule here makes that a one-line change.
    return functional.split("-")[0]

function most_populous_submodel(entries):
    # Group by (basis, functional); return the largest group.
    groups = group_by(entries,
        key = lambda e: (e.context.basis, e.context.functional))
    if len(groups) == 0:
        return []
    return max(groups.values(), key = len)
```

The two k-NN stages share one inverse-distance-weighted
helper.  Weights are `1 / (d + distance_floor)`, normalized to
sum to 1.0, so an exact match dominates without dividing by
zero.

```
function knn_weights(entries, distance_of):
    # distance_of(entry) -> float >= 0.  Returns the
    # neighbor_count nearest as (entry, weight) pairs.
    scored = sort(entries, key = distance_of)        # ascending
    nearest = scored[: min(neighbor_count, len(scored))]
    raw = [1.0 / (distance_of(e) + distance_floor)
           for e in nearest]
    total = sum(raw)
    return [(e, r / total) for e, r in zip(nearest, raw)]
```

```
function stage1(query, entries):
    # Chemistry -> electronic character.  d1 combines the
    # composition L2 distance with the lattice-family one-hot
    # term, the latter halved so a full mismatch maps to the
    # same [0,1] range as composition (DESIGN 7.6 step 3).
    function d1(e):
        composition_sq = sum_sq(
            sub(query.composition_vector,
                e.signature.composition_vector))
        lattice_sq = sum_sq(
            sub(query.lattice_onehot,
                e.signature.lattice_onehot))
        return sqrt(composition_weight * composition_sq
                    + lattice_family_weight * lattice_sq / 2.0)

    nbrs = knn_weights(entries, d1)
    pgap = sum(w * e.measured.gap_ev   for e, w in nbrs)
    pmag = sum(w * intensive_mag(e)    for e, w in nbrs)
    # Confidence_1 from the weighted gap variance (7.6).
    var = sum(w * (e.measured.gap_ev - pgap) ** 2
              for e, w in nbrs)
    conf1 = exp(-sqrt(var) / gap_confidence_scale)
    return pgap, pmag, conf1, [e.entry_id for e, _ in nbrs]
```

`intensive_mag(e)` is the entry's net moment per atom --
`abs(e.measured.total_magnetization) / e.context.cell_atom_count`
(Bohr magnetons per atom) -- the predictor's spin-character
feature (DESIGN 7.6).  It is per-atom (intensive, so a cell and
its supercell compare equal), taken in magnitude (the up/down
labeling is arbitrary), and built from the measured moment
rather than `spin_polarization` (which imago never surfaces, so
it is always 0.0).

```
function stage2(pgap, pmag, entries):
    # Electronic character -> k-density.  d2 uses the
    # PREDICTED character (from stage 1), not the query's
    # chemistry: "find calcs whose gap-and-magnetization look
    # like what this query is likely to produce"
    # (DESIGN 7.6 step 4).
    function d2(e):
        gap_term = (gap_weight
                    * (pgap - e.measured.gap_ev) ** 2
                    / gap_distance_scale ** 2)
        mag_term = (magnetization_weight
                    * (pmag - intensive_mag(e)) ** 2
                    / magnetization_distance_scale ** 2)
        return sqrt(gap_term + mag_term)

    nbrs = knn_weights(entries, d2)
    pkpd = sum(w * e.measured.kpoint_density for e, w in nbrs)
    var = sum(w * (e.measured.kpoint_density - pkpd) ** 2
              for e, w in nbrs)
    conf2 = exp(-sqrt(var) / kpoint_density_confidence_scale)
    return pkpd, conf2, [e.entry_id for e, _ in nbrs]
```

### 15.6 Builder record and shared primitives (DESIGN 7.7)

These live in `src/scripts/kaleidoscope/` (domain-aware
convenience; the dispatch core stays dumb, Principle 12).  They
are the pieces the producer's convergence builders (`predict_-
kpoint_density` and `build_mesh_unit`, both 4e.7) and the harvest
share: the `PredictionRecord` the predictor fills and the harvest
hook later recovers, the calc-tag encoding a unit's directory name
uses, and the cache-identity key fields.  The mesh-climb search
that drives the builders is 4e; the prediction itself is 4e.7.

```
dataclass PredictionRecord:    # 7.7-derived; serialized as
                               # [flight.predictions.<id>]
    policy                  : str   # trust_no_verify |
                                    #   wide_grid_no_prior |
                                    #   verify_around_prediction |
                                    #   curator_override
    predicted_kpoint_density : float | None
    confidence              : float
    is_under_trained        : bool
    neighbor_entry_ids      : tuple[str]
    predicted_gap           : float | None
    predicted_magnetization : float | None
    system_type             : str
    feature_vector          : Signature
    basis                   : str   # the (basis, functional,
    functional              : str   #   kpoint_integration)
    kpoint_integration      : str   #   sub-model this run used;
                                    #   per-record (not flight-
                                    #   level fixed_axes) so a
                                    #   combined multi-structure
                                    #   flight whose structures
                                    #   differ in sub-model is
                                    #   still harvestable (DESIGN
                                    #   6.2.9 / 7.8 step 3f)
    kpoint_convergence_threshold : float  # the two RESOLVED
    metal_gap_threshold          : float  #   manifest knobs the
                                    #   harvest needs and cannot
                                    #   look up: it is a standalone
                                    #   tool pointed at a finished
                                    #   workspace and never sees a
                                    #   manifest, so the producer
                                    #   stamps both here when it
                                    #   builds the flight (DESIGN
                                    #   7.8 steps 3c / 3d').  The
                                    #   first is the grid-flatness
                                    #   tolerance, eV per atom, from
                                    #   the solid's own value else
                                    #   [harvest] else 2e-3 (5.7);
                                    #   the second is the absolute
                                    #   band gap in eV below which a
                                    #   run counts as metallic, the
                                    #   database-wide
                                    #   [kpoint_climb] knob the
                                    #   climb reads as
                                    #   config.metal_gap_threshold
                                    #   (4e.3), so both harvest
                                    #   paths and the climb apply
                                    #   ONE resolved value
```


```
function encode_axis_value(v):
    # DESIGN 6.2.4 rule 3: '.' -> 'p', leading '-' -> 'm'.
    # The flight builder rounds k-density to an integer
    # first (below), so in v1 this is just the decimal int;
    # the general encoder is kept for future axes.
    if v == round(v):
        text = str(int(round(v)))
    else:
        text = trim_trailing_zeros(repr_compact(v))
    negative = text.startswith("-")
    if negative:
        text = text[1:]
    text = text.replace(".", "p")
    return ("m" + text) if negative else text

function build_calc_tag(calc_axes):
    # calc_axes is an ORDERED mapping {axis: value}; the
    # order must match Flight.sweep.varied_axes.  Returns a
    # tuple of "<axis>-<encoded-value>" directory components
    # (DESIGN 6.2.1 / 6.2.4).  In v1 it is one component.
    components = []
    for axis, value in calc_axes.items():       # insertion order
        require(matches(axis, "^[a-z0-9-]+$"),
            "calc axis name not a slug: " + axis)
        components.append(axis + "-" + encode_axis_value(value))
    return tuple(components)
```


```
# The scalar option keys that define the producer's run identity
#   (DESIGN 6.2.1/6.2.10): just `converg`, the SCF convergence
#   limit, a makeinput option and the concrete name for DESIGN's
#   "scf_threshold".  Taken from a unit's options when present.
#
#   The engine build identity is NOT here.  The key asks whether
#   this is the same CALCULATION, not whether its result is still
#   good (DESIGN 6.2.5): a rebuilt engine does not make a stored
#   potential wrong, since that potential is a starting point every
#   later SCF re-converges, while comparing the build would miss
#   the cache on every ordinary development commit.  It travels on
#   `CalcUnit.record` instead -- recorded per run, printed in the
#   reuse plan, never compared (13.1 / 13.5).
KEY_SCALAR_NAMES = ("converg",)


# The makeinput outputs byte-compared as the producer's key
#   (DESIGN 6.2.5).  BOTH are needed, and the second is not
#   optional polish: `structure.dat` bakes in the type/species
#   assignment, the basis, the functional, and the potential, but
#   NOT the k-point integration scheme, which reaches `kp-scf.dat`
#   as KPOINT_INTG_CODE.  With `structure.dat` alone, one solid at
#   one mesh under two different integration schemes shares a run
#   directory and hits -- returning the other scheme's answer under
#   the name of the one asked for.  `kp-scf.dat` also carries the
#   point operations, so a run that suppresses the mesh reduction
#   while keeping its atomic symmetry is likewise distinguished.
#
#   Adding the scheme to KEY_SCALAR_NAMES instead would invalidate
#   every stored cache_key.toml at once (the scalars are compared
#   as a whole table), which is a mass false miss.  A key FILE
#   costs nothing PROVIDED the path names a file every unit has --
#   which is why both are declared under `inputs/`, where
#   makeinput writes them for every unit whatever its job reads.
#   Naming them at the run-directory root instead reaches them
#   only for jobs that run an SCF, so every fingerprint unit --
#   which runs none -- misses forever and in silence (TODO D23).
#   The paths are relative to the unit's directory and are joined
#   the same way on both sides of the compare (13.4 / 11.4).
KEY_FILE_PATHS = ("inputs/structure.dat", "inputs/kp-scf.dat")


function standard_key_fields(structure, options):
    # DESIGN 6.2.5: the producer's cache identity -- the scalars
    # taken from `options` (the SCF threshold, and nothing else)
    # plus the key files above, byte-compared.  The key files are
    # makeinput's OUTPUTS, not the raw skeleton, so an input that
    # changes the result misses the cache on its own, with no
    # hand-listed "options that matter" to fall stale.  Each
    # KeyFile `source` is provisional here (the skeleton
    # `structure`); the driver's prepare step (11.4) re-points
    # EVERY key file at its built copy once those files exist.
    return KeyFields(
        scalars = { name : options[name]
                    for name in KEY_SCALAR_NAMES
                    if name in options },
        files   = [KeyFile(path = path, source = structure)
                   for path in KEY_FILE_PATHS])
```


This rests on the generic `metadata` field on `Flight` plus
the `serialize_flight` loop that emits each `metadata[key]`
as a `[flight.<key>]` block -- both now defined canonically
in DESIGN 6.2.1 / PSEUDOCODE 13.1.  It keeps the core
domain-agnostic while letting the §7 helper persist the
prediction provenance the harvest needs.


### 15.7 Harvest and promote (DESIGN 7.8)

**`guidance_harvest.py`** turns a finished flight into
staged guidance entries.  It reads the flight workspace
(it is the producer that has workspace access; the curator
later works on staging files alone).

**The three-source rule (Model 1, settled 2026-05-30).**  Each
entry field is filled from exactly ONE of three inputs, so the
information flow stays simple:

- **`flight.toml`** -- the *plan*: the unit list (id, structure,
  calc tags), the `[flight.sweep]` block (which axis varied; in
  v1 nothing is held fixed), and the
  `[flight.predictions.<id>]` tables (one per structure, each
  carrying that structure's prediction AND the (basis,
  functional, kpoint_integration) sub-model it ran under).  Each
  grid point's swept k-density is read out of its
  calc tag (`kpt-density-<int>`) using the sweep's ordered
  `varied_axes` -- the makeinput `options` are deliberately NOT
  persisted in `flight.toml`, so the calc tag is the on-disk
  source of the swept value.
- **each run's `result.toml`** -- the *per-run facts*: the final
  SCF `total_energy`, the measured `gap_ev` / `gap_kind` /
  `total_magnetization`, and the `scf_threshold` the run used
  (imago.py writes all of these; the per-run record is
  self-contained, DESIGN 6.1).
- **the structure `.skl`** -- the *structural facts*: the harvest
  loads it anyway for `compute_signature`, and the same load
  yields `cell_atom_count` (`num_atoms`) and
  `cell_volume_per_formula_unit` (the cell volume in Bohr^3,
  formula-unit count Z = 1 in v1).

Three v1 conventions.  The grid-flatness `metric_threshold` is
the solid's resolved `kpoint_convergence_threshold` -- energy per
atom, in eV (DESIGN 7.8 / 5.7) -- and rides on the structure's
prediction record, a manifest/resolved fact absent from any run's
`result.toml`; it is distinct from the SCF's own `scf_threshold`,
which stays a per-run context fact.  `grid_energies` are stored
RAW -- total-cell energies in hartree, exactly as the runs
produced them -- and every site that compares them against the
per-atom eV `metric_threshold` normalizes at the point of use
(`pick_converged` here, `auto_promote_ok` in the promoter); the
physical values keep the record honest, and `cell_atom_count` is
recorded alongside for the conversion.  `imago_commit` is a
*recorded* fact rather than a measured one: the producer hangs it on
each unit's `record` (11.4), the driver stamps it into `status.toml`
and the wingbeat echoes it into `result.toml` (13.2/13.5), and this
harvest reads it there with the other per-run facts, falling back to
`"unknown"` for a run that recorded none.  `spin_polarization`
is recorded as `0.0` -- imago surfaces the magnetic *moment*, not
a polarization, so the predictor keys its spin character on
`total_magnetization` instead (DESIGN 7.6).

```
function harvest_flight(workspace_root, db_root, dataspace):
    flight = read_flight_toml(
        join(workspace_root, "flight.toml"))
    # Per-structure predictions, keyed by structure id (DESIGN
    # 6.2.9); a single-structure flight carries a one-entry map.
    predictions = flight.metadata.get("predictions", {})  # 15.6

    # The swept axis names which calc-tag component carries each
    #   grid point's value.  A flight that declared no sweep
    #   (sweep is None) cannot feed THIS harvester -- it is the
    #   GUIDANCE harvester, and a guidance entry IS the claim
    #   "this k-density is converged," which only a grid can
    #   establish; with no varied axis there is nothing to read a
    #   swept value along.  A single one-off calculation is not
    #   blocked by this: it is a length-1 SWEEP (trust mode or a
    #   pinned kpoint_spec), still harvested for the producer's
    #   potential deliverable, and merely skipped for guidance
    #   staging at step (c) below (one point is not convergence
    #   evidence); a known-good k-density can instead be seeded by
    #   hand as a source="manual" entry (DESIGN 7.9).
    #   The sub-model (basis/functional/kpoint_integration) is NOT
    #   read here from sweep.fixed_axes: it rides on each
    #   structure's prediction record so a combined
    #   mixed-sub-model flight is harvestable (DESIGN 6.2.9 / 7.8
    #   step 3f); the per-structure read is in step (g) below.
    axis = flight.sweep.varied_axes[0]            # v1: single axis

    # Only convergence-sweep runs are grid points; other kinds
    # (e.g. "fingerprint" loen runs) share a structure id but
    # belong to a different harvester (DESIGN 6.2.9).
    convergence_units = [u for u in flight.units
                         if u.kind == "convergence"]

    for unit_id, units in group_by_id(convergence_units):
        # The prediction this structure was launched under.  It is
        #   the SOLE source of system_type (step f) and the
        #   sub-model (step g), so a structure with no record cannot
        #   be staged -- skip it.  The helper (15.6) always attaches
        #   one; a record-less sweep is a hand-built flight outside
        #   the predict-then-verify path (DESIGN 7.8 / 7.9).
        prediction = predictions.get(unit_id)
        if prediction is None:
            log(unit_id + ": no prediction record (not staged)")
            continue
        # a. The verification sub-grid for this structure, sorted
        #    by swept k-density (read out of each unit's calc tag).
        grid = sort(units, key = lambda u: swept_value_of(u, axis))

        # b. Parse each converged run's result.toml.  meshes[i]
        #    is the run's resolved kpoint_mesh -- the axial counts
        #    imago integrated over (6.1.2 / PSEUDOCODE 4c.6) --
        #    and feeds the duplicate-mesh guard in step (d).  It
        #    is None when result.toml carries no kpoint_mesh (an
        #    older run, or imago not yet emitting it), which the
        #    guard treats as "cannot collapse" (see collapse_by_mesh).
        #    A point whose SCF did not converge is DROPPED here, so
        #    it never reaches the flatness test (DESIGN 7.8 step 3b).
        #    Its energy is wherever the iteration stopped, for a
        #    reason unrelated to the mesh, so it can read flat by
        #    coincidence and can break a plateau that was real.
        #    Dropping (rather than stopping the structure, as the
        #    climb does) is right here because this grid is a fixed
        #    set of points, not a sequence that chooses its next
        #    member -- removing one cannot stall anything.  If too
        #    few survive, step (e)'s "too few distinct meshes" arm
        #    already covers it.  Only an EXPLICIT false drops a
        #    point; a result.toml with no such field cannot be judged
        #    and is kept.  Never silent: the drops are reported.
        kpoint_densities, energies, meshes, rts = [], [], [], []
        dropped = []
        for u in grid:
            rt = read_result_toml(
                join(workspace_root, "wingbeats", u.id, *u.calc,
                     "result.toml"))
            if rt.get("converged") is false:
                dropped.append(swept_value_of(u, axis))
                continue
            kpoint_densities.append(swept_value_of(u, axis))
            energies.append(rt["total_energy"])
            meshes.append(rt.get("kpoint_mesh"))
            rts.append(rt)
        if dropped:
            warn(unit_id + ": dropped " + str(len(dropped))
                 + " grid point(s) whose SCF did not converge: "
                 + str(dropped))

        # c. A single-point grid harvests deliverables but does
        #    NOT auto-stage a guidance entry (DESIGN 6.2.1 / 7.7):
        #    one converged calc is weaker evidence than a grid.
        #    This covers both trust_no_verify and a single-point
        #    curator_override, and MUST precede pick_converged:
        #    the two-sided convergence test below needs >= 3
        #    points and would otherwise misreport one as "energy
        #    still moving."
        #
        #    Counted over the SURVIVING points, not the requested
        #    ones, since everything downstream reads the survivors --
        #    and the zero case has to be caught before scf_threshold
        #    is read off rts[0].
        if len(rts) == 0:
            log(unit_id + ": every grid point failed to converge"
                + " (not staged)")
            continue
        if len(rts) == 1:
            log(unit_id + ": single point (not staged)")
            continue

        # The k-point flatness tolerance rode in on this
        #   structure's prediction record (15.6): per atom, in eV,
        #   the solid's resolved kpoint_convergence_threshold
        #   (DESIGN 7.8 / 5.7).  The SCF's own criterion is a
        #   separate per-run fact kept for context below.
        kpoint_threshold = prediction[
            "kpoint_convergence_threshold"]
        scf_threshold    = rts[0]["scf_threshold"]

        # The per-atom comparison needs the cell size, so load the
        #   structure once here; step (f) reuses it for the
        #   signature and the cell facts.
        sc = load_structure(grid[0].structure)   # read_input_file
        cell_atom_count = sc.num_atoms

        # d. Collapse duplicate-mesh rungs, then pick the
        #    converged grid point (DESIGN 7.8 step 3c).  Two rungs
        #    that resolved to the same mesh are one calculation run
        #    twice; their zero energy delta would fool the two-
        #    sided test, so collapse_by_mesh reduces the grid to one
        #    rung per distinct mesh (keeping the lowest-density
        #    member) before the test.  `kept[j]` maps a collapsed
        #    index back to its original grid position.  The energies
        #    are raw total-cell hartree, so pick_converged
        #    normalizes each delta to eV per atom before comparing.
        c_kpoint_densities, c_energies, kept = collapse_by_mesh(
            kpoint_densities, energies, meshes)
        j = pick_converged(c_energies, cell_atom_count,
                           kpoint_threshold)

        # e. No flat interior point -- energy still moving at the
        #    top of the range, or the grid collapsed below the
        #    three distinct meshes the interior test needs.  Tag
        #    and SKIP: a non-converged sweep earns no entry.
        if j is None:
            warn(unit_id + ": no converged point (energy still"
                 + " moving, or too few distinct meshes) -- skipped")
            tag_prediction_mismatch(workspace_root, unit_id)
            continue
        idx = kept[j]              # original index of chosen rung

        # f/g. Build the rich entry from the already-chosen facts:
        #    the COLLAPSED distinct-mesh ladder (step d), the chosen
        #    rung's k-density, and its result.toml.  ONE entry
        #    builder feeds both harvests (the Q1-Q2 shared core,
        #    DESIGN 5.7): this density sweep -- which PICKED via
        #    collapse_by_mesh + pick_converged above -- and the
        #    producer's in-memory climb (11.4, which picks via the
        #    climb and record_converged) hand build_entry the SAME
        #    shape, so they stage identical entries and cannot drift.
        #
        #    This path has no climb to take a verdict FROM, so it
        #    makes the multi-rung reading itself, over the gaps of
        #    every point in its own grid -- the same any-rung rule
        #    the climb applies to its ladder (4e.3), on the evidence
        #    this path happens to hold.  The gaps were parsed in step
        #    (b) and simply never looked at before.
        ladder_is_metal = any(
            is_gapless_value(rt.get("gap_ev"),
                             prediction["metal_gap_threshold"])
            for rt in rts)
        entry = build_entry(
            workspace_root, grid[0].structure, prediction,
            dataspace, sc, kpoint_threshold,
            c_kpoint_densities, c_energies,
            kpoint_densities[idx], rts[idx],
            ladder_is_metal = ladder_is_metal)

        # g'. A metal builds no entry (DESIGN 7.8).  build_entry
        #    returns None and BOTH harvests skip on it, so the one
        #    place the rule lives is the one builder they share.
        if entry is None:
            log(unit_id + ": metal -- no guidance entry staged")
            continue

        # h. Stage it.  save_entry fills entry_id = slug.
        path = save_entry(entry, db_root)
        log(unit_id + ": staged " + path)
```

`swept_value_of(unit, axis)` reads the value out of the calc
component whose prefix matches `axis` (each component is
`"<axis>-<encoded-value>"`, DESIGN 6.2.4), inverting the
builder's `encode_axis_value` (`"p"` -> decimal point, leading
`"m"` -> minus).  `flight_id_of(workspace_root)` is the
workspace root's basename; both live in
`kaleidoscope.workspace` alongside `read_flight_toml`.

```
function build_entry(workspace_root, source_structure, prediction,
                     dataspace, structure, kpoint_threshold,
                     grid_values, grid_energies, converged_density,
                     chosen_result, ladder_is_metal = false,
                     ladder_gaps = None):
    # Assemble one converged structure's GuidanceEntry (DESIGN 7.8
    # step 3f).  The ONE entry builder both harvests feed (the
    # Q1-Q2 shared core, DESIGN 5.7): the density sweep
    # (harvest_flight, picking via collapse_by_mesh + pick_converged)
    # and the producer's in-memory climb (11.4, picking via the
    # climb + record_converged) each PICK their converged rung, then
    # hand the identical already-chosen facts here --
    #   grid_values / grid_energies  the distinct-mesh flatness
    #                                ladder (ascending), raw total-
    #                                cell hartree (Option B; the
    #                                consumer normalizes per atom)
    #   converged_density            the chosen rung's k-density
    #   chosen_result                the chosen run's result.toml
    # -- so the two paths stage identical entries and cannot drift.
    # Everything measured, the SCF threshold, the exact mesh, and
    # the commit come from chosen_result; the sub-model and
    # system_type from the prediction record (its sole home,
    # 6.2.9); the cell facts and signature from the loaded structure.
    #
    # A METAL BUILDS NO ENTRY -- return None (DESIGN 7.8).  An
    # entry's whole content is "for a structure like this, this
    # k-density is converged," and a metal cannot make that claim:
    # its energy does not converge in k-points, and the climb
    # short-circuits at the FIRST gapless rung and settles there
    # as a deliberately rough potential (DESIGN 3.12.3).  That
    # settled rung is a stopping point, not a converged density,
    # and the predictor would read it as evidence -- often as the
    # only member of its lattice family, and so as the dominant
    # neighbor for every later query in that family (7.6).  The
    # guard sits HERE, in the shared builder, so neither harvest
    # path can grow its own version of the rule.
    #
    # TWO readings, EITHER sufficient (DESIGN 7.8 d').
    #
    # `ladder_is_metal` is the CALLER's multi-rung reading, passed in
    # rather than re-derived here, because only the caller has the
    # ladder: the producer passes the climb's verdict (11.4), the
    # standalone harvest passes the any-rung test over its own grid
    # (below).  It defaults to false, which means "no multi-rung
    # evidence offered" and leaves the chosen-rung test as the only
    # one -- never "known not to be a metal".
    #
    # The chosen rung's own gap is the second reading.  Neither is
    # redundant.  A metal on a discrete mesh shows an artificial gap
    # whose size depends on where the mesh points fall (DESIGN 1.6),
    # so the single chosen rung is weak evidence -- fcc Al reads zero
    # at several meshes and 0.124 eV at another -- while the ladder
    # taken whole is strong.  Taking either as sufficient means the
    # stronger reading cannot be overruled by the weaker one, and
    # nothing the chosen-rung test already caught stops being caught.
    #
    # The cut is read off the prediction record, which is how a
    # manifest knob reaches a standalone tool that never sees a
    # manifest (15.6, the same channel kpoint_threshold uses).  The
    # test itself is `is_gapless_value`, the scalar core the climb's
    # rung-shaped `is_gapless` also calls (4e.2), so one rule serves
    # both -- including its side on missing data: an UNKNOWN gap is
    # not metallic, because a missing reading must never silently
    # suppress an entry a real insulator earned.
    if ladder_is_metal:
        return None
    if is_gapless_value(chosen_result.get("gap_ev"),
                        prediction["metal_gap_threshold"]):
        return None
    system_type = prediction["system_type"]
    sig = compute_signature(structure, system_type,
                            dataspace.group_table)
    return GuidanceEntry(
        entry_id     = "",                    # set by save_entry
        generated_at = now_iso8601_utc(),
        source       = "flight",
        signature    = sig,
        measured     = Measured(
            gap_ev              = chosen_result["gap_ev"],
            gap_kind            = chosen_result["gap_kind"],
            spin_polarization   = 0.0,        # not measured (7.6)
            total_magnetization = chosen_result.get(
                                    "total_magnetization", 0.0),
            kpoint_density      = converged_density),
        context      = Context(
            # sub-model from THIS structure's record (the sole home;
            #   DESIGN 7.8 step 3f / 6.2.9), never sweep.fixed_axes.
            basis      = prediction["basis"],
            functional = prediction["functional"],
            kpoint_integration = prediction["kpoint_integration"],
            scf_threshold   = chosen_result["scf_threshold"],
            cell_atom_count = structure.num_atoms,
            cell_volume_per_formula_unit =
                structure.real_cell_volume * ANGSTROM3_TO_BOHR3),
        verification = Verification(
            # The distinct-mesh ladder the picker judged, so the
            #   curator's auto_promote_ok re-judges flatness on the
            #   same calculations (DESIGN 7.8 step 3f).
            grid_values   = tuple(grid_values),
            grid_energies = tuple(grid_energies),
            converged_at  = converged_density,
            # The chosen rung's resolved mesh, stored exact beside
            #   the density (DESIGN 3.12.4 / 7.2); read from its
            #   result.toml (6.1.2), absent -> None.  Both harvests
            #   supply chosen_result, so both fill it identically.
            converged_mesh = (tuple(chosen_result["kpoint_mesh"])
                              if "kpoint_mesh" in chosen_result
                              else None),
            metric           = "total_energy",
            metric_threshold = kpoint_threshold,
            predictor_confidence   = prediction["confidence"],
            predictor_neighbor_ids =
                tuple(prediction["neighbor_entry_ids"])),
        provenance   = Provenance(
            flight_id        = flight_id_of(workspace_root),
            source_structure = source_structure,
            # The build behind the run, read from result.toml like
            #   every other per-run fact -- the wingbeat echoed it
            #   there out of the unit's `record` (13.2), so this
            #   harvest stays on its three sources and never opens
            #   the dispatch core's status.toml.  "unknown" remains
            #   the floor for a run that recorded nothing; it is
            #   non-empty, so the schema's rule-11 check passes and
            #   a curator can spot it on review.
            imago_commit     = chosen_result.get("imago_commit")
                               or "unknown",
            curator          = "guidance_harvest.py"))


# How many ladder positions away the gap's flatness is measured.
# TWO, not one, and forced by measurement rather than chosen: a
# k-point ladder carries a strong parity sawtooth in the gap.  On
# diamond silicon adjacent rungs disagree by 19% -- [11,11,11]
# reads 0.9572 eV against [12,12,12]'s 0.8046 -- even where the gap
# has settled to about 1% within one parity family.  Odd and even
# meshes sample the zone differently near the band edges and reach
# the same limit at different rates, so comparing a rung to its
# immediate neighbours calls every ladder unsettled and
# discriminates nothing (DESIGN 7.2).
GAP_SPREAD_STRIDE = 2


function measure_gap_spread(ladder_gaps, chosen_index):
    # How far the gap still moves with the mesh at the chosen rung:
    # the largest RELATIVE change between its gap and the rungs
    # GAP_SPREAD_STRIDE positions either side, as a fraction of the
    # chosen gap.  None when it cannot be measured (DESIGN 7.2).
    #
    # Relative rather than absolute, for a reason the seed ladders
    # supply: near the top of its ladder si_ia-3's gap moves by
    # 0.010-0.014 eV per two rungs, SMALLER in absolute terms than
    # diamond silicon's mid-ladder movement -- yet si_ia-3's gap is
    # collapsing to zero while silicon's has settled.  As fractions
    # they separate cleanly, ~20% against ~1%.
    #
    # None means NOT MEASURED, never "settled": too short a ladder,
    # a missing reading, or a zero gap (a metal, whose relative
    # change is undefined and which stages no entry anyway).
    if ladder_gaps is None or chosen_index is None:
        return None
    if chosen_index outside range(ladder_gaps):
        return None
    chosen_gap = ladder_gaps[chosen_index]
    if chosen_gap is None or chosen_gap <= 0:
        return None
    # Either side alone suffices.  The converged rung often sits
    #   near the top of a ladder whose upper neighbours were never
    #   computed, and one side still answers the question asked:
    #   is the gap still moving?
    spreads = []
    for offset in (-GAP_SPREAD_STRIDE, +GAP_SPREAD_STRIDE):
        j = chosen_index + offset
        if j inside range(ladder_gaps) and ladder_gaps[j] is not None:
            spreads.append(abs(ladder_gaps[j] - chosen_gap)
                           / chosen_gap)
    return max(spreads) if spreads else None


function is_gapless_value(gap_ev, gap_threshold):
    # The metal test, on a bare gap reading (DESIGN 3.12.3 / 7.8).
    # `gap_threshold` is an ABSOLUTE band gap in eV -- not a per-atom
    # energy -- low enough that no real insulator crosses it and high
    # enough to catch a true metal's near-zero reading.
    #
    # A MISSING gap (None) is NOT metallic.  Both callers depend on
    # that side of it, for the same reason from opposite ends: in the
    # climb an absent reading must not stop a search that was
    # converging (4e.2), and in the harvest it must not suppress an
    # entry a genuine insulator earned.  Defaulting the other way
    # would make an unwired gap look like a collection with no
    # insulators in it, which is the failure nobody would question.
    #
    # The scalar core, so the rung-shaped is_gapless (4e.2) and the
    # result-dict-shaped call in build_entry share ONE rule; neither
    # caller has to build a shape it does not have.
    return gap_ev is not None and gap_ev <= gap_threshold


function per_atom_ev(total_energy_hartree, cell_atom_count):
    # A raw total-cell energy (hartree) expressed as eV per atom
    # -- the basis the k-point threshold is stated in (DESIGN 7.8
    # / 5.7, Option B).  HARTREE_TO_EV is the shared hartree->eV
    # constant.  Single-sourced here so pick_converged and
    # auto_promote_ok normalize identically and cannot drift.
    return total_energy_hartree * HARTREE_TO_EV / cell_atom_count


# Two runs of the same resolved mesh are the same calculation and
# must give the same total energy; ENERGY_MATCH_EPS is the gap
# below which they count as equal.  It is tight (near float noise,
# in hartree): a single-threaded imago run is deterministic, so
# genuine duplicates agree to the last digits and a wider gap
# means the runs were not in fact identical.
ENERGY_MATCH_EPS = 1e-9


function collapse_by_mesh(kpoint_densities, energies, meshes):
    # DESIGN 7.8 step 3c guard.  Reduce a density-sorted grid to
    # one rung per distinct resolved mesh, keeping the lowest-
    # density member.  Inputs are parallel arrays ordered by
    # ascending requested density; meshes[i] is rung i's resolved
    # kpoint_mesh (the axial counts, 6.1.2), or None.
    #
    # Returns (collapsed_densities, collapsed_energies, kept),
    # where kept[j] is the ORIGINAL index of the j-th surviving
    # rung, so a caller can map a collapsed index back to the run
    # it names.
    #
    # If any mesh is None -- an older result.toml, or imago not
    # yet emitting kpoint_mesh -- the guard cannot act: return the
    # grid unchanged with identity indices.  The guard is thus
    # INERT until result.toml carries the mesh (DESIGN 6.1.2 /
    # 3.11), and behavior matches the pre-guard code.
    if any(m is None for m in meshes):
        return kpoint_densities, energies, list(range(len(meshes)))

    kept = []
    for i in range(len(meshes)):
        # Same mesh as the last surviving rung?  Then rung i is
        # that calculation run again (the map is monotone in
        # density, 3.7, so equal meshes are contiguous).  An equal
        # mesh MUST give an equal energy; a mismatch means the
        # runs were not identical and is surfaced, not averaged.
        if kept and meshes[i] == meshes[kept[-1]]:
            if abs(energies[i] - energies[kept[-1]]) > ENERGY_MATCH_EPS:
                error("k-density rungs "
                      + kpoint_densities[kept[-1]] + " and "
                      + kpoint_densities[i]
                      + " resolved to the same mesh "
                      + str(meshes[i]) + " but disagree in total"
                      + " energy -- not the same calculation")
            continue                        # drop the duplicate
        kept.append(i)

    return ([kpoint_densities[i] for i in kept],
            [energies[i] for i in kept],
            kept)


function pick_converged(energies, cell_atom_count, threshold):
    # DESIGN 7.8 step 3c: the smallest interior grid index i
    # at which BOTH consecutive-pair energy deltas fall below
    # `threshold`.  Two-sided so a single-point fluke does not
    # masquerade as convergence.  Returns None if the energy
    # is still moving (no flat interior point).
    #
    # `energies` are raw total-cell hartree (Option B) and
    # `threshold` is per atom in eV, so normalize the whole
    # ladder to that basis once (per_atom_ev) before comparing.
    # The per-atom scale keeps a big cell from being held to a
    # tighter bound than a small one (DESIGN 7.8).
    e = [per_atom_ev(x, cell_atom_count) for x in energies]
    for i in range(1, len(e) - 1):
        below_up   = abs(e[i] - e[i + 1]) < threshold
        below_down = abs(e[i] - e[i - 1]) < threshold
        if below_up and below_down:
            return i
    return None
```

**`guidance_promote.py`** is the curator helper.  Four
modes; promotion is a `mv` of the file from
`staging/<system_type>/` to `entries/<system_type>/`, so the
contents (and provenance) never change.

```
function dedup_key(entry):
    # The identity of the CLAIM an entry makes (DESIGN 7.8):
    # same system, same settings, same structure.  The three
    # settings fields are the predictor's sub-model partition,
    # so a gaussian and a gaussian-0.1 run of one solid are
    # NOT re-runs of each other.  The structure's basename, not
    # its path, because the path records only where the
    # structure cache sat and that has moved (ARCH 8.1).
    # imago_commit is deliberately absent: it is what the
    # comparison examines, not what the key partitions on.
    return (entry.signature.system_type,
            entry.context.basis,
            entry.context.functional,
            entry.context.kpoint_integration,
            basename(entry.provenance.source_structure))


function mesh_of(entry):
    # The entry's converged mesh, or None when it has no
    # verification block at all (a manual entry, 7.9).  REPORTING
    # ONLY: no branch tests this value.  The occupied case acts
    # the same whether the meshes agree or not (DESIGN 7.8), so
    # the mesh is shown to the person and nothing else.
    if entry.verification is None:
        return None
    return entry.verification.converged_mesh
```

```
function promote(db_root, mode):
    # Load the promoted corpus ONCE, already keyed by claim, and
    # keep each entry's PATH beside it so REPLACE can retire the
    # file it names.  This is the one judgment the acceptance rule
    # cannot make from the staged file alone (DESIGN 7.8);
    # entries/ is small and local, so the cost is a directory
    # read, not a workspace re-read.  staging/ and superseded/ are
    # NOT loaded: only a promoted entry can hold a claim.
    #
    # The index is LIVE, not a snapshot: every branch below that
    # fills a slot updates it, which is what makes a within-batch
    # duplicate fall out of the ordinary rule (DESIGN 7.8).
    #
    # ONE shape throughout: key -> (path, entry), where `path` is
    # where that entry's file lives NOW, under entries/.  Every
    # slot-filling branch stores the pair, never a bare entry,
    # because REPLACE has to retire the file the promoted entry
    # occupies and therefore needs its path.  Store bare entries and
    # a SECOND replace against the same claim has nothing to retire.
    # load_promoted_entries walks entries/*/*.toml, so the paths cost
    # nothing to carry; deriving them from the entry_id instead would
    # work today only because save_entry happens to name each file
    # <entry_id>.toml, and that is an invariant not worth leaning on.
    promoted = load_promoted_entries(db_root)   # key -> (path,
                                                #         entry)

    # dry-run evaluates every decision below and moves nothing,
    # so a curator sees the whole outcome -- promotions,
    # retirements, and conflicts -- before a file is touched.
    dry = (mode == "dry-run")

    for system_type in VALID_SYSTEM_TYPES:
        # Sorted for determinism only.  There is NO separate
        # batch-resolution pass and no generated_at tie-break: the
        # index below is updated as entries are promoted, so a
        # second staged file making a claim the first just filled
        # simply finds it occupied and takes the ordinary branch
        # (DESIGN 7.8).  One rule, applied uniformly, is why
        # staging IS a uniqueness namespace without extra
        # machinery to make it one.
        staged = sorted(glob(
            join(db_root, "staging", system_type), "*.toml"))

        for path in staged:
            entry = load_entry(path, system_type, {})
            key = dedup_key(entry)
            prior = promoted.get(key)

            if prior is not None:
                prior_path, prior_entry = prior
                # OCCUPIED.  An existence test, not a comparison:
                #   the collection holds one record per structure
                #   per settings, and this claim is already held.
                #   Report both -- meshes, dates, builds -- and
                #   retire the newcomer.  The promoted entry is
                #   left byte-identical; nothing here can retract
                #   it.  The meshes are PRINTED, never tested: the
                #   action is the same whether they agree or not,
                #   and a disagreement is a person's question
                #   (DESIGN 7.8, VISION 16).
                report_occupied(entry, prior_entry)
                                        # meshes + dates + builds

                # Interactive review is the ONLY mode offering
                #   REPLACE, and REPLACE is the only route by
                #   which anything leaves entries/.  Unattended
                #   modes never retract.
                if mode == "interactive":
                    choice = ask("REPLACE / SKIP / DELETE")
                    if choice == "REPLACE":
                        # The promoted entry's own file is retired,
                        #   which is why its path had to be carried.
                        #   `destination` is where the newcomer now
                        #   lives, so the index keeps pointing at a
                        #   real file and a later REPLACE against
                        #   this same claim has something to retire.
                        destination = entries_path(
                            db_root, system_type, path)
                        if not dry:
                            move_to_superseded(
                                prior_path, db_root, system_type)
                            destination = move_to_entries(
                                path, db_root, system_type)
                        promoted[key] = (destination, entry)
                        record(entry, "replaced",
                               "replaced " + prior_entry.entry_id)
                        continue
                    if choice == "DELETE":
                        if not dry:
                            remove_file(path)
                        record(entry, "deleted", "")
                        continue
                    record(entry, "skipped", "left in staging")
                    continue

                if not dry:
                    move_to_superseded(path, db_root, system_type)
                record(entry, "superseded",
                       "already promoted as " + prior_entry.entry_id)
                continue

            # FREE: the ordinary path, decided by mode.  `all`
            # still passed the existence test above -- it waives
            # the quality rule, not the correctness guard.
            # Each branch that fills the slot stores (path, entry),
            # the one shape the index has; move_to_entries returns
            # the destination it renamed to, so the path recorded is
            # the file that now exists rather than a guess at it.
            if mode == "dry-run":
                print_would_promote(entry,
                    auto_promote_ok(entry))
                # A dry run must model the index too, or a second
                # staged file claiming this slot would be reported
                # as free when a real run would find it taken.  It
                # moves nothing, so it records where the file WOULD
                # land.
                if auto_promote_ok(entry):
                    promoted[key] = (entries_path(db_root,
                                         system_type, path), entry)
            elif mode == "all":
                promoted[key] = (
                    move_to_entries(path, db_root, system_type),
                    entry)
            elif mode == "auto-promote":
                if auto_promote_ok(entry):
                    promoted[key] = (
                        move_to_entries(path, db_root, system_type),
                        entry)
                # else: leave in staging for review.
            else:                                # interactive
                print_summary(entry)             # sig+measured
                                                 #   +verif+prov
                choice = ask("PROMOTE / SKIP / DELETE")
                if choice == "PROMOTE":
                    promoted[key] = (
                        move_to_entries(path, db_root, system_type),
                        entry)
                elif choice == "DELETE":
                    remove_file(path)
                # SKIP: leave in staging.
```

```
function load_promoted_entries(db_root):
    # The promoted corpus, keyed by claim: key -> (path, entry).  It
    # yields PAIRS, not bare entries, because REPLACE retires the
    # file a promoted entry occupies and so needs its path (above).
    # Only entries/ is read -- staging/ and superseded/ hold no
    # claims -- and it is a directory read of small local files, not
    # the flight-workspace re-read the single-file acceptance rule
    # exists to avoid (DESIGN 7.8).
    #
    # On two promoted files sharing one claim the later read wins,
    # silently.  That is a pre-existing state this pass cannot fix
    # (it has no verb for removing a promoted entry unasked), and
    # refusing to run would leave the curator unable to promote
    # anything else either.
    promoted = {}
    for system_type in VALID_SYSTEM_TYPES:
        for path in sorted(glob(
                join(db_root, "entries", system_type), "*.toml")):
            entry = load_entry(path, system_type, {})
            promoted[dedup_key(entry)] = (path, entry)
    return promoted


function entries_path(db_root, system_type, staged_path):
    # Where a staged file lands once promoted.  Promotion is a pure
    # rename that keeps the basename, so this is the one formula for
    # that destination; move_to_entries computes the same thing and
    # returns it, and dry-run -- which renames nothing but must still
    # model the index (above) -- asks here instead.
    return join(db_root, "entries", system_type,
                basename(staged_path))


function move_to_superseded(path, db_root, system_type):
    # Retire a staged entry that a promoted one already claims.
    # Mirrors move_to_entries: a pure rename into
    # superseded/<system_type>/, refusing a pre-existing
    # destination rather than overwriting.  Retired and not
    # deleted so the record of what a re-run produced stays
    # recoverable, and so staging/ does not accrete files that
    # every later promotion pass re-examines (ARCH 10.1).
    destination = join(db_root, "superseded", system_type,
                       basename(path))
    if exists(destination):
        raise "superseded collision; resolve by hand"
    rename(path, destination)
    return destination
```

`promote` returns one `(entry_id, action)` record per staged
file, so a caller or a test sees the outcome without parsing
printed text.  The uniqueness rule adds two actions to the
existing `promoted` / `skipped` / `deleted` (and the `would-`
forms): `superseded` for a retired re-run, and `replaced` for
the curator-driven swap that is the only way an entry ever
leaves `entries/`.

```
function auto_promote_ok(entry):
    # DESIGN 7.8 objective acceptance test, evaluated from
    # the staging file alone (this is why harvest records
    # grid_energies).
    v = entry.verification
    if v is None or v.grid_energies is None:
        return False                      # manual / no sweep

    # 1. Converged density in the middle 60% of the grid
    #    (not at either endpoint -- a grid that may have been
    #    too narrow).
    lo, hi = v.grid_values[0], v.grid_values[-1]
    if hi == lo:
        return False
    position = (v.converged_at - lo) / (hi - lo)
    if not (0.2 <= position <= 0.8):
        return False

    # 2. Top-three grid points convincingly flat: their SPREAD
    #    (max - min), per atom in eV, below metric_threshold * 10.
    #    grid_energies are raw total-cell hartree (Option B), so
    #    normalize each via the same per_atom_ev helper
    #    pick_converged uses; cell_atom_count rides on the entry's
    #    context.  A spread is a like-for-like linear quantity; a
    #    variance would be an energy squared against a linear
    #    threshold, so it is deliberately not used (DESIGN 7.8).
    n = entry.context.cell_atom_count
    top3 = [per_atom_ev(x, n) for x in v.grid_energies[-3:]]
    if (max(top3) - min(top3)) >= v.metric_threshold * 10.0:
        return False

    # 3. gap_ev / gap_kind consistent.
    is_metal = (entry.measured.gap_ev == 0.0)
    if (entry.measured.gap_kind == "none") != is_metal:
        return False

    return True
```

In practice the rule auto-promotes ~80% of a seed flight
(TODO C75) and leaves the ~20% endpoint-converged or
not-yet-flat outliers for the curator's interactive review
-- the friction that makes a 250-entry seed tractable
without rubber-stamping every entry (DESIGN 7.8).

## 16. Resource & Cost Guidance Dataspace (DESIGN 8)

The cost prong, and a deliberate *sibling* of section 15
rather than an extension of it (ARCHITECTURE 11.1).  Where
the historical-guidance dataspace records what operating
point is **accurate**, this one records what a run **costs**:
the problem-size signature, the parallel execution
configuration, the build the binary was compiled with, and
the measured peak memory, disk, and walltime.  A
physics-informed regressor learns the cost surface, and the
near-term consumer turns a prediction into a scheduler
request that neither overflows memory nor exceeds the
walltime limit (DESIGN 8.1).

The two stay separate because a converged k-density
transfers across machines and a walltime does not.  So this
dataspace is partitioned by a **hardware fingerprint**, and
its atomic unit is one **execution observation** -- a single
run under a single configuration, never collapsed to a
per-system summary.  That is what lets the same artifact
serve provisioning now and configuration optimization, build
comparison, and scaling studies later, with no schema change.

Blocks, helpers first then drivers: the library
`resource_db.py` (16.1 shapes and registries, 16.2 hardware
fingerprint, 16.3 reader, 16.4 emitter, 16.5 predictor); the
dispatch-time capture hooks (16.6); the producers
`resource_harvest.py` / `resource_promote.py` (16.7); and the
provisioning consumer in the flight layer (16.8).  All Python
under `src/scripts/`.

**Read 16.9 before implementing 16.5.**  DESIGN 8.9 leaves
four questions open, and two of them reach into the
predictor.  This section specifies everything the design
pins and stops -- visibly -- where it does not.  A pseudocode
that guessed past those points would read as settled and
would be nothing of the kind.

### 16.1 Constants, registries, and shapes (DESIGN 8.4)

The constants and the three registries live in one place so a
post-seed recalibration or a new knob is a one-file change.
A key absent from its registry is rejected at load (rule 11):
extensibility goes through the registry, never through silent
key drift.

```
SCHEMA_VERSION          = 1
VALID_SOURCES           = ("flight", "manual")
VALID_OUTCOMES          = ("completed", "oom", "timeout",
                           "failed")
WAVEFUNCTION_COMPONENTS = (1, 4)     # Schrodinger | Dirac
SPIN_CHANNELS           = (1, 2)

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

# The execution knobs that must be present and > 0 (rule 7).
# `binding` is validated against VALID_BINDINGS instead.
REQUIRED_POSITIVE_EXECUTION = ("node_count", "cores_per_node",
    "total_cores", "mpi_ranks", "omp_threads_per_rank")

# The metrics a completed run must carry, all > 0 (rule 9).
REQUIRED_COMPLETED_METRICS = ("peak_memory_bytes",
    "disk_bytes", "walltime_seconds")

# Calibration knobs.  DESIGN 8.6 states plainly that the
# thin-group threshold and the safety margin are tuned AFTER
# the seed flight (TODO C82) and are deliberately not required
# for the artifact to begin accumulating data.  They are named
# here, in one place, so calibration is a one-line edit -- this
# document does not invent their values.  Section 15 fixes its
# analogous knobs because its seed exists; ours does not yet.
MIN_GROUP_SAMPLES = <calibrated by the seed; 8.6/8.9>
SAFETY_MARGIN     = <calibrated by the seed; 8.6>

# k-NN fallback knobs, used only when a group is too thin to
# fit (16.5).  Same rule: named here, valued by the seed.
neighbor_count        = <calibrated by the seed; 8.6>
distance_floor        = 1e-6   # keeps 1/distance finite when a
                               #   neighbour coincides with the
                               #   query point
secular_dim_weight    = 1.0    # the dominant cost axis, so the
                               #   other two are scaled against
                               #   it and it holds the unit
kpoint_count_weight   = <calibrated by the seed; 8.6>
spin_channel_weight   = <calibrated by the seed; 8.6>
```

The dataclasses mirror the schema block for block.
`ExecutionConfig` and `BuildConfig` hold open dictionaries
rather than fixed fields *precisely so* the registries -- not
the dataclass definitions -- are the single source of truth
for which knobs exist.  Promoting a studied compiler flag to a
first-class feature is then a registry edit (ARCHITECTURE
11.3).

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
    knobs : dict          # registry-validated key -> value

dataclass BuildConfig:
    knobs : dict          # registry-validated COARSE knobs;
                          #   the verbatim compile string is
                          #   Provenance.compile_string

dataclass MeasuredResources:
    metrics  : dict       # registry-validated metric -> value
    censored : bool       # True when outcome != completed:
                          #   a BOUND, not a point measurement

dataclass Provenance:
    flight_id        : str
    source_structure : str
    imago_commit     : str
    hostname         : str
    compile_string   : str   # build fidelity layer, verbatim
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
    schema_version              : int
    observations_by_fingerprint : dict   # fp -> list[Obs]
    hardware_registry           : dict   # fp -> attributes
```

### 16.2 Hardware fingerprint and registry (DESIGN 8.5)

The fingerprint is the coarse partition *within which cost is
comparable*.  Its whole job is to be stable: the CPU string is
normalized to vendor plus microarchitecture family, and the
stepping, base clock, and exact model number are dropped, so
routine BIOS or microcode churn does not fragment the data.

```
function hardware_fingerprint(attrs):
    # v1 recipe (DESIGN 8.5):
    #   <cpu_vendor>-<cpu_microarch>-<cores>c-<mem_gb>gb
    # e.g. "intel-haswell-24c-128gb".  Normalization is what
    # makes the slug stable, so it is done here and nowhere
    # else: lowercase, spaces to hyphens, and NOTHING of the
    # stepping / clock / model number survives.
    vendor    = normalize_token(attrs.cpu_vendor)
    microarch = normalize_token(attrs.cpu_microarch)
    return (vendor + "-" + microarch + "-"
            + str(attrs.cores_per_node) + "c-"
            + str(attrs.memory_per_node_gb) + "gb")
```

`hardware_registry.toml` maps each fingerprint to its full
probed attributes (exact CPU model, socket count, memory,
interconnect) for diagnostics.  The observation files carry
only the fingerprint and never the repeated attributes,
mirroring how section 15 keeps the element group table out of
individual entries.

```
function load_hardware_registry(path):
    raw = tomllib.load(path)
    registry = {}
    for fingerprint, attrs in raw.items():
        # The slug must be REDERIVABLE from the attributes it
        # maps to.  A registry whose key disagrees with its own
        # body would silently partition two machines together.
        require(hardware_fingerprint(attrs) == fingerprint,
            path, fingerprint,
            "fingerprint disagrees with its probed attributes")
        registry[fingerprint] = attrs
    return registry
```

The three registry validators are one function, because the
three registries differ only in their name and contents.  The
error names the file, the block, and the offending key -- the
same discipline as section 15 and DESIGN 5.2.

```
function require_registered(table, registry, path, block):
    # Rule 11.  An unknown key is a HARD error: it is nearly
    # always a typo or a knob someone added without extending
    # the registry, and either way silently dropping it loses
    # the very fact the observation exists to record.
    for key in table:
        require(key in registry, path, block,
            "unknown key " + key + " (extend the registry and "
            "bump SCHEMA_VERSION to add a knob)")
```

### 16.3 TOML reader `load()` (DESIGN 8.2, rules 1-12)

`load(root)` walks `entries/<fingerprint>/` and returns a
validated `ResourceDataspace`.  Staging is read the same way
by the promoter (16.7) and is never mixed into a prediction
pool.

```
function load(root):
    marker = read_file(join(root, "SCHEMA_VERSION")).strip()
    require(marker == str(SCHEMA_VERSION),               # 1
        join(root, "SCHEMA_VERSION"),
        "marker " + marker + " != " + str(SCHEMA_VERSION))

    registry = load_hardware_registry(
        join(root, "hardware_registry.toml"))

    by_fingerprint = {}
    seen_ids = {}                    # observation_id -> path
    for fingerprint in sorted(registry):
        subdir = join(root, "entries", fingerprint)
        if not exists(subdir):
            continue                 # a registered, unseeded
                                     #   machine: legal, empty
        by_fingerprint[fingerprint] = [
            load_observation(path, fingerprint, seen_ids,
                             registry)
            for path in sorted(glob(subdir, "*.toml"))]

    return ResourceDataspace(
        schema_version              = SCHEMA_VERSION,
        observations_by_fingerprint = by_fingerprint,
        hardware_registry           = registry)
```

`load_observation` checks the schema BEFORE building the
dataclass (rule 12), so an omission surfaces as a validation
failure naming the field rather than a bare constructor
error.

```
function load_observation(path, fingerprint_dir, seen_ids,
                          registry):
    raw = tomllib.load(path)

    for f in ("schema_version", "observation_id",
              "generated_at", "source", "outcome",
              "hardware_fingerprint"):
        require(f in raw, path, "top-level", "missing: " + f)

    require(raw["schema_version"] == SCHEMA_VERSION,     # 1
        path, "top-level", "schema_version must be "
        + str(SCHEMA_VERSION))

    obs_id = raw["observation_id"]                       # 2
    require(obs_id == file_stem(path), path, "top-level",
        "observation_id must equal the file stem")
    require(obs_id not in seen_ids, path, "top-level",
        "duplicate observation_id, also in "
        + seen_ids.get(obs_id, "?"))
    seen_ids[obs_id] = path

    fingerprint = raw["hardware_fingerprint"]            # 3
    require(fingerprint in registry, path, "top-level",
        "unregistered hardware_fingerprint " + fingerprint)
    require(fingerprint == fingerprint_dir, path, "top-level",
        "hardware_fingerprint disagrees with its directory")

    require(raw["source"] in VALID_SOURCES,              # 4
        path, "top-level", "source must be one of "
        + str(VALID_SOURCES))
    outcome = raw["outcome"]                             # 5
    require(outcome in VALID_OUTCOMES, path, "top-level",
        "outcome must be one of " + str(VALID_OUTCOMES))

    signature  = load_signature(raw, path)               # 6
    execution  = load_execution(raw, path)               # 7
    build      = load_build(raw, path)                   # 8
    resources  = load_resources(raw, path, outcome)      # 9
    provenance = load_provenance(raw, path, raw["source"])  # 10

    return Observation(
        observation_id       = obs_id,
        generated_at         = raw["generated_at"],
        source               = raw["source"],
        outcome              = outcome,
        hardware_fingerprint = fingerprint,
        signature            = signature,
        execution            = execution,
        build                = build,
        resources            = resources,
        provenance           = provenance)
```

```
function load_signature(raw, path):
    # Rule 6.  Every count is a physical quantity, so every
    # bound below is a statement about physics, not taste:
    # a cell has atoms, an all-electron run has electrons, the
    # valence space is a subset of them, and the spinor
    # structure is 1-component or 4-component, nothing between.
    block = require_block(raw, "observation.signature", path)
    for f in ("atom_count", "electron_count",
              "basis_function_count", "secular_dimension",
              "kpoint_count"):
        require(block[f] > 0, path, "signature", f + " must be > 0")
    require(0 <= block["valence_electron_count"]
              <= block["electron_count"],
        path, "signature",
        "valence_electron_count must lie in [0, electron_count]")
    require(block["wavefunction_components"]
              in WAVEFUNCTION_COMPONENTS,
        path, "signature",
        "wavefunction_components must be 1 or 4")
    require(block["spin_channels"] in SPIN_CHANNELS,
        path, "signature", "spin_channels must be 1 or 2")
    return SizeSignature(**block)
```

```
function load_execution(raw, path):
    # Rule 7.  total_cores is recorded, never derived from
    # node_count x cores_per_node, so a partially-packed node
    # stays faithful -- which is exactly the case a derived
    # value would quietly misreport.
    block = require_block(raw, "observation.execution", path)
    require_registered(block, EXECUTION_KNOB_REGISTRY,
                       path, "execution")
    for f in REQUIRED_POSITIVE_EXECUTION:
        require(f in block and block[f] > 0, path, "execution",
            f + " must be present and > 0")
    require(block["binding"] in VALID_BINDINGS, path,
        "execution", "binding must be one of "
        + str(VALID_BINDINGS))
    return ExecutionConfig(knobs = block)
```

```
function load_build(raw, path):
    # Rule 8.  The COARSE layer only: bucketed values that make
    # a build a comparable feature rather than a fragmenting
    # one (an optimization LEVEL, not a flag string; a MAJOR
    # version, not a patch).  The verbatim compile string is
    # provenance, below.
    block = require_block(raw, "observation.build", path)
    require_registered(block, BUILD_KNOB_REGISTRY, path, "build")
    for f in ("compiler_family", "optimization_level"):
        require(f in block and block[f] != "", path, "build",
            f + " must be present and non-empty")
    return BuildConfig(knobs = block)
```

```
function load_resources(raw, path, outcome):
    # Rule 9, and the censoring rule of DESIGN 8.7.  A
    # non-completed run is NOT discarded: an OOM kill is
    # positive evidence that this configuration is insufficient
    # at this size, and a timeout bounds the walltime from
    # above.  Marking `censored` here is what stops the
    # predictor (16.5) from ever reading a bound as a point.
    block = require_block(raw, "observation.resources", path)
    require_registered(block, RESOURCE_METRIC_REGISTRY,
                       path, "resources")

    if outcome == "completed":
        for f in REQUIRED_COMPLETED_METRICS:
            require(f in block and block[f] > 0, path,
                "resources", f + " must be present and > 0 "
                "for outcome=completed")
    else if outcome == "oom":
        require("peak_memory_bytes" in block, path, "resources",
            "outcome=oom must record the memory limit it hit "
            "(a LOWER bound on the true need)")
    else if outcome == "timeout":
        require("walltime_seconds" in block, path, "resources",
            "outcome=timeout must record the walltime limit "
            "(an UPPER bound on the true need)")
    # outcome=failed carries no usable cost signal at all; it
    #   is kept for diagnostics and never promoted (16.7).

    return MeasuredResources(
        metrics  = block,
        censored = (outcome != "completed"))
```

```
function load_provenance(raw, path, source):
    # Rule 10.  compile_string is required for EVERY source --
    # it is the build fidelity layer, and recording it verbatim
    # is what lets a flag that is not a coarse knob today be
    # recovered post-hoc when it turns out to matter.
    block = require_block(raw, "observation.provenance", path)
    require(block["compile_string"] != "", path, "provenance",
        "compile_string must be present for every source")
    if source == "flight":
        for f in ("flight_id", "source_structure",
                  "imago_commit"):
            require(block.get(f, "") != "", path, "provenance",
                f + " must be non-empty for source=flight")
    return Provenance(**block)
```

### 16.4 Hand-formatted emitter `save_observation()` (8.3)

The same deterministic discipline as 15.4 and DESIGN 7.5:
fixed block sequence, fixed key order within each block,
`%.16e` for every real, byte-identical output for a given
in-memory observation.  The observation_id doubles as the file
stem (rule 2), and its shape -- fingerprint plus a short hash
-- makes a file self-identifying on sight.

```
function observation_slug(obs):
    # DESIGN 8.3: "<fingerprint>-<6 hex>".  The hash is taken
    # over the provenance that distinguishes two runs which
    # agree on everything else, so two harvests of the same
    # flight cannot collide unless they are the same run.
    blob = (obs.provenance.flight_id
            + obs.provenance.source_structure
            + obs.generated_at)
    return obs.hardware_fingerprint + "-" + sha256_hex(blob)[:6]

function save_observation(obs, root, area = "staging"):
    # Harvest writes to staging/; the curator's promote (16.7)
    # is what moves a file to entries/.  Refuse a collision
    # rather than overwrite: rule 2 makes the id unique across
    # BOTH trees, so an existing path is a real conflict.
    subdir = join(root, area, obs.hardware_fingerprint)
    make_dirs(subdir)
    path = join(subdir, obs.observation_id + ".toml")
    require(not exists(path),
        "save_observation: " + path + " already exists")
    write_file(path, format_observation(obs))
    return path
```

```
function format_observation(obs):
    # Block order is FIXED, and matches DESIGN 8.3's sketch, so
    # a diff between two observations is a diff of values and
    # never of layout.  Optional metrics are emitted only when
    # present; phase_timings is a sub-table and therefore last,
    # as TOML requires of a table's own keys before sub-tables.
    lines = [emit_scalar("schema_version", SCHEMA_VERSION),
             emit_scalar("observation_id", obs.observation_id),
             emit_scalar("generated_at", obs.generated_at),
             emit_scalar("source", obs.source),
             emit_scalar("outcome", obs.outcome),
             emit_scalar("hardware_fingerprint",
                         obs.hardware_fingerprint)]

    emit_block(lines, "[observation.signature]",
               obs.signature, SIGNATURE_FIELD_ORDER)
    emit_block(lines, "[observation.execution]",
               obs.execution.knobs, EXECUTION_KNOB_REGISTRY)
    emit_block(lines, "[observation.build]",
               obs.build.knobs, BUILD_KNOB_REGISTRY)

    # Resources, minus phase_timings, then phase_timings.
    emit_block(lines, "[observation.resources]",
               obs.resources.metrics,
               [m for m in RESOURCE_METRIC_REGISTRY
                if m != "phase_timings"])
    if "phase_timings" in obs.resources.metrics:
        emit_block(lines, "[observation.resources.phase_timings]",
                   obs.resources.metrics["phase_timings"],
                   PHASE_ORDER)   # setup, scf, eigensolve, ...

    emit_block(lines, "[observation.provenance]",
               obs.provenance, PROVENANCE_FIELD_ORDER)
    return join(lines, "\n") + "\n"
```

Round-trip is the emitter's contract, exactly as it is the
manifest writer's (11.6): whatever `load()` reads,
`format_observation` writes back.  Every block the reader
parses therefore has a writing counterpart above, and a
metric the reader accepts but the writer drops would silently
delete a measurement on the next curation pass.

### 16.5 Predictor `predict()` (DESIGN 8.6)

Within a fixed `(hardware_fingerprint, build_bucket)` group,
cost is a smooth, physics-grounded function of size and
parallel configuration.  That is why this predictor is a
**physics-informed regression** and not the pure k-NN of 15.5:
we know the shape of the curve before we see the data.  Peak
memory scales roughly as the square of `secular_dimension` and
the eigensolve as its cube, so the fit is a power law in log
space and the exponent is *recovered from the data* rather
than assumed:

```
log(resource) = log(A) + p * log(secular_dimension)
                + (parallel and spin correction terms)
```

Expected `p` is near 2 for memory and near 3 for walltime.
Recovering it rather than fixing it is the point: an exponent
that comes back far from its expectation is telling you
something true about the code, and a fixed exponent could
never say it.

```
function build_bucket(build):
    # The group key's second half.  A build is a COMPARABLE
    # feature because its knobs are bucketed (8.2); the tuple
    # is taken in registry order so the key is deterministic.
    return tuple(build.knobs[k] for k in BUILD_KNOB_REGISTRY)

function group_key(obs):
    return (obs.hardware_fingerprint, build_bucket(obs.build))
```

```
function feature_row(signature, execution):
    # The design matrix row.  DESIGN 8.6 names the quantities
    # the correction terms capture -- the speedup from
    # `mpi_ranks` and `omp_threads_per_rank`, the memory split
    # across ranks, and the spin channels -- but says plainly
    # that "the exact functional form ... [is a] tuning knob
    # calibrated after the seed flight."
    #
    # So this is the v1 form, in ONE place, and it is expected
    # to be re-derived once the seed exists (16.9, TODO C82).
    # It is a named function, not an expression inlined into
    # the fit, precisely so recalibration touches one thing.
    return [1.0,
            log(signature.secular_dimension),   # the exponent p
            log(signature.kpoint_count),
            log(signature.spin_channels),
            log(execution.knobs["mpi_ranks"]),
            log(execution.knobs["omp_threads_per_rank"])]
```

```
function fit_group(observations, metric):
    # Least squares in log space over ONE (fingerprint,
    # build_bucket) group.  Returns coefficients, or None when
    # the group is too thin to fit a stable exponent -- in
    # which case predict() falls back to k-NN, below.
    #
    # !! CENSORED OBSERVATIONS DO NOT ENTER HERE. !!  DESIGN 8.7
    # is explicit that an OOM memory figure is a LOWER bound and
    # a timeout walltime an UPPER bound, never a point; and
    # DESIGN 8.9 records, as an OPEN QUESTION, how such bounds
    # should enter the least-squares fit (a censored/Tobit-style
    # regression, or a weighting scheme).  That question is not
    # settled, so this pseudocode does not settle it: see 16.9.
    # Feeding a bound in as though it were a measurement would
    # bias every exponent the artifact ever reports.
    usable = [o for o in observations
              if not o.resources.censored
              and metric in o.resources.metrics]
    if len(usable) < MIN_GROUP_SAMPLES:
        return None

    rows    = [feature_row(o.signature, o.execution)
               for o in usable]
    targets = [log(o.resources.metrics[metric]) for o in usable]
    return least_squares(rows, targets)
```

```
function predict_metric(coefficients, signature, execution):
    row = feature_row(signature, execution)
    return exp(dot(coefficients, row))
```

The fallback chain has two rungs, and DESIGN 8.5 constrains
both with one sentence: the predictor "never silently predicts
from one machine for another."  So every result names the
fingerprint it was fitted on, and a cold start says so.

```
dataclass ResourcePrediction:
    peak_memory_bytes : real | None
    disk_bytes        : real | None
    walltime_seconds  : real | None
    source_fingerprint : str   # the group actually fitted --
                               #   may differ from the query's
    is_borrowed        : bool  # fitted on a NEIGHBOURING machine
    is_cold_start      : bool  # no usable data anywhere
    sample_count       : int
```

```
function predict(dataspace, query_fingerprint, build, signature,
                 execution):
    # Rung 1: this machine, this build bucket.
    pool = dataspace.observations_by_fingerprint.get(
        query_fingerprint, [])
    key  = (query_fingerprint, build_bucket(build))
    group = [o for o in pool if group_key(o) == key]

    fitted = {m: fit_group(group, m)
              for m in ("peak_memory_bytes", "disk_bytes",
                        "walltime_seconds")}
    if all(c is not None for c in fitted.values()):
        return prediction_from(fitted, signature, execution,
            source_fingerprint = query_fingerprint,
            is_borrowed = False, is_cold_start = False,
            sample_count = len(group))

    # Rung 2: a thin group.  k-NN over the SAME machine, on the
    # size and parallel axes, before ever leaving the machine --
    # a same-machine neighbour is better evidence than a fitted
    # curve from a different one.
    if len(pool) > 0:
        return knn_predict(pool, signature, execution,
            source_fingerprint = query_fingerprint)

    # Rung 3: this machine has NO observations.  Borrow from the
    # nearest registered relative, and SAY SO -- the consumer
    # (16.8) widens its margin on a borrowed prediction.
    neighbour = nearest_fingerprint(query_fingerprint,
                                    dataspace.hardware_registry)
    if neighbour is not None and \
            len(dataspace.observations_by_fingerprint
                .get(neighbour, [])) > 0:
        result = predict(dataspace, neighbour, build, signature,
                         execution)
        return with_flags(result, source_fingerprint = neighbour,
                          is_borrowed = True)

    # Rung 4: nothing anywhere.  Day-1 on a fresh cluster.
    return ResourcePrediction(None, None, None,
        source_fingerprint = query_fingerprint,
        is_borrowed = False, is_cold_start = True,
        sample_count = 0)
```

```
function knn_predict(pool, signature, execution,
                     source_fingerprint):
    # DESIGN 8.6's thin-group fallback: k-NN over
    # secular_dimension, kpoint_count, and spin_channels,
    # scaled by the parallel config.
    #
    # Distances are taken in LOG space on the two count axes,
    # because cost spans orders of magnitude: a 4000-dimension
    # secular problem and a 400-dimension one are a factor of
    # ten apart, not "3600 apart," and a linear distance would
    # let every large system swamp the neighbourhood of every
    # small one.  Censored observations are excluded for the
    # same reason they are excluded from the fit (16.5): a
    # bound is not a measurement.
    usable = [o for o in pool if not o.resources.censored]
    if len(usable) == 0:
        return cold_start_prediction(source_fingerprint)

    function distance_to(observation):
        other = observation.signature
        size_gap = log(signature.secular_dimension
                       / other.secular_dimension)
        kpt_gap  = log(signature.kpoint_count
                       / other.kpoint_count)
        # Spin is a 1-or-2 label, not a magnitude, so it enters
        # as a mismatch penalty rather than a ratio.
        spin_gap = (0.0 if signature.spin_channels
                            == other.spin_channels else 1.0)
        return sqrt(secular_dim_weight  * size_gap ** 2
                    + kpoint_count_weight  * kpt_gap ** 2
                    + spin_channel_weight  * spin_gap ** 2)

    # The parallel config scales the ANSWER, not the distance:
    # two runs of the same problem on different core counts are
    # equally near in problem space, and it is their costs that
    # differ.  So each neighbour's cost is rescaled to the
    # query's configuration before it is averaged in.
    neighbors = nearest(usable, distance_to, neighbor_count)
    weights   = [1.0 / max(distance_to(o), distance_floor)
                 for o in neighbors]

    return weighted_average_of_rescaled_costs(
        neighbors, weights, signature, execution,
        source_fingerprint)
```

```
function nearest_fingerprint(fingerprint, registry):
    # "The nearest related fingerprint by probed attributes"
    # (DESIGN 8.5).  The concrete ordering is this pseudocode's
    # reading of that phrase, and it is deliberately narrow:
    # only a machine of the SAME vendor and microarchitecture is
    # a candidate at all, because a walltime does not carry
    # across instruction sets.  Among those, the closest by core
    # count and memory, in log space.
    attrs = registry[fingerprint]
    candidates = [f for f, a in registry.items()
                  if f != fingerprint
                  and a.cpu_vendor    == attrs.cpu_vendor
                  and a.cpu_microarch == attrs.cpu_microarch]
    if len(candidates) == 0:
        return None
    return min(candidates, key = lambda f: (
        abs(log(registry[f].cores_per_node
                / attrs.cores_per_node))
        + abs(log(registry[f].memory_per_node_gb
                  / attrs.memory_per_node_gb))))
```

### 16.6 Capture: what the run records (DESIGN 8.7)

An observation is assembled at harvest from the four sources
of ARCHITECTURE 11.4.  Two of them must be written *while the
run happens* and cannot be reconstructed afterwards, which is
why capture is its own step rather than part of the harvest.

```
# 1. Dispatch time, written by the wingbeat into the run dir.
#    The size signature is known BEFORE the run (it derives from
#    the makeinput inputs and the structure), and the execution
#    config is what the wingbeat is about to launch.  Neither
#    survives in any output file, so both are recorded here.
function capture_dispatch(unit, wingbeat_dir, execution):
    write_toml(join(wingbeat_dir, "resource_capture.toml"),
        { "signature": size_signature_of(unit),
          "execution": execution.knobs })

# 2. Build time, emitted once by CMake (TODO C78).
#    build_info.toml carries BOTH layers: the coarse bucketed
#    knobs, and the verbatim compile string and library detail.

# 3. After the run, from the scheduler's accounting:
#    sacct -> MaxRSS (peak memory), disk high-water, Elapsed.

# 4. Optionally, imago's own per-phase timings, self-reported.
```

`secular_dimension` is recorded **directly from imago** where
imago reports it, with `basis_function_count` and
`wavefunction_components` kept alongside so the derived value
can be cross-checked against the authoritative one.  DESIGN
8.9 leaves the choice of provenance open; the schema records
it directly and keeps the primitives, which is what makes the
cross-check possible either way.

### 16.7 Harvest and promote (DESIGN 8.7)

```
function harvest_flight(workspace_root, db_root, registry):
    # Walk a finished flight, build one Observation per run
    # directory, write each to staging/<fingerprint>/.
    build_info = load_build_info(workspace_root)  # both layers
    for wingbeat_dir in run_directories(workspace_root):
        capture = read_toml(
            join(wingbeat_dir, "resource_capture.toml"))
        accounting = read_scheduler_accounting(wingbeat_dir)

        # The outcome is what decides whether the resources are
        # measurements or bounds (16.3), so it is resolved FIRST
        # and everything downstream reads it.
        outcome = classify_outcome(accounting)
        resources = assemble_resources(accounting, outcome,
            phase_timings = read_self_report(wingbeat_dir))

        obs = Observation(
            observation_id = <filled by observation_slug>,
            generated_at   = now_utc_iso8601(),
            source         = "flight",
            outcome        = outcome,
            hardware_fingerprint = hardware_fingerprint(
                probe_attributes(accounting.hostname)),
            signature  = SizeSignature(**capture["signature"]),
            execution  = ExecutionConfig(capture["execution"]),
            build      = BuildConfig(build_info.coarse),
            resources  = resources,
            provenance = Provenance(
                compile_string = build_info.compile_string,
                library_detail = build_info.library_detail,
                curator        = "resource_harvest.py", ...))
        save_observation(obs, db_root, area = "staging")
```

```
function classify_outcome(accounting):
    # The scheduler already knows.  An OOM kill and a timeout
    # are DATA (16.3), not failures to be dropped: the first
    # says this configuration is insufficient at this size, the
    # second bounds the walltime from above.  Only a Fortran
    # abort unrelated to resources carries no cost signal.
    if accounting.state == "OUT_OF_MEMORY":  return "oom"
    if accounting.state == "TIMEOUT":        return "timeout"
    if accounting.exit_code != 0:            return "failed"
    return "completed"
```

The curator's promoter mirrors 15.7's four modes (`dry-run`,
`all`, `auto-promote`, interactive), and promotion is a move
of the file from `staging/<fingerprint>/` to
`entries/<fingerprint>/`, so contents and provenance never
change.  One rule DESIGN 8.7 states outright:

```
function may_promote(obs):
    # A `failed` run is a Fortran abort unrelated to resources.
    # It carries no usable cost signal and is staged only for
    # diagnostics -- NEVER promoted (DESIGN 8.7).
    if obs.outcome == "failed":
        return False
    # `completed`, `oom`, and `timeout` all carry signal: a
    # measurement, a lower bound, an upper bound.
    return True
```

DESIGN 8.7 says the curator promotes "with the same discipline
as section 7.8," but 7.8's discipline includes an *objective
auto-promote rule* (a flatness test on the staged file alone),
and no such rule is stated for a cost observation.  What would
one even test?  A cost has no analogue of convergence.  So the
`auto-promote` mode above has no criterion beyond
`may_promote`, and that is a gap in DESIGN, not an omission
here: see 16.9.

### 16.8 Provisioning consumer and day-1 (DESIGN 8.6, 8.8)

The near-term consumer.  The flight layer asks the predictor
what a proposed configuration will cost, applies a margin, and
emits the scheduler request -- the resource-and-cost sibling of
what 15.6 does for k-density.

```
function provision(dataspace, fingerprint, build, signature,
                   execution):
    p = predict(dataspace, fingerprint, build, signature,
                execution)

    if p.is_cold_start:
        # DESIGN 8.8: a fresh fingerprint cannot be predicted
        # for.  Ask generously, run, harvest -- and the run
        # SEEDS the fingerprint.  The artifact improves
        # monotonically from here; the first job pays for it.
        return conservative_cold_start_request(signature)

    margin = SAFETY_MARGIN
    if p.is_borrowed:
        # Fitted on a neighbouring machine (16.5 rung 3).  The
        # prediction is not wrong so much as untested here, so
        # widen rather than trust.  It is surfaced, never
        # silent (DESIGN 8.5).
        margin = widened(margin)

    return SchedulerRequest(
        memory   = p.peak_memory_bytes * margin,
        disk     = p.disk_bytes        * margin,
        walltime = p.walltime_seconds  * margin)
```

A small `manual` seed -- a handful of hand-entered
observations spanning the size range on the local machine --
accelerates day 1, exactly as the section-15 seed flight (C75)
bootstraps convergence guidance.  This is TODO C82.

### 16.9 What this section deliberately does not specify

DESIGN 8.9 leaves four questions open.  Two of them reach into
the predictor, so the pseudocode above stops at them rather
than guessing past them.  Writing a plausible-looking
algorithm here would make the chain *look* complete while
leaving the code with nothing to be checked against -- which
is the failure the chain exists to prevent.

1. **How censored observations enter the fit.**  `fit_group`
   (16.5) excludes them and says so.  That is a placeholder,
   not an answer: an OOM lower bound is real evidence about
   the cost surface and throwing it away wastes it, while
   feeding it in as a point measurement biases every exponent.
   The honest options are a censored (Tobit-style) regression
   or a bound-weighting scheme, and DESIGN 8.9 says the choice
   is calibrated after the seed.  **C77 cannot implement
   `fit_group` faithfully until DESIGN settles this.**  The
   library's schema, loader, emitter, and fingerprint (16.1 -
   16.4) do not depend on it and can be built now.

2. **The exact correction-term form.**  `feature_row` (16.5)
   gives the v1 terms DESIGN 8.6 names by quantity, in one
   place, to be re-derived from the seed.  The pseudocode
   fixes *where* the form lives, not *what* it is.

3. **Aggregate versus per-rank memory.**  Whether
   `peak_memory_bytes` is a job aggregate, a per-node, or a
   per-rank figure changes how the parallel correction is
   modelled and how `sacct` / `time` / self-report are
   reconciled (ARCHITECTURE 11.8).  The schema records one
   number; which number it is must be settled before the
   parallel term means anything.

4. **Build effects on numerics, not just cost.**  Whether the
   build block is ever referenced from the section-7
   convergence side -- against the no-cross-reference boundary
   -- is to be settled in DESIGN, not assumed here.

Two further gaps this pseudocode pass surfaced, which DESIGN 8
does not currently cover:

5. **No objective auto-promote rule for a cost observation.**
   DESIGN 8.7 borrows 7.8's curation discipline, but 7.8's
   `auto_promote_ok` rests on a flatness test that has no cost
   analogue (16.7).  Either a criterion is designed, or
   `auto-promote` should be dropped from the promoter's modes
   and the curator reviews every observation.

6. **`MIN_GROUP_SAMPLES` and `SAFETY_MARGIN` have no values.**
   DESIGN 8.6 says both are calibrated after the seed, which
   is right, and 16.1 names them in one place so calibration is
   a one-line edit.  But a margin is a *safety* parameter: the
   consumer must not ship with it unset, so the seed (C82) is a
   hard prerequisite for the provisioning consumer (C81), not
   merely a source of accuracy.

## 17. Output Control and Citation Banner (DESIGN 10, ARCH 12)

Two facilities, specified together because one gates the other.
`O_Verboseness` decides what the engine is willing to say; the
banner modules are its first and so far only client.

Everything below is Fortran-side.  The line and file references
are to the engine as it stands, and the call sites are named
exactly, because DESIGN 10.6 makes the ordering load-bearing and
an ordering error here is silent -- the banner would test an
uninitialized mask and simply not print.

### 17.1 The module split, and why it is three modules not two

ARCHITECTURE 12.4 maps this onto two modules, `O_Verboseness`
and a single `O_Banner` holding both the identity block and the
method-citation registry.  That map cannot be compiled.  The
identity block prints from `parseCommandLine`, so `O_CommandLine`
must `use` it; the registry's predicates read `kPointIntgCode`,
so it must `use O_KPoints`; and `kpoints.f90:1093` already reads
`use O_CommandLine, only: doSYBD_SCF, ...`.  One module holding
both halves therefore closes a Fortran module cycle:

```
O_CommandLine -> O_Banner -> O_KPoints -> O_CommandLine
```

which no compilation order resolves.  The constraint is factual
and forced, not a preference.

The split that removes it falls on a seam DESIGN 10.2 has
already drawn.  That section separates the two blocks by when
their information becomes available: the identity block knows
everything it needs before any work begins, and the methods
block cannot be written until the run is over.  The same fact
governs what each may depend on.  A module that prints before
the work starts can depend on nothing the work produces, and a
module that prints after it is finished may depend on all of it.
So:

- **`O_Verboseness`** (`verboseness.f90`) -- the mask, the
  category table, and the query.  Depends on nothing.
- **`O_Banner`** (`banner.f90`) -- the identity block: locate
  `banner.txt`, echo it, gated on the `banner` category.
  Depends on `O_Verboseness` only.
- **`O_MethodCitations`** (`methodCitations.f90`) -- the
  registry of DESIGN 10.5 and its predicates.  Depends on
  `O_KPoints` and whatever later state a future predicate
  needs; placed late in the build so it may.

This keeps the property DESIGN 10.5 says matters -- that adding
a method and adding its citation are one task -- which the
alternative loses.  The alternative is to leave the registry
inside `O_Banner` and pass its predicate inputs down as
arguments from a call site that already holds them.  That also
breaks the cycle, but it makes every new reference an edit to a
subroutine signature and to its call in `imago.F90`, so the
citation stops being a local addition.

**ARCHITECTURE 12.4 needs its module list corrected to match.**
This is the one legitimate upward edit: the constraint was
discovered at this level and the higher level is factually
wrong, not merely less detailed.

### 17.2 O_Verboseness: the category table

The table is the single source of truth of ARCHITECTURE 12.2.
Bit positions do not appear in it, because they do not need to:
a category's bit *is* its row index minus one.  There is then no
second column that can drift out of step with the first, and the
translation from name to bit exists in exactly one expression.

```
module O_Verboseness

    numCategories = 1

    # The public contract.  Callers pass these, never a number.
    # The value is a row index, not a bit position; nothing
    # outside this module ever sees a bit position at all.
    integer, parameter :: VERB_BANNER = 1

    character(len=20), dimension(numCategories) :: categoryName
    logical,           dimension(numCategories) :: categoryInNormal

    # Private.  Bit i-1 of the mask is row i of the table.
    integer, private :: verbosenessMask

subroutine initCategoryTable
    categoryName(VERB_BANNER)     = "banner"
    categoryInNormal(VERB_BANNER) = .true.
```

Two names are reserved and may never be used for a category:
`normal` and `none`.  They are not categories but set-valued
aliases -- `normal` expands to every row whose
`categoryInNormal` is true, and `none` expands to the empty set.
Adding a category is one row of `initCategoryTable`, one
`parameter`, and a decision about whether the default includes
it; nothing else moves.

The mask is a default `integer`, so bit 30 is the last usable
position and the table has a ceiling of 31 rows.  That is far
beyond any plausible category count, but it is a real limit and
the parser should not pretend otherwise (17.3).

### 17.3 initVerboseness: parsing IMAGO_VERBOSENESS

```
subroutine initVerboseness

    character(len=1024) :: requestString
    integer :: readStatus, valueLength

    call initCategoryTable
    verbosenessMask = 0

    call get_environment_variable (NAME="IMAGO_VERBOSENESS", &
          & VALUE=requestString, LENGTH=valueLength, &
          & STATUS=readStatus)

    # Unset is not silent (ARCHITECTURE 12.2).  Status 1 means
    # the variable does not exist; a variable that exists but
    # holds only blanks is the same request.
    if ((readStatus == 1) .or. (len_trim(requestString) == 0)) then
        call applyNormalDefault
        return
    endif

    # Status -1 means the value was longer than the buffer and
    # has been truncated, which would silently drop trailing
    # categories.  Say so, then parse what did arrive.
    if (readStatus == -1) then
        warn to unit 20: "IMAGO_VERBOSENESS is", valueLength,
              "characters and was truncated to 1024; trailing"
              " categories were ignored."
    elseif (readStatus /= 0) then
        warn to unit 20: "could not read IMAGO_VERBOSENESS"
              " (status", readStatus, "); using the normal"
              " default."
        call applyNormalDefault
        return
    endif

    # Split on commas and set one bit per recognized token.
    for each comma-delimited token in requestString:
        candidate = lowercase(trim(adjustl(token)))
        if (len_trim(candidate) == 0) cycle    # ",," or trailing
        if (candidate == "none") cycle         # contributes nothing
        if (candidate == "normal") then
            call applyNormalDefault             # OR-ed in, not assigned
            cycle
        endif
        matched = .false.
        do row = 1, numCategories
            if (candidate == trim(categoryName(row))) then
                verbosenessMask = ibset(verbosenessMask, row - 1)
                matched = .true.
                exit
            endif
        enddo
        if (.not. matched) then
            warn to unit 20: "unrecognized IMAGO_VERBOSENESS"
                  " category '", trim(candidate), "' ignored."
        endif

subroutine applyNormalDefault
    do row = 1, numCategories
        if (categoryInNormal(row)) then
            verbosenessMask = ibset(verbosenessMask, row - 1)
        endif
    enddo
```

Four parsing decisions are settled here rather than left to the
code, because each is the kind of thing that gets decided by
accident otherwise.

**Tokens combine by union, and nothing subtracts.**  `none` is
the empty set, so `none,banner` is exactly `banner` -- it is not
a mute switch that a later token overrides.  A caller that wants
silence writes `none` alone.  Likewise `normal,banner` is
`normal`, since `banner` is already in it.

**Matching is case-insensitive.**  `BANNER` and `Banner` are the
same request.  ARCHITECTURE 12.2 requires an unrecognized name
to be reported, and a case difference is not the kind of mistake
that report is for; treating it as one would produce a warning
that names a category the user can plainly see in their own
environment.

**Empty tokens are skipped in silence.**  A trailing comma is a
typo with no ambiguity about intent, and warning about it would
train the reader to ignore the warnings that matter.

**Nothing here stops the run.**  Every failure path warns to
unit 20 and continues, per ARCHITECTURE 12.2: a mistyped
environment variable must not kill a queued cluster job hours
after it was submitted.  This is the opposite of how
`elementData.f90:76` treats a missing `IMAGO_DATA`, and the
difference is deliberate -- missing element data makes the run
impossible, whereas a bad verboseness request only makes the log
wrong.

### 17.4 isVerbose: the query

```
logical function isVerbose (category)

    integer, intent(in) :: category    # a table row, e.g. VERB_BANNER

    if ((category < 1) .or. (category > numCategories)) then
        warn to unit 20: "isVerbose called with out-of-range"
              " category", category, "-- treating as off."
        isVerbose = .false.
        return
    endif

    isVerbose = btest (verbosenessMask, category - 1)
```

The bounds test catches an out-of-range literal, and that is all
it can catch.  ARCHITECTURE 12.2's rule that callers never pass
a bare number is not mechanically enforceable -- a hardcoded `1`
is indistinguishable from `VERB_BANNER` at run time -- so it
remains a review discipline.  Worth stating plainly, since a
guard that catches part of a problem invites the belief that it
catches all of it.

### 17.5 O_Banner: echoing banner.txt

```
module O_Banner

    use O_Verboseness, only: isVerbose, VERB_BANNER

subroutine echoBannerFile

    character(len=100) :: dataDirectory
    character(len=100) :: bannerFileName
    character(len=132) :: lineBuffer
    integer, parameter :: bannerUnit = 314
    integer :: openStatus, readStatus

    call get_environment_variable (NAME="IMAGO_DATA", &
          & VALUE=dataDirectory, STATUS=openStatus)
    if (openStatus /= 0) return          # see the note below
    bannerFileName = trim(dataDirectory)//"/banner.txt"

    open (unit=bannerUnit, file=bannerFileName, &
          & form='formatted', status='old', IOSTAT=openStatus)
    if (openStatus /= 0) then
        write (20,*) "Could not open ", trim(bannerFileName)
        write (20,*) "Continuing without the identity block."
        return
    endif

    do
        read (bannerUnit, fmt='(a)', IOSTAT=readStatus) lineBuffer
        if (readStatus /= 0) exit
        write (20, fmt='(a)') trim(lineBuffer)
    enddo

    close (bannerUnit)
```

Four details carry the whole section.

**The read and write formats are `'(a)'`, never list-directed.**
This is the concrete Fortran form of the warning in DESIGN 10.3
that a reader must not "fix" the whitespace.  `read` under
`fmt='(a)'` preserves leading blanks, which carry the kerning
and the centring; `write` under `fmt='(a)'` starts at column one.
A list-directed `write (20,*)` would insert a leading blank on
every line, shifting the whole butterfly one column right of the
`fmt='(a51)'` rules that `O_TimeStamps` prints directly beneath
it.  The result would look almost right, which is worse than
looking wrong.

**`trim` is safe and `adjustl` is not.**  `trim` removes trailing
blanks only, which are insignificant.  `adjustl` would left-
justify and destroy the art.  The distinction is easy to lose
when tidying this routine later.

**The buffer is 132 characters, not 51.**  DESIGN 10.3 fixes the
*artwork* at 51 columns to match `opLabels`, and the art obeys
that -- its widest line is 50.  The citation lines beneath it do
not and need not: the longest today is 64 characters, the one
naming the repository URL.  A `character(len=51)` buffer would
silently truncate the citation the block exists to deliver, and
the DOI line of TODO A12 will be longer still.

**A missing `IMAGO_DATA` returns quietly.**  It cannot happen in
practice: `Imago` calls `initElementData` before
`parseCommandLine` (`imago.F90`), and that routine stops the run
outright when `IMAGO_DATA` is unreadable, so a run that reaches
the banner has already proved the variable good.  The test is
there so the routine is honest on its own terms rather than
relying on a caller two files away, and it is silent because
anything it could say has already been said by the code that
actually failed.

The unit number is 314 because 313 belongs to `elementData.f90`,
9 to `potential.f90`, and 20 to the log.  Any free number does;
the point is that it was checked.

### 17.6 printIdentityBlock

```
subroutine printIdentityBlock

    if (.not. isVerbose(VERB_BANNER)) return

    call echoBannerFile
    write (20,*)
    call flush (20)
```

The gate is the only thing this routine adds, and it is the
reason `initVerboseness` must already have run (DESIGN 10.6).
The trailing blank line and the `flush` follow what
`timeStampStart` already does, so the banner and the first
timestamp rule are separated the same way every later pair is.

No version string is composed here.  `banner.txt` is echoed
whole, and the version and DOI live in it as text (DESIGN 10.4),
to be filled in by TODO A12.

### 17.7 O_MethodCitations: the registry and its predicates

DESIGN 10.5 specifies a registry pairing each reference with a
predicate answering "did this run use it."  The predicates read
state the engine already holds.

```
module O_MethodCitations

    numMethodRefs = 2
    maxRefLines   = 4

    integer, parameter :: METHOD_MONKHORST_PACK = 1
    integer, parameter :: METHOD_BLOECHL_LAT    = 2

    character(len=76), dimension(maxRefLines,numMethodRefs) :: refText
    integer,           dimension(numMethodRefs) :: refNumLines

subroutine initMethodCitations
    refNumLines(METHOD_MONKHORST_PACK) = 3
    refText(1,METHOD_MONKHORST_PACK) = &
      & "  H. J. Monkhorst, J. D. Pack, ""Special points for"
    refText(2,METHOD_MONKHORST_PACK) = &
      & "  Brillouin-zone integrations,"" Phys. Rev. B 13, 5188"
    refText(3,METHOD_MONKHORST_PACK) = &
      & "  (1976).  DOI: 10.1103/PhysRevB.13.5188"
    ... likewise Bloechl, Jepsen and Andersen, Phys. Rev. B 49,
    ... 16223 (1994), DOI 10.1103/PhysRevB.49.16223

logical function methodWasUsed (methodRef)

    use O_KPoints, only: kPointStyleCode, kPointIntgCode

    select case (methodRef)
    case (METHOD_MONKHORST_PACK)
        # Style 1 is an explicit mesh plus shift and style 2 a
        # minimum density; both build a Monkhorst-Pack mesh.
        # Style 0 is a bare list of k-points the user supplied,
        # which is not necessarily one (kpoints.f90:20-30).
        methodWasUsed = ((kPointStyleCode == 1) .or. &
                       & (kPointStyleCode == 2))
    case (METHOD_BLOECHL_LAT)
        methodWasUsed = (kPointIntgCode == 1)
    case default
        methodWasUsed = .false.
    end select

subroutine printMethodsBlock

    call initMethodCitations

    # The log must still be open.  Writing to a closed unit does
    # not fail -- Fortran reconnects it to fort.20 and truncates,
    # taking the run's whole log with it (DESIGN 10.6).  Say so
    # where it can still be seen, and write nothing.
    inquire (unit=20, opened=logIsOpen)
    if (.not. logIsOpen) then
        write (6,*) "The log was closed before the citations"
              " could be written; they are omitted."
        return
    endif

    # Say nothing at all rather than print an empty invitation.
    if (no methodWasUsed(ref) is true) return

    write (20,*)
    write (20,fmt='(a51)') &
      & '***************************************************'
    write (20,*) "Methods exercised by this run.  Please cite:"
    write (20,*)
    do methodRef = 1, numMethodRefs
        if (.not. methodWasUsed(methodRef)) cycle
        do line = 1, refNumLines(methodRef)
            write (20,fmt='(a)') trim(refText(line,methodRef))
        enddo
        write (20,*)
    enddo
    call flush (20)
```

The registry is not gated on verboseness (DESIGN 10.7): it is a
handful of lines, and it is the part a reader is meant to copy
into a paper.

**Two entries, not three.**  DESIGN 10.5 names Rappe and
co-workers alongside Monkhorst-Pack and Bloechl, but the UFF
parameters have no presence in the engine -- `UFF` does not
appear anywhere in `src/imago/`.  That reference belongs to the
force-field path of DESIGN 4.8, which is Python, and so do
Cornell, Jorgensen, and the LAMMPS reference beside it in the
`## References` list.  Adding a registry entry no predicate can
ever select would create exactly the dead weight DESIGN 10.8
worries about, in the first version of the file.  If those
methods should announce themselves, it is `make_reactions.py`
and `create_lammps_files` that must do it, and that is a
separate task on the Python side.

Both entries restate a citation that `## References` in DESIGN
already carries.  The duplication is accepted for the same
reason 10.4 accepts it for the citation text -- the engine
cannot read a design document -- and the two must be edited
together.

### 17.8 Call sites

**The identity block, in `parseCommandLine`
(`commandLine.f90:91`).**  The two new calls go between the
`open` and the `timeStampStart` that currently follows it on the
next line:

```
    open (20,file='fort.20',status='unknown',form='formatted')

    call initVerboseness        # MUST come first: it sets the
    call printIdentityBlock     #   mask the banner tests.

    call timeStampStart (24)
```

The slot is forced from both sides.  Nothing may print earlier
because unit 20 does not exist earlier, and nothing may print
later because `timeStampStart(24)` writes the first rule of the
log, which the banner must sit above.

The order within the slot is forced too, and this is the one
place in the whole section where a mistake is invisible.
Reversed, `printIdentityBlock` tests a mask that is still zero,
`isVerbose` returns false, and the routine returns without
printing or complaining.  The run succeeds and the banner is
simply absent.  Nothing in the compiler or the output will point
at the cause.

**The methods block, at the end of `Imago` (`imago.F90`).**  It
goes after the `doLoEn` branch and immediately before
`end subroutine Imago`, which is the last point at which every
branch that could have exercised a method has run.

**Relocating the close of the log, in `imago.F90`.**  That last
point is, as the engine stands, after the log has already been
closed, so the call above cannot be added on its own.  Three
`close (20)` statements are removed and replaced by one:

```
    # cleanUpSCF:  close (20) under "if (doPSCF < 0)"  -- REMOVE
    # cleanUpPSCF: close (20)                          -- REMOVE
    # loen:        close (20) before opening fort.2    -- REMOVE

    # The fort.2 completion signal moves for the same reason
    # (DESIGN 10.6).  All three of these certify success from a
    # routine that knows about one stage only, and cleanUpSCF's
    # fires before the post-SCF stage has even begun.
    # cleanUpSCF:  open (unit=2,file='fort.2',...)      -- REMOVE
    # cleanUpPSCF: open (unit=2,file='fort.2',...)      -- REMOVE
    # loen:        open (unit=2,file='fort.2',...)      -- REMOVE

    # end of Imago, after every branch above.  The order is the
    # meaning: nothing may follow the signal.
    call printMethodsBlock
    close (20)
    open (unit=2,file='fort.2',status='unknown')
```

DESIGN 10.6 gives the reasoning; what matters here is that this
is a correction to existing code, not an addition beside it.
Leaving those closes in place and having the methods block
reopen the file would also work, and is the wrong repair: the
log would still be opened at one level and closed at three, and
the next thing appended to the end of a run would meet the same
trap.

Nothing else depends on either being early.  The conditional on
`cleanUpSCF`'s close exists only to keep the log alive for a
post-SCF stage, which closing once at the end does
unconditionally and better; and a run that reaches the end of
`Imago` at all is exactly the run whose success `fort.2` is
meant to certify.

The relocation is not a design change.  DESIGN 6.1.2 already
states that `fort.2` certifies the binary ran without an
abortive error, and `imago.py` already treats it that way. The
code does not implement that, which makes this the ordinary
case of code disagreeing with a specification, repaired in the
direction the chain requires.

### 17.9 Build wiring

`verboseness.f90`, `banner.f90`, and `methodCitations.f90` are
added to the two source lists that build an engine --
`src/imago/real/CMakeLists.txt` (target `imagoG`) and
`src/imago/complex/CMakeLists.txt` (target `imago`).  The
`auxiliary` targets do not use `O_CommandLine` and need none of
this.

Those lists are ordered by module dependency, and the ordering
is what 17.1 is about, so it is not free to choose:

- `verboseness.f90` before `commandLine.f90`, and it may go
  immediately after `kinds.f90` since it depends on nothing.
- `banner.f90` after `verboseness.f90` and before
  `commandLine.f90`.
- `methodCitations.f90` after `kpoints.f90`, and before
  `imago.F90` which calls it.  Anywhere in that range works.

`src/data/banner.txt` is already in the `DATABASES` list of
`src/data/CMakeLists.txt` and installs to `share` with the rest,
so no build change is needed for the artwork itself.

### 17.10 What this section does not specify

**The dead `banner` variable.**  `O_TimeStamps` declares
`character(len=51) :: banner` at `timeStamps.f90:24` and nothing
reads it.  It predates this work and should be deleted in the
same change that adds `O_Banner`, before the two can be confused
for each other.  It is noted here rather than specified because
deleting an unread variable needs no design.

**Retrofitting the existing writes.**  ARCHITECTURE 12.3 puts
the existing unconditional `write (20,...)` calls out of scope,
and they stay out of scope here.  One category exists and one
call site consults it.

**The categories themselves.**  Everything the debugging and
parallelization campaign will want -- SCF iteration detail,
integral diagnostics, resource reporting, developer trace -- is
deliberately absent, per ARCHITECTURE 12.3.  The shape above is
what makes that cheap later: a row, a parameter, and a gate.

Two questions this pass could not close.

1. **Whether `printMethodsBlock` should fire when the run
   aborts.**  Every `stop` in the engine bypasses the end of
   `Imago`, so a run that dies in the secular equation prints no
   methods block, having already exercised the mesh it would
   have cited.  That is arguably correct -- a failed run has no
   results to cite -- but it is a consequence of where the call
   sits rather than a decision anyone made, and DESIGN 10.6 does
   not address it.

   The first draft of this section had a worse version of the
   same blind spot, and it is worth recording how it was caught.
   It placed the call at the end of `Imago` without noticing
   that the log is closed before that point, so the first real
   run wrote its citations into a freshly truncated `fort.20`
   and destroyed every line the calculation had produced.  The
   code implemented the pseudocode faithfully; the pseudocode
   was wrong.  Nothing in the chain would have found it, because
   every level had been checked against the level above and the
   error was introduced at the bottom of the stack of documents
   and inherited upward from the code's actual behaviour.  Only
   running it found it.  That is the argument for the guard now
   in 17.7 and for treating "it compiles and the section reads
   correctly" as the weakest kind of evidence.

2. **How a method added without a predicate is caught.**  This
   is DESIGN 10.8's open question and the pseudocode does not
   answer it.  Worth recording what this section did *not* do
   about it: the split of 17.1 puts `O_MethodCitations` next to
   the engine state its predicates read, which makes a missing
   predicate easier for a reader to notice, and nothing more.
   It remains true that a method can be added to Imago and go
   uncited with no symptom at all.

---

## 18. POPTC Decomposition Index (DESIGN 11)

Before any transition is computed, a partial optical run must
decide which partial each basis function belongs to.  This
section specifies that assignment for all four offered cells.
The IBZ correction applied afterwards is section 7a.

### 18.1 What is built

Five things, all sized from the decomposition request.  The
first three are what the accumulation and the output need; the
last two describe the *layout* and exist so that section 7a can
carry a partial through a symmetry operation.

- `sumNumPartials` -- how many partials this run produces.
  The stored pair matrix is this squared, so it is also the
  cost driver of DESIGN 11.4.  Note that a type-grouped
  request is not automatically the small one: `numAtomTypes`
  is fixed and tiny for a crystal or a defect supercell, but
  in an amorphous cell a type is an environment bin and the
  count grows with the cell (DESIGN 11.4).
- `pOptcIndex(valeDim)` -- for each basis function, the
  partial it contributes to.  This is the whole of the
  decomposition; everything downstream just accumulates
  through it.
- `partialsIndex(sumNumPartials)` -- how many basis functions
  feed each partial.  The Kramers-Kronig consumer needs it to
  normalize, since a complete set of partials must carry the
  additive constant of `eps1 = 1 + (2/pi) Int[...]` exactly
  once between them rather than once each.
- `segmentBase(numSegments + 1)` -- where each segment's block
  of partials begins, as a zero-based offset, with the extra
  final entry holding `sumNumPartials`.  A segment is a type
  for the type-grouped cells and a site for the atom-grouped
  ones.
- `slotsPerSegment(numSegments)` -- how many partials each
  segment owns, which is `segmentBase(s+1) - segmentBase(s)`.
  Stored rather than recomputed because section 7a walks it
  directly.

**The last two must outlive this routine**, which is the only
reason they are listed here at all.  Section 7a builds its
`partialPerm` table by taking each site's block of partials and
re-basing it onto the block of the site's image under a
symmetry operation, and that re-basing is expressible only in
terms of these two arrays.  A version of this walk that treated
them as scratch would leave 7a with nothing to build from.

### 18.2 One walk, two parameters

DESIGN 11.2 defines a cell by two independent choices, and
the assignment mirrors that directly rather than branching
per code.  The two parameters are:

- **segment key** -- what a partial belongs to.  The atom's
  TYPE index for the type-grouped cells, the atom's SITE
  index for the atom-grouped ones.
- **slots per segment** -- ONE for the total-resolved cells,
  and one per radial function (per QN_nl pair) for the
  nl-resolved ones.

Writing it as one parameterized walk rather than four
branches is worth doing for a reason beyond tidiness: it
makes the grid of DESIGN 11 visible in the code, so that a
cell added later is a new parameter value rather than a new
branch that must remember to do everything the others do.

The walk itself is indifferent to what a type means -- it
reads `atomTypeAssn` and groups by it.  What that grouping
is *worth* to the reader of the output is not indifferent,
and DESIGN 11.6 is where that is set out: for a defect
supercell in particular, a type-grouped partial averages the
defect neighbourhood into the bulk, because the types were
assigned from the pre-defect symmetry on purpose.

```
# Resolve the request into the two parameters.
#   detailCodePOPTC 0 does not reach here at all.
if detailCodePOPTC in (1, 2):  grouping = TYPE
else:                          grouping = ATOM
if detailCodePOPTC in (1, 3):  resolution = TOTAL
else:                          resolution = NL

if grouping == TYPE: numSegments = numAtomTypes
else:                numSegments = numAtomSites

# Lay out the partials: each segment owns a contiguous block,
#   segmentBase records where each block starts, and
#   slotsPerSegment records how long it is.  Both are kept
#   past the end of this routine for section 7a.
segmentBase(1) = 0
for s = 1, numSegments:
   typeOfSegment = (grouping == TYPE) ? s
                                      : atomTypeAssn(s)
   if resolution == TOTAL:
      slots = 1
   else:
      slots = sum over l of
              atomTypes(typeOfSegment)%numQN_lValeRadialFns(l)
   slotsPerSegment(s) = slots
   segmentBase(s+1)   = segmentBase(s) + slots

sumNumPartials = segmentBase(numSegments + 1)

# Assign every basis function.  The walk is over SITES in all
#   four cells, because that is the order the basis functions
#   are laid out in; only the destination differs.
partialsIndex(:) = 0
valeDimIndex = 0
for site = 1, numAtomSites:
   currentType = atomTypeAssn(site)
   segment = (grouping == TYPE) ? currentType : site

   slot = 0
   for l = 1, lAngMomCount:                 # 1=s 2=p 3=d 4=f
      for radialFn = 1, atomTypes(currentType)%
                        numQN_lValeRadialFns(l):

         # A total-resolved segment has a single slot that
         #   every radial function shares.  An nl-resolved
         #   one advances, so the slots run s, p, d, f in
         #   the order the basis is laid out.
         if resolution == TOTAL: slot = 1
         else:                   slot = slot + 1

         partial = segmentBase(segment) + slot

         for m = 1, (l-1)*2 + 1
            valeDimIndex = valeDimIndex + 1
            pOptcIndex(valeDimIndex) = partial
            partialsIndex(partial) = partialsIndex(partial) + 1
```

**Note what this removes.**  Counting `partialsIndex` by
incrementing once per assigned basis function is correct for
every cell, including the type-grouped ones where several
sites feed one partial.  A type-grouped count obtained
instead by measuring runs of equal `pOptcIndex` would be
correct only while every atom of a type is contiguous in the
basis ordering.  That happens to hold -- Imago sorts sites by
type for reasons that have nothing to do with this code --
but depending on it means an unrelated change to site
ordering breaks the normalization silently.  The increment
does not care.

**What the segment index means, and why 7a can assume it.**
The walk assigns `segment = currentType` for the type-grouped
cells and `segment = site` for the atom-grouped ones, so for
codes 3 and 4 the segment index *is* the atomic site index and
`segmentBase` and `slotsPerSegment` are indexed by site
directly.  Section 7a relies on exactly that when it writes
`segmentBase(site)`, and it is entitled to: the correction runs
only for `detailCodePOPTC >= 3`, which is the atom-grouped
half.  For codes 1 and 2 those two arrays are indexed by type
instead and no correction is applied, so the two readings never
have to coexist.

### 18.3 The withdrawn cell

There is no branch for a QN_nlm resolution.  Assigning each
basis function its own partial is a one-line change to the
walk above, and it must not be made: DESIGN 11.2 gives two
independent reasons, the binding one being that the pair
matrix would be `valeDim` squared.  The `QN_mLetter` naming
tables that only an nlm-resolved printer would need are
correspondingly absent.

### 18.4 What the printer walks

`printSpectrumPOPTC` must visit the partials in the same
order this section lays them out, since the file's sequence
numbers are the only thing tying a spectrum to its label.
The nesting is therefore the same doubly-parameterized walk,
run over the pair (initial partial, final partial):

```
for each initial partial p in layout order:
   for each final partial q in layout order:
      write SEQUENCE_NUM, the two labels, then the
        TOTAL/x/y/z columns for that pair
```

Labels follow the grouping.  A type-grouped partial names its
type; an atom-grouped one names its site.  An nl-resolved
partial appends the QN_l letter and the radial function
number.  The total spectrum is always written first, as
sequence number 1, before any pair.

---

## 19. LAT Optical Integration (DESIGN 12)

The optical properties accumulate their broadened spectrum
by Gaussian smearing over IBZ k-points.  This section
specifies the tetrahedron alternative that DESIGN 12 adds
alongside it.  Both pathways remain, selected by
`kPointIntgCode`, and neither replaces the other.

The quantity integrated is a joint density of states over
band PAIRS: the surface where `e_j(k) - e_i(k)` equals the
output energy, weighted by the squared momentum matrix
element.  Section 8 solves the same shape of problem for the
partial DOS, and this section follows its two-pass
structure, with the band-energy DIFFERENCE where section 8
has a band energy and the squared matrix element where it
has a Mulliken projection.

### 19.1 What is dispatched, and what is not

DESIGN 12.2 restructures the optical path to match the DOS
path rather than imitate it.  `optcCond` and
`optcCondPOPTC` move to module scope in O_OptcSpectra,
`printOptcResults` splits in two, and `subroutine optc`
calls the halves in turn -- which is where and how
`subroutine dos` already makes the same choice.

```
# subroutine optc, in imago.F90.  Compare subroutine dos
#   immediately above it: same shape, same place.
call computeOptcSpectra(doOPTC)
call printOptcSpectra(doOPTC)


# computeOptcSpectra, in O_OptcSpectra.  Setup, then the
#   pathway choice.  Both branches leave the same module
#   arrays filled, so the printer never learns which ran.
build energyScale
allocate optcCond (and optcCondPOPTC if decomposing)

if kPointIntgCode == 1:
   latFactor = the DESIGN 12.4 normalization
   call accumulateOptcCond_LAT(latFactor)
   if detailCodePOPTC /= 0:
      call accumulateOptcCondPOPTC_LAT(latFactor)
else:
   kPointFactor = kPointWeight * 0.5 / sigmaSqrtPi
                  / hartree / spin
   call accumulateOptcCond(kPointFactor, sigma)
   if detailCodePOPTC /= 0:
      call accumulateOptcCondPOPTC(kPointFactor, sigma)
```

**On the names.**  `getOptcCond` does not get anything; it
accumulates into an array handed to it, so it and its
partner are renamed to say so.  The `_LAT` suffix marks the
integration method and the unsuffixed name is the Gaussian
one, following section 3's convention for
`electronPopulation`.  Note that the Gaussian meaning of an
unsuffixed name is a convention rather than something the
name states.

**What is shared and what is not.**  The momentum matrix
elements and the decomposition index of section 18 are
shared: both pathways need the same physics from the same
eigenvectors.  What diverges is how the resulting transition
strengths are indexed, filtered, ranged and occupied --
section 19.2 -- so the producers differ even though the
physics inside them does not.  Downstream nothing changes:
the accumulators leave the same arrays the printer already
consumes, and the `1/E` scaling and unit conversion stay in
`printSpectrum` where they are.

### 19.2 Pass 1: re-index the matrix elements by band pair

**This pass exists because the Gaussian path's storage
cannot be reused.**  `transitionProb` is indexed by a
position in a list sorted by transition energy, and the
band identity is discarded by that sort; the pair count and
the band ranges also vary per k-point.  A tetrahedron needs
the same band pair at all four corners, so the LAT path
stores under a band-pair index instead (DESIGN 12.4).

**The physics is shared; only the bookkeeping differs.**
Five things separate this producer from `computePairs`, and
only one of them is the sort: the storage slot, the
cutoff filtering, the band range, the occupation source, and
the ordering.  None of them touches the momentum matrix
element itself.  So the construction of `conjWaveMomSum` --
the sum over basis functions of the conjugated wave function
against the momentum matrix, which is the expensive part and
the part where an error would be hardest to see -- is
extracted into a routine both producers call.

### 19.2.1 What this pass consumes, and from where

This producer attaches to a running program rather than
standing alone, so its inputs are listed before its
algorithm.  Each row is the answer to "who fills this, and
when": the rows are what the code below is entitled to
assume, and a later reader checks them against the source
rather than trusting the prose.

```
quantity            supplied by                 when
--------------------------------------------------------
valeVale,           readDataSCF / readDataPSCF, once per
valeValeMM            called in                 (spin,
                      computeTransitions          k-point)
energyEigenValues   the secular solution,       resident
                      dimensioned (numStates,     for ALL
                      numKPoints, spin)           k-points
firstOccupiedState  getEnergyStatistics,        before the
  and companions      per (k-point, spin)         loop
electronPopulation  computeElectronPopulation   NOT on the
  _LAT                _LAT, called only from      optical
                      subroutine bond             path
tetrahedra,         generateTetrahedra, from    before optc
  tetraVol            initializeKPoints when      runs
                      kPointIntgCode == 1
fullKPToIBZKPMap,   the IBZ fold in             before optc
  fullKPToIBZOpMap    initializeKPointMesh        runs
pOptcIndex,         buildPOPTCIndex, called    before the
  segmentBase,        once from                   k-point
  slotsPerSegment,    computeTransitions          loop
  sumNumPartials,
  partialPerm
```

**A row's "when" has two halves, and the second one is where
this section went wrong.**  Every row above says when a
quantity is FILLED.  For anything the code releases, when it
is FREED is equally part of the contract, and a consumer
that runs in a later phase can be starved by a release that
looked correct beside the loop that built it.  The rule that
holds is that **a lifetime spans every consumer**, not that
the routine owning a loop owns what the loop used.

Two quantities here run past the k-point loop and must not
be freed with it:

  - `pairIsWanted`, read by both accumulators to skip band
    pairs that were never filled.
  - `partialPerm`, read by the decomposed accumulator of
    19.4 to permute the two partial indices at each corner.

Both are therefore released by `subroutine optc` after the
spectra have been written, alongside the banded stores they
describe.  `cleanUpPOPTCIndex` is called from there for the
same reason, rather than at the end of `computeTransitions`.

Two of these rows decide the structure, and they are the
reason this section can no longer be written as a routine
that owns a loop over k-points.

**The momentum matrix is loaded per k-point, inside a loop
this pass does not own.** `computeTransitions` reads it and
then calls a producer once, and `buildConjWaveMomSum`
consumes what was read rather than reading anything itself.
A producer owning its own k-point loop would have to
duplicate that read, including its SCF and post-SCF variants
and its PACS matrix codes -- two copies of the logic that
decides what to read, free to drift apart.  **So the banded
producer is called once per k-point from inside the existing
loop, exactly where `computePairs` is called today**, and
the array it fills is allocated before the loop because it
spans every k-point.  This also satisfies for free the
one-call-per-read rule that the Gamma build needs, since the
in-place Hermiticity fix runs once per read either way.

**`electronPopulation_LAT` is absent on this path.** It is
built inside `subroutine bond`, which `subroutine optc` never
reaches, so an optics-only run finds it unallocated. The
optical path calls it once, before the k-point loop, after
the eigenvalues have been shifted to put the Fermi level at
zero -- which is why the energy argument is zero.

```
function buildConjWaveMomSum(firstFin, lastFin,
                             initComponent, finComponent,
                             conjWaveMomSum):
    # Exactly the loop that computePairs runs today, lifted
    #   unchanged.  Both pathways call it; neither owns it.
    #
    # It takes NO k-point or spin argument. The wave
    #   functions and the momentum matrix are module arrays
    #   that the caller has already loaded for the k-point
    #   in hand, so passing an index would imply a lookup
    #   this routine does not perform.
    #
    # The POPTC counterpart, buildConjWaveMomSumPOPTC,
    #   carries the extra partial index of section 18 and is
    #   extracted the same way.
    ...


function computeTransProbBanded(h, numKPoints,
        firstInitAll, lastInitAll,
        firstFinAll,  lastFinAll):

    # PRECONDITION, and it is not met by the existing call
    #   sequence.  computeElectronPopulation_LAT is invoked
    #   from subroutine bond, which subroutine optc never
    #   reaches; an optics-only run leaves the array
    #   unallocated.  The optical path calls it here, after
    #   the eigenvalues have been shifted so the Fermi level
    #   sits at zero, which is why the argument is zero.
    if kPointIntgCode == 1 and
            .not. allocated(electronPopulation_LAT):
       call computeElectronPopulation_LAT(0.0)

    # The band ranges are per k-point arrays.  The banded
    #   store must span the UNION over k-points, so that a
    #   band pair present at any corner has a slot at every
    #   corner.  Taking the per-k range here is the error
    #   that would leave holes at exactly the k-points a
    #   metal's Fermi surface passes through.
    initLo = min over k of firstInitAll(k, h)
    initHi = max over k of lastInitAll(k, h)
    finLo  = min over k of firstFinAll(k, h)
    finHi  = max over k of lastFinAll(k, h)

    # INDEX ORDER IS DELIBERATE (DESIGN 12.4).  Fortran
    #   stores the leftmost index fastest.  Section 19.3
    #   holds the band pair fixed and walks tetrahedra,
    #   fetching four different k-points per tetrahedron,
    #   so kIBZ belongs immediately after the component and
    #   NOT last.  The whole block (:, :, i, j) is then
    #   contiguous -- three components by the k-point count,
    #   tens of kilobytes -- and stays cache resident across
    #   every tetrahedron for that band pair.  Writing the
    #   more obvious (dim3, i, j, kIBZ) instead strides by
    #   dim3*nOcc*nUnocc on each corner and misses on every
    #   one.  Do not reorder this without reordering 19.3.
    allocate transProbBanded(dim3, numKPoints,
                             initLo:initHi, finLo:finHi)
    transProbBanded = 0.0

    for kIBZ = 1 to numKPoints:
       conjWaveMomSum = buildConjWaveMomSum(kIBZ, h, ...)

       for i = initLo to initHi:
          for j = finLo to finHi:

             # The squared momentum matrix element times the
             #   initial state's occupancy and the final
             #   state's vacancy.  Occupancies come from
             #   electronPopulation_LAT on this path, not
             #   from the Gaussian electronPopulation, so
             #   that one scheme sets both the geometry and
             #   the filling (DESIGN 12.4c).
             #
             # No cutoff test here.  computePairs drops
             #   pairs failing the transition-energy cutoff,
             #   which is safe when each k-point stands
             #   alone and unsafe here: a pair may fail at
             #   one corner of a tetrahedron and pass at the
             #   other three.  The pruning is done instead by
             #   the pre-pass of 19.2.2, over all k-points at
             #   once, and reaches this loop as pairIsWanted.
             for c = 1 to dim3:
                transProbBanded(c, kIBZ, i, j) =
                      |M_ij^c(kIBZ)|^2
                      * occ(i, kIBZ) * (1 - occ(j, kIBZ))

    return transProbBanded
```

### 19.2.1a The decomposed producer

Everything above describes the undecomposed store.  A partial run
needs both: the total spectra are wanted whether or not a
decomposition was requested, exactly as on the Gaussian side, so
the undecomposed producer runs in both cases and the decomposed
one runs in addition.

```
function computeTransProbPOPTCBanded(kIBZ, h):
    # Same loops, same occupancies, same pruning mask.  What
    #   changes is that the matrix element is resolved by partial
    #   instead of collapsed, so conjWaveMomSum comes from
    #   buildConjWaveMomSumPOPTC and carries the final-state
    #   partial as an extra index.
    for i, j in the wanted band pairs:
       for c = 1 to dim3:

          # Resolve by the partial the INITIAL basis function
          #   belongs to.  The final-state partial is already
          #   carried by conjWaveMomSum.
          pairMatrix(:,:) = 0
          for basisFn = 1, valeDim:
             for n = 1, sumNumPartials:
                pairMatrix(pOptcIndex(basisFn), n) +=
                      wavefn(basisFn, i)
                      * conjWaveMomSum(basisFn, n, jSlot, c)

          # The totals that turn a sum of squares back into a
          #   squared sum, which is what makes the partials add
          #   up to the undecomposed transition probability.
          sumRe = sum of the real parts of pairMatrix
          sumIm = sum of the imaginary parts of pairMatrix

          # n OUTER, o INNER.  o is the leftmost index of both
          #   the scratch matrix and the store, so it runs
          #   innermost.
          for n = 1, sumNumPartials:
             for o = 1, sumNumPartials:
                transProbPOPTCBanded(o, n, c, kIBZ, i, j) =
                      (Re(pairMatrix(o,n)) * sumRe
                       + Im(pairMatrix(o,n)) * sumIm)
                      * occ(i, kIBZ) * (1 - occ(j, kIBZ))
```

**Store layout, and why its rule is weaker than 19.2's.**  The
order is `(o, n, dim3, kIBZ, i, j, spin)`.  The same reasoning
applies -- leftmost is fastest, choose against the consuming loop
-- but it cannot reach the same conclusion, because this array is
far too large to hold one band pair's slice in cache.  So the aim
narrows to keeping the pair matrix for one component at one
corner contiguous, and running the partial loops leftmost-index
innermost.  The accumulation's destination stays scattered
whatever is done, since permuting the two partial indices is the
entire purpose of the operation (DESIGN 12.4).

**The decomposition index is built once, before the loop.**
Section 18 describes it as a standalone construction and it now
is one: `buildPOPTCIndex` fills `pOptcIndex`, `segmentBase`,
`slotsPerSegment`, `sumNumPartials` and `partialPerm`.  Two
reasons it cannot stay inside a producer.  The tetrahedron
store is sized from `sumNumPartials`, which is therefore needed
before the loop begins; and both producers read the index while
neither owns it, so a producer that freed it would be guessing
whether it was the last to run.

**It is released after the SPECTRA are written, not after the
loop ends.**  `cleanUpPOPTCIndex` is called from `subroutine
optc`, which is the routine that owns the whole optical phase.
Releasing it at the end of `computeTransitions` instead is
wrong, and wrong silently: the accumulation of 19.4 runs later
and reads `partialPerm`, so the table it needs would already be
gone, its `allocated` guard would read false, and the corner
permutation would be skipped without any complaint.  A guard
written to protect a legitimately absent table cannot tell that
case apart from a table freed too early -- which is why the
lifetime has to be settled here rather than left to the guard.

The tempting shorter rule, that the routine owning the loop
owns the lifetime, is what produced that error.  It is right
about who must not free the index (the producers inside the
loop) and wrong about who must (the phase, not the loop).

### 19.2.2 Cost, measured against the Gaussian producer

**Pass 1 costs almost exactly what the Gaussian pass costs,
which is not obvious and is worth showing.**  The dominant
term is building `conjWaveMomSum`: three components by the
final-state count by `valeDim` squared, per k-point.  That
term is IDENTICAL in both pathways, because the array is
built for the whole final-state range BEFORE any cutoff is
applied.  `computePairs` prunes with two `exit` statements
inside the pair loop, so what the cutoff actually saves is
the per-pair dot product, three components by `valeDim`
each -- not the expensive part.

Putting numbers on it, for the KNbO3 case where `valeDim` is
90 with roughly 20 occupied and 50 unoccupied states in
range: the shared build is about 1.2 million multiply-adds
per k-point, while the complete unpruned set of pairs costs
about 270 thousand.  **So a producer that prunes nothing at
all pays at most about twenty percent more than one that
prunes perfectly.**

Where pruning does matter is STORAGE, and there it matters a
great deal for the decomposed case.  Each retained pair
costs three doubles per k-point for the total spectra, but
three times `sumNumPartials` squared for the pair matrix --
about 31 kilobytes per pair per k-point at 36 partials.

So the pruning is worth keeping, and it can be recovered
without touching a single matrix element:

```
# Cheap pre-pass, eigenvalues only.  Cost is one subtraction
#   per (i, j, k), which is nothing beside a valeDim squared
#   build, and it runs before any matrix element is formed.
#
# A band pair is worth storing if its transition energy
#   falls within the cutoff at ANY k-point.  Taking the
#   minimum over k is what makes this safe for tetrahedra: a
#   pair that is out of range at one corner and in range at
#   another must still be present at both.
for i, j in the union band ranges:
   pairIsWanted(i,j) = ( min over k of
                         (eigen(j,k,h) - eigen(i,k,h)) )
                       <= maxTransEnergy
```

**Pass 2 is where the real cost sits, and it hinges on one
implementation choice.**  The Gaussian accumulation is
`numKPoints x pairs x W_gauss`, where `W_gauss` is the
number of energy bins a single transition actually reaches
before the exponent cuts it off.  The LAT accumulation is
`numTetrahedra x pairs x W_tetra`, where `numTetrahedra` is
six times the FULL mesh count and `W_tetra` is the number of
bins spanned by the four corner values of the energy
difference.

Those widths are only real if the energy loop is BOUNDED.
As written today, `getOptcCond` loops over every energy
point and tests inside the loop, so with the 5001-point
grid these runs use and a Gaussian reaching about 140 bins,
roughly thirty-five evaluations in thirty-six are discarded.
That is tolerable at `numKPoints` of four.  It would not be
tolerable here: the same mesh gives 384 tetrahedra against
those four k-points, so the unbounded form would do about
ninety-six times the outer iterations, each sweeping all
5001 points.  **Compute the first and last energy index from
`sortedDiff(1)` and `sortedDiff(4)` and loop only that
range.  This is not an optimization to add later; the
unbounded form is not viable.**

With the loop bounded, the ratio of LAT to Gaussian work is
six times the IBZ reduction factor times `W_tetra` over
`W_gauss`.  For the KNbO3 mesh that is roughly seventy.

**But the scaling favours LAT, and that is the argument for
paying the constant.**  `W_gauss` is fixed by the broadening
width, so the Gaussian cost grows as the k-point count, or
as the cube of the points per axis.  `W_tetra` is set by how
much the energy difference varies between ADJACENT mesh
points, which falls off as the mesh is refined, so the
tetrahedron count rising as the cube is partly cancelled and
the LAT cost grows roughly as the square.  LAT is the more
expensive method on a coarse mesh and the cheaper one on a
fine mesh, with a crossover that is a measurement rather
than a prediction.  The caveat is that this assumes the
energy difference is roughly linear between adjacent points,
which is the same assumption the tetrahedron method itself
rests on -- so where the assumption fails, the accuracy
argument fails with the cost argument.

### 19.3 Pass 2: tetrahedron accumulation, total spectra

```
function accumulateOptcCond_LAT(transProbBanded, eigenValues,
        tetrahedra, numTetrahedra, tetraVol,
        fullKPToIBZKPMap, fullKPToIBZOpMap,
        energyScale, numEnergyPoints, h):

    optcCond = 0.0

    # The LAT replacement for kPointFactor, derived in
    #   DESIGN 12.4.  tetraVol sums to 1 while kPointWeight
    #   sums to 2, so sum(kPointWeight) restores the scale.
    #   The Gaussian normalization 1/(sigma*sqrt(pi)) is
    #   absent because the corner density weight IS the
    #   delta function, evaluated exactly.  The other three
    #   factors belong to the quantity rather than to the
    #   integration and are carried over unchanged -- the
    #   0.5 especially, which stops the k-point weights from
    #   counting two electrons per state a second time.
    #   Dropping it would halve every spectrum.
    latFactor = sum(kPointWeight) * 0.5
                / hartree / real(spin)

    for i = initLo to initHi:
       for j = finLo to finHi:
          for T = 1 to numTetrahedra:

             for c = 1 to 4:
                kFull(c) = tetrahedra(c, T)
                kIBZ(c)  = fullKPToIBZKPMap(kFull(c))
                opIdx(c) = fullKPToIBZOpMap(kFull(c))

                # The DIFFERENCE band is what the corner
                #   weights are built from.  Both
                #   eigenvalues are read at the same IBZ
                #   representative, and e(Rk) = e(k), so no
                #   permutation enters here.
                epsDiff(c) = eigenValues(j, kIBZ(c), h)
                           - eigenValues(i, kIBZ(c), h)

             sigma      = argsort(epsDiff)
             sortedDiff = epsDiff(sigma)

             for iE = 1 to numEnergyPoints:
                E = energyScale(iE)
                if E < sortedDiff(1) or E >= sortedDiff(4):
                   cycle

                cornerDOSWt_LAT(1:4) =
                      bloechlCornerDOSWt(E, sortedDiff)

                for c = 1 to 4:
                   # The weights come back in SORTED corner
                   #   order.  sigma(c) carries each one
                   #   back to the corner it belongs to, and
                   #   the matrix element must be fetched
                   #   for THAT corner.  Fetching for corner
                   #   c instead pairs a weight with the
                   #   wrong k-point and yields a plausible
                   #   spectrum rather than a broken one
                   #   (DESIGN 12.4a).
                   orig = sigma(c)
                   kIc  = kIBZ(orig)

                   # kIc is the SECOND index, so this slice
                   #   sits inside the block already held
                   #   for this band pair (19.2).
                   optcCond(:, iE, h) +=
                         cornerDOSWt_LAT(c) * tetraVol
                         * latFactor
                         * transProbBanded(:, kIc, i, j)

    return optcCond
```

**What is deliberately absent: the Cartesian rotation.**
The momentum operator is a vector, so the components at a
full-mesh corner are mixed by the operation `opIdx(orig)`
relative to the IBZ representative.  The loop above does not
apply that mixing, exactly as the Gaussian path does not.
The isotropic column is unaffected, because summing the
three components is a trace and the mixing cancels
(PSEUDOCODE 7a proves this for the pair matrix and the
argument is the same here); the per-axis columns remain
unverified.  **This is TODO O3, and this loop is where its
fix belongs** -- rotate `transProbBanded(:, i, j, kIc)` by
`opIdx(orig)` at the fetch.  It is not attempted here so
that the integration change can be validated on its own.

### 19.4 The partial counterpart

Identical in structure, with the pair matrix carried through
and both of its indices permuted at the corner.

```
function accumulateOptcCondPOPTC_LAT(
        transProbPOPTCBanded, ...):

    optcCondPOPTC = 0.0
    latFactor     = sum(kPointWeight) * 0.5
                    / hartree / real(spin)

    ... same band, tetrahedron, corner and energy loops ...

                for c = 1 to 4:
                   orig = sigma(c)
                   kIc  = kIBZ(orig)
                   R    = opIdx(orig)

                   # b OUTER, a INNER.  The partial store is
                   #   laid out (a, b, dim3, kIBZ, i, j), so
                   #   a is the fastest index and must be
                   #   the innermost loop.  This array is
                   #   far too large to hold a band pair's
                   #   slice in cache, so the aim here is
                   #   narrower than in 19.3: keep the walk
                   #   over the pair matrix sequential.  The
                   #   destination stays scattered whatever
                   #   is done, because the permutation is
                   #   the whole point (DESIGN 12.4).
                   for b = 1 to sumNumPartials:
                      bRot = partialPerm(R, b)
                      for a = 1 to sumNumPartials:
                         aRot = partialPerm(R, a)

                         optcCondPOPTC(aRot, bRot, :, iE, h)
                            += cornerDOSWt_LAT(c) * tetraVol
                               * latFactor
                               * transProbPOPTCBanded(
                                    a, b, :, kIc, i, j)
```

`partialPerm` is built by section 18 and is unchanged: the
same table that the star average used serves here.  The
deposit-forward direction matches section 7a -- the partial
that IBZ index `a` represents AT THIS CORNER is
`partialPerm(R, a)` -- so the two blocks express one map and
a reader can check them against each other.

**This loop is the reason the index outlives the k-point
loop** (19.2.1a).  It runs from `computeOptcSpectra`, well
after the producers have finished, and it is the last reader
of `partialPerm`.  A release placed beside the construction,
or beside the loop that filled the stores, leaves this line
indexing a table that is gone.

For the type-grouped detail codes 1 and 2 no permutation is
needed, exactly as in 7a, and `partialPerm` is not built.
Guard the two inner lines with the same
`detailCodePOPTC >= 3` threshold and use `a` and `b`
directly otherwise.

### 19.5 What happens to the star average

**The section 7a block is not called on this path.**  It is
not moved and not modified; it simply does not run.

The reason is worth stating where a reader will look for the
block and fail to find it.  7a exists because the Gaussian
path visits only IBZ representatives and must spread each
one's contribution over the members of its star.  The loop
above visits full-mesh corners directly, and applies the
operation once per corner as the matrix element is fetched.
That is the same arithmetic reaching the same answer by a
shorter route, so applying both would double-count the
symmetry.

7a remains live and necessary on the Gaussian path.  Both
pathways must be correct under IBZ reduction; they reach it
differently.

### 19.6 Guards and checks

**Degenerate corners.**  When the four `epsDiff` values
coincide the analytic denominators vanish.  Section 2a's
guards apply unchanged, but the case is not rare here:
parallel bands make `epsDiff` flat by construction, and
parallel bands are what produce the sharp structure an
optical spectrum is computed to show (DESIGN 12.4d).  Test
this branch directly rather than assuming it is inherited.

**The band-pair union, not the per-k range.**  Stated in
19.2 and repeated because it fails silently: a store built
from one k-point's `firstOccupiedState` leaves holes at
other corners.

**Cross-checks available.**

- The partials must still sum to the total, since the
  accumulation is linear in the pair matrix.  This holds
  under any permutation of the two indices and so proves
  nothing about the unfolding -- the same blindness section
  7a records.
- The two pathways must converge to the same spectrum as
  the mesh is refined.  This is the real check, and DESIGN
  12.5 governs how to run it: `sigmaOPTC` means different
  things on the two paths, so hold its MEANING fixed rather
  than its value.
- A gapped system must give a spectrum that is identically
  zero below the gap on both paths.  Under LAT this is
  exact rather than approximate, since no broadening leaks
  weight below `sortedDiff(1)`.

---

## 20. Symmetrize Atom-Resolved Results (DESIGN 1.7)

Section 1's decomposition is carried onto itself by the point
group for lattices whose operations permute the mesh axes up
to sign, and not for hexagonal or rhombohedral ones.  This
closes the remainder by averaging the finished result over
the group.  It runs on the tetrahedron pathway only; the
Gaussian star average already distributes each irreducible
point's contribution evenly over its star.

**Where it does NOT go.**  Effective charge and bond order
need nothing.  `computeElectronPopulation_LAT` pools every
corner's weight onto the corner's IBZ representative and
`computeBond` spreads it back over the star divided by the
star size, which IS the star average.  Adding this there
would average an already-averaged quantity.  The quantities
that need it are the energy-resolved ones, which attach a
weight to an individual corner: the LAT PDOS of section 8.3
and the LAT optical accumulators of section 19.

```
function symmetrizeChannels(values, permTable, numOps,
                            numChannels, numEnergyPoints):

    # values(channel, energyPoint), modified in place.
    #
    # No orbit is enumerated. Summing over EVERY operation
    #   is the orbit average, because the operations
    #   carrying a channel onto a given orbit member form a
    #   coset and so contribute that member equally often.
    #   This is also why the direction of the permutation
    #   does not matter: channelPermTbl is built from
    #   invAtomPerm and partialPerm from atomPerm, and both
    #   give the same average when R runs over the group.

    # The averaged values must be built from the
    #   UNSYMMETRIZED ones, so the sum cannot go in place.
    #   One energy point at a time keeps the scratch to a
    #   single channel vector.
    allocate scratch(numChannels)

    for iE = 1 to numEnergyPoints:
        scratch = 0.0
        for R = 1 to numOps:
            for alpha = 1 to numChannels:
                scratch(alpha) += values(permTable(R, alpha),
                                         iE)
        values(:, iE) = scratch / real(numOps)


function symmetrizePairs(values, permTable, numOps,
                         numPartials, numEnergyPoints):

    # values(a, b, dim3, energyPoint), modified in place.
    #   Both indices are permuted by the same operation,
    #   which is what preserves the meaning of a pair.
    allocate scratch(numPartials, numPartials, dim3)

    for iE = 1 to numEnergyPoints:
        scratch = 0.0
        for R = 1 to numOps:
            # b OUTER, a INNER: a is the leftmost index of
            #   both the scratch slab and the store.
            for b = 1 to numPartials:
                bRot = permTable(R, b)
                for a = 1 to numPartials:
                    aRot = permTable(R, a)
                    scratch(a, b, :) +=
                          values(aRot, bRot, :, iE)
        values(:, :, :, iE) = scratch / real(numOps)
```

**Call sites, and the window each must sit in.**

- LAT PDOS: after `integratePDOS_LAT` fills `pdosComplete`
  and BEFORE the shared output phase writes it, while
  `channelPermTbl` is still allocated -- `computeDOS`
  releases that table at the end of the routine.  Applies
  to `detailCodePDOS` 1 and 2 only; mode 0 is a type-level
  sum and already invariant, and mode 3 is refused on this
  pathway (DESIGN 1.4).
- LAT optical: after the accumulators fill
  `optcCondPOPTC` and before `printOptcSpectra`, which is
  inside the window where `cleanUpPOPTCIndex` keeps
  `partialPerm` alive (19.2.1).  Applies to
  `detailCodePOPTC` 3 and 4 only; codes 1 and 2 are type
  grouped and already invariant.

**Guard, and say so out loud.**  For k-point style code 0 the
symmetry maps are never built, so `atomPerm` and therefore
both permutation tables are unallocated.  Guard on the
table's allocation, skip the averaging, and write a line to
fort.20 saying it was skipped and why.  Silence here would be
indistinguishable from having done the work.

**Report the spread before averaging.**  For each group of
channels the averaging will merge, compute the largest
deviation among them first and write the maximum over all
groups to fort.20.  An imposed equality that leaves no trace
cannot be told from an earned one, and this also gives every
run a free measurement of the residual asymmetry -- which is
what makes `SYMMETRIZE_LAT_PARTIALS 0` a rarely-needed
diagnostic rather than the only way to see the number.

**Checks.**

- Totals must not move at all.  Summing the averaged
  weights over k-points returns the original sum, so the
  TDOS and the total spectra are unchanged to round-off.
  A total that shifts means the permutation table is not a
  permutation, and that is the first thing to test.
- Symmetry-equivalent atoms must agree to round-off
  afterwards, on any lattice, including the hexagonal case
  the decomposition alone cannot fix.
- The partials-sum-to-total identity proves nothing here.
  It holds under any permutation of the indices, which is
  the same blindness sections 7a and 19 record.

---

## 21. Optical Cartesian Components (DESIGN 13)

The momentum operator is a vector, so unfolding a matrix
element onto a star member has to ROTATE its three Cartesian
components and not merely permute atoms.  The code squares
the matrix element before anything could rotate it, and a
squared modulus cannot be rotated, so this section is a
change of what is stored rather than an extra step at the
deposit.

### 21.1 What this consumes, and from where

Each row answers "who fills this, and when", and -- because
that omission is what broke section 19 -- also when it is
freed.

```
quantity            supplied by            lifetime
------------------------------------------------------------
abcRealPointOps     computeRealPointOps,   built in
                      already conjugated     initializeKPoints
                      into the LOADED cell   for style codes
                      basis (DESIGN 2.7)     0, 1, 2; never on
                                             the SYBD branch;
                                             freed in
                                             cleanUpKPoints
realVectors         O_Lattice, read from   resident for the
                      the structure file     whole run.  Its
                                             COLUMNS are the
                                             lattice vectors --
                                             the k-point file's
                                             CONV_LATTICE block
                                             uses ROWS, and the
                                             two must not be
                                             confused
invert3x3           already in kpoints.f90 no lifetime; a
                      beside                 helper
                      computeRealPointOps
xyzRealPointOps     NEW, built here        same branches and
                                             same lifetime as
                                             abcRealPointOps
fullKPToIBZOpMap    the IBZ fold in        before optc runs;
                      initializeKPointMesh   absent for style
                                             code 0
fullKPToIBZKPMap    same                   same
cartesianCodeOPTC,  readOptcControl, from  read once, BEFORE
  cartesianCode-      fort.5's               any store is
  POPTC               OPTC_INPUT_DATA        sized, because the
                                             stores are
                                             dimensioned from
                                             them
```

**One precondition is not met by the existing call
sequence** and must be arranged rather than assumed. The
stores are allocated in `getEnergyStatistics`, which runs
before `computeTransitions`; the direction codes are read in
`parseInput`, which runs before both. So the ordering is
already correct and needs only to be stated -- but a later
reader must not move the allocation earlier.

### 21.2 The Cartesian rotations

```
function computeXyzPointOps():
    # abcRealPointOps acts on LOADED-cell fractional
    #   coordinates. A fractional vector f sits at Cartesian
    #   r = L f, with L the lattice vectors as columns, so
    #   the same operation in Cartesian form is
    #   R_xyz = L R_abc L^-1.
    #
    # This starts from abcRealPointOps rather than from
    #   convAbcPointOps, so it inherits the full/prim
    #   cell-mode handling already done and does not repeat
    #   that branch.
    allocate xyzRealPointOps(3, 3, numPointOps)

    call invert3x3(realVectors, realVectorsInv)

    for opIdx = 1 to numPointOps:
        xyzRealPointOps(:,:,opIdx) =
            matmul(realVectors,
                   matmul(abcRealPointOps(:,:,opIdx),
                          realVectorsInv))
```

**Check this before trusting anything built on it.** Every
`xyzRealPointOps` must be orthogonal to round-off:
`R R^T = I`. That is not automatic -- it holds only because
the operation really is a symmetry of the lattice -- so a
failure means the conjugation is wrong or the cell is not
what the space group says. Test it once at construction and
stop loudly, because every downstream symptom of getting
this wrong looks like a physics result.

### 21.3 What the stores hold

The direction code decides, and the two codes are read
independently for totals and partials (DESIGN 13.7).

```
  code 0   isotropic only.  Store ONE real per transition:
             sum over c of |M^c|^2.  No rotation is ever
             needed, because that sum is invariant.  A third
             of today's storage.
  code 1   diagonal.  Store THREE COMPLEX: M^c itself.
  code 2   full symmetric tensor.  Store the same three
             complex numbers.  Levels 1 and 2 differ only in
             what is formed at the deposit and how wide the
             output record is.
```

So the producers change like this, with `computePairs`
standing for all four:

```
    valeValeXMom(k) = sum(valeVale(:,i,1)
                          * conjWaveMomSum(:,finalStateIndex,k))

    if cartesianCode == 0:
        # Collapse immediately. The occupancy factors are
        #   real and fold in here as they do today.
        store(pair) = sum over k of
                      (real(valeValeXMom(k))**2
                       + aimag(valeValeXMom(k))**2)
                      * initStateFactor * finStateFactor
    else:
        # Keep the complex element. Occupancies are NOT
        #   folded in: see below.
        store(k, pair) = valeValeXMom(k)
```

**Occupancy factors stay out of the stored element.**
Folding them in as a square root would be wrong, because the
final-state vacancy `1 - occupancy` can go slightly negative
through round-off near a Fermi surface and the square root
of a negative number is a crash rather than a rounding
error. They are real scalars per (initial band, final band,
k-point), so carry them in their own small array and apply
them at the deposit where those indices are already in hand.

### 21.4 The tetrahedron pathway deposit

This pathway already visits full-mesh corners and applies
that corner's operation as it fetches, so the rotation joins
the fetch and nothing moves.

```
    for c = 1 to 4:                      # tetrahedron corner
       orig = sigma(c)
       kIc  = kIBZ(orig)
       R    = xyzRealPointOps(:,:, opIdx(orig))

       # Rotate, THEN square. The other order is the defect.
       rotated(:) = matmul(R, storedElement(:, kIc, i, j))

       occupancy = initFactor(i, kIc) * finFactor(j, kIc)

       if cartesianCode == 1:
          for d = 1 to 3:
             optcCond(d, iE) += weight * occupancy
                                * abs(rotated(d))**2
       else:                              # code 2
          for (d, e) in the six upper-triangle pairs:
             optcCond(pairIndex(d,e), iE) += weight
                * occupancy
                * real(rotated(d) * conjg(rotated(e)))
```

Code 0 needs no rotation at all: the stored scalar is
already the invariant sum, and the deposit is what it is
today with one component instead of three.

### 21.5 The Gaussian pathway and the star-summed rotation

**This pathway has no star loop, and that is the whole
difficulty.** It accumulates over irreducible points with
the multiplicity carried inside `kPointWeight`, so there is
nowhere for a per-member rotation to live.

Rather than add the loop -- which would multiply the
accumulation by the reduction factor, 4 to 48 -- precompute
the star AVERAGE of the rotation product, once per
irreducible k-point:

```
function buildStarRotationAverage():
    # For each irreducible k-point, average R_c1,d R_c2,e
    #   over the operations that carry it to its star
    #   members. This depends on the star alone: not on the
    #   band, not on the energy, not on the matrix element.
    allocate starRotationAvg(3, 3, 3, 3, numKPoints)
    starRotationAvg = 0.0

    for kFull = 1 to numFullMeshKP:
        kIBZ = fullKPToIBZKPMap(kFull)
        R    = xyzRealPointOps(:,:, fullKPToIBZOpMap(kFull))
        starSize(kIBZ) += 1

        for c1, c2, d, e in 1..3:
            starRotationAvg(c1,c2,d,e, kIBZ) +=
                R(c1,d) * R(c2,e)

    for kIBZ = 1 to numKPoints:
        starRotationAvg(:,:,:,:, kIBZ) /= starSize(kIBZ)
```

**Why four component indices and not three** (DESIGN 13.5).
A diagonal entry has both factors of the squared modulus
carrying the same component, and for direction codes 0 and 1
that is all that is ever asked for. Direction code 2 asks
for the off diagonal tensor entries too, where the two
factors carry different components, so the second index `c2`
cannot be tied to `c1`. Setting `c2 = c1` recovers the
diagonal form, so codes 0 and 1 read their entries out of
this same array unchanged.

If there is no folded mesh -- an explicitly listed k-point
set, where each point stands for itself -- there is no star
to average over. Fill the identity entries and return; the
accumulation is then exactly what it was before this change.

**It must be the AVERAGE and not the sum.** `kPointWeight`
already carries the star multiplicity, so a star SUM here
would count the symmetry twice -- the same trap section 19.5
records for applying both the star average and the
per-corner permutation. An unreduced mesh gives every star
one member, the average is the identity operation's own
product, and the accumulation reduces to what it does today.

The deposit then becomes

```
    for d = 1 to 3:
       for e = 1 to 3:
          T(d,e) = storedElement(d, ...) 
                   * conjg(storedElement(e, ...))

    # One slot per accumulated spectrum: one at code 0,
    #   three at code 1, six at code 2. The slot names the
    #   pair of components (c1, c2) it is built from, in the
    #   order 21.7 declares: xx, yy, zz, xy, xz, yz.
    for slot = 1 to numAccumComp:
       c1 = firstOfPair(slot)
       c2 = secondOfPair(slot)
       optcCond(slot, iE) += kPointFactor(kIBZ) * broadening
          * occupancy
          * sum over d, e of
            starRotationAvg(c1,c2,d,e, kIBZ) * real(T(d,e))
```

Form the strength ONCE per transition, outside the energy
loop. It does not depend on the energy point, and the energy
grid runs to thousands of points, so building it inside
would repeat the whole contraction for every one of them.

### 21.6 The decomposed case, and the loop that moves

The partial store holds a cross-term between one pair's
matrix element and the total over all pairs:

  Re[ M_c(o,n) conjg(S_c) ],  S_c = sum over (o',n') of M_c

Both factors carry the component index, so both are rotated.
The store therefore keeps the complex `M_c(o,n)` AND the
complex total `S_c`, which is one vector per (initial band,
final band, k-point, spin) rather than one per pair and so
costs nothing beside the pair matrix.

```
    rotatedPair(:) = matmul(R, storedPair(:, o, n, ...))
    rotatedTotal(:) = matmul(R, storedTotal(:, ...))

    for d = 1 to 3:
       optcCondPOPTC(o, n, d, iE) += weight * occupancy
          * real(rotatedPair(d) * conjg(rotatedTotal(d)))
```

**The section 7a star average must move INSIDE the component
loop.**  As written, `computePOPTCPairs` walks the star with
the component index held fixed outside it, which is correct
only while an operation cannot mix components:

```
    for pairIndex, for component:          # TODAY
       for each star member:
          permute the two partial indices

    for pairIndex:                         # REQUIRED
       for each star member:
          rotate the component index AND
          permute the two partial indices
```

This is the only existing loop nest the change reorders
rather than extends, and it is the first place to look if
the sum rule breaks afterwards.

### 21.7 Output width

`printSpectrum` writes a fixed five-column record and a
five-field header; `printSpectrumPOPTC` writes
`COL_LABELS 4` and `TOTAL x y z` per unit.  Both become
width-aware:

```
  code 0   COL_LABELS 1   TOTAL
  code 1   COL_LABELS 4   TOTAL x y z
  code 2   COL_LABELS 7   TOTAL xx yy zz xy xz yz
```

Those are the widths of epsilon-2 as the engine writes it,
and of epsilon-1 as imagoKKc writes it. The other five
derived files follow `numScalarCols` instead, so at code 2
they keep the code 1 width of four rather than widening to
seven:

```
  code 2   eps2, eps1              7   TOTAL xx yy zz xy xz yz
  code 2   ELF, n, k, R, alpha     4   TOTAL x y z
```

A consumer must therefore take each file's width from that
file, and never from the run's direction code.

The isotropic column stays FIRST in every case, so a
consumer that only wants it needs no knowledge of the code.
`imagoKKc` and `processPOPTC.py` read these files and must
be checked against the declared width rather than assuming
four.

**How each consumer learns the width.** They differ, and the
difference is not cosmetic. `processPOPTC.py` parses the raw
decomposed file, which carries `COL_LABELS`, so it reads the
count and the label text from there and passes both on to
every file it writes. `imagoKKc` is handed the TOTAL
spectrum files as well, and those name their columns in a
header without declaring a count, so it counts the numeric
fields on the first data line instead. Counting the data is
the one method that serves both file kinds, and it cannot
disagree with the data it describes.

```
  numDirCols = (fields on a data line) - 2
                 # drop the energy and the isotropic column
  numWorkCols = 1 if numDirCols == 0 else numDirCols
  numScalarCols = min(3, numWorkCols)
```

`numWorkCols` sizes epsilon-2 and epsilon-1, which carry
every component the input holds: 1, 3 or 6.

`numScalarCols` sizes the five quantities that are functions
of a scalar dielectric function -- the energy loss function,
the refractive index, the extinction coefficient, the
reflectivity and the absorption. It is the DIAGONAL count,
so it is 1, 3 and 3 for the three direction codes (DESIGN
13.7). At code 2 those five are computed from the first
three columns and the off-diagonal three are carried through
epsilon-1 only.

**The isotropic average runs over `numScalarCols`, never
over `numWorkCols`.** At code 2 they differ, and averaging
epsilon-1 over all six would fold the off-diagonal entries
into a quantity that is supposed to be one third of the
TRACE. The result stays finite and plausible and is wrong.

```
  totalEps1(i) = sum(eps1(1:numScalarCols, i))
                 / numScalarCols
```

At code 1 that is the division by three this has always
done. At code 0 it is the identity on the single column;
dividing by a literal three there would shrink every derived
spectrum to a third of its value with no symptom but the
numbers.

**The list-directed read is the trap.** At code 0 a record
holds two values. Reading three would not fail -- Fortran's
list-directed input runs on into the NEXT record to satisfy
the list, shifting every subsequent line by one and yielding
a complete, plausible, wrong spectrum. The read is therefore
branched on the width rather than left to consume what it
needs.

### 21.8 Checks, in the order they are worth running

1. **Orthogonality of the rotations**, at construction
   (21.2).  Everything else is meaningless until it passes.
2. **An unreduced mesh must be unchanged.**  With every star
   of size one the rotations are the identity and the star
   average is trivial, so the whole change must be a no-op
   there.  Any movement means the rotation is being applied
   where it should not be.
3. **A reduced run must now reproduce the unreduced one,
   per axis.**  This is the defect's own test and it has a
   known target: on cubic KNbO3 the unreduced run gives
   48.470 in all three columns while the reduced run today
   gives 37.607 / 72.701 / 85.237.
4. **The isotropic column must not move at all**, on any
   mesh, for any direction code.  It is correct today, so a
   shift means the rotation has broken something that was
   working -- and it is the column every user reads.
5. **Partials must still sum to the total** (section 11).
   Weak, as always: it holds under any permutation of the
   pair indices and so cannot see a wrong rotation.  It can
   see a reordered loop nest, which is exactly what 21.6
   introduces.

## 22. imagoKKc's Output File Lifecycle (TODO O14)

`imagoKKc` opens all eight of its outputs with
`status="new"`, which fails if the file already exists. A
run that dies partway therefore leaves them behind and
BLOCKS every later run, and the failure surfaces as a bare
non-zero exit that names nothing.

### 22.1 The seam: who creates, consumes and removes each file

Read from the callers rather than assumed, because that is
where this project's specification errors have actually
happened.

```
file                created by     consumed by
------------------------------------------------------------
fort.50, fort.51    the optc run   READ by imagoKKc as
  (total eps2)                       status="old"
fort.450            makePDOS.py,   READ by imagoKKc as
  (one pair's         once per       status="old"
   eps2)              pair
fort.100-170        imagoKKc       imago.py safe_move's each
fort.101-171          (spin 1        one into the job
  (total spectra)     and 2)         directory, which
                                     removes the source
fort.500-570        imagoKKc       processPOPTC.copy_data
  (one pair's         (POPTC         appends each into the
   spectra)           path)          accumulated raw file,
                                     then os.remove()s it
```

**So in a normal run nothing is ever left behind**, which is
why `status="new"` has held up for as long as it has. The
window is narrow: a crash BETWEEN `imagoKKc` writing and the
caller consuming. The POPTC path widens it by running
`imagoKKc` once per pair, so an abort at pair seventeen
leaves a full set of eight.

### 22.2 What changes

```
  the three INPUTS   stay status="old"
  the 24 OUTPUTS     become status="replace"
```

**The inputs must not be touched.** Opening a missing input
as anything but `"old"` would CREATE it empty, and imagoKKc
would then read a spectrum of zeros and write derived
spectra of zeros without failing. That is a worse defect
than the one being fixed, and it is one edit away.

Re-running a job is meant to be idempotent. These files are
transient scratch whose consumers remove them, so replacing
a stale one restores the intended behaviour rather than
discarding anything a user wanted.

### 22.3 Say when it happened

Replacing silently would trade a loud wrong behaviour for a
quiet one, and stale output is evidence that an earlier run
died. So `inquire` the outputs before opening them and, if
any existed, write ONE line naming how many were replaced.

```
  count = number of the 8 output files already present
  if count > 0:
      write "imagoKKc: replaced <count> stale output
             file(s) from an earlier run that did not
             finish."
```

One line and not one per file: eight lines per pair across
twenty-six pairs is noise, and the count is what a reader
needs. A clean run prints nothing.

### 22.4 Checks

1. **A normal run is unchanged**, and prints no replacement
   message. Nothing is left behind in a run that completes,
   so the count is zero.
2. **A run after a crash now succeeds** where it previously
   failed. Reproduce by writing a stale `fort.500` by hand
   before invoking `processPOPTC.py`.
3. **The message appears exactly once** for that run, and
   names a non-zero count.
4. **A missing INPUT still fails**, and fails at the open
   rather than by producing zeros. This is the check that
   guards against the one-edit-away defect in 22.2.

## 23. Gaussian Per-Atom PDOS Unfolding (TODO C148)

The Gaussian accumulation deposits straight into
`pdosAccum(pdosIndex(valeDimIndex))` and permutes nothing,
so on a symmetry-reduced mesh every member of a k-point's
star is credited with the REPRESENTATIVE's atom-by-atom
breakdown. `buildChannelPermTable` is called only inside
`if (kPointIntgCode == 1)`, so the table the tetrahedron
path uses was never built for this one.

Live in shipped output, on the default integration setting.
Measured on cubic KNbO3: three symmetry-equivalent oxygens
that must be identical differ by 7.1e-01.

### 23.1 Why a single call at the end is exact

The obvious fix is to permute inside the k-point loop, per
star member. That is not needed, and the reason is worth
writing down because it looks too convenient to be true.

Let `P(k)` be the per-channel vector at irreducible k-point
`k`, and `w_k` its star multiplicity. What is accumulated
today is

```
  A = sum over k of  w_k * P(k)
```

What is WANTED is each representative spread over its own
star,

```
  correct = sum over k of sum over R in the star of k
                            of  R . P(k)
```

`P(k)` is invariant under the little group of `k` -- the
operations that leave `k` fixed also carry the crystal onto
itself, so they can only permute equal projections among
equivalent atoms. Therefore summing over star members is
the same as averaging over the WHOLE point group and
multiplying by the multiplicity:

```
  sum over R in star of R . P(k)  =  w_k * groupAvg(P(k))
```

and because `groupAvg` is linear,

```
  correct = sum over k of w_k * groupAvg(P(k))
          = groupAvg( sum over k of w_k * P(k) )
          = groupAvg(A)
```

**So one group average of the finished accumulation is
exactly the per-k-point star average.** No loop moves and
nothing inside the accumulation changes.

This holds for the accumulated PDOS even though individual
band projections inside a degenerate subspace are NOT
invariant -- the diagonalizer's basis there is arbitrary.
Degenerate bands sit at one energy and so carry one
broadening weight, which makes the sum over a degenerate
multiplet invariant, and that sum is what the PDOS holds.

**This is NOT the same as the LAT path's symmetrization.**
There, correctness comes from permuting per tetrahedron
corner, and `symmetrizePDOS_LAT` is an OPTIONAL cleanup
(DESIGN 1.7) that a user can switch off. Here the group
average IS the correctness mechanism and is not optional.
Wiring it to `symmetrizeLATPartials` would let a user turn
off a fix rather than a polish.

### 23.2 Seam inventory

Every quantity the change consumes, and when it lives.

```
quantity        supplied by          lifetime
------------------------------------------------------------
invAtomPerm     buildInvAtomPerm,    built for EVERY run
                  called from          except SYBD, right
                  imago.py's setup     after buildAtomPerm;
                  and pscf paths       independent of
                                       kPointIntgCode, so it
                                       is already available
                                       here. Guard on
                                       allocated() anyway --
                                       SYBD does not build it
channelPermTbl  buildChannelPerm-    TODAY built and freed
                  Table                only under
                                       kPointIntgCode == 1.
                                       BOTH guards must widen,
                                       or the Gaussian path
                                       leaks it or reads a
                                       table that was never
                                       filled
cumulNumDOS,    allocated before     present on both paths
  cumulDOSTotal   the spin loop        already
pdosComplete    the accumulation     the average must run
                  in the spin loop     AFTER the k-point loop
                                       for that spin closes
                                       and BEFORE the output
                                       phase reads it
```

**The lifetime row is the one that matters.** Freeing a
permutation table before its consumer runs is exactly the
defect that made O9's fix silently do nothing for weeks
(PSEUDOCODE 19.2.1a). Widening the build guard without
widening the release guard leaks; widening the release
without the build reads an unallocated array.

### 23.3 What runs, and for which detail codes

```
  mode 0  type-resolved. NOTHING to do: a sum over all
            atoms of a type is already invariant, and an
            operation carries an atom onto one of the same
            type (DESIGN 2.3)
  mode 1  per-atom total.    group-average the channels
  mode 2  per-atom, per QN_nl. group-average the channels
  mode 3  per-atom, per QN_nlm. NOT COVERED -- see below
```

```
if (kPointIntgCode /= 1) and (detailCodePDOS is 1 or 2)
      and allocated(invAtomPerm):
    build channelPermTbl before the spin loop
    ...
    after the k-point loop for spin h closes:
        symmetrizeChannels(pdosComplete(:,:), channelPermTbl)
    ...
    release channelPermTbl after the spin loop
```

`symmetrizeChannels` already exists in `O_MathSubs` and
already performs exactly this group average; it was written
for the LAT path in C149 and needs no change.

### 23.4 Mode 3 is the same defect and this does NOT fix it

`buildChannelPermTable` has branches for detail codes 0, 1
and 2 and none for 3. That is not an oversight: the m
components of an orbital MIX under a rotation, so permuting
channel indices cannot express what happens to them. The
`D^l(R)` representation matrices are needed, which is
exactly why DESIGN 1.4 REFUSES mode 3 on the tetrahedron
path.

The Gaussian path does not refuse it. So mode 3 on a reduced
mesh is silently wrong TODAY, on both paths, and one of them
says so while the other does not. Making that consistent is
a decision with user-visible consequences -- refusing it
would stop runs that people are doing now -- and it is
deliberately left out of this section rather than settled by
whoever writes the code.

### 23.5 Checks

1. **The three oxygens of cubic KNbO3 must agree**, on a
   reduced mesh, to the precision the effective charge from
   `-bond` already reaches. That control found the defect
   and is in `jobs/knbo3/cubic`; a fresh reproduction is in
   `jobs/knbo3/o9_pdos/gauss`. Spread today is 7.1e-01.
2. **An unreduced mesh must not move at all.** Every star
   has one member, the group average is over the identity
   alone, and the whole change is a no-op there.
3. **The total DOS must not move**, on any mesh. It is a
   sum over channels, which no permutation can change, so
   any shift means the average is being applied where it
   should not be.
4. **Mode 0 must not move**, for the same reason.
5. The partials-sum-to-total identity CANNOT see this
   defect and is not a check. It holds under any permutation
   of the channel indices, which is the whole point.

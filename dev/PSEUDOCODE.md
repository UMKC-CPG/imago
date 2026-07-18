# Pseudocode

> **Document hierarchy:** VISION -> ARCHITECTURE -> DESIGN
> -> **PSEUDOCODE** -> Code. For the design rationale behind
> these algorithms, see `DESIGN.md`.

---

## 1. Generate Tetrahedra (DESIGN 1.2)

```
function generateTetrahedra(nA, nB, nC):
    numTetrahedra = 6 * nA * nB * nC
    allocate tetrahedra(4, numTetrahedra)
    t = 0
    for a = 1 to nA:
        for b = 1 to nB:
            for c = 1 to nC:
                # 8 corners with periodic wrapping
                M1 = idx(a,     b,     c    )
                M2 = idx(a+1,   b,     c    )
                M3 = idx(a,     b+1,   c    )
                M4 = idx(a,     b,     c+1  )
                M5 = idx(a+1,   b+1,   c    )
                M6 = idx(a+1,   b,     c+1  )
                M7 = idx(a,     b+1,   c+1  )
                M8 = idx(a+1,   b+1,   c+1  )

                # 6 tetrahedra sharing diagonal M1-M8
                tetrahedra(:, t+1) = [M1, M2, M5, M8]
                tetrahedra(:, t+2) = [M1, M3, M5, M8]
                tetrahedra(:, t+3) = [M1, M3, M7, M8]
                tetrahedra(:, t+4) = [M1, M4, M7, M8]
                tetrahedra(:, t+5) = [M1, M4, M6, M8]
                tetrahedra(:, t+6) = [M1, M2, M6, M8]
                t = t + 6

function idx(a, b, c):
    # Periodic wrapping, 1-based indexing
    return getIndexFromIndices(
        mod(a-1, nA) + 1,
        mod(b-1, nB) + 1,
        mod(c-1, nC) + 1)
```

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

                for i = 1 to 4:
                    ki = corners(sigma(i))
                    electronPopulation_LAT(
                        n, ki, spin) +=
                        cornerIntgWt_LAT(i)
                            * tetraVol

    return electronPopulation_LAT
```

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
            # position. R preserves species, so only
            # atoms of the same type can match.
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
    #   h = (prod|b_i| / (recipCellVolume * D))^(1/3),
    #   x_i = |b_i| / h.
    h = (recipMag(1) * recipMag(2) * recipMag(3)
         / (recipCellVolume * density)) ^ (1/3)
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
    n = [max(1, round(x(i))) for i in 1..3]

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

The pure Monkhorst-Pack grid in fractional abc coordinates
with equal base weight.  This is the primary product of the
section; 4c.5 optionally folds it.

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
                mesh(:,p) = -1/2
                          + ([i, j, k] - 1 + shift) * delta
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
# In initializeKPoints, after the mesh is built (styles 1
# and 2), write to the main output unit:
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
4c.7 mirror, not for the run itself.  Only the mesh style codes
(1 and 2) emit these; an explicit-list run (style 0) builds no
axial mesh and emits none, so imago.py records them as absent
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
    # its two endpoints' per-atom energies are within the SAME
    # `threshold` the two-sided test uses.  Because a stride adds many
    # k-points, a small change across it is strong evidence the energy
    # has settled -- but only evidence: the refine phase VERIFIES with
    # the full two-sided test, so a coincidentally flat stride (an
    # oscillating near-metal energy) is caught there, not trusted here.
    lo = per_atom_ev(lo_rung.energy, cell_atom_count)
    hi = per_atom_ev(hi_rung.energy, cell_atom_count)
    return abs(hi - lo) < threshold


function at_ceiling(mesh, max_count):
    # The fixed per-axis backstop (DESIGN 3.12.3): the LARGEST axial
    # count reaching max_count, not the product or the stride.  A cost
    # ceiling from the resource dataspace (16) layers on later; the
    # climb stops at whichever bites first.
    return max(mesh) >= max_count
```

### 4e.3 One material's next action

Three search shapes (DESIGN 3.12.3 / 3.12.5) share the two-sided
stop test (4e.2) and the rung rule (4e.1); they differ only in
which rungs they compute.  `climb_next` dispatches on the mode,
threading a per-material search `state` for the stateful
bracket-refine shape; the grid and the unit-step climb ignore it.
Every shape returns one of: `RUN(mesh)` -- run one more mesh;
`CONVERGED(rung)` -- done, that rung; `CEILING` -- stop,
non-converged.

```
function climb_next(rungs, state, config):
    # rungs: the material's computed {mesh, energy}, sorted ascending.
    # state: the per-material search state (bracket-refine only).
    # Returns (action, state').
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
    # at the flat interior rung, or CEILING.
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
    if state.phase == BRACKET:
        top = state.endpoints[-1]
        if length(state.endpoints) == 1:
            # Only the seed is computed; launch the first stride
            #   (stride 1).  No flatness to test yet.
            return strideUp(state, 1, config)
        # Two or more endpoints: test the last stride's flatness.
        prev = state.endpoints[-2]
        if stride_is_flat(rung_at(rungs, prev), rung_at(rungs, top),
                          config.cell_atom_count, config.threshold):
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

    else:  # REFINE -- fill [lo, hi] one position at a time, then test
        gap = nextFillMesh(rungs, state.lo, state.hi,
                           config.classes, config.recipMag)
        if gap is not None:
            return RUN(gap), state              # keep filling
        # Filled: judge the now-consecutive block with the two-sided
        #   test.  Only the bracket's own rungs are passed, so the
        #   sparse bracket endpoints outside it never mislead the
        #   neighbour comparison (DESIGN 3.12.3).
        block = [r for r in rungs
                 if meshSize(state.lo) <= meshSize(r.mesh)
                    and meshSize(r.mesh) <= meshSize(state.hi)]
        idx = pick_converged_climb(block, config.cell_atom_count,
                                   config.threshold, config.flat_needed)
        if idx is not None:
            return CONVERGED(block[idx]), state
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

The numeric knobs that policy reads -- the `confidence_high`
threshold, the two `flat_needed` counts, `grid_width`, the two
`start_offset` values, `max_stride`, the climb-shape choice, and
the per-axis `max_count` ceiling -- are config, not constants
(Principle 11).  They are sourced from the
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
than falling through to a default shape.  The provisional default
values themselves are still to be fixed by the seed experiment
(3.12.6).

### 4e.5 Concurrent orchestration across materials

```
function converge_by_climb(materials, configs, seed_densities,
                           dispatcher, on_non_converged):
    # Drive every material through the climb to a verdict --
    # converged, or non-converged (a ceiling, or a rung that failed
    # to run).  Serial within a material, concurrent across, and NO
    # material waits on another: a chain climbs on the instant its
    # own rung lands (DESIGN 3.12.5).  The injected `dispatcher`
    # owns the in-flight set so this loop tracks only its per-
    # material ladders (Principle 12); it exposes two calls (4e.7):
    #   dispatcher.send(mesh_lists) -- launch one calc per (material,
    #     mesh) WITHOUT waiting (send_off, 13.5).
    #   dispatcher.next_rung() -- block until the next rung lands and
    #     return (material, result), where result is a {mesh, energy}
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
    active   = set(materials)
    in_air   = {}                  # rungs still in flight, per m
    opening  = set(materials)      # still in the opening (grid) phase

    # retire m with a verdict and drop it from the active set; a
    #   non-converged stop tags the mismatch (7.8 3d).
    function retire(m, verdict):
        outcomes[m] = verdict
        if verdict is NON_CONVERGED:
            on_non_converged(m)
        active.discard(m)

    # judge m's ladder (climb_next, 4e.3) and either retire it or
    #   launch its single next rung.  climb_next threads the per-
    #   material search state, so the bracket-refine phase persists
    #   across landings; the grid and unit-step climbs ignore it.
    function judge(m):
        (action, search[m]) = climb_next(rungs[m], search[m],
                                          configs[m])
        if action is CONVERGED:
            outcomes[m] = action.rung                # the Rung, not a
            active.discard(m)                        #   mismatch
        elif action is CEILING:
            retire(m, NON_CONVERGED)                 # 7.8 3d
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
                retire(m, NON_CONVERGED)             # run failure
            else:
                judge(m)
        else:
            # A continuation is exactly one rung.  If it failed to
            #   run the climb cannot advance, so stop the material
            #   rather than re-dispatch it forever (7.7); otherwise
            #   judge the extended ladder.
            if result is FAILED:
                retire(m, NON_CONVERGED)             # run failure
            else:
                judge(m)

    return outcomes, rungs
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


function build_mesh_unit(structure, options, mesh, id):
    # One explicit-mesh convergence unit (DESIGN 7.7 / 6.2.1).
    # `scfkp` is the makeinput key for an explicit axial-count mesh
    # (a style-code-1 k-point file); `kpt-mesh` is its calc-tag axis.
    # The cache identity is the same one the density units used
    # (6.2.1), so a mesh re-run in a later round is a cache hit and
    # costs nothing.
    unit_options = copy(options)
    unit_options["scfkp"] = mesh                 # [a, b, c]
    calc = buildCalcTag({ "kpt-mesh": encodeMeshValue(mesh) })
    return CalcUnit(id=id, calc=calc, structure=structure,
                    options=unit_options, wingbeat="imago",
                    key_fields=standardKeyFields(structure, options))


function predict_kpoint_density(structure, dataspace, system_type,
                                submodel, center):
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
    record = buildPredictionRecord(policy, density, confidence,
                                   under_trained, result, sig,
                                   system_type, submodel)
    return density, confidence, under_trained, record


function make_climb_dispatcher(structures, options_by_material,
                               workspace, parsl_config, executor,
                               force):
    # Build the dispatcher converge_by_climb (4e.5) drives, closing
    # over each material's structure and options, the workspace root,
    # the resolved Config (13.7), and the ONE shared executor every
    # send runs under (make_executor, 13.5): the pre-flight loen batch
    # and every climb rung run beneath the SAME executor, so the whole
    # run rides one warm pool (DESIGN 6.2.11) and lands in one tree.
    # `force` bypasses the run-reuse cache exactly as the pre-flight
    # dispatch does.  The material key doubles as the unit id
    # (materials ARE the reference ids the producer already uses).
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
                                       id = m)
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
        if entry.status != "done":
            return (m, FAILED)
        res = readResult(unit)                   # result.toml (6.1.2)
        assert res.kpoint_mesh == mesh           # honoured exactly
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
                            * tetraVol / hartree
                            * projArray(
                                permA, n, kIc)

    return pdosComplete
```

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

### 11.1 TOML Reader (DESIGN 5.2, 5.4)

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

```
function load(path, known_methods = None):
    raw = tomllib.load(path)

    # Rule 3 (top-level half): every required
    # top-level key must be present.  Check before
    # any value-level rule so the error message
    # names the missing field rather than failing
    # later with a value mismatch.
    for f in ("schema_version", "element_symbol",
              "nuclear_z", "nuclear_alpha",
              "covalent_radius"):
        require(f in raw, path,
            "missing top-level field: " + f)

    # Rule 1: schema version must equal 2.
    require(raw["schema_version"] == 2,
        path, "unsupported schema_version "
              + str(raw["schema_version"])
              + " (expected 2)")
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
to consider here -- the reader rejects any
`schema_version != 2` (DESIGN 5.2), and the producer
(11.4) writes every on-disk file as v2, so a loaded
database is always v2 and always carries a `default` tag.

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
        # Wrap the existing reduce algorithm: for
        # EVERY atom in the structure (not just one
        # element), compute a shell-code vector from
        # sub_spec's (level, thick, cutoff)
        # parameters.  Returned in site-index order
        # so 11.3.c can index by full-structure atom
        # index.  The species-pass filter on
        # method.element handles per-element
        # selection at the call site; computing all
        # atoms keeps the matcher contract uniform
        # across Python-side and loen-side families
        # (loen naturally writes one row per site of
        # the whole structure).
        return run_reduce_in_python(structure,
            sub_spec)

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
function buildInitialPotentials(manifest_path,
        force, single_element, dispatch_shape,
        partition, nodes, walltime, profile,
        save_config):
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
        # -> converg, shift -> kpshift -- and adds the imago_commit
        # cache identity.  Kept for the whole climb: every round's
        # mesh unit copies them (build_mesh_unit, 4e.7).
        options = make_producer_options(ref, imago_commit)
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
        density, confidence, under_trained, record = \
            predict_kpoint_density(
                struct, dataspace, ref.system_type,
                submodel, center = ref.kpoint_spec.density)
        seed_densities[ref.reference_id] = density

        # Everything the climb needs for THIS solid, gathered once
        # (build_climb_config below): the reciprocal-cell geometry
        # the rung mechanics read, the confidence-derived mode /
        # persistence / grid from the resolved policy, and the
        # energy / ceiling knobs.
        configs[ref.reference_id] = build_climb_config(
            ref, struct, confidence, under_trained,
            thresholds, max_count)

        # Store the record as a plain dict (metadata must be TOML-
        # serializable), and stamp the resolved per-atom k-point
        # flatness tolerance onto it: the guidance harvest reads
        # both, and the tolerance is a manifest/resolved fact
        # absent from any run's result.toml (DESIGN 7.8 / 5.7).
        predictions[ref.reference_id] = as_dict(record)
        predictions[ref.reference_id][
            "kpoint_convergence_threshold"] = (
                ref.kpoint_convergence_threshold)

        # Geometry-only fingerprint units: one structure-only
        # `-loen -scf no` unit per Fortran-side declaration, tagged
        # kind = "fingerprint" (DESIGN 6.2.9).  The bispectrum
        # fingerprint depends on geometry alone, so these need not
        # wait for a converged mesh; they dispatch in the pre-flight
        # below, and their run dirs persist for the harvest.
        loen_units.extend(build_loen_units(ref, struct, workspace))

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
    try:
        if loen_units is not empty:
            loen_flight = Flight(
                units        = loen_units,
                root         = workspace,
                parsl_config = parsl_config,
                sweep        = SweepRecord(
                    varied_axes = (), fixed_axes = {}))
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
        dispatcher = make_climb_dispatcher(
            struct_of, options_of, workspace,
            parsl_config = parsl_config, executor = executor,
            force = force)
        materials = [ref.reference_id
                     for ref in manifest.reference_solids]
        outcomes, rungs = converge_by_climb(
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
            id = ref.reference_id)
        converged_result = read_result_toml(workspace, converged)
        # The run log records the converged mesh AND its k-density
        # (both from record_converged) and the SCF iteration count.
        log.append(make_run_log_entry(
            ref, harvest, converged_result))

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
        # climb hands only the density and ladder here.  A converged
        # climb always carries at least the three distinct rungs the
        # stop test required (4e.3), so every converged solid
        # contributes an entry.
        entry = build_entry(
            workspace, struct, predictions[ref.reference_id],
            dataspace, load_structure(struct),
            ref.kpoint_convergence_threshold,
            harvest.grid_values, harvest.grid_energies,
            harvest.converged_kpoint_density,
            converged_result)
        save_entry(entry, "share/historicalGuidanceDB/")

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

    return ClimbConfig(
        classes           = classes,
        recip_mag         = recip_mag,
        recip_cell_volume = recip_cell_volume,
        mode              = policy.mode,
        flat_needed       = policy.flat_needed,
        grid_width        = policy.grid_width,
        start_offset      = policy.start_offset,
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
            require(f in ref or f in raw_defaults, path,
                "manifest rule 2: [[reference_solid "
                + ref.get("reference_id", "?")
                + "]] run setting " + f + " not resolvable"
                + " (absent here and from [defaults])")
        # The harvest setting kpoint_convergence_threshold is
        # EXEMPT from this resolvability rule: it carries a
        # built-in default (5e-4 eV/atom; DESIGN 5.7 / 7.8), so a
        # solid naming neither it nor a [harvest] block is
        # accepted -- apply_manifest_defaults supplies the default.

        rid = ref["reference_id"]

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


# The five run settings that may live in [defaults] and be
#   inherited per solid (DESIGN 5.7).  system_type is NOT among
#   them -- it is structure metadata, always named per solid.
RUN_SETTING_KEYS = ("basis", "functional",
    "kpoint_integration", "kpoint_spec", "scf_threshold")

# The producer's built-in k-point flatness tolerance, used when a
#   solid names neither its own kpoint_convergence_threshold nor a
#   [harvest] block (DESIGN 5.7 / 7.8).  Per atom, in eV.
DEFAULT_KPOINT_CONVERGENCE_THRESHOLD = 5.0e-4    # 0.5 meV/atom

# The adaptive-climb tuning knobs that may live in the optional
#   [harvest.kpoint_climb] sub-table (DESIGN 5.7 / 3.12.6).  Six
#   name the confidence-to-policy PolicyThresholds (4e.4); max_count
#   is the per-axis ceiling (4e.2).  Each carries a provisional
#   built-in default (mesh_climb), so the sub-table -- and any knob
#   -- may be omitted.  Database-wide: no per-solid override.
KPOINT_CLIMB_KEYS = ("confidence_high", "grid_width",
    "start_offset_moderate", "start_offset_cold",
    "flat_needed_confident", "flat_needed_cold", "max_count")


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
        # No network.
        return ref.structure_path

    # cod_id ref: fetch the pinned revision once to a
    # plain local location.  Strict on failure (network
    # down / COD outage / pinned revision missing) --
    # never falls back to another revision, because a
    # silent fallback would desync the build from the
    # pinned manifest (DESIGN 5.7).
    local = ("share/atomicBDB/cache/structures/"
             + ref.reference_id + cod_extension(ref))
    if not file_exists(local):
        fetch_cod_structure(
            cod_id       = ref.cod_id,
            cod_revision = ref.cod_revision,
            dest         = local)
    return local


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
    return report          # the CLI prints it plus a tally


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


function harvestFingerprints(flight, ref, env,
        result_toml, characterization):
    # Every environment harvests the database-wide
    # [characterization] preferred recipe (one sub_spec per
    # method, marked preferred=true), plus any rare per-
    # entry override the customization added (extra NON-preferred
    # sub_specs).  DESIGN 5.7.
    decls = []
    for fp in characterization:
        decls.append({"method": fp["method"],
                      "sub_spec": fp["sub_spec"],
                      "preferred": true})
    for fp in env.overrides:
        decls.append({"method": fp["method"],
                      "sub_spec": fp["sub_spec"],
                      "preferred": false})
    if decls is empty:
        return []

    # INTERIM (until C55/C58): the Fortran-side bispectrum
    # harvest is not built yet, so refuse any loen-side
    # declaration up front rather than silently dropping a
    # fingerprint the recipe asked for.  When C55/C58 land,
    # this guard is replaced by a per-declaration dispatch to
    # harvestLoenFingerprint (below) for every matcher whose
    # needs_loen_run is true.
    for d in decls:
        if MATCHERS[d["method"]]().needs_loen_run:
            raise NotImplementedError(
                "method needs a loen run; the Fortran-side "
                "harvest is C55/C58 and is not built yet")

    # Every remaining declaration is Python-side (reduce).
    # Read the run's EXPANDED full-cell structure
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
    # applies.
    max_cutoff = max(d["sub_spec"]["cutoff"]
                     for d in decls)
    structure = read_structure(
        result_toml.outputs["structure"])
    build_min_dist_matrix(structure, max_cutoff)

    # The expanded skeleton is ordered by the run's sorted
    # (dat) numbering, but atom_site is a skeleton index, so
    # map it to the structure row through datSkl.map (the same
    # map step i reads); the map yields both the row and that
    # row's element symbol.
    (dat_index, map_element) = skeleton_to_dat(
        result_toml.outputs["datSkl_map"])[env.atom_site]
    # Guard the numbering assumption: the structure row and
    # the map must name the same element, or the expansion and
    # the map have desynced and the fingerprint would describe
    # the wrong atom.  Strict refusal beats a silent mismatch.
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
        # In-process compute against the shared structure;
        # compute_query trims to this declaration's own
        # sub_spec cutoff.  Wrap the chosen vector via
        # matcher.build_payload so the per-matcher payload
        # field name (DESIGN 5.2: bispec uses `values`, reduce
        # uses `shell_code`) flows through the same accessor
        # the loen-side branch uses.
        vectors = matcher.compute_query(structure, sub_spec)
        payload = matcher.build_payload(vectors[dat_index])
        fingerprints.append(FingerprintRecord(
            method    = method,
            sub_spec  = sub_spec,
            preferred = d["preferred"],
            payload   = payload))
    return fingerprints


function harvestLoenFingerprint(flight, ref,
        atom_site, matcher, sub_spec):
    # FINISHED-STATE path (C55/C58): the per-declaration
    # dispatch the interim guard in harvestFingerprints stands
    # in for.  Read the fort.21 of the `-loen -scf no` unit
    # that kaleidoscope already dispatched for this
    # (solid, method, sub_spec) back in step 1b.  No loen run
    # happens here and there is no separate loen cache --
    # kaleidoscope's run-reuse cache (DESIGN 6.2.5) already
    # owns recompute avoidance.  The unit's run directory
    # follows the calc-tag convention (DESIGN 6.2.4): id =
    # reference_id, calc = "loen-<method>-<slug>"; the slug
    # encodes the sub_spec so two declarations differing in
    # any key or value land in different run directories by
    # construction.
    slug     = sub_spec_slug(sub_spec)
    calc_tag = "loen-" + matcher.name + "-" + slug
    run_dir  = unit_run_dir(flight.root,
        ref.reference_id, calc_tag)
    out_path = path_join(run_dir, "fort.21")

    rows = matcher.parse_loen_output(out_path,
        sub_spec)
    # atom_site is a skeleton index per the manifest contract;
    # like the Python-side branch it must be mapped to the
    # run's row numbering through datSkl.map before indexing
    # rows (left naive here -- wire the same skeleton_to_dat
    # step when C55/C58 land).  The matcher's build_payload
    # accessor wraps the vector in the per-matcher payload
    # shape (DESIGN 5.2: bispec uses `values`, reduce uses
    # `shell_code`) so producer and consumer stay symmetric on
    # field naming.
    return matcher.build_payload(
        rows[atom_site - 1])


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
    name   : str     # the staged copy's path, relative to the
                     #   run dir (what a prior run left there)
    source : str     # the current input byte-compared against
                     #   that staged copy

dataclass KeyFields:
    scalars : dict   # verbatim-compared identity fields,
                     #   e.g. {scf_threshold, imago_commit}
    files   : list   # KeyFile entries to byte-compare (name +
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
        # the "SCF after loen" the seed run hit.)  imago_commit is
        # a cache-only scalar, not an imago.OPTION_KEYS member.
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
        write_toml(join(wingbeat_dir, "result.toml"),
                   as_dict(result))

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
    if unit.prepared_dir is not None:
        commit_prepared_inputs(unit.prepared_dir, wingbeat_dir)
    else if not is_prepared(wingbeat_dir):
        mk_opts = { k : v for k, v in unit.options.items()
                    if k not in imago.OPTION_KEYS
                    and k not in CACHE_ONLY_KEYS }
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
    # Copy the driver-staged inputs (structure.dat, imago.dat,
    # scfV, kp files -- DESIGN 6.2.5) into the run dir so
    # run_prepared finds them.  The staging area is transient
    # (the prepare pass, 11.4, rebuilds it each producer run),
    # so the commit simply copies from it.
    for name in list_files(prepared_dir):
        copy_file(join(prepared_dir, name),
                  join(wingbeat_dir, name))
```

An ASE wingbeat (D12) and future adapters implement the
same protocol; the dispatch core (13.5) never changes
when one is added.  `commit_prepared_inputs` is
ImagoWingbeat's own step; another wingbeat stages its
inputs however its toolchain requires.

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
    # Key files: byte-compare each declared key file's current
    # source against the copy already staged in the run dir under
    # its name.  No hashing -- a developer can diff the files to
    # see why a cache missed (DESIGN 6.2.5 / 5.7).
    for key_file in unit.key_fields.files:
        staged = join(wingbeat_dir, key_file.name)
        if not exists(staged) \
           or not files_byte_equal(key_file.source, staged):
            return False
    return True
```

```
function write_cache_key(wingbeat_dir, unit):
    # The identity snapshot, written on launch (13.5).  Only the
    # key-file NAMES are recorded (for inspection); the byte
    # compare reads the staged files themselves, not this list.
    write_toml(join(wingbeat_dir, "cache_key.toml"),
        { scalars = unit.key_fields.scalars,
          files   = [kf.name for kf in unit.key_fields.files] })
```

Each `KeyFile` names both halves the compare needs -- the
`source` (the current input) and the `name` (the staged copy's
run-dir path) -- so the core stays oblivious to how a client's
inputs map onto staged files (DESIGN 6.2.5).  For the producer,
the driver's prepare step (11.4, Phase 1b) points the
`structure.dat` KeyFile's `source` at the staged copy it builds.

### 13.5 Dispatch driver (DESIGN 6.2.3)

One task per unit; per-future exception capture so a single
failure never aborts the flight (Principle 10).  Resuming a
flight is just re-running it: the hit-test skips the `done`
units and re-dispatches the rest -- unless the re-run switch
(`force`) is set, which bypasses the cache so every unit
re-launches (DESIGN 6.2.5).

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
    outstanding = []                   # list of (unit, future)
    for unit in units:
        outstanding.append(
            (unit, dispatch_unit(flight, unit, executor, force)))
    return outstanding
```

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
function dispatch(flight, executor=None, force=False):
    # The one-shot convenience form: send every unit off, then
    # collect them all in unit order (DESIGN 6.2.3).  Behaviour is
    # identical to the pre-split driver, so every existing caller
    # is unchanged; the climb (4e.5) uses send_off + collect_next
    # directly instead.  We tear the executor down at the end only
    # if we built it here (a caller-supplied executor is the
    # caller's to close).
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
    # and hand the unit to the executor.
    makedirs(wingbeat_dir, exist_ok=True)
    write_cache_key(wingbeat_dir, unit)          # 13.4
    write_status(wingbeat_dir, id=unit.id, calc=unit.calc,
        status="queued",
        wingbeat=(unit.wingbeat or flight.default_wingbeat),
        submitted_at=now())
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
# rest sizes a job, and there are exactly two job CLASSES: a
# worker (one calculation, the per-worker keys) and an
# orchestrator (a driver that prepares units and fans them out,
# the grouped block at the end).  The file grows by job class,
# never by builder (ARCHITECTURE 9.4).
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
    }
```

```
function merge_settings(base, overlay):
    # The ONE merge every overlay uses -- profile, queue, and the
    # per-run flags alike (DESIGN 6.2.11, decision 1).  Per key, and
    # one level down: when a setting is itself a BLOCK of settings
    # (`orchestrator` is the only one today), the overlay names only
    # the keys it means to change and the rest keep the value the
    # layer beneath gave them.
    #
    # Replacing the whole block instead would silently discard facts
    # the curator never mentioned -- "the driver needs 2G on the
    # debug queue" would also drop its cores and walltime, and they
    # would reappear as plausible-looking fallbacks rather than as
    # an error.  A block holds plain values, never further blocks,
    # so one level is the whole of the descent.
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

    lines = ["#!/bin/bash",
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
    settings = _starter_schema()
    for key in ["partitions", "cores_per_node",
                "memory_per_node", "gpus_per_node"]:
        if key in facts: settings[key] = facts[key]
    return render_full_dict(settings, notes = _SETTINGS,
                            options = facts, account_hint = facts,
                            blanks = ["partitions", "worker_init"])
```

## 14. makeinput Callable Build API (DESIGN 6.3)

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
    return Verification(
        grid_values            = tuple(grid),
        grid_energies          = (tuple(energies)
                                  if energies is not None
                                  else None),
        converged_at           = v["converged_at"],
        converged_mesh         = (tuple(mesh)
                                  if mesh is not None else None),
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
#   (DESIGN 6.2.1/6.2.10): `converg` (the SCF convergence limit,
#   a makeinput option -- the concrete name for DESIGN's
#   "scf_threshold") and `imago_commit` (the build identity,
#   producer-injected).  Taken from a unit's options when present.
KEY_SCALAR_NAMES = ("converg", "imago_commit")


function standard_key_fields(structure, options):
    # DESIGN 6.2.5: the producer's cache identity -- the scalars
    # taken from `options` (the SCF threshold and imago_commit)
    # plus one key file, `structure.dat`, byte-compared.  The key
    # file is makeinput's OUTPUT, not the raw skeleton: it bakes
    # in every input that changes the result (the type/species
    # assignment, basis, functional, potential), so any of those
    # changing misses the cache on its own, with no hand-listed
    # "options that matter" to fall stale.  The KeyFile `source`
    # is provisional here (the skeleton `structure`); the driver's
    # prepare step (11.4, Phase 1b) re-points it at the built
    # structure.dat once that file exists.
    return KeyFields(
        scalars = { name : options[name]
                    for name in KEY_SCALAR_NAMES
                    if name in options },
        files   = [KeyFile(name = "structure.dat",
                           source = structure)])
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
recorded alongside for the conversion.  `imago_commit` falls back
to `"unknown"` when the producer injected none.  `spin_polarization`
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
        kpoint_densities, energies, meshes, rts = [], [], [], []
        for u in grid:
            rt = read_result_toml(
                join(workspace_root, "wingbeats", u.id, *u.calc,
                     "result.toml"))
            kpoint_densities.append(swept_value_of(u, axis))
            energies.append(rt["total_energy"])
            meshes.append(rt.get("kpoint_mesh"))
            rts.append(rt)

        # c. A single-point grid harvests deliverables but does
        #    NOT auto-stage a guidance entry (DESIGN 6.2.1 / 7.7):
        #    one converged calc is weaker evidence than a grid.
        #    This covers both trust_no_verify and a single-point
        #    curator_override, and MUST precede pick_converged:
        #    the two-sided convergence test below needs >= 3
        #    points and would otherwise misreport one as "energy
        #    still moving."
        if len(grid) == 1:
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
        entry = build_entry(
            workspace_root, grid[0].structure, prediction,
            dataspace, sc, kpoint_threshold,
            c_kpoint_densities, c_energies,
            kpoint_densities[idx], rts[idx])

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
                     chosen_result):
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
            imago_commit     = chosen_result.get("imago_commit")
                               or "unknown",
            curator          = "guidance_harvest.py"))


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
function promote(db_root, mode):
    for system_type in VALID_SYSTEM_TYPES:
        staged = sorted(glob(
            join(db_root, "staging", system_type), "*.toml"))
        for path in staged:
            entry = load_entry(path, system_type, {})

            if mode == "dry-run":
                print_would_promote(entry,
                    auto_promote_ok(entry))
            elif mode == "all":
                move_to_entries(path, db_root, system_type)
            elif mode == "auto-promote":
                if auto_promote_ok(entry):
                    move_to_entries(path, db_root, system_type)
                # else: leave in staging for review.
            else:                                # interactive
                print_summary(entry)             # sig+measured
                                                 #   +verif+prov
                choice = ask("PROMOTE / SKIP / DELETE")
                if choice == "PROMOTE":
                    move_to_entries(path, db_root, system_type)
                elif choice == "DELETE":
                    remove_file(path)
                # SKIP: leave in staging.
```

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

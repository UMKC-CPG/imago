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
(section 4b below) from the Cartesian xyz operations
that arrive on disk.  Atom Cartesian positions are
converted to fractional using invRealVectors (=
recipVectors / 2*pi) of the same loaded lattice.

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
    for A = 1 to numAtomSites:
        for i = 1 to 3:
            abcAtomPos(i, A) =
                sum(invRealVectors(i,:)
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
                             realVectors,
                             invRealVectors):
    # Conjugate conv-abc operations into the basis of
    # the loaded real lattice for use by buildAtomPerm.
    #
    #   C          = M_loaded^{-1} * M_conv
    #              = invRealVectors^T * convLattice
    #   R_real_abc = C * R_conv_abc * C^{-1}
    #   t_real_abc = C * t_conv_abc
    #
    # invRealVectors^T equals realVectors^{-1} by the
    # orthogonality identity
    #   realVectors^T * recipVectors = 2*pi*I,
    # so the formula uses pre-computed O_Lattice matrices
    # and never re-inverts realVectors at runtime.
    allocate abcRealPointOps(3, 3, numPointOps)
    allocate abcRealFracTrans(3, numPointOps)

    if cellMode == 'full':
        # Identity shortcut: M_loaded == M_conv, so
        # C = I and the conjugation collapses to a copy
        # from the on-disk arrays into the runtime
        # arrays.  Saves the inverse and the per-op
        # similarity kernel for the most common mode.
        for i = 1 to numPointOps:
            abcRealPointOps(:,:,i) =
                convAbcPointOps(:,:,i)
            abcRealFracTrans(:,i) =
                convAbcFracTrans(:,i)
    else:
        # Full conjugation path: form C once, then
        # apply the similarity transform per operation.
        C     = invRealVectors^T * convLattice
        C_inv = inverse_3x3(C)
        for i = 1 to numPointOps:
            abcRealPointOps(:,:,i) =
                C * convAbcPointOps(:,:,i) * C_inv
            abcRealFracTrans(:,i) =
                C * convAbcFracTrans(:,i)

    return (abcRealPointOps, abcRealFracTrans)
```

```
function computeRecipPointOps(numPointOps,
                              convAbcPointOps,
                              cellMode, convLattice,
                              realVectors,
                              invRealVectors):
    # Conjugate conv-abc operations into the basis of
    # the loaded reciprocal lattice for use by
    # initializeKPointMesh (IBZ folding).
    #
    # Point group rotations transform identically under
    # the dual abc bases, so the similarity uses the
    # same change-of-basis matrix C as the real-space
    # sibling.  No translation field -- reciprocal-space
    # operations have no translation.
    allocate abcRecipPointOps(3, 3, numPointOps)

    if cellMode == 'full':
        # Same identity shortcut as the real-space side.
        for i = 1 to numPointOps:
            abcRecipPointOps(:,:,i) =
                convAbcPointOps(:,:,i)
    else:
        C     = invRealVectors^T * convLattice
        C_inv = inverse_3x3(C)
        for i = 1 to numPointOps:
            abcRecipPointOps(:,:,i) =
                C * convAbcPointOps(:,:,i) * C_inv

    return abcRecipPointOps
```

Both routines run unconditionally in every style-code
branch (style 0 sets up trivial identity operations as
before; styles 1 and 2 read real symmetry from the kp
file).  The `cellMode` flag selects the identity
shortcut versus the full conjugation path -- no other
`full`-vs-`prim` branching exists outside these two
routines.

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
`makeinput.py` (11.3), the regeneration pipeline
(11.4), and the validation harness (11.5). All live in
Python under `src/scripts/`. The Fortran side does not
change.

### 11.1 TOML Reader (DESIGN 5.2, 5.4)

Parses a per-element `s_gaussian_pot.toml` file and
applies the nine validation rules from DESIGN 5.2
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

        # Length consistency (rule 4)
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

        # Fingerprint sub-blocks (rules 8 and 9).
        # Each [[potential.fingerprint]] record
        # contributes one FingerprintRecord.  The
        # method/sub_spec pair is unique per entry
        # (rule 8); the method must be known if a
        # registry was supplied (rule 9).
        fingerprints = []
        seen_method_subspec = set()
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

            # Payload = all keys other than method
            # and sub_spec.  Matchers validate their
            # own payload shape at lookup time.
            payload = {k: v for (k, v) in fp_dict
                       if k not in ("method",
                                    "sub_spec")}
            fingerprints.append(FingerprintRecord(
                method   = method,
                sub_spec = sub_spec,
                payload  = payload))

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
                  "convergence_threshold",
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
themselves.  Any real diff between two regenerations
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

    # method and sub_spec come first, in that order;
    # the payload keys follow in the payload dict's
    # iteration order.  Width alignment spans
    # method/sub_spec plus the payload's scalar and
    # multi-line keys, so all `=` signs align within
    # the block.
    fixed_keys   = ["method", "sub_spec"]
    payload_keys = list(fp.payload.keys())
    align_keys   = fixed_keys + payload_keys
    width        = max(len(k) for k in align_keys)

    out.append(format_kv("method", fp.method, width))
    out.append(format_kv(
        "sub_spec",
        format_inline_table(fp.sub_spec),
        width))

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
              "kpoint_spec", "convergence_threshold",
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

Most runs exercise only a reduced subset of this flow.
Sub-section 11.3.0 pins that subset -- the path active
when no environment-based matcher (`-reduce` / `-bispec`)
is selected.  It is the path the first consumer milestone
(C47) implements before C54 onward layer in the matcher
machinery; reading it first makes the seven full-flow
sub-sections easier to place.

The mapping from DESIGN to PSEUDOCODE here is one-to-
one:

  PSEUDOCODE   DESIGN
  -----------------------------------------------------
  11.3.a       8.9 + 5.10.5 (matcher protocol and
               LOEN parameter contract)
  11.3.b       5.6.3        (per-element preflight,
                             including the rule-4
                             coverage check)
  11.3.c       5.6.4        (species pass and scope
                             resolution)
  11.3.d       5.6.5        (manifest-entry pick per
                             species)
  11.3.e       5.6.6        (type pass and electronic-
                             state perturbation)
  11.3.f       5.10         (nested-makeinput
                             bootstrap)
  11.3.g       5.6.7        (driver and on-the-wire
                             emit)

---

#### 11.3.0 Reduced Flow (no environment matcher; DESIGN 5.6)

The seven sub-sections above specify the full Phase-2
selection flow.  Most runs -- and the first consumer
milestone (C47) -- travel only the reduced path taken
when the CLI selects no environment-based matcher
(`-reduce` / `-bispec`).  This sub-section names which
branches are live in that path so the reduced consumer
has a self-contained spec, separate from the matcher
machinery that C54 and onward layer in.

When `first_environment_matcher` (11.3.g) returns None:

  - The nested-makeinput bootstrap (11.3.f) never fires;
    it is gated on a matcher whose `needs_loen_run` is
    true, and here there is no matcher at all.
  - The preflight (11.3.b) still loads each element's
    database and still marks missing ones for the legacy
    path, but skips `require_coverage` -- that check is
    meaningful only while a matcher is active.
  - The species pass (11.3.c) returns an empty
    `env_species_ids` set; only the position-based flags
    (`-target`, `-block`, `-xanes`) and the default
    one-species-per-element grouping remain.
  - The entry pick (11.3.d) loses its precedence-2
    fingerprint branch: that branch is gated on both
    `active_matcher is not None` and `species_id in
    env_species_ids`, and the latter set is now empty.
    Only precedence 1 (`-pot LABEL`) and precedence 3
    (`default_entry`) survive.

An element whose augmented file is absent is handled the
same way as in the full flow: the preflight marks
`databases[elem] = None` and the driver (11.3.g, steps 5
and 7) emits it via the legacy `pot1`/`coeff1` path, with
no library entry consulted.  There is no schema-v1 case
to consider here -- the reader rejects any
`schema_version != 2` (DESIGN 5.2), and the producer
(11.4) regenerates every on-disk file as v2, so a loaded
database is always v2 and always carries a `default` tag.

The reduced entry pick is the `active_matcher = None`
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
    # no active environment matcher there is no
    # precedence-2 fingerprint match to attempt first.
    return default_entry(db)
```

**Carry-forward to the full flow.**  When C54 onward turn
the matcher machinery on, this reduced pick is subsumed by
the full 11.3.d: precedence 2 slots back in between the
override and the default, gated on the now-non-empty
`env_species_ids`.  No branch written for the reduced flow
is rewritten -- the full flow only inserts the middle
precedence.

---

#### 11.3.a Matcher Protocol (DESIGN 5.6, 5.10.5; ARCH 8.9)

Each matcher knows one descriptor family.  The species
pass calls `compute_query` and `distance` to bucket
atoms (11.3.c); the entry pick calls `representative`
and `distance` against per-entry fingerprints (11.3.d);
the producer (11.4) and the consumer bootstrap (11.3.f)
both call `to_loen_input` and `parse_loen_output` on
matchers whose `needs_loen_run` is true.  The protocol
isolates Imago's Fortran side from manifest-schema
growth: a new descriptor family is a new class plus a
new `MATCHERS` registry entry.

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
        # Python-side matchers compute in-process;
        # loen-side matchers fire the bootstrap
        # (11.3.f) and parse its fort.21, which
        # naturally produces one row per potential
        # site of the whole structure.
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

    function build_payload(query_vector):
        return {"shell_code": query_vector}

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

    function compute_query(structure, sub_spec):
        # Outer entry point.  Triggers the bootstrap
        # of 11.3.f on `structure` with this matcher's
        # loen-side parameters, then parses fort.21
        # via parse_loen_output.  Element-aware mode
        # is deferred to TODO C62 / D10.
        if sub_spec.get("by_element", False):
            error("element-aware bispectrum"
                + " (sub_spec by_element=true) is"
                + " not yet implemented; see TODO"
                + " D10 / C62")
        fort21 = runLoenBootstrap(structure, this,
            sub_spec)
        return this.parse_loen_output(fort21,
            sub_spec)

    function distance(a, b):
        # Euclidean distance between bispectrum
        # vectors of length 2*twoj2 + 1.  Symmetric,
        # cheap to compute, and consistent with the
        # element-wise mean used by `representative`
        # below.
        return l2_norm(vector_subtract(a, b))

    function representative(members):
        # Element-wise arithmetic mean of the member
        # vectors.  All members share the same length
        # (2*twoj2 + 1) because they come from one
        # bootstrap run with one sub_spec, so the
        # mean is well-defined slot by slot.
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
            "max_neigh"    : sub_spec.get(
                "max_neigh", 20),
            "cutoff"       : sub_spec.get(
                "cutoff", 5.0),
            "angleSqueeze" : sub_spec.get(
                "angle_squeeze", 0.85)}

    function parse_loen_output(path, sub_spec):
        # fort.21: one row per potential site of the
        # whole structure; 2*twoj2 + 1 real values
        # per row, followed by a sum column the
        # matcher ignores.  Returns a list of
        # vectors in site-index order, length =
        # structure atom count.
        n_slots = 2 * sub_spec["twoj2"] + 1
        rows = read_text_rows(path)
        return [row_first_n_floats(r, n_slots)
                for r in rows]

    function extract_query_vector(payload):
        # Bispectrum records carry the vector in the
        # `values` field (DESIGN 5.2 / 5.3 sketch).
        return payload["values"]

    function build_payload(query_vector):
        return {"values": query_vector}


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
fallback path, and confirms (when an environment-based
scheme is active) that every element's database covers
the requested `(method, sub_spec)`.  Failing fast here
keeps the expensive bootstrap (11.3.f) out of the
"nothing to harvest" case.

```
function perElementPreflight(structure,
        active_matcher):
    # `active_matcher` is None when no environment-
    # based scheme is in play.  Returns a mapping
    # elem -> ElementDatabase | None; None marks the
    # legacy-fallback path that 11.3.g consumes.
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

        # Coverage check (5.6.3 step 4).  Only
        # meaningful when an environment-based
        # matcher is active; with only spatial flags
        # or a manual -pot override, the default-tag
        # fallback in 11.3.d step 3 always succeeds
        # regardless of fingerprint records.
        if active_matcher is not None:
            require_coverage(databases[elem],
                active_matcher, elem, path)

    return databases


function require_coverage(db, matcher, elem, path):
    # At least one entry in `db` must carry a
    # FingerprintRecord whose (method, sub_spec)
    # matches the active matcher.  We use the
    # library's find_fingerprint to keep the sub_spec
    # comparison canonical (matching rule 8's
    # equality semantics).
    method   = matcher.name
    sub_spec = matcher.active_sub_spec
    for entry in db.potentials:
        try:
            find_fingerprint(entry, method, sub_spec)
            return       # coverage exists; done
        except KeyError:
            continue
    error("preflight coverage check: no entry in "
        + path + " carries a fingerprint for ("
        + method + ", " + str(sub_spec) + ").  Add"
        + " the [[reference_solid.entry.fingerprint]]"
        + " declaration to the curation manifest and"
        + " regenerate the database (DESIGN 5.7).")
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

```
function speciesPass(structure, settings, databases,
        atom_fingerprints):
    # `atom_fingerprints` is None when no environment-
    # based matcher is active.  When non-None it
    # holds one vector per atom of the whole
    # structure (in site-index order); the bucketing
    # step indexes into it by full-structure index.
    #
    # The third return value, env_species_ids, is the
    # set of species IDs that were produced by an
    # environment-based scheme.  Per DESIGN 5.6.5
    # step 2, fingerprint matching applies only to
    # those species; the entry-pick pass (11.3.d)
    # gates on this set rather than on the presence
    # of fingerprint vectors.
    n_atoms          = len(structure.atoms)
    atom_species_id  = [1] * n_atoms
    named_regions    = {}     # name -> atom-index set
    env_species_ids  = set()  # species made by an
                              # environment scheme

    for method in settings.methods:
        if method.op == "spatial":
            # Position-based (-target, -block).
            # Existing geometric grouping; the new
            # `name=` keyword is consumed here.
            # Spatially-grouped species do NOT join
            # env_species_ids -- DESIGN 5.6.5 step 2
            # forbids fingerprint matching on them.
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
            # loop.  The new species IDs the
            # bucketing produces are recorded in
            # env_species_ids so 11.3.d knows to
            # attempt fingerprint matching on them.
            scope = resolve_scope(method.scope,
                named_regions, n_atoms)
            atoms_in_scope = \
                atoms_of_element_in_scope(structure,
                    method.element, scope)
            new_ids = assign_species_by_bucketing(
                atom_species_id,
                atoms_in_scope,
                atom_fingerprints,
                method.matcher)
            env_species_ids.update(new_ids)

        elif method.op == "electronic":
            # -xanes etc. defer to the type pass
            # (11.3.e); no species-level effect.
            continue

    return (atom_species_id, named_regions,
            env_species_ids)


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

    # Assign a fresh species ID to each bucket.
    # Atoms outside `atoms_in_scope` keep whatever
    # species ID earlier flags produced.  Return the
    # set of newly-minted species IDs so the caller
    # can record them in env_species_ids (consumed
    # by 11.3.d to gate the fingerprint-match step).
    new_species_ids = set()
    next_id = max(atom_species_id) + 1
    for b in buckets:
        for i in b.atom_indices:
            atom_species_id[i] = next_id
        new_species_ids.add(next_id)
        next_id += 1
    return new_species_ids
```

---

#### 11.3.d Manifest-Entry Pick per Species (DESIGN 5.6.5)

For each `(element, species)` pair, chooses exactly
one `PotentialEntry` from the element's database.
Three-step precedence: `-pot LABEL` manual override,
fingerprint match against entry-attached records, then
default-tag fallback.  The default fallback is the
single point that always succeeds, guaranteed by
per-element-database rule 7.

```
function pickManifestEntry(species_id, species_atoms,
        element, db, active_matcher, pot_override,
        atom_fingerprints, env_species_ids):
    # species_id        -- the species being picked
    # species_atoms     -- atom indices in this
    #                      (element, species) bucket
    # db                -- the element's
    #                      ElementDatabase
    # active_matcher    -- the active environment-
    #                      based matcher, or None
    # pot_override      -- the -pot LABEL value, or
    #                      None
    # atom_fingerprints -- per-atom fingerprint
    #                      vectors (length = N), or
    #                      None
    # env_species_ids   -- set of species IDs
    #                      produced by an environment-
    #                      based scheme.  Per DESIGN
    #                      5.6.5 step 2, fingerprint
    #                      matching applies only to
    #                      members of this set.

    # Precedence 1: manual override.
    # KeyError here is fatal; -pot is a deliberate
    # user choice and silently falling back would
    # mask the intent.
    if pot_override is not None:
        try:
            return lookup(db, pot_override)
        except KeyError:
            error("-pot " + pot_override
                + " not found in "
                + "share/atomicPDB/" + lower(element)
                + "/s_gaussian_pot.toml; manual"
                + " override must match an existing"
                + " label")

    # Precedence 2: fingerprint match.  Gated on
    # species_id in env_species_ids per DESIGN 5.6.5
    # step 2 -- spatially-grouped species (-target,
    # -block) skip this branch even when an
    # environment-based matcher is active elsewhere
    # in the run.
    if (species_id in env_species_ids
            and active_matcher is not None
            and atom_fingerprints is not None):
        member_vectors = [atom_fingerprints[i]
                          for i in species_atoms]
        rep = active_matcher.representative(
            member_vectors)
        best_entry    = None
        best_distance = +infinity
        method   = active_matcher.name
        sub_spec = active_matcher.active_sub_spec
        for entry in db.potentials:
            try:
                fp = find_fingerprint(entry, method,
                    sub_spec)
            except KeyError:
                continue      # entry has no matching
                              # fingerprint; skip
            # The payload's vector field name is
            # matcher-specific (DESIGN 5.2: bispec
            # uses `values`, reduce uses
            # `shell_code`); the matcher's
            # extract_query_vector accessor knows
            # which key to read so this code stays
            # descriptor-agnostic.
            entry_vector = \
                active_matcher.extract_query_vector(
                    fp.payload)
            d = active_matcher.distance(rep,
                entry_vector)
            if d < best_distance:
                best_distance = d
                best_entry    = entry
        if (best_entry is not None
                and best_distance
                    <= active_matcher
                       .default_similarity_floor):
            return best_entry
        if best_entry is not None:
            warn("species in " + element
                + " has best fingerprint match "
                + best_entry.label + " at distance "
                + str(best_distance) + " (> floor "
                + str(active_matcher
                      .default_similarity_floor)
                + "); falling back to default tag")
        # If best_entry is None, the preflight
        # already confirmed at least one entry has
        # the fingerprint for this element -- the
        # per-species absence here is normal when
        # the curator added the fingerprint only to
        # the most-relevant entries.  Fall through.

    # Precedence 3: default tag.  Guaranteed to
    # succeed by per-element-database rule 7.
    return default_entry(db)
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

#### 11.3.f Nested-makeinput Bootstrap (DESIGN 5.10)

Triggered only when the active matcher's
`needs_loen_run` is true AND the preflight (11.3.b) has
already confirmed coverage.  Spawns a stripped-down
nested `makeinput.py` to build a throwaway `imago.dat`,
runs `imago.py -loen -scf no` against it, and parses
the resulting `fort.21`.  The recursion guard --
explicit `--no-loen-bootstrap` flag plus an argument
list that contains no environment-based grouping flag
-- is the sole mechanism preventing the nested call
from triggering its own bootstrap.

```
function runLoenBootstrap(structure, matcher,
        sub_spec):
    # Recursion guard: the running invocation must
    # not itself be a nested bootstrap call.  If it
    # is, the trigger logic in 11.3.g misfired; abort
    # with a clear message rather than silently
    # recursing.
    require(not settings.no_loen_bootstrap,
        "runLoenBootstrap called inside a"
        + " --no-loen-bootstrap invocation; the"
        + " trigger condition in 11.3.g (only fires"
        + " when matcher.needs_loen_run AND no"
        + " --no-loen-bootstrap is set) misfired")

    # 1. Scratch directory.  Re-create from scratch
    #    every run so stale fort.21 contents cannot
    #    be silently reused.
    scratch = ".inputTemp/loen_bootstrap/"
    fresh_directory(scratch)

    # 2. Build a minimal imago.dat via a stripped-
    #    down nested makeinput.py.  The argument list
    #    deliberately omits:
    #      * -pot LABEL  (per-element preflight in
    #        the nested call falls through to the
    #        default-tagged entry on the no-scheme
    #        branch -- exactly the path needed; no
    #        per-element -pot machinery on the CLI).
    #      * any grouping flag (-target, -block,
    #        -reduce, -bispec, -xanes).  Every atom
    #        in the nested call becomes its own
    #        species and type, which is fine for
    #        loen because loen iterates over
    #        potential sites regardless of grouping.
    #      * a non-trivial k-point mesh.  Γ-only is
    #        sufficient; loen reads the k-point
    #        block but does not use it for
    #        fingerprint output.
    nested_argv = [
        "makeinput.py",
        structure.skeleton_path,
        "-kpd",    "1",            # Γ-only
        "-scfkpd", "1",
        "--no-loen-bootstrap",     # recursion guard
        "--scratch",     scratch,
        "--loen-params",
            serialize_loen_params(
                matcher.to_loen_input(sub_spec))]
    run_subprocess(nested_argv,
        cwd   = current_working_directory(),
        check = True)

    # 3. Run loen against the freshly produced
    #    imago.dat.  -scf no skips the SCF run; loen
    #    only needs the structure and the LOEN_INPUT_
    #    DATA block the nested call just wrote.
    imago_dat = path_join(scratch, "imago.dat")
    fort21    = path_join(scratch, "fort.21")
    run_subprocess(
        ["imago.py", "-loen", "-scf", "no",
         "--input",   imago_dat,
         "--workdir", scratch],
        check = True)

    # 4. Hand the fort.21 path back to the caller.
    #    The matcher's parse_loen_output (11.3.a)
    #    knows the row width and column semantics;
    #    the bootstrap is descriptor-agnostic.
    return fort21


function serialize_loen_params(param_dict):
    # The nested makeinput call writes the
    # LOEN_INPUT_DATA block from whatever parameter
    # dict it receives via --loen-params.  Today the
    # block is emitted at lines 4659-4665 of
    # makeinput.py with hardcoded defaults; TODO C58
    # wires the dict produced by
    # matcher.to_loen_input(sub_spec) through that
    # emission instead.  JSON is sufficient as the
    # wire format.
    return json_dumps(param_dict)
```

---

#### 11.3.g Driver (DESIGN 5.6.7)

Top-level orchestrator.  Chains preflight, optional
bootstrap, species pass, manifest-entry pick, type
pass, and emit.  The driver is matcher-agnostic; all
descriptor-family knowledge lives in the matcher
classes (11.3.a) and the `MATCHERS` registry.

```
function emitInitialPotentials(structure, settings,
        imago_input):
    # 1. Identify the active environment-based
    #    matcher.  DESIGN 5.6.2 mutual exclusion is
    #    enforced at argparse time, so at most one
    #    matcher is active here -- often none, when
    #    only spatial flags and/or -pot are used.
    active_matcher = first_environment_matcher(
        settings.methods)

    # 2. Per-element preflight.  Loads each element's
    #    database, runs the coverage check when a
    #    matcher is active, and marks elements
    #    without a database for the legacy fallback.
    databases = perElementPreflight(structure,
        active_matcher)

    # 3. Compute per-atom fingerprints once (if
    #    needed) so both the species-pass bucketing
    #    and the entry-pick representative comparison
    #    see the same vectors.  For loen-side
    #    matchers, compute_query fires the bootstrap
    #    of 11.3.f and parses its fort.21.  For
    #    Python-side matchers, the call stays in
    #    process.
    atom_fingerprints = None
    if active_matcher is not None:
        atom_fingerprints = \
            active_matcher.compute_query(structure,
                active_matcher.active_sub_spec)

    # 4. Species pass.  Position-based and
    #    environment-based flags compose in CLI
    #    order; output is a per-atom species ID
    #    array, the dict of named regions used in
    #    scope resolution, and the set of species
    #    IDs produced by an environment-based scheme
    #    (consumed by 11.3.d to gate fingerprint
    #    matching per DESIGN 5.6.5 step 2).
    atom_species_id, named_regions, env_species_ids \
        = speciesPass(structure, settings, databases,
            atom_fingerprints)

    # 5. Manifest-entry pick per species.  Atoms in
    #    elements with no database file
    #    (databases[elem] is None) bypass the
    #    manifest machinery and emit via the legacy
    #    pot1/coeff1 path in step 7.
    species_potentials = {}
    for species in unique_species(atom_species_id):
        atoms   = atoms_of_species(species,
            atom_species_id)
        element = structure.atoms[atoms[0]].element
        if databases[element] is None:
            species_potentials[species] = \
                LegacyEntry(element)    # marker
                                        # consumed in
                                        # step 7
            continue
        species_potentials[species] = \
            pickManifestEntry(
                species_id        = species,
                species_atoms     = atoms,
                element           = element,
                db                = databases[element],
                active_matcher    = active_matcher,
                pot_override      = settings
                                     .pot_override,
                atom_fingerprints =
                    atom_fingerprints,
                env_species_ids   = env_species_ids)

    # 6. Type pass.  Inherit from the species pass
    #    and apply electronic-state flags (XANES
    #    today; future flags layer in unchanged).
    atom_type_id, type_potential = typePass(
        structure, atom_species_id,
        species_potentials, settings)

    # 7. Emit per-Imago-type blocks in today's on-
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

### 11.4 Regeneration Pipeline (DESIGN 5.7, schema v2)

Reproducible rebuild of all affected
`s_gaussian_pot.toml` files (see DESIGN 5.7 for the
layered reproducibility contract: bit-level emitter,
precision-level numerics, free metadata).  Step 1
always rebuilds the "isolated" entries from current
`pot1`/`coeff1` files (so atomSCF changes propagate);
step 2 layers curated solid-state entries on top from
each cached or freshly-run reference SCF, with one
`FingerprintRecord` harvested per declared
`[[reference_solid.entry.fingerprint]]`.

The v2 manifest reader (`load_manifest_v2` below)
enforces the nine validation rules of DESIGN 5.7 --
the v1-era three-rule reader has been retired with
the rest of schema v1.  Producer-side fingerprint
harvest splits Python-side matchers (in-process,
descriptor-agnostic) from Fortran-side matchers
(cached `imago.py -loen -scf no` runs producing
`loen-<m>-<s>.out` files alongside the SCF cache).

```
function buildInitialPotentials(manifest_path,
        force, single_element):
    manifest     = load_manifest_v2(manifest_path)
    imago_commit = git_sha("HEAD")
    timestamp    = iso8601_now_utc()

    # Step 1: refresh "isolated" entries.  atomSCF
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
        # Drop and rebuild the isolated entry so
        # atomSCF refreshes always propagate.
        db.potentials = [e for e in db.potentials
                         if e.label != "isolated"]
        db.potentials.insert(0,
            build_isolated_entry(elem, imago_commit,
                timestamp, manifest))
        databases[elem] = db

    # Step 2: process each reference solid.
    log = []
    for ref in manifest.reference_solids:
        if not force and is_cached_v2(ref,
                imago_commit):
            scf = load_cache(ref)
        else:
            # Materialize the structure (COD fetch
            # on cache miss for cod_id refs; disk
            # read for structure_path refs).  Cache
            # the bytes before running SCF so the
            # cache directory is self-describing.
            structure_bytes = materialize_structure(
                ref)
            write_cache_structure(ref,
                structure_bytes)
            scf = run_imago_scf(ref, structure_bytes)
            write_cache_inputs(ref, imago_commit)
            store_cache(ref, scf)
        log.append(make_run_log_entry(ref, scf))

        for spec in ref.entries:
            elem = spec.element
            if elem not in databases:
                continue   # filtered out by --element
            coeffs, alphas = extract_potential(scf,
                spec.atom_site)
            # Harvest one FingerprintRecord per
            # declared [[reference_solid.entry.
            # fingerprint]] declaration.  Python-side
            # matchers compute in-process; Fortran-
            # side matchers run cached loen jobs.
            fingerprints = harvestFingerprints(ref,
                spec, imago_commit)
            new = PotentialEntry(
                label         = spec.label,
                default       = spec.default,
                description   = spec.description,
                num_gaussians = len(coeffs),
                alpha_min     = min(alphas),
                alpha_max     = max(alphas),
                coefficients  = coeffs,
                alphas        = alphas,
                provenance    = make_imago_provenance(
                    imago_commit, timestamp,
                    ref, spec.atom_site,
                    scf.iterations),
                fingerprints  = fingerprints)
            db = databases[elem]
            # Replace any prior entry with the same
            # label.  Manifest rule 6 enforces
            # cross-solid (element, label)
            # uniqueness, so the only possible prior
            # is the same entry from a previous run.
            db.potentials = [e for e in db.potentials
                             if e.label != spec.label]
            db.potentials.append(new)

    # Step 3: write all affected element files via
    # the deterministic emitter (5.5).
    for elem, db in databases.items():
        save(db, element_path(elem))

    # Step 4: run log capturing the manifest
    # snapshot, per-run iteration counts, and the
    # Imago commit.  The validation harness (11.5)
    # reads this log.
    write_run_log(
        "share/curation/run_log.toml",
        manifest_snapshot = manifest,
        imago_commit      = imago_commit,
        per_run_log       = log)


function load_manifest_v2(path):
    # Strict-refusal validator implementing manifest
    # rules 1-9 of DESIGN 5.7.  Every failure names
    # the failing rule and the offending entry; no
    # warning-and-continue path exists.
    raw = tomllib.load(path)

    # Rule 1: schema_version must equal 2.
    require(raw.get("schema_version") == 2, path,
        "manifest rule 1: schema_version must equal"
        + " 2 (found "
        + str(raw.get("schema_version")) + ")")

    solids               = raw.get("reference_solid",
                                   [])
    seen_ref_ids         = set()
    seen_element_label   = set()
    default_per_element  = {}    # elem -> count of
                                 #   default=true
                                 #   entries
    seen_elements        = set()

    for ref in solids:
        # Rule 2: required per-solid fields.
        for f in ("reference_id", "kpoint_spec",
                  "convergence_threshold"):
            require(f in ref, path,
                "manifest rule 2: [[reference_solid"
                + "]] missing field: " + f)

        rid = ref["reference_id"]

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

        # Per-entry checks (rules 3, 6, 7, 8, 9).
        for entry in ref.get("entry", []):
            # Rule 3: required entry fields.
            for f in ("element", "atom_site",
                      "label", "default",
                      "description"):
                require(f in entry, path,
                    "manifest rule 3:"
                    + " [[reference_solid.entry]] in"
                    + " " + rid + " missing field: "
                    + f)

            elem  = entry["element"]
            label = entry["label"]
            seen_elements.add(elem)

            # Rule 6: (element, label) uniqueness
            # across the entire manifest.
            key = (elem, label)
            require(key not in seen_element_label,
                path,
                "manifest rule 6: duplicate"
                + " (element, label): "
                + str(key))
            seen_element_label.add(key)

            # Rule 7 tally: count default=true
            # entries per element.  Final check
            # happens after the loop.
            if entry["default"]:
                default_per_element[elem] = (
                    default_per_element.get(elem, 0)
                    + 1)

            # Rule 8 / 9 on fingerprint declarations.
            seen_fp = set()
            for fp in entry.get("fingerprint", []):
                require("method" in fp
                        and "sub_spec" in fp, path,
                    "manifest rule 8: fingerprint"
                    + " declaration must carry both"
                    + " method and sub_spec ("
                    + rid + ", label=" + label
                    + ")")
                # Rule 9: method must be a known
                # matcher.
                require(fp["method"] in MATCHERS,
                    path,
                    "manifest rule 9: unknown"
                    + " matcher method '"
                    + fp["method"] + "' (" + rid
                    + ", label=" + label + ")")
                # Rule 8: (method, sub_spec) unique
                # within this entry.  Use the same
                # canonicalization the per-element-
                # database reader uses (rule 8 in
                # 11.1).
                canon = canonicalize_sub_spec(
                    fp["sub_spec"])
                k2 = (fp["method"], canon)
                require(k2 not in seen_fp, path,
                    "manifest rule 8: duplicate"
                    + " (method, sub_spec) in entry "
                    + label + " of " + rid)
                seen_fp.add(k2)

    # Rule 7 post-loop: exactly one default-tagged
    # entry per element that appears in the manifest.
    for elem in seen_elements:
        count = default_per_element.get(elem, 0)
        require(count == 1, path,
            "manifest rule 7: element " + elem
            + " has " + str(count)
            + " default-tagged entries across the"
            + " manifest (need exactly one)")

    return parse_manifest_object(raw)


function is_cached_v2(ref, imago_commit):
    # Cache hit / miss per DESIGN 5.7's contract:
    # field-by-field compare on the snapshotted
    # inputs, then byte-for-byte compare of the
    # structure file.  Reports the specific field
    # that drove a miss for debuggability.
    cache_dir = ("share/atomicBDB/cache/scf/"
                 + ref.reference_id)
    inputs    = path_join(cache_dir, "inputs.toml")
    if not file_exists(inputs):
        return False

    cached = tomllib.load(inputs)

    # Field-by-field compare on scalar SCF inputs.
    # The cache deliberately excludes `entries`
    # (so adding a new harvest target reuses the
    # cached SCF) and `reference_id` (used only as
    # the cache directory name).
    for f in ("kpoint_spec", "convergence_threshold",
              "imago_commit"):
        expected = current_field_value(ref, f,
            imago_commit)
        if cached.get(f) != expected:
            log("cache miss for " + ref.reference_id
                + ": field " + f + " differs ("
                + str(cached.get(f)) + " -> "
                + str(expected) + ")")
            return False

    # Structure-file byte compare.  COD entries
    # fetch at the pinned cod_revision; structure_
    # path entries read from disk.  Cached bytes
    # are the source of truth for the cached SCF
    # result by construction.
    cached_struct  = find_cached_structure(cache_dir)
    current_bytes  = materialize_structure(ref)
    if read_bytes(cached_struct) != current_bytes:
        log("cache miss for " + ref.reference_id
            + ": structure file changed")
        return False

    return True


function harvestFingerprints(ref, spec,
        imago_commit):
    # Build one FingerprintRecord per declared
    # [[reference_solid.entry.fingerprint]].  The
    # Python-side / Fortran-side split mirrors the
    # consumer's matcher dispatch (11.3.a): the
    # matcher knows whether it computes in process
    # or needs a loen run.
    fingerprints = []
    for fp_decl in spec.fingerprints:
        method   = fp_decl["method"]
        sub_spec = fp_decl["sub_spec"]
        matcher  = MATCHERS[method]()
        if matcher.needs_loen_run:
            payload = harvestLoenFingerprint(ref,
                spec.atom_site, matcher, sub_spec,
                imago_commit)
        else:
            payload = harvestPythonFingerprint(ref,
                spec.atom_site, matcher, sub_spec)
        fingerprints.append(FingerprintRecord(
            method   = method,
            sub_spec = sub_spec,
            payload  = payload))
    return fingerprints


function harvestLoenFingerprint(ref, atom_site,
        matcher, sub_spec, imago_commit):
    # Cache one fort.21 per (method, sub_spec)
    # alongside the SCF cache.  The sub_spec is
    # slugged deterministically so the file name
    # encodes the parameters; two declarations whose
    # sub_spec differs in any key (or value)
    # produce different slugs and different cache
    # files by construction.
    cache_dir = ("share/atomicBDB/cache/scf/"
                 + ref.reference_id)
    slug      = sub_spec_slug(sub_spec)
    out_path  = path_join(cache_dir,
        "loen-" + matcher.name + "-" + slug
        + ".out")

    if not is_cached_loen(out_path, ref,
            imago_commit):
        # Run imago.py -loen -scf no against the
        # cached SCF directory.  The matcher's
        # loen-side parameters flow into the
        # LOEN_INPUT_DATA block; the SCF result
        # itself is reused from the cache.
        run_imago_loen(
            workdir  = cache_dir,
            sub_spec = sub_spec,
            matcher  = matcher,
            output   = out_path)

    rows = matcher.parse_loen_output(out_path,
        sub_spec)
    # atom_site is 1-based per the manifest
    # contract; the row list is 0-indexed.  The
    # matcher's build_payload accessor wraps the
    # vector in the per-matcher payload shape
    # (DESIGN 5.2: bispec uses `values`, reduce
    # uses `shell_code`) so the producer and
    # consumer stay symmetric on field naming.
    return matcher.build_payload(
        rows[atom_site - 1])


function is_cached_loen(out_path, ref, imago_commit):
    # Loen output is valid iff (a) the file exists
    # and (b) the parent SCF cache it depends on is
    # also still a hit.  The (method, sub_spec) are
    # already encoded in the file name (via
    # sub_spec_slug, see above) so no separate
    # parameter comparison is needed -- a sub_spec
    # change produces a different slug and a
    # different cache path by construction.
    if not file_exists(out_path):
        return False
    return is_cached_v2(ref, imago_commit)


function harvestPythonFingerprint(ref, atom_site,
        matcher, sub_spec):
    # In-process matcher call against the cached
    # structure file.  No subprocess, no on-disk
    # cache file -- recomputing in-process is
    # cheaper than the cache bookkeeping would be.
    # Wraps the result via matcher.build_payload so
    # the per-matcher payload field name (DESIGN
    # 5.2) flows through the same accessor the
    # loen-side branch uses.
    cached_struct = find_cached_structure(
        "share/atomicBDB/cache/scf/"
        + ref.reference_id)
    structure = read_structure(cached_struct)
    vectors = matcher.compute_query(structure,
        sub_spec)
    return matcher.build_payload(
        vectors[atom_site - 1])


function sub_spec_slug(sub_spec):
    # Deterministic file-name slug.  Keys in
    # alphabetical order, joined as "key_value"
    # segments, hyphen-separated.  Floats format
    # as "%.6g" -- long enough to disambiguate the
    # parameters humans actually pick, short
    # enough to fit on a file-name line.
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
    # is set to true iff the manifest contributes
    # no other entry for `elem` -- so the per-
    # element database always has exactly one
    # default-tagged entry (rule 7 of 5.2).
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
        coefficients  = coeff1.coefficients,
        alphas        = coeff1.alphas,
        provenance    = {
            "source"       : "atomSCF",
            "commit"       : commit,
            "generated_at" : ts},
        fingerprints  = [])


function is_isolated_default_for(elem, manifest):
    # True iff the manifest declares no other
    # default-tagged entry for `elem`.  Manifest
    # rule 7 forbids zero defaults for any element
    # that *does* appear in the manifest, so a
    # manifest entry with default=true always wins
    # over the isolated baseline.  For elements
    # with no manifest entry at all, the isolated
    # baseline is the database file's only entry
    # and must therefore carry the default tag.
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

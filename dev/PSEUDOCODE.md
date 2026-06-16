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
  11.3.f       5.10         (makegroups.py bispectrum
                             grouping)
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

  - No bispectrum grouping is involved; that happens
    earlier, in makegroups (11.3.f), and never inside this
    makeinput driver.
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
the producer (11.4) and the `makegroups` bispectrum flow
(11.3.f) both call `to_loen_input` and `parse_loen_output`
on matchers whose `needs_loen_run` is true.  The protocol
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
fallback path, and confirms (when an environment-based
scheme is active) that every element's database covers
the requested `(method, sub_spec)`.  Failing fast here
keeps a doomed run out of the "nothing to harvest" case.
(For bispectrum the equivalent coverage check lives in
makegroups, 11.3.f, since its grouping runs before
makeinput.)

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

```
function group_by_bispectrum(skeleton_path, sub_spec,
        similarity_floor):
    matcher = MATCHERS["bispectrum"]()

    # 1. First makeinput: a provisional imago.dat with no
    #    grouping.  The LOEN_INPUT_DATA block carries the
    #    sub_spec via matcher.to_loen_input; the potential
    #    is irrelevant (bispectrum is geometric), so the
    #    default-tagged entry is used and no -pot is given.
    run_makeinput(skeleton_path,
        loen_params = matcher.to_loen_input(sub_spec))

    # 2. Run loen.  -scf no skips the SCF; loen needs only
    #    the structure and the LOEN block.  Produces a
    #    self-describing fort.21 (DESIGN 5.10.3).
    run_imago(flags = ["-loen", "-scf", "no"])

    # 3. Read fort.21.  Each row carries its own identity
    #    (site#, element, species, type_in_species,
    #    type_flat) and the bispectrum vector, so the
    #    row -> atom mapping is read off the file -- no
    #    separate datSkl.map lookup (DESIGN 5.10.3).
    rows = matcher.parse_loen_output("fort.21", sub_spec)

    # 4. Bucket atoms by fingerprint distance within the
    #    floor, per element, refreshing each bucket's
    #    representative as it grows (the same bucketing as
    #    11.3.c, but run here in the orchestrator rather
    #    than inside makeinput).
    species_of = bucketByFingerprint(rows, matcher,
        similarity_floor)

    # 5. Rewrite the skeleton with explicit per-element
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
    #    see the same vectors.  Only Python-side
    #    matchers (reduce) reach this driver:
    #    compute_query stays in process.  Loen-side
    #    (bispectrum) grouping never runs here -- it is
    #    done ahead of makeinput by makegroups (11.3.f),
    #    so those atoms arrive already typed.
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
precision-level numerics, free metadata).  The
producer is a **kaleidoscope client**: it runs no SCF
itself.  It works in three phases -- *build*,
*dispatch*, *harvest*.  The build phase always rebuilds
the "isolated" entries from current `pot1`/`coeff1`
files (so atomSCF changes propagate), then, per
reference solid, materializes the structure, asks the
guidance predictor for a verification grid
(`build_kpoint_convergence`, 15.6), and emits one `CalcUnit`
per k-density grid point plus one structure-only
`imago.py -loen -scf no` unit per Fortran-side
fingerprint.  The dispatch phase hands every collected
unit to kaleidoscope as one flat batch.  The harvest
phase picks each solid's converged grid point, extracts
the potential, harvests one `FingerprintRecord` per
declared `[[reference_solid.entry.fingerprint]]`, and
contributes the same converged point back to the
guidance dataspace.

The v2 manifest reader (`load_manifest_v2` below)
enforces the nine validation rules of DESIGN 5.7 --
the v1-era three-rule reader has been retired with
the rest of schema v1.  Producer-side fingerprint
harvest splits Python-side matchers (in-process,
descriptor-agnostic) from Fortran-side matchers, which
read the `fort.21` of the `-loen` unit kaleidoscope
already dispatched (no separate loen cache; the
kaleidoscope run-reuse cache of DESIGN 6.2.5 subsumes
it).

```
function buildInitialPotentials(manifest_path,
        force, single_element):
    manifest     = load_manifest_v2(manifest_path)
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
        # Drop and rebuild the isolated entry so
        # atomSCF refreshes always propagate.
        db.potentials = [e for e in db.potentials
                         if e.label != "isolated"]
        db.potentials.insert(0,
            build_isolated_entry(elem, imago_commit,
                timestamp, manifest))
        databases[elem] = db

    # Step 1b: per reference solid, materialize the
    # structure and build its verification grid.  Every
    # solid's units accumulate into ONE flight so the
    # whole producer run dispatches as a single flat
    # parallel batch (DESIGN 5.7 / 6.2.1).
    all_units   = []
    predictions = {}    # reference_id -> PredictionRecord
    struct_of   = {}    # reference_id -> local struct path
    for ref in manifest.reference_solids:
        struct = materialize_structure(ref)
        struct_of[ref.reference_id] = struct

        # Fixed (non-swept) RUN SETTINGS, already translated
        # into the tools' own coded vocabulary (DESIGN 6.2.10):
        # make_producer_options maps the manifest's
        # human-readable physics -- functional -> xccode,
        # kpoint_integration -> scfkpint, basis -> the imago
        # scf_basis, scf_threshold -> converg, shift ->
        # kpshift -- and adds the imago_commit cache identity.
        # These are tool-facing only; the wingbeat (§13.2)
        # routes each key to the tool that recognises it.
        options = make_producer_options(ref, imago_commit)

        # The predictor and the PredictionRecord speak the human
        # physics names, not the codes, so the sub-model travels
        # to the builder in its OWN dict -- never mixed into the
        # tool-facing options (DESIGN 6.2.8 / 6.2.10), which would
        # both duplicate the basis and make makeinput reject
        # "functional" / "kpoint_integration".
        submodel = {"basis":              ref.basis,
                    "functional":         ref.functional,
                    "kpoint_integration": ref.kpoint_integration}

        # One builder call for every mode (DESIGN 6.2.9).
        # `center` carries the curator override when the
        # manifest pins kpoint_spec.density (the builder then
        # bypasses the predictor); otherwise it is None and the
        # builder runs full predict-then-verify, falling back to
        # the wide default grid inside itself when the dataspace
        # cannot predict.  Only the units are kept here; the
        # combined flight below supplies the shared root.
        flight_i, record_i = build_kpoint_convergence(
            struct, options, dataspace, ref.system_type, submodel,
            id     = ref.reference_id,
            center = ref.kpoint_spec.density)

        # Domain-specific loen units: the generic builder
        # knows nothing about fingerprints, so the
        # producer appends one structure-only
        # `-loen -scf no` unit per Fortran-side
        # declaration.  The bispectrum fingerprint
        # depends on geometry alone, so these need not
        # wait for the converged grid point.  Each is
        # tagged kind = "fingerprint" (DESIGN 6.2.9) so the
        # convergence harvest skips it and only the
        # fingerprint harvest reads its fort.21.
        loen_units = build_loen_units(ref, struct,
            workspace)

        all_units.extend(flight_i.units)
        all_units.extend(loen_units)
        # Store the plain dict (as attach_prediction_record
        # does), since metadata must be TOML-serializable.
        predictions[ref.reference_id] = as_dict(record_i)

    # The combined flight carries every solid's units and
    # stashes the per-solid prediction records in the
    # opaque metadata dict so the harvest recovers them
    # without re-reading the dataspace.
    flight = Flight(
        units    = all_units,
        root     = workspace,
        sweep    = SweepRecord(
            varied_axes = ("kpt-density",),
            fixed_axes  = {}),
        metadata = {"predictions": predictions})

    # ===== Phase 2: dispatch ==========================
    # Kaleidoscope runs and tracks every unit through the
    # wingbeat seam and its run-reuse cache (DESIGN
    # 6.2.5).  `force` bypasses cached runs so they
    # re-run; fresh results still repopulate the cache.
    dispatch(flight,
        executor = curation_executor(force))

    # ===== Phase 3: harvest ===========================
    log = []
    for ref in manifest.reference_solids:
        struct = struct_of[ref.reference_id]
        # Two-sided delta-below-threshold rule (DESIGN
        # 7.8): the converged grid point's run, or None
        # if nothing converged.
        converged = pick_converged_unit(flight, ref,
            scf_threshold = ref.scf_threshold)
        if converged is None:
            # Non-convergence (e.g. at the top of the
            # range): flag it, harvest no potential.
            log.append(
                make_nonconverged_log_entry(ref))
            continue
        log.append(make_run_log_entry(ref, converged))

        for spec in ref.entries:
            elem = spec.element
            if elem not in databases:
                continue   # filtered out by --element
            # extract_potential reads the site's TYPE block
            # from the converged scfV (NUM_TYPES + per-type;
            # the type number comes from datSkl.map, 9.7).
            coeffs, alphas = extract_potential(
                converged, spec.atom_site)
            # One FingerprintRecord per declared
            # [[reference_solid.entry.fingerprint]].
            # Python-side matchers compute in-process
            # from `struct`; Fortran-side matchers read
            # the loen unit kaleidoscope already ran.
            fingerprints = harvestFingerprints(flight,
                ref, spec, struct)
            new = PotentialEntry(
                label         = spec.label,
                default       = spec.default,
                description   = spec.description,
                num_gaussians = len(coeffs),
                alpha_min     = min(alphas),
                alpha_max     = max(alphas),
                coefficients  = coeffs,
                alphas        = alphas,
                # Provenance records ref.system_type for
                # forensics (DESIGN 5.7 rule 2).
                provenance    = make_imago_provenance(
                    imago_commit, timestamp,
                    ref, spec.atom_site,
                    converged.iterations),
                fingerprints  = fingerprints)
            db = databases[elem]
            # Replace any prior entry with the same
            # label.  Manifest rule 6 enforces
            # cross-solid (element, label) uniqueness,
            # so the only possible prior is the same
            # entry from a previous run.
            db.potentials = [e for e in db.potentials
                             if e.label != spec.label]
            db.potentials.append(new)

    # Step 3b: guidance contribution.  The same
    # converged grid points feed the historical guidance
    # dataspace staging, so every solid the producer
    # converges sharpens the predictor.  harvest_flight
    # recovers each solid's prediction from
    # flight.metadata["predictions"] and skips trust-mode
    # (length-1) solids -- a single point is weak
    # evidence (DESIGN 7).
    harvest_flight(workspace,
        "share/historicalGuidanceDB/", dataspace)

    # ===== Write outputs ==============================
    # All affected element files via the deterministic
    # emitter (5.5).
    for elem, db in databases.items():
        save(db, element_path(elem))

    # Run log capturing the manifest snapshot, per-run
    # iteration counts, the converged k-density per
    # solid, and the Imago commit.  The validation
    # harness (11.5) reads this log.
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
        # Rule 2: required per-solid fields.  basis,
        # functional, and kpoint_integration are required
        # so that nothing the producer emits depends on an
        # implicit default (VISION Principle 5); together
        # with system_type they select the guidance
        # predictor's sub-model (DESIGN 7.6).
        for f in ("reference_id", "system_type", "basis",
                  "functional", "kpoint_integration",
                  "kpoint_spec", "scf_threshold"):
            require(f in ref, path,
                "manifest rule 2: [[reference_solid"
                + "]] missing field: " + f)

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


function harvestFingerprints(flight, ref, spec,
        result_toml):
    # Build one FingerprintRecord per declared
    # [[reference_solid.entry.fingerprint]].  An entry
    # with no declarations harvests nothing and never
    # touches the structure -- the common case so far.
    if spec.fingerprints is empty:
        return []

    # INTERIM (until C55/C58): the Fortran-side bispectrum
    # harvest is not built yet, so refuse any loen-side
    # declaration up front rather than silently dropping a
    # fingerprint the curator asked for.  When C55/C58 land,
    # this guard is replaced by a per-declaration dispatch to
    # harvestLoenFingerprint (specified below) for every
    # matcher whose needs_loen_run is true.
    for fp_decl in spec.fingerprints:
        if MATCHERS[fp_decl["method"]]().needs_loen_run:
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
    # conditions enter only here).  Today only one declaration
    # ever appears, so "max" is just that cutoff; the build-
    # once form matters only once differing cutoffs coexist.
    max_cutoff = max(fp_decl["sub_spec"]["cutoff"]
                     for fp_decl in spec.fingerprints)
    structure = read_structure(
        result_toml.outputs["structure"])
    build_min_dist_matrix(structure, max_cutoff)

    # The expanded skeleton is ordered by the run's sorted
    # (dat) numbering, but atom_site is a skeleton index, so
    # map it to the structure row through datSkl.map (the same
    # map step i reads); the map yields both the row and that
    # row's element symbol.
    (dat_index, map_element) = skeleton_to_dat(
        result_toml.outputs["datSkl_map"])[spec.atom_site]
    # Guard the numbering assumption: the structure row and
    # the map must name the same element, or the expansion and
    # the map have desynced and the fingerprint would describe
    # the wrong atom.  Strict refusal beats a silent mismatch.
    if lower(structure.atom_element_name[dat_index])
            != lower(map_element):
        raise ValueError(
            f"site {spec.atom_site}: datSkl.map names "
            f"{map_element} but the expanded structure row "
            f"{dat_index} is a different element; "
            f"numbering desync")

    fingerprints = []
    for fp_decl in spec.fingerprints:
        method   = fp_decl["method"]
        sub_spec = fp_decl["sub_spec"]
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
            method   = method,
            sub_spec = sub_spec,
            payload  = payload))
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

    iters  = None
    energy = None
    conv   = False
    if "iteration" in outputs:
        # One read of the last data row yields the
        # convergence metric, the total energy, and the
        # iteration count together (all 1-based columns).
        row       = last_data_row(outputs["iteration"])
        threshold = read_scf_threshold(
                        join(run_dir, "imago.dat"))
        # Column 4 is the SCF convergence metric;
        # converged iff it is below the imago.dat criterion.
        conv   = (column(row, 4) < threshold)
        # Column 5 is the last iteration's total energy
        # (converged or not).
        energy = column(row, 5)
        # Column 1 is a per-run cycle counter that resets
        # to 1 each SCF invocation, so the last row's value
        # is THIS run's iteration count -- reruns append
        # rows but never inflate it.
        iters  = int(column(row, 1))

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
        outputs = outputs, job = job_identity(settings),
        runtime_seconds = seconds,
        message = status.name)
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
dataclass KeyFields:
    scalars : dict   # verbatim-compared identity fields,
                     #   e.g. {kpoint_spec, scf_threshold,
                     #   imago_commit}
    files   : list   # logical names of key files to
                     #   byte-compare, e.g. ["structure"]

dataclass CalcUnit:
    id          : str          # stable per-structure key
    calc        : tuple[str,...]  # per-axis directory
                               #   components (DESIGN
                               #   6.2.1); () = no second
                               #   level; one element per
                               #   varied sweep axis
    structure   : str          # path to an imago.skl
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
        # Default wingbeat: drive the §12 API.  The wingbeat
        # owns the makeinput/imago option split (DESIGN
        # 6.2.10): route each option to the tool that
        # recognises it, drop the cache-only build identity,
        # then build the run dir and run it.
        mk_opts    = {}
        imago_opts = {}
        for key, value in unit.options.items():
            if key in imago.OPTION_KEYS:      # job, edge,
                imago_opts[key] = value       #   scf_basis...
            else if key in CACHE_ONLY_KEYS:   # imago_commit:
                continue                      #   dropped, 6.2.5
            else:
                mk_opts[key] = value          # strict makeinput
        makeinput.build_run_dir(
            unit.structure, mk_opts, wingbeat_dir)
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

An ASE wingbeat (D12) and future adapters implement the
same protocol; the dispatch core (13.5) never changes
when one is added.

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
    # Key files: byte-compare each declared key file's
    # current source against the copy already staged in the
    # run dir.  No hashing -- a developer can diff the files
    # to see why a cache missed (DESIGN 6.2.5 / 5.7).
    for logical in unit.key_fields.files:
        current = key_file_source(unit, logical)
        staged  = key_file_staged(wingbeat_dir, logical)
        if not exists(staged) \
           or not files_byte_equal(current, staged):
            return False
    return True
```

```
function write_cache_key(wingbeat_dir, unit):
    # The identity snapshot, written on launch (13.5).
    write_toml(join(wingbeat_dir, "cache_key.toml"),
        { scalars = unit.key_fields.scalars,
          files   = unit.key_fields.files })
```

### 13.5 Dispatch driver (DESIGN 6.2.3)

One Parsl app per unit; per-future exception capture so a
single failure never aborts the flight (Principle 10).
Resuming a flight is just re-running it: the hit-test
skips the `done` units and re-dispatches the rest.

```
function dispatch(flight):
    validate_flight(flight)            # 13.3
    makedirs(flight.root, exist_ok=True)
    serialize_flight(flight)           # 13.1

    with parsl_loaded(flight.parsl_config):
        pending = []                  # list of (unit, future)
        for unit in flight.units:
            pending.append(
                (unit, dispatch_unit(flight, unit)))

        # Gather.  Catch PER future; never let one failure
        # propagate out of the loop.
        entries = []
        for (unit, fut) in pending:
            entry = collect_future(flight, unit, fut)
            entries.append(entry)
            if flight.on_outcome is not None:
                flight.on_outcome(entry)    # stream hook

    return FlightReport(entries = entries)
```

```
function dispatch_unit(flight, unit):
    wingbeat_dir = unit_run_dir(flight, unit)
    if is_cache_hit(unit, wingbeat_dir):         # 13.4
        # Hit: no Parsl app; resolve immediately from the
        # existing status.toml.
        return completed_future(
            report_entry_from_status(unit, wingbeat_dir))
    # Miss: prepare the dir, snapshot the key, mark queued,
    # and submit the wingbeat as a python_app.
    makedirs(wingbeat_dir, exist_ok=True)
    write_cache_key(wingbeat_dir, unit)          # 13.4
    write_status(wingbeat_dir, id=unit.id, calc=unit.calc,
        status="queued",
        wingbeat=(unit.wingbeat or flight.default_wingbeat),
        submitted_at=now())
    return submit_app(execute_wingbeat_task, flight, unit, wingbeat_dir)
```

```
@parsl_python_app
function execute_wingbeat_task(flight, unit, wingbeat_dir):
    # Runs on a worker.  Returns the WingbeatOutcome; raising
    # here surfaces to collect_future as a worker-side
    # failure.
    write_status(wingbeat_dir, id=unit.id, calc=unit.calc,
        status="running", started_at=now())
    wingbeat  = resolve_wingbeat(flight, unit)   # name->Wingbeat
    outcome = wingbeat.run(unit, wingbeat_dir)         # 13.2
    write_status(wingbeat_dir, id=unit.id, calc=unit.calc,
        status=("done" if outcome.ok else "failed"),
        detail=outcome.detail,
        finished_at=now(),
        runtime_seconds=outcome.runtime_seconds,
        message=outcome.message)
    return outcome
```

```
function collect_future(flight, unit, fut):
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
    # In every case status.toml is now terminal; build the
    # report entry from it (single source of truth).
    return report_entry_from_status(unit, wingbeat_dir)
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
k_min          = 3        # below this, refuse the sub-model
k_neighbors    = 5        # neighbors used at each k-NN stage
epsilon        = 1e-6     # numerical floor on distance
w_comp         = 1.0      # composition weight in d1
w_latt         = 0.25     # lattice-family weight in d1
w_gap          = 1.0      # gap weight in d2
w_spin         = 0.5      # magnetization weight in d2
sigma_gap      = 1.0      # gap normalization (eV) in d2
sigma_spin     = 0.5      # magnetization normalization in d2
#                           (Bohr magnetons per atom)
sigma_gap_ref  = 1.0      # gap-spread -> confidence_1 (eV)
sigma_kpd_ref  = 50.0     # kpd-spread -> confidence_2
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
    return Verification(
        grid_values            = tuple(grid),
        grid_energies          = (tuple(energies)
                                  if energies is not None
                                  else None),
        converged_at           = v["converged_at"],
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
    if len(exact) >= k_min:
        return exact, False

    # 2. Most-populous (basis, functional) sub-model within
    #    the same functional family (DESIGN 7.6: (mb,gga-pbe)
    #    -> (fb,gga-pbe)).
    fam = functional_family(functional)
    family = [e for e in pool
              if functional_family(e.context.functional) == fam]
    best = most_populous_submodel(family)
    if len(best) >= k_min:
        return best, False

    # 3. The whole system_type pool, context ignored.
    if len(pool) >= k_min:
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
helper.  Weights are `1/(d+epsilon)` normalized to sum 1.0.

```
function knn_weights(entries, distance_of):
    # distance_of(entry) -> float >= 0.  Returns the
    # k_neighbors nearest as (entry, weight) pairs.
    scored = sort(entries, key = distance_of)        # ascending
    nearest = scored[: min(k_neighbors, len(scored))]
    raw = [1.0 / (distance_of(e) + epsilon) for e in nearest]
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
        comp_sq = sum_sq(sub(query.composition_vector,
                             e.signature.composition_vector))
        latt_sq = sum_sq(sub(query.lattice_onehot,
                             e.signature.lattice_onehot))
        return sqrt(w_comp * comp_sq
                    + w_latt * latt_sq / 2.0)

    nbrs = knn_weights(entries, d1)
    pgap = sum(w * e.measured.gap_ev   for e, w in nbrs)
    pmag = sum(w * intensive_mag(e)    for e, w in nbrs)
    # Confidence_1 from the weighted gap variance (7.6).
    var = sum(w * (e.measured.gap_ev - pgap) ** 2
              for e, w in nbrs)
    conf1 = exp(-sqrt(var) / sigma_gap_ref)
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
        gap_term = w_gap  * (pgap - e.measured.gap_ev) ** 2 \
                   / sigma_gap ** 2
        mag_term = w_spin * (pmag - intensive_mag(e)) ** 2 \
                   / sigma_spin ** 2
        return sqrt(gap_term + mag_term)

    nbrs = knn_weights(entries, d2)
    pkpd = sum(w * e.measured.kpoint_density for e, w in nbrs)
    var = sum(w * (e.measured.kpoint_density - pkpd) ** 2
              for e, w in nbrs)
    conf2 = exp(-sqrt(var) / sigma_kpd_ref)
    return pkpd, conf2, [e.entry_id for e, _ in nbrs]
```

### 15.6 Flight-builder helper (DESIGN 6.2.8 / 7.7)

Lives in `src/scripts/kaleidoscope/` (a domain-aware
optional convenience; the dispatch core stays dumb,
Principle 12).  It turns a structure + tool-facing options
+ a separate `submodel` dict (the physics names the predictor
reads) + a loaded Dataspace into a verification-grid `Flight`
plus a `PredictionRecord` the harvest hook later recovers.

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
function default_wide_kpoint_density_grid():
    # DESIGN 7.9: 8-point bracket spanning a factor of 16,
    # used when no usable predictor exists.  Lives HERE, not
    # in the dataspace (an empty crystalline subtree has no
    # content to consult -- the chicken-and-egg of 7.9).
    return [25.0, 50.0, 100.0, 150.0,
            200.0, 250.0, 300.0, 400.0]
```

```
function logspace(lo, hi, n):
    # n geometrically-spaced points, endpoints inclusive.
    if n <= 1:
        return [sqrt(lo * hi)]            # geometric midpoint
    step = (log(hi) - log(lo)) / (n - 1)
    return [exp(log(lo) + i * step) for i in range(n)]

function build_verification_grid(center, confidence):
    # DESIGN 7.7: width and point count scale inversely with
    # predictor confidence.  conf=1 -> tight 3-point span
    # [center/1.2, center*1.2]; conf=0 -> wide 7-point span.
    width    = 1.2 + 1.5 * (1.0 - confidence)
    n_points = round(3 + 4 * (1.0 - confidence))
    return logspace(center / width, center * width, n_points)
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
function build_kpoint_convergence(structure, options, dataspace,
                                  system_type, submodel,
                                  verify = True,
                                  id = None, center = None,
                                  root = ""):
    # The k-point-density convergence flight builder (DESIGN
    # 6.2.8/6.2.9).  `options` holds the dest-keyed coded run
    # settings forwarded verbatim to the tools and carries NO
    # physics-name keys.  `submodel` MUST carry "basis",
    # "functional", and "kpoint_integration" -- the human names
    # that select the predictor sub-model (DESIGN 7.6 step 2) and
    # are stamped on the record; they are kept out of `options`
    # because makeinput would reject them (DESIGN 6.2.10).  When
    # `center` is given it is a curator-pinned k-density (the 5.7
    # kpoint_spec override): the predictor is bypassed and the
    # grid is built around that value instead.  `root` is the
    # workspace root for the returned flight; "" lets a
    # multi-structure producer supply the shared root when it
    # merges the per-structure flights.
    sc = (structure if is_structure_control(structure)
          else load_skl(structure))
    query_sig = compute_signature(
        sc, system_type, dataspace.group_table)

    # Decide the grid + policy (DESIGN 7.7 step 3 / 6.2.9).
    if center is not None:
        # Curator override: a tight verify grid centred on the
        # pinned density, or a single point when not verifying.
        # No prediction is consulted.
        result = None
        grid_values = (build_verification_grid(center, 1.0)
                       if verify else [center])
        policy = "curator_override"
    else:
        result = predict(dataspace, query_sig,
                         submodel["basis"], submodel["functional"],
                         submodel["kpoint_integration"])
        if not verify:
            grid_values = [result.predicted_kpoint_density]
            policy = "trust_no_verify"
        elif result.is_under_trained:
            grid_values = default_wide_kpoint_density_grid()
            policy = "wide_grid_no_prior"
        else:
            grid_values = build_verification_grid(
                result.predicted_kpoint_density, result.confidence)
            policy = "verify_around_prediction"

    # Round + dedupe to integer k-densities so the on-disk
    # tag parses back to exactly the swept value (6.2.4).
    kpd_grid = sorted(set(round(v) for v in grid_values))

    unit_id = id if id is not None else slug_from_path(structure)
    units = []
    for kpd_int in kpd_grid:
        unit_options = dict(options)           # tool-only copy
        unit_options["kpd"] = kpd_int          # makeinput key
        calc_axes = ordered_map({"kpt-density": kpd_int})
        units.append(CalcUnit(
            id         = unit_id,
            calc       = build_calc_tag(calc_axes),
            structure  = structure,
            options    = unit_options,
            wingbeat   = "imago",
            kind       = "convergence",        # default role
            key_fields = standard_key_fields()))

    # Record the sweep shape (DESIGN 6.2.8 step 6) so
    # serialize_flight emits [flight.sweep] and harvest
    # recovers the varied axis without path-parsing.
    # fixed_axes is empty: the (basis/functional/
    # kpoint_integration) sub-model is carried on the
    # per-structure record below, not duplicated here, so the
    # same fact never lives in two places (DESIGN 6.2.9).
    sweep = SweepRecord(
        varied_axes = ("kpt-density",),
        fixed_axes  = {})

    # The prediction record.  In override mode it documents the
    # curator-pinned value (full confidence, no neighbors, no
    # predicted character); otherwise it mirrors the predictor.
    if center is not None:
        record = PredictionRecord(
            policy                   = policy,
            predicted_kpoint_density = float(center),
            confidence               = 1.0,
            is_under_trained         = False,
            neighbor_entry_ids       = (),
            predicted_gap            = None,
            predicted_magnetization  = None,
            system_type              = system_type,
            feature_vector           = query_sig,
            basis                    = submodel["basis"],
            functional               = submodel["functional"],
            kpoint_integration       = submodel["kpoint_integration"])
    else:
        record = PredictionRecord(
            policy                   = policy,
            predicted_kpoint_density = result.predicted_kpoint_density,
            confidence               = result.confidence,
            is_under_trained         = result.is_under_trained,
            neighbor_entry_ids       = result.neighbor_entry_ids,
            predicted_gap            = result.predicted_gap,
            predicted_magnetization  = result.predicted_magnetization,
            system_type              = system_type,
            feature_vector           = query_sig,
            basis                    = submodel["basis"],
            functional               = submodel["functional"],
            kpoint_integration       = submodel["kpoint_integration"])

    # One flight, one structure: stash the record in the
    # predictions mapping under this structure's id (DESIGN
    # 6.2.9).  A multi-structure producer merges these one-entry
    # mappings into a combined flight.
    flight = Flight(units = units, sweep = sweep, ...)
    attach_prediction_record(flight, unit_id, record)
    return flight, record
```

```
function standard_key_fields():
    # DESIGN 6.2.1: the producer's cache identity -- scalar
    # scf_threshold + imago_commit, with the
    # structure file byte-compared.
    return KeyFields(
        scalars = {"scf_threshold", "imago_commit"},
        files   = ["structure"])
```

**Attaching the PredictionRecord without teaching the core
about it.**  DESIGN 6.2.8 step 6 / 6.2.9 / 7.7 step 6 say the
record is attached to "the Flight's metadata" and serialized
under `[flight.predictions.<id>]`, but the §13.1 `Flight` has
no field for it, and the dispatch core must not interpret
domain data (Principle 9).  The pseudocode pins this the way
the opaque `WingbeatOutcome.detail` is handled: `Flight`
carries a generic `metadata: dict[str, dict]` (default
empty) that `serialize_flight` (§13.1) emits VERBATIM as
`[flight.<key>]` tables, never reading the contents.

```
function attach_prediction_record(flight, structure_id, record):
    # Domain-aware helper stashes a plain dict keyed by the
    # structure id (DESIGN 6.2.9), so one combined flight can
    # carry many structures' predictions.  The core only
    # round-trips it; harvest (15.7) reads it back by id.
    flight.metadata.setdefault("predictions", {})
    flight.metadata["predictions"][structure_id] = as_dict(record)
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

Two v1 conventions: `metric_threshold = scf_threshold` (the
same criterion the SCF converged to is reused as the
grid-flatness threshold), and `imago_commit` falls back to
`"unknown"` when the producer injected none.  `spin_polarization`
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

        # b. Parse each converged run's result.toml.
        kpds, energies, rts = [], [], []
        for u in grid:
            rt = read_result_toml(
                join(workspace_root, "wingbeats", u.id, *u.calc,
                     "result.toml"))
            kpds.append(swept_value_of(u, axis))
            energies.append(rt["total_energy"])
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

        # metric_threshold = the SCF criterion the runs used (the
        #   v1 convention); it is a per-run fact in result.toml.
        thr = rts[0]["scf_threshold"]

        # d. Pick the converged grid point (DESIGN 7.8 3c).
        idx = pick_converged(energies, thr)

        # e. Non-convergence at the top of the range: tag and
        #    SKIP -- a non-converged sweep earns no entry.
        if idx is None:
            warn(unit_id + ": energy still moving at top of"
                 + " grid -- skipped")
            tag_prediction_mismatch(workspace_root, unit_id)
            continue

        # f. Signature: system_type from the prediction record (it
        #    carries it from 7.7; the record is guaranteed present
        #    -- record-less structures were skipped at the loop
        #    top); composition + lattice via compute_signature.
        #    The SAME loaded structure also supplies the cell facts
        #    in (g).
        st  = prediction["system_type"]
        sc  = load_skl(grid[0].structure)        # read_input_file
        sig = compute_signature(sc, st, dataspace.group_table)

        # g. Build the rich GuidanceEntry from the three sources.
        chosen_rt  = rts[idx]
        chosen_kpd = kpds[idx]
        entry = GuidanceEntry(
            entry_id     = "",                # set by save_entry
            generated_at = now_iso8601_utc(),
            source       = "flight",
            signature    = sig,
            measured     = Measured(
                gap_ev              = chosen_rt["gap_ev"],
                gap_kind            = chosen_rt["gap_kind"],
                spin_polarization   = 0.0,    # not measured; see
                #                               DESIGN 7.6
                total_magnetization = chosen_rt.get(
                                        "total_magnetization", 0.0),
                kpoint_density      = chosen_kpd),
            context      = Context(
                # sub-model from THIS structure's record (the sole
                #   home; DESIGN 7.8 step 3f / 6.2.9), never from
                #   sweep.fixed_axes.
                basis      = prediction["basis"],
                functional = prediction["functional"],
                kpoint_integration =
                    prediction["kpoint_integration"],
                scf_threshold = thr,          # result.toml
                cell_atom_count = sc.num_atoms,
                cell_volume_per_formula_unit =
                    sc.real_cell_volume * ANGSTROM3_TO_BOHR3),
            verification = Verification(
                grid_values   = tuple(kpds),
                grid_energies = tuple(energies),  # 7.8 / 7.2
                converged_at  = chosen_kpd,
                metric        = "total_energy",
                metric_threshold = thr,
                # prediction is guaranteed present (record-less
                #   structures were skipped at the loop top).
                predictor_confidence   = prediction["confidence"],
                predictor_neighbor_ids =
                    tuple(prediction["neighbor_entry_ids"])),
            provenance   = Provenance(
                flight_id        = flight_id_of(workspace_root),
                source_structure = grid[0].structure,
                imago_commit     = chosen_rt.get("imago_commit")
                                   or "unknown",
                curator          = "guidance_harvest.py"))

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
function pick_converged(energies, threshold):
    # DESIGN 7.8 step 3c: the smallest interior grid index i
    # at which BOTH consecutive-pair energy deltas fall below
    # threshold.  Two-sided so a single-point fluke does not
    # masquerade as convergence.  Returns None if the energy
    # is still moving (no flat interior point).
    for i in range(1, len(energies) - 1):
        below_up   = abs(energies[i] - energies[i + 1]) < threshold
        below_down = abs(energies[i] - energies[i - 1]) < threshold
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

    # 2. Top-three grid points' total-energy variance below
    #    metric_threshold * 10 (convincingly flat).
    top3 = v.grid_energies[-3:]
    if variance(top3) >= v.metric_threshold * 10.0:
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

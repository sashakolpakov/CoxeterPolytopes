# Exact exclusion of H9/14 record 34

## Result

The second-layer H9/14 survivor with index 34 and graph6 string

```text
M?????????O?O?_@?
```

is not the nerve of a compact hyperbolic Coxeter 9-polytope.  The proof is
solver-independent.  It combines exact local `g`-vector arithmetic, forced
minimal-nonface shadows, blocker duality, and a complete bicoloured-graph
census.  The old bounded Z3 timeout is corroborative only and contributes no
logical step.

The missing-edge graph is

\[
K_3\{0,11,13\}\;\sqcup\;K_2\{1,12\}\;\sqcup\;9I.
\]

The resulting theorem-facing partition of the 142 graph shadows is

```text
121 independently impossible + 6 ordinary-polytopal + 15 undecided.
```

This is a theorem about one graph shadow, not a general dimension-nine
non-existence theorem.

## Dehn--Sommerville branches

Write `a` for the number of missing four-faces.  The exact schedules are:

| `a` | minimal A5 count | facets | summed vertex-link `g4` | disposition |
|---:|---:|---:|---:|---|
| 0 | 6 | 114 | 36 | exhaustive six-row census; all target blockers fail ridges |
| 1 | 12 | 112 | 30 | exact residual two-A5 blocker census |
| 2 | 18 | 110 | 24 | forced A5 shadow has at least 19 members |
| 3 | 24 | 108 | 18 | forced A5 shadow has at least 27 members |
| 4 | 30 | 106 | 12 | contradicts strict positivity in 14 vertex links |
| 5 | 36 | 104 | 6 | contradicts strict positivity in 14 vertex links |

## The positive-`a` proof

The exact local A4 interval gate leaves 546 rows: 126 four-subsets of the
nine isolates and 420 rows containing exactly one active vertex.  Every one
of the 420 active rows is impossible in a family:

- an A4 through a `K3` vertex forces eight incident A5s, giving link `g4=0`;
- an A4 through a `K2` vertex forces nine, again giving `g4=0`; and
- if a `K2` vertex lies in two A4s, the two incident shadows have union 17 or
  18, giving link `g4=-2` or `-3`.

These statements are recomputed over all 420 single rows and all 6,972
same-`K2` pairs; they are not inferred from selected representatives.
Consequently every A4 lies in the isolate block.

Each such A4 has ten independent A5 extensions.  The complete shadow-union
censuses are

```text
a=2: 19:1260, 20:6615
a=3: 27:1680, 28:16380, 29:117180, 30:190260.
```

They exceed the scheduled A5 counts 18 and 24.

For `a=1`, isolate symmetry fixes the A4 to `{2,3,4,5}`.  Its ten forced
A5s leave exactly two additional minimal A5s.  Necessary vertex-link
capacities leave 371 rows and 15,645 labelled pairs.  Exactly 1,645 pairs
have no transversal of size at most four.  Their order-five blocker counts
are

```text
90:480, 96:600, 102:490, 108:75.
```

The required count 112 never occurs.

## The `a=0` census

Represent the six missing A5s as rows.  Let

```text
A = {0,11,13},  B = {1,12},  I = {2,...,10}.
```

The complete local link intervals say:

- each row has five vertices and contains at most one active vertex;
- each `A` vertex occurs in at most one row;
- each `B` vertex occurs in at most two rows; and
- each isolate occurs in at most five rows.

Thus the isolate incidence is a bicoloured graph with six row vertices and
nine column vertices.  A row has degree four when it carries an active label
and degree five otherwise; every column has degree at most five.  `genbg`
enumerates the isomorphism classes for edge counts 24 through 30, and the
C++ certificate enumerates every locally admissible `S3(A) x S2(B)` active
colouring.

The complete exact funnel is

```text
1,288,153 bicoloured representatives
26,598,204 active colourings
25,501,701 states with six distinct rows
 3,101,794 states with no transversal of order <=4
        18 states with the required 114 order-five transversals
        18 pure order-five blockers
         0 states satisfying all ridge equations.
```

For each of the 18 target states, an implementation-independent Python
replay reconstructs all `2^14` subsets, the complete blocker and face vector,
and every six-set ridge degree.  Every target has a ridge of degree 1, 3, or
4.  Hence none is even a pseudomanifold, so none can be a simplicial sphere
or a Coxeter nerve.

The six-set test has a useful conceptual form.  If `X` is the active part of
a six-set and `Y` its isolate part, its ridge degree is

\[
d(X,Y)=\sum_{y\in Y}h(X,Y\setminus\{y\})
       +\sum_{x\in X}h(X\setminus\{x\},Y),
\]

where `h(U,Z)=1` exactly when `U union Z` hits the missing-edge graph and all
six A5 rows.  Since a graph cover needs at least three active vertices, only
the unary, binary, and ternary isolate-column cases occur.  The exhaustive
program evaluates the equivalent direct six-set equations.

### Independent 64-pattern ridge CSP

A second, smaller enumeration checks the ridge obstruction without reading
the `genbg` output or any of the 18 target states.  There is one important
logical proviso: for an actual candidate, absence of a transversal of order
at most four first identifies each five-transversal `h=1` with a facet
complement.  Only then does the displayed formula compute ridge degree.  The
formula is not claimed for arbitrary row systems having smaller
transversals.

After deriving the ridge equations for a valid candidate, the independent
audit deliberately relaxes away the no-small-transversal condition, the
isolate column cap, distinctness of the six rows, and the target count 114.
It retains only the exact row sums and ridge equations.  Each isolate column
is one of the 64 six-bit words, ordered nondecreasingly to quotient `S9`.

There are 23 formal active profiles `(a,b0,b1)`; the main census removes
`(2,2,2)` immediately and checks the other 22.  Fix the active four-set
consisting of all three `A` vertices and the second `B` vertex.  Its 36
binary equations, one for each isolate pair, are already inconsistent in 22
of the 23 formal profiles.  The unique binary survivor is `(3,0,0)`.  For
that profile, adding the six minimal-active-cover families of ternary
equations, `6*C(9,3)=504` equations, is inconsistent as well.  Thus the
relaxed CSP has zero survivors.

For comparison, the unreduced binary/ternary formulation has

```text
5*C(9,2) = 180 binary equations
6*C(9,3) = 504 ternary equations
                 684 total equations.
```

The nine unary equations are unnecessary.  The self-contained C++20 audit
uses direct backtracking over the 64 column words, with no SAT/SMT library.
Its frozen transcript gives `22` binary-UNSAT profiles, the single
binary-SAT profile `(3,0,0)`, and `23` final-UNSAT profiles.  This is an
independent corroboration of the larger solver-independent census.

## Reproduction

The fast frozen verifier is

```bash
PYTHONPATH=src python3 scripts/certify_h9_14_record34.py --verify-only
PYTHONPATH=src python3 -m pytest -q tests/test_h9_14_record34_certificate.py
```

The full `a=0` census uses nauty 2.9.3:

```bash
c++ -std=c++20 -O3 scripts/certify_h9_14_record34_a0.cpp -o /tmp/record34_a0
for e in {24..30}; do
  genbg -q -d4:0 -D5:5 6 9 $e:$e
done | /tmp/record34_a0
```

The independent ridge CSP is replayed by

```bash
c++ -std=c++20 -O3 scripts/audit_h9_14_record34_ridge_csp.cpp \
  -o /tmp/record34_ridge_csp
/tmp/record34_ridge_csp
PYTHONPATH=src python3 -m pytest -q \
  tests/test_h9_14_record34_ridge_csp.py
```

The frozen data are
`data/enumerations/h9_14_record34_a0_census.json` and
`data/enumerations/h9_14_record34_certificate.json`.

## Nonclaims

No Coxeter labels or Lorentzian Gram matrices are classified here.  The
argument excludes record 34 at the necessary nerve layer; it does not produce
a new compact polytope and does not exclude the remaining 15 graph shadows.

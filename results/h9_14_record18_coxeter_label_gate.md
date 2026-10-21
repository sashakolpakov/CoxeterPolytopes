# Record 18: selected-nerve compact-Coxeter no-go

## Theorem and scope

The explicit 110-vertex abstract nerve candidate for record 18 is not the
nerve of a compact hyperbolic Coxeter 9-polytope.

This is an unrestricted compact-label theorem for that one selected
incidence completion.  It does not prove that the candidate is a PL sphere
or an ordinary convex polytope, and it does not need either positive result:
the conclusion is conditional only in the direction needed for a no-go.  It
does not exclude another completion of record 18 or a finite-volume
noncompact realization.

## 1. The compact conditions

In a compact Coxeter chamber, every inclusion-minimal missing face supports
a Lanner diagram and the nine walls through every actual vertex form an
elliptic diagram.  A fourteen-wall chamber in `H9` has normalized Gram
inertia `(9,1,4)`.

Among the eight minimal five-faces of the selected record-18 candidate are
the disjoint pair

\[
L=\{0,5,6,7,8\},\qquad
R=\{1,4,9,10,11\}.
\]

Every one of their 25 cross pairs occurs in an actual vertex.  Thus every
cross label is finite, and every affected actual-vertex restriction must be
elliptic in a compact realization.

## 2. Exhaustive symmetry quotient

There are 420 labelled order-five Lanner placements on each support.  The
certificate constructs an eight-element subgroup of the selected nerve
symmetry from

```text
(5 6),  (7 8),  (2 3)(4 9)(5 7)(6 8)(10 11).
```

Every generated permutation is checked on all thirteen minimal nonfaces and
all 110 actual vertices, and it stabilizes `L` and `R`.  The 420 left
placements split into 57 orbits, with size histogram `4^9,8^48`.  All 420
right placements are retained literally, so no independent symmetry on the
right is assumed.

The exact local census is

```text
representative m=3 single-edge trials      598,500
representative m=3 survivors                 3,780
full labelled single-edge trials          4,410,000
full labelled m=3 survivors                  28,800
m=4 survivors                                      0
m=5 survivors                                      0
representative raw two-edge trials         7,182,000
full raw two-edge trials                  52,920,000
state pairs with two individually viable edges            0
```

The 3,780 representative survivors reduce to 3,660 distinct triple orbits:
120 have size four and 3,540 have size eight, accounting exactly for all
28,800 full-space survivors.  The actual-vertex restriction sizes are
`6^16,7^44,8^25`; the first failure histogram is
`6^310276,7^281322,8^3122`.

For normalized Gram `G=I-A`, increasing a finite label or adding a non-right
edge increases the symmetric nonnegative matrix `A` entrywise.
Perron--Frobenius monotonicity promotes failure at the minimal non-right
label `m=3` to arbitrary stronger labels.  The explicit zero counts at
`m=4,5` and the absence of a state with two viable cross edges therefore
make the census unrestricted.

## 3. Exact principal-block obstruction

The residual 28,800 `m=3` cases consist of two exact weighted-block
isomorphism classes, each of full multiplicity 14,400.  For twice the
normalized ten-node Gram matrix, both classes have exact leading-minor sign
string

```text
+++-----+
```

after adjoining `Delta_0=1`.  Direct arithmetic in
`Q(sqrt(2),sqrt(5))` gives inertia `(8,2,0)` in both cases.  Every leading
minor is frozen and independently recomputed; no floating-point eigensign is
used.

If all cross pairs are right-angled, the principal matrix is instead the
orthogonal direct sum of two Lanner blocks and also has negative index two.
Thus every locally possible compact labeling contains a principal Gram block
with at least two negative directions.  Cauchy interlacing contradicts the
single negative direction in ambient inertia `(9,1,4)`, proving the theorem.

## 4. Vinberg finite-volume warning

This theorem does not imply infinite volume.  A finite-volume noncompact
chamber can have parabolic ideal-vertex diagrams, whereas the obstruction
above uses ellipticity of every actual vertex in the compact case.

A noncompact decision requires an actual Coxeter diagram and Lorentzian Gram
matrix followed by Vinberg's full audit: the relevant elliptic edge
subdiagrams and all finite/ideal endpoint extensions must be checked, and
the maximal parabolic cusp subdiagrams must be enumerated globally.  Neither
that finite-volume test nor a nonexistence proof for it is claimed here.

## 5. Frozen artifacts

```text
certificate:
data/enumerations/h9_14_record18_coxeter_labels.json

certificate digest:
0b168250c951ed814f0678b01122f67de8da8943a3d37bd4f5f1d881aed6c88a

certificate file SHA-256:
a485d2b186a9c0c50b7e881d270353f0bbadf3a8a562d3012dc25477da1248e9

independent hostile-replay digest:
dafa1cb1ee54efc83d256b84571d61a4a58e5701bf339b0ec95028d084368032
```

Reproduction commands:

```bash
PYTHONPATH=src:. python3 scripts/independent_h9_14_record18_coxeter_audit.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record18_coxeter_labels.py --verify
PYTHONPATH=src:. python3 -m pytest -q tests/test_h9_14_record18_coxeter_labels.py
```

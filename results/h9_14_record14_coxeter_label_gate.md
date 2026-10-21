# Record 14: selected-nerve compact-Coxeter no-go

## Theorem and exact scope

The explicit ordinary-polytopal record-14 completion certified in
`results/h9_14_record14_polytopality.md` is **not** the nerve of a compact
hyperbolic Coxeter 9-polytope.

This is an unrestricted compact-label theorem for that one selected nerve.
It does not exclude another completion of the same missing-edge graph
shadow.  It also does not exclude a finite-volume noncompact realization of
the selected nerve.

## 1. The compact Vinberg conditions used

For a compact Coxeter chamber in `H9` with fourteen walls:

1. every inclusion-minimal missing face supports a Lanner diagram;
2. the diagram on the nine walls through every actual vertex is elliptic;
3. the normalized full wall Gram matrix has inertia `(9,1,4)`.

The selected nerve has the two disjoint minimal five-faces

\[
L=\{3,4,5,6,7\},\qquad
R=\{2,8,9,10,11\}.
\]

Thus a putative compact labeling restricts to a labelled order-five Lanner
diagram on each of `L` and `R`.  Every cross pair occurs together in an
actual vertex, so every cross label is finite and that actual-vertex
restriction must remain elliptic.

## 2. Exhaustive reduced-symmetry census

There are 420 labelled order-five Lanner placements on each support.  Record
14 does **not** have the full `S5` symmetry used for some earlier multiwedge
examples, so reducing each support to five abstract Lanner types would be
invalid here.

Instead, the certificate verifies directly on both the minimal nonfaces and
all 106 actual vertices the nerve symmetry

\[
S_{\{3,4,5\}}\times S_{\{6,7\}}\times S_{\{9,10,11\}},
\]

of order 72.  It gives the following exact orbit quotients:

| support | all labelled states | orbit representatives | orbit sizes |
|---|---:|---:|---|
| `L` | 420 | 39 | `6^8, 12^31` |
| `R` | 420 | 73 | `3^6, 6^67` |

For every one of the `39*73` state pairs, all 25 cross pairs are tested at
the smallest non-right label `m=3`.  The result is:

```text
representative single-edge trials       71,175
representative m=3 survivors                420
full labelled single-edge trials      4,410,000
full labelled m=3 survivors              28,800
m=4 survivors                                  0
m=5 survivors                                  0
full raw two-edge trials                52,920,000
state pairs with two individually viable edges        0
```

The 420 surviving representative triples form 420 symmetry orbits: 40 have
size 36 and 380 have size 72, accounting exactly for the 28,800 full-space
survivors.  Hence the quotient has not silently lost lower-stabilizer cases.
The actual-vertex restriction sizes are `7^47,8^25`; the first failures in
the 71,175 single-edge trials have histogram `7^70454,8^301`.

Write a normalized finite Coxeter Gram matrix as `G=I-A`, with
`A_ij=cos(pi/m_ij)>=0`.  Increasing a label beyond three, or adding another
non-right cross edge, increases `A` entrywise.  Perron--Frobenius
monotonicity therefore cannot restore positive definiteness after a local
failure at `m=3`.  The explicit `m=4,5` tests leave zero one-edge survivors,
and the absence of a state pair with two individually viable edges excludes
every multiple-cross-edge labeling.

## 3. Exact inertia of the residual blocks

Up to weighted-graph isomorphism, the 28,800 `m=3` exceptions fall into two
classes, each occurring 14,400 times.  For twice the normalized Gram matrix,
the exact leading-principal-minor sign string of both classes is

```text
+++-----+
```

with the initial `Delta_0=1`.  Exact arithmetic in
`Q(sqrt(2),sqrt(5))` gives inertia `(8,2,0)` for each ten-node block.  The
certificate freezes every leading minor; no floating-point eigensign is
used.

If all 25 cross pairs are right-angled, the same principal matrix is an
orthogonal direct sum of two Lanner Gram matrices and again has negative
index two.  Thus every locally possible compact labeling contains a
principal Gram submatrix with two negative directions.  Cauchy interlacing
contradicts the one negative direction of ambient inertia `(9,1,4)`.  This
proves the theorem.

## 4. Finite-volume warning

The conclusion is **not** “infinite volume.”  In a noncompact finite-volume
chamber an ideal vertex has a parabolic, rather than elliptic, local diagram,
so the compact actual-vertex gate above no longer applies.

A finite-volume decision would require a candidate Coxeter labeling and
Lorentzian Gram matrix, followed by the full Vinberg audit: the relevant
elliptic edge subdiagrams must have the required finite or ideal endpoint
extensions, and the maximal parabolic cusp subdiagrams must be enumerated
and checked globally.  Neither that audit nor its negation is contained in
this certificate.

## 5. Frozen artifacts

```text
certificate:
data/enumerations/h9_14_record14_coxeter_labels.json

certificate digest:
c31d43853470a12c6b88fa95e01bbb7b8078ce314a00b9342cc1fae2c1ad8f4e

certificate file SHA-256:
340577ee357f23a214fd3ac54c94bdde9fd0847e22601d848fed26a66dc371a3

independent hostile-replay digest:
b754de10cfd74822fb6e8708546b872d3b3e3cf9892a31e8af320db955567eec
```

Reproduction commands:

```bash
PYTHONPATH=src:. python3 scripts/independent_h9_14_record14_coxeter_audit.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record14_coxeter_labels.py --verify
PYTHONPATH=src:. python3 -m pytest -q tests/test_h9_14_record14_coxeter_labels.py
```

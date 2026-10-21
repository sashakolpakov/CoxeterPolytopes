# Record 35: an exact rational multiwedge realization

## Result

H9/14 second-layer survivor record 35 is the nerve of an ordinary rational
convex simple 9-polytope with fourteen facets.  This upgrades it from a
bounded timeout to an exact ordinary-polytopal object.  This construction
does not supply Coxeter labels or a Lorentzian Gram matrix.  A subsequent
unrestricted actual-vertex obstruction proves that the same nerve cannot be
the nerve of a compact Coxeter polytope; see
`results/h9_14_record35_coxeter_label_gate.md`.
The later `results/h9_14_record35_whole_shadow_compact_no_go.md` upgrades that
selected-nerve obstruction to an exhaustive compact-only theorem for every
incidence completion over the record-35 missing-edge graph.

The graph6 string is

```text
M?????????O?O?w??
```

and the complete minimal-nonface clutter is

\[
\begin{aligned}
E={}&\{\{0,11\},\{1,12\},\{0,13\},\{1,13\},\{2,13\}\},\\
Q={}&\{\{2,3,4,6,8\},\{0,5,7,9,10\},
\{1,5,7,9,10\},\\
&\hspace{35mm}\{4,5,7,9,10\},\{3,6,8,11,12\}\}.
\end{aligned}
\]

The complete blocker has exactly 102 members, all of order five.  Their
complements generate a pure eight-dimensional complex with face vector

\[
(14,86,309,726,1171,1309,986,459,102).
\]

The fixed complex passes ridge degree two, dual connectivity, all enforced
strict local \(g\)- and Macaulay inequalities, spherical Euler characteristic
for every link, and exact mod-two sphere homology for the complex and every
face link.  Polytopality, however, comes from the independent construction
below rather than from these necessary topological tests.

The subsequent record-34 exclusion gives the current theorem-facing partition

```text
121 independently impossible + 6 ordinary-polytopal + 15 undecided.
```

## Inverse multiwedge

There are two clone blocks

\[
C=\{3,6,8\},\qquad D=\{5,7,9,10\}.
\]

Collapse them to carrier vertices \(c=14\) and \(d=15\).  The resulting
nine-vertex base has vertex order

\[
(0,1,2,4,11,12,13,c,d)
\]

and minimal nonfaces

\[
\begin{gathered}
\{0,11\},\{1,12\},\{0,13\},\{1,13\},\{2,13\},
\{0,d\},\{1,d\},\{4,d\},\\
\{2,4,c\},\qquad \{11,12,c\}.
\end{gathered}
\]

Expanding every occurrence of \(c\) to the full block \(C\), and every
occurrence of \(d\) to the full block \(D\), reproduces the ten displayed
record-35 minimal nonfaces exactly.  Thus the target is the simplicial
multiwedge with multiplicities three and four over this base.

## Integral realization of the base

In the vertex order above, take the following points in \(\mathbb Z^4\):

\[
\begin{array}{c|rrrr}
0&0&3&0&-5\\
1&2&5&-1&-1\\
2&2&4&3&-3\\
4&0&5&-4&0\\
11&5&-2&0&3\\
12&-4&-4&-2&1\\
13&3&-5&-1&0\\
c&2&-1&0&-6\\
d&3&-5&1&1
\end{array}
\]

The verifier evaluates the homogeneous \(5\times5\) determinant for every
one of the \(\binom95=126\) five-point subsets.  Every determinant is
nonzero, and the minimum absolute value is 9.  It then checks every
four-subset as a possible supporting tetrahedron.  Exactly 19 support the
configuration, and they are exactly the tetrahedra allowed by the base
minimal nonfaces.  The minimum absolute determinant against a point off one
of these supporting hyperplanes is 21.  Consequently the convex hull is a
uniform integral simplicial 4-polytope with precisely the required boundary.

A simplicial wedge preserves convex polytopality.  Two wedges at \(c\) and
three at \(d\) raise the dimension from four to nine and the vertex count
from nine to fourteen.  The resulting simplicial 9-polytope has exactly the
record-35 nerve as its boundary; its polar is the required simple 9-polytope.

## Algebro-geometric carrier

Translate a rational interior point of the integral base polytope to the
origin and take its polar.  This gives a rational simple 4-polytope.  The five
facet wedges can be made rational, producing a rational simple 9-polytope.
Its rational simplicial normal fan defines a projective \(\mathbb Q\)-factorial
toric 9-fold whose toric boundary-stratum dual complex is the record-35
nerve.

No unimodularity or Delzant calculation has been made, so smoothness and an
SNC toric boundary are not claimed.  This toric carrier records the incidence
design; it does not manufacture Coxeter angles.

## Reproduction and subsequent compact-label theorem

The frozen proof is implemented by

- `src/coxeter_polytopes/combinatorics/h9_14_record35.py`;
- `scripts/certify_h9_14_record35.py`;
- `data/enumerations/h9_14_record35_polytopality.json`; and
- `tests/test_h9_14_record35.py`.

Replay it with

```bash
PYTHONPATH=.:src python3 scripts/certify_h9_14_record35.py
PYTHONPATH=.:src pytest -q tests/test_h9_14_record35.py
```

Its five order-five supports overlap through the two multiwedge carriers
rather than decomposing as orthogonal factors.  The complete order-five
Lanner join has since been checked against all 102 actual vertices.  Every
assignment contains an affine principal subdiagram in an actual rank-nine
vertex, and Perron--Frobenius monotonicity proves that adding arbitrary
non-right finite Coxeter labels cannot repair ellipticity.  Hence record 35
is excluded from the compact Coxeter track, not merely from a sparse label
ansatz.  The theorem does **not** exclude a noncompact finite-volume chamber
with a changed ideal-face structure, for which the finite/ideal vertex data
and Vinberg finite-volume criterion would have to be reformulated and checked
separately.

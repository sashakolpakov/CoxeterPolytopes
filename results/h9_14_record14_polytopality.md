# Record 14: one exact ordinary-polytopal completion

## Theorem and scope

For the H9/14 second-layer graph-shadow record 14

```text
M???????????o?w??
```

there is an explicit completion whose nerve is the boundary complex of a
convex simplicial 9-polytope on 14 vertices.  The construction is a
five-fold simplicial multiwedge of an exactly Gale-realized simplicial
4-polytope on nine vertices.

This is a theorem about **one displayed incidence completion**.  It is not a
classification of every completion of record 14.  In particular, it is not
yet a Coxeter-polytope theorem: no Coxeter labels, Lorentzian Gram matrix, or
hyperbolic realization are supplied.

## 1. The completion

The minimal nonfaces are the following five pairs, two four-sets, and four
five-sets:

```text
{0,12} {1,12} {0,13} {1,13} {2,13}

{1,3,4,5} {1,9,10,11}

{0,2,6,7,8}
{2,8,9,10,11}
{3,4,5,6,7}
{8,9,10,11,12}
```

The missing-edge graph has precisely the three inclusion-minimal covers

\[
\{12,13\},\qquad \{0,1,13\},\qquad \{0,1,2\}.
\]

Taking the complete blocker of the displayed minimal-nonface clutter gives
exactly 106 inclusion-minimal transversals, all of order five.  Their
complements are therefore 106 eight-dimensional facets.  Direct face
enumeration gives

\[
(f_0,\ldots,f_8)=
(14,86,310,733,1192,1344,1020,477,106).
\]

All 477 ridges have degree two, the facet-adjacency graph is connected, every
face link has the spherical Euler characteristic, and every face link has
the mod-2 homology of the appropriate sphere.  The strict local `g` and
Macaulay checks also pass.  These tests verify the selected sphere; ridge
degree two alone is not being used as a finite-volume test.

## 2. Exact Gale realization

There are three clone blocks

\[
\{3,4,5\},\qquad \{6,7\},\qquad \{9,10,11\}.
\]

Collapsing each block to its first member leaves the ordered base set

\[
(0,1,2,3,6,8,9,12,13).
\]

The collapsed complex has dimension four and

\[
(f_0,f_1,f_2,f_3)=(9,28,38,19).
\]

In the displayed base-vertex order, an integral Gale matrix is

\[
G=\begin{pmatrix}
1&0&0&0&31&-53&-52&-31&47\\
0&1&0&0&21&-39&-32&-64&60\\
0&0&1&0&23&-21&-22&-50&45\\
0&0&0&1&50&-20&-44&-10&15
\end{pmatrix}.
\]

All \(\binom94=126\) maximal minors are nonzero.  Its 19 positive circuits
are exactly the complements of the 19 base facets.  Moreover

\[
G\begin{pmatrix}
107156\\652142\\43267\\470205\\32655\\42227\\58450\\128091\\172953
\end{pmatrix}=0,
\]

and every coordinate of this kernel vector is positive.  Every Gale element
is both contained in and avoided by a positive circuit.  Thus Gale duality
produces a convex simplicial 4-polytope with exactly the collapsed boundary
complex.  Expanding the three clone blocks performs

\[
(3-1)+(2-1)+(3-1)=5
\]

simplicial wedges.  Wedges preserve convex polytopality, raise the dimension
from four to nine, and reproduce the eleven displayed minimal nonfaces
exactly.  This proves the theorem.

## 3. Vinberg finite-volume boundary

No hyperbolic volume conclusion follows from the preceding incidence data.
In particular,

- two facets around each ridge is the pseudomanifold condition, not
  Vinberg's finite-volume criterion;
- spherical or homology-spherical links are topological, not assertions that
  the associated Coxeter subdiagrams are elliptic or parabolic;
- an ordinary convex realization does not supply a Lorentzian Gram matrix.

Any positive hyperbolic-volume claim must first label the 1-skeleton and
solve the exact Gram constraints with signature \((9,1)\) and the required
corank.  Only then can one run the full Vinberg audit: enumerate the relevant
elliptic edge subdiagrams and their finite or ideal endpoint extensions, and
separately enumerate and test the maximal parabolic cusp subdiagrams.

The compact branch has since been settled negatively without constructing a
full Gram matrix: the two disjoint five-node Lanner supports force an exact
negative-index-two principal block for every locally possible compact
labeling.  See `results/h9_14_record14_coxeter_label_gate.md`.  That later
theorem still does not decide finite-volume noncompact realizability.

## 4. Reproducibility

The frozen certificate is
`data/enumerations/h9_14_record14_polytopality.json`.  It records the complete
106-facet transcript, the full link/homology audit, the Gale certificate, and
an independently implemented hostile replay.

```text
certificate digest:
0fba773d1449780048625c0243eb4740efc0ccd23d87b266e5f3ca06deb7edb1

certificate file SHA-256:
01ac2e8ca83a72f40dd813ee0051081911bab203be01966902882eb61d751dc1

independent replay digest:
149d05f32c19e056e33d750c228784c28d2d61d06952c3a61018b509195e2ff4

facet-complement transcript SHA-256:
db6456207158b9c5640dca0e3930719547d08ef59c993e7e6dc3d0f635702d3f
```

Reproduction commands:

```bash
python3 scripts/independent_h9_14_record14_audit.py
PYTHONPATH=src:. python3 scripts/certify_h9_14_record14.py --verify-only
python3 -m pytest -q tests/test_h9_14_record14.py
```

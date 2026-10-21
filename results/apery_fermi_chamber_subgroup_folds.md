# Exact subgroup-fold census of Shimada's Apéry--Fermi chamber

## Outcome

Shimada's 80-wall Apéry--Fermi chamber supplies a complete, finite fishing
stock, but it does **not** produce a Coxeter polytope above dimension 9 by any
subgroup fold of its weighted wall graph.

More precisely, the exact weighted wall Gram has an abstract automorphism
group of order 32.  This is exactly the wall action of Shimada's published
group `OSXD0`; it strictly contains the order-16 geometric subgroup
`AutXD0`.  All 58 subgroups of the order-32 group were enumerated.  The
geometric subgroup was not used as a gate.

Among nontrivial subgroup fixed sections above dimension 9 there are exactly:

- six sections in \(\mathbb H^{11}\), each with 22 genuine facets; four fail
  the Coxeter-angle test at exactly two facet pairs and two fail at eight
  pairs;
- five sections in \(\mathbb H^{10}\), each with 28 genuine facets; the best
  fails at four pairs and the other four fail at ten pairs.

Thus the high-dimensional result is a finite relative no-go theorem, not an
unfinished search.  The closest \(\mathbb H^{11}\) sections have only two bad
pairs, both with normalized squared product \(3/8\).  The closest
\(\mathbb H^{10}\) section has four bad pairs, all with normalized squared
product \(1/16\).  Neither number is \(\cos^2(\pi/m)\) for an integral Coxeter
label \(m\geq 2\).

## Positive lower-dimensional controls

The same exhaustive census contains 33 genuine Coxeter sections in dimensions
2 through 8.  Vinberg's complete edge-endpoint criterion was run on every one
of them, and all 33 are finite-volume and noncompact.  There are no merely
Coxeter-compatible but infinite-volume survivors in these dimensions.

The largest positive controls are:

- four \(\mathbb H^8\) sections with 14 facets, 102 finite vertices, 15 ideal
  vertices, and 594 elliptic edge diagrams; every edge has exactly two
  endpoints and there are no bad edges;
- six \(\mathbb H^7\) sections with 14 facets, again finite-volume and
  noncompact.

There are also eight one-dimensional, one-facet restrictions.  They are
recorded separately as degenerate \(\mathbb H^1\) Coxeter-compatible sections;
the Vinberg test was not invoked for them.  No novelty claim is made for any
of the lower-dimensional polytopes.

The four \(\mathbb H^8\) outputs are not four distinct diagrams: their exact
normalized Coxeter diagrams are edge-label-preserving isomorphic.  They should
therefore be treated as one explicit diagram candidate with four realizations
inside this chamber stock.  A scan of the currently vendored CoxIter corpus
finds no dimension-8 file with 14 facets (its dimension-8 files have 9 through
12 facets).  This makes the 14-facet diagram a strong **new-diagram
candidate**, but it is not a global literature novelty proof.

### Exact invariants of the \(\mathbb H^8\) diagram

For each of the four realizations, the saturated integral fixed lattice has
rank 9, determinant \(512=2^9\), inertia \((1,8)\), and Smith invariants

\[
  (2,2,2,2,2,2,2,2,2).
\]

For the representative `subgroup_012` (the involution `OSXD0[30]`), GAP's
exact integral nullspace basis gives the fixed-lattice Gram

\[
\begin{pmatrix}
-208&-274&0&-354&30&-324&12&4&2\\
-274&-420&0&-510&50&-460&18&0&0\\
0&0&-4&-6&6&0&0&0&0\\
-354&-510&-6&-662&78&-580&20&4&2\\
30&50&6&78&-22&54&0&0&0\\
-324&-460&0&-580&54&-526&20&4&2\\
12&18&0&20&0&20&-2&0&0\\
4&0&0&4&0&4&0&-4&2\\
2&0&0&2&0&2&0&2&-2
\end{pmatrix}.
\]

Half this matrix is integral, odd, unimodular, and indefinite of signature
\((1,8)\).  The classification of indefinite odd unimodular lattices therefore
identifies the fixed lattice as

\[
  S_X^H \cong 2I_{1,8}
  =\langle2\rangle\oplus\langle-2\rangle^{\oplus8}.
\]

The 14 facet normals have square multiset \((-4)^8,(-2)^6\).  An independent
CoxIter computation gives

\[
  f=(117,594,1186,1268,825,342,90,14,1),
\]

15 ideal vertices, orbifold Euler characteristic \(17/28800\), and covolume

\[
  \operatorname{vol}=\frac{17\pi^4}{189000}.
\]

CoxIter also marks the group arithmetic.  Its numerical signature display is
unreliable for this input; the exact positive-diagonal normalized Gram has
rank 9 and inertia \((8,1,5)\).  Equivalently, our negative-diagonal normal
Gram has inertia \((1,8,5)\).

### Relation to the standard unimodular reflection group

The lattice identification also explains why this positive example should
not yet be advertised as a globally new arithmetic group.  After dividing
the form by two, its six norm-\(-2\) and eight norm-\(-4\) facet normals are
respectively norm-\(-1\) and norm-\(-2\) integral roots of the standard odd
unimodular lattice \(I_{1,8}\).  Hence the 14 reflections generate a subgroup
of the classical full root-reflection group of \(I_{1,8}\).

The standard nine-facet simplex for that group is the CoxIter input
`group units qf/8-1,1,1,1,1,1,1,1,1.coxiter`; its exact Euler
characteristic is \(17/1393459200\).  Comparing with \(17/28800\) shows that
the Apéry--Fermi reflection subgroup has index

\[
  \frac{17/28800}{17/1393459200}=48384=2^8\,3^3\,7
\]

in the standard simplex reflection group.  Thus the algebraic construction
produces a nontrivial 14-facet finite-index Coxeter subgroup and a new
divisor-theoretic realization; whether this particular subgroup diagram was
previously printed remains a separate literature question.

## Exact reduction and completeness

The Néron--Severi lattice parsed from the source is

\[
  S_X \cong U\oplus E_8(-1)\oplus E_8(-1)\oplus\langle-12\rangle,
\]

with rank 19, determinant 12, and inertia \((1,18)\).  The 80 published wall
normals consist of 24 outer and 56 inner walls.  Their squares are
\(-2,-4/3,-1,-3/4,-1/3,-1/12\).  The class `a32` has square 32, is fixed by
all 32 wall symmetries, and pairs strictly positively with every wall.

For a subgroup \(H\), each wall inequality restricts to the fixed space through
the sum of its normal over the corresponding \(H\)-orbit.  Every resulting
inequality is disposed of by one of two exact certificates:

1. a positive rational dependence expressing it in the cone generated by the
   other orbit inequalities; or
2. an exact timelike point on that wall which pairs strictly positively with
   every other orbit inequality.

Consequently the retained wall lists are complete and irredundant.  Floating
linear programming is used only to discover small supports for the rational
dependences; every dependence and every facet witness is replayed exactly.
The Coxeter gate then uses exact normalized squared products.  Vinberg's
finite-volume conclusion is invoked only after that gate passes and with the
complete simple-wall list.

## Reproduction

- Exact imported stock:
  [`data/known_examples/apery_fermi_chamber_stock.json`](../data/known_examples/apery_fermi_chamber_stock.json)
- Full 58-subgroup certificate:
  [`data/enumerations/apery_fermi_chamber_subgroup_folds.json`](../data/enumerations/apery_fermi_chamber_subgroup_folds.json)
- Parser and verifier:
  [`scripts/search_apery_fermi_chamber_folds.py`](../scripts/search_apery_fermi_chamber_folds.py)
- Regression tests:
  [`tests/test_apery_fermi_chamber_folds.py`](../tests/test_apery_fermi_chamber_folds.py)
- Download provenance:
  [`data/external/apery_fermi_shimada_2025/README.md`](../data/external/apery_fermi_shimada_2025/README.md)

The primary source is Ichiro Shimada's open dataset,
<https://doi.org/10.5281/zenodo.14842583>.  Only the 549 KB exact chamber file
was extracted; the archive's 369 MB auxiliary face table was not needed.

# Strict lower-bound and ridge-support filters for H9/14

## Result and exact scope

The first graph-level archive contains 444 possible missing-edge graphs for
the nerve \(K=\partial Q\) of the polar \(Q=P^*\), where \(Q\) would be a
simplicial 9-polytope on 14 vertices.  Two independent exact filters give:

- the strict generalized lower-bound filter leaves **263** graphs;
- the closed-pseudomanifold ridge-support filter leaves **197** graphs; and
- their intersection contains **142** graphs.

The 444-to-197 count is independent of the generalized lower-bound argument.
The 142 graphs are the frontier at this layer, not completed nerves. Later
global facet work, including the record-23/34/37/54/74/92/93/46 certificates and
the four later multiwedge realizations, gives the current theorem-facing
partition \(121+6+15\); see
[`h9_14_facet_completion.md`](h9_14_facet_completion.md). No homology-sphere,
PL-sphere, oriented-matroid, polytopality, Coxeter-labeling, or Gram-matrix
claim is made for the 142 at this layer.

The complete per-graph support cores and rejection certificates are archived
in
[`data/enumerations/h9_14_nerve_support.json`](../data/enumerations/h9_14_nerve_support.json).
They are produced by
[`scripts/h9_14_nerve_support_filter.py`](../scripts/h9_14_nerve_support_filter.py)
and recomputed without a graph or SAT library by
[`coxeter_polytopes.combinatorics.nerve_support`](../src/coxeter_polytopes/combinatorics/nerve_support.py).

## A strict generalized lower-bound obstruction

For a simplicial \(d\)-polytope, Murai--Nevo's generalized lower-bound theorem
states that, if

\[
 h_{r-1}=h_r,\qquad 1\le r\le d/2,
\]

then the polytope has a geometric \((r-1)\)-stacked triangulation: every face
of the triangulation of dimension at most \(d-r\) is already a boundary face.
This is Theorem 1.2 of their paper.  Its hypotheses apply here because the
object being excluded is the simplicial convex polar \(Q\), not merely an
abstract homology sphere.

Take \(d=9\) and \(r=4\).  The equality \(g_4=h_4-h_3=0\) would give a
3-stacked triangulation with no interior face of dimension at most five.  The
triangulation is nontrivial: \(Q\) has 14 vertices, whereas a single
9-simplex has ten.  A nontrivial triangulation of a ball has an interior face.
Choose an inclusion-minimal one.  It has dimension at least six, hence at
least seven vertices, and all its proper faces lie on the boundary.  Its
vertex set is therefore a missing face of \(K\) of order at least seven.

For a compact hyperbolic Coxeter polytope, every missing face supports a
Lannér diagram and consequently has order at most five.  This contradiction
proves the strict condition

\[
 g_4>0.
\]

If \(a\) denotes the number of missing four-faces and
\(g_4=\beta(H)-a\), the interval from the first-stage \(g\)-theorem filter is
therefore tightened by replacing its upper endpoint with
\(\min\{a_{\max},\beta(H)-1\}\).  This rejects 181 graphs and leaves 263.

This argument generalizes the earlier \(g_2=0\) stacked obstruction.  The new
step uses the equality characterization, not just nonnegativity of the
\(g\)-vector.

## Exact facet-complement formulation

Let \(V\) be the 14 vertices and let \(H\) be the missing-edge graph.  A facet
\(F\) of \(K\) has nine vertices.  Its complement

\[
 B=V\setminus F
\]

has five vertices and must meet every edge of \(H\).  Thus every actual facet
complement is among the five-vertex covers of \(H\).  The converse is not
asserted: a cover is only a candidate because missing four- and five-faces may
exclude its complementary nine-set.

For a six-set \(C\subset V\), consider

\[
 N_C=\#\{B:\ B\text{ is a selected facet complement and }B\subset C\}.
\]

The complement \(R=V\setminus C\) has eight vertices.  If \(N_C>0\), then
\(R\) is a ridge of the pure 8-complex and the closed-pseudomanifold condition
requires exactly two containing facets.  Hence

\[
 N_C\in\{0,2\}.
\]

Crucially, \(N_C=0\) is allowed.  The algorithm never assumes that every
\(H\)-allowed eight-set is a ridge or even a face.

## The greatest ridge-support core

Start with all five-vertex covers of \(H\).  Delete a candidate \(B\) if one
of its nine six-set supersets \(B\cup\{v\}\) contains fewer than two current
candidates, and repeat to a fixed point.

This deletion is complete as a negative test.  Let \(X\) be the selected
facet complements of any actual closed pseudomanifold completion.  Initially
\(X\) lies in the cover set.  Inductively, if \(X\) is contained in a current
set and \(B\in X\), then every six-set containing \(B\) contains exactly two
members of \(X\), so it contains at least two current candidates.  Thus \(B\)
is never deleted.  Every possible completion is contained in the final core.

The no-missing-triples theorem now supplies the rejection certificate.  Every
independent triple \(T\) of \(H\) must lie in a facet, equivalently there must
be a core complement \(B\) with \(B\cap T=\varnothing\).  If none exists, no
global selection can work.

It is enough to test independent triples.  Since \(|E(H)|\le9\), every
independent pair extends to an independent triple: otherwise the other twelve
vertices would require at least twelve distinct edges incident to the pair.
Every vertex first extends to an independent pair and hence to a triple.
Consequently triple coverage also ensures that all vertices and all nonmissing
pairs occur in the putative pure complex.  Purity itself comes from defining
the complex as the downward closure of selected nine-vertex facets.

The archive stores the entire fixed-point core as sorted 14-bit masks, the
number removed in every iteration, and the lexicographically first uncovered
independent triple when one exists.  Recomputing these certificates gives the
independent count **444 to 197**.

## Exact counts

| missing edges \(e\) | first-layer input | pass strict GLBT | pass ridge support | pass both |
|---:|---:|---:|---:|---:|
| 2 | 2 | 2 | 2 | 2 |
| 3 | 5 | 5 | 5 | 5 |
| 4 | 11 | 11 | 11 | 11 |
| 5 | 23 | 21 | 19 | 17 |
| 6 | 50 | 42 | 35 | 29 |
| 7 | 89 | 60 | 46 | 33 |
| 8 | 118 | 63 | 44 | 28 |
| 9 | 146 | 59 | 35 | 17 |
| **total** | **444** | **263** | **197** | **142** |

Across all 444 inputs there are 58,521 initial facet complements and 58,191
members in the fixed-point cores.  The joint 142 survivors have between 82
and 880 core candidates, 26,859 in total.  This explains why enumerating all
global subsets without further propagation remains unattractive even after
the exact local filters.

On the current machine, verification of the 1.45 MB archive takes about 6.5
seconds; generation plus repeated verification takes about 19.6 seconds.  In
addition to the production bitmask implementation, the verifier independently
rebuilds the critical GLBT intervals and ridge cores with frozensets.  It also
checks all 444 Dehn--Sommerville schedules, frozen graph hashes, certificates,
and count vectors.

## Dehn--Sommerville and the forced number of missing five-faces

Let \(I_k(H)\) be the number of independent \(k\)-sets, let
\(\mathcal A_4\) be the selected missing four-faces, and put
\(a=|\mathcal A_4|\).  Then

\[
 f_3=I_4(H)-a.
\]

The relation \(h_5=h_4\) fixes

\[
 f_4=6f_3-21f_2+56f_1-126f_0+252.
\]

If \(c(\mathcal A_4)\) is the number of independent five-sets containing at
least one member of \(\mathcal A_4\), minimality gives

\[
 |\mathcal A_5|=I_5(H)-c(\mathcal A_4)-f_4.
\]

Thus the number of missing five-faces is forced after \(\mathcal A_4\) is
chosen; it is not an independent search parameter.  The manifest archives,
for every still-admissible \(a\), the complete DS-forced \(f\)-vector and the
affine offset \(I_5(H)-f_4\), so that the eventual generator need only subtract
\(c(\mathcal A_4)\).  The higher DS relations similarly fix
\(f_5,f_6,f_7,f_8\).  At this count-only level all 263 strict-GLBT survivors
obey the obvious independent-set upper bounds, so this check causes no
additional graph rejection.

## Positive fixtures

Two ordinary polytopes survive and test every incidence convention:

1. \(\Delta_1^3\times\Delta_3^2\) has 128 nerve facets and 576 ridges.  Its
   missing faces have orders \(2,2,2,4,4\), and its
   \(f\)-vector is
   \((14,88,328,804,1352,1568,1216,576,128)\).
2. \(C_5\times\Delta_3\times\Delta_4\) has 100 nerve facets and 450 ridges.
   Its missing faces have orders five times 2, then 4 and 5, and its
   \(f\)-vector is
   \((14,86,309,725,1165,1295,970,450,100)\).

Their complete facet lists are constructed in the tests.  The verifier checks
ridge degree two, connected facet adjacency, all independent triples, all
minimal nonfaces, the DS vectors, and the forced \(|\mathcal A_5|\) formula.
They are not claimed to be hyperbolic Coxeter examples.

## What remains for an exhaustive nerve enumeration

For each of the 142 graph shadows, introduce a Boolean variable \(x_B\) for
each core complement.  A complete finite search must impose

\[
 \sum_{B\subset C}x_B\in\{0,2\}
 \quad\text{for every six-set }C,
\]

coverage of every independent triple, the exact number and placement of
missing four-faces, the DS-forced missing five-faces, connected facet
adjacency, and all local homology-sphere conditions.  Every model must then be
checked for a PL-sphere certificate.

Even that would enumerate abstract nerves, not polytopes.  The Gale transform
has rank

\[
 14-(9+1)=4,
\]

so a realizable rank-four oriented-matroid or another exact realization
certificate is still required before claiming polytopality for an arbitrary
completion. Records 43, 44, and 47 were subsequently proved polytopal by an
independent multiwedge theorem, without realizing their archived Gale sign
words. A Coxeter labeling and Lorentzian Gram certificate come only after the
ordinary-polytopality gate.

## References

- Satoshi Murai and Eran Nevo, *On the generalized lower bound conjecture for
  polytopes and spheres*, Acta Mathematica 210 (2013), 185--202,
  [DOI 10.1007/s11511-013-0093-y](https://doi.org/10.1007/s11511-013-0093-y),
  especially Theorem 1.2.
- Amanda Burcroff, *Near Classification of Compact Hyperbolic Coxeter
  \(d\)-Polytopes with \(d+4\) Facets and Related Dimension Bounds*,
  [arXiv:2201.03437](https://arxiv.org/abs/2201.03437), for the no-missing-
  triples reduction in this target range.

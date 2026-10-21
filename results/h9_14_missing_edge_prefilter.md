# Exact missing-edge prefilter for the H9/14 theorem track

## Result and scope

Let \(P\) be a hypothetical compact Coxeter 9-polytope with 14 facets, let
\(Q=P^*\) be its simplicial polar, and let \(K=\partial Q\).  Define \(H\) on
the 14 vertices of \(K\) by declaring \(ij\in E(H)\) exactly when
\(\{i,j\}\notin K\), equivalently when the corresponding facets of \(P\)
are disjoint.

The command

```text
geng -q 14 2:10
```

produces 6,789 unlabeled graphs.  The exact \(g_3\) test leaves 592 graphs;
allowing an arbitrary number of missing four-faces and applying the exact
\(g_4\) feasibility test leaves 575.  Of these, the 131 graphs with ten edges
are excluded by the equality case of the lower-bound theorem.  Thus **444
missing-edge graphs** survive this first prefilter.

This is a coverage result for missing-edge graphs, not an enumeration of
nerves.  A surviving graph still needs missing faces of orders four and five,
a simplicial-sphere realization, polytopality, a Coxeter labeling, and a Gram
matrix.

The complete transcript and all per-graph invariants are archived in
[`data/enumerations/h9_14_missing_edges.json`](../data/enumerations/h9_14_missing_edges.json).
They are generated and independently checked by
[`scripts/h9_14_missing_edge_prefilter.py`](../scripts/h9_14_missing_edge_prefilter.py)
and the pure-Python verifier in
[`coxeter_polytopes.combinatorics.missing_edges`](../src/coxeter_polytopes/combinatorics/missing_edges.py).

## Mathematical reduction

Burcroff's Lemma 10.10 proves that the target has no missing face of order
three.  Missing faces of a compact Coxeter polytope correspond to Lannér
subdiagrams, hence have order at most five.  Results of Felikson--Tumarkin
also force at least two pairs of disjoint facets in this facet range.  Write

\[
 e=|E(H)|,\qquad
 w=\sum_{v\in V(H)}\binom{\deg_H(v)}2,\qquad
 t=\#\{\text{triangles of }H\}.
\]

Since there are no missing triples, a triple is a face of \(K\) exactly when
it contains no edge of \(H\).  Inclusion--exclusion gives

\[
 f_1=\binom{14}{2}-e,
 \qquad
 f_2=\binom{14}{3}-12e+w-t.
\]

For the boundary of a simplicial 9-polytope,

\[
 h_1=5,\qquad h_2=15-e,
\]

and direct substitution in the \(h_3\) formula gives

\[
 g_2=h_2-h_1=10-e,
 \qquad
 g_3=h_3-h_2=20-4e+w-t.
\]

The necessity half of the \(g\)-theorem says that the \(g\)-vector is an
M-sequence.  Consequently

\[
 0\le e\le10,
 \qquad
 0\le g_3\le g_2^{\langle2\rangle},
\]

where the upper bound is Macaulay raising.  Combined with the two-disjoint-
pair theorem, only \(2\le e\le10\) need be generated.

The graph also determines the possible range of \(g_4\), even though it does
not determine the missing four-faces themselves.  Let \(I_4(H)\) denote the
number of independent four-sets of \(H\), and let \(a\) be the number of
missing faces of order four.  Every independent four-set has all its pairs,
and the no-missing-triples result gives all its triples.  It is therefore
either a 3-face or exactly a missing four-face.  Hence

\[
 f_3=I_4(H)-a.
\]

Put

\[
 \beta(H)=I_4(H)-6f_2+21f_1-56\cdot14+126-h_3.
\]

The \(h_4\) formula then gives \(g_4=\beta(H)-a\).  The next Macaulay
condition, \(0\le g_4\le g_3^{\langle3\rangle}\), can hold for some possible
number of missing four-faces only if the integer interval

\[
 \max\{0,\beta(H)-g_3^{\langle3\rangle}\}
 \;\le a\le\;
 \min\{I_4(H),\beta(H)\}
\]

is nonempty.  This is an exact necessary graph-level test: it rejects a graph
only when no choice of the number \(a\) can satisfy the \(g_4\) condition.  It
does not assert that every integer in the interval is realized by a compatible
family of missing four-faces.

If \(e=10\), then \(g_2=0\).  The equality case of the lower-bound theorem
implies that \(Q\) is stacked.  Since \(Q\) has 14 vertices rather than the
10 vertices of a 9-simplex, at least one stacking occurred.  The removed
9-vertex facet remains a missing face: stacking removes the facet but none
of its proper faces, and later stackings do not remove proper faces either.
This missing face has order nine, contradicting the Lannér bound five.
Thus the whole \(e=10\) stratum is impossible.

## Exact counts

| \(e\) | all unlabeled graphs | pass through \(g_3\) | admit \(a\) through \(g_4\) | survive stacked test |
|---:|---:|---:|---:|---:|
| 2 | 2 | 2 | 2 | 2 |
| 3 | 5 | 5 | 5 | 5 |
| 4 | 11 | 11 | 11 | 11 |
| 5 | 26 | 25 | 23 | 23 |
| 6 | 68 | 54 | 50 | 50 |
| 7 | 177 | 89 | 89 | 89 |
| 8 | 495 | 118 | 118 | 118 |
| 9 | 1,464 | 146 | 146 | 146 |
| 10 | 4,541 | 142 | 131 | 0 |
| **total** | **6,789** | **592** | **575** | **444** |

The archived nauty provenance is:

- nauty/Traces 2.9301 (Homebrew package 2.9.3, 32-bit graph words);
- executable SHA-256
  `ad2f68adf733dbed7cad543841cfa329740596ee5cd136f17a0a17f6e744f5ad`;
- graph6 stdout SHA-256
  `780c577cc7344d270d4551f1ba6db7834e0eef97e96a7e35da1c8abe8728d64f`.

The independent verifier decodes graph6 itself, recomputes \(e,w,t,I_4(H)\),
all listed \(f,h,g\) entries, both Macaulay bounds, the feasible interval for
\(a\), and every disposition; it also checks the transcript hash and frozen
count vector.  It does not purport to re-prove nauty's canonical-generation
algorithm.

## Two explicit surviving graph shadows

There are elementary simple 9-polytopes satisfying the same missing-face
orders, so the graph filter is non-vacuous:

1. \(\Delta_1^3\times\Delta_3^2\) has 14 facets.  Its minimal nonfaces have
   orders \((2,2,2,4,4)\); its missing-edge graph is a matching of size three,
   and \((g_2,g_3)=(7,8)\).
2. \(C_5\times\Delta_3\times\Delta_4\) has 14 facets.  Its minimal nonfaces
   consist of the five diagonals of the pentagon, one 4-set, and one 5-set;
   its missing-edge graph is a 5-cycle, and \((g_2,g_3)=(5,5)\).

These are ordinary combinatorial polytopes, not claimed Coxeter examples.
They are useful positive fixtures for the next nerve-extension stage.

## Correct Gale rank and dataset audit

The appropriate oriented-matroid rank is **four**, not five.  A homogeneous
vertex matrix of \(Q\) has rank \(9+1=10\), so its Gale kernel has dimension

\[
 14-10=4.
\]

After choosing a generic affine chart one may reorient individual Gale
vectors and draw a signed/colored configuration in \(\mathbb R^3\).  The
signs are essential: because the Gale oriented matroid is totally cyclic,
forgetting them does not give an equivalent ordinary affine point order type.

No complete public dataset at the required size was located.  The Homepage
of Oriented Matroids reaches uniform rank-four oriented matroids on nine
elements; that class already has 9,276,601 reorientation classes.  The
complete realizability work of Fukuda--Miyata--Moriyama treats rank four on
eight elements, not fourteen.  Firsching's public polytope data consist of
other small or neighborly classes, and polyDB does not advertise an exhaustive
collection of simplicial 9-polytopes with 14 vertices.

Consequently, downloading a complete nerve list is not a viable route.

## Recommended complete next-stage generator

There are two rigorous routes, best used as cross-checks.

### 1. Missing-face hypergraph extension

For each of the 444 surviving \(H\), enumerate antichains

\[
 \mathcal M=E(H)\cup\mathcal A_4\cup\mathcal A_5,
\]

where every member of \(\mathcal A_4\) and \(\mathcal A_5\) is independent
in \(H\), and no member contains another missing face.  The complex consists
of subsets containing no member of \(\mathcal M\).  During canonical
backtracking impose:

- every maximal face has size nine;
- every ridge lies in exactly two facets;
- all links have the integral homology of the appropriate sphere;
- the whole complex is a PL 8-sphere; and
- isomorphic hypergraphs are rejected under \(\operatorname{Aut}(H)\).

This finite generator is complete for abstract candidate nerves under the
proved missing-face restrictions.  Polytopality is a separate final filter.

### 2. Rank-four chirotope generator

A uniform rank-four chirotope on 14 elements has
\(\binom{14}{4}=1001\) sign variables.  Its Grassmann--Plücker axioms can be
encoded by 45,045 six-element constraints.  Generate single-element
extensions with canonical augmentation as in Finschi--Fukuda, retain the
appropriate totally cyclic Gale orientations, reconstruct positive circuits
and their nerves, and apply the missing-edge manifest as an early filter.
Finally decide realizability exactly; abstract oriented-matroid feasibility
alone does not prove polytopality.

The hypergraph route should be attempted first because it begins with only
444 graph shadows and applies the strongest Coxeter-specific information
before oriented-matroid realizability.  The chirotope route is the independent
coverage check required before a nonexistence claim.

## Primary and computational references

- Amanda Burcroff, *Near Classification of Compact Hyperbolic Coxeter
  \(d\)-Polytopes with \(d+4\) Facets and Related Dimension Bounds*,
  [arXiv:2201.03437](https://arxiv.org/abs/2201.03437), especially
  Proposition 3.5, Proposition 4.1, and Lemma 10.10.
- Richard Stanley, *The number of faces of a simplicial convex polytope*,
  Advances in Mathematics 35 (1980), 236--238, for necessity in the
  \(g\)-theorem.
- David Barnette and Gil Kalai's lower-bound/equality theory; a convenient
  exact statement is Theorem 1.1 of Hailun Zheng,
  [arXiv:1607.05757](https://arxiv.org/abs/1607.05757).
- Lukas Finschi and Komei Fukuda, *Generation of Oriented Matroids---A Graph
  Theoretical Approach*, Discrete & Computational Geometry 27 (2002),
  117--136, [DOI 10.1007/s00454-001-0056-5](https://doi.org/10.1007/s00454-001-0056-5).
- Komei Fukuda, Hiroyuki Miyata, and Sonoko Moriyama, *Complete Enumeration
  of Small Realizable Oriented Matroids*,
  [arXiv:1204.0645](https://arxiv.org/abs/1204.0645).
- Brendan McKay and Adolfo Piperno, nauty and Traces; the archived generation
  used nauty 2.9.3's `geng`.

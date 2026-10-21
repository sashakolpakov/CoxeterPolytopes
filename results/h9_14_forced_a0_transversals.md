# Forced-\(a=0\) minimal-transversal obstruction for H9/14

## Result and exact scope

Among the 142 missing-edge graph shadows surviving the strict generalized
lower-bound and ridge-support layer, exactly 36 force

\[
a=|\mathcal A_4|=0.
\]

Dehn--Sommerville then forces the number \(b=|\mathcal A_5|\) of missing
five-faces:

- 35 graphs have \(b=1\);
- one graph has \(b=0\).

There are exactly 30,934 admissible \(\mathcal A_5\) choices at this stage.
For every one, the archived certificate supplies an inclusion-minimal
transversal of the resulting minimal-nonface hypergraph whose size is three
or four rather than five. Hence none can be the boundary of a pure
simplicial \(9\)-polytope.

This exactly excludes all 36 forced-\(a=0\) graph shadows. Only one of them
lies among the corrected 30 local-link survivors, so the combined graph
frontier is

\[
30-1=29.
\]

At this stage the 29 were unresolved graph shadows. The first global
completion snapshot found two ordinary-polytopal product selections and left
27 theorem-undecided. Subsequent exact certificates, including the record-34
exclusion, and four multiwedge realizations give the current global partition
(121) impossible, (6) ordinary-polytopal, and (15) theorem-undecided; see
[`h9_14_facet_completion.md`](h9_14_facet_completion.md).

The implementation and certificate are:

- coxeter_polytopes.combinatorics.forced_a0_transversals;
- scripts/h9_14_forced_a0_transversals.py;
- data/enumerations/h9_14_forced_a0_transversals.json; and
- tests/test_h9_14_forced_a0_transversals.py.

All archived source paths are repository-relative.

## 1. Maximal faces and minimal transversals

Let \(\mathcal M\) be the hypergraph of minimal nonfaces of a simplicial
complex \(K\) on vertex set \(V\). A set \(F\subseteq V\) is a face exactly
when it contains no member of \(\mathcal M\). Equivalently, its complement

\[
B=V\setminus F
\]

meets every member of \(\mathcal M\), so \(B\) is a transversal.

**Lemma 1.** The map \(F\mapsto V\setminus F\) is a bijection between
maximal faces of \(K\) and inclusion-minimal transversals of \(\mathcal M\).

**Proof.** Suppose \(F\) is a face and \(B=V\setminus F\). If \(B\) is not
an inclusion-minimal transversal, then \(B\setminus\{v\}\) is still a
transversal for some \(v\in B\). Its complement \(F\cup\{v\}\) is a face,
so \(F\) is not maximal. The converse is the same argument reversed.
\(\square\)

The target boundary is pure of dimension eight, so every maximal face has
nine vertices. Since \(|V|=14\), Lemma 1 gives the exact necessary and
sufficient purity condition

\[
\boxed{\text{every inclusion-minimal transversal of }\mathcal M
       \text{ has size }14-9=5.}
\]

One inclusion-minimal transversal of any other size therefore certifies
failure of purity in the target dimension.

## 2. Re-deriving the forced stratum

For each of the 142 source graphs, the verifier recomputes the strict
missing-four interval directly from \(H\):

\[
\max\{0,\beta(H)-g_3^{\langle3\rangle}\}
\leq a\leq
\min\{I_4(H),\beta(H)-1\}.
\]

Exactly 36 intervals are the singleton \(\{0\}\). Their distribution is:

| missing edges | forced-\(a=0\) graphs |
|---:|---:|
| 4 | 1 |
| 6 | 2 |
| 7 | 5 |
| 8 | 11 |
| 9 | 17 |
| **total** | **36** |

This count is derived from the graph6 strings and is not read from the new
certificate.

## 3. Dehn--Sommerville forces \(b\)

When \(a=0\), every independent four-set is a face:

\[
f_3=I_4(H).
\]

The no-missing-triples reduction also gives

\[
f_1=\binom{14}{2}-|E(H)|,
\qquad
f_2=I_3(H).
\]

The relation \(h_5=h_4\) fixes

\[
f_4=6f_3-21f_2+56f_1-126\cdot14+252.
\]

With no missing four-face, an independent five-set is either a four-face
of \(K\) or a minimal missing five-face. Therefore

\[
\boxed{b=I_5(H)-f_4.}
\]

Direct evaluation gives \(b=1\) for 35 graphs and \(b=0\) for the graph

\[
\texttt{M???????C?G?G?C??}.
\]

## 4. The complete \(\mathcal A_5\) choice space

Under the proved restrictions, all minimal nonfaces have orders two, four,
or five. In the forced-\(a=0\) stratum,

\[
\mathcal M=E(H)\cup\mathcal A_5.
\]

If \(b=0\), there is exactly one choice, \(\mathcal A_5=\varnothing\).

If \(b=1\), its unique member must be an independent five-set of \(H\);
otherwise it would contain a missing edge and would not be minimal. Conversely,
every independent five-set is retained at this stage. Hence the complete
choice space is exactly

\[
\{\,\{A\}:A\subseteq V,\ |A|=5,\ A\text{ independent in }H\,\}.
\]

Across the 35 graphs this contains 30,933 choices. Including the one empty
choice gives 30,934.

## 5. Explicit certificates

For every choice, the archive stores a pair

\[
(A,T),
\]

where \(A=0\) is the sentinel for the unique empty \(\mathcal A_5\) choice
and \(T\) is a vertex mask. The independent verifier checks directly:

1. the archived \(A\)-masks equal the complete independently enumerated
   choice space in canonical order;
2. \(T\) meets every edge of \(H\) and also \(A\), when \(A\neq0\);
3. for every \(t\in T\), the set \(T\setminus\{t\}\) fails to meet at least
   one hyperedge, proving inclusion-minimality; and
4. \(|T|\neq5\).

No SAT result or search trace is trusted by the verifier.

The exact witness-size distribution is:

| bad transversal size | choices |
|---:|---:|
| 3 | 3,192 |
| 4 | 27,742 |
| **total** | **30,934** |

A size-three transversal has an eleven-vertex complementary maximal face,
and a size-four transversal has a ten-vertex complementary maximal face.
Both contradict a pure eight-dimensional boundary, whose facets must have
nine vertices.

## 6. Combination with the corrected local-link layer

The corrected local-link theorem, including the fact that simplex face
figures are forbidden, leaves 30 graph shadows. Exactly one is in the
forced-\(a=0\) set, namely the unique \(b=0\) graph above. Therefore the
two independent necessary filters leave

\[
\boxed{29}
\]

graph shadows.

This is not yet an exhaustive nerve nonexistence result. The remaining
graphs still require simultaneous missing-four and missing-five placement,
facet-support selection, spherical-link and polytopality tests, Coxeter
labels, and exact Lorentzian Gram feasibility.

# Polytopal-link constraints for the H9/14 nerve search

## Result and exact scope

Let \(K=\partial Q\), where \(Q\) would be the simplicial polar of a compact
hyperbolic Coxeter \(9\)-polytope with \(14\) facets, and let \(H\) be the graph
of missing edges of \(K\). The preceding exact filters leave 142 possible
graphs \(H\).

The Lanner order bound is inherited by every link. Applying it together with
the generalized lower-bound equality theorem to the polytopal face figures
gives the exact reduction

\[
142\longrightarrow72\longrightarrow39\longrightarrow36\longrightarrow30.
\]

| missing edges | input | vertex \(g_1>0\) | strict vertex \(g_2\) | vertex \(g_3\) interval | edge/triangle feasibility |
|---:|---:|---:|---:|---:|---:|
| 2 | 2 | 2 | 2 | 2 | 2 |
| 3 | 5 | 5 | 5 | 5 | 5 |
| 4 | 11 | 10 | 9 | 8 | 7 |
| 5 | 17 | 15 | 10 | 9 | 8 |
| 6 | 29 | 22 | 9 | 8 | 6 |
| 7 | 33 | 14 | 3 | 3 | 2 |
| 8 | 28 | 4 | 1 | 1 | 0 |
| 9 | 17 | 0 | 0 | 0 | 0 |
| **total** | **142** | **72** | **39** | **36** | **30** |

The last two stages mean only that every displayed linear incidence interval
is individually nonempty. They do not prove that one simultaneous missing-four
family realizes all vertex, edge, and triangle incidences. Thus the 30
survivors remain graph shadows at this local-link layer, not yet nerves or
polytopes. Later global work produces two product nerves and four exact
multiwedge nerves, while the record-23/34/37/54/74/92/93/46 certificates give the
current partition (121) impossible, (6) ordinary-polytopal, and (15)
theorem-undecided; see
[`h9_14_facet_completion.md`](h9_14_facet_completion.md).

The exact implementation, archive, and tests are:

- coxeter_polytopes.combinatorics.local_link;
- scripts/h9_14_local_link_filter.py;
- data/enumerations/h9_14_local_links.json; and
- tests/test_h9_14_local_links.py.

The archive freezes the ordered graph transcript, every per-vertex invariant,
the count vectors above, and its canonical SHA-256 digest.

## 1. Inheritance of the Lanner bound by links

**Lemma 1.** Let \(K\) be a simplicial complex, let \(F\in K\), and let \(T\)
be a minimal nonface of \(\operatorname{lk}_K(F)\). There is a minimal
nonface \(M\) of \(K\) such that

\[
M\subseteq F\cup T,
\qquad
M\setminus F=T.
\]

In particular, \(|T|\leq|M|\). If \(F=\{v\}\), then \(M=T\) or
\(M=T\cup\{v\}\).

**Proof.** Since \(F\cup T\notin K\), choose an inclusion-minimal nonface
\(M\subseteq F\cup T\). Because \(F\in K\), the set \(M\setminus F\) is
nonempty. If some \(t\in T\) did not belong to \(M\), then

\[
M\subseteq F\cup(T\setminus\{t\}),
\]

contradicting the minimality of \(T\) in the link. Hence \(T\subseteq M\),
so \(M\setminus F=T\). \(\square\)

Every minimal nonface of the target Coxeter nerve supports a Lanner diagram
and therefore has order at most five. Lemma 1 proves

\[
\boxed{\text{every minimal nonface of every link has order at most five.}}
\]

## 2. There is no simplex exception

If \(F\) has \(s\) vertices, its link is the boundary of a simplicial
face figure of dimension

\[
D=9-s.
\]

A \(D\)-simplex boundary has one minimal nonface: its full set of \(D+1\)
vertices. For each relevant link dimension \(D=8,7,6\), this order exceeds
five. Consequently no such face figure can be a simplex, and

\[
\boxed{g_1=f_0-D-1>0.}
\]

Equivalently, in the one-simplex triangulation of a \(D\)-simplex, the top
\(D\)-simplex is an interior face: it is not a boundary face, while all its
proper faces are boundary faces. Thus the generalized-lower-bound
minimal-interior-face argument also includes this case.

This corrects the earlier draft, which incorrectly exempted simplex face
figures.

## 3. Strict higher local generalized lower bounds

Suppose \(g_r=h_r-h_{r-1}=0\), where \(1\leq r\leq D/2\). Murai--Nevo's
equality theorem gives a geometric \((r-1)\)-stacked triangulation with no
interior face of dimension at most \(D-r\). Choose an inclusion-minimal
interior simplex. All its proper faces lie in the boundary, so its vertex set
is a missing face of the link. Its order is at least

\[
D-r+2.
\]

Whenever \(D-r+2>5\), this contradicts Lemma 1. Hence

\[
\boxed{
g_r>0
\quad\text{for}\quad
1\leq r\leq \min(\lfloor D/2\rfloor,D-4).
}
\]

For \(r=1\), the direct simplex-boundary argument above is the clearest proof.

| linked face size \(s\) | dimension \(D\) | strict entries | missing orders forced by equality |
|---:|---:|---|---|
| 0 | 9 | \(g_1,g_2,g_3,g_4\) | \(10,9,8,7\) |
| 1 | 8 | \(g_1,g_2,g_3,g_4\) | \(9,8,7,6\) |
| 2 | 7 | \(g_1,g_2,g_3\) | \(8,7,6\) |
| 3 | 6 | \(g_1,g_2\) | \(7,6\) |

For the new vertex filter:

- \(D=8,r=1\) excludes a nine-vertex simplex figure;
- \(D=8,r=2\) forces a missing face of order at least eight;
- \(D=8,r=3\) forces a missing face of order at least seven.

## 4. Exact vertex-figure formulas

Fix a vertex \(v\). Let

\[
W_v=\{u\neq v:uv\notin E(H)\},
\qquad n_v=|W_v|,
\]

let \(e_v=|E(H[W_v])|\), and let \(I_{3,v}\) count independent triples of
\(H[W_v]\). The no-missing-triples theorem gives

\[
f_0(v)=n_v,
\qquad
f_1(v)=\binom{n_v}{2}-e_v.
\]

If \(a_v\) is the number of missing four-faces containing \(v\), then

\[
f_2(v)=I_{3,v}-a_v.
\]

For the simplicial \(8\)-dimensional figure,

\[
\begin{aligned}
g_1(v)&=n_v-9,\\
g_2(v)&=f_1(v)-8n_v+36
       =\binom{n_v-8}{2}-e_v,\\
g_3(v)&=f_2(v)-7f_1(v)+28n_v-84
       =\beta_v-a_v,
\end{aligned}
\]

where

\[
\beta_v=I_{3,v}-7f_1(v)+28n_v-84.
\]

The exact conditions are

\[
g_1(v)>0,
\qquad
1\leq g_2(v)\leq g_1(v)^{\langle1\rangle},
\]

and

\[
\boxed{
\max\{0,\beta_v-g_2(v)^{\langle2\rangle}\}
\leq a_v\leq\beta_v-1.
}
\]

Since \(n_v=13-\deg_H(v)\), strict \(g_1\) immediately gives

\[
\deg_H(v)\leq3.
\]

The vertex \(g_1\), \(g_2\), and \(g_3\)-interval stages leave 72, 39, and
36 graphs, respectively.

## 5. Edge and triangle constraints

### Edge figures

Let \(E=\{p,q\}\) be an independent pair of \(H\), hence an edge of \(K\), and
put

\[
W_E=\{u\notin E:up,uq\notin E(H)\},
\qquad n_E=|W_E|.
\]

Let \(I_{2,E}\) count independent pairs in \(H[W_E]\), and let \(a_E\) count
missing four-faces containing \(E\). The figure has dimension seven and

\[
g_1(E)=n_E-8,
\qquad
g_2(E)=I_{2,E}-a_E-7n_E+28=\gamma_E-a_E.
\]

It must satisfy

\[
g_1(E)>0
\]

and

\[
\boxed{
\max\{0,\gamma_E-g_1(E)^{\langle1\rangle}\}
\leq a_E\leq\gamma_E-1.
}
\]

There is no \(g_1=0\) branch.

### Triangle figures

Every independent triple \(T\) is a face. Let \(q_T\) be the number of
vertices extending it to an independent four-set, and let \(a_T\) count
missing four-faces containing it. Its six-dimensional figure has
\(q_T-a_T\) vertices. Strict \(g_1>0\) gives

\[
\boxed{a_T\leq q_T-8.}
\]

Applying all edge intervals and triangle capacities to the 36 vertex survivors
leaves 30 graph shadows.

The incidence checks

\[
\sum_v a_v=4a,
\qquad
\sum_{\substack{E\in K\\|E|=2}}a_E=6a,
\qquad
\sum_{\substack{T\in K\\|T|=3}}a_T=4a,
\]

hold for \(a=|\mathcal A_4|\).

## 6. Constraints once all face indicators are available

The same theorem should be imposed inside the facet solver:

- at every vertex,
  \[
  1\leq g_4(v)
  =f_3(v)-6f_2(v)+21f_1(v)-56f_0(v)+126
  \leq g_3(v)^{\langle3\rangle};
  \]
- at every edge,
  \[
  1\leq g_3(E)
  =f_2(E)-6f_1(E)+21f_0(E)-56
  \leq g_2(E)^{\langle2\rangle};
  \]
- at every independent triple,
  \[
  1\leq g_2(T)=f_1(T)-6f_0(T)+21.
  \]

There are no simplex equality branches in these dimensions.

## 7. Further exact topological and polytopal tests

For every face \(F\in K\),

\[
\sum_{\substack{U\cap F=\varnothing\\F\cup U\in K}}(-1)^{|U|}
=(-1)^{9-|F|}.
\]

This local Euler equation is necessary but not sufficient for spherical links.

Facets containing \(F\) give the facet-adjacency graph of its link. If
\(D=9-|F|\), this graph is \(D\)-regular by the ridge equations and
\(D\)-vertex-connected by Balinski's theorem. These local cut constraints are
stronger than one global connectivity test.

All links can be checked over \(\mathbb Z\) by Smith normal forms, or first
over finite fields by exact boundary-matrix ranks. Boundary complexes of
convex polytopes and all their links are shellable; only an exhaustive
no-shelling certificate may be used for rejection.

Over any field, a homology \(8\)-sphere on 14 vertices has a
Stanley--Reisner ring of Gorenstein codimension five, with

\[
\beta_{i,j}=\beta_{5-i,14-j},
\qquad
\beta_{5,14}=1.
\]

Hochster's formula and the identity between the alternating Betti polynomial
and \(h_K(t)(1-t)^5\) provide an independent missing-face and homology
checksum.

## 8. Rank-four Gale route

The linear Gale rank is

\[
14-(9+1)=4.
\]

After a general-position perturbation preserving a simplicial face lattice,
one may use a uniform rank-four chirotope. It has
\(\binom{14}{4}=1001\) signs. The three-term Grassmann--Plücker sign axioms
give

\[
\binom{14}{2}\binom{12}{4}=45{,}045
\]

small constraints. For a five-set \(B=\{b_0<\cdots<b_4\}\), its circuit signs
are

\[
c_i=(-1)^i\chi(B\setminus\{b_i\}).
\]

The complement of \(B\) is a facet exactly when all five \(c_i\) agree.
Oriented-matroid feasibility remains weaker than exact realization.

Existing complete oriented-matroid catalogues are far below rank four on
fourteen elements, so no public dataset found in this audit replaces the
required generator.

## References

- Satoshi Murai and Eran Nevo, *On the generalized lower bound conjecture for
  polytopes and spheres*, Acta Math. 210 (2013), Theorem 1.2,
  <https://arxiv.org/abs/1203.1720>.
- Amanda Burcroff, *Near Classification of Compact Hyperbolic Coxeter
  \(d\)-Polytopes with \(d+4\) Facets and Related Dimension Bounds*,
  <https://arxiv.org/abs/2201.03437>.
- Komei Fukuda, Hiroyuki Miyata, and Sonoko Moriyama, *Complete Enumeration
  of Small Realizable Oriented Matroids*,
  <https://arxiv.org/abs/1204.0645>.

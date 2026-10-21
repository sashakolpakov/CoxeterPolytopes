# One-ideal simplex blocks at a product cusp

There is a cusp counterpart to the compact-simplex disjoint-port theorem,
but its capacity is larger and its scope excludes triangles with two or
three ideal vertices. The following gives a proved incidence constraint
for triangles with **exactly one** ideal vertex. It supplies a possible
ingredient for the original-face budget in equation(11) of the
[finite-volume report](fresh_finite_volume_simplex_truncation_20261017.md),
not a completed budget or a new dimension bound.

Let P be a finite-sided, finite-volume acute Coxeter n-polytope. Fix an
ideal vertex v with link

\[
 L_v=\prod_{i=1}^{s}\Delta^{r_i},\qquad
 r_i\ge1,\quad \sum_i r_i=n-1,\quad
 V=\prod_i(r_i+1).
\]

A port at v means an incident old edge of P, hence a vertex of L_v.
A **factor fiber** is a maximal simplex face of L_v obtained by varying
one factor and fixing a vertex in every other factor. A fiber of factor i
has r_i+1 vertices. There are V/(r_i+1) such fibers.

Let A range over maximal actual simplex faces of P of dimension at least
two whose sole ideal vertex is v. Maximality is within this specified
class. Then every A belongs to a unique factor fiber F, and

\[
 \boxed{\quad\sum_{A\text{ assigned to }F}\dim A\le r_i+1.\quad}
\tag{1}
\]

The port sets of the summands in (1) are disjoint. Consequently

\[
 \sum_A\dim A\le sV,
 \qquad
 N_{3,1}(v)=\sum_A\binom{\dim A}{2},
\tag{2}
\]

where N_{3,1}(v) counts triangular two-faces with v as their only ideal
vertex. More precisely, a fixed old edge through v belongs to at most s
of these maximal simplex faces, at most one for each factor.

## Proof and the finite-vertex condition

A simplex face of a product of simplices, of dimension at least one,
varies in exactly one factor. A face varying in two factors contains a
square and cannot be a simplex. Thus every triangle through v, and every
simplex of dimension at least two containing it, determines a unique
factor fiber of the link.

The fiber F = Delta^{r_i} is the link at v of an actual (r_i+1)-face P_F
of P. Although v may be nonsimple in P, it is simple in P_F. Its ports
within P_F are exactly the r_i+1 vertices of F. All finite vertices of
P_F are simple as well.

Relate two ports only if their two-face is a triangle with v as its sole
ideal vertex. In two such triangles sharing an edge from v, both the
other endpoint of the shared edge and the two remaining triangle vertices
are finite. The generated three-face has degree three at v, because v is
simple in P_F, and degree three at the other endpoint, because it is
finite. Deleting the other two triangle vertices would separate those
endpoints from any further vertex. Balinski's three-connectivity excludes
such a vertex. This three-face is a tetrahedron with v as its only ideal
vertex. The port relation is transitive.

For a class of k >= 2 ports, their neighbors are finite. The triangles
give the complete graph on v and these k neighbors in the generated
k-face. Every displayed vertex is simple in that face and already has
degree k. Connectivity excludes any further vertex, finite or ideal.
The face is a k-simplex whose only ideal vertex is v. Conversely, every
such simplex gives a clique. Maximal classes therefore give the unique
maximal simplex extensions within P_F, and their port sets are disjoint.

A larger one-ideal simplex in P would have a simplex link containing the
positive-dimensional link of the original face. Its factor fiber must
therefore be the same F. Maximality within P_F is also maximality in P.
This proves (1) and the unique assignment of A. Every one-ideal triangle
belongs to exactly one block, and a k-simplex has binom(k,2) triangles
through v. This proves the second identity in (2). Summing (1) over the
V/(r_i+1) fibers for each i proves its first inequality.

The proof only uses simplicity at v inside P_F and at the displayed
finite neighbors. Other vertices of P_F may be nonsimple ideals. A port
ending at another cusp cannot participate in these triangles and is not
counted as a one-ideal simplex aggregate. Isolated ports are omitted.

The graph-connectivity input is
[Balinski's original theorem, printed p.432](https://msp.org/pjm/1961/11-2/pjm-v11-n2-p03-s.pdf).
The product-link and finite-vertex simplicity hypotheses are recorded in
[Potyagailo–Vinberg, printed p.64](https://ems.press/content/serial-article-files/43032).

## Weighted capacity and its limitation

For arbitrary nonnegative port weights u(F,e), disjointness gives the
stronger inequality

\[
 \sum_{A\text{ assigned to }F}\ \sum_{e\subset A,\ e\ni v}u(F,e)
 \le\sum_{e\in\operatorname{Vert}(F)}u(F,e).
\tag{3}
\]

It is thus legitimate to credit each one-ideal simplex block by its
weighted port incidence and debit the corresponding cusp factor fibers.
If weights depend only on factor i and are constant on its ports, say
gamma_i >= 0, the total debit is exactly V sum_i gamma_i. The unweighted
special case is sV. These are actual incidence constraints; no optimization
or assumption about realizability of an abstract block distribution is
being made.

The number of original two-faces through v is
E_v=(n-1)V/2. In the dimension23–26 target range the necessary benchmark in
the earlier report is alpha_n E_v, with alpha_n equal to9/46 or1/5.
Equation(3) does not prove that the one-ideal triangle demand can be paid
with gamma_i satisfying V sum_i gamma_i <= alpha_n E_v. A complete charge
must also account for the finite corners of those triangles, compact
faces, and triangles or quadrilaterals with several ideal corners.

For example, assigning the entire demand2 of each one-ideal triangle to
its cusp gives only the general upper estimate

\[
 2N_{3,1}(v)
 =\sum_A (\dim A)(\dim A-1)
 \le V\sum_i r_i=(n-1)V=2E_v.
\tag{4}
\]

This is much larger than alpha_n E_v and cannot close the global budget.
The gain in (1)–(3) is the constraint on **simultaneous block incidence**,
which can be combined with other payments; it is not a free replacement
of2 by alpha_n.

## Sharp local example and scope boundary

The exact Coxeter square pyramid in the
[previous certificate](fresh_finite_volume_simplex_truncation_20261017.json)
has link Delta1 x Delta1 at its apex. Its four side triangles each have
exactly one ideal vertex. They are the four maximal blocks, each of
dimension2. Thus sum dim A=8=sV and N_{3,1}=4. Each edge from the apex is
in two blocks, one for each factor. A cusp disjoint-port capacity1 would
already be false in this literal finite-volume example.

Combinatorially, a pyramid over any product link saturates (1) in every
fiber, but no Coxeter realization of these higher-dimensional pyramids
is asserted. The exact square pyramid is sufficient to show why the
factor multiplicity cannot simply be discarded.

**FINAL/FROZEN.** This is a written geometric proof and a precise candidate
incidence constraint. Exact_hub_counts independently checked the unique
factor-fiber assignment, relative simplicity, completion steps, counts and
finite-ended-port scope. This is an internal conceptual review, not a new
arithmetic replay, external review or dimension exclusion. All predecessor
files remain unchanged.

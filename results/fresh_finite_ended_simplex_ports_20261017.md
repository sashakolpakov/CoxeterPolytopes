# Finite-ended ports in a chamber with ideal vertices

Let P be a finite-sided, finite-volume acute Coxeter polytope. At a finite
vertex v, call an incident edge a finite-ended port if its other endpoint
is finite. Triangular adjacency need not be transitive on all ports.
However, the closed neighborhood of every finite-ended port is a clique,
and that clique generates a whole simplex face, possibly with several
ideal vertices.

Consequently distinct maximal simplex faces having at least two finite
vertices cannot share a finite-ended port at v. Every triangular face with
at least two finite vertices belongs to a unique such maximal simplex.
This extends the compact-face result while keeping a different scope from
the [one-ideal cusp blocks](fresh_cusp_one_ideal_simplex_capacity_20261017.md).

## Proof

Suppose the triangles v,a,b and v,a,c share the finite edge va; b and c
may be ideal. The three ports at v generate a three-face because v is
simple. In this three-face v and a each have degree three, and their
neighbors are already exactly the other three displayed vertices. If any
additional vertex existed, deleting b and c would disconnect v,a from it.
Balinski's three-connectivity excludes this. The three-face is a
tetrahedron, so the ports to b and c are also triangularly adjacent.
Thus all neighbors of the finite-ended port va are mutually adjacent.

More generally, let a clique of d ports contain a finite-ended port. In
their generated d-face F, the displayed vertices are v and its d
neighbors. They form a complete graph. Let q of these neighbors be ideal.
Since v and at least one neighbor are finite, q <= d-1. Every displayed
finite vertex is simple in F and already has all its d edges within the
displayed graph. If an additional vertex existed, deleting the q displayed
ideal vertices would disconnect the remaining displayed finite vertices
from it. This contradicts d-connectivity. For q=0 ordinary connectivity
already suffices. Hence F has exactly d+1 vertices and is a simplex.

The closed neighborhood of a finite-ended port is therefore a simplex.
Every simplex containing that port has its port set inside this closed
neighborhood. It is the unique maximal simplex containing the edge. If
two maximal simplex faces through v shared any finite-ended port, both
would be this same face. A triangle with at least two finite vertices has
a finite edge, so it has a unique maximal extension as asserted.

Only the displayed finite vertices are assumed simple. The ideal vertices
may be nonsimple in P or F before graph exhaustion. The connectivity input
is [Balinski's original theorem, printed p.432](https://msp.org/pjm/1961/11-2/pjm-v11-n2-p03-s.pdf),
applied to the ordinary convex projective closure of F.

## Valid capacities and face counts

Let f(A) and h(A) be the finite and ideal vertex counts of a selected
maximal simplex A, with f(A) >= 2. At each finite vertex v,

\[
 \sum_{A\ni v}\bigl(f(A)-1\bigr)\le\deg_{ff}(v),
\tag{1}
\]

where deg_ff(v) is the number of edges from v to other finite vertices.
Arbitrary nonnegative weights on these finite-ended ports give the
corresponding weighted inequality. No capacity1 is asserted for an
ideal-ended port; distinct blocks may share those ports.

If isolated finite edges are included as dimension-one maximal simplex
blocks, every finite edge belongs to exactly one block. Writing E_ff for
the total number of finite-to-finite edges gives

\[
 \sum_A\binom{f(A)}2=E_{ff}.
\tag{2}
\]

Without the dimension-one blocks the same formula is an inequality <=.
Every triangle with at least two finite vertices is counted once, so

\[
 N_{3,0}+N_{3,1}
 =\sum_A\left[\binom{f(A)}3+h(A)\binom{f(A)}2\right].
\tag{3}
\]

Here N_{3,s} counts triangles with exactly s ideal vertices. The terms
with two or three ideal vertices are deliberately absent: they do not
have a finite edge, and the uniqueness argument does not cover them.
Even if such a triangle lies in one of the displayed maximal simplexes,
it must not be counted in every containing block without a separate
uniqueness or allocation argument.

These are combinatorial constraints on actual faces. They do not identify
an induced simplex face with a Coxeter simplex, do not supply a dimension
bound, and do not establish the complete noncompact-face charge needed
for an ambient dimension theorem.

**FINAL/FROZEN.** Root proposed the finite-ended-port refinement;
frontier_literature independently checked both graph-connectivity steps,
the ideal-vertex deletion count q<d, and the capacity scope. This is a
written internal proof review, with no arithmetic replay, optimization,
external review or formal proof-assistant claim. Predecessor files were
not edited.

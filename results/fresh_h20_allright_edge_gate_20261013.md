# Every section crossing an all-right source edge is non-Coxeter

Every hyperplane section of the frozen Borcherds H21 chamber that crosses
an open edge incident to either all-right finite vertex is **non-Coxeter**.
The normal can have arbitrary real coordinates. This extends the earlier
[all-right-vertex obstruction](fresh_h20_borcherds_sections_20261013.md):
the section need not contain either old vertex, the center, or a cusp.
It does not exclude sections avoiding these 42 open edges.

The [standard-library exact checker](../scripts/fresh_h20_allright_edge_gate_20261013.py)
and its [certificate](fresh_h20_allright_edge_gate_20261013.json) verify
both vertices, all 42 first edges, all 420 auxiliary cusp-to-cusp edges,
and 8,862 support-avoidance cases. The proof below supplies the continuous
argument, with no coefficient grid or bound on Coxeter labels. A separate
agent reviewed the geometric argument, including facet survival and
transversality at the second edge. A
[separate Fraction reconstruction](../scripts/fresh_h20_allright_edge_gate_20261013_review.py)
also [passes](fresh_h20_allright_edge_gate_20261013_review.json); it imports
no producer or repository geometry helper and rechecks all saved incidence
inputs directly against the integer source.

## Sparse pairings at a new finite vertex

More generally, let an interior section cross an ambient open edge with
positive-definite Gram matrix G of its 20 incident normals. The section
has a finite vertex there. If its normal is e, the projected Gram is

    H = G - aa^T,       a_i = (e,r_i)/sqrt(e^2).

It is positive definite: the section is transverse to the edge, so its
20 incident projected normals remain independent. Every one is a true
facet near this vertex. A finite Coxeter diagram is a forest.

Three a_i of the same nonzero sign would give a negative triangle in H.
If all three angles were Coxeter, its normalized off-diagonal entries
would be at most -1/2, contradicting positive definiteness. Opposite-sign
nonzero entries must lie at adjacent vertices of the old diagram, since
an old zero entry would project to a positive, obtuse one. Two entries
of each sign would therefore require an old four-cycle. Consequently
at most three of the 20 pairings can be nonzero. A support of size three
has two signs and forms a two-edge star in the old forest; at least one
old star edge must become zero to avoid a new negative triangle.

This reduces a general search through an ambient edge to at most five
linear dimensions for e: three active pairing coordinates and the
two-dimensional Lorentzian edge plane. It is a reduction, not a complete
classification of such sections.

For an all-right edge, the old graph is empty. Opposite signs are already
impossible, so **at most two** incident pairings can be nonzero.

## Exact source geometry

The 42 short roots split into the two point/line parts of PG(2,4). Fix
one part of 21 mutually orthogonal norm-two roots and write v for its
all-right vertex, scaled so that

    v^2 = -2,       (v,r_p) = 0 for every point p,
    (v,r_l) = -1 for every line l.

Each line meets exactly five points, with pairing -1, and is orthogonal
to the other point roots. Put delta_p=v-r_p. These are future null rays
satisfying all 210 inequalities. The open edge from v to delta_p has
exactly the other 20 point walls as its vanishing inequalities.

For distinct points p,k, their unique common line l gives another actual
ambient edge. Its endpoints are delta_p and delta_k. Their common zero
set is exactly the other 19 point walls together with l. Their Gram is
D4 plus 16 A1 components, and

    (delta_p,delta_k) = -2.

The checker verifies these statements separately for both point/line
choices and every indicated pair. Thus for every 0<u<1,

    x = u*delta_p+(1-u)*delta_k
      = v-u*r_p-(1-u)*r_k,       x^2=-4u(1-u)<0

is in the relative interior of the asserted ambient edge. There are no
unexamined inequalities hidden in this incidence claim.

## The unavoidable forbidden ridge

Suppose e-perpendicular crosses the first edge at

    P = v-u*r_p,       0<u<1.

If the section were Coxeter, the first local argument leaves at most two
active point roots S among the other 20. The 21 point roots and v are an
orthogonal basis of the ambient vector space. Therefore

    e = c*r_p-c*u*v + sum_{a in S} b_a*r_a.

Here c is nonzero, since otherwise the entire first edge is contained in
the section. That case is covered by the earlier all-right-vertex theorem.

There are five lines through p. At most two meet S, so choose a line l
through p meeting none of S and a point k on l distinct from p. Then
(e,r_k)=0 and the second edge point x above also belongs to the section.
The crossing is transverse because

    (e,delta_p-delta_k) = -2c != 0.

This is another genuine finite section vertex, with all 20 projected
incident facets present. In particular l and the three other points of
l, apart from p,k, give a positive-definite projected D4 block. Those
three point roots are all orthogonal to e. But

    (e,r_l) = -c*(1-u) != 0.

If N is the squared norm of the projected line normal, its block has
three diagonal entries 2, all three central pairings -1, and center N.
Its Schur complement and the nonzero pairing give

    3/2 < N = 2-(e,r_l)^2/e^2 < 2.

Each central edge is an actual dihedral ridge of the section, with

    1/4 < cos^2(theta) = 1/(2*N) < 1/3.

No Coxeter angle lies in that interval: the values for labels 3 and 4
are 1/4 and 1/2, with higher labels larger still. This contradiction
excludes every real transverse section in the declared class.

## Replay and scope

    python scripts/fresh_h20_allright_edge_gate_20261013.py
    python scripts/fresh_h20_allright_edge_gate_20261013_review.py

All genuine interior sections inherit finite volume from the source
Klein polytope; the obstruction here is a forbidden actual angle. No
new finite-volume Coxeter example or dimension is established. The
sparse-pairing reduction leaves other source-edge orbits as a concrete
construction route, rather than another enlargement of an unbounded
wall dictionary.

# Global affine exclusion and the all-simple dimension bound61

Let P be a finite-sided finite-volume hyperbolic Coxeter polytope in H^n.
If every ideal vertex of P is simple, then **n<=61**. Finite vertices of
an acute-angled polytope are automatically simple. The statement permits
both finite and ideal vertices and requires no arithmetic hypothesis.

The key structural fact is stronger than a restriction on black deletions:
**every connected affine subdiagram of an all-simple P has exactly n
nodes**. Consequently, when n>=10, every minimal indefinite subdiagram
of order at most10 is Lanner, and every actual simplex face has dimension
at most4. Supporting white roots, including affine subdiagrams containing
every black root, are covered by this global argument.

The [checker](../scripts/fresh_simple_cusp_affine_exclusion_20261017.py)
and [saved certificate](fresh_simple_cusp_affine_exclusion_20261017.json)
give the small exact radius, graph-cost, and threshold checks. This is a
fresh successor to the frozen [distance15 proof](fresh_simple_cusp_distance_bound_20261017.md).
Its compact-facet conclusion n<=127 remains unchanged. The unrestricted
finite-volume bound163 remains unchanged: arbitrary ideal vertices need
not be simple. No novelty or current-best published-bound claim is made.

## An affine diagram must be a whole simple cusp

Use outward spacelike normals e_i and a future timelike strict interior
point H, so H.e_i<0. Distinct normals have nonpositive products.
Let C be a connected affine principal subdiagram. Its Gram matrix has a
strictly positive kernel vector (c_i). Define

    z = sum_{i in C} c_i e_i.

Then z^2=0, z.e_i=0 for i in C, and z.e_j<=0 for every other wall j.
Moreover H.z<0, so z is nonzero and lies on the future null sheet.
Thus [z] is a feasible point at infinity of P.

The projective closure of a finite-volume polytope is a finite convex
polytope contained in the closed Klein ball. Every point of this closure
on the boundary sphere is extreme: strict convexity of the ball forbids
expressing it as a nontrivial convex combination of two distinct points
of the closed ball. Hence [z] is an actual ideal vertex, not merely a
feasible direction unrelated to the vertex set. Every wall of C contains
this vertex.

At an ideal vertex the Euclidean link is a product of simplices. A simple
vertex has exactly n incident walls, so its (n-1)-dimensional link is one
simplex; its entire Coxeter diagram is connected affine of order n.
Every proper principal subdiagram of this connected affine diagram is
positive definite. Since it contains C, which is singular, C must be the
whole diagram. In particular,

    order(C) = n.                                      (1)

The primary [Potyagailo--Vinberg paper, printed pp.64 and69](https://ems.press/content/serial-article-files/43032)
states the finite/ideal local face structure and the finite-volume
projective convex-hull convention used here. The positive-kernel and
strict-convexity argument above supplies the particular global conclusion.
There is also direct primary confirmation in
[Felikson--Tumarkin, Section1, PDF pp.3--4](https://arxiv.org/pdf/math/0502413):
every parabolic subdiagram lies in an ideal-vertex diagram, and a simple
ideal vertex has a connected parabolic diagram. Their displayed
condition(iii), that parabolic diagrams have order n, uses these simplicity
facts; it does not use the all-ideal assumption needed for their separate
condition(ii) excluding elliptic diagrams of order n. We use only the
former implication in the mixed-vertex class.

## Consequences for small indefinite stars

An inclusion-minimal indefinite principal diagram L is Lanner or
quasi-Lanner: if a proper principal subdiagram is affine, take one
connected affine component C. Any outside node coupled to its positive
kernel makes C plus that node indefinite. Minimality forces L to be
exactly C plus that node. Thus a quasi-Lanner diagram of order at most10
contains a connected affine diagram of order at most9.

For n>=10 this contradicts (1). Therefore every minimal indefinite
diagram of order at most10 is Lanner and has order at most5. These
classification inputs are the same ones used in the frozen
[original-star simplex proof](fresh_multi_ideal_triangle_scope_20261017.md).
The checker additionally reads all26 high-rank quasi-Lanner diagrams and
confirms that their listed affine deletions are connected and of order
at most9; it does not repeat their classification.

In an actual k-simplex face, the original star has n-k supporting white
roots and k+1 black roots. Every indefinite principal subdiagram contains
all black roots, because deleting any black leaves normals at an actual
finite or ideal vertex. The original-star proof gives a minimal L of
order at most10. Hence for n>=10,

    k+1 <= order(L) <= 5,  so k<=4.                    (2)

In particular actual simplex faces of dimensions7,8,9 are all excluded,
not just the dimension9 case. No assumption that an induced face is a
Coxeter simplex is made. Equation (1) also handles the formerly delicate
white affine deletion containing every black node; that affine diagram
produces its own actual cusp elsewhere in P and is too small.

## Triangle and quadrilateral payments at distance7

Retain the original-wall-star and graph-distance conventions of the
distance15 proof. Put weight1 at a polygon corner when its two side
normals are joined by a path of length at most7 within that corner's
incident diagram, and weight0 otherwise.

For a triangle, a minimal indefinite star contains its three black nodes.
It is now Lanner of order at most5. The tree joining the three marked
nodes has at least two marked-to-marked paths whose interiors avoid the
third black node. Each has length at most order(L)-2<=3. These paths lie
in actual corner diagrams, so every triangle receives at least2.

For a quadrilateral, each opposite pair together with its white support
is indefinite. Otherwise the positive-definite common-face argument or
the nonzero feasible affine-kernel argument would put a finite or ideal
point in two disjoint CLOSED opposite edges. This remains impossible for
quadrilaterals with ideal vertices. Each resulting minimal indefinite
diagram L or M contains both opposite black nodes and has order at most5.

In a connected N-node graph with two marked nodes, the distance of any
node to the nearer mark is at most N-2: a shortest path to their set has
no marked interior nodes. Hence this radius is at most3 for L and M.
The two indefinite diagrams cannot be disjoint and orthogonal, because
that would give two negative directions in a Lorentz space. They either
overlap or have a joining edge. Joining their paths and taking a shortest
path between their two black sets gives length at most

    3+1+3 = 7.                                        (3)

Its interior has only white nodes; its endpoints are consecutive sides,
so it lies in an actual corner diagram. Every quadrilateral receives
at least1. The checker confirms the radius bound for every connected
simple graph of orders2 through5 and every marked pair, rather than
assuming a particular list of low-order Coxeter diagrams.

## The vertex cost is at most7n, in every rank

At a finite vertex the n incident normals form a finite Coxeter forest.
At a simple ideal vertex they form a connected affine diagram on n nodes.
For paths, cycles, single forks and double forks, the all-rank formulas
proved in the preceding distance15 report give at most N unordered pairs
at each positive fixed graph distance. Therefore the total number at
distance at most7 is at most7N. Summing over finite components gives7n.

The exceptional connected finite/affine diagrams have at most9 nodes,
so even all their unordered pairs number at most4N<=7N. This includes
the affine D4 four-leaf star, separated from the double-fork formula.
No finite rank cutoff is used. The checker reconfirms the new cost for
all245 independently counted classical graphs and the seven saved
exceptional graphs; the inherited analytic formulas establish all ranks.

Thus the total polygon-corner weight is at most7n f_0(P), while the
triangle/quadrilateral payments give at least2f_2^{(3)}+f_2^{(4)}.

## Strict h-vector threshold

The projective closure of P is simple. The same strict simple-polytope
h-vector inequality used in the frozen distance proof gives

    2 f_2^{(3)} + f_2^{(4)} > b_n f_0(P),
    b_n = n(n-6)/8                 if n is even,
          (n-1)(n-5)/8             if n is odd.

Consequently b_n<7n. For n=62+2t and n=63+2t, respectively,

    8(b_n-7n) = 4t^2+124t,
                 4t^2+128t+68.

Both are nonnegative for every integer t>=0. Equality at n=62 still
contradicts the STRICT lower inequality. Therefore n<=61. Dimensions
n<=9, where the preceding Lanner reduction was not invoked, already
satisfy this conclusion.

## Verification and limits

Generation and a fresh saved replay pass: 7,521 marked low-order graph
checks, all26 high-rank affine-deletion checks, and245 inherited classical
graph-cost checks. Syntax, nine direct source hashes, and report links
pass. The [independent mathematical/source review](fresh_simple_cusp_affine_exclusion_20261017_review.md)
finds no gap and explicitly does not claim a duplicate arithmetic replay.
The checker uses exact integers/Fractions and frozen source manifests.
The polynomial identities and graph formulas, not the
finite confirmation ranges, justify the unbounded dimension and rank
claims. No optimizer, arithmetic field restriction, or induced-face
Coxeter assumption enters the proof.

This lemma supplies neither an unrestricted bound61 nor a construction.
The mixed-vertex all-simple class differs from the simple ALL-IDEAL
class, for which [Felikson--Tumarkin](https://arxiv.org/pdf/math/0502413)
prove the stronger dimension bound8. No transfer of that all-ideal result
to the mixed case is used here.

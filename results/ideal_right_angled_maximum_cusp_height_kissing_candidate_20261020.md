# Maximum cusp height, kissing configurations, and a candidate ideal-dimension obstruction

20 October 2026. Written mathematical candidate, awaiting independent
adversarial review. No computation, Lean, census or build was run.
Do not promote a dimension bound from this note before that review.
The proposed argument is short enough that each geometric reduction
and its equality case are included below.

## 1. Precise proposed conclusion and external inputs

**Candidate theorem.** There is no entirely ideal finite-volume
right-angled hyperbolic five-polytope. Consequently entirely ideal
right-angled finite-volume polytopes have dimension at most four.

The two substantive external theorems used are:

1. A subset of the unit sphere in Euclidean four-space whose distinct
   vectors have pairwise product at most one half has at most 24 points.
   This is the main theorem on printed page 4 of O. R. Musin,
   [The kissing number in four dimensions, Annals of Mathematics 168
   (2008), 1–32](https://annals.math.princeton.edu/wp-content/uploads/annals-v168-n1-p01.pdf).
   Only its cardinality bound is used; uniqueness of the kissing
   configuration is neither asserted nor needed.
2. An entirely ideal right-angled hyperbolic four-polytope has at least
   24 facets, and equality characterizes the regular ideal right-angled
   24-cell. This is Theorem 2 of A. Kolpakov,
   [On the optimality of the ideal right-angled 24-cell, Algebraic &
   Geometric Topology 12 (2012), 1941–1960](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-s.pdf).

The other inputs are the usual all-wall acuteness statement for an
acute-angled polytope, inheritance of right angles and finite volume
by faces of dimension at least three, and connectedness of the facet
adjacency graph of a convex polytope. Inheritance is used only for
four-facets, not for ideal polygon angles.

## 2. Height and an actual point of the maximizing facet

Use Lorentz signature `(-,+,+,+,+,+)`. All facet normals are outward
unit spacelike vectors, with inequalities `<x,m> <= 0`. Assume that
`P` is an entirely ideal right-angled finite-volume five-polytope.
Choose one of its ideal vertices and a future null representative `v`.
For each actual facet normal `m`, put

    H(m) = -<v,m> >= 0.

The facets incident to the chosen cusp have height zero. At least one
facet has positive height: otherwise all facet hyperplanes contain
that same ideal point, which is incompatible with a finite-volume
full-dimensional polytope. Equivalently, an ideal vertex of a
full-dimensional polytope is not incident to every facet.

There are finitely many facets, so choose a unit normal `n` with
maximum height `H0 > 0`. Define

    q = (v + H0*n)/H0.

Its square is `-1`, and `<q,n>=0`. It is future because its product
with the future null `v` is negative. For ANY actual facet normal `m`,
all-wall acuteness gives `<n,m> <= 0` when `m != n`; hence

    <q,m> = -H(m)/H0 + <n,m> <= 0.

For `m=n` the product is zero. Thus `q` is an ACTUAL point of the
four-facet `Q=P intersect n-perp`. No inference from visible pairwise
wall intersections is used. The point may initially lie on the
boundary of `Q`; the argument below allows height-zero neighboring
walls and does not assume that `q` is interior.

## 3. The spherical code of neighboring walls

Let `m_1,...,m_f` be exactly the actual facets adjacent to `n` across
its ridges. Their restrictions are precisely the facets of `Q`.
Each has `<m_i,n>=0`, so it is also an intrinsic unit normal in the
Lorentz five-space `n-perp`. In coordinates with `q=(1,0,0,0,0)`,
write

    m_i = (t_i, sqrt(1+t_i^2)*u_i),
    t_i = -<q,m_i> = H(m_i)/H0 in [0,1],
    u_i in S^3.

For two distinct neighbors, their ambient inner product equals their
intrinsic inner product and is at most zero. Therefore

    u_i.u_j <= t_i*t_j / sqrt((1+t_i^2)*(1+t_j^2)) <= 1/2.

The last bound follows from `0<=t_i,t_j<=1`. The directions are
distinct: equal directions would give positive Lorentz product for
the two unit normals. Thus Musin's theorem applies to an honest
`f`-point spherical code and gives `f<=24`.

The facet `Q` is an entirely ideal finite-volume right-angled
four-polytope. Kolpakov's theorem gives `f>=24`, and equality makes
`Q` isometric to the regular ideal 24-cell.

## 4. The equality point is forced to be the center

This step uses the METRIC equality conclusion for `Q`, not any
unproved uniqueness of a spherical kissing arrangement.

In the standard Lorentz model of the regular ideal right-angled
24-cell, its 24 outward unit normals are

    (1, epsilon*e_i + delta*e_j),
    1<=i<j<=4, epsilon,delta in {+1,-1}.

Each has square one, and the cell center is `c=(1,0,0,0,0)`.
Opposite normals have the form `(1,r)` and `(1,-r)`. Transport the
point `q` to this model and write it `(a,w)`, where `a>0` and
`a^2-|w|^2=1`. The already proved inequalities `t_i<=1` apply to
both members of every opposite pair:

    a-r.w <= 1,    a+r.w <= 1.

Adding gives `a<=1`. Future unit timelikeness gives `a>=1`, so
`a=1` and `w=0`. The point is the center, and consequently

    t_i=1,    H(m_i)=H0

for EVERY facet adjacent to `n`.

## 5. Global propagation and the proposed contradiction

Let `M` be the set of actual facets with height `H0`. It is nonempty.
The preceding proof can be repeated at any member of `M`, since each
still has globally maximum positive height. It proves that all of
its facet-adjacency neighbors also lie in `M`.

The facet-adjacency graph is connected, so every facet of `P` lies
in `M`. But the facets incident to the originally chosen cusp have
height zero, whereas `H0>0`. This is the proposed contradiction.

If the five-dimensional conclusion survives independent review,
face inheritance would exclude every higher entirely ideal dimension:
an entirely ideal right-angled polytope of dimension at least six
has an entirely ideal finite-volume right-angled five-face.

## 6. Adversarial checklist and status

The cardinality bound uses `S^3` in Euclidean FOUR-space; it is not
the three-dimensional kissing number. The regular-cell equality
uses Kolpakov's metric characterization, not a uniqueness claim for
four-dimensional kissing configurations. Zero-height neighbors are
allowed until the equality-center argument rules them out. Unit
normal scaling is fixed throughout, so the maximum-height choice
and all ratios are comparable. Every facet of `Q` comes from an
actual orthogonal ridge neighbor, and the point `q` is checked
against ALL walls. Neither a small set of feasible vertices nor a
bounded local search is being used to infer finite volume.

The conclusion remains a CANDIDATE pending independent checking of
these steps. No repository status table, manuscript, frozen proof
packet or verification claim was changed.

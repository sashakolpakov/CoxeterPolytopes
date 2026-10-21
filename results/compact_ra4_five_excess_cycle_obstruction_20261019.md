# Compact right-angled H4: the five-excess cycle is impossible

19 October2026. New written proof, awaiting independent review. No
scientific calculation, census, optimizer, build or Lean process ran.
Earlier source-bound evidence and the compact-capacity manuscript remain
unchanged. The extra metric input is stated explicitly below.

## Statement and inherited inputs

For a compact right-angled four-polytope P put

    W=sum_(ridges R)(number_of_sides(R)-5).

Then

    W=0 or W>=6.                                    (1)

Together with the reviewed zero-excess recognition, W=0 is exactly the
regular right-angled120-cell. This improves the positive-excess gap5
to6. It does not prove the120-cell minimum or double-volume gap, or
change a dimension ceiling. The previous46/46/96 facet-size cutoffs are
unchanged.

Use the already reviewed facts that a compact right-angled3-polyhedron
with a k-gon has at least2k+2 faces, and equality has Loebell combinatorics
with two disjoint k-gons and2k pentagons. The positive-excess graph has
one vertex per positive-excess facet and one edge of weight k-5 per
nonpentagonal shared ridge. Its weighted degree at F is d(F), its total
edge weight is W, it has no degree-one vertex, and every degree-two
vertex has equal edge weights and the two disjoint Loebell bases.
Consequently its two neighbors do not meet even in the full facet graph.
The earlier proof already excludes0<W<5.

Source pins:

* `compact_ra4_neighbor_excess_and_zero_excess_rigidity_20261019.md`:
  `9c670a7b3c25dda48997c3f70a6d1be9be88bcdefb475b7710f5011068aebcbe`.
* `review_compact_ra4_neighbor_excess_root_20261019.md`:
  `292187dbf4fb5ea327baed6086177ec1b0fac174116ca3c579a6d75bc43f0baa`.

One further external input is the UNIQUENESS clause of the compact
nonobtuse three-dimensional Andreev theorem: a fixed marked combinatorial
type and fixed dihedral angles have at most one compact hyperbolic
realization up to isometry. See Roeder--Hubbard--Dunbar,
[Theorem1.4, p832](https://aif.centre-mersenne.org/article/AIF_2007__57_3_825_0.pdf#page=9).
Only this metric-rigidity clause is needed; no list of four-polytopes or
published minimum-volume claim is used. The Loebell combinatorial
identification is supplied by the earlier own closing proof.

## 1. W=5 forces exactly five Loebell L6 facets

Suppose W=5. Since there is no degree-one vertex, the positive graph has
minimum degree two. It cannot have at most four edges, by the previous
argument, so it has exactly five edges, all of weight1. It is connected:
two components with minimum degree two would already require six edges.

There are at most five vertices. With four vertices the graph would be
K4 minus an edge, and each degree-two vertex has two adjacent neighbors,
contrary to their proved disjointness in P. Three vertices cannot carry
five distinct edges. Thus there are five vertices, all of degree two,
forming C5. Denote the corresponding actual facets by F_i, i modulo5.

Each F_i has exactly two nonpentagonal ridges, both hexagons. Its excess
is2 and it has14 faces. Equality in the3D capacity lemma makes it the
Loebell L6 combinatorial type: two opposite hexagons, with two belts of
six pentagons. Its hexagonal bases are precisely F_i intersect F_(i-1)
and F_i intersect F_(i+1). No assertion about an arbitrary abstract
five-cycle is substituted for these actual facet incidences.

## 2. Own normal calculation of the L6 base distance

Let Q be ANY compact right-angled realization of this L6 type, with
outward unit base normals a,b. We prove

    a.b=-sqrt(3).                                    (2)

The combinatorial rotation taking each belt pentagon to the next and
preserving both bases has order six. By the stated marked metric
uniqueness it is realized by an isometry of Q. It is a genuine60-degree
rotation around a geodesic axis: it preserves the inward side of each
base and acts on the base hexagon by a one-step orientation-preserving
rotation. Both base planes are perpendicular to that axis. The axis
meets the interior of each base, since it fixes the center of its
rotationally symmetric compact hexagon. Their center segment lies in Q.

In Lorentz coordinates(-t^2+x^2+y^2+z^2), place the midpoint of this
segment at(1,0,0,0). Write its half-length as h>0, and put
c=cosh(h),s=sinh(h). Then

    a=(s,0,0,c),       b=(s,0,0,-c).

Label upper pentagon normals u_i and lower normals v_i so that u_i is
orthogonal to a,u_(i-1),u_(i+1),v_(i-1),v_i. These are the Loebell
incidences. Choose angles theta_i=i*pi/3. Rotational symmetry gives

    u_i=(T,R*cos(theta_i),R*sin(theta_i),Z),
    v_i=(T1,R1*cos(theta_i+psi),R1*sin(theta_i+psi),Z1).

The positive radial coefficients R,R1 are nonzero: otherwise a belt
normal would be fixed by the rotation. Orthogonality to the base gives
Z=T*s/c and Z1=-T1*s/c. The midpoint is strict on belt walls, giving
T,T1>0. The square-one and adjacent-upper orthogonality equations are

    R^2-T^2/c^2=1,
    R^2/2-T^2/c^2=0.

Therefore R=sqrt(2),T=c,Z=s. Applying the identical two equations to
the lower belt gives R1=sqrt(2),T1=c,Z1=-s. No assumed equality of
upper/lower edge lengths or arbitrary prescribed metric has been used.

Now u_0.v_0=u_1.v_0=0 gives

    2*cos(psi)=c^2+s^2=2*cos(psi-pi/3).

Thus psi=pi/6 modulo pi. The common value c^2+s^2 is positive, so the
negative-cosine branch is excluded. It follows that

    c^2+s^2=2*cos(pi/6)=sqrt(3),
    a.b=-(c^2+s^2)=-sqrt(3),

proving(2). In particular, the metric assertion follows from an actual
rigidity premise and explicit equations, rather than from the name L6
alone. The same equations provide a consistent isolated symmetric L6
normal model; the obstruction below occurs when five such facets must
fit in a single four-dimensional Lorentz space.

## 3. The actual ambient Gram matrix has three negative directions

Let n_i be the outward unit normal of F_i in P. Consecutive facets share
a hexagonal ridge, so n_i.n_(i+1)=0. Inside F_i, the neighboring normals
n_(i-1),n_(i+1) already lie in n_i-perp and are the outward unit normals
of its two bases. There is no rescaling or extra projection. Equation(2)
therefore gives

    n_(i-1).n_(i+1)=-sqrt(3).

All ten off-diagonal entries of the five-normal Gram matrix are now
specified: adjacent products0, nonadjacent products-sqrt(3). The
complement of C5 is again C5, whose adjacency eigenvalues are
2,1/phi,1/phi,-phi,-phi, where phi=(1+sqrt(5))/2. Hence the Gram
eigenvalues are

    1-2sqrt(3),
    1-sqrt(3)/phi       (multiplicity2),
    1+phi*sqrt(3)       (multiplicity2).

Since phi^2=(3+sqrt(5))/2<3, the first three eigenvalues are negative.
An actual normal Gram matrix in Lorentz signature(4,1) has negative
index at most one. This contradiction excludes W=5 and proves(1).

The same proof forbids a connected C5 component of five unit-weight
edges in the positive-excess graph even if other positive components
are present. It does not exclude all degree-two cycles: their lengths
or equal edge weights can differ, and the L6 metric computation cannot
be transferred to larger Loebell bases.

## What remains open in this packet

A distinct polytope at or below the120-cell volume must have W>=6 and
facets with at most46 vertices. Excluding W=5 does not force f+W>=120.
No rigidity for larger positive-excess graphs, global120-cell
tessellation theorem, compact minimum/equality theorem beyond W=0,
or compact second-volume theorem follows. Independent review of this
successor proof and later formalization are separate obligations.

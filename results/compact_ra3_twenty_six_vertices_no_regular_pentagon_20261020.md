# Twenty-six vertices exclude a regular pentagonal face

20 October2026. Written proof awaiting independent review. No enumeration,
scientific execution, numerical realization, manuscript or Lean change.
The classification of the underlying fullerene is an explicit published
input; its hyperbolic edge calculation below is derived here.

**Theorem.** A compact right-angled hyperbolic three-polyhedron with26
vertices has no entire regular right-angled pentagonal face.

Consequently a positive facet adjacent to a dodecahedral facet in a
compact RA4 polytope has excess at least four. In particular every
compact RA4 polytope has at least61 facets and volume at least47pi^2/6.
These are lower bounds, not a universal120-cell minimum or second-volume
claim. No worldwide novelty claim is made.

## 1. The precise classification input and an explicit model

A compact RA3 polyhedron is trivalent, so26 vertices give39 edges and15
faces. Its polygon capacity is F>=2k+2 for a k-gonal face. Thus it has
no face of size at least seven. Euler then forces twelve pentagons and
three hexagons: its one-skeleton is a fullerene in the usual graph sense.

The only classification imported is that there is ONE fullerene graph
on26 vertices. The primary account is Brinkmann, Goedgebeur and McKay,
*The Generation of Fullerenes*, arXiv:1207.7010v2 (16October2012),
Table2, PDF page15, row(nv,nf)=(26,15), column“fullerenes”=1:

https://arxiv.org/pdf/1207.7010v2 .

The paper defines fullerenes as cubic plane graphs with only pentagonal
and hexagonal faces; it also explains why abstract and embedded graph
isomorphism coincide here. We do not import a chemical bond geometry
or infer any edge lengths from the abstract symmetry group.

Here is a full dual triangulation, so no external graph diagram or
symmetry identification is needed. All indices are modulo three. Its
vertices are the face labels

    A_i, A'_i, B_i, B'_i, H_i,       i=0,1,2.

Its triangles are

    A0 A1 A2,
    A_i A_(i+1) B_i,
    A_i B_i H_i,
    A_i H_i B_(i-1),

    A'0 A'1 A'2,
    A'_i A'_(i+1) B'_i,
    A'_i B'_i H_i,
    A'_i H_i B'_(i-1),

    B_i B'_i H_i,
    B_i B'_i H_(i+1).

There are26 triangles and39 edges. Every edge lies in two triangles;
the complex is connected. The links of A_i,B_i,H_i are respectively
the cycles

    A_(i-1), A_(i+1), B_i, H_i, B_(i-1);
    A_i, A_(i+1), H_(i+1), B'_i, H_i;
    A_i, B_i, B'_i, A'_i, B'_(i-1), B_(i-1),

and the primed A/B links follow by interchange of primes. Thus this
is a connected closed triangulated2-manifold of Euler characteristic
15-39+26=2, hence a sphere. Its face-adjacency lists are exactly

    A_i: A_(i-1), A_(i+1), B_i, B_(i-1), H_i;
    B_i: A_i, A_(i+1), H_i, H_(i+1), B'_i;
    A'_i and B'_i: the primed versions of these lists;
    H_i: A_i, A'_i, B_i, B_(i-1), B'_i, B'_(i-1).

Thus the twelve A/B vertices have degree five, and the three H vertices
have degree six. This is a simplicial sphere, whose dual is a cubic
polyhedral graph on26 vertices with the stated face sizes. The imported
uniqueness identifies the dual of any polyhedron in the theorem with
this explicit triangulation.

## 2. Why the geometric symmetry normalization is legitimate

The displayed triangulation has the automorphisms

    rho: every index i maps to i+1;
    sigma: interchange A_i with A'_i, B_i with B'_i, fix H_i;
    tau: A_i maps to A_-i, H_i to H_-i, B_i to B_(-i-1),
         with the same rule on primed labels.

They generate the group D3 times C2. The usual compact3D Andreev metric
uniqueness realizes every such combinatorial automorphism isometrically
on an ASSUMED actual polyhedron Q. This is a conditional use of
uniqueness: we do not obtain compactness or actual incidence from a
formal list of normals. If no such Q existed the theorem would already
hold. Given Q, its actual edge and vertex incidences are exactly the
listed dual edges and triangles.

This finite group fixes an interior point O. For example sum the future
unit representatives of the orbit of an interior point and normalize;
the result is future timelike and still satisfies every wall strictly.
Use Lorentz coordinates(-t^2+x^2+y^2+z^2) centered at O. The order-three
rho acts by a120-degree rotation around the z-axis. Since sigma fixes
each of the three distinct H faces, their horizontal projections are
nonzero and sigma acts identically on the horizontal plane. It must
therefore negate z, interchanging the upper and lower faces. The
involution tau is a vertical reflection fixing the horizontal direction
of H0. Its z action is positive: A0 is fixed by tau and is distinct
from its sigma-image A'0. Set theta_i=2pi*i/3 and e(theta)=(cos theta,
sin theta). The actual OUTWARD UNIT normals therefore have the form

    H_i  =(h,R e(theta_i),0),
    A_i  =(a,U e(theta_i),A),    A'_i=(a,U e(theta_i),-A),
    B_i  =(b,V e(theta_i+pi/3),B),
    B'_i =(b,V e(theta_i+pi/3),-B).

Here h,a,b are strictly positive because O is interior. Choose R>0.
Orthogonality of A_i,H_i and of B_i,H_i forces U,V>0. Choose A>0
by interchanging upper and lower; the A_i,B_i orthogonality below
then forces B>0. The phases of B_i follow from the vertical reflection
rho*tau fixing B_0, followed by the rho rotations. This accounts for possible signs, rather than
silently choosing an inward normal or a wrong orbit alignment.

## 3. Exact orthogonality equations

Put X=h^2. The norm equation for H_i gives R^2=1+X. The edges A_i H_i
and B_i H_i give

    U R=a h,             V R=2b h.                    (1)

The top triangle contains the edge A_i A_(i+1), so those normals are
orthogonal. Combined with their unit norm this gives

    A^2=a^2+U^2/2,       U^2=2/3,
    a^2=2(1+X)/(3X),     A^2=(2+3X)/(3X).              (2)

The edge B_i B'_i similarly gives

    B^2=V^2-b^2,         B^2=1/2,
    b^2=(1+X)/(2(3X-1)).                               (3)

In particular X>1/3. For the orthogonal pair A_i,B_i, equation(1)
implies U V/2=ab*X/(1+X), hence

    A B=ab/(1+X).

Square this positive equality and substitute(2),(3). It becomes

    (2+3X)/(6X)=1/(3X(3X-1)),
    (2+3X)(3X-1)=2,
    9X^2+3X-4=0.

The positive solution is

    X=(sqrt(17)-1)/6.                                  (4)

No numerical approximation or square-root sign choice is needed.

## 4. One nonregular edge in each pentagon orbit

For an actual compact RA3 edge, its two defining unit facet normals
are orthogonal, and the two third normals at its endpoints belong to
the edge's Lorentz two-space. Their product is -cosh(edge length).
This is the usual direct Lorentz formula for the distance between the
two endpoint hyperplanes in that two-space.

The dual edge A_i A_(i+1) has endpoint triangles A0 A1 A2 and
A_i A_(i+1) B_i. Its endpoint third normals are A_(i-1) and B_i,
whose horizontal directions are opposite. Equations(1)--(3) give

    <A_(i-1),B_i>=-ab-U V+A B=-3ab*X/(1+X),
    cosh^2(length(A_i intersect A_(i+1)))
        =3X/(3X-1)=(7+sqrt(17))/4.                    (5)

This is greater than11/4 since sqrt(17)>4, whereas
phi^2=(3+sqrt5)/2<11/4 since sqrt5<5/2. Thus every A pentagon, upper
or lower, has an edge STRICTLY LONGER than the regular pentagon edge.

The dual edge B_i B'_i has endpoint third normals H_i,H_(i+1).
Their product gives

    cosh(length(B_i intersect B'_i))
        =X+(1+X)/2=(1+sqrt(17))/4.                    (6)

This is less than3/2 since sqrt(17)<5, whereas phi>3/2. Thus every B
pentagon has an edge STRICTLY SHORTER than the regular pentagon edge.
The two orbits exhaust all twelve pentagons. This proves the theorem.

## 5. The unconditional RA4 consequence and the remaining target

Let Q be a positive facet adjacent to an actual D facet in compact RA4.
Their common ridge is an ENTIRE regular pentagon in Q. Excess two is
impossible by the reviewed L6 metric, whose lateral pentagons contain
hexagonal base edges b with cosh(b)=2 rather than phi. Excess three would give
V(Q)=26 and is impossible by the theorem. Therefore d(Q)>=4 and
V(Q)>=28. A D exists by the reviewed farthest-facet theorem. Unless
the parent is the all-D120-cell, its connected facet graph contains
such a D--positive adjacency. Facet capacity gives f>=2*28+5=61.
The120-cell has120 facets and satisfies that bound as well. Combining
f>=61, the nonzero-excess gap W>=17 and the exact volume formula gives

    Vol(P)>=pi^2*(61+17+16)/12=47pi^2/6.

In the global charge notation, this replaces d(Q)>=3 by d(Q)>=4 for
every boundary-positive facet. The old source is held unchanged pending
review; no new120-cell minimum, recognition of all positive facets, or
classification of all compact RA4 polytopes follows.

## Other unchanged proof inputs and verification boundary

* `manuscripts/compact_ra4_volume_constraints.tex`, SHA256
  `afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f`:
  compact incidence and simplicity; polygon and facet capacities;
  compact3D metric uniqueness; regular pentagon length cosh L=phi;
  L6 recognition; W=0 or W>=17 and all-D recognition; volume formula.
* `results/compact_ra4_farthest_dodecahedron_wall_radius_20261020.md`,
  SHA256 `d082e4a918d4504a59ac200f32bbb141f1f2cb52d701006eaea7ba5b0e5bcc22`:
  existence of a dodecahedral facet for every compact RA4 polytope.

The sole finite classification input is the cited primary fullerene
count. The graph, actual symmetry normalization and both edge formulas
are proved above. No local or global numerical checker has been run,
and no claim of a self-contained fullerene enumeration is made.

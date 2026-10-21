# A global low-excess charge in compact right-angled dimension four

20 October 2026. Written proof, awaiting independent review. No searches,
execution, manuscript change or formal verification. This gives unconditional
count and volume inequalities, followed by a precise remaining homology
branch. It does not prove the 120-cell minimum or its second volume.

Let P be compact and right-angled in H4. Write f,v for its facet and vertex
counts, W=sum_R(k(R)-5)=v-5f, and

    d(Q)=(V(Q)-20)/2,       sum_Q d(Q)=2W.

A zero-excess facet is the regular dodecahedron D. Write N2 for the number
of facets of excess two, N+ for the number of positive facets, and Nu,Ni
for the numbers of uniform and irregular D facets. Uniform means that all
twenty outgoing edges have the regular length L, cosh L=phi.

## 1. The excess-two facets have no dodecahedral neighbors

Every excess-two facet is the right-angled Loebell L6. Its two nonpentagonal
faces are regular hexagons, each with edge length b satisfying cosh b=2.
Its twelve lateral pentagons each contain such a base edge. These lateral
pentagons are not regular pentagons of length L, since 2 is not phi.
Consequently all fourteen parent neighbors of an L6 are positive: the
two hexagonal neighbors because their common ridge is nonpentagonal,
and the twelve others because a face of D would be a regular pentagon.
These are properties of the entire actual faces, not formal face types.

There is always a D facet, by the reviewed farthest-dodecahedron theorem.
If P is not the all-D 120-cell, connectedness of its facet graph therefore
gives a positive facet Q adjacent to D. This Q cannot have excess two by
the preceding paragraph. Since every positive facet has excess at least
two, it has d(Q)>=3 and V(Q)>=26. The reviewed facet capacity gives

    f>=2V(Q)+5>=57.                                      (1)

The 120-cell has120 facets, so (1) holds for every compact RA4 polytope.
If P is not the 120-cell, the reviewed W gap also gives W>=17. The exact
volume formula consequently yields

    Vol(P)=pi^2(f+W+16)/12 >=15pi^2/2.                   (2)

The regular120-cell also satisfies (2). No equality recognition or claim
of worldwide novelty is attached to this lower bound.

## 2. Charge the whole interior positive core and its interface

Partition the positive facets into I, those with no D neighbor, and B,
those with at least one D neighbor. These names refer to facet adjacency,
not membership in the interior or boundary of the ambient polytope.
Put m=|I|, b=|B| and

    C=sum_(Q in I)(d(Q)-1),       EB=sum_(Q in B)(d(Q)-3).

Every excess-two facet belongs to I by Section1. All facets in B have
excess at least three; hence C>=N2 and EB>=0.

Both parent facets of a nonpentagonal ridge are positive. Write WII,
WIB,WBB for the total ridge weights according to these two parent
types. Thus W=WII+WIB+WBB, and exact incidence counting gives

    TI=sum_I d=2WII+WIB,
    TB=sum_B d=2WBB+WIB=3b+EB.                          (3)

A nonpentagonal ridge incident to I has NO dodecahedral edge-carrier.
Such a carrier would meet its parent I facet along a face, making D
a neighbor of that facet. Actual simplicity ensures this adjacency;
a mere formal orthogonality is not being substituted for incidence.

The reviewed irregular-edge argument says that every irregular D has
at least three nonpentagonal perpendicular ridges. A k-gonal ridge has
at most floor(2k/3)<=4(k-5) dodecahedral edge-carriers, by the actual
three-consecutive-L-edge closure lemma. Only B--B ridges can contribute
here, so

    3Ni<=4WBB.                                         (4)

Use f=Nu+Ni+m+b and (3) to obtain

    Nu >= f-m-b-(4/3)WBB
        = f-2W+C+(2/3)WIB+EB/3.                       (5)

This is an unconditional global charge: the excess of every interior
positive facet, the nonpentagonal interface weight, and the excess
above three of boundary positive facets enter the SAME inequality.
It does not add overlapping incidence allowances as separate budgets.
Dropping its nonnegative terms and using C>=N2 retains the useful
corollary

    Nu >= f-2W+N2.                                      (6)

In particular f>=57 implies Nu>=1 whenever W<=28. Both (5) and (6)
are valid at every W; they do not assume W<20 or the absence of
excess-two facets. The case I empty and all B facets of excess three
has C=WIB=EB=0: the refined charge does not penalize this possible
boundary-only separator beyond its original W contribution.

## 3. What this charges at or below the120-cell volume

The following uses the reviewed center-cut theorem, solely here. Let K
be the simplicial dual boundary3-sphere of P, and let K+ be its full
subcomplex on positive facets. Its complementary full subcomplex KD
is on the dodecahedral facets. For finite full subcomplexes with these
complementary vertex sets, S3 minus |K+| deformation retracts to |KD|:
in barycentric coordinates discard the positive coordinates and
renormalize the remaining ones. Alexander duality therefore gives

    reduced H2(K+;F2) is isomorphic to
    reduced H0(KD;F2).                                  (7)

Only connected components of KD enter this identity. No assertion that
K+ is two-dimensional or is a surface is made.

Adjacent D facets have the same matched120-cell center. Here is the
explicit calculation. Let nF,nG be their outward normals, and p the
center of their common regular pentagon. With

    h=sqrt(phi/2), B=sqrt(1+h^2),
    T0=phi^(3/2), B0=sqrt(1+T0^2)=sqrt(2)*phi,

their intrinsic centers satisfy oF=Bp-h*nG and oG=Bp-h*nF. The matched
centers qF=B0*oF-T0*nF and qG=B0*oG-T0*nG coincide because B0*h=T0.
Every wall adjacent to either D has height T0 at this common point.

If KD is connected and its common center q lies in P's interior, then
every D wall has height T0 at q. By the farthest-dodecahedron theorem,
the maximum of ALL actual wall heights is consequently T0. A positive
facet adjacent to D would also have height T0, contradicting the same
theorem, which makes EVERY maximizing facet D. Connectedness of the
whole facet graph now forces P to be the all-D120-cell.

Suppose Vol(P)<=Vol(120-cell), and P is different from that cell. A
uniform D has its matched center in P's interior: otherwise an actual
wall reaching or cutting that center would force Vol(P)>=23pi^2/2,
strictly greater than Vol(120-cell)=34pi^2/3. Thus if H2(K+)=0, (7)
and the preceding argument imply Nu=0. The full charge (5) gives

    f+C+(2/3)WIB+EB/3<=2W.

Its corollary (6), combined with (1), gives

    f+N2<=2W,       W>=29,       f+W>=86.                (8)

In particular a non120-cell at or below its volume with W<=28 must
have disconnected KD and nonzero H2(K+).

This remaining homology is carried even after deleting the excess-two
positive vertices. Indeed the complementary full complex then has
vertex set D together with those excess-two vertices. Section1 proves
there is no edge from the latter to ANY D vertex. Adding them cannot
merge distinct D components. Apply (7) again to obtain nonzero H2 in
the full subcomplex on positive facets of excess at least three.

The unresolved global obligation is to charge the separating boundary
facets, particularly those of excess three, strongly enough to give
f+W>=120, or to obtain another whole-volume comparison. Nonzero H2 alone
supplies no such weighted inequality. Neither treating its support as a
surface nor assuming garland structure is justified by this argument.

## Unchanged proof inputs

* `manuscripts/compact_ra4_volume_constraints.tex`, SHA256
  `afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f`:
  actual wall and incidence conventions; capacity; excess identities;
  Loebell degree-two recognition and regular-base lengths; three-edge
  closure; W=0 or W>=17 and all-D recognition; volume formula.
* `results/compact_ra4_farthest_dodecahedron_wall_radius_20261020.md`,
  SHA256 `d082e4a918d4504a59ac200f32bbb141f1f2cb52d701006eaea7ba5b0e5bcc22`:
  existence of D and the assertion that EVERY maximizing facet at any
  interior point is D.
* `results/compact_ra4_irregular_outgoing_facet_excess_count_20261020.md`,
  SHA256 `f57760ea1bf0e023857ee8798bbdfba68ae3c31a2c3dcfccd13ee6589f297795`:
  the at-least-three irregular incidences and ridge edge-carrier bound.
* `results/compact_ra4_uniform_outgoing_center_cut_volume_cost_20261020.md`,
  SHA256 `78135c4cca5d64c583727d612d5f22cf8f99f50a6a9b0aedbae17780e01c2b2c`:
  used only for Section3's low-volume uniform-center implication.

The first two global inequalities and weighted identity (3)--(6) are paper
proofs. No new executable certificate or formal-geometry coverage is
claimed. Central statuses and primary manuscripts remain unchanged.

# Compact right-angled H4: neighboring excess and zero-excess rigidity

19 October 2026. New paper proof, awaiting independent review. No scientific
calculation, classification census, optimizer, Lean, or typesetting was run.
The earlier capacity packet and its manuscript are unchanged. The universal
minimum-volume and second-volume questions are not settled here.

## Inherited inputs and normalization

Let P be a compact right-angled hyperbolic four-polytope, with f facets and
v vertices. Every facet is a compact right-angled three-polyhedron. For a
ridge R write k(R) for its number of sides, and put

    W = sum_R (k(R)-5) = v-5f,
    d(F) = sum_(R contained in F)(k(R)-5)
         = (V(F)-20)/2 = number_of_faces(F)-12.

All summands are nonnegative and

    sum_F d(F) = 2W,
    Vol(P) = pi^2(f+W+16)/12.                         (1)

Use the already reviewed capacity proof for these two statements:

* A compact right-angled three-polyhedron containing a k-gon has at least
  2k+2 faces. At equality its two nonpentagonal faces, when present, are the
  disjoint bases of the usual Loebell combinatorics.
* A facet F with V vertices in P forces f >= 2V+5.

The first statement is the actual perpendicular-edge opposite construction
and the triangulated-sphere ear argument, not an assumed finite catalogue.
The second includes the singleton-footprint closing argument.

Source pins:

* `compact_ra4_perpendicular_facet_capacity_research_20261018.md`:
  `8e6db87b001b59a16c1b9b55f03fa9373630bcfb20c1697796f4fd88fbaff8ff`.
* `review_compact_ra4_facet_capacity_root_20261018.md`:
  `44c89c904208005db2e5111cbb62bf470a54f6dbb6f58c24f91b94f77f9e64ef`.
* `compact_ra4_minimum_and_pentagon_excess_scope_audit_20261018.md`:
  `20c19e0946a8fa2eb94ff0ab345bd77553396773209f1c55b3c4733fb2de12e7`.

## 1. A neighbor forces twice its shared-ridge excess

Fix a facet F, with distinct neighboring facets F_1,...,F_M, and let the
ridge F intersect F_i have k_i sides. Applying the three-dimensional
capacity theorem inside F_i gives

    d(F_i) = number_of_faces(F_i)-12 >= 2(k_i-5).

The sum of these k_i-5 is exactly d(F). Since F and its neighbors are
distinct facets and every remaining excess is nonnegative, (1) gives

    2W >= d(F)+sum_i d(F_i) >= 3d(F).

Consequently every facet with V vertices satisfies the universal bounds

    W >= ceil(3(V-20)/4),
    f+W >= 2V+5+ceil(3(V-20)/4)
          >= 11V/4-10.                              (2)

No neighbor excess is required to be disjoint from another neighbor's
excess. It is sum_F d(F)=2W which permits this counting: a shared ridge
may occur twice on the right, just as it occurs twice in that identity.
This addresses the possible overlap rather than treating ridge charges
as separate globally available units.

The regular right-angled 120-cell has volume 34pi^2/3, corresponding to
f+W=120. Its double across one facet has volume 68pi^2/3, corresponding to
f+W=256. Since V is even, (2) implies:

| Hypothesis on P | Necessary maximum V in every facet |
|---|---:|
| Vol(P) < 34pi^2/3 | 46 |
| Vol(P) <= 34pi^2/3 | 46 |
| Vol(P) < 68pi^2/3 | 96 |

Indeed V=48 already forces f+W>=101+21=122. V=98 forces
f+W>=201+59=260. At V=46 the lower bound is117; at V=96 it is254,
so neither endpoint is removed by this argument. These improve the earlier
48/50/104 local-size cutoffs, and do not prove a universal volume gap.

## 2. A small positive-excess obstruction

Make a simple graph whose vertices are the positive-excess facets and
whose edges are the nonpentagonal ridges. Give an edge of a k-gonal ridge
the positive integral weight k-5. Its total edge weight is W, and the
weighted degree at a facet is d(F).

There is no degree-one vertex: its sole incident edge of weight w would
give d(F)=w, whereas the capacity theorem gives d(F)>=2w. If a vertex has
exactly two incident edges, of weights w_1,w_2, then

    w_1+w_2 >= 2 max(w_1,w_2),

so the weights are equal. The facet has exactly 2k+2 faces, where
k=w_1+5. The equality clause of the capacity proof makes those two
nonpentagonal faces disjoint. Its two neighboring positive-excess facets
therefore do not meet in P: if they met, their normals together with F's
normal would be orthogonal and acute projection would give an actual
common edge in F.

It follows that

    W=0 or W>=5.                                    (3)

For if 0<W<=4, the positive graph has at most four edges, minimum degree
two, and hence is a triangle or a four-cycle. A triangle contradicts the
degree-two conclusion just proved. In a four-cycle that conclusion rules
out both diagonals even in the full facet graph. Such an induced square
cannot occur: its unit-normal Gram matrix splits, after reordering, into
two 2-by-2 blocks with off-diagonal products below-1. Each block has a
negative eigenvalue, contradicting ambient Lorentz index one. The strict
gap uses compactness, exactly as in the pinned capacity proof.

This elementary obstruction is not a classification of the W=5 case and
does not by itself approach f+W>=120.

## 3. Zero excess determines the whole 120-cell

Here use one explicit external three-dimensional classification input,
already recorded in the pinned scope audit: a compact right-angled
three-polyhedron with12 faces (equivalently20 vertices) is the regular
right-angled dodecahedron. This includes its metric uniqueness; a purely
combinatorial identification would not suffice below. The known regular
right-angled 120-cell T supplies the comparison object.

**Claim.** W=0 if and only if P is isometric to T.

The reverse implication is its pentagonal ridges. If W=0, every facet of
P has d(F)=0 and is the regular dodecahedron by the stated classification.
Isometrically match one whole facet F to a facet F_T of T, with the two
four-polytopes on the same inward side of its supporting hyperplane.

Take a pentagonal ridge R of this matched facet. There is a unique other
hyperplane through the hyperbolic plane spanned by R which is orthogonal
to F's supporting hyperplane. In normals: the two-dimensional positive
normal plane to R contains F's unit normal a, and the line perpendicular
to a in that plane is unique. Thus the adjacent facets of P and T have
the same supporting hyperplane, with the same halfspace orientation.

Inside that common hyperplane, a regular dodecahedron with the fixed whole
regular pentagonal face R is unique on each side of R. This is immediate
from regularity: any two placements differ by an isometry of the regular
pentagon, and its dihedral symmetries extend to symmetries of the regular
dodecahedron. The halfspace of F fixes which side is inward. Hence the
whole adjacent facets of P and T agree, not merely their shared ridge.

Repeat across ridges. The facet adjacency graph of a convex polytope is
connected. Every facet of P is therefore a facet of T. Conversely the
nonempty set of matched T facets is closed under all their neighbors,
since every ridge in a matched whole dodecahedron is also a ridge of P.
Connectedness of T's facet graph makes this set all120 facets. The
supporting halfspaces agree, so P=T.

This argument uses the hypothesis that *every* facet is dodecahedral.
It gives no extension of a lone dodecahedral facet when other facets may
be larger, and no tessellation or garland classification of arbitrary P.

## Remaining gap and verification status

A counterexample below the 120-cell volume must now have W>=5 and every
facet at most46 vertices. A distinct equal-volume example has the same
necessary restrictions. A counterexample to the double gap has facets
at most96 vertices. These statements leave many positive-excess facet
compatibilities unresolved. No minimum theorem, second-volume theorem,
new dimension ceiling, or new example is established.

Sections1--2 are new written consequences of the internally reviewed
capacity proof and elementary convex-polytope facts. Section3 additionally
uses the explicit3D metric classification just stated. No end-to-end
formalization or independent review of this successor note is claimed.

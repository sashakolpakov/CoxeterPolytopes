# Raw23: exact exclusion of adjacent-spacelike two-ray interiors

No norm2 normal satisfying the model's native Coxeter contact conditions
over Q(√5) can have Lorentz projection in the relative
interior of any of the336 adjacent-spacelike two-ray faces of the three
original retained raw23 cones. This remains true with arbitrary passive
components. The complete finite reduction tests8,400 native anchor pairs:
4,864 have a nonpositive ray coefficient and3,536 have projected norm
greater than2. None reaches further native or arithmetic tests.

This excludes a precise class of non-extreme projected directions. It
does **not** exclude the raw23 pairs or arbitrary acute completions.
Faces containing a timelike generator, higher-dimensional projected
faces and newly formed two-dimensional faces after adding a wall remain
open. No vertex capacity, face inequality, model column, optimization or
arithmetic bound changes.

The [producer](../scripts/fresh_raw23_two_ray_completion_preflight_20261018.py),
[certificate](fresh_raw23_two_ray_completion_preflight_20261018.json) and
[standalone checker](../scripts/check_raw23_two_ray_saved_20261018.py)
bind the frozen [extreme-projection polar packet](fresh_raw23_union_polar_completion_preflight_20261018.md).
Its actual supporting-wall interpretation, complete original cones and
arithmetic ground field remain mathematical inputs.

## Complete domain and adjacency

Each original twelve-wall cone is full-dimensional and pointed in a
ten-dimensional Lorentz space. Its frozen complete extreme-ray list has
18 timelike and32 spacelike rays. There are496 unordered pairs of
spacelike rays per cone. The producer checks all1,488 pairs, without
selecting only a positive-mass source class or a visually promising pair.

Let x and y be two distinct spacelike extreme rays. All known-wall
products are nonpositive. Consequently the active walls of x+y are
exactly the common active walls of x and y. The smallest cone face
containing both has dimension10 minus the rank of those wall rows.
Thus the pair spans a two-dimensional face exactly when this rank is8.
A pointed two-dimensional face has exactly two extreme rays.

Exactly112 pairs in each cone meet this criterion. For each, the
two-dimensional Gram matrix

    H = [[x², x·y], [x·y, y²]]

is checked positive definite. This is a computed property of every
selected face, not an assumption about arbitrary pairs in a Lorentz
space. The relative interior consists precisely of z=a*x+b*y with
a>0 and b>0. Its boundary rays were treated by the frozen polar packet;
zero coefficients are not counted as new interior candidates.

## A sharp bound on each known-wall contact

Write a proposed full norm2 normal as n=z+w, where w lies in the
positive-definite orthogonal passive space. Its four passive directions
are unrestricted here. Orthogonality gives

    z² <=2,    w²=2−z².

For known scaled wall j, let s_j be its squared scale and

    p_j = (−x·wall_j, −y·wall_j).

Both entries are nonnegative. Its physical contact square with n is

    (p_j·(a,b))²/s_j.

The passive part has zero product with these twelve walls because they
lie in the known Lorentz span. This does not assert zero contacts with
the passive supporting facets; all such additional constraints are
omitted from the relaxation.

For p_j nonzero, maximize the displayed square over the closed positive
quadrant subject to (a,b)H(a,b)^T<=2. If H^−1*p_j has both entries
positive, the maximum is

    2*p_j^T*H^−1*p_j/s_j.

Otherwise the maximum is the larger boundary value

    2*p_j[0]²/(H[0,0]*s_j),
    2*p_j[1]²/(H[1,1]*s_j).

This follows from the positive-definite Cauchy inequality and the
stationary point on the compact quadrant ellipse. An interior maximum
would have to be proportional to H^−1*p_j. If that vector is not
strictly positive, the maximum occurs only on the boundary; every
relative-interior value is then strictly smaller. The checker verifies
the inverse, stationary vector and both boundary values exactly.

For every face the sharp positive-quadrant certificate identifies two walls whose
contact forms are linearly independent and whose squared contacts are
strictly below4 throughout the relative interior. A bound equal to4
is used only when equality occurs exclusively on a boundary ray.

Every chosen p_j is nonzero and has nonnegative entries. Hence a,b>0
makes both anchor contacts strictly negative, not right-angled. This
justifies five positive native values per anchor, rather than adding a
zero-contact alternative. Boundary cases remain separate.

## The finite native reduction is complete

In the fixed arithmetic field K=Q(√5), a positive finite norm2 Coxeter
contact square is one of

    1, 2, (3+√5)/2, 3, (5+√5)/2.

Indeed4cos²(π/m)=2+2cos(2π/m) lies in K only for finite labels
3,4,5,6,10: real cyclotomic degree at most2 also permits8 and12,
but their values involve√2 and√3. Label2 gives zero, already excluded
for these anchors. A strictly positive square below4 cannot be parallel
or dotted either. The fixed-field/arithmetic hypothesis is inherited;
no varying-ground-field assertion is made.

For anchor walls i,j and chosen squares d,e, the coefficients satisfy

    p_i·(a,b) = sqrt(d*s_i),
    p_j·(a,b) = sqrt(e*s_j).

The right-hand sides have positive signs because the outward products
are nonpositive. The anchor determinant is nonzero, so these equations
have exactly one real solution. There are exactly25 choices of(d,e)
per face. This is an exhaustive finite algebraic reduction of the
continuous interior, not a rational grid or a bounded denominator search.

The producer represents the two square roots in their exact
multiquadratic extension of K and tests the resulting coefficients and
norm. The saved certificate contains every choice, field basis,
coefficient vector and rejection reason. Across all336 faces:

|Exact rejection|Number|
|---|---:|
|At least one coefficient is nonpositive|4,864|
|Both coefficients are positive but z²>2|3,536|
|Remaining candidates|0|

The second rejection is independent of all unknown passive contacts:
their positive squared norm cannot reduce z². Therefore no candidate
survives to the remaining known-wall angle, algebraic-integrality,
conjugate Schur, finite-vertex retention or future-wall tests. These
further tests were not run or asserted sufficient in this packet;
they could only remove more candidates. In particular no normal
preserving the18 known finite vertices exists in the tested interiors,
because no norm2 native-anchor candidate exists there at all.

## Independent saved verification and limits

Final source-bound generation78228/040a02 passed. It compared every
case and rejection count with the initial bounded generation24238/ce8e3e
before rewriting only the unfrozen packet. Separate producer replay34278/bdddf5
and standalone checker70740/5cc055 also terminated exit0/PASS. After the
helper's docstring-only arithmetic-scope precision, final checker50141/7f4244
again terminated exit0/PASS. Metadata preflightb0a91b checked all5 direct
bindings, both Python ASTs and all4 local links.

Root completed the full producer/report/standalone-helper mathematical
read (5c4271/7373a7/75e233): PASS. It covered the complete adjacency
domain, positive-definite two-plane property, sharp quadrant bound
including boundary-only equality4, independent finite anchors and all25
choices, coefficient uniqueness/signs, passive norm relaxation and
exact8,400-candidate exclusion. It executed no scripts and did not
reprove the original cone, catalogue or field inputs. The requested
opening precision explicitly restricts the conclusion to the model's
native Coxeter contacts; arbitrary norm2 vectors need not obey those
anchors.

Construction completed the full producer/report/standalone-helper
mathematical read (946eee/f97a6b/b2dd1d): PASS. It checked the complete
1,488-pair/336-face domain, rank8 adjacency, positive planes and sharp
quadrant maxima, strictly nonzero independent anchors, all25 native
square pairs, unique coefficient recovery and passive-norm rejection.
The helper's reverse ranks, anchor equations and completed-square norms
were read with the explicitly shared field-arithmetic design. That review
executed no scripts and did not replay the original cone, field,
catalogue or actual-support inputs. Both reviews retain the relative-
interior and three-original-cone scope. Source, certificate, report and
helper are final/frozen after these completed reviews.

The standalone checker imports no producer. Its self-contained rational K
arithmetic, exact square roots and multiquadratic signs follow the frozen
independent-checker arithmetic design also used by the producer; a second
unrelated field implementation is not claimed. It reconstructs every candidate
field and independently computes ranks with reverse-column elimination,
rechecks every pair and cone bound, and verifies the saved coefficients
by the invertible anchor equations. It computes the norm using the
completed-square identity

    H00*(a+(H01/H00)*b)² + (det(H)/H00)*b²,

which differs from the producer's expanded quadratic evaluation.
The original extreme-ray enumeration, native union geometry, actual
finite-vertex argument and source catalogue remain frozen inputs; this
is not a replay of those historical computations.

The three original raw23 pair configurations remain retained unknowns.
An actual completion could use an untested projected direction, and
necessity of additional facets is not itself an exclusion. This packet
produces no extra occurrence quota and does not weaken or replace any
old ordinary, marker, signed or multi-owner model constraint.

    python3 scripts/fresh_raw23_two_ray_completion_preflight_20261018.py --verify
    python3 scripts/check_raw23_two_ray_saved_20261018.py

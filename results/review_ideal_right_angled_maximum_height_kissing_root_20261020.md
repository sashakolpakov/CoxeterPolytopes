# Root review: entirely ideal dimension obstruction

20 October 2026. Complete written mathematical review PASS, conditional
on the explicitly named published inputs. No computation or Lean
compilation is claimed in this review.

I read the complete proposed proof
`ideal_right_angled_maximum_cusp_height_kissing_candidate_20261020.md`,
SHA256 `95ef1ec3edef6083992aad06dd56e2608b1f120ec098a0d197a2320b55823537`,
and independently rederived its geometric argument. I also read the full
independent review
`review_ideal_right_angled_maximum_height_kissing_batch_guard_20261020.md`,
SHA256 `053fff76e71ea51a0be89bc2938b59e3504142a8ae697d0a303e37bf0f9e5550`.

## External statements and scope

The primary journal PDF of Kolpakov, AGT 12 (2012), printed page 1945,
**Theorem 4**, supplies both the lower bound of 24 facets for an entirely
ideal right-angled four-polytope and metric rigidity at equality. The
candidate's Theorem 2 numbering must be corrected in the current paper.
The primary Musin Annals 168 (2008) PDF, printed page 4, supplies the
cardinality bound of 24 for a spherical pi/3-code in S3. Neither theorem
is being reproved or counted as repository novelty.

Standard all-wall acuteness is a genuine geometric input: distinct
outward unit normals of an acute polytope have nonpositive Lorentz
product, including nonadjacent facets. Applying the same inequality
only at visible ridges would be insufficient. Intrinsic face inheritance
is used in dimensions at least three. Connectedness is for the actual
facet graph, with edges corresponding to ridges. These are valid in
the stated finite-volume convex setting.

## Independent calculation

Let v be a future null cusp vector and H(m)=-v.m. Every H is
nonnegative, cusp-incident walls have height zero, and some wall has
positive height. An elementary way to check the last point is that
otherwise all walls would be vertical in the cusp upper half-space
model; a nonempty full-dimensional intersection would have infinite
volume toward its finite boundary. Maximize H over the finite wall
set, obtaining H0>0 at n.

For q=v/H0+n, direct expansion gives q²=-1, q.n=0 and q.v=-H0<0.
For every distinct wall m, q.m=-H(m)/H0+n.m<=0. Thus q is an actual
point of the selected facet Q; no guessed local completion or
interior-point assumption enters. Zero-height neighbors are allowed.

Every intrinsic wall of Q is an actual ridge neighbor m_i perpendicular
to n. It is consequently still unit in n-perp, which has signature
(4,1). Its time coordinate relative to q is t_i=H(m_i)/H0 in [0,1].
The spatial unit directions satisfy

    u_i.u_j <= t_i*t_j/sqrt((1+t_i²)(1+t_j²)) <= 1/2.

Repeated directions are impossible. The space of directions is S3,
so Musin bounds the number of intrinsic walls by 24, with the correct
Euclidean dimension four. Idealness, right angles and finite volume
are inherited by Q. Kolpakov now gives equality and identifies its
metric with the regular ideal 24-cell.

In its standard unit-normal model, opposite wall normals are (1,r)
and (1,-r), with |r|²=2. Write the transported point q=(a,w).
The two previously proved t<=1 inequalities imply 2a<=2. Future
unit timelikeness gives a>=1, so a=1 and w=0. All intrinsic time
coordinates are therefore exactly one. This step uses the metric
24-cell equality theorem, not uniqueness of maximal kissing codes.

Thus every neighbor of a positive globally maximal-height facet has
the same height. Repeating the local argument at each such neighbor
is legitimate because every facet is entirely ideal. Connectivity
propagates the positive height to the original cusp's height-zero
walls, a contradiction.

Higher-dimensional entirely ideal right-angled polytopes have an
actual five-face satisfying the same hypotheses, so the exclusion
extends to all dimensions at least five. The regular ideal 24-cell
supplies dimension four. The written argument therefore proves a
sharp ceiling of four under the named published inputs.

## Limitations of this review

This is a complete ordinary mathematical proof review, with no
worldwide priority assertion. It is not external refereeing. A later
conditional Lean module may check the cardinality squeeze, path
propagation and final dimension implication; the Lorentz reduction,
face inheritance and the two external theorems must remain explicit
interfaces. A finite model script cannot replace any of those facts.

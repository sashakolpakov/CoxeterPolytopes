# Transported H22 cusp completions and the next incidence obstruction

A coherent exact construction now repairs 483,840 specified cusp links of the preceding H22 chamber by adding6,289,920 walls. The result has **11,688,030 actual Coxeter facets, but still infinite volume**. The exact strict escape has square59356314529/250000000000. No finite-volume example in a new dimension has been found.

The construction also exposes why this type of repair alone does not close the chamber: every new facet has all its intersecting neighbors through one cusp. Each therefore has an explicitly certified uncapped ideal direction and no finite vertex. Further constructive work must introduce finite intersections away from that cusp, or replace the facet; simply appending more local cusp walls with the same incidence pattern cannot produce finite volume.

## Construction

The starting chamber is the [incidence-profile extension](fresh_h22_incidence_binary_20261013.json) with 5,398,110 facets. Its original210 walls and complete previously certified symmetric banks remain. The old pair of center walls has already been replaced by four compatible incidence walls, followed by560 binary-profile walls.

A lower-placement singleton repair wall and an original norm-two wall form a parallel pair. Their sum is a null raydelta. At a representative such ray,17 original walls and that repair wall are incident. The original finite components have sizes 5,3,4,3,1. Five negative highest roots complete those components. Four explicitly rational orthogonal directions in the remaining quotient contribute their positive and negative versions, giving 13 proposed bases in all.

In coordinates whose final entry denotes the coefficient of sqrt(11), the form is rational`diag(Q,11)`. For each baseb use

`r_i = b_i + t_i*delta`, with

`t = [3,4,3,4,4,4,4,4,5,8,7,9,6]`.

The five highest-root norms are 2. The remaining four pairs have norms 5,26/5,168/13,110/7. Every wall is spacelike and negative at the fixed timelike center. Exact generator closure transports the full frame through 241,920 positive-height cusp rays; height reflection supplies their negative-height partners. Thus13 walls at each signed cusp give`13*2*241920 = 6,289,920` additions.

The [primary verifier](../scripts/fresh_h22_transported_cusp_completion_20261013_verify.py) checks the full frame orbit and all new relations by symmetry reduction, including both height signs. Every nonincident retained relation and every nonlocal new–new relation is strictly dotted. All same-cusp relations are Coxeter. Nonpositive pairings and the strict center prove actual, nonredundant facets by projecting the center onto each wall.

[Primary exact certificate](fresh_h22_transported_cusp_completion_20261013_verify.json) · [explicit local inputs](../data/enumerations/fresh_h22_transported_cusp_completion_20261013_local_proposal.json) · [frame metadata](../data/enumerations/fresh_h22_transported_cusp_completion_20261013_frames.json).

## Bounded links, with the necessary positive relation

At each specified cusp there are 31 incident normals:17 original walls, one previous repair wall, and13 new ones. Their exact Gram matrix has rank 21; their actual normal span has rank 22. More importantly, **every coefficient is positive** in

`sum(all 31 incident normals) = 66*delta`.

In the 21-dimensional positive Euclidean quotientdelta-perp/moddelta, this is a full-support positive dependence. Any recession direction has nonpositive pairing with each normal. The positive relation forces every pairing to vanish, and full quotient span then forces the direction to vanish. The link is consequently bounded. The group transport and height reflection prove this for all 483,840 specified rays.

A [separate stdlib Fraction reconstruction](../scripts/fresh_h22_transported_cusp_completion_20261013_local_independent.py) independently checks signature, all 13 norms, all 2,730 original-wall pairings, all 465 local pairings, the null relation and both ranks. Its [certificate](fresh_h22_transported_cusp_completion_20261013_local_independent.json) passed saved replay. This is not a complete cusp census of the enormous resulting chamber.

## The global volume gate fails exactly

A numerical constraint-generation scout supplied a point, which was moved slightly toward the strict center and rounded to rational coordinates. The saved rational point itself is the certificate; no numerical tolerance is used in its verification.

The [escape script](../scripts/fresh_h22_transported_cusp_completion_20261013_escape.py) replays all 5,398,110 retained inequalities and computes exact maxima across each of the 13 complete new positive-height families. A family maximum has the form`(A+B*sqrt11)/D`, with both height signs covered by taking`B>=0`. Every family satisfies`A<0` and`A^2>11*B^2`. All11,688,030 inequalities are therefore strict.

The point is future-oriented and has exact square

`59356314529/250000000000 > 0`.

This proves infinite volume. The [full escape certificate](fresh_h22_transported_cusp_completion_20261013_escape.json) records its coordinates, inherited replay, integer maxima and positive square margins. Completing many cusp links, even equivariantly and exactly, does not itself control every remaining end.

## A constructive condition for the next step

The [facet diagnostic](../scripts/fresh_h22_transported_cusp_completion_20261013_facet_obstruction.py) gives a more specific target than another large wall bank. Fix one of the new normals r, its cusp delta, and the strict center h. Define

`x0 = h - (h,r)*r/r^2`,

`eta = x0 - x0^2*delta/(2*(x0,delta))`.

The point x0 is a strict timelike point of the actual facet. The vector eta is future null, lies inr-perp, and is distinct fromdelta. All other local inequalities are strict there because their normals annihilate delta. Every nonlocal wall is dotted to r; its projected normal inr-perp is future timelike and thus strictly negative on eta. The archive checks the local identities and all 210 original-wall signs independently for each of 13 representatives.

Since`(eta,delta)<0`, the point`eta-epsilon*delta` is spacelike for positive epsilon. For sufficiently small epsilon every other inequality remains strict. A final small move toward h makes the inequalityr strict too. This supplies a geometric infinite-volume proof independent of the remote rational escape.

The same incidence condition excludes a finite completion consisting only of walls whose hyperplanes meeting r also contain delta. Parallel and dotted walls do not cut the intrinsic hyperbolic interior of the facet. For finitely many added parallel tangencies, a generic projected interior point avoids their ideal tangency points. Thus a finite-volume completion retaining r must add a wall meeting r at finite hyperbolic points and not containing delta, or must replace r. The claim assumes the retained walls remain acute-compatible facets; it is not about arbitrary inequalities that may make them redundant.

The [facet certificate](fresh_h22_transported_cusp_completion_20261013_facet_obstruction.json) contains13 exact actionable null targets. Profile 7 is marked as the preferred next target because its transported facet occurs in the numerical remaining-gap support. No new search along this direction is claimed here.

## Validation and scope

The complete primary replay, separate local Fraction replay, exact global escape replay, and structural facet replay all passed. [Review and dependency scope](fresh_h22_transported_cusp_completion_20261013_review.md) explains which arithmetic was separately reconstructed and which global arrays were replayed using the primary implementation.

```sh
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_h22_transported_cusp_completion_20261013_verify.py --verify
python3 scripts/fresh_h22_transported_cusp_completion_20261013_local_independent.py --verify
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_h22_transported_cusp_completion_20261013_escape.py --verify
python3 scripts/fresh_h22_transported_cusp_completion_20261013_facet_obstruction.py --verify
python3 scripts/fresh_h22_transported_cusp_completion_20261013_review.py --verify
```

The positive result is a coherent exact geometric construction with bounded specified links and explicit new completion targets. The negative result is infinite volume for this chamber and a scoped obstruction to preserving its cusp-only facet incidence pattern. Neither result is an unrestricted obstruction to finite-volume H22 polytopes or a literature-priority claim.

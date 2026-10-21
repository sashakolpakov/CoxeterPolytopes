# A global obstruction to perpendicular H2 completion

The current 98-wall H20 chamber cannot acquire finite volume through **any finite coherent Coxeter extension in which every added normal is orthogonal to S61**. This includes arbitrary coefficient fields, root lengths and finite labels. Even if such additions close S61 to a compact polygon, the alternative S75 face remains an infinite sector with an ambient open ideal region.

The [exact geometry certificate](fresh_dimension_global_routes_20261014_orthogonal_gate.json) and [checker](../scripts/fresh_dimension_global_routes_20261014_orthogonal_gate.py) have passed generation and saved replay. The result uses the frozen [98-wall independent review](fresh_h20_bridge_edge_closure_20261014_review.md) and the [complete S75 support-contact gate](fresh_h20_s75_projection_gate_20261014.md). No frozen predecessor was changed. This is a global obstruction to the stated extension family, not to arbitrary completions of the retained core.

## Exact geometry of the two faces

The norm-two positive-definite supports are

```text
S61 = [1,2,3,4,5,6,7,8,9,10,12,13,14,15,17,31,61,87],
S75 = [1,2,3,4,5,6,7,8,9,10,12,13,14,15,17,31,75,87].
```

They share 17 roots. Put `V = span(S75)^perp`, a Lorentzian three-space, and let `J` be the inverse S75 Gram. The exact certificate verifies

`J[75,75] = 5`,

and reconstructs `U`, the projection of old wall 61 into `V`. It verifies that `U` is strictly future timelike. Therefore

`span(S61,75) = span(S75,61)`

has dimension 19 and signature `(18,1)`. Its perpendicular space `W` is a positive-definite two-plane. Inside `V`, it is exactly the tangent plane perpendicular to `U`.

The current 98-wall chamber restricts to the **entire** hyperbolic plane in `V`: every nonsupport old wall projects strictly future timelike. This includes walls 95, 96 and 97 and is a bound input of the independent 98-wall review, not a numerical assumption.

## Every perpendicular addition is redundant or concurrent on S75

Consider a finite proposed extension retaining all 98 walls, with pairwise acute Coxeter normals and a common strict future timelike center `C`. Let an added normal `r` have square `q>0` and be orthogonal to S61.

It is right to the 17 common support roots. Its only possible nonzero S75 pairing is with wall 75. If that pairing is nonright, the Coxeter condition gives

`(r·N75)^2 >= q/2`.

The Schur projection formula then gives

`(projection_V r)^2 = q - 5(r·N75)^2 <= -3q/2 < 0`.

This projection is future timelike, not past timelike. To see the orientation without assuming an old center remains feasible, project the extension's own `C` off S75. The support inverse is entrywise nonnegative, so its projection coefficients are negative. Since every wall is acute to the support, subtracting those components makes all nonsupport pairings strictly more negative or leaves them unchanged. The projected center is future timelike and strict for every proposed wall. Hence any timelike projected normal is future. Such a wall is negative everywhere on the future H2 plane and supplies no boundary there.

If `r·N75=0`, then `r` is orthogonal to both supports and lies in the positive plane `W`. Its S75 boundary is a geodesic through the same finite point `[U]`. Thus **all effective added boundaries on S75 are concurrent at `[U]`**.

The projected strict center ensures their intersection has positive angular width. A finite intersection of their retained half-planes is a nonempty geodesic sector based at `[U]`, with no restriction on radial distance. All other projected inequalities are strict everywhere. Its area is infinite: an angular interval of positive length has area proportional to `integral_0^infinity sinh(rho) d rho`.

## An ambient strict escape, without a face-volume assumption

There is also a direct ambient argument. Choose a future null direction `nu` strictly inside that angular sector. Every nonsupport old and new inequality is strict at `nu`, while all 18 S75 support pairings vanish.

Let `Z` be minus the sum of the S75 fundamental dual vectors. Then `Z·N_i=-1` for every support root, `Z^2>0`, and `nu·Z=0`. Choose a positive number `tau` small enough that

`Z^2 + 2 tau (nu·C) > 0`;

for example `tau = Z^2/(-4 nu·C)`. For sufficiently small positive `epsilon`,

`X_epsilon = nu + epsilon Z + epsilon^2 tau C`

is future-oriented spacelike and strict for every old and new wall. Indeed its leading square is `epsilon^2 [Z^2+2 tau(nu·C)]>0`, its support inequalities are negative, and every other inequality remains negative by continuity. The family is finite, so one positive epsilon works for all of them.

Together with the strict timelike center this gives an open ideal-boundary patch and hence infinite hyperbolic volume. The proof does not require S75 itself to be a Coxeter polygon, nor does it rely on an unproved assertion about volumes of lower-dimensional faces.

## What the local doubling construction does and does not provide

The [frozen H2 doubling note](fresh_h20_h2_doubling_20261014.md) supplies a genuine finite local construction: on the 96-wall predecessor, it retains the cutter and eventually closes S61 to a compact right-angled polygon with all new normals perpendicular to S61. It does not compute that polygon, does not automatically retain the later two endpoint caps, and makes no ambient finite-volume claim. The theorem above shows that any compatible continuation of this perpendicular kind from the current 98-wall chamber still cannot finish the ambient problem.

Ordinary doubling preserves dimension and replaces a chamber by its union with a reflected copy. The even-angle doubling criterion is [Allcock, *Infinitely many hyperbolic Coxeter groups through dimension 19*, Lemma 1.3](https://web.ma.utexas.edu/users/allcock/research/doubling.pdf). A finite doubling or gluing operation that retains an isometric copy of the present infinite-volume chamber necessarily still has infinite volume. It cannot repair the gap by itself.

Nor should the H2 face be assumed Coxeter simply because the ambient chamber is Coxeter: the sufficient face theorem in [Allcock, Theorem 2.1](https://web.ma.utexas.edu/users/allcock/research/doubling.pdf) excludes spherical components of type `A_n` and `D5`. Both supports here have A-components. The local perpendicular construction avoids that issue by giving its polygon normals and angles explicitly.

## A second exact restriction and the remaining construction problem

The separate [S75 escape dual](fresh_h20_s75_completion_20261014.md), with its [independent physical identity/rank review](fresh_h20_s75_completion_20261014_gate_review.md), proves that **every** acute Coxeter cap inducing a boundary geodesic on S75 strictly misses the current rational escape `E`. It uses an exact 21-root spanning identity for `E`, not the bounded numerical profile scout.

In fact only the identity coefficients at walls 8 and 15 are negative. Every escape cutter must be nonright to at least one of those roots. Their inverse diagonals `24/11` and `4` force that cutter's S75 projection to be future timelike. Consequently an escape cutter leaves S75 completely strict. A cap that bounds S75 and a cap that cuts this escape must be different members of any successful retained-core completion.

These results supply concrete requirements for a new route:

- At least one added wall must both bound S75 and be nonorthogonal to S61; otherwise the concurrent-sector theorem applies.
- At least one further wall must cut the saved escape while leaving the whole S75 plane strict.
- These different kinds of walls must have compatible mutual angles and jointly close the other ideal regions. Neither separate local face completion nor a finite list of completed faces proves that condition.

There is an exact asymmetry worth preserving for a changed-contact construction. S75-boundary caps must be right to wall 75. A cap orthogonal to S75 with a finite contact to wall 61 is not excluded by the reverse diagonal argument: the S61 inverse diagonal at 61 is only `2/3`. This leaves a possible contact-coupling direction; no compatible globally enclosing family has been constructed here. The stronger S75 escape dual still applies to every wall that actually bounds S75.

Embedding a complete lower-dimensional source also remains only a necessary subproblem. The earlier [complete H12 face transfer](fresh_h20_contact_closure_20261014_e8_face_review.md) already established that, after all source facets are retained, no further acute spacelike normal wholly in that source subspace can be added. The [complete H13 face transfer](fresh_h20_e7_intrinsic_caps_20261014.md) likewise left an ambient strict escape. A fresh construction therefore needs transverse contacts or a changed retained core; rational similarity or a complete face alone is not an ambient completion certificate.

No theorem here excludes arbitrary finite-volume nonarithmetic H20 polytopes, or all possible extensions of these 98 walls. No new finite-volume example is claimed.

## Reproduction and scope of arithmetic

```sh
python3 scripts/fresh_dimension_global_routes_20261014_orthogonal_gate.py --verify
python3 scripts/fresh_h20_s75_projection_gate_20261014.py --verify
python3 scripts/fresh_h20_s75_completion_20261014_gate_review.py --verify
```

All three saved replays pass. The first reconstructs the common timelike point, the support inverse and the relevant signatures using the previously independent field implementation. The last independently reconstructs the dual identity and signs and proves rank by finite-field specialization. Old chamber all-pair and vertex enumerations are explicit frozen dependencies, not duplicated here.

An independent conceptual review by `/root/compact_search` found the concurrent-sector and ambient perturbation arguments sound under the stated finite-extension and strict-center premises. It specifically checked the proposed center's role in orientation, the diagonal-five dichotomy, strictness of the angular null direction, and the positive leading square of `X_epsilon`. This was a source/proof review, not another field replay.

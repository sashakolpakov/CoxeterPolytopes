# Strict-escape cutters and the coverage condition for an H8 attachment

There cannot be a global theorem excluding every coherent Coxeter cap that cuts a strict spacelike escape: **a cap dotted to every old wall always exists**. This applies to arbitrary finite acute hyperbolic chambers, including the mixed-D4 H9 seed. It does not provide finite-volume progress. The new facet is the entire H8 plane and has no finite contact with an old wall.

The [exact checker](../scripts/fresh_h9_mixed_d4_global_caps_20261014.py) and [certificate](fresh_h9_mixed_d4_global_caps_20261014.json) give a small explicit instance for the frozen [12-wall mixed-D4 seed](fresh_h9_mixed_d4_coupling_20261014.md). They retain all twelve old walls, cut the rational escape proposed in `/private/tmp/mixed_d4_robust_gap.json`, and certify all twelve new pairings as strictly dotted. The decimal escape proposal is embedded literally as rational original coordinates in the certificate; the temporary file is not a dependency.

This is a search-triage lemma and exact example, not a proposed larger finite-volume chamber. No further cap bank is generated.

## A general dotted cutter

Let the old inequalities be `(x,n_i)≤0`, with finitely many spacelike normals. Let `H` be a strict future timelike point and `E` a strict spacelike feasible point. Along

`z(t)=(1−t)H+tE`, `0≤t≤1`,

the first light-cone crossing `ν=z(t₀)` is future null and strict on every old inequality. Since `(ν,H)<0`, the identity

`(ν,E)=−(1−t₀)(ν,H)/t₀>0`

shows that the crossing already has positive pairing with the escape.

For `t>t₀` sufficiently close to `t₀`, the vector `z(t)` is spacelike, retains negative pairing with H and every old wall, and has positive pairing with E. Normalize it to

`r=√(2/z(t)²) z(t)`.

As `t` decreases to `t₀`, every `(r,n_i)` tends to negative infinity. Finiteness of the old wall set gives a common interval in which every new pair is strictly dotted. Thus r is a norm-two Coxeter cap, cuts E and shares the strict center H. Standard acute center projection proves that the new and retained facets are actual.

In the exact H9 example, rational bisection succeeds after 16 steps at

`t=65323/65536`.

The checker uses exact quartic arithmetic to verify the raw positive norm, positive escape pairing, negative center pairing, and all twelve inequalities

`(z,n_i)² > z² n_i²`, with `(z,n_i)<0`.

No expansion of the additional normalization radical is needed. All old facets remain strict at their center-projection points.

## Why this alone cannot help finite volume

The projection of an old wall to the new hyperplane is

`n_i′=n_i−(n_i,r)r/2`.

Strict dottedness makes `n_i′` timelike. Its orientation is fixed by the projected strict center, at which its pairing is negative. It therefore imposes no boundary anywhere on the future H8 hyperplane `r⊥`. The checker verifies these negative projected norms and the orientation directly for all twelve old walls.

Consequently the new facet is the entire H8 plane. Each old wall is a positive distance away; the minimum of finitely many positive distances gives an inward collar of uniform positive width. That collar has infinite volume. Any finite-volume completion retaining this facet needs at least **nine additional actual facets meeting it**, by the elementary minimum of `d+1` facets for a finite-volume hyperbolic d-polytope. This is only a necessary lower bound, not an H8 construction or a sufficiency claim.

## A precise coverage test for a perpendicular attachment

Suppose walls perpendicular to r are added so that their restriction to `r⊥` is a compact polytope `P⊂H8`. Assume the full ambient walls have the intended Coxeter pairings and a common interior point; the following is a geometric volume test, not a replacement for those angle checks. Positively normalize r and each old n to square two. Write

`g=(n,r)<−2`, `γ=−g/2>1`, `b=√(γ²−1)`,

`n′=n−gr/2=√2 b v`, where `v²=−1` is future directed.

The old retained side of r is parameterized by

`X(τ)=cosh(τ)x−sinh(τ)r/√2`, with `x∈P`, `τ≥0`.

The minus sign is essential: it makes `(X(τ),r)≤0`. The old inequality becomes

`(X(τ),n)=√2[−b cosh(τ) cosh d(x,v)+γ sinh(τ)]≤0`.

It truncates this inward ray at finite τ precisely when

`cosh d(x,v)<γ/b=|g|/√(g²−4)`.

Thus the truncating region in P is the open hyperbolic ball centered at v with radius

`ρ=arcosh(|g|/√(g²−4))=asinh(2/√(g²−4))`.

If an open subset of P misses all closed balls, its inward rays are unbounded, proving infinite volume. Boundary-only uncovered sets require separate ideal-link or integrability checks; this note makes no claim that mere closed-ball coverage suffices.

Conversely, if compact P is covered by the **interiors** of the finitely many balls, the continuous function

`f(x)=min_i [(b_i/γ_i) cosh d(x,v_i)]`

has maximum strictly below one on P. Every inward height is at most `artanh(max_P f)`. The resulting chamber is compact, subject to the previously stated compatibility and actuality hypotheses. This gives a concrete coverage objective for placing a lower-dimensional source.

## Why an all-dotted attachment fails

Let u be a norm-two side normal of P, perpendicular to r. Acuteness puts each old projected center v in its retained side. If u is strictly dotted to n, then

`(u,n)=−2b sinh d(v,u⊥)<−2`,

so `d(v,u⊥)>asinh(1/b)=ρ`.

If every sidewall of compact P is strictly dotted to every old wall, every truncating ball lies strictly inside P. Their finite union is a compact subset of the interior, leaving an open boundary collar outside all balls. That collar has unbounded inward height and infinite volume.

Parallel contacts do not remove this obstruction. If every old/side pairing is parallel or dotted, the closed balls are contained in P. Equality makes a ball tangent to a side hyperplane at one point. The relative interior of an H7 facet cannot be covered by finitely many such tangency points. Choose a facet point away from all of them; it has positive distance from the finite closed-ball union. Its neighbourhood contains an open interior patch of P outside every ball, again giving infinite volume.

Therefore attaching a compact H8 chamber with only dotted old/new interactions cannot repair this end. A productive attachment must create finite contacts to old walls so that the old truncating balls can reach the base boundary. Placement, exact contacts and coverage must be solved together; cutting one escape or closing the H8 face alone is insufficient.

## A source-independent placement condition

Acuteness of a perpendicular sidewall u with an old wall n is equivalent to `(u,v)≤0`. Every normalized projected center v must therefore lie in P. Any finite-volume attachment of this type necessarily satisfies

`conv_h{v_i} ⊂ P ⊂ ⋃ closed B(v_i,ρ_i)`.

The second inclusion follows because a point of P strictly outside the finite closed-ball union has a neighbourhood containing an uncovered open interior patch of P. This remains true when the chosen point lies on the boundary of P. Thus a single point of the hyperbolic convex hull of the centers, strictly outside all the balls, would rule out every acute perpendicular source attachment for this literal cap, independently of the source's facet count or incidence. Absence of such a witness would not prove that a suitable Coxeter source exists.

There is a convenient unnormalized test. For a raw cap z of square q, put `g_i=(z,n_i)` and `p_i=n_i−g_i z/q`. Any nonzero positive combination X of the future-timelike p_i represents a point in the center hull. It is strictly outside every closed ball exactly when

`q(X,n_i)²+g_i² X² > 0` for every i.

This avoids normalization radicals in an eventual exact certificate.

For the explicit cap above, the separate [bounded diagnostic](../scripts/fresh_h9_mixed_d4_global_caps_20261014_coverage_scout.py) and [saved output](fresh_h9_mixed_d4_global_caps_20261014_coverage_scout.json) test all 66 sums of two unit centers, all 220 sums of three, and the full barycenter. Thirteen local optimization starts are also recorded, including their failure statuses. No point strictly outside all balls was found. The coefficients are evaluated at 80 digits before numerical optimization. This is only a bounded no-find result, not a proof that the center hull is covered and not an all-source obstruction for this cap. No further search was performed.

Exact generation and saved replay pass. The scalar arithmetic reuses the previously independent quartic checker; the general interpolation and warped-coordinate arguments are proved above rather than inferred from the example.

```sh
python scripts/fresh_h9_mixed_d4_global_caps_20261014.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h9_mixed_d4_global_caps_20261014_coverage_scout.py --verify
```

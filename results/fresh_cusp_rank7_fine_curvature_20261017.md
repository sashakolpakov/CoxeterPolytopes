# Fine-budget cusp curvature from minimum factor rank7

For a product of simplex cusp links with factor ranks r_i>=7 and total
rank d=sum r_i>=29, let s count the factors, V=product(r_i+1), and let
kappa be the curvature in the [exact shape inequality](fresh_cusp_link_shape_curvature_20261017.md).
Then

    kappa/V >= 3s.

Equality occurs only for the single factor of rank29. This improves the
rank17 hypothesis of the preceding scalar comparison. It is an all-rank
inequality, not an ambient dimension theorem: the geometric fine-payment
rule must separately be proved for the chosen cusp class. In particular
it supplies the scalar side of a possible minimum-rank9 dimension29 proof.

## Monotonicity and reduction to minimal profiles

The [pointwise centered-uniform calculation](fresh_cusp_link_rank_monotonicity_20261017.md)
gives, for X_r uniform on {-r/2,...,r/2},

    sup_y [E|y+X_(r+1)|-E|y+X_r|]
      = (r+2)/(4(r+1)) for even r,
      = (r+3)/(4(r+2)) for odd r.

These identities were proved for every r by the four parity/breakpoint
formulas; they do not assume r>=17. For r=7+2t or8+2t the common cap is
(10+2t)/(4(9+2t)). The quantity5/18 minus that cap is

    t/[18(9+2t)] >= 0.

The bound is uniform in y and can be averaged over the independent sum
of the other factors. Since

    kappa/V = 3d/4-(5/2)E|sum_i X_(r_i)|,

increasing any rank already at least7 raises normalized curvature by
at least3/4-(5/2)(5/18)=1/18>0. The comparison with3s preserves s.

For s>=7, decrease every rank to7. The variance of each X_7 is21/4,
so Cauchy--Schwarz gives

    kappa/V >= 21s/4-(5/2)sqrt(21s/4) > 3s.

Both terms in the last comparison are positive. Squaring it gives
16[(9s/4)^2-(25/4)(21s/4)]=s(81s-525)>0; at s=7 the linear factor
is42 and grows thereafter. This handles arbitrarily many factors.

For 2<=s<=6, decrease ranks, never below7, until their total is
max(29,7s). This is always possible: if the total exceeds that target,
some rank exceeds7. After sorting, the remaining profiles are the21
partitions listed in the [certificate](fresh_cusp_rank7_fine_curvature_20261017.json).
For s=2,3,4 there are respectively8,10,1 profiles of total29. For s=5,6
there is one profile each, all ranks7. Their exact curvature values all
exceed3s. The smallest margin is41/256 at(7,7,7,8). These are a complete
set of minimal profiles under a proved monotonicity reduction, not a
sampled upper rank cutoff.

The coefficients of H(t)=product_i(1+...+t^(r_i)) determine each value
by kappa/V=[5 sum_j h_j min(j,d-j)-dV]/(2V). The checker computes
these coefficients twice: by polynomial convolution and independently
by inclusion-exclusion for bounded integer compositions. Two partition
recursions check that all21 minimal profiles are included.

For s=1, kappa/V=beta_n/n, n=d+1, and the existing parity identities
give beta_n>=3n for n>=30. Their numerators are4t^2+60t at n=30+2t
and4t^2+64t+36 at n=31+2t. Only n=30 gives equality.

## Why rank6 does not support this same comparison

The local Euclidean link (Delta^6)^5 has d=30 and s=5. Its exact
normalized curvature is454415/33614, so

    kappa/V-3s = -49795/33614 < 0.

Thus minimum rank7 is the smallest integer minimum-rank hypothesis
for which this particular3s comparison holds uniformly in n>=30.
The example is a valid Euclidean product-link type, and the all-pairs
fine rule on an affine A6 factor has cost3 per fiber vertex. No
finite-volume hyperbolic completion of this local type is asserted.
Different geometric weights or global incidence estimates might still
handle rank6 and smaller factors.

## Verification and geometric scope

The [checker](../scripts/fresh_cusp_rank7_fine_curvature_20261017.py)
binds the existing shape and all-rank pointwise proofs, verifies the
new exact cap and positive polynomial coefficients, and independently
computes every minimal-profile coefficient and the rank6 obstruction.
Generation, fresh saved replay, syntax and local links pass. Separate
construction-agent and frontier-agent mathematical/source reads pass the
pointwise cap, monotone reduction, all-rank variance tail, equality case
and rank6 obstruction. Those reads do not claim a duplicate arithmetic
replay. The source, certificate and this report are final and frozen. No finite dimension tail, optimization, or large
partition enumeration is used.

If all triangles and quadrilaterals receive the fine payments with
costs3n at finite vertices and3sV at ideal vertices in a class with
minimum cusp-factor rank at least7, then the shape inequality's
positive5(n-1) slack rules out every n>=30 in that class. This remains
a conditional geometric implication. The completed rank17 dimension29
theorem is unchanged by this scalar result alone.

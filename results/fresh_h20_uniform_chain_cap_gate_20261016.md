# An all-label obstruction to a uniform cap over the H20 chain

No additional wall forming Coxeter pairs with all 133 retained walls can meet all 26 added chain walls at one common finite Coxeter angle. This covers arbitrary real normal coefficients, positive normalizations, every finite label, and arbitrary right, finite Coxeter, parallel or dotted contacts to the 18 retained S61 support walls. As usual, a coherent addition must admit a common strict interior point with the retained chamber.

This is a construction-family obstruction, not a nonexistence theorem for H20 chambers. Unequal angles to the chain, walls dotted to the entire chain, and changed retained-wall sets remain outside its scope. No new wall or finite-volume chamber is constructed.

The [checker](../scripts/fresh_h20_uniform_chain_cap_gate_20261016.py) and [exact certificate](fresh_h20_uniform_chain_cap_gate_20261016.json) bind the frozen [133-wall construction](fresh_h20_s61_multicap_20261014.md) and its [independent review](fresh_h20_s61_multicap_20261014_exact_review.md). That chamber has a compact 29-sided S61 face, a second compact tracked face, at least 168 known finite vertices, and a strict spacelike escape. Its ambient finite-vertex census is not claimed complete.

## Equal chain angles force a single residual direction

Normalize the proposed wall to square 2. The 26 saved chain normals also have square 2 and are perpendicular to the S61 support. In its orthonormal Lorentz frame they are

`n_j = a x + b_j e + c_j f`, with `x²=−1`, `e²=f²=1`,

where `a≈16.766732018447345` is the frozen isolated algebraic root. The spatial pairs lie on a circle and consecutive pairs differ by a rotation with

`C=a²/(a²+2)`, `S=2√(a²+1)/(a²+2)`.

Three consecutive points are noncollinear: after removing the positive radius factor, their determinant is `2S(1−C)>0`. Equal pairings with the chain therefore force the spatial part of the proposed wall's S61-perpendicular projection to vanish.

Let `G` be the positive-definite support Gram matrix, `J=G⁻¹`, and `ω_i` its dual weights. With support contacts `d_i=−(r,N_i)≥0`, write

`r = t x − Σ d_i ω_i`, so `dᵀJd=2+t²`.

A common angle `π/m`, `m≥3`, gives `t=2cos(π/m)/a`. The exact isolating interval satisfies `a²>280`, hence

`1/a ≤ t < 2/a`, and `2 < dᵀJd < 2+4/a² < 141/70`.

For the common right angle, `t=0`. Then `r=−Σ(Jd)_i N_i`, with `Jd≥0`. A common strict interior point `h` has `(N_i,h)<0`, so `(r,h)>0` for nonzero `r`, contradicting the proposed chamber inequality. Thus the right-angle case is impossible independently of the cost gap below.

## Complete support-contact reduction

The exact support type is `A10+A4+A1+A1+A2`. Its inverse Gram is entrywise nonnegative. The two isolated A1 nodes are old walls 17 and 31, with inverse diagonal `1/2`; every other inverse diagonal is at least `2/3`. The sum of its four smallest inverse diagonals is `7/3`.

Every nonzero finite contact has value `2cos(π/k)≥1`. Thus the cost ceiling `141/70` allows at most three nonzero support contacts. High labels can also be removed without a cutoff assumption. Set `z=2cos(2π/7)>1`; it satisfies

`z³+z²−2z−1=0`.

The polynomial is increasing for `z≥1` and negative at `31/25`. Consequently every `k≥7` has `(2cos(π/k))²>81/25`. Such a contact on a non-A1 node costs more than `54/25`; on an A1 node together with any other contact it costs more than `53/25`. Both exceed `141/70`. A lone finite A1 contact has cost less than 2 and is also impossible.

It remains to enumerate one, two or three finite contacts with labels 3,4,5,6. Standalone Fraction arithmetic in `Q(√2,√3,√5)` gives:

| Support cost | Number of profiles |
|---|---:|
| less than 2 |97|
| equal to 2 |8|
| at least `137/66` |54,639|
| total |54,744|

The minimum strictly above 2 is exactly `137/66`, which is larger than `141/70`. Its eight minimizing profiles have all three labels 3, choosing one node from each of `{1,10}`, `{17,31}`, `{61,87}`. Thus no all-finite support profile can lie in the required interval.

A parallel or dotted support contact has value at least 2. On a non-A1 node it costs at least `8/3`; on an A1 node together with any other nonzero contact it costs at least `5/2`. Therefore the only remaining possibility is a single nonfinite contact at 17 or 31. The norm equation fixes it to `√(4+2t²)>2`, so it is strictly dotted.

## The two remaining continuous axes are blocked

Since both nodes are isolated A1 components, their dual weights are literally `ω17=N17/2` and `ω31=N31/2`. The remaining normal is

`r_j(t)=t x−√(4+2t²)ω_j`, for `j∈{17,31}`.

The checker reconstructs its physical pairings with retained roots from the frozen radical-tower data. Outward rational intervals at 448 bits cover the whole closed interval `[1/a_hi,2/a_lo]`, not just integer angle values:

| Single dotted support contact | Retained blocking wall | Certified pairing interval, decimal summary |
|---|---:|---:|
|17|83|`[0.88227,1.15264]`|
|31|11|`[0.95007,1.48126]`|

Both intervals are strictly positive. They contradict acuteness with the retained wall, completing the exclusion. The certificate stores exact rational endpoints; the displayed decimals are outward summaries.

## Construction consequence and provenance

The earlier [nonright S61 bank](fresh_h20_nonright_s61_20261014.md) closed 847 particular recipes and established the future-timelike residual condition. The present theorem instead removes an entire continuous uniform-chain-angle family. Any cap with finite contacts to all 26 chain walls must use at least two different angle values; allowing a nonfinite chain contact is another way outside this theorem. No claim is made that either remaining route succeeds.

The alternative fixed-source routes were checked before selecting this test. The frozen [H8 side-tilt theorem](fresh_h9_h8_side_tilts_20261014.md) and [H8 height-family theorem](fresh_h9_h8_cap_height_20261014.md) already cover those particular variations. The primary literature also limits small source replacements: the compact dimension 8 polytope with 11 facets is unique, while no compact dimension 8 Coxeter polytope has 12 facets. See [Tumarkin, *Compact hyperbolic Coxeter n-polytopes with n+3 facets*, introduction](https://arxiv.org/pdf/math/0406226) and [Felikson–Tumarkin, *On compact hyperbolic Coxeter d-polytopes with d+4 facets*, Theorem 2](https://arxiv.org/pdf/math/0510238). These facts do not exclude larger H8 sources or other construction mechanisms. [Allcock's doubling lemma](https://arxiv.org/pdf/0903.0138) applies along a wall whose finite meeting angles have even denominators; it preserves dimension. A finite union containing an infinite-volume chamber cannot repair that chamber's volume gap.

The new finite cost census uses a standalone multiquadratic Fraction implementation. The physical blocker calculation deliberately reuses the frozen field and interval helpers; this report is not an independent reimplementation of the 133-wall geometry. No external novelty or optimality claim is made.

## Reproduction

```sh
python scripts/fresh_h20_uniform_chain_cap_gate_20261016.py --verify
```

Primary generation and fresh saved replay pass. The JSON binds the checker and all exact input sources. A separate read-only mathematical and scope audit by `construction_audit` found no gap in the reductions or the two-axis conclusion. That audit did not independently replay the 54,744-profile census or the physical blocker arithmetic; a separate full arithmetic review of this new theorem has not been performed.

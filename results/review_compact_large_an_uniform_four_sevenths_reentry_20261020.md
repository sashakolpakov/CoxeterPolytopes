# Outside-author source audit: uniform actual A-rank-eight four-sevenths bound

20 October 2026. **Mathematical source-review PASS** for
`results/compact_large_an_actual_four_sevenths_uniform_transfer_candidate_20261020.md`,
SHA-256 `390bbf35828e490fe05918af404c14ba946d1cddd0deed325c93d0f36522d1a5`
(289 lines). I found no gap in the claimed actual-projection theorem:
every A_r component, r>=8, at the chosen farthest vertex has squared
projection norm >4/7, with arbitrary finite spectators. This is a
written-proof review, with no mathematical execution, Lean or TeX.
It does not establish a new ambient dimension bound by itself.

## Source bindings and independence

All three direct pins match the bytes read:

| Input under results/ | SHA-256 |
| --- | --- |
| compact_universal_a8_three_foot_actual_projection_20261020.md | `0f1b3bbe8bfcaba6f53d717664286a7b19dbaaba0911bba8626a50fe91a8800f` |
| compact_large_an_actual_half_uniform_transfer_20261020.md | `7351b471ba3b0edaf0a524ccf082010e458b7fedd8d063491a75632223bd70e7` |
| compact_universal_a7_two_end_actual_projection_20261020.md | `5ae61f571761f517222bef38f35468462ec8fb9ad9d879850c419a6c316828a5` |

The relevant inherited explicit pins also match:

| Input under results/ | SHA-256 |
| --- | --- |
| compact_universal_a6_actual_projection_three_sevenths_20261020.md | `022c0cc77512576dfc7bed74be372a94956618498c9ed74508aff53c40a50ae2` |
| compact_farthest_vertex_mixed_replacement_inequalities_20261020.md | `3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84` |
| compact_actual_wall_reciprocal_shadow_20261020.md | `c913ec3bd8e0f26a8c679a1ede2b98b251ad3f36efeae60486213ddfa3e8badc` |
| compact_farthest_vertex_orthogonal_component_opposite_obstruction_20261020.md | `4216d06c2e03fd336898b8dac289b9969fa096c8a1d08bd67ddeb8ab5ec2f23d` |

I independently derived the transfer, inventories, scalar seals and
actual-point signs below. The author is resource_cleanup; the original
A8/core-foot work and actual-wall interfaces have shared contributions
credited in their sources. My separate N17/A16 candidate also uses those
inherited feet and inventories; it corroborates the new-rank conclusion,
but is not independent evidence for those shared inputs. No existing
PASS label was used as a substitute for the mathematical checks here.
The predecessor's separate forest-composition appendix is not imported.

## Forcing a common pair, lines 48--81

For i<=j<=4, the same-side and opposite-side dual correlations add to

`(u_i+u_(N-i)).(u_j+u_(N-j))`
`=2N sqrt(i/[j(N-i)(N-j)])`.

This includes i=j. Its logarithmic derivative is strictly negative,
since both reciprocal terms exceed 1/(2N). Thus the positively weighted
paired-end Gram decreases, although each individual near-end correlation
increases. The A8 calibrations therefore transfer in the stated direction.
I checked its N9 inverse-Gram sum and rational bound
`W9^2<2614/35<747/10`, the pair demand 3269/500, and the positive
104527/1750000 squared Cauchy margin. This forces one of the four
selected common pairs. Reversing indices is only a relabeling.

## Complete profiles at cosh L=11/3, lines 83--147

I reconstructed both retained-white inventories using finite spherical
classification, including E8 at rank eight. A nonadjacent central double
contact is retained in one connected rank-r white diagram with both new
bonds internal; neither a terminal B bond nor an E/D weighted exception
is possible. Thus both labels are 3, producing the forbidden affine cycle.
For the adjacent 3/4 profile, the two displayed endpoint bounds increase
with N. At N9 they are respectively
`(6+5sqrt(2))/sqrt(6)>26/5` and
`(3sqrt(2)+4)sqrt(5/18)>41/10`, both >11/3.
These elementary bounds use sqrt(2)>7/5, sqrt(6)<5/2 and
sqrt(5/18)>1/2. Larger labels or foreign cost strengthen the endpoint
ratio. Single5@3 gives `phi^3 sqrt(5)>8` already at N9.
Consequently no central C foot is needed at this target; the only
surviving central wall is F.

The remaining single3 Q bounds follow from the two retained A-path
inverse Grams exactly as printed. In particular the full height at
positions 2,r-1 is `A^2=Q-2/N`, while Q<1/2 or Q<2/(N-3)<=1/3
permits at most one foreign component. Its cost q=Q contradicts q<A^2.
The end heights are <1/4. Single4 and single5 restrictions, the internal
affine E7/E8 trees, and zero-own-contact exclusion all persist at this
stronger target without using an old projection-energy contradiction.

For the end common pair, the only nonadjacent weighted finite-white
possibilities are 4@1+3@3 (H) and 4@1+3@r (affine F4). Adjacent
contacts are L. The single3-at2 two-foreign-component branch has both
first bonds 3, since any first bond 4 costs at least 1/2; its six-node
two-fork is the exact affine D5 principal. All remaining single-contact
costs give precisely L/H/F. In particular H and F are retained rather
than incorrectly discarded using the predecessor's 1/2 target.

Spectators are whole finite components throughout. A retained component
has at most one contact; inverse positivity gives q>=1/4 and the stated
cost gap. Replacing its contacted node has Schur remainder q-A^2:
if positive, actual acute projection makes this same wall an opposite
in an orthogonal component; if zero, the positive corrected wall is a
feasible future null vector, impossible in a compact polytope. Thus
q<A^2 is valid for arbitrary spectator ranks, with no full-star-finiteness
assumption.

## Actual feet and uniform scalar comparisons, lines 149--230

For L, positive spectator clearing uses the inverse of each entire
contacted finite block. The displayed A2 coefficients obey
`lambda2=2(1-b^2-Q)/B>0` and
`B lambda1=2bd+1-Q>0`. Direct products are zero against c,n1 and
all cleared spectator normals, and `-3h/(2B)<0` against n2.
Its negative norm follows from 0<z0<1/3. H and F have the pinned
positive coefficients and selected nonpositive products. Every other
actual wall has nonpositive product with each positive combination by
all-wall acuteness. Their unchanged positive time coefficient makes
each normalized negative vector an actual point of the full polytope.
Farthest comparison therefore concerns the same actual x_A, not merely
an edge-cover vector.

An independent raw-coordinate check gives, before centering, H entries
`(3sqrt(2)/4,3sqrt(2)/4,1/2,1/2)` and F entries
`((2-phi)/sqrt(2),...,(2-phi)/sqrt(2),(5phi-8)/sqrt(2))`,
with four copies of the first F entry. Their sums are sqrt(2)B as stated.
This checks both general-N A^2 formulas and the correlation identity (3).
For fixed foot parameters A^2 increases while k,B stay fixed. Comparing
the algebraic parameter formulas with N9 does not assume the existence
of matching polytopes at different ranks.

I rechecked all N9 table constants from their explicit feet: L's positive
quadratic has discriminant -72976/3; H has A^2/k>54/5>23409/2209;
F has A^2/k>961/81 with integer margin 5752125. Their density bounds
and the H/F far correlations follow from the raw-sum identity. At N10
the L numerator square is exactly 5/36 and
`7/15+64/2187-49/100=3891/656100>0`, giving
`rho_L^2<125/441<64/225`, with cross-product margin 99.

The L perpendicular-Cauchy expression decreases for
v>=sqrt(9/20), since rho^2 E0<v^2 J throughout. Its stated coarse
upper bound is <51/50, below 41/40 by 1/200. Both radical rational
seals and the separate H/F bounds <97/100 and <98/100 check. Thus
every far-distinct case is genuinely excluded.

## Far-common closure, lines 232--271

Normal reversal has the required minus sign. Disjoint raw supports give
the POSITIVE centered correlation `2/(N sqrt(dd'))`, which decreases
as N increases. At N9 the supports already are disjoint for every pair
in the five-row table, including H/F (4+5 coordinates). All pressure
and correlation seals transfer, and I checked every table result.
The two closest rational results, L/F and H/F, exceed 4/7 with
integer cross-product margins 14169 and 1909 respectively.

F/F is handled by actual normals, independently of any overlapping-foot
formula. With beta=phi/2 and no foreign contact,
`A^2=phi^2(N-1)/(2N)-1`, so
`c.c'=-A^2+(phi^2/4)(2/N)=1-phi^2/2+phi^2/N`.
It is >-1/2 because phi^2<3 and is negative already at N9 because
phi^2=phi+1>13/5>18/7. It decreases with N. This lies in the
forbidden actual native gap (-1/2,0); the two isolated end walls are
distinct. This seals the final alternative.

The uniform r>=8 theorem is therefore the stronger usable lemma;
the independently reviewed N17/A16 argument can remain a supporting
fallback. Classical spherical classification, all-wall acuteness,
actual compact edges and the stated actual-point interfaces remain
the geometric boundary of this written review. No new geometric
premise, arithmeticity, machine certificate or global composition
claim is introduced.

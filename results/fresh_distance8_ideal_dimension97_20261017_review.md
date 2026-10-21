# Independent review of the finite-nine / ideal-eight dimension97 proof

**PASS: independent exact scalar reconstruction and mathematical/source
composition review.** The final [primary proof](fresh_distance8_ideal_dimension97_20261017.md)
has SHA256 `7922b4dce4448026057c072540f12f47ab582d9fa055ece647c80d35e13df28c`.
No missing mathematical case was identified in its deduction of dimension
at most97 for finite-sided finite-volume real hyperbolic Coxeter polytopes.
The compact-facet ambient ceiling79 remains unchanged.

The [independent checker](../scripts/fresh_distance8_ideal_dimension97_20261017_review.py)
and [saved data](fresh_distance8_ideal_dimension97_20261017_review.json)
import no producer or other project arithmetic code. They enumerate the
reduced shapes by first choosing the exceptional smaller/larger ranks and
then forcing the multiplicity of16. This reproduces exactly957 shapes for
link dimensions97 through107, with no missing or repeated profile.

For each profile, a sliding-window integer convolution reconstructs every
coefficient of `H(t)=product_i(1+...+t^r_i)`. Independently, write

    product_i(1-t^(r_i+1)) = sum_b a_b t^b,
    s = number of ranks, D = floor(d/2), V = product_i(r_i+1).

The cumulative-tail identity gives

    G = D*V - 2*sum_b a_b binom(D-b+s,s+1).

This follows by summing `G=sum_(t=1)^D [V-2 Pr_count(S<t)]` and the
hockey-stick identity for bounded compositions. It is distinct from the
producer's truncated first-moment formula. Both calculations agree with
`G=sum_j h_j min(j,d-j)` on every profile and the boundary profile. The
checker additionally verifies the centered absolute-deviation expression
for curvature. All operations are exact integers and Fractions.

Every scalar row and per-dimension minimum agrees with the primary data.
The minimum curvature-minus-cost is `112101/1419857`, and the six-rank16
boundary at ambient97 has margin `-185352/1419857`. The finite slack
`beta_n-9n` is positive for every `n=98,...,108`. Generation and a fresh
saved replay pass all957 profiles and the boundary comparison.

The written composition was checked separately from these scalars:

* The inherited108 proof pays triangles by two paths of length at most8,
  and quadrilaterals by a distance-at-most9 corner or a component of order
  at most17. Finite corner weights remain unchanged.
* At an ideal corner the varying affine factor contains the marked pair
  and any connecting corner path. Its pair multiplicity is exactly
  `V/(r_i+1)`, including nonsimple cusps. The all-rank affine graph bound
  therefore gives normalized cost `sum_i min(8,r_i/2)`.
* The only lost old payment could be a distant ideal distance9 pair in a
  large factor. The independently reviewed
  [universal replacement](fresh_distance8_ideal_corner_replacement_20261017_review.md)
  supplies another paid corner. Neither a simple-cusp assumption nor a
  global exclusion of small affine diagrams is inserted.
* The inherited convex-order rank transfer is valid for an arbitrary cap.
  With cap16, concentrating two ranks below16 or two above16 preserves
  the stated cost and cannot increase curvature. The sum of squared ranks
  strictly increases, so reduction terminates with only copies of16 and
  at most one positive rank on either side. This proves coverage of all
  cusp partitions, not merely the957 displayed shapes.
* The unrestricted shape inequality has positive global slack. The
  scalar comparisons exclude dimensions98 through108, and the completed
  inherited108 theorem excludes the entire higher tail.

The published finite/affine classifications, actual face/cusp geometry,
all-rank convex-order argument and inherited108 proof remain explicit
mathematical inputs. This review does not rerun a high-core census or claim
formal verification. A separate frontier-agent mathematical/source read
also passes the payment, cap16 and dimension-tail composition, without
duplicating the957-profile arithmetic. External refereeing and publication
priority are not asserted.

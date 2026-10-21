# Independent review of the continuous T8 rank gate

**PASS.** The [primary certificate](fresh_h9_t8_metric_20261014.json) excludes compact H9 realizations of the saved mandatory-label assignment of T-shaped incidence case 8, for **every finite integer label `m >= 7`** at the edge `(1,2)`. It does not exclude other assignments or the whole T-shaped incidence.

The [separate checker](../scripts/fresh_h9_t8_metric_20261014_review.py) and [review certificate](fresh_h9_t8_metric_20261014_review.json) use Fraction arithmetic and an independent polynomial determinant routine. They import no primary producer or primary arithmetic. The frozen independent square-root-basis engine is reused; only its `Q(sqrt(2))` subfield is needed. All source manifests, including those of the geometry and saved assignment, are checked.

The contact domain is verified before the determinant calculation:

- The literal eight-root support `[1,2,4,5,7,8,12,13]` occurs at an actual vertex. Cross contacts between its blocks are forced right by twelve actual three-root paths with determinant `3/4-t < 0`, and twelve actual `4–3–4` paths with explicit positive null vectors.
- The mutual contact of anchors `11` and `17` is right: testing label 3 creates an actual affine triangle. Larger labels cannot restore positive definiteness.
- Four actual eight-root supports give the anchors' projection budgets. The `I2(m)` inverse has diagonal `1/(1-t) > 4` and positive entries; any nonright contact to that block alone costs more than one. Each anchor is therefore right to both parameter roots. Of the other two B2 blocks, it may contact at most one root, necessarily at label 3. This gives exactly five necessary projection patterns per anchor.
- Two actual seven-root supports contain leaf `14`, its two retained own contacts, and the foreign B2 blocks. The two known label-3 contacts cost `1/4+1/4`; any foreign nonright contact adds at least `1/2`. Thus every foreign B2 contact is right. The third projection has only the three specified nonzero pairings, at roots `2,4,5`.
- Both upper-contact witnesses for `(11,14)` and `(14,17)` are independently reconstructed as exact nonpositive leading minors at label 4. Their actual vertices do not contain both parameter roots, so these witnesses are independent of `m`. Both contacts are limited to labels 2 or 3.

Unknown acute contacts are relaxed to right in these tests. Restoring them can only increase the largest eigenvalue of the nonnegative matrix in a normalized Gram matrix `I-A`; it cannot repair a failed positive-definiteness test.

The parameter is **`t = cos(pi/m)^2`**, not `cos(pi/m)`. For every finite `m >= 7`, it lies in `(3/4,1)`. The review verifies the symbolic inverse identity for the `I2(m)` block and the full positive-block determinant `(1-t)/8`. Because both anchors are right to this block, the parameter enters only the third residual norm:

`q14 = -sqrt(2)/2 - 1/[4(1-t)]`.

All 25 anchor residual matrices have determinant `A < 0`. For each of their four possible contacts to leaf 14, the independent polynomial determinant gives

`D(t) = D0 - A/[4(1-t)]`.

The denominator is positive and `D'(t) = -A/[4(1-t)^2] > 0`. Every one of the **100** cases is already positive at `t=3/4`; the exact minimum there is **`3/4 + 5sqrt(2)/8`**. Consequently all residual determinants are strictly positive throughout the whole interval. Combined with the PD rank-eight block, this forces principal rank eleven, impossible in the ten-dimensional ambient space for H9. There are no residual rank-zero cases and no numerical label cutoff.

Generation and saved-result replay pass. The review uses no numerical acceptance predicate. The later two-parameter extension is a separate claim and is not covered by this certificate.

```sh
python3 scripts/fresh_h9_t8_metric_20261014_review.py --verify
```

# Independent review of the49th nonpolar H9 cap

**PASS:** the new norm2 wall is compatible with all48 retained walls, cuts
their common strict escape, and gives49 actual facets with a common strict
future center. No discrepancy was found in the [primary construction](fresh_nonarithmetic_common_escape_20261016.md).

The [independent source](../scripts/fresh_nonarithmetic_common_escape_20261016_cap_review.py)
and [certificate](fresh_nonarithmetic_common_escape_20261016_cap_review.json)
import only frozen independent arithmetic and geometry helpers. They solve
the ten contact equations in reverse root and column order over a separately
instantiated octic field, rather than reading the producer's inverse or
numerical proposal coordinates. Full rank10 and the exact square2 follow
from this direct reconstruction.

The defining contacts are right to0,1,2,4,5,6,27, pairing-1 to8 and10,
and pairing-phi to17. The unique reconstructed vector agrees exactly with
the primary. Independent nested-radical arithmetic verifies every remaining
old-wall contact: the cap is additionally right to polar wall40, and all37
remaining relations are strictly dotted. Thus its full contact census is
eight right angles, two pi/3 angles, one pi/5 angle and37 dotted pairs.

The inherited center stays strict. For every wall i, the reviewer checks
the exact projected point `q_i H-(H.n_i)n_i`, its negative square, zero
own-wall pairing and future orientation. Acuteness gives strict negativity
against every other wall, so all49 listed walls are actual facets.

The previous escape's pairing with the new wall is independently proved
positive, approximately0.00824633 in the saved normalization. This is an
actual escape-cutting construction, not only a compatible Gram row.

Generation and fresh saved replay pass. The preceding48-wall geometry is
bound to its [independent batch review](fresh_nonarithmetic_joint_regions_20261016_batch_review.md).
The [separate complete-cone review](fresh_nonarithmetic_common_escape_20261016_cone_review.md)
checks the finite incidence and remaining infinite-volume obstruction.
This cap certificate alone supplies no complete cone census or new-dimensional
finite-volume example.

```sh
python3 scripts/fresh_nonarithmetic_common_escape_20261016_cap_review.py --verify
```

# Independent review of the finite-angle parallel partner

**PASS.** The [review source](../scripts/fresh_h20_h4_parallel_partner_20261016_review.py)
and [certificate](fresh_h20_h4_parallel_partner_20261016_review.json)
reconstruct the [146-wall chamber](fresh_h20_h4_parallel_partner_20261016.md)
from the earlier independent physical geometry, without producer arithmetic
imports. Generation and a fresh saved replay both pass.

The new normal is defined by an exact null-line recipe. Orthogonal-frame
identities prove its squared norm is 2, its eighteen right contacts, and
its parallel contact with wall 141. The quotient defining the remaining
parameter gives the angle pi/3 with wall 140. Rational outward intervals
certify the other 125 contacts are strictly dotted, the common center is
strict, and all 146 facets are actual.

The positive relation

    r/sqrt(2) + n141/sqrt(q141) = z*nu2,  z > 0

identifies an actual ideal vertex. Its twenty incident physical normals
have rank 20 and Gram rank 19. Its link is unbounded: the affine A1 factor
still has spherical A10, A4, A2 and two A1 factors beside it. Thus this
stage alone does not complete a cusp. All possible new finite vertices
would use the nineteen finite old neighbors, but that support contains
the dotted pair 140,142; there are none.

The new wall cuts E4. All four old proper cusps and at least 171 finite
vertices persist. A different strict escape, E5, has exact square
31500061/5000000000 and proves infinite volume. The review does not give a
complete ambient cone or ideal-vertex census.

The preceding parameter isolation and old chamber acuteness remain explicit
inherited inputs. Interval overlaps check consistency of algebraic equalities;
the equalities follow from the displayed recipes and exact frame identities,
not from overlap alone. No new finite-volume example or dimension is claimed.

```sh
python3 scripts/fresh_h20_h4_parallel_partner_20261016_review.py --verify
```

# Scope and proof review of the mixed-D4 endpoint profiles

The frozen [source](../scripts/fresh_h9_mixed_d4_endpoint_profiles_20261014.py) and [certificate](fresh_h9_mixed_d4_endpoint_profiles_20261014.json) pass this narrow mathematical/source review. The result is complete for **local spherical extensions of the 52 saved open rank-eight supports**: 1,825 finite profiles and 72 isolated dihedral tails. It is not a compatible-cap census or a chamber completion.

The cycle argument forbids two attachments into one old connected component. The four-contact star and the three-contact star with a label at least4 have nonpositive principal quadratic forms. A label at least6 adjacent to any other nonright edge likewise gives a nonpositive three-node path. Thus the finite labels3/4/5 and the separate single-attachment I2(m), for every finite integer m≥6, exhaust the possibilities. These are all-label arguments; the finite table has no hidden maximum-label assumption.

Contacts in different old components have zero inverse cross terms. The stated positive Schur complement is therefore exactly the necessary and sufficient local positive-definiteness test. The source separately checks the inherited physical labels, support positivity, inverse identities, and outward frame orientation. The hyperbolic parameterization and its quadratic inequality after setting T=exp(t) have the stated signs.

For precision, t=0 must be excluded because the proposed wall is an actual additional facet, not merely because the finite start satisfies a weak inequality. At t=0 the new normal is perpendicular to the old finite vertex. Its pairings with all nine old incident simple roots are nonpositive. The inverse of their positive definite acute Gram matrix is nonnegative, so this normal is a negative linear combination of those roots. It then pairs positively with every strict interior point of the old chamber and cannot define a compatible new facet. Consequently a genuine added endpoint cap retaining the old sector has t>0; the negative residual branch meets the backward half-edge.

The [small review checker](../scripts/fresh_h9_mixed_d4_endpoint_profiles_20261014_scope_review.py) independently verifies the elementary star/path bounds, the distinct-support and count summaries, and source hashes. Generation and saved replay pass. **It does not repeat the 10,753 exact profile trials or the physical ray census.** The primary remains responsible for those arithmetic enumerations.

```sh
python scripts/fresh_h9_mixed_d4_endpoint_profiles_20261014_scope_review.py --verify
```

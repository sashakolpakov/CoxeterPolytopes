# Independent review of the complete H22 cusp repair

The uniform, closer-shift and lowered-singleton constructions pass this review. Each has 5,397,548 actual facets, completes the 241,920 specified signed cusp links, and still has an exact strict spacelike escape. The lowered version cuts the earlier escape, but a different one survives. These are separately certified chambers. There is no finite-volume discovery claim.

I freshly ran both primary checks:

```sh
python scripts/fresh_h22_new_cusp_completion_20261013_verify.py --verify
python scripts/fresh_h22_new_cusp_completion_20261013_escape.py --verify
```

Both passed. Their certificates are the [complete orbit repair](fresh_h22_new_cusp_completion_20261013_verify.json) and [strict escape](fresh_h22_new_cusp_completion_20261013_escape.json). The escape has square 366040322341/10^12>0 and satisfies all 5,397,548 inequalities strictly, uniformly over the retained center labels m≥3.

A separate [Fraction implementation](../scripts/fresh_h22_new_cusp_completion_20261013_local_review.py) reconstructs the representative cusp in the rational form diag(Q,11), representing the added coordinate as z√11. Its [certificate](fresh_h22_new_cusp_completion_20261013_local_review.json) passed creation and deterministic replay. It checks all 31 incident normals and their 465 distinct pairings, all 11 new norms, every positive affine relation, normal-span rank 22 and Gram rank 21. No producer field implementation is imported by that local checker. I also recomputed the 11 representative old-coordinate covectors and Gram contractions using unbounded Python integers; they agree exactly with the native arrays.

The main proof has the required ingredients:

- The 120,960 old-coordinate cusp vectors form the complete orbit under the archived generators. Exact closure and distinctness are checked. The connectivity test supplies a lower-index neighbor for each noninitial entry, so induction on the index proves reachability. Both tangent frames transform equivariantly, including stabilizers.
- The seven negative highest roots are rebuilt from every cusp's complete old incidence set. The connected components are chains of sizes 5,4,3,3 and three singletons. This reconstruction makes their covariance under the same group explicit; it is not assumed from the proposal arrays.
- For each cusp, old component roots plus their added negative highest root sum to 10δ. The old parallel pair sums to δ, and each of the two added residual pairs sums to 20δ. Their mutually orthogonal positive quotient ranks add to 18+1+2=21. A recession vector must pair nonpositively with every link normal; each positive relation forces those pairings to vanish, and full quotient rank forces the vector to vanish. Thus these links are bounded, rather than merely having a large Gram rank.
- All 11 representative types are checked against the entire positive new orbit. Equivariance transfers this to every same-sign new pair. Opposite height signs decrease pairings; the uniform shift and the norm bounds make them strictly dotted. The bounds against every retained family use rational upper bounds for square roots, including the full range m≥3 for the old center caps.
- The common timelike center is strict. Once complete pair acuteness is known, its orthogonal projection onto any wall remains timelike and strict on every other wall. Hence the listed distinct planes are actual facets; the facet count is not merely a count of inequalities.

The review retains the dependencies on the prior exact Lorentz signature, original chamber and earlier complete-orbit certificates. Their hashes are checked by the primary replay. The proposal NPZ arrays are not accepted as geometry without the exact closure, incidence and pairing tests.

## Closer-shift placement

Fresh replays of the cleaned sources also passed:

```sh
python scripts/fresh_h22_new_cusp_completion_20261013_tight.py --verify
python scripts/fresh_h22_new_cusp_completion_20261013_tight_escape.py --verify
```

The [closer-placement certificate](fresh_h22_new_cusp_completion_20261013_tight.json) uses shifts

```
[3, 3, 4, 4, 4, 4, 4, 10, 5, 1/20, 1/20].
```

The highest-root shifts depend only on component size, so the archived old symmetry group preserves them. The oriented tangent frames were already proved equivariant, permitting different shifts for their two signs. Shifting a base normal by a multiple of its orthogonal null vector leaves its norm unchanged. The primary check recomputes the full positive-height pair bounds for this placement; it does not reuse the uniform-shift margins.

There is a separate reason that opposite-height pairs are dotted even with these unequal shifts. Every new normal has a future timelike old-coordinate projection, checked by its negative square and negative pairing with the old center. Thus two such old projections have negative pairing. Their opposite added coordinates make the full pairing still smaller, with magnitude greater than the product of the added-coordinate magnitudes; each of those exceeds the square root of the respective full positive norm.

I checked the constants in all three rational cross-family estimates: the timelike family uses √(55/2)/10<21/40; the normalized null family uses √11/10<1/3; and the center family uses 4√11/10<4/3 together with center norm at most 2. The new k=6 family comparison is exact. The squared inequalities always also check the pairing's negative sign.

The separate Fraction checker now accepts `--tight` and produces a [second local certificate](fresh_h22_new_cusp_completion_20261013_tight_local_review.json), while preserving the previous uniform certificate unchanged. It reconstructs all 31 incident normals, 465 pairs, 11 norms and the two ranks again. The seven highest-root affine multipliers are their positive shifts; the two residual-pair multipliers are now 15 and 1/10. The resulting positive relations still give zero recession and bounded links.

The [closer-placement escape certificate](fresh_h22_new_cusp_completion_20261013_tight_escape.json) reuses the same rational spacelike point and verifies every new inequality separately. For a group of walls with twenty-times-shift t, the maximum old-coordinate pairing m is negative and satisfies 100m²>11t²DEN²z². This is precisely the strict inequality against both signs of the added coordinate. The square remains 366040322341/10^12>0, so the closer placement also has infinite volume.

## Lowered singleton placement

The final cleaned source and new escape checker both passed fresh replay:

```sh
python scripts/fresh_h22_new_cusp_completion_20261013_tight.py --lower --verify
python scripts/fresh_h22_new_cusp_completion_20261013_lower_escape.py --verify
```

This [third placement](fresh_h22_new_cusp_completion_20261013_lower.json) uses shifts

```
[3, 3, 4, 4, 1, 1, 2, 10, 5, 1/20, 1/20].
```

Its new rule depends on component size and, for singleton components, norm 2 or 4; these invariants are preserved by the archived symmetry group. The full new/new orbit gates and all retained-family bounds were recomputed. The separate Fraction checker also passed with `--lower`, producing a [third local certificate](fresh_h22_new_cusp_completion_20261013_lower_local_review.json). It now independently checks all 2,310 relations between the 11 new representative types and the 210 original walls for all three placements; the older saved local certificates remain unchanged.

Lowering the norm 2 singleton shifts from 4 to 1 changes each representative's original-wall relations to 18 right, one π/3, two parallel and 189 dotted. Lowering the norm 4 singleton shift from 4 to 2 gives 17 right, three π/4, one parallel and 189 dotted. These are exact Fraction checks. In particular, new parallel relations occur: the proof completes the specified 241,920 cusp links, **not a complete census of every cusp of the lowered chamber**.

The [new global escape certificate](fresh_h22_new_cusp_completion_20261013_lower_escape.json) proves both the constructive cut and the remaining failure of finite volume. Highest-root orbit entry 613840, at cusp 87691, has strictly positive pairing with the preceding remote escape. A different rational point has square 3511216341/10^10>0 and is strict on all 5,397,548 inequalities for every old center label m≥3.

A further [separate Fraction check](../scripts/fresh_h22_new_cusp_completion_20261013_lower_escape_review.py), with its [saved certificate](fresh_h22_new_cusp_completion_20261013_lower_escape_review.json), reconstructs that cutting wall directly from its archived highest root and cusp frame. It verifies norm 2, the exact positive pairing A+B√11, and the new escape's square, future orientation and 210 original strict inequalities. The full remaining inequality coverage comes from the fresh primary global replay.

Scope: all three placements above have been source-reviewed and freshly replayed. The separate arithmetic does not independently repeat the millions of distant pair checks. No claim is made here about combining the placements or about all real choices of shifts.

Local arithmetic replay:

```sh
python scripts/fresh_h22_new_cusp_completion_20261013_local_review.py --verify
python scripts/fresh_h22_new_cusp_completion_20261013_local_review.py --tight --verify
python scripts/fresh_h22_new_cusp_completion_20261013_local_review.py --lower --verify
python scripts/fresh_h22_new_cusp_completion_20261013_lower_escape_review.py --verify
```

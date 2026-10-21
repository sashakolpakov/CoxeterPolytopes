# Independent review of the cubic H9 fixed-point symmetry

**PASS.** The nearby fixed-point argument correctly certifies the proposed involution as a symmetry of the full reflection chamber of `⟨−2θ⟩+E8+⟨1⟩`, conditional on the inherited complete height-100 Vinberg prefix. It does not require first-contact certification at the far transported corner.

The [separate checker](../scripts/fresh_h9_cubic_corner_20261014_review.py) uses the project's `Cubic7` arithmetic, rather than the producer's independent Fraction-polynomial engine, and imports no producer. Its [saved result](fresh_h9_cubic_corner_20261014_review.json) binds the [extension certificate](fresh_h9_cubic_corner_20261014_extension.json), its sources, the original prefix and the field implementation. Generation and saved `--verify` both pass.

The replay independently verifies:

- `T²=I`, exact form preservation and integrality of all ten actual lattice-basis images. Since T is its own inverse, both module directions are integral. Its trace is two, so this is a finite-order symmetry.
- Future orientation, the literal fixed point `H=u+Tu`, its negative norm, all seventeen old-wall inequalities, and its seven incident simple facets.
- The exact positive normal permutation `[0,1,2,3,4,7,8] -> [4,3,2,1,0,7,8]` at H. There is no sign reversal or numerical normal matching.
- The distinct union of the seventeen original walls and their images: twenty-two normals, with all 231 pair gates independently recomputed (128 right, 19 finite nonright, 84 dotted), and a strict timelike center.
- The transported initial `E8+A1` corner, its nine incident walls and all twenty-two inequalities. The former prefix escape is cut, as reported.

## Why the bounded prefix suffices at H

Write `a=2θ` and `u=(1,0,…,0)`, so `u²=−a`. For a future timelike point X, its squared hyperbolic radial sine is

```text
sinh² d(u,X) = a X₀²/(−X²) − 1.
```

For a spacelike reflecting normal r, the corresponding squared sine of the distance from u to its mirror is `a r₀²/r²`, exactly the search's Vinberg height. Every mirror separating u from H meets the connecting segment, so its distance from u is at most `d(u,H)`. The same bound holds for every mirror containing H. Height-zero mirrors are handled by the initial `E8+A1` sector.

The independent calculation gives

```text
sinh² d(u,H) = 2+7θ+3θ² < 16 < 100.
```

Consequently no omitted full-chamber wall can separate H or be incident to it if the inherited height-100 prefix is complete. The exact prefix inequalities place H in the full chamber, and its incident simple facets are precisely the seven recorded ones. This uses the actual distance formula, not an informal comparison of coordinate heights.

Completeness of the norm alphabet, height-zero chamber and bounded Vinberg enumeration is an inherited premise, with its [previous scope audit](fresh_cubic7_h9_independent_scope_audit_20261013.md). The present review checked how that premise is used and binds the search sources; it does **not** repeat the positive-height candidate census. In particular, the earlier saved-root audit alone is not being relabeled as an independent exhaustion proof.

## From a fixed sector to the full chamber

An integral lattice isometry with integral inverse conjugates every integral reflection to another such reflection. T therefore preserves the full mirror arrangement. It is future-preserving, fixes H, and permutes all simple facets at H with their inward half-spaces. Thus it preserves the local chamber sector at H. Two reflection chambers with the same sector there have the same interior germ and are the same chamber. Hence `T(C)=C`.

This reasoning applies at a face point; H need not be a full finite corner or a strict interior point. It also avoids any assumption that the distant proposed corner had already been reached by a complete edge walk. Once T is known to stabilize C, the image of the genuine initial finite corner and every image of an actual chamber wall are genuine.

The extension establishes a full-chamber **involution**, not an infinite stabilizer. The twenty-two-wall intersection still requires a global volume analysis, and loss of the old escape is not a compactness certificate. The independent review makes no reflectivity or nonreflectivity conclusion and does not audit the separate complete elliptic-edge census.

```sh
python scripts/fresh_h9_cubic_corner_20261014_review.py --verify
```

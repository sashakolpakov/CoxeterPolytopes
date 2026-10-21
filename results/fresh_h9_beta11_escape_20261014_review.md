# Independent prefix escapes and finite-vertex persistence

**PASS.** Independent arithmetic confirms strict future spacelike escapes
for the 21-wall height-100 and 35-wall height-300 Vinberg prefixes. Both
finite prefixes therefore have infinite volume. This does not imply
nonreflectivity or infinite volume of the eventual full chamber.

The [checker](../scripts/fresh_h9_beta11_escape_20261014_review.py) and
[certificate](fresh_h9_beta11_escape_20261014_review.json) use standalone
Fraction pairs in \(1,\phi\), with no producer or project arithmetic
imports. They bind both Vinberg inputs and both geometry certificates:
[height 100](fresh_h9_beta11_geometry_20261014.json),
[height 300](fresh_h9_beta11_geometry_h300_20261014.json).

From the literal form, the review independently checks the E8 leading
principal minors and binary block, obtaining signatures \((9,1)\) and
\((10,0)\) at the two embeddings. It reconstructs all root norms, both
geometry centers, both saved Vinberg centers, and the two escape vectors.
All 112 center-wall products are strictly negative; all 56 corresponding
escape-wall products are strictly negative. Each center is future
timelike, each escape is spacelike, and each escape has negative product
with its center. The certificate saves the exact products and norms.

The segment from a strict timelike center to its strict spacelike escape
crosses the future light cone while remaining strictly inside every wall
halfspace. It therefore produces an open ideal-boundary region, proving
infinite volume. This argument needs no complete extreme-ray census.

The old height-100 escape is also checked against all 35 height-300 roots:
**only root 25 has positive pairing; none has zero pairing**. The new
height-300 witness survives all 35. This explicitly separates a prefix
escape from a permanent obstruction.

## Finite vertices do persist under later simple walls

The proposed preservation argument is correct for a true acute Vinberg
prefix. Let \(v\) be a finite vertex, and let \(n_1,\ldots,n_9\) be
independent incident outward roots. Their Gram matrix B is positive
definite with nonpositive off-diagonal entries, hence \(B^{-1}\) is
entrywise nonnegative. Let H be the future timelike Vinberg basepoint,
already on the retained side of each prefix wall. Decompose

\[
 H=a v+\sum_i b_i n_i,\qquad
 a=\frac{(H,v)}{v^2}>0,\qquad
 b=B^{-1}\bigl((H,n_i)\bigr)_i\le0.
\]

For a later actual positive-height simple root r,
\((H,r)<0\) and \((n_i,r)\le0\). Thus

\[
 0>(H,r)=a(v,r)+\sum_i b_i(n_i,r)\ge a(v,r),
\]

so \((v,r)<0\). Every such later wall strictly retains the vertex and its
incident walls. The complete zero-height sector must already have been
chosen, as it has here; an unprocessed zero-height wall would not satisfy
the strict premise. Local finiteness of the integral mirror arrangement
then gives the same finite vertex and local chamber sector in the full
chamber. No distance-to-basepoint cutoff is required.

This proof applies to the finite vertices certified by the primary cone
calculation. The present arithmetic review does **not** repeat its full
ray census or verify each of the reported 880 height-300 vertices; those
data remain an inherited prerequisite for that numerical count.

There is a useful precise map corollary. Suppose a future-preserving
integral isometry maps one such finite vertex to another and maps the
**complete incident outward root rays with positive scales** onto the
other incident rays. It preserves the full integral mirror arrangement
and maps the first local chamber sector to the second. These are sectors
of the same full chamber, so the isometry stabilizes that chamber. Equal
abstract incident diagrams alone do not establish integrality, future
orientation, root-ray orientation or this conclusion. An infinite-order
test is still required before using such a map as a nonreflectivity
certificate.

Generation and saved replay, source manifests, syntax and links pass.
No frozen predecessor was edited.

```sh
python scripts/fresh_h9_beta11_escape_20261014_review.py --verify
```

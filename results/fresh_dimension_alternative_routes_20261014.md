# A coupled H3-face lift: an exact first obstruction

The bounded new experiment was to use the known compact six-dimensional Bugaenko chamber with ten facets as an `H3`-supported face of a compact H9 chamber. This differs from the previously tested H5+H4 seed: H3 admits a spherical H4 extension, apparently allowing a nonorthogonal facet to meet the source face. The exact test below removes that apparent advantage **for this particular intrinsic source**. Every source-face facet is forced orthogonal to H3. No compact H9 example or general nonexistence theorem is claimed.

The [standalone exact script](../scripts/fresh_dimension_alternative_routes_20261014.py) and [certificate](fresh_dimension_alternative_routes_20261014.json) check the H3 inverse, all inequalities used in the all-label argument, and a nonright-neighbor witness for each of the ten source facets. Generation and saved replay pass. The original H6 source's compactness is the published input, not a new full geometric census.

## Why this was a distinct candidate

The old H5+H4 experiment starts with a complete compact H5 face and a positive H4 transverse block. Every nonright H4 contact consumes more than the entire positive norm budget. H3 instead has a terminal inverse diagonal below two. The known H6 source is explicitly exhibited in Burcroff's Figure 5 and §8, attributed there to Bugaenko. [Burcroff, *Near classification of compact hyperbolic Coxeter polytopes with few facets*](https://arxiv.org/pdf/2201.03437).

This experiment reverses a legitimate Coxeter-face operation: a spherical H3 support satisfies the Borcherds face theorem's restriction on spherical components. It does not assume that arbitrary facets of Coxeter polytopes are themselves Coxeter. [Allcock, Theorem 2.1](https://arxiv.org/pdf/0903.0138).

The August 2026 H6 preprint's v1 makes the unrestricted facet assertion in its reductions. Its current official code repository describes a larger search that avoids those reductions. Neither the v1 uniqueness claim nor the expanded search is needed here; only the previously published explicit H6 example and the literal diagram are used. [Preprint v1, §4](https://arxiv.org/html/2608.02894v1), [official implementation and diagram](https://github.com/johnamcleod/10-facets-6-d-polytopes).

## The all-label local gate

Normalize the H3 roots to square two and write `phi=(1+sqrt(5))/2`. Their Gram matrix and inverse diagonal are

```text
B = [[2,-phi,0],[-phi,2,-1],[0,-1,2]],
diag(B^-1) = ((9+3sqrt(5))/4, 3+sqrt(5), (5+sqrt(5))/4).
```

The entire inverse is strictly positive. A norm-two ambient wall inducing a facet of the H3-supported face has positive projected square `2-c^T B^-1 c`, where its support pairings are `-c_i`. Each nonzero Coxeter contact has `c_i>=1`. The positive inverse eliminates two or more contacts, and either of the first two single contacts. At the last node, label four or larger has `c_i^2>=2` and is also impossible. Thus an inducing wall is either orthogonal to H3 or has exactly one terminal label-three contact. In the latter case its projected square is

```text
kappa = (3-sqrt(5))/4.
```

Call this latter kind special. For a special wall and an orthogonal wall, an ambient nonright contact gives intrinsic squared cosine at least `1/(2*kappa)>1`. They therefore cannot meet in the source face. If they do meet, their ambient and intrinsic angles are right. Two special walls have projected pairing at most `-B^-1[2,2]`; its magnitude divided by `kappa` exceeds one, so they cannot meet each other either. Consequently **every special source facet must meet all of its intrinsic neighboring facets at right angles**.

This conclusion allows arbitrary ambient finite labels, arbitrary divergent distances, and arbitrary choices of which source facets might be special. It is a local positive-Gram argument, not a finite label scan.

## Application to the explicit H6 source

In one-based numbering, the source has the following finite nonright edges:

```text
(1,2):5, (2,3):3, (3,4):3, (4,5):3, (5,6):3,
(9,10):5, (2,7):4, (2,8):4, (5,9):3, (6,10):5.
```

Every one of its ten facets occurs in this list. Hence no facet can be special; all ten must be orthogonal to H3. The list is also given in the author's machine-readable reconstruction of Burcroff's diagram. [Explicit diagram source](https://raw.githubusercontent.com/johnamcleod/10-facets-6-d-polytopes/main/verify_polytope.py).

If these ten orthogonal walls are retained, acuteness forces every additional old-space projection into the complete compact source cone, hence into its future causal cone. In particular the seemingly favorable positive residual `kappa` cannot produce a new compatible wall. This explains why simply choosing H3 instead of H4 does not open the desired contact mechanism.

## What remains executable

A coupled H3-face construction must change the intrinsic H6 metric/source so that at least one facet is right-angled to every neighbor, or change the spherical support. A source with such a facet can be tested by replacing its normal with the unique H3-to-H4 profile above, then checking all other pairings and the full projective ray cone. The current Bugaenko source fails the necessary test before any cap bank is generated. Outer walls disjoint from the H3 face remain possible, but they return to the already difficult global covering problem; this note does not exclude them.

The literature operations examined here do not supply another automatic dimension increase: doubling and gluing keep the ambient dimension. The earlier [H18 two-coordinate experiment](fresh_dimension_construction_routes_20261014.md) records their exact scope, and the new [arbitrary-height obstruction](fresh_h18_general_height_lift_20261014_review.md) now closes its entire bare fixed-source-direction lift family. A fresh successful construction must change projected directions or the source incidence/metric, not only choose more general finite height sets.

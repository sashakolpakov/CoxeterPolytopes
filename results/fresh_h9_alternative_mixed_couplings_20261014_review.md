# Independent rank-four coupling review

**PASS:** the [separate checker](../scripts/fresh_h9_alternative_mixed_couplings_20261014_review.py) reproduces the [fixed-source classification](fresh_h9_alternative_mixed_couplings_20261014.md), with a [saved exact review](fresh_h9_alternative_mixed_couplings_20261014_review.json). It imports no coupling producer. Its arithmetic uses rational coefficients in the basis `1, √2, √5, √10`, permutation determinants, and adjugate inverses.

All 4,096 labeled support diagrams are tested: 182 are positive definite and reduce to 18 types. Their 216 positive-projection contact profiles and all 5,168 pairs at the source π/5 edge agree exactly. Only the ordinary/ordinary pair survives that edge. The original source's finite-edge graph is reconstructed from its exact normals: the π/3 chain forces facets 0–5 ordinary, and the remaining π/4 attachment allows only zero or unit support-weight square. Only irreducible B4 and D4 support the latter. Every reducible allowed profile leaves an entire support component unused.

The unbounded-label reductions are analytic. A support edge with label at least six must form an isolated dihedral factor; its inverse diagonal is at least two and prevents a nonzero contact with positive remaining projection square. A contact label at least six can occur only on an isolated A1, with no second support contact. The acute Coxeter comparison with any finite source neighbor then excludes it. Thus the finite census is not an assumed label cutoff.

This is a review of the fixed seven-facet H5 source, four spherical support facets, and exactly one ambient lifted wall per source facet. It neither classifies arbitrary H9 chambers nor excludes additional caps on the resulting seeds.

```sh
python scripts/fresh_h9_alternative_mixed_couplings_20261014_review.py --verify
```

# Complete geometry of the 30-wall chamber

The [nine-contact wall](fresh_h9_nine_contact_cap_20261015.md) extends the separate geometric chamber to 30 actual facets. A complete exact cone census gives:

| Quantity | Count |
|---|---:|
| Finite vertices | 89 |
| Ideal vertices | 0 |
| Hyperideal extreme rays | 805 |
| Finite edges | 250 |
| Open finite-corner germs | 301 |

All 88 old finite vertices persist, and one new vertex appears. The identity \(9\cdot89=2\cdot250+301\) checks the finite edge incidences. All 435 wall pairs satisfy the Coxeter conditions.

The [primary checker](../scripts/fresh_h9_nine_contact_extension_20261015.py) produces the [full certificate](fresh_h9_nine_contact_extension_20261015.json), including wall normals, the Gram matrix, actual-facet witnesses, the complete ray census, edge incidence and an exact strict escape. The escape is spacelike, lies on the future side and pairs strictly negatively with every wall. Hence the chamber has **infinite volume**.

The [independent checker](../scripts/fresh_h9_nine_contact_extension_20261015_review.py) and [review certificate](fresh_h9_nine_contact_extension_20261015_review.json) reconstruct the entire cone with separate field arithmetic and complete-face clipping. They reproduce every ray support, vertex, finite edge, open germ and escape sign. Primary and independent generation and saved replay pass.

A [bounded numerical rescore](fresh_h9_nine_contact_rescore_20261015.md) finds no cutter of this new escape among compatible entries in three saved proposal banks. New frontier geometry or changed retained walls remain outside that diagnostic. No new finite-volume example or dimension is established.

```sh
python scripts/fresh_h9_nine_contact_extension_20261015.py --verify
python scripts/fresh_h9_nine_contact_extension_20261015_review.py --verify
```

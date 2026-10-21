# The fixed H8 source fails along the entire dotted-cap segment

**No finite-volume attachment of the specified Bugaenko H8 source works anywhere on the all-dotted part of the frozen H–E cap segment.** This covers every isometric placement of that source perpendicular to the cap, with all twelve old walls retained. It does not exclude other cap directions, other source chambers, additional source facets, nonperpendicular attachments or changed retained walls.

The [exact checker](../scripts/fresh_h9_h8_facet_cover_20261014_line_exclusion.py) and [certificate](fresh_h9_h8_facet_cover_20261014_line_exclusion.json) combine the [diameter theorem](fresh_h9_h8_facet_cover_20261014_parameter_gate.md) with the [source contact-gap theorem](fresh_h9_h8_facet_cover_20261014.md). No numerical placement search is a premise.

Set `a=(√5−1)/4`. For a norm-two old wall, its projection p to the cap hyperplane has `κ²=−p²`. The gap theorem says that if `κ²<4a`, every allowable projection is a positive multiple of one of three exceptional source-vertex rays, at source vertices 0, 32 and 33. Each ray's minimal contact vector has norm square a and has a nonzero contact equal to one. Its multiple by s therefore has a contact equal to s, which must belong to the Coxeter set `0,1,√2,φ,…` or be at least two. In particular, no allowable projection has `a<κ²<2a`.

The diameter theorem restricts a possible finite-volume attachment to

`0.99674949 < t < 0.996752404`.

Split this interval at the exact rational value `T=0.9967505`. The same theorem proves that every κ² decreases as t increases.

* For `t≤T`, wall 7 has `a<κ_7²<2a`. The lower inequality is checked exactly at T, and the upper at the interval's lower endpoint. Its exceptional-ray multiplier would lie strictly between one and √2, which is impossible for a Coxeter contact.
* For `t≥T`, walls 0, 1, 2 and 3 all have `κ_i²<4a`, checked exactly at T. Their four centers are distinct for every all-dotted cap: the original mutual products are nonpositive, and the projected center-distance formula gives cosh distance strictly larger than one. Four distinct centers cannot occupy only three exceptional source vertices.

Both parts include T. Together with the exact outer diameter exclusions they cover the entire all-dotted segment. This is a stopping rule for one concrete attachment construction, not a nonexistence theorem for compact H9 polytopes.

```sh
python scripts/fresh_h9_h8_facet_cover_20261014_line_exclusion.py --verify
```

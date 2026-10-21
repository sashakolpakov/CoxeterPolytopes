# Compact arithmetic dimension21: one field remains

**Internally reviewed research reduction:** if a compact arithmetic Coxeter
polytope of dimension21 exists, its ground field must be **Q(√5)**.
This does not exclude dimension21 itself or lower the overall arithmetic
bound21 candidate. External validation and priority remain open.

The [binder](../scripts/fresh_compact21_single_field_20261016.py) and
[certificate](fresh_compact21_single_field_20261016.json) combine the frozen
[two-field reduction](fresh_compact21_field_reduction_20261016.md) with the
new [F25 proof](fresh_sqrt2sqrt5_lanner_catalogue_20261016.md) and its
[independent review](fresh_sqrt2sqrt5_lanner_catalogue_20261016_review.md).
Generation and a fresh saved replay pass. All recursive source manifests
are checked; this composition does not rerun inherited large enumerations.

The thirteen dimension21 branches are retained explicitly: the rational
field is excluded by the published isotropy argument; eleven nonrational
fields now have reviewed bound20 candidates; Q(√5) remains open. The larger
dimension14 field list is used here, never the three-field restriction
valid only from dimension22. The new exclusion concerns exact field
Q(√2,√5), and is not applied to its proper subfield Q(√5).

The arithmetic bound21 candidate, the published noncompact arithmetic
bound21, the unrestricted compact bound29, and the unrestricted finite-volume
163 candidate are unchanged. No existence result follows from the remaining
field branch.

```sh
python3 scripts/fresh_compact21_single_field_20261016.py
python3 scripts/fresh_compact21_single_field_20261016.py --verify
```

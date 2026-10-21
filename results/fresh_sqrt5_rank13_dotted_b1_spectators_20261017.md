# 119,698 rank13 quadrilateral contexts from the dotted cores

Every one of the1,399 certified dotted B1 cores has been extended by
**every finite orthogonal spectator forest of the required rank**. This
produces119,698 locally realizable compact quadrilateral contexts in H13.
After identifying equal marked corner types and duplicate face sums,
the data contain:

| Quantity | Count |
|---|---:|
| Local quadrilateral cores | 1,399 |
| Spectator contexts | 119,698 |
| Distinct four-corner weight inequalities | 82,263 |
| Marked corner types | 63,618 |
| Whole rank13 vertex forests appearing | 6,716 |

The [source](../scripts/fresh_sqrt5_rank13_dotted_b1_spectators_20261017.py)
and [certificate](fresh_sqrt5_rank13_dotted_b1_spectators_20261017.json)
use the [exact core package](fresh_sqrt5_rank13_dotted_b1_quad_cores_20261017.md).
Spectator ranks0 through11 have respectively
1,1,6,9,29,47,117,196,419,704,1366,2276 finite component multisets.
All finite Q(sqrt5)-compatible component types are included.

Adding the spectators as orthogonal white roots leaves the common-white
face cone unchanged. The admissible form, stable integral reflection
module, and strict interior point extend by orthogonal positive blocks,
as proved in the [realization lemma](fresh_singular_quad_local_realization_20261017.md).
Each resulting Gram has15 roots,11 white normals, physical inertia(13,1,1)
and conjugate inertia(14,0,1). Its four corner diagrams have rank13.

The JSON interns marked corners in `corner_types`. Each `face_rows` entry
has required payment1 and a list of `[corner_type_index,multiplicity]`
pairs; its `contexts` list retains every originating core and spectator
choice. No weight vector is assumed. This format can be joined directly
to the whole-forest corner-orbit model by `forest` and `marked_forest_code`.

These data bring the local constraints into6,716 of the6,789 finite
rank13 forest types. That is coverage of vertex types appearing in this
family, not completeness of all triangles, quadrilaterals, or possible
ambient chambers. A new dimension theorem still requires a complete
face-payment argument and the appropriate global vertex budget.

Final generation and fresh saved replay passed after normalizing the
core's JSON representation. All119,698 contexts and82,263 distinct rows
replay exactly. Syntax, source manifests and local report links are
checked separately. The package is frozen.

```sh
python3 scripts/fresh_sqrt5_rank13_dotted_b1_spectators_20261017.py --verify
```

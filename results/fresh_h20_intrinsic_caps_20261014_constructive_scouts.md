# Intrinsic H20 cap construction and replacement scouts

This round produced two exactly verified changes to the golden H20 chamber: a coherent sixty-cap layer, followed by a radius-four replacement cutting a different gap. Both final chambers have **587 actual facets and an exact strict escape**, hence infinite volume. The remaining searches below are numerical proposal and finite-bank tests, not new finite-volume examples.

The exact outcomes are documented separately:

- [Sixty quarter-face caps](fresh_h20_intrinsic_caps_20261014_quarter_face_exact.md), with all 33,390 new pair relations checked exactly, 252 retained finite vertices, and a new eighth-vertex escape.
- [Radius-four replacement](fresh_h20_intrinsic_caps_20261014_radius4_replacement.md), which replaces one old radius-one cap, creates an actual finite endpoint, and cuts the annular escape. The eighth-vertex escape remains.
- [Uniform annular obstruction](fresh_h20_intrinsic_caps_20261014_annular_gate.md), covering all finite selections from the specified small-radius families, including the continuous central-contact family below.

## Quarter-face construction

The seed alcove is the product of affine A11, A2, and E6 simplexes with inradii `1,phi,1`. For each factor, choose a set of omitted-vertex indices and take its affine-mark-weighted barycenter. The numerical bank uses:

| Factor | Subsets |
| --- | --- |
| A11 | 48 contiguous cycle subsets of sizes one through four |
| A2 | All seven nonempty subsets |
| E6 | All 35 nonempty subsets whose mark sum is at most four |

Divide the product face center by four and use radius `3/2`. Of 11,760 recipes, 2,388 cover at least one quarter-vertex probe. All 2,388 pass the numerical seed and frozen 527-wall compatibility gates. Seven lazy MILP rounds select **60 mutually compatible caps** covering all 252 quarter probes; all are strictly covered numerically. There are twelve right pairs and 1,758 dotted pairs among the selected caps.

The [scout certificate](fresh_h20_intrinsic_caps_20261014_quarter_face_scout.json) stores every selected recipe. [Its source](../scripts/fresh_h20_intrinsic_caps_20261014_quarter_face_scout.py) reconstructs factor-local omitted-vertex indices and weights explicitly. The separate exact certificate verifies only the selected sixty, not the whole numerical bank or solver optimality. The [numerical geometry replay](fresh_h20_intrinsic_caps_20261014_quarter_face_geometry_scout.json) also identifies the first old vertex divided by eight as a remaining gap; the exact construction uses a square-`1/2` spacelike lift there.

## Central contact and the annulus stop

For a face center p, take center `gamma*p` and radius `R=2(1-gamma)`. This preserves the retained seed contacts at labels three and five. A prescribed finite or parallel contact with the central radius-two sphere gives the equation

```text
(p²-4) gamma² + 8(1+c) gamma - 8(1+c) = 0,
```

where `c=cos(pi/m)` or `c=1`. The positive root is solved directly, without a parameter grid.

The [central-contact scout](fresh_h20_intrinsic_caps_20261014_central_contact_contiguous_scout.json) tests central labels 2 through 20 and parallel contact on the same 11,760 face rays. It finds 127,800 seed-compatible recipes, of which 8,262 retain all 587 walls numerically. **None of those 8,262 cuts any of the 252 eighth-vertex probes.** Finite central contact alone was therefore not counted as closure progress.

A [joint replacement search](fresh_h20_intrinsic_caps_20261014_joint_contact_scout.json) instead fixes the first 527 walls and allows the quarter layer to change. It uses 5,304 candidates: 2,388 quarter caps and 2,916 useful central-contact caps with central labels 2, 3, or 5. It requires coverage of all 504 quarter- and eighth-vertex probes. After forty rounds and 6,642 accumulated pair conflicts, the search is **incomplete**, not infeasible. Individual solver calls were limited to ten seconds; no mutually compatible cover was found.

The exact annular theorem then supplied a stronger geometric stopping rule. The point `3*z_first/4` has Euclidean norm between 13 and 14. Vertex-centered balls of radius at most one lie beyond it; half-alcove radius-one balls, quarter-alcove radius-`3/2` balls, and all central-contact balls with radius at most two lie inside it. Thus every finite selection from those families misses a common open region, regardless of replacement or compatibility choices. No further MILP was run on that closed route. This theorem does not change the archived solver status.

## Larger endpoint caps and coherent replacement tests

The [earlier endpoint bank](fresh_h20_intrinsic_caps_20261014_completion_scout.json) has 221,592 numerical entries. Its radii come from finite endpoint-profile and omitted-wall contact equations. Its finite-label and isolated-dihedral-tail sampling remain bounded; this is not a complete all-radius endpoint bank.

Against the first three-quarter-vertex gap, [79 entries](fresh_h20_intrinsic_caps_20261014_annulus_endpoint_scout.json) cut the target and retain the central, half-vertex, and selected quarter layers. Some conflict with original radius-one vertex caps, which may be replaced. The simple recipe `[transport 0, profile 1, radius 4, finite]`, endpoint-bank index 12, conflicts only with the first original vertex cap. It became the exact replacement linked above. Transport 0 has endpoint **type 2**, not type 0.

Testing all 252 three-quarter probes yields 25,092 cutters; [24,984 retain all interior layers](fresh_h20_intrinsic_caps_20261014_annulus_endpoint_all_scout.json). Two restricted radius sets were then tested with the full numerical compatibility graph, allowing all 252 old radius-one vertex caps to be optional:

| Radius set for new caps | New options | Total options | Conflicting pairs | Numerical model result |
| --- | ---: | ---: | ---: | --- |
| `4,6`, endpoint contacts only labels 3 and 5 | 2,448 | 2,700 | 129,636 | Infeasible |
| `4,2*sqrt(3)*phi`, endpoint contacts labels 3,4,5,6 | 2,016 | 2,268 | 48,528 | Infeasible |

The corresponding [R4/R6 artifact](fresh_h20_intrinsic_caps_20261014_endpoint_replacement_scout.json) and [R4/sqrt3 artifact](fresh_h20_intrinsic_caps_20261014_endpoint_sqrt3_scout.json) record solver status, sources, counts, and scope. These are floating finite-bank infeasibility results, not exact universal obstructions. The second radius was chosen from a coherent three-cap cover of three triangle-related three-quarter probes; it was not an arbitrary finer radius grid.

Both models demand **strictly negative ball power** at all 252 old vertex positions as well as the 252 three-quarter positions. Thus a cap whose endpoint profile is semidefinite and merely touches an old vertex is not counted as covering it. This is a search for strict removal of those old ideal directions; infeasibility does **not** exclude a finite-volume completion that retains some as ideal vertices and separately bounds their cusp links. No such bounded-link claim was made for the sampled semidefinite profiles.

No broader endpoint bank or additional radius sweep was launched after these two tests. All jobs are complete. The exactly verified positive changes remain useful construction checkpoints; the finite-volume objective is unresolved.

## Reproduction and validation scope

The scouts use NumPy and SciPy/HiGHS, with numerical angle tolerances. Recipe and dependency hashes are archived. Their source files are [quarter faces](../scripts/fresh_h20_intrinsic_caps_20261014_quarter_face_scout.py), [central contacts](../scripts/fresh_h20_intrinsic_caps_20261014_central_contact_scout.py), [joint replacement](../scripts/fresh_h20_intrinsic_caps_20261014_joint_contact_scout.py), [annular endpoint cutters](../scripts/fresh_h20_intrinsic_caps_20261014_annulus_endpoint_scout.py), [R4/R6 selection](../scripts/fresh_h20_intrinsic_caps_20261014_endpoint_replacement_scout.py), and [R4/sqrt3 selection](../scripts/fresh_h20_intrinsic_caps_20261014_endpoint_sqrt3_scout.py). The numerical geometry replay passes with:

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_intrinsic_caps_20261014_quarter_face_geometry_scout.py --verify
```

The exact quarter-face, annular, and radius-four reports provide their separate arithmetic replay commands and independent-review scopes. No solver success, failure, or finite probe cover in this report is substituted for those exact checks or for global volume certification.

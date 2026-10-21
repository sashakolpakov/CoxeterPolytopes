# Exact26-wall execution and pencil review

19 October2026. Root review PASS for the stated local algebra and written
pencil obstructions, with the inherited norm alphabet and H8 facet facts
explicit. This is not a finite-volume construction or a nonreflectivity
theorem. No Lean or TeX batch was needed for this review.

## What actually ran

The complete source of
`scripts/check_h19_twenty_six_wall_seed_20261019.py` was read before its
separate bounded generation and replay. Source SHA256:
`f0ac1d28225a37be17e9b61f2e436ea2fc3b1511efb32b02212b33539fecb363`.
The10,261-byte saved certificate is
`results/h19_twenty_six_wall_seed_20261019.json`, SHA256
`e310723817ce78bce5c8ed1c11f2a210925d6edc3c6fcb8faf443d5458fc011b`.
Both source and certificate remain unchanged.

The two execution ledgers are
`results/h19_twenty_six_wall_seed_generation_20261019.json` and
`results/h19_twenty_six_wall_seed_replay_20261019.json`.
They record PASS in0.229/0.225 seconds with process-tree peaks
56,836,096/57,147,392 bytes. Disk was183.05GiB free and available RAM
6.38/6.69GiB at launch. Each run held the shared lock under the
59-second/1GiB wrapper; no concurrent scientific process was started.
There were no generated build intermediates to remove.

The exact checks cover:

* Physical signature(9,1), positive-definite conjugate,26 integral roots
  and all their integral reflection covectors.
* All325 wall products:174 right,32 label3,23 label5 and96 strictly
  ultraparallel. All26 inequalities have actual-facet witnesses obtained
  by acute projection of a common strict timelike centre.
* The original finite vertex p, the new finite vertex p_s, their full
  active-wall lists and all other inequalities. Each active9-root Gram
  is positive definite.
* The former25-wall spacelike escape and its positive product with the
  new mirror s, which cuts that escape.
* All nine initial edges from p: seven first intersections are finite
  vertices; the first intersections on the old3 and old4 directions are
  spacelike and satisfy every retained inequality.

The two complete41-vertex H8 facets and their disjointness are inherited
from the explicit facet-identification proofs, not re-enumerated here.
Together with p and p_s, which are on neither facet, they give at least
84 finite vertices. A complete ambient vertex census is not claimed.

## Why the volume is nevertheless infinite

Take either checked feasible spacelike endpoint X and the checked strict
timelike centre c. The segment from c to X crosses the future light cone
before reaching X. Every point before X satisfies all26 wall inequalities
strictly, because c does. The crossing therefore has an open ideal-boundary
neighbourhood inside the chamber. Its hyperbolic volume is infinite.
Thus the two edge diagnostics already disprove finite volume; a global
cone census is unnecessary for this negative conclusion.

## Review of the subsequent all-height arguments

Root read the complete independent report
`results/compact_h9_discriminant19_edge_pencil_obstructions_20261019.md`,
SHA256 `41e78966ecb8137714347e79ddc5e434044409a85d5d98bcecd437350e8ace0a`.
Its three pure-pencil exclusions, mixed finite-angle obstruction and52
Schur-positive labeled profiles were checked directly. The profile count
is42 zero/label3 plus10 profiles with a label5 bond; removing the obstructed
zero profile leaves51. The displayed inverse-A contributions reproduce
every allowed and excluded combination.

Root also checked the subsequent report
`results/h19_omega3_short_pencil_completion_20261019.md`: modulo(D),
O/(D)=F19 and10t²=1 forces t²=2, impossible. This removes the height
cutoff for the remaining short omega3 pencil. The actual negative bounded
search and replay are recorded there separately from this paper proof.

The complete canonical reflection norm alphabet{2,2b} uses the primitive
determinant argument in the unchanged seed report: d divides the prime
determinant ideal, O is a PID and totally positive units are squares.
These hypotheses are explicit. Neither orthogonal-pencil exclusion
rules out a short root with nonzero finite contact to the fixed edge
walls. Such profiles remain a construction research problem.

## Reproduce the26-wall check

```
python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock -- scripts/check_h19_twenty_six_wall_seed_20261019.py --verify
```

Use the existing resource gates. This command rechecks the local algebra;
it is not an exhaustive root search, an H8 census, or a Lean proof.

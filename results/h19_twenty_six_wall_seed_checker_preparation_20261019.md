# Source preparation for the26-wall H19 exact check

19 October 2026. **SOURCE PREPARED, NOT EXECUTED.** Only an AST syntax
parse was run. No checker body, matrix calculation, scientific producer,
optimizer, Lean process or build was executed by this lane.

Source: `scripts/check_h19_twenty_six_wall_seed_20261019.py`,272lines,
12,402bytes, SHA256
`f0ac1d28225a37be17e9b61f2e436ea2fc3b1511efb32b02212b33539fecb363`.
The coordinator owns source review and any gated generation/replay.

## Mathematical/source review before execution

The pinned26th-wall report is
`results/compact_h9_discriminant19_nonpolar_pell_mirror_20261019.md`,
SHA256 `b817f3a362144053d8d0ff8ca4f71de2c56faecf8b357ed1751fb802b4b130c4`.
Its new exact pencil identities were independently checked on paper:

    u=10+17phi, t=12+19phi,
    7t^2-Du^2=1,
    Du-Mt=3-phi,
    7t-Mu=2(phi-1).

They give the norm-two root s, its sparse covector, all25 negative/zero
old products and its strict positive product with the old escape X.
The center pairing uses C^(-1)1/10 correctly and equals
-(808+1311phi)/10. The projected point p_s has the stated nine active
normals and diagram3A1+A6: deleting native E8 node1 leaves A1+A6.
There is no source geometry flaw in those statements. The two full H8
facets and their complete41-vertex censuses remain inherited, as in
the earlier independent polar-cap review.

## Exact coverage

Only standard-library arithmetic and the pinned pure field/inverse/inertia
helpers are used. Native a+bphi is mapped to the helper's sqrt5 basis by
Q(a+b/2,b/2). The body checks:

* native form inertia at both embeddings, all22 archived norm-two roots,
  the19 retained roots and all7 later roots r,k,w,a,d,z,s;
* all26 canonical norms and integral coordinates, every reflection
  covector divisibility criterion, and all7 displayed sparse covectors;
* ALL325 wall pairs, including exact short/short labels2,3,5 or
  parallel/ultraparallel alternatives, and strict mixed/long dotted gaps;
* the common center's norm, future orientation and strict26 products,
  followed by26 timelike future facet projections strict on all other
  wall inequalities;
* complete26-wall feasibility of p and p_s, their exact nine active
  normals, positive-definite active Grams, and the p_s square identity;
* full25-wall feasibility and positive square of X, the26th-wall cut,
  the exact norm identity and its rational field norm;
* the nine fixed initial-p directions, their first other-wall parameters
  and exact endpoint norm signs, as a local diagnostic.

The last item DOES NOT assume all nine edges close at finite vertices.
It records finite/ideal/spacelike first intersections. An open direction
is not silently relabeled a closed edge, and this nine-direction check
is not a final-cone, full-vertex or cusp census. Any resulting geometric
volume consequence requires the separate written interpretation.

Seven exact source/input hashes are checked, including the12KB native
prefix, the pure helper and all constructive reports. The7.5KB complete
H8 archive is hash-bound but its root completion and41-vertex census are
not rerun. No scientific producer or solver is imported. The certificate
records this inherited-only scope explicitly.

## Generation and separate replay

The checker owns NO inner flock. Use the already reviewed bounded wrapper
with `--parent-lock`; it supplies the single shared lock,6GiB/40GiB start
gates,3GiB/8GiB abort gates,59second deadline,1GiB process-tree RSS limit,
parent-death protection and bounded output. Expected source-level work is
one8x8 inverse and a few fixed10-dimensional rational loops, so execution
should take seconds rather than require a scientific batch; no timing or
RSS has yet been measured.

Generation argv, after the coordinator authorizes the slot:

```sh
python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock -- \
  scripts/check_h19_twenty_six_wall_seed_20261019.py --write
```

Separate replay argv:

```sh
python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock -- \
  scripts/check_h19_twenty_six_wall_seed_20261019.py --verify
```

The certificate path is `results/h19_twenty_six_wall_seed_20261019.json`.
Generation uses exclusive creation and refuses replacement. Replay
recomputes the complete same fixed check and compares canonical JSON
structures; there are no tuple/list or timestamp differences. Neither
mode invokes an optimizer, searches roots, edits prior evidence, or claims
finite volume, compactness, nonreflectivity, a complete census or a new
dimension ceiling.

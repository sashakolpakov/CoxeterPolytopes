# Compact RA4: source ledger for the consolidation draft

24 October 2026. This is the source ledger for the **same compact-volume
paper family**, not a seventh paper or a new numerical increment. The new
TeX is a complete translation of accepted source arguments; it is **not
yet independently reviewed, accepted as a new edition, or compiled**.
Its author also wrote the 24 October fresh source audit. That audit is
not a different-author review of this new translation.

## Draft and preservation

| Item | Bytes | SHA256 |
|---|---:|---|
| `manuscripts/revisions/compact_ra4_volume_constraints_consolidated.tex` | 78302 | `f51b8c9f5752a7c218b5191589e032333be31fe5017300a8f24d2872b04495ab` |
| Preserved canonical `manuscripts/compact_ra4_volume_constraints.tex` | 53508 | `afe807672f3a6bbf4bb2dd95d2861b60391d7c2d7c879c6d6382e749e7ad976f` |
| Preserved canonical `manuscripts/build/compact_ra4_volume_constraints.pdf` | 405039 | `920edadbbd7a973d70931e6ced9b8d0d5cfa12922d3b5dc0def3dc19ed2e8bf2` |
| Preserved C28 successor freeze, binding 27 scientific/control/evidence files | 5674 | `67344acf21c0e5424a0d21de62766baaac54d6283077e28ae07ba6496ba5ce0d` |

The last row is
`results/compact_ra3_c28_successor_20261020_freeze.json`. All 27 members
were freshly hash/size matched in the source audit; the full list and
four authenticated historical executions are in
`results/audit_compact_ra4_c28_fresh_source_inventory_20261024.json`,
SHA256 `7dcdc8e65fe6dada79d0dedf4820c814e554778ff5f241ae0c71859b5f2a137f`.
This draft neither rewrites those files nor changes the old checkers'
canonical-source pin. The old C28 two-unit-normal-fiber proof is retained
in its original source and certificate; its shorter actual-adjacency
consequence is used in the draft.

## Result and proof correspondence

The labels below are stable TeX proof homes, rather than references to
Markdown fragments a reader must assemble.

| Result or internal bridge | Draft proof home | Source translated |
|---|---|---|
| Actual normal gap, orthogonal intersections, exact outgoing opposites | `lem:gap`, `lem:footprint`, `lem:opposite` | Canonical TeX, normal geometry |
| Polygon capacity and its full Loebell equality lattice | `lem:polygon` and following equality paragraph | Canonical TeX, neighborhood closing |
| Facet capacity `f >= 2V+5` | `thm:capacity` | Canonical TeX, capacity proof |
| Exact excess/volume identities; neighboring excess; local 46/96 cutoffs | `eq:facet-excess`, `thm:neighbor-excess`, `cor:local-size` | Canonical TeX, excess and volume |
| Zero excess recognizes the regular 120-cell | `prop:zero` | Canonical TeX, whole-facet propagation |
| Complete positive-excess gap `W >= 17` | `lem:small-degree` through the proof of `thm:gap` | Canonical TeX, geometric and small-graph dispatch; proved once |
| Every maximum-height facet is dodecahedral; sharp maximum wall height | `thm:farthest-dodecahedron` | Farthest-dodecahedron source below, all five proof steps |
| Any two regular-length pentagon edges; actual finite symmetry realization | `lem:any-two-pentagon`, `lem:actual-symmetry` | Elementary identities and the C26/C28 symmetry bridges |
| C26 has no entire regular pentagon | `prop:c26` | C26 source below: full dual sphere, symmetry normalization and two edge witnesses |
| C28 has no entire regular pentagon | `prop:c28` | C28 source below: Loebell L7, both fullerene types, propagation, Lorentz coordinates and the shorter `C.H1=0` contradiction |
| Universal `f >= 65`, `Vol >= 49 pi^2/6` | `thm:global-lower` and its proof | C28 consequence, using the shared core and farthest D theorem |
| Fixed graph/arithmetic algorithms, halt and coverage statements | `sec:verification` | The two fixed checkers and their retained certificates |

The common method is actual-wall geometry, with all original inequalities
checked at projection points. The method branches only after it supplies
a dodecahedral facet: its small positive neighbor is then handled using
explicit three-dimensional combinatorics and metric uniqueness. No
arithmeticity, general four-dimensional classification, 120-cell minimum
or compact second-volume theorem enters the argument.

## Exact translated sources

| Source | SHA256 |
|---|---|
| `results/compact_ra4_farthest_dodecahedron_wall_radius_20261020.md` | `d082e4a918d4504a59ac200f32bbb141f1f2cb52d701006eaea7ba5b0e5bcc22` |
| `results/compact_ra3_twenty_six_vertices_no_regular_pentagon_20261020.md` | `1bbe6adc630ff2f44424f6ef3c573b3c71dcef4f25496def9b13fd35b2e4a03f` |
| `results/compact_ra3_twenty_eight_vertices_no_regular_pentagon_20261020.md` | `076864a5d836456e56ac377dffbcaa581a4738ec977dbf58a7905580d9838afc` |
| `scripts/check_compact_ra3_c26_graph_metric_20261020.py` | `47108a5ff3f732bf4b54d7893cc925f8e2a9423106c0ac812f690d355a9976e4` |
| `scripts/check_compact_ra3_c28_graph_metric_20261020.py` | `7356fa22ba8095e480c81cd2800f62832764925ffc6c0706be904ec8c0a6bf4e` |
| `results/compact_ra3_c26_graph_metric_20261020.json` | `3f5567187bdaa7f2f42c6d4064a8cb46d09f30c99cdf12f667be050b3261cd1b` |
| `results/compact_ra3_c28_graph_metric_20261020.json` | `2f9fe47902fd08af4e7627cc92dd87e4ea0dfc48a4400548c4b2c120557d9415` |
| Fresh written source audit: `results/audit_compact_ra4_primary_and_c28_fresh_20261024.md` | `29880f5978cdb91b95cb95b8fb7fc58357078f41ebd11557981fac7b168b7902` |

The fresh audit reconstructed these source deductions and authenticated
the frozen packet without executing scientific code. The older root and
different-author reviews are preserved in that packet. Their acceptance
supports the source results, not this new TeX translation or a worldwide
novelty claim.

## External statements and exact scope

The paper names the acute-chamber wall property, Lorentz geometry,
convex-polytope boundary and Euler facts, compact reflection-orbifold
Gauss--Bonnet, twelve-face dodecahedron recognition, and marked compact
nonobtuse three-dimensional metric uniqueness. The latter is used only
on an assumed **actual** compact polyhedron, never to construct geometry
from a partial Gram matrix. Its primary location is
[Roeder--Hubbard--Dunbar, Theorem 1.4, printed p. 832](https://aif.centre-mersenne.org/item/10.5802/aif.2279.pdf).

The sole finite-classification input in the small-neighbor branch is
[Brinkmann--Goedgebeur--McKay, Table 2, PDF p. 15](https://arxiv.org/pdf/1207.7010v2):
one fullerene graph at 26 vertices and two at 28 vertices. Its page 2
specifies cubic plane graphs with pentagonal/hexagonal faces and explains
the equivalence of abstract and embedded isomorphism here. The manuscript
constructs one and two explicit dual spheres, respectively. Comparing
their count with that external result is what proves exhaustiveness;
the repository scripts do not re-enumerate fullerenes or assume chemical
edge lengths. The official primary PDF text was retrieved in the fresh
audit. No visual-inspection claim rests on the failed screenshot request.

## Mathematical meaning of the fixed computations

The inputs are the displayed finite triangle lists and permutations,
and rational-pair field expressions for the displayed metric identities.
The graph procedure constructs every edge incidence and every vertex
link. Two triangles per edge, connected cyclic links, connectedness and
Euler characteristic two verify each stated sphere. It then checks every
triangle image for each claimed group element. These tests are complete
for the **given** complexes and actions. Their finite loops halt; generator
closure for C26 is a subset of the finite permutation group. No completeness
claim for a general graph search follows.

The C28 propagation procedure starts from each indicated orbit, adds
pentagons having two regular neighbors, and closes under the four-element
action. It is monotone on a twelve-element set and halts after at most
twelve strict enlargements. The two outer initial orbits reach all twelve
pentagons; the middle orbit remains `DEJK` and is resolved geometrically.

Field operations are exact on rational pairs `a+b sqrt(d)`, with `d=17`
or `5` and the positive real embedding. Rational comparison after squaring
settles opposite-sign terms, inversion uses the nonzero norm, exponent
halving terminates powers, and decreasing order terminates the determinant
recursion. This verifies listed identities, not universal geometry.

The retained outputs give C26 dual counts `(15,39,26)`, generated action
order 12, and witness values `(7+sqrt(17))/4` for squared cosh and
`(1+sqrt(17))/4` for cosh. C28 gives two dual counts `(16,42,28)`, actions
of orders 24 and 4, the stated orbit closures and the exact quartic roots.
The old longer fiber contradiction is also retained. The newer short
ending is a written consequence of the same normal coordinates and
checked nonadjacency, and is not represented as a newly executed worker.

## Retained execution and new execution status

Generation and separate full-certificate fresh-process replay are already
recorded for both fixed checkers. They compare complete certificates;
`--verify-code` alone is not such a replay. Those are historical executions,
not executions of this consolidation task.

| Historical record | Elapsed seconds | Peak process-tree RSS bytes |
|---|---:|---:|
| `results/compact_ra3_c26_graph_metric_generation_20261020.json` | 0.066 | 56033280 |
| `results/compact_ra3_c26_graph_metric_replay_20261020.json` | 0.065 | 55640064 |
| `results/compact_ra3_c28_graph_metric_generation_20261020.json` | 0.068 | 56754176 |
| `results/compact_ra3_c28_graph_metric_replay_20261020.json` | 0.067 | 56721408 |

If a substantive change requires replay, the coordinating agent can use
the assertion-safe public entry, one command at a time:

```sh
python3 -B scripts/run_proof_batch_20261024.py check -- scripts/check_compact_ra3_c26_graph_metric_20261020.py --verify
python3 -B scripts/run_proof_batch_20261024.py check -- scripts/check_compact_ra3_c28_graph_metric_20261020.py --verify
```

These are documented commands, **not new execution claims**. They retain
the shared exclusive lock, 40 GiB disk / 6 GiB available-memory admission,
59-second / 1 GiB / one-thread limits, process-tree monitoring and abort
floors. Existing successes need no rerun merely for a new edition.

This task launched **zero** scientific workers/imports, enumerations,
tests, Lean or TeX jobs. Only sources and small metadata were read or
written. No Lean module is represented as verifying this consolidated
geometric manuscript. Conditional internal kernels and script assumptions
must retain their named external and geometric premises.

## Excluded material and next gate

The optional global low-excess charge, uniform outgoing center cut,
homology, inball, and conditional central-collar refinements are excluded.
Their sources, certificates and history remain preserved outside this
draft. The longer C28 fiber proof likewise remains in its original source
and certificate. The compact minimum target `f+W >= 120`, a second-volume
gap and a double-of-the-120-cell equality statement are not proved here.

The new edition needs another author's full mathematical review and root
acceptance, then a bounded build and a receipt tied to its final hash.
The intended build route is a byte-exact promotion to a **new** publication
edition filename directly under `manuscripts/`, using the unchanged bounded
generic builder. The pinned old path remains untouched. A draft under
`manuscripts/revisions/` is intentionally not buildable by that frozen
builder; no extra builder or proof route is required. Build, publication
catalogue and any conditional Lean interface updates remain root gates.

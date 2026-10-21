# Changed incidences in the 15-facet H9 source

No compact H9 example was found. This round changes the incidence structure of the noncompact source `A16_14_1_C01`, rather than adding caps to the previous partial chambers. Two entire prescribed families are excluded exactly. A further bounded search releases the old parallel incidences and finds abstract incidence proposals with 216 candidate vertices, but none of the four proposals taken into the metric stage produces a rank-ten Gram matrix.

The source is the exact 15-wall Gram and affine-circuit data in the [previous surgery archive](fresh_h9_incidence_surgery_20261013.json). The [previous report](fresh_h9_incidence_surgery_20261013.md) explains its source identification. The present work does not revisit the old H9/14 multiwedge schedules or the T1/T2 four-angle deformation family.

| Stage | Scope | Result and status |
|---|---|---|
| Flag obstruction | All 72 previous minimum right/dotted patterns, with arbitrary additional right-to-dotted changes and arbitrary dotted distances | Exact exclusion |
| New finite contacts | All 11 minimum seven-pair right-contact covers, each labeled independently 3 or 4; arbitrary further right-to-dotted changes | All 1,408 profiles excluded exactly |
| Joint missing-face proposals | First 64 recipes in each of 12 starting patterns; three disjoint Lanner triangles, at most four added finite contacts | 768 numerical incidence scouts |
| Fixed old disjointness | Twelve selected joint recipes; old designated dotted incidences retained | Twelve solver-reported UNSAT cases; no exported solver proof |
| Released old disjointness | Ten selected recipes; old parallel and designated dotted pairs may become right or newly dotted | Nine recipes produce 16 incidence proposals, each with 216 candidate vertices and three dotted pairs |
| Bounded metric constraint generation | First model of the first four successful released recipes; additional finite contacts with labels 2 through 6 | 55 exact rank-11 rejections, 145 numerical vertex rejections; no model reaches dotted-distance solving |

The [exact family certificate](fresh_h9_incidence_repair_20261014_verify.json) and its [stdlib checker](../scripts/fresh_h9_incidence_repair_20261014_verify.py) passed generation and a saved-result replay. A separate root-agent [source/proof review and fresh replay](fresh_h9_incidence_repair_20261014_review.md) also passed. That review explicitly excludes the numerical metric scout.

For every one of the 72 first-stage patterns, replace all designated dotted entries by zero. The resulting 15-node finite-only Gram matrix is positive definite; the certificate records its rational LDL pivots. Hence every clique avoiding dotted pairs has a positive-definite Gram matrix. Further right-to-dotted changes only remove cliques. Any compact realization would therefore have a flag simplicial 8-sphere as its dual nerve. Such a sphere needs at least 18 vertices: the nine ridges of one nine-vertex facet require nine distinct outside vertices, since a repeated outside vertex would form a ten-clique. Fifteen facets cannot suffice. This argument allows arbitrary dotted distances and does not use a rank cutoff.

The second exact family keeps all original nonright finite labels and changes originally right contacts instead. Seventeen surviving affine circuits require at least seven such changed pairs. There are exactly eleven minimum covers; `(1,8)`, `(7,13)`, and `(13,14)` occur in every cover. The checker enumerates both labels 3 and 4 on each of the seven contacts. For each of the 1,408 profiles it examines all 631 nine-node subsets avoiding the original six parallel pairs. Positive definiteness is exact: a nonright cycle has nonpositive quadratic form on its cycle vertices, and each forest is diagonally congruent to a rational Gram with root norms 2 or 4, tested by Fraction LDL.

Every resulting potential vertex is removed by unique-ridge peeling. If an eight-node face is contained in only one remaining potential nine-node vertex, that vertex cannot belong to a closed compact nerve; delete it and repeat. All profiles peel completely. Additional right-to-dotted changes can only remove potential vertices, so cannot restore a compact realization. This theorem does not cover further finite contacts, other labels, or changing an original parallel pair into a finite contact.

The positive-search step therefore changes those excluded assumptions. It combines weakened old finite contacts with new finite contacts forming three disjoint Lanner triangles, then permits the old parallel incidences themselves to change. The [joint numerical recipes](fresh_h9_incidence_repair_20261014_lanner.json), [fixed-incidence SAT trials](fresh_h9_incidence_repair_20261014_joint_sat.json), and [released-incidence proposals](fresh_h9_incidence_repair_20261014_released_sat.json) preserve the literal choices. Nine of ten released recipes yield saved models; two second-model queries time out, and one recipe is reported UNSAT. These solver outcomes are exploratory, not completeness certificates.

For all sixteen saved models, the new [solver-free replay](../scripts/fresh_h9_incidence_repair_20261014_scout_verify.py) independently reconstructs the selected incidences from the original source and reuses the exact spherical-block predicate. It verifies precisely 216 spherical nine-node blocks, no spherical ten-node block, extension of every positive-definite face to a nine-node block, and exactly two such blocks above every positive-definite eight-node face. The [saved replay certificate](fresh_h9_incidence_repair_20261014_scout_verify.json) binds the final helper and input hashes. These are necessary combinatorial conditions only. No spherical-nerve classification, rank-ten realization, Lorentzian signature, or compact hyperbolic polytope is inferred from them.

The [metric-stage script](../scripts/fresh_h9_incidence_repair_20261014_metric.py) then retains each of four proposed 216-vertex incidences while allowing further currently right pairs to take labels 2 through 6. Individual numerical vertex tests select the allowed contact pool. Joint failures generate monotone vertex constraints. Each surviving assignment is tested for a nonzero eleven-node principal minor that avoids the three free dotted entries. The minors are reduced modulo 241 using roots of 2, 3, and 5; a nonzero residue proves the exact determinant is nonzero in `Z[1/2,sqrt(2),sqrt(3),sqrt(5)]`. Thus that assignment cannot have Gram rank ten, regardless of the dotted distances.

The [metric archive](fresh_h9_incidence_repair_20261014_metric.json) records respectively 13, 3, 29, and 10 exact rank failures, together with 18, 6, 107, and 14 numerical joint-vertex failures. The solver-free replay checks all 55 modular determinants afresh. All four bounded solver spaces become UNSAT before the 500-round budget; no nonlinear distance trial is needed. The contact-pool restriction and the 145 vertex decisions remain numerical, so these four UNSAT outcomes are not promoted to an exact unrestricted-family exclusion.

Reproduce the two exact family conclusions and the saved combinatorial/rank certificates with:

```sh
python scripts/fresh_h9_incidence_repair_20261014_verify.py --verify
python scripts/fresh_h9_incidence_repair_20261014_scout_verify.py --verify
```

The bounded metric stage can be rerun, with solver choices potentially differing, using:

```sh
python scripts/fresh_h9_incidence_repair_20261014_metric.py --cases 4 --rounds 500
```

The earlier numeric scout files are historical proposal-generation records; the new replay certificate is authoritative for the saved released models and their rank witnesses. In particular it does not rely on the intermediate helper hash recorded before the released-incidence option was added. The first two exact family results are independent of those scouts altogether.

All jobs for this checkpoint have finished. No new Gram realization or example remains pending. Further work would have to change the finite-contact choices or the proposed vertex incidence structure again, rather than claiming that these combinatorial models already define compact polytopes. Compact H9 remains an appropriate target: the arithmetic noncompact dimension-21 restriction does not rule out compact H9. This search imposes Coxeter angles and geometric constraints, not a prescribed arithmetic lattice.

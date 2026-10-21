# Complete weighted corner cones at the new face masses

**Completed exact certificate; generation and fresh replay both passed.** The complete nonnegative orbit-weight cones of the saved conservative compatibility graphs still leave a witness above12. The full bound changes from12.368567257872337 to

`5064933545222915 / 410583979722408 = 12.335925889381436`.

The earlier restriction to the named248-class geometry union, of which209 classes have positive current mass, changes from12.089357918752444 to

`18670750401222348 / 1548357527926405 = 12.058423241708674`.

There is no global optimization or dimension theorem here. The [source](../scripts/fresh_current_historical_weighted_corners_20261017.py) and [certificate](fresh_current_historical_weighted_corners_20261017.json) retain all5,053 current face masses and only increase vertex multipliers. The inputs are the frozen [full model](fresh_weighted_raw110_full_model_20261017.md), its [support diagnostic](fresh_weighted_raw110_full_model_20261017_support.md), and the saved physical-corner ancestry of the [raw110 package](fresh_raw110_mixed_actual_corners_20261017.md).

## Domain and provenance

The eight geometry packages in that ancestry contain350 distinct historical classes. Of these,248 have positive current mass:209 in the latest named248-class union and39 additional classes available in older packages. The two occurrences of248 refer to different sets. The selected positive domain consists of61 raw23,23 raw24,42 raw25,4 raw33,5 raw39,19 raw87,22 raw88,21 raw89,3 raw110,16 raw178,8 raw295,14 raw296 and10 raw297 classes.

Every selected face is checked against the latest explicit-provenance or unique-context table: face index, original provenance when present, source context, canonical core and raw family must match. No merged multi-context row is assigned an arbitrary representative. The raw core order plus the spectator rank is15, and every saved literal prism corner is retained. A d-prism has2d vertices, whereas its face demand is `d(d−1)/2`; the source checks these as separate quantities. The older raw24 package lacks explicit rank fields, so those fields are reconstructed from the consistent frozen raw-family order tables and the actual spectator names.

At each canonical finite forest, a full star is identified by its unordered native external contact pair. Repeated presentations must also agree in occupied port pairs, both squared timelike heights and the native omitted product. Complete orbit node sets may coincide, in which case their identical actual occurrence records are merged. Partial overlaps are rejected. Actual corner origins are counted once, and their complete union must equal the selected whole-face corner list. Coefficients count actual prism vertices, never the number of orientation embeddings.

The result contains550 forests,2,072 complete orbits and35,064 literal stars. The2,072 actual corner occurrences all pass the complete occurrence accounting.

## Pair ownership and conservative compatibility

The source order is a deterministic preorder of the raw110 ancestry, preserving each producer's recorded predecessor order. A source owns a pair whenever it contains both full native star keys. Ownership is determined before reading the pair's outcome. Thus an earlier unknown result cannot be replaced by a later exclusion. Recursive predecessor transport checks all native star data, injectivity, reused pair totals and reused exclusions. Direct exclusions must lie outside every inherited mapped-pair domain.

Across the expanded domain there are4,856,393 literal pairs:

- 4,369,518 excluded pairs, justified by frozen physical certificates;
- 100,189 inherited unknown pairs, which remain compatible;
- 386,686 pairs with no saved common source, which also remain compatible.

No new contact projection, algebraic field test or physical pair census is performed. In particular, a compatible set in this graph need not be realizable by actual walls. This direction is deliberate: every actual simultaneous set remains graph-compatible, which is sufficient for every exported upper capacity.

## Complete cone certificates

Let the complete corner orbits at a forest be `O_i`, and let `d_i` be their current actual occurrence demands. Each literal compatible clique C gives an orbit-count vector `N(C)`. For any nonnegative orbit charge w, its graph capacity is `max_C N(C)·w`. It suffices to enumerate maximal cliques and discard coordinatewise dominated count vectors.

The exact covering problem is

`min Σ λ_C`, subject to `λ_C≥0` and `Σ λ_C N_i(C)≥d_i`.

Its dual maximizes `d·w` subject to `w≥0` and `N(C)·w≤1` for every clique. Every saved forest has both a feasible rational charge and a nonnegative rational covering mixture of exactly the same value. Their equality proves the optimum for the entire nonnegative orbit-weight cone at those fixed demands, not merely for the exported charge list. Averaging whole cliques under the common forest automorphism group gives the corresponding uniform literal orbit marginals when desired; independent choices of orientations are not used.

The full and restricted demands are recomputed separately. In the restricted calculation the39 additional historical classes have zero mass. The same full clique vectors remain sufficient: every clique in the induced retained domain extends to a maximal full clique, and unused orbit charges may be zero.

Enumeration found1,429,288 maximal literal cliques,1,954 distinct count vectors and1,922 nondominated vectors. Exact certificates use1,561 full covering terms and1,405 restricted covering terms. Every forest completed; no resource guard or recovery case remains open. Local LPs only propose rational certificates during generation. Saved verification enumerates the graph cliques again and checks every exact inequality, marginal coverage and equality without an optimizer.

## Repairs and preserved constraints

The full certificate requires75 vertex increases; the restricted certificate requires53. There are78 distinct exported weighted columns, including columns needed by only one of the two witnesses. Each column has vertex cost1 and exact nonnegative face coefficients obtained from its orbit charges and actual occurrence multiplicities. Its load is recomputed at both sets of face masses.

The repaired normalizers are respectively

`17107665821767 / 17062500000000`

and

`309671505585281 / 308880000000000`.

Each repaired multiplier is the maximum of its frozen baseline and the corresponding exact cone value. All439,410 inherited finite-model columns remain valid, including96,176 zero-budget signed columns and52 positive-vertex-cost signed columns. Their variables are nonnegative; signed face coefficients do not make the dual constraints equalities. With face masses fixed, increasing a vertex multiplier can only relax a positive-cost dual inequality and leaves a zero-budget inequality unchanged. The independently replayed model/support certificates supply the baseline validity.

No old-mass implicit simplex-packing or joint-port cover is imported as a new-mass certificate. The conclusion is restricted to the represented corner graph cone just proved. At the fixed current face masses, adding any further nonnegative orbit charge from this same graph cannot invalidate the repaired witnesses: the saved covering mixtures bound every such charge simultaneously. Stronger physical pair information, additional actual-face coverage or constraints coupling different incidence domains can still strengthen the model.

## Validation scope

Generation and a separate fresh saved replay both terminated exit0/PASS, covering all550 forests with zero open cases. The replay performed no optimization. A separate standard-library Fraction check, without importing any producer, reconstructed every demand, all charge inequalities against the1,922 saved maximal vectors, all2,966 covering terms, all78 column loads and both multiplier repairs and ratios. This saved-data check did not re-enumerate the compatibility graphs or physical wall pairs. All34 direct source hashes, Python syntax and local report links pass. The inherited exact local charge/cover recovery is reused from the [earlier weighted-cone producer](fresh_weighted_actual_corner_capacities_20261017.md); all demands and certificates here are new.

An independent peer gave a narrow full source/mathematical read of the complete-orbit merge, ordered ownership including unknowns, separate full/restricted demands and exact cone accounting. Its final delta read also confirmed the two corrected preparation guards described above. This is not an independent graph, physical-pair or LP census. Frozen upstream geometry and its recorded reviews remain explicit inputs.

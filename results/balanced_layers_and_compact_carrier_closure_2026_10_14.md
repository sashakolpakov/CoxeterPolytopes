# Balanced cap layers and a complete compact-carrier exclusion

**14 October 2026. No new finite-volume or compact Coxeter polytope is
established.** This checkpoint records exact constructions, exact reasons
they still have infinite volume, and a complete exclusion of one specified
compact-H9 incidence family. It continues the
[preceding balanced-cusp checkpoint](balanced_cusp_and_bugaenko_surgery_2026_10_14.md).

## Literature determines the targets

The [fresh primary-literature check](fresh_h20_balanced_boundary_20261014_literature.md)
confirms the distinction requested by the user:

- Noncompact arithmetic reflection groups occur only in dimensions at most
  19 or in dimension 21. Arithmetic H20 and arithmetic H22+ are excluded.
- General arithmetic reflection groups and all compact Coxeter polytopes
  have upper bound 29. The dimension-21 theorem does not exclude compact
  arithmetic examples in dimensions 22 through 29.
- The sources checked retain compact dimension 8 and noncompact dimension
  21 as the known records; dimension 20 remains absent from the known
  finite-volume dimension list.

Allcock's published edgewalking algorithm takes a rational Lorentzian
`Z`-lattice as input. It is not a ready-made arbitrary nonarithmetic
combinatorial enumerator. A separate incidence-and-Gram search can allow
nonarithmetic realizations, but must establish exact metric realization and
global finite volume. The active targets remain compact H9 and nonarithmetic
finite-volume H20. No arithmetic lattice search is warranted in the excluded
noncompact dimensions.

## H20 constructions and their obstructions

All three branches use the intrinsic affine `A6 × A6 × A7` cusp rather than
an arithmetic host. Their pair certificates prove coherent actual facets;
their strict spacelike escapes prove infinite volume.

| branch | actual facets | useful construction result | failure certificate |
|---|---:|---|---|
| rational integer-gap layer | 810 | preserves the old 418 facets and cuts all 392 `2z/5` probes | strict `3z/5` escape; stronger whole-family obstruction |
| unit-inradius contact chains | 1,206 | two full outer layers, seven inner caps, all 392 `2z/5` probes cut | common escape even for all three complete chain levels |
| golden third-factor chains | 815 | all 392 `7z/20` probes cut with the old endpoint layer retained | strict inner escape and explicit deficient cusp links |

The [integer-gap report](fresh_h20_balanced_integer_20261014.md) includes
all 327,645 exact pair checks, 810 facet witnesses, 392 retained finite
endpoints and the literal rational wall list. A
[separate full chamber reconstruction](fresh_h20_balanced_integer_exact_20261014_review.md)
verifies these results without producer imports. A separate exhaustive integer
calculation excludes finite-volume completion by **any** radius-two cap with
nonnegative integer seed gaps while retaining the full old endpoint layer.
It considers every weak composition, includes tangencies, and checks all
possible finite Coxeter labels against the old endpoint walls. Every
compatible cap leaves the same norm-`1/1400` strict escape.
The [independent review](fresh_h20_balanced_integer_obstruction_20261014_review.md)
checks the proof and reconstructs all 1,110 archived cutter recipes with
separate rational matrix arithmetic.

The [unit-chain report](fresh_h20_balanced_boundary_20261014.md) uses
`k=(33−sqrt(65))/32` and caps centered at `(1−k^j)z` with radii `2k^j`.
Seven first-level caps and every second- and third-level cap give 1,206
actual facets. Exact field inequalities handle all 640,710 new wall pairs.
Even the relaxed bank containing every first-level cap leaves a norm-`1/8`
strict escape at `3z/5`. The fourth same-ray cap conflicts with the retained
endpoint facet. This closes that fixed three-level dictionary, not all
possible chain variants.

The [golden-chain report](fresh_h20_balanced_golden_20261014.md) also checks
417- and 479-wall intermediate configurations. The 479-wall branch replaces
the uniform endpoint radius and permits eight-step chains; it is explicitly
distinct from simply adding walls to the old baseline. The
[separate symbolic review](fresh_h20_balanced_golden_20261014_review.md)
verifies the construction identities, exact escape and cusp-rank analysis.
It does not independently rerun every pair classification.

Deficient cusp links show that a **current** chamber has infinite volume.
They are not automatically obstructions to adding walls: crystallographic
finite components can have local affine completions. Any proposed global
completion still needs compatible new walls, all new edges and cusps checked,
and the absence of every spacelike escape. Covering finitely many sampled
probes or increasing the facet count does not establish that condition.

## Compact H9: a complete incidence family is excluded

The [compact-search report](fresh_h9_bugaenko_closure_20261014.md) begins with
the known compact Bugaenko H8 nerve, wedges a facet to dimension nine, and
makes two ordered face truncations of codimensions two through five.
For the central-singleton wedge carrier, all 15,790 operation orbits are
excluded: 15,582 have a missing face too large for a compact Coxeter
polytope, and the remaining 208 recipes represent 128 distinct incidence
types, all exactly excluded.

This has no field, dotted-distance or finite-label cutoff. Necessary Lannér
assignments and spherical vertex diagrams are propagated exactly; disjoint
orthogonal negative blocks contradict Lorentzian index one. Selected fixed
minors give additional exact rank obstructions. A changed assignment order
settles the last central-carrier type in 126 states. Node-limited searches
in other incidence families remain unresolved.

Across the frozen runs, 642 of 650 tested catalogue types are exactly
excluded. The eight unresolved cases are `22,23,318,582,629,633,636,637`;
2,303 catalogue types have not been tested. The final first-carrier sweep
alone excludes 502 of its 505 types, leaving its three limited cases open.
These counts describe this 2,953-type surgery catalogue, not all compact
H9 incidence structures.

The formerly unresolved case zero has a separate short certificate: its
864 necessary partial labelings form six symmetry orbits of size 144, and
every representative forces two orthogonal negative directions. An
independent arithmetic checker reconstructs the orbits and verifies the
negative blocks and forced cross-pair orthogonality. The aggregate central
closure also has a separate reconstruction of its 208 actual surgery
recipes and relabeling witnesses.

## Frozen verification and next direction

The individual reports state which checks are producer replays, independent
arithmetic checks, or numerical proposals. The
[archive auditor](../scripts/audit_balanced_layers_checkpoint_20261014.py)
checks local report links and recursively declared source hashes; its
[saved result](balanced_layers_and_compact_carrier_closure_2026_10_14_audit.json)
is a provenance check, not a substitute for the mathematical certificates.

Do not enlarge the excluded integer-gap or fixed three-level chain banks
while retaining the same baseline. A useful next H20 construction must
change endpoint walls, cap geometry or the seed, with global closure as the
objective. In compact H9, skip the fully excluded central carrier and use
the sharper assignment order on genuinely unresolved incidence types.
The mathematical existence problem remains open in this work.

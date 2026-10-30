# Compact arithmetic dimension21: proof and reproduction map

For current reading, see the [current dimension-bounds manuscript](../arxiv/hyperbolic_coxeter_dimension_bounds/main.tex) and the
[submission guide](arxiv_submission_editions_20261028.md). The supporting
arguments and recorded verification scopes below remain distinct from those
current submission builds. Standalone draft typesetting routes have been retired.

This accompanies supporting argument. It translates the completed internal computer-assisted proof of compact arithmetic dimension at most21. It is not a new catalogue run or a proof that compact dimension21 is attained. The broader arithmetic bound21 is attained by a published **noncompact** example. External refereeing, priority, and complete proof-assistant verification are separate matters.

Frozen reports retain their historical “candidate” language. The final internal mathematical review found no missing global hypothesis or unresolved mathematical step. The new manuscript keeps the published assumptions, internal geometric proofs, and finite computational obligations distinct.

## Dependency chain

```text
published compact ceiling29 + fields at n>=22
             |
             +-- Q(sqrt2): complete Lanner/union/all-white proof
             |              + whole-component weights <=2n-2 -> bound20
             |
             +-- cubic-seven: complete high-edge/all-white proof
             |                 + whole D5 halves <=2n-2 -> bound20
             |
             +-- Q(sqrt5): complete arithmetic stopping domain through31 roots
                            + exact dimension-specific angle rules
                            + all common spectators + all finite forests
                            -> exclude every n=22,...,29 -> bound21
             |
             +-- exactly24 compact field/dimension branches excluded
```

The common analytic implication is a strict lower bound on average vertex weight:

`beta(n)=n(n−6)/8` for even n and `(n−1)(n−5)/8` for odd n.

It follows from the nonnegative symmetric h-vector and the face inequalities triangle≥2, quadrilateral≥1, k-gon≥5−k. The strictness comes from h0=hn=1. It remains strict when an upper certificate equals beta, as at n29.

## Published assumptions suitable for an explicit formal interface

| Suggested interface name | Precise content and use | Source |
|---|---|---|
| `CompactCoxeterSimple` | A compact acute Coxeter chamber is simple; finite corners have positive-definite Gram matrices. | Vinberg, §§3.1–3.2/Propositions7–9, as pinpointed in the [global audit](../results/fresh_arithmetic_dimension21_hostile_20261016.md). |
| `FiniteAndAffineCoxeterClassification` | Complete finite connected diagrams, positive affine kernels, and the Lanner order bound5. | Vinberg §2.4 and the standard finite/affine classification. |
| `ArithmeticCyclicCriterion` | Integral cyclic products in the actual ground field, totally real entry field, conjugate positive-semidefinite Gram forms. | [Nikulin §2](https://arxiv.org/pdf/0708.3991), recalling Vinberg's criterion. |
| `LannerGroundField` | A principal arithmetic Lanner block has the actual ground field used by the high-field reduction. | Vinberg/Nikulin input named in the [cubic proof](../results/fresh_compact_global_bound_20261015.md). |
| `SimplePolytopeHVector` | Symmetry, nonnegativity, h0=hn=1 and f/h identities for the simplicial dual of a compact simple polytope. | Standard h-vector theorem; weighted-face application in Vinberg and the [independent weight proof](../results/fresh_sqrt2_weight_bound_20261015_review.md). |
| `CompactCeiling29` | No compact hyperbolic Coxeter polytope in dimensions≥30. | [Vinberg Theorem1](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt). |
| `ArithmeticFieldsAbove21` | At dimensions≥22, the possible arithmetic ground fields are Q, Q(sqrt2), Q(sqrt5), Q(cos(2pi/7)). This restriction is **not** applied at21. | Vinberg Theorem2; [Nikulin introduction](https://arxiv.org/pdf/0708.3991). |
| `ReflectionLatticesSimplestType` and `MeyerIsotropy` | Arithmetic reflection lattices are of simplest type; a rational indefinite form in at least five variables is isotropic, excluding the rational field in the compact high-dimensional range. | [Belolipetsky Lemma2.1](https://arxiv.org/pdf/1506.03111), with the standard Meyer theorem. |
| `NoncompactArithmetic21` | Noncompact arithmetic reflection groups have dimension≤21 (and not20). Used only for the optional all-arithmetic corollary. | [Esselmann](https://doi.org/10.1006/jnth.1996.0141). |
| `BorcherdsExample21` | A noncompact arithmetic reflection example exists in dimension21. Used only for sharpness of the combined class. | Borcherds example, described in [Belolipetsky §7](https://arxiv.org/pdf/1506.03111). |

This inventory is explicit about the standard theorems admitted from the literature. It is not a claim that all these theorems already have Lean implementations in this workspace.

## Internal arguments that must not be hidden as published axioms

| Internal statement | Mathematical proof / finite obligation |
|---|---|
| Elliptic subsets are actual faces | The manuscript proves this by projection from a strict interior timelike vector and nonnegativity of the inverse positive-definite Z-Gram. It also explains affine compactness exclusion. |
| Complete quadrilateral star has corank1 | A corner supplies k−2 positive directions, a Lanner subdiagram supplies one negative direction, and n+2 roots lie in dimension n+1. |
| Arbitrary white support is covered | A first touching white node enters a complete contact domain; subsequent nodes can be ordered along their white component. Orthogonal spectators are retained. Every singular state can stop, even if it can also extend. |
| Whole-component use of bonuses | Conditional sealing or an all-size equality argument proves the marked factor is entire. An embedded proper subdiagram is insufficient. |
| Dotted parameter finiteness | A common Lorentz block forces a positive-semidefinite endpoint Schur complement. Physical and conjugate bounds plus integral cyclic products give a finite integer box. This is not sampled distance coverage. |
| Cross-component additivity | Exact source multiplicities a,2b,3c cancel their actual denominators. The three destination exclusions make the families disjoint. Capped counts select coefficients only. |
| Finite-forest optimization | Fix exact hub counts, remove those types from the residual alphabet, and prove the unbounded-coin recurrence by removing/appending a last component. |
| Strict global contradiction | Double count planar angles and pair h-vector coefficients. No compact-face argument is substituted for whole-chamber compactness. |
| Field-specific bounds | They follow from the preceding arguments **and** complete exact finite checks below. Importing these bounds as axioms proves only the final composition, not the entire internal proof. |

## Exact finite certificate interfaces

All links below are frozen proof inputs. The manuscript is a translation of these inputs, not an alteration of them.

### Q(sqrt2)

- [Complete Lanner catalogue](../results/fresh_sqrt2_lanner_catalogue_20261015.md): finite alphabet2,3,4,6,8, nine triangles and two cycles, plus the full dotted parameter family.
- [Complete quadrilateral/all-white proof](../results/fresh_sqrt2_quad_exhaustion_20261015.md) and [independent reviewer](../scripts/fresh_sqrt2_quad_exhaustion_20261015_review.py): 51 marked pieces,6,010 identifications,115 literal supports/52 classes, complete one-white domains, and the all-size pendant/connecting-component proof.
- [Independent weight proof](../results/fresh_sqrt2_weight_bound_20261015_review.md) and [checker](../scripts/fresh_sqrt2_weight_bound_20261015_review.py): literal whole-component corrections, all-rank upper2n−2 and strict exclusion at n≥21.
- [Dotted-tail review](../results/fresh_sqrt2_dotted_tail_20261015_review.md): the unbounded pendant classification used by the full quadrilateral proof.

### Cubic-seven

- [Full proof](../results/fresh_compact_global_bound_20261015.md) and [independent reviewer](../scripts/fresh_compact_global_bound_20261015_review.py): all35,904 binary union supports, nine structural survivors, arbitrary-white reduction, exact cubic determinant coefficient identities, and entire D5 half payments.
- The inequality `det H≤1−delta²−epsilon²≤1/2≤8L` meets the singularity equality `det H=8L`, so no eliminated external branch can survive. This is the geometric basis for the whole-D5 payment.

### Q(sqrt5)

- [Complete Lanner catalogue](../results/fresh_sqrt5_lanner_catalogue_20261015.md) and [independent catalogue review](../results/fresh_sqrt5_lanner_catalogue_20261015_review.md).
- [No-short-corner minimal unions](../results/fresh_sqrt5_quad_exhaustion_20261015.md) and [independent complete traversal](../scripts/fresh_sqrt5_quad_exhaustion_20261015_incidence_review.py):5,440 literal supports/1,180 marked classes. The later unrestricted8,298-class census is not required as a substitute: the bound21 rule never removed the persistent short-path baseline.
- [Complete original coverage argument](../results/fresh_sqrt5_dimension27_coverage_20261015.md): initial arithmetic filters, conditional seals,147 singular first attachments,51 frontier states and17 closure states, with all literal continuation maps.
- [Exact dotted parameter review](../scripts/fresh_sqrt5_dotted_parameters_20261015_review.py):24 families,532 integral candidates,32 seeds.
- [Complete dotted queue](../scripts/fresh_sqrt5_dotted_rank29_exhaustion_20261015_review.py):438 states,264 singular obligations, all423 below-boundary domains. The31-node ceiling is geometric.
- [Complete finite queue](../scripts/fresh_sqrt5_nonsingular_rank29_fast_20261015_review.py):5,172 states,51,599 contact profiles,2,070 singular obligations, including expanded states. The separate arithmetic implementation and its fresh replay are recorded in its [report](../results/fresh_sqrt5_nonsingular_rank29_fast_20261015_review.md).
- Rank22 rational weight table: `variable_order` indexes `certificate.exact_primal`; the last coordinate is the bound, preceding coordinates are bonuses. All30,169 bonus variables are nonnegative rational numbers. The literal orbit records define the meaning of every coefficient.
- [Independent rank22 rule/forest checker](../scripts/fresh_sqrt5_dimension22_h3_cross_20261016_review.py):71 finite types,410 exact hub contexts, explicit graph-orbit multiplicities and reverse-type recurrence; upper21977/500.
- [Final golden-field coverage](../scripts/fresh_sqrt5_h3_hub_coverage_20261016.py): all2,661 obligation identities,1,634 feasible rows,91,726 common spectator triples; minimum quadrilateral payment560003/560000.
- Earlier dimension-specific complete rules are [28–29](../results/fresh_sqrt5_dimension27_coverage_20261015.md), [25–27](../results/fresh_sqrt5_dimension24_coverage_20261015.md), [24](../results/fresh_sqrt5_dimension23_coverage_20261015.md), and [23](../results/fresh_sqrt5_hub_count_coverage_20261016.md). The rank22 coefficients are not extrapolated to these dimensions.

### Final composition

- [All-field composition](../scripts/fresh_arithmetic_dimension21_composition_20261016.py): checks all24 field/dimension branches22–29 against the complete field bounds and recursively binds the finite certificates.
- [Standalone scope checker](../scripts/fresh_arithmetic_dimension21_scope_check_20261016.py): standard-library exact bookkeeping, all2,661 identities,91,726 triples,410 contexts, and the24 branches. It does not replay every Gram/queue calculation.
- [Hostile mathematical audit](../results/fresh_arithmetic_dimension21_hostile_20261016.md): actual-face geometry, rank bridge, kernel/stopping completeness, lower-rank reuse, shared spectator quantifiers, strictness and field scope. This is a written proof audit, not a second whole-catalogue execution.

## Reproduction without mutating frozen certificates

Use Python3 with the repository's numerical/symbolic dependencies installed. The two quick final checks below use only standard-library code; the full field/graph replays may require SymPy and NetworkX through their explicitly bound helper modules. No optimizer is needed for saved `--verify` modes. Do not run generation mode over frozen output files.

Quick final identity/scope reproduction, from the repository root:

```sh
python3 scripts/fresh_arithmetic_dimension21_composition_20261016.py --bound 21 --verify
python3 scripts/fresh_arithmetic_dimension21_scope_check_20261016.py --verify
```

Full independent finite-check entry points, also in read-only verification mode:

```sh
python3 scripts/fresh_sqrt2_quad_exhaustion_20261015_review.py --verify
python3 scripts/fresh_sqrt2_weight_bound_20261015_review.py --verify
python3 scripts/fresh_compact_global_bound_20261015_review.py --verify
python3 scripts/fresh_sqrt5_quad_exhaustion_20261015_incidence_review.py --verify
python3 scripts/fresh_sqrt5_dotted_parameters_20261015_review.py --verify
python3 scripts/fresh_sqrt5_dotted_rank29_exhaustion_20261015_review.py --verify
python3 scripts/fresh_sqrt5_nonsingular_rank29_fast_20261015_review.py --verify
python3 scripts/fresh_sqrt5_dimension22_h3_cross_20261016_review.py --verify
python3 scripts/fresh_sqrt5_h3_hub_coverage_20261016.py --verify
```

These full replays can be substantial. They verify their stated layers; they do not by themselves replace the written all-size geometric proofs or published inputs. Earlier higher-dimensional coverage is recursively bound by the final golden-field certificate. To replay each of those changed rules directly, use its linked producer's documented `--verify` command.

Build the manuscript without placing auxiliary files among frozen results:

```sh
mkdir -p manuscripts/build
```

## Concrete next reductions, not current theorems

1. **Compact arithmetic21→20.** The later [single-field reduction](../results/fresh_compact21_single_field_20261016.md) leaves Q(sqrt5) at dimension21. The [current exact incidence model](../results/fresh_raw357_coface_full_model_20261018.md) still has dual lower12.331381486554207>12. A successful next step requires new valid incidence/coface coupling or stronger geometric exclusion, followed by complete shared-budget verification; changing masses in the same deficient matrix is insufficient.
2. **General compact29→28.** At n29, beta=84 while the classical vertex allowance3n is87. Extract the exact maximizing rank29 forests under the existing rule, find a geometrically justified saving of at least3, and check all current triangle/quadrilateral outputs under it. This is separate from the arithmetic field restriction.
3. **Compact H9 construction.** The exact66-wall seed retains all original facets and205 finite vertices but is proved infinite-volume. Its new four-finite-duty quadrilateral-cap family fails the necessary norm equation for all seven profiles and all finite labels. The appropriate next construction changes the incidence, for example a split cap; it must retain the full native, center, facet, finite-vertex and escape obligations. No new compact or finite-volume example follows from those local seeds.

## Translation validation record

The TeX and this map are newly authored from the listed frozen proofs. On 2026-10-18 the manuscript compiled successfully with `latexmk` (terminal run `6c4e78`). The read-only composition check passed with 93 recorded interfaces (`1c954d`), and the independent scope check passed with 91,726 common spectator triples and 24 field/dimension branches (`d0ed37`). A separate map audit checked all 33 local links and the documented verification entry points (`69621c`). These checks did not rerun the full historical geometry queues during translation. The Lean formalization status is owned by the root formalization lane and must be stated at its actual checked scope.

# Heavy-centre forest: conditional scalar and incidence kernel

23 October 2026. **PREPARED SOURCE ONLY: uncompiled, zero actual axiom
outputs.** Eleven declarations request separate prints. Full different-author
source/map review and root preclear are required before any guarded compile
and separate actual audit. Preparation proves no geometric instance.

Source: `formal/CoxeterProofs/HeavyCentreForestScalarKernel.lean`, SHA256
`a839befdbd0454d70fd67cc458641469524dce1e69bdb1b4ba89e777c4fa2060`,
12329 bytes/263 lines. Only `Std` is imported. Existing modules, objects,
papers, guides, runners, D18 sources and all verification evidence are
unchanged. The A18 scalar module supplied an implementation style; it is
neither imported nor rerun.

Root authored the new forest theorem and requested the precise two-local-
bridge chain interface. I authored its inherited finite-centre classification
and derived this Lean implementation, including the strong induction. The
present author source check is not an independent review of the new
implementation or of my own input.

Before extending the scalar preparation, BOTH frozen uncompiled predecessors
were byte-preserved and gzip-roundtrip verified. The small preservation
manifest is
`results/source_archive/heavy_centre_forest_before_local_chain_induction_20261023/manifest.json`,
SHA256 `ff4b711b610797450027d14cca79378325e160ef1b3122f8d396ebd64b311a3a`.
It binds source828a17b9 and map e61b5f5c with their exact archived bytes.
Neither predecessor received compilation or axiom-output credit.

## Exact written bindings

| File | SHA256 | Restricted role |
|---|---|---|
| results/finite57_heavy_centre_forest_and_linear_exit_expansion_candidate_20261023.md | 059e3fcf6c80700a0b11778638cbb47e758ca245f6cd9573c747e1ff8caa158f | Complete written actual-edge, metric forest and incidence theorem |
| results/review_finite57_heavy_centre_forest_linear_exit_expansion_full_scope_compact_reviewer_20261023.md | 2510181660697bbc329f3c34762891d46f8fcc96f4a089aa818f96b809358b06 | Full different-author mathematical source review |
| results/review_finite57_heavy_centre_forest_and_linear_exit_expansion_root_acceptance_20261023.md | 97ff5c89037fc7e70e6b1393568750e52c95d50a6900d70b136c287ef43cefb5 | Root acceptance of that written theorem, without capacity or dimension56 |
| manuscripts/finite_volume57.tex | c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0 | Actual all-wall, finite-simplicity and incidence setting |
| results/finite57_complete_finite_centre_recipients_and_heavy_component_escape_candidate_20261023.md | 994d3b9ed36a32e8698610140ee656b40a88507c7f0a7b6b1978faef21b00186 | Complete finite-centre classification, full retained support and reverse ports |
| results/review_finite57_finite_centre_safe_union_and_heavy_escape_root_acceptance_20261023.md | 1c2e887918d4277efffd9c568aed5d9e203a61dce2ca7c98eca96665831792f9 | Accepted scope of that inherited classification |

All six bindings were freshly metadata-checked. Hashes authenticate the
sources; they do not turn their geometric conclusions into Lean instances.
The historical candidate/review labels in immutable headers are retained.

## Primitive scalar data and derived bounds

All metric quantities range over one generic linearly ordered commutative
ring `R`, using the existing `Lean.Grind.OrderedRing` interface. There is
no real square-root, hyperbolic function or polytope object in this module.
Every theorem uses explicit parameters, avoiding body-only section-variable
capture. The primitive long-edge data are exactly

    16 <= c, 0 <= s, s^2=c^2-1, 0 <= t, 2*t^2=1.

In the intended application these are `c=cosh L`, `s=sinh L` and
`t=1/sqrt(2)`. Existence and interpretation remain external. The kernel
derives, rather than assumes, all of

    c^2>=256, 4*t<3, 16*s>15*c, s<c,
    c-s*t>4, s-c*t>3.

The first four bounds use differences of squares. In particular assuming
`16*s<=15*c` gives `256-31*c^2>=0`, contradicting `c^2>=256`.
The two gap inequalities then follow from the positive product
`c*(3-4*t)` and the nonnegative product `(c-s)*t`. No supplied chord,
numerical radical approximation or final cone gap is a premise.

| Declaration in namespace `CoxeterProofs.HeavyCentreForestScalarKernel` | Internally derived conclusion |
|---|---|
| `cone_primitive_bounds` | The four bounds immediately above, from the primitive square data |
| `long_edge_cone_gaps` | Both strict cone gaps, from the same primitive data |
| `cosine_rule_growth` | `cD>4*cR` from `cR>0`, `0<=sR<=cR`, `cosB<=t` and the exact cleared cosine rule |
| `cotangent_rule_growth` | `cotC>3` from `cothR>=1`, `cosB<=t`, `0<sinB<=1` and the cleared cotangent rule |
| `triangle_cone_step` | Both growth conclusions together with the SAME long-edge and angular data |
| `finite_chain_cone_invariant` | Every prefix `1<=k<=n` satisfies its cone and `path k != path 0`, derived by strong Nat induction from two LOCAL bridges |
| `finite_chain_no_return` | `path n != path 0` for `n>=1`, by the preceding invariant; no global no-return or no-cycle premise |
| `projected_centre_length_lower` | `coshL>=rank+2>=16` from the primitive projected-normal equations |
| `forest_boundary_identity` | `boundary=H+2*components` from degree-three incidence and forest Euler identities |
| `connected_component_exits` | `boundary=H+2` when the supplied Euler identity has one component |
| `finite_receiver_count` | `H+2*components<=3*receivers` given the ACTUAL three-port inequality |

The cosine-rule proof explicitly permits negative `cosB`. It derives

    cR*t-sR*cosB=(cR-sR)*t+sR*(t-cosB)>=0,
    cD-4*cR=cR*(c-s*t-4)+s*(cR*t-sR*cosB)>0.

Thus it does not multiply an inequality by a potentially negative cosine.
The cotangent rule is supplied only in the denominator-cleared form

    sinB*cotC=cothR*s-c*cosB.

Its numerator exceeds3 by the derived gap and two nonnegative products.
If `cotC<=3`, positivity and the upper bound on `sinB` force the numerator
at most3, a contradiction. There is no inverse operation or division.

For centre normalization the primitive data are

    rank>=14, delta>0, (rank+1)*delta=1,
    shadow>=1+delta, wall<=0, delta*coshL=shadow-wall.

The actual ranks are14/20. The kernel treats `rank` as an ordered-ring
quantity and proves the stronger continuous statement, without casting
heights or lengths to naturals. The key identity is

    delta*(rank+2-coshL)=1+delta-shadow+wall<=0.

## Conditional finite-chain induction

Let `V` be an arbitrary type, `path:Nat->V`, `n>=1`, and
`cone:Nat->Prop`. Primitive base facts are `cone 1` and
`path 1 != path 0`. The two explicit LOCAL external geometric bridges,
each restricted to `1<=i<n`, have exactly the types

    cone i -> path (i+1) != path 0,
    cone i -> path (i+1) != path 0 -> cone (i+1).

The first comes from tangent-angle control. The second comes from the
nondegenerate hyperbolic triangle and cotangent bound AFTER the new point
is known distinct from the start. They are not a global no-return premise.
Strong induction on `k` obtains the previous cone at `k-1`, then uses the
first bridge to exclude the next return and the second to propagate the
cone. It proves `cone k` and `path k != path 0` for EVERY `1<=k<=n`.
The second declaration specializes to `k=n`. Thus the conditional global
no-return logic is internal; real geometry realizing the two local bridges
and converting a graph cycle into such a path remains external.

No distance is represented by a natural number. These Nat arguments are
only finite path indices. Quantitative iteration of the strict factor4
distance growth, beyond the one-step scalar theorem, remains in the written
proof rather than an extra formal claim.

## Geometric and counting boundary

The written theorem, not this implementation, must supply all of:

* the ORIGINAL finite-volume Coxeter57 polytope, all heavy vertices of
  whole type `2D17+D23`, the three actual centre edges, and their actual
  finite/ideal endpoints;
* the complete heavy-completion classification and positive full-white
  inverse-Gram calculation, including the projected norms, white shadow,
  actual wall product and normalized `coshL` equation;
* orthogonality of distinct centre directions, tangent-sphere angle
  control and the existence/nondegeneracy of each hyperbolic triangle;
* the cosine and cotangent rules, signs of the hyperbolic/sine quantities,
  and the interpretation `cotC>3` as the strict backward-angle cone;
* the interpretation of `cone`, the two LOCAL geometric bridges just
  listed, and realization of any graph cycle as an eligible finite path;
* finite graph incidence, the forest Euler identity, and the actual
  reverse-port bound when all relevant targets are finite.

The Nat variables `H,E,components,boundary,receivers` are finite counts
only. The three incidence declarations prove linear arithmetic consequences;
they do not create graphs, realize the local bridges or identify counts.
There is no supplied global `noCycle` premise in the chain declarations.
The receiver inequality is the denominator-cleared count corresponding to
a ceiling bound, not a receiver-capacity, debt or matching assertion.

For a subset of heavy vertices, boundary may include edges to heavy
vertices outside that subset. Such incidences supply no automatic fresh
capacity. No ideal reverse-port bound, full8 financing, Hall condition,
universal364 cap, dimension56, example or global ceiling follows here.
The original C0 rule and every existing allocation remain unchanged.

## Requested verification

The footer requests exactly the eleven fully qualified declarations in the
table. There are presently zero actual outputs and no successful object
or compiler binding. Root alone may admit one module compile and a
SEPARATE axiom audit after full static review. Failed/refused attempts
would receive no proof credit. Even successful checks would certify only
the conditional ordered-ring algebra, local-to-global no-return induction
and incidence algebra. Actual geometric instances and graph realization
remain at the explicit written boundary above.

# Compact17: complete recovery from the retained local energy proofs

This formulation reuses the accepted local inequalities and complete local
contradiction proofs in
`manuscripts/revisions/compact_farthest_vertex_bounds_consolidated.tex`
(SHA256 `6222f5c6cee8b30b2c1019d648bc31cde8bfcddf08f031a0355a3e7137c107ac`).
It does not attempt to rederive them from a new Gram or SOS relaxation.
The accompanying pure supplier is
`scripts/compact17_proof_recovery_energy_catalogue_20261118.py`.
Neither its generator nor an LP has been run during this preparation.

The output is a complete finite **family catalogue**, with the infinite
ordinary dihedral parameter kept as a uniformly covered symbolic domain.
It is ready for root's bounded materialization and exact LP check. The
source statements remain mathematical inputs; their labels do not replace
their proofs. No final theorem excluding a dimension is an input row.

## One energy account and a small complete cover

Choose the original farthest vertex and original vector as in the retained
manuscript. The orthogonal incident-component energies are nonnegative and
sum to `||x||²<1`. The external published compact ceiling29 reduces the
remaining actual dimensions to18–29; no face inheritance or rank decrement
is used. The compact problem has no arithmetic field restriction.

Call `A_r,r>=7` and `B_r,D_r,r>=9` dominant. The following is a direct
uniform cover of **all** spherical forests, before using any historical
pruned frontier:

1. No dominant: one aggregate energy slot of rank `n`, with lower bound
   `(n−1)/14`. Every ordinary nondominant rank costs at least rank/14;
   the maximum permitted defect debit is1/14.
2. At least two dominants: two actual dominant slots, each with floor1/2,
   plus the orthogonal remainder. Their ranks need not be fixed; their
   energies partition the original budget and the remainder may be empty.
3. Exactly one dominant of rank `r`, with spectator rank `S=n−r>=8`:
   retain a dominant slot with floor1/2 and one aggregate spectator slot
   with floor `(S−1)/14`. These rows already imply total energy≥1.
4. Exactly one dominant and `0<=S<=7`: generate every nondominant spherical
   partition of `S`, every allowed defect profile, and all applicable rows
   below, unless the uniform component cap already gives an explicit exit
   or its strongest applicable spectator-uniform floor `d` satisfies
   `d+(S−1)/14>=1`. In that latter case retain one complete aggregate branch
   with those two rows and a unit-coefficient witness. This is proof-backed
   coarsening of the whole context, not unrecorded energy pruning.

Thus the catalogue does not enumerate all coloured forests up to29. It
keeps the first three families as explicit coarse LP branches, with their
unit-coefficient dual witnesses, rather than silently deleting them.
Exceptional types and components of rank8 in the first three families are
covered by the same complete ordinary theorem. They cannot fit in a
detailed spectator partition of rank≤7.

For the detailed recursion use `A1`, `A2=I2(3)`, the full `I2(m),m>=4`
bucket, `A3..A6`, `B3..B7`, `D4..D7`, `E6,E7,F4,H3,H4`.
No `B2` or `D3` alias is separately counted. The unique dominant is an
actual `A`, `B` or `D` at its own rank. Component permutations do not change
any applicable row. An actual defect on equal-type slots may be relabelled
to the first slot (first two for `J`); each ordinary dihedral slot still has
its independent full integer parameter range. An `O7` mark fixes only its
owned old dihedral to7. This covers every actual marked forest without
asserting that any formal forest is realizable.

The permitted marks are none, `I`, `J`, `N7`, `O7`, at most one total:

| Mark | Actual old support and retained meaning | Replacement floor |
|---|---|---:|
| I | one old A1, actual isolated dotted short opposite | 0 |
| J | two old A1 factors in the actual common short-wall profile | 17/50 jointly |
| N7 | old A2 with its actual isolated native7 opposite | 9/100 |
| O7 | old I2(7) with its actual isolated native3 opposite | 9/100 |

Marked supports do not also receive ordinary floors. The source proves
compatibility and incompatibility of these marks (lines380–410 and
`app:ordinary`), rather than treating their names as sufficient geometry.
The generic spectator debit is0 for none/J,1/14 for I, and37/700 for either
typed defect; J occupies two ranks. The coarse maximum debit1/14 is valid
for every profile. All finer rows retain the actual profile identity.

## Complete scalar and joint rows used by the supplier

Every entry below is a strict proved lower bound, safely weakened to≥ in
the master. The contradiction still uses the actual strict total `<1`.
All decimal strings denote their exact terminating rational values.

| Ordinary irreducible | Energy floor |
|---|---:|
| A1 | 1/10 |
| A2 | 319/2000 |
| I2(m), every integer m≥4 | 1/7 |
| A3 | 1/4 |
| B3,H3 | 2/7 |
| A4 | 2/7 |
| B4,D4,F4 | 1/3 |
| H4 | 11/36 |
| A5 | 3/8 |
| A6,B5,B6,D6 | 3/7 |
| D5 | 5/14 |
| B7,D7 | 1/2 |
| B8,D8 | 4/7 |
| E6,E7,E8 | 9/20,5/9,2/3 respectively |

These are the source's complete common table (lines362–378), proved in
`app:ordinary`. The B8/D8/E8 entries participate in the coarse proof; rank8
does not appear among detailed spectators.

For A-components retain all applicable thresholds `r>=7:1/2`,
`r>=8:4/7`, `r>=13:5/7`, and the exact-rank refinements:

| Actual rank | Stronger floor | Complete local source |
|---:|---:|---|
| 12 | 13/20 | `app:A12`,4150 |
| 15 | 413/500 | `app:A15`,3939 |
| 16 | 43/50 | `app:prop:A16`,3466 |
| 17 | 6/7 | `app:Achain` |
| 18 | 91/100 | `app:Achain`, common table430 |
| 19 | 563/625 | `app:tail:lem:A19ordinary`,16801 |
| 20 | .92274 | `app:tail:lem:A20A1`,15995 |

These permit arbitrary finite spectators. The old A16 floor199/250 is
dominated by43/50 and need not be separately attached. No arithmetic-only
improvement is imported.

Every actual B/D component of rank≥9 has floor9/14
(`app:prop:BDnine`,2208). For each `M=7..16` the genuine generic-prefix
hypotheses are `B_r:r>=M+5` and `D_r:r>=M+6`. Its exact conservative floor
`P_M`, and the stronger replacement `P_M+3969/800000`, are:

| M | P_M |
|---:|---:|
| 7 | .7349975 |
| 8 | .7683275 |
| 9 | .797425 |
| 10 | .823735 |
| 11 | .8472225 |
| 12 | .8689525 |
| 13 | .88864 |
| 14 | .90715 |
| 15 | .9240775 |
| 16 | .9401975 |

The formula and rank hypotheses are `prop:prefix`,453–504; rational
evaluations are displayed at551–569,4730–4745,15160–15180. The extra gain
is `app:tail:lem:prefixgain`,15774–15811. It replaces a floor on the same
support; it is never charged as an extra orthogonal energy. All applicable
prefix rows may coexist. The supplier does not extrapolate the generic
rank domain or use a numerical approximation to an irrational expression.

Three separately supplied actual-rank continuations add component floors:
entire B16/D16+ordinary T2 gives .8689525; entire B17+ordinary A1 gives
.90715 (`app:exact18layers`,4769–4977). Entire D19+ordinary A1 gives
.91211125 after the actual D4 continuation through15
(`app:tail:lem:D4proper`,16026 onward, the prefix-gain calculation, and the
explicit exact20 consumption at17930). This last item is a local
continuation/cost, not the exact20 exclusion proposition as a premise.

The complete same-vector joint rows are I+A1≥1/5, I+two ordinary A1≥3/10,
I+ordinary T2≥21/100, and I+ordinary T3≥2/7 (`app:Ijoint`,2034–2114).
Retain a row for every matching support. Beside N7 or O7 every *other*
ordinary rank-two component alone costs≥1/5. Every other A1 alone costs
≥1/5 beside N7 and≥2/15 beside O7 (`app:typedpropagation`,2117–2196).
The defect's own energy is not included a second time in these latter rows.

## Complete local contradiction leaves

The following are named local theorems with full native, continuous,
foreign-contact and original-wall hypotheses. They are retained wholesale;
the master does not repeat their geometric proof or invent a scalar floor
for a null/native-angle/acuteness exit. An exact matching branch receives
`exit={theorem,hypotheses}`, meaning its mass is0. Component cap exits allow
arbitrary spectators; every other use below is restricted to the indicated
**entire** forest even if its source proves a stronger spectator statement.

| Actual forest family | Local source |
|---|---|
| any A21..29, B23..29, D24..29 component | `app:tail:prop:caps`,14306 |
| B_r+I,r≥20; D_r+I,r≥21 | `app:tail:prop:classicalI`,15189, no extra spectators |
| A16+N7/O7 | `app:A16typed18`,11583 |
| A17+ordinary A1 or I | `app:A17rankone18`,11941 |
| B16/D16+N7/O7 | `app:BD16typed18`,12320 |
| B17+I | `app:B17I18`,12925 |
| D17+ordinary A1; D17+I | `app:D17ordinary18`,12735; `app:D17I18`,13814 |
| whole A18 | `app:prop:wholeA18`,5895, complete proof4980–5916 |
| whole B18 | `sec:b18`,638, complete `app:wholeB18`,5918 onward |
| whole D18 | `sec:d18-setup`,662; `app:D18foundationcomposition`,11533; `sec:composition`,1255–1352 |
| A17+N7/O7; A18+I | `app:tail:lem:A17typed`,17951; `app:tail:lem:A18I`,18054 |
| D17+ordinary T2; B17+N7/O7; B18/D18+ordinary A1 | `app:tail:lem:proper19standard`,18691 |
| D17+N7/O7; B18+I; D18+I | `app:tail:lem:D17typed`,18265; `...:B18I`,18188; `...:D18I`,18452 |
| whole A19/B19/D19 | `app:exact19:prop:wholeA19`,19733; `...:wholeB19`,21873; `...:wholeD19`,25061 |
| A19+I | `app:tail:lem:A19I`,16885 |
| B18/D18+N7/O7; B19/D19+I | `app:tail:lem:proper20classical`,16687 |
| whole A20/B20/D20 | `app:tail:lem:wholeA20`,17042; `...:wholeB20`,17207; `...:wholeD20`,17527 |
| A19+N7/O7; A20+I | `app:tail:lem:A19typed`,15970; `...:A20A1`,15995 |
| D19+N7/O7; D20+I | `app:tail:lem:proper21`,16082 |
| whole B21/D21 | `app:tail:prop:wholeB21`,16194; `...:wholeD21`,16270 |
| whole B22/D22/D23 | `app:tail:prop:whole58`,15556 |

In particular the whole-D18 leaf includes *all* head and fork branches,
not just the last near/far calculation. The actual adjacent-vertex local
theorems are already included within that complete local contradiction;
their new vector's energy is never summed with the original vector's energy.

## Master LP and exact certificate

The supplied schema is exactly
`{family,dimensions,sources,premises,coverage,branches}` plus a scope status.
Each branch has `key,dimension,domain,slots,rows,exit`; rows have exact
rational `lower`, their support slot keys and theorem id. Aggregate slots
are explicitly marked; unknown constituent ranks are `null`, not fictitious
rank0. Detailed component partitions have integer ranks summing to `n`.
Coarse branches include their evident unit-coefficient `hand_witness`.

For each nonexit branch b introduce mass x_b≥0 and energies e_bj≥0,
and impose Σ_b x_b=1 and every row

    Σ_(j in support_l) e_bj >= lower_l * x_b.

Minimize Σ_bj e_bj. Exit branches have x_b=0. Any actual polytope maps to
one branch with mass1 and its original energies, giving objective<1.
Overlapping joint cuts remain simultaneous constraints on the same energies.
Alternative floors on a single component give their maximum, not their sum.

Root's smaller dual-feasibility form uses a free common η with η≥1 and
nonnegative α_bl for every nonexit branch row:

    Σ_(l covering j) α_bl <= 1       for each branch slot j,
    Σ_l lower_l * α_bl - η >= 0     for every nonexit branch b.

An exact feasible rational point proves the original total energy≥1 in
every nonexit branch. Every multiplication and capacity inequality must be
checked against the complete generated catalogue. A solver status is not
the certificate. A normalized slack-cone primal alternative is given in
`existing_proof_energy_lp_master_20261118.md`; it handles the strict budget
without choosing an arbitrary fixed epsilon.

The logical split is precise. Original-vector geometry, classification,
defect compatibility, all named local provider proofs and the external29
ceiling are retained mathematical inputs. Hierarchical coverage, all row
applicability matches, exact finite catalogue construction, and the final
capacity inequalities are the new finite LP recovery work. The catalogue
does not use `thm:main`, `thm:tail`, or the exact-rank exclusion propositions
as premises. Their accounting paragraphs are only source cross-checks of
the displayed local floor evaluations and terminal-family inventory.

## Source-derived construction bounds, without running the generator

All generated theorem ids and branch ids have the `compact.` prefix, so
arithmetic local theorems cannot collide with them in the combined catalogue.
The single-dominant `(n,family,r)` loop has582 contexts: over n18–29 the
A counts sum to210 and the B/D counts to186 each. Add24 no-/multi-dominant
coarse branches. Every context initially accounts for one branch. Component
caps and sufficient aggregate costs remain one branch each.

Only these contexts can enter a nontrivial spectator partition recursion:

| S | Possible dominant types/ranks | Number of contexts |
|---:|---|---:|
| 4 | A14,D14 | 2 |
| 3 | A15; B15,B16; D15,D16,D17 | 6 |
| 2 | A16..A20; B16..B19; D16..D20 | 14 |
| 1 | A17..A20; B17..B22; D17..D23 | 17 |
| 0 | any remaining whole component | one branch per context |

There is no detailed S≥5 context. At S5 the weakest A bound is5/7 and
the weaker D13 prefix already exceeds5/7. At S6 the A12 bound13/20 and
the B/D floor9/14 suffice. At S7 the A11 bound4/7 and the B/D floor9/14
suffice. The S≥8 branch uses the half bound directly. These are rational
comparisons of displayed theorem constants, not a generated-case count.

The ordinary spherical partition coefficients through rank4 are
`F0..F4=1,1,3,6,14`. They follow from the formal truncated product

    (1-t)^(-1) (1-t²)^(-2) (1-t³)^(-3) (1-t⁴)^(-5) modulo t⁵.

Marking one I removes one A1; marking J removes two; marking N7 or O7
removes its one specified rank-two factor. Thus the complete marked count
is `G_s=F_s+F_(s−1)+3F_(s−2)`, giving `1,2,7,12,29` for s0..4.
These are elementary source-derived bounds; no catalogue was enumerated.
Consequently the total number of branches is at most

    582+24 + 2(29-1)+6(12-1)+14(7-1)+17(2-1) = 829.

Every branch has at most5 slots. A detailed branch has at most21 dominant
rows (the BD floor plus ten prefix/replacement pairs). With S≤4, all
ordinary/defect/joint/propagation rows together require at most9 further
rows; at most one exceptional exact-rank continuation adds a row. Each row
has at most3 support slots. Coarse branches have at most3 slots and2 rows.
Safe complete catalogue bounds are therefore:

    branches <= 829
    sum slot counts <= 4,145
    sum cost-row counts <= 25,699
    sum row-support sizes <= 77,097.

These include exit branches and deliberately overcount incompatible row
maxima. They are below root's proposed preallocation guards. The dual has
at most25,700 variables including eta, at most4,975 constraints including
eta>=1, and at most103,626 matrix entries when both cost/eta coefficients
are counted. No scientific allocation, resource-fit measurement or solve
has occurred. Root must still enforce the actual joined arithmetic+compact
catalogue limits before allocating numerical arrays.

No executed enumeration count, numerical solve, exact candidate, or recovered global
bound is credited here. The next concrete step is bounded materialization,
source-bound catalogue review and one exact certificate check by root.

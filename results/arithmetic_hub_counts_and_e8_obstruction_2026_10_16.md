# Arithmetic hub counts and the E8 host identification

**Internally reviewed computer-assisted proof candidate: every finite-volume
arithmetic hyperbolic Coxeter polytope has dimension at most 21.** The bound
would be sharp, since Borcherds' established example attains 21. Independent
weight checks, complete face coverage, and the all-field composition pass
generation and fresh saved replay. External mathematical validation and
literature priority remain open. No new finite-volume example has been
constructed.

## Two successive weight improvements

The [preceding checkpoint](arithmetic_dimension23_and_complete_face_coverage_2026_10_15.md)
gave an internally reviewed arithmetic bound-23 candidate. The next rule
distinguishes zero, one and at least two entire I2(5) components, together
with the exact number of A1 components. Its
[independent review](fresh_sqrt5_hub_count_weights_20261016_review.md) checks
all 54,993 face cases and 144 exact vertex contexts. Its
[complete coverage](fresh_sqrt5_hub_count_coverage_20261016.md) and
[all-field composition](fresh_arithmetic_hub_dimension_composition_20261016.md)
give the completed intermediate arithmetic bound-22 candidate.

Dimension 22 requires a genuinely additional pair-weight family. The
[new H3 rule](fresh_sqrt5_dimension22_h3_cross_20261016.md) lets an entire H3
component pay across to a root in another component. If the vertex contains
`a` A1, `b` I2(5), and `c` H3 components, the three cross payments divide by
the actual source-root counts `a`, `2b`, and `3c`. Each later family excludes
the earlier source types, so no pair is counted twice. Coefficients also
depend on exact `a`, capped `b`, and whether `c` is positive.

| Dimension excluded | Exact vertex maximum | Strict lower threshold |
|---:|---:|---:|
| 22 | 21977/500 = 43.954 | 44 |
| 23 | 485211/10000 = 48.5211 | 49.5 |
| 24, preceding checkpoint | 267527/5000 = 53.5054 | 54 |

At fixed exact `(a,b,c)`, summing over all source roots cancels each
denominator, leaving additive component costs. An independent type-first
dynamic program checks all 410 exact vertex contexts and 71 finite component
types in dimension 22. The independent face check covers 2,151 rows and
89,255 common spectator triples, including all 1,113 rank-feasible stops
from the complete list of 2,070 finite-queue singular obligations. Its
minimum face slack is 3/560000. Floating-point optimization supplied only
the proposal; all displayed bounds and inequalities are exact rational
checks.

The fixed triangle and short-path/H4 payments are preserved. The
[complete final coverage](fresh_sqrt5_h3_hub_coverage_20261016.md) directly
reevaluates every older, dotted and complete finite-queue branch under the
new rule. It retains all 2,661 stopping rows, checks all 1,634 rows fitting
dimension 22 in 91,726 common spectator cases, and verifies 471 literal maps.
The binder separately reconstructs the complete external-triple domain and
requires each exact four-corner sum to be at least one. The strict h-vector
lower threshold is `5·55−binom(22,2)=44`. Both the
[independent weight review](fresh_sqrt5_dimension22_h3_cross_20261016_review.md)
and the final coverage pass generation and saved replay. A
[separate standard-library scope checker](fresh_arithmetic_dimension21_scope_check_20261016.md)
then verifies all 91,726 saved payment sums, every external-triple domain,
all 2,661 originating row identities, the 410 vertex contexts, the strict
threshold, and all 24 field/dimension branches. Its generation and fresh
replay also pass.

## Scope and the limit of the arithmetic programme

Vinberg's published field restriction applies from dimension 22 upward.
The [completed all-field composition](fresh_arithmetic_dimension21_composition_20261016.md)
combines the Q(√5) bound 21 with the preceding Q(√2) and cubic-seven bound-20
candidates and excludes all 24 field/dimension branches for dimensions 22–29. Noncompact arithmetic reflection groups already satisfy
the published bound 21, attained by Borcherds. Thus the arithmetic bound-21 candidate is sharp across compact and
noncompact cases together. It does not establish that a compact example in
dimension 21 exists or exclude one. Every arithmetic reflection lattice is of simplest
type; arbitrary arithmetic lattices are outside this statement.

The unrestricted compact ceiling remains the published 29. The separate
unrestricted finite-volume candidate 163 is unchanged. External mathematical
validation and literature priority remain open for all repository proof
candidates. A [fresh scope audit](fresh_arithmetic_dimension23_scope_audit_20261016.md)
found no gap in the earlier geometric, arithmetic, attachment, spectator
and h-vector transitions. It also checked 112,764 saved payment sums and
their spectator quantifiers independently. It explicitly records which
large enumerations are inherited and corrects a frozen prose overstatement
about the primary fast queue's saved replay.

The audit's bounded literature check found that a 2026 author survey still
reports the published compact/arithmetic exclusion from dimension 30 and
no known examples above dimension 21. This corroborates the comparison,
without establishing an exhaustive priority result.
[Mühlherr–Weiss, Section 7](https://londmathsoc.onlinelibrary.wiley.com/doi/pdf/10.1112/jlms.70379).

## Constructive result: a duplicate host, now identified exactly

The [height-100 E8 continuation](fresh_silver_e8_h9_continuation_20261016.md)
has 48 walls, 100 finite vertices and 182 finite edges. An
[independent Fraction-only calculation](fresh_silver_e8_full_chamber_20261016_review.md)
repeats all 29 controls and the complete three-norm prefix using a different
covector recursion. It checks the three complete finite corners and two
integral sector maps needed for a full-chamber stabilizer with trace
`10+4√2>10`. The full reflection chamber therefore has infinite volume.
The separate [conceptual audit](fresh_silver_e8_scope_audit_20261016.md)
checks the persistence, local-finiteness and chamber-stabilization argument.

There is an important correction to the host search: the literal lattice
`〈−1−√2〉+E8+〈1〉` is integrally isometric to the previously excluded
`〈−1−√2〉+F4+F4+〈1〉` host. The
[explicit identification](fresh_silver_e8_f4_identification_20261016.md)
and [independent matrix check](fresh_silver_e8_f4_identification_20261016_review.md)
give an integral future change of basis and integral inverse with determinant
one. It mixes time and spatial coordinates. Different positive blocks did
not make different Lorentz lattice classes. The new symmetry proof remains
valid but excludes no new integral isometry class. Future hosts should be
screened for these identifications before long root searches.

## Diagnostics and preservation

Several smaller weight models have exact dual obstructions: the dimension-22
[two-count family](fresh_sqrt5_dimension22_i25_counts_20261016.md), including
exact I2(5) counts, the [E8 context](fresh_sqrt5_dimension22_e8_context_20261016.md),
and the [A2 context](fresh_sqrt5_dimension22_a2_context_20261016.md) and
[A2 cross family](fresh_sqrt5_dimension22_a2_cross_20261016.md).
These are restricted-model obstructions on conservative face catalogues,
not existence results or barriers to all weight rules. The H3 family is
the successful enlargement.

The [archive checker](../scripts/audit_hub_counts_e8_checkpoint_20261016.py)
and [integrity certificate](arithmetic_hub_counts_and_e8_obstruction_2026_10_16_audit.json)
bind the new artifacts while preserving the previous frozen source and
report hashes. Archive integrity and mathematical validity remain distinct.
The original goal of finding a new finite-volume example or construction
dimension remains open.

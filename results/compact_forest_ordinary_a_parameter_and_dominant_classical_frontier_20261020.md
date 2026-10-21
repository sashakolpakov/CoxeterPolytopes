# Ordinary A ranks and the unique large B/D forest frontier

20 October 2026. Complete paper aggregation, awaiting root and outside
review. No scientific execution, census, optimizer, Lean, TeX or
accepted-source edit. The new B5--B8 rates are used as accepted inputs;
the exception accounting is checked here independently of the type
of the dominant component.

## Statement, with the A premise explicit

Fix an integer M>=6. Assume that every WHOLE A_j component with
5<=j<=M at an actual farthest vertex pays energy strictly greater
than j/14, with arbitrary finite spectators. Let p be a vertex
farthest from a strict interior point of a compact hyperbolic Coxeter
n-polytope, n>=15. If its finite diagram has NO A component of rank
greater than M, then it has exactly one B_r or D_r component of
rank at least9, and

    r>=max(10,n-5),    spectator total rank S=n-r<=5.       (1)

The other components and their contacts are not assumed ordinary
without exception: the proof keeps the possible one short A1 or one
isolated rank-two7 defect explicitly.

The premise is presently accepted for M6. No A7/A8 ordinary rate is
asserted here; an independently proved extension of that premise can
be substituted without changing this aggregation. In dimension29,
the accepted large-component exclusions then imply that EVERY
farthest vertex must contain an A_j with

    M+1<=j<=20.                                         (2)

In particular, the current unconditional-on-components case M6
leaves a necessary A7--A20 component. This is not an unconditional
improvement of the compact or arithmetic dimension ceiling.

## 1. Pinned rates and their arbitrary-spectator scope

All cited files are under results/. At a farthest vertex the common
vector satisfies ||x||^2<1, and its orthogonal whole-component
projections have additive squared norms. This is the convention in
`compact_farthest_vertex_energy_obstruction_20261020.md`, SHA256
`50c64e5a521d3e4e0f277a8ece02d3a29aa4ae8b71485cf3088211e587f45bff`.

The complete rank-at-most-four list and its exception interfaces are
in `compact_farthest_rank_four_forest_exclusion_20261020.md`, SHA256
`68ad62a619603eb7e6c2acc5689c901dde3f9d0897a7c1580b3deb948f5832b6`,
using `compact_universal_small_a_component_energy_20261020.md`, SHA256
`801581c4b855289aca7982b2d0f4bb6f5849f151a78f40b86c9d403e84f415e6`.
Their component conclusions, not just their whole-small-forest
corollaries, hold with arbitrary finite spectators. They give:

* nonshort A1 directions pay more than1/10>1/14;
* every other rank<=4 component pays more than rank/14, except the
  two isolated rank-two defects recalled in Section2.

The remaining accepted ordinary components are supplied by:

| Components and sufficient rate | Source | SHA256 |
| --- | --- | --- |
| A5>3/8>5/14 | `compact_universal_a5_actual_projection_three_eighths_20261020.md` | `fd6f024d975e525184de05961d4d6943f1b0bc8eef4acc279dde98cbc1bc4456` |
| A6>3/7 | `compact_universal_a6_actual_projection_three_sevenths_20261020.md` | `022c0cc77512576dfc7bed74be372a94956618498c9ed74508aff53c40a50ae2` |
| B5/B6>3/7, B7>1/2, B8>4/7 | `compact_b5_b8_ordinary_component_energy_20261020.md` | `079f6fca7185f001eedea069505b8276d9bf7a8ff0beb27cd6f57327234dfe56` |
| D5>5/14, D6>3/7; E6>9/20, E7>5/9, E8>2/3 | `compact_universal_e6_e7_e8_component_energy_20261020.md` | `acf261d97f6eb1a11a64502e17e079ef7f418b9746fbf85dbdf6915750fe9a18` |
| D7>1/2 | `compact_dn_center_neighbor_common_opposite_energy_20261020.md` | `e371611c4903ad362654ee0d4656457140cdc5e4374b4e19a2697718c2d679b1` |
| D8>4/7 | `compact_dn_uniform_energy_and_chain_relaxation_20261020.md` | `1bd7c4eb372083973ca059b34121cdf0e3f575945a1552d771d7b09c84617804` |

Every rate in this table strictly exceeds or equals rank/14, with a
strict energy inequality in every case. The finite irreducible
classification has no further types of rank>=5 beyond A,B,D,E.
The M premise covers the further A types allowed by the statement.

The TWO possible large-component types both have the uniform actual
projection bound

    energy>9/14 for B_r or D_r, r>=9.                    (3)

For B use `compact_bn_actual_projection_nine_fourteenths_20261020.md`,
SHA256 `08056719c47080f233c42d4b92c7d3537696956d01dd847639373adf9569c8cd`;
for D use `compact_dn_actual_projection_nine_fourteenths_20261020.md`,
SHA256 `99bde3780f9d720d9809f5d8bea6795e4a77fcaacf8ede5871f0f831fb7823e9`.
These are actual projection estimates with arbitrary spectators,
not an extrapolated rank/14 cost or an edge-only relaxation. The
proof below uses no component-specific contact inventory from them.

## 2. Defect compatibility is independent of the large component

The pinned universal small-A Section5 proves that there is at most
ONE short A1 direction, where short means L<asinh(1/sqrt2).
Its actual opposite has only A1 spatial support, of type I or J:

    I: one isolated dotted A1 contact, squared height h in(0,1/2);
    J: dotted A1 plus one A1_3 contact, h in(1/3,3/8).

All other A1 directions individually pay more than1/10. For the
coarse aggregation here the short direction may be discarded
entirely, even for J; its light A1 direction is then nonshort and
still pays its ordinary cost. No borrowing from a B/D component is
used in this step.

The rank-four source Sections2--4 gives exactly two possible
nonordinary rank-two BLOCKS:

    OLD A2 / NEW isolated7: h_A in(0,1/12);
    OLD I2(7) / NEW isolated3: h_D in(0,1/3).

Each such block nevertheless pays energy greater than9/100.
Both named opposite walls have spatial support ONLY in that old
rank-two block. Their isolation statements are universal: they
already retain arbitrary other old finite components and use the
minimum positive foreign cost1/4, not a classification of a large
D component. Every other rank-two block has its full ordinary rate.

Consequently two defect blocks cannot coexist, and neither can
coexist with a short A1 of either type I or J. Indeed their spatial
supports are orthogonal. For the corresponding actual unit walls
the full product is minus the product of their positive heights,

    c.d=-sqrt(h_c*h_d) in(-1/2,0),

because every displayed squared height is below1/2. This contradicts
the actual Coxeter wall gap: two different wall products are0 or at
most-1/2. The isolated walls are distinct since their nonzero spatial
supports lie in different old blocks. Multiplicity of walls within
one defect block does not create a second rank debit.

This reasoning refers only to the small defect supports and their
heights. It does not refer to the dominant component, assume its
walls are disjoint from defect walls, impose a native D cost, or
extend any of the old whole-small-forest propagation estimates.
The only exceptional accounting needed below is therefore:

    either no rank-two defect and at most one short A1;
    or one rank-two defect of cost>.09 and no short A1.      (4)

## 3. Existence, uniqueness, and the five-rank remainder

First suppose that there is no B/D component of rank>=9. By the
finite classification, the M hypothesis and the list in Section1,
all components are ordinary apart from(4). With no rank-two defect,
discard the possible one short A1 and obtain

    ||x||^2>(n-1)/14>=1.

If no short A1 actually occurs the stronger n/14 lower bound holds.
With one rank-two defect the other n-2 ranks are ordinary, giving

    ||x||^2>9/100+(n-2)/14>=713/700>1.

These contradictions prove existence of a large B/D component.
Two such components would pay more than18/14>1 by(3), so it is
unique. Its type may be B or D; denote its rank by r and write
S=n-r. Every spectator block is now ordinary except the named
possibilities(4), since a second large B/D has already been excluded.

The complete spectator accounting is:

| Spectator case | Contradiction threshold | Consequence |
| --- | --- | --- |
| no exception | S>=5 gives energy>9/14+5/14=1 | S<=4 |
| one short A1, no rank-two defect | S>=6 leaves ordinary rank>=5 after discarding that A1 | S<=5 |
| one rank-two defect, no short A1 | S>=6 gives energy>9/14+4/14+9/100=713/700 | S<=5 |

Thus always S<=5. Since n>=15, r=n-S>=n-5>=10. This proves(1).
Every strict lower bound concerns orthogonal projections of the
SAME actual farthest vector. No spectral energy of one formal star
is substituted for that common vector, and no rank is counted twice.

## 4. The dimension29 frontier and the current A obligation

For this final consequence pin the accepted arbitrary-spectator
large-component exclusions:

* A21--A29: `compact_an_large_components_arbitrary_spectator_transfer_20261020.md`,
  SHA256 `8d3b999761be7dd46e7990973b3ec4d4a746d4d682b044d53695f4aae405e3e5`;
* B23--B29: `compact_bn_large_components_arbitrary_spectator_transfer_20261020.md`,
  SHA256 `9852635a5271c1a733999fb9bd49cd4170d019ae3f37f2ae8bb4af246ebc79ed`;
* D24--D29: `compact_dn_large_components_short_tail_transfer_20261020.md`,
  SHA256 `651a9d047d904f1b6496ba2664df5f313878466f67d011dcda144edef6a5620a`.

If n29 and there were no A component of rank>M, Section3 would
force a unique B/D component of rank>=24. The B23 and D24 exclusions
both contradict this. Therefore an A component of rank at least M+1
must occur. The A21 exclusion bounds its rank by20, proving(2).

For the accepted M6 premise, every farthest vertex of a hypothetical
compact29-polytope must consequently contain A7,...,A20. This is
stronger than a necessary large A/B alternative, but the remaining
A components have NOT been excluded. If an ordinary A7 or A8 rate
is proved independently, the same lemma conditionally advances the
lower endpoint to8 or9; no such pending premise is promoted here.
A coupled charge for an A component and some spectators would
require its own debit accounting before substitution.

After importing the existing ambient compact ceiling29, the restricted
class admitting a farthest forest with no A rank>=7 has ceiling28.
The unrestricted compact29 and arithmetic20 ceilings remain unchanged.

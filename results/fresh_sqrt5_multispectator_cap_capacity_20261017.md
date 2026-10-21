# Arithmetic cap capacities with two isolated A1 spectators

The selected two-spectator family gives **555 exact capacity constraints**.
After including every full colored-core orientation, 23 constraints are
violated by the frozen dual of ratio12.25714671503906. This rejects that
particular hypothetical incidence distribution after the new constraints
are imposed. It is not an optimized bound, a feasibility claim, or a
new dimension theorem. No optimizer or witness repair is run here.

The [source](../scripts/fresh_sqrt5_multispectator_cap_capacity_20261017.py)
and [data](fresh_sqrt5_multispectator_cap_capacity_20261017.json) reconstruct
the selected domain directly from frozen actual prism/base provenance.
They do not use the temporary discovery scouts as proof inputs. The
input face masses are the final
[joint-port certificate after 851 cap constraints](fresh_sqrt5_profile_clique_joint_ports_20261017.md).

## Exact geometric domain and normalization

Fix an actual maximal simplex S, its full supporting star, and an
identification of its unique connected Lorentz core with a canonical
weighted core. The full S star consists of this core and an orthogonal
finite spectator forest and has fourteen independent roots in dimension
thirteen. Its normalized roots therefore span the entire ambient space.

The selected actual prism-base incidences have the following literal
property. Deleting the additional wall b from the prism star leaves the
S core and **two internal isolated A1 roots**, q and p. The coface root is
q. The native products are

`b.q = -sqrt(tau), b.p = -1`,

where tau is the squared norm-two product for label5, label10 or a dotted
bond, and tau is never one. All remaining spectator products are zero.
The contact vector with the complete core, together with the subset B
of S vertices forming this particular prism base, is retained as part of
the profile. The complete core contacts in this selected family are right
or label3; no omitted contact is assigned a default value.

Write v for the projection of b to the Lorentz core. With all roots
normalized to square two, orthogonality gives the exact decomposition

`b = v - sqrt(tau) q/2 - p/2`,

and consequently

`v^2 = h = (3-tau)/2`.

No sign assumption on h is made: the finite label5 case has positive h.
The producer independently reconstructs all 555 literal base maps and
all 25 core inverses and verifies this identity. It checks that q,p are
actually isolated after deletion, that b.p has label3, that the complete
core graph maps exactly, and that the saved base vertices are precisely
the other omitted roots in the actual PD two-root deletions. This last
check fixes the convention for `canonical_vertices`.

For cores admitting a physical Gram over Q(sqrt5), v² is computed directly
from its inverse. Otherwise an independently constructed tree gauge is
used: Hii=2s_i, and the gauge contacts are -sqrt(s_i) at the simple-contact
roots. Exact positive square roots of the scale ratios are found in
Q(sqrt5), and the quadratic form is evaluated with those factors. This
includes the two-contact gauge cases; it does not assume every cap has
only one core contact.

All external spectator components of the selected actual star are
orthogonal to b. On fixing an actual S, any other selected occurrence with
the same profile can therefore use only two of the k actual isolated A1
copies in its full spectator forest. It has no unspecified residual
component. Other spectator/contact families are simply not assigned a
new charge.

## Pair products and exact arithmetic exclusions

For a fixed core orientation, represent a cap by the ordered pair q→p
of distinct A1 copies. The roles are distinguished because tau differs
from the simple label3 squared product one. Two distinct caps have the
following products:

| Relation between ordered pairs | Native product |
|---|---:|
| disjoint supports | h |
| common head p | h+1/2 |
| one head equals the other tail, without swapping | h+sqrt(tau)/2 |
| swapped ordered pairs | h+sqrt(tau) |
| same tail q, distinct heads | h+tau/2=3/2 |

The last case is impossible by acuteness. The same ordered pair and the
same complete core contacts determine the same normalized ambient wall,
because the full S star spans. With B also fixed, it is the same base
incidence, not a second occurrence. This prevents duplicate flags from
artificially enlarging a capacity.

For each of eleven exact tau values, the first four products are tested
for two necessary conditions: nonpositivity and algebraic integrality of
the squared normalized product. The latter is a necessary arithmetic
Vinberg two-cycle condition. Passing these tests does **not** assert a
valid Coxeter angle, field compatibility, or a realizable pair; keeping
such possibilities only weakens the resulting upper capacity safely.

The square-root-in-Q(sqrt5) test is complete. If tau has such a square
root, the square of the product is evaluated directly in that field and
its trace/norm integrality is checked against the integral-basis test.
Otherwise, writing the squared product as A+B sqrt(tau), its relative
polynomial is

`X^2 - 2 A X + A^2 - B^2 tau`.

Taking its norm from Q(sqrt5) to Q gives a monic rational polynomial.
Every coefficient is checked exactly, and the polynomial is evaluated
back at the formal radical to verify zero. A coefficient is nonintegral
exactly when the algebraic number fails integrality: the norm polynomial
is a product of all conjugates, with multiplicity if necessary; integral
conjugates have rational integral elementary symmetric functions. When
B=0 the smaller quadratic field polynomial is used instead. This avoids
assuming that every displayed quartic is a minimal polynomial.

## The five directed-pair capacities

Order the four allowed flags as disjoint, common-head, chain, swapped.
Only the following patterns occur. Each bound holds for every k>=2.

| Pattern | Maximum per fixed orientation | Reason |
|---|---:|---|
| 0000 | 1 | No pair of distinct arrows is allowed. |
| 1001 | 2 floor(k/2) | Supports split into disjoint two-vertex components; each permits at most the two opposite arrows. |
| 1000 | floor(k/2) | All arrow supports are pairwise disjoint. |
| 0010 | 3 if k>=3, otherwise1 | Two arrows form a directed path; a third compatible with both must close its directed triangle. No fourth can join. |
| 0100 | k-1 | Every pair has the same head, so all arrows share that head and have distinct tails. |

For the fourth row, same-tail pairs and common-head pairs are forbidden,
so every vertex has in-degree and out-degree at most one. Swapping is also
forbidden. A compatible third arrow is forced to close the three-cycle;
this proves the all-k assertion without a bounded clique census. The
producer additionally exhausts the small k=2,3,4,5 cases as a separate
sanity check, not as the proof for larger k.

## Orientations, incidence multiplicity and signed columns

The complete weighted core automorphism group is reconstructed by exact
backtracking and checked against the inherited colored group. The
spherical-deletion set is recomputed from principal positive-definite
matrices. Every permutation preserves that set, hence the black/white
coloring used for the maximal S. Finite spectator components cannot be
exchanged with the connected Lorentz core.

An orientation is the simultaneous image of the **entire ordered core
contact vector and B**, under this one automorphism group. Fixing any
single identification of an actual S with its canonical star assigns each
actual incidence to exactly one of these orientations. For each
orientation, the directed-pair bound above applies to the k actual A1
copies. Therefore the safe total capacity is

`number of distinct orientations × single-orientation capacity`.

Different orientations are not assumed independent; counting each with
its own upper bound remains valid. In the selected export, 452 columns
have one orientation and 103 have two. The largest resulting capacity is
twelve.

Identical profiles are grouped by base type, simultaneous orientation
orbit and tau. Each actual source base is counted once; repeated face or
provenance entries receive their true integer multiplicity. A column has
positive coefficients on its actual prism-face rows and the negative
coefficient `-capacity` on the corresponding maximal-simplex row. Thus
the required inequality is `sum(actual incidence counts) <= capacity *
number of maximal S`. This has the same signed orientation as the frozen
[shared-base theorem](fresh_prism_simplex_shared_capacity_20261017.md).

## Result and verification scope

The producer examines 1,827 actual A1-coface occurrences and selects all
555 with the stated two-internal-A1 pattern. These give 555 profile
columns, of which 66 have positive load in the final input dual. After
including the orientation factors, **23** have negative slack: ten from
pattern0000, nine from1001, three from1000, and one from0100.

For example, source base2439 has tau=(3+sqrt5)/2, k=3, one orientation,
and capacity one. Its constraint is

`alpha[41621] <= alpha[79450]`.

The saved exact masses are 0.0005982455 and 0.0001994152, respectively,
so the slack is -0.0003988303. Here h>0 already excludes every distinct
same-orientation cap; no subtle quartic test is needed for this example.
The data exports every column, exact load/slack, source base and
provenance, and the required individual simplex-mass lower bounds.

Generation72363 and a fresh saved replay33442 both terminate with PASS.
The checker repeats all raw-map and projection checks, all 25 core inverse,
automorphism and deletion checks, all 44 pair tests, complete grouping and
all exact dual loads. It reuses the frozen independent Fraction arithmetic
helper; it is a new producer, not a second implementation of itself.
The separate [independent review](fresh_sqrt5_multispectator_cap_capacity_20261017_review.md)
passes generation and a fresh saved replay: 44 independent minimal-polynomial
and sign tests, 30 small directed-graph checks, and all 555 saved
occurrence/column capacity identities. Its source/proof read passes the
normalization, orientations and incidence counting; it explicitly does not
duplicate the 555 raw projection reconstructions.
No model repair, new global solve or improved dimension bound is claimed.

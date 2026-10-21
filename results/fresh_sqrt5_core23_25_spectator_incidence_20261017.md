# The raw23/raw25 incidence rule across all represented spectator contexts

The [forced-quadrilateral rule](fresh_sqrt5_core23_25_forced_quad_20261017.md)
extends to all52 matching common-simplex types in the frozen actual-provenance
domain. It yields one further violated inequality after the first
[one-coordinate repair](fresh_sqrt5_core23_25_incidence_repair_20261017.md).
Raising only one additional vertex mass satisfies the whole family and
leaves the exact partial-model obstruction ratio **12.11531601432504 > 12**.
No optimizer is used, and no dimension theorem follows.

The [source](../scripts/fresh_sqrt5_core23_25_spectator_incidence_20261017.py)
and [certificate](fresh_sqrt5_core23_25_spectator_incidence_20261017.json)
save every literal profile map, every family inequality and target forest,
the exact violated row, and the repaired masses.

## The general spectator formula

Let S have the same five-node Lorentz core1 as in the frozen proof and
spectator forest

\[
 F=kA_1+T,\qquad k\ge2,\qquad \operatorname{rank}F=9,
\]

where T has no A1 component. Use the same raw23 profile A and raw25 profile
B, attached to one of the k isolated A1 supporting roots. Their core
contacts and norms do not depend on F. The core has trivial automorphism
group, while spectator permutations give exactly the k interchangeable
coface choices. Other spectator automorphisms fix the complete contacts.

The full S-star is still a nondegenerate14-root basis. The same exact
products prove that there is at most one A cap overall, at most one B
cap per coface, and no mixed A/B pair on the same coface. Therefore

\[
 A_i\in\{0,1\},\quad B_i\le k-A_i,\quad
 X_i=A_iB_i\ge kA_i+B_i-k.
\tag{1}
\]

This counts actual selected prism occurrences. It introduces no factor k
for a diagram automorphism: the k slots already occur in the per-S count.
Each selected prism contributes exactly one base of the specified S type;
the source checks this multiplicity directly in the full base-incidence
records.

For every distinct-slot mixed pair, the identical core projection creates
the compact quadrilateral from the frozen proof. The unused spectators
are orthogonal and positive definite. The minimal five-root quadrilateral
has white spectator forest

\[
 (k-2)A_1+T+H_3,
\]

and its unique native5 corner has full vertex forest

\[
 Y=(k-2)A_1+T+2H_3.
\]

The other three corner forests are obtained by adding respectively
`3 A1`, `A1+A2`, and `3 A1` to those white spectators. Thus the baseline
total remains2. All four corners are actual finite vertices by the same
strict-center projection argument; no further distance parameter appears.

The native5 H3 component again fixes the roles A=center, B=5-leaf, and
w=3-leaf. At the opposite old corner, every originating maximal simplex
contains the same w port. The disjoint maximal-simplex port theorem gives
at most one origin. Hence mixed pairs inject into these pointed actual
quadrilaterals. If h is the number of H3 components of Y, each Y vertex
has exactly h eligible native5 pairs. Other finite components, including
H4 when present, do not create another eligible H3-component corner.
Summing(1) gives the necessary inequality

\[
 k a+b-k s\le h v_Y.
\tag{2}
\]

This is a condition on actual whole-face incidences. Independent local
realizability of each saved row does not supply it automatically.

## Literal coverage and conservative missing profiles

The scan traverses every frozen actual base occurrence whose containing
simplex has core1 and at least two A1 spectators, retaining raw source23
or25. It finds65 profile occurrences across52 S types:

- 13 S types have both selected actual profiles;
- 39 have only a selected A profile;
- 49 have positive source-face mass in the current witness.

Every retained profile has the same literal core isomorphism, omitted root,
coface root, and canonical base-vertex set as its independently frozen
raw23/raw25 geometry input. Each contributes once to its S base type.
The provenance-to-current-face map is injective, and all52 target Y forests
are distinct and already occur in the6558-forest budget universe.

For a missing profile its coefficient is zero. This makes(2) weaker and
remains valid: deleting selected occurrences cannot increase a positive
face-count load. In particular the39 A-only rows reduce to an inequality
already implied by A≤S. The scan does not assert that an unselected profile
cannot occur in an ambient chamber. Its completeness is with respect to
the saved actual-provenance domain, not every possible geometric profile.

Some S rows were absent from the early incidence catalogue's LP and were
added by the later forced-base patch. The current labelled base maps and
row indices are used; the historical `represented_in_lp` flag is not
mistaken for their present status.

## The one additional violation and exact repair

The new violated row has base type55293, S face80185, k=2, and profiles
Q26406/Q26876. Its target is

`Y = A3 + 2 H3 + 2 I2(5)`, vertex index5413.

It requires

\[
 2\alpha_{26406}+\alpha_{26876}-2\alpha_{80185}
 \le 2\mu_{5413}.
\]

The saved left side is `215857/10000000000`; its old right side is
`15573/10000000000`. Thus the exact gap is `50071/2500000000`.
The simultaneous family repair changes only

\[
 \mu_{5413}:\frac{15573}{20000000000}
 \longmapsto\frac{215857}{20000000000},
\]

an increment of `50071/5000000000`. The earlier S580 inequality remains
an equality after its already frozen repair. Every row in(2) is checked
under the final masses.

All face masses are unchanged. Consequently the inherited signed and cap
loads, complete local simplex-packing marginals and their3712 literal
mixture terms, and all17 joint3/4/5-family loads are unchanged. Only a
nonnegative vertex allowance increases. The same monotonicity argument
as in the preceding repair preserves all their checks without rerunning
the2,711,432 joint inequalities.

The objective numerator stays `3842394855727/309375000000`. The exact
normalizer and rational ratio are saved in the certificate and compared
with12 using fractions. The resulting12.11531601432504 is a certified
repaired witness, not an optimization result. Further global compatibility
conditions and complete face coverage remain necessary for a dimension
bound.

**FINAL/FROZEN:** generation and fresh saved replay passed. Direct input
hashes, local report links, and source syntax were checked. This package
inherits the frozen forced-quad proof and its stated internal review scope;
it does not claim a separate arithmetic implementation of the65 mappings.

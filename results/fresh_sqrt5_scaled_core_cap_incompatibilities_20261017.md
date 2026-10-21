# Two scaled-core cap incompatibility families

Two new capacity-one inequalities violate the current exact
[12.11531601432504 witness](fresh_sqrt5_core23_25_spectator_incidence_20261017.md):

| Containing maximal-edge row | Selected quadrilateral rows | Current excess over one per edge |
|---:|---|---:|
| 82048 | 64334 | 2991/200000000 |
| 80010 | 82719 + 64876 | 9/2500000000 |

The first obstruction is acuteness alone. The second uses acuteness, a
forbidden finite-angle interval, and the necessary integrality of arithmetic
two-cycles. These invalidate the stated witness; this package does not
repair it or infer a dimension bound.

The [source](../scripts/fresh_sqrt5_scaled_core_cap_incompatibilities_20261017.py)
and [certificate](fresh_sqrt5_scaled_core_cap_incompatibilities_20261017.json)
check all 17 represented actual base occurrences of the three literal
profiles. They yield 15 capacity-one columns, three with positive load and
two violated. The coverage is this selected family, not all cap profiles.

## Why a complete contact vector determines the cap

The full star of the containing maximal simplex is a basis of the ambient
Lorentz space. Here the maximal simplex is an edge: its core has exactly
two positive-definite one-root deletions, labelled 0 and 1. The rest of its
star consists of orthogonal finite spectator components.

Each selected cap has contacts with the core and one isolated A1 root q,
and is right to all other spectators. Fix a K-valued core Gram H, with
diagonal entries 2s_i, and write its inverse as J. A cap scaled to square
2t has core covector p. Its core projection v is unique and satisfies

\[
v^2=p^T Jp.
\]

If the normalized squared contact with q is tau, the remaining component
has square t tau/2. The source checks the exact norm identity

\[
p^T Jp+t\tau/2=2t.
\]

Two caps on distinct orthogonal A1 copies therefore have product
`h=p^T Jp'`. Their normalized two-cycle is

\[
\frac{h^2}{tt'}.
\]

On the same A1 copy their scaled product instead equals

\[
h+\tfrac12\sqrt{tt'\tau\tau'}.
\]

All gauges, norm identities, contact squares and signs in these formulas
are checked exactly over Q(sqrt(5)), with exact square-root tests when
needed. The cap scales and core diagonal scales are checked totally
positive. The source also checks the physical Lorentz signature through a
positive-definite deletion and negative determinant, the positive-definite
conjugate core, and the complete core automorphism group preserving the
two deletion roots. No assumed numerical normalization enters the gate.

## Raw254: a native label5 coface gives capacity one

The seven-root core is core555. Its dotted normalized squared product is
`144+64sqrt(5)`, whose physical square root is `8+4sqrt(5)`. Thus this
cyclic core has a literal norm2 Gram over K. Its only nontrivial
automorphism exchanges roots 2 and 4 and fixes the two edge vertices.

The cap has native label5 contact with its isolated A1 coface, so
`tau=(3+sqrt(5))/2`. This finite contact was outside the earlier
single-A1 family restricted to dotted cofaces. There are exactly two core
orientations. They contact core root 1 and respectively root 2 or root 4
by label3, and are right to the other core roots.

The two core projections satisfy

\[
v_0^2=v_1^2=(5-\sqrt5)/4>0,
\qquad
v_0\cdot v_1=(3-\sqrt5)/4>0.
\]

Consequently any two distinct caps in these orientations, on different
A1 copies, have positive mutual product. On the same A1 copy the extra
term is positive as well. Neither is possible for distinct walls of an
acute Coxeter polytope. Repeating the identical orientation at the
identical copy reproduces the identical normal, not another incidence.
There is therefore at most one such cap across all A1 copies and both
orientations of a fixed containing edge.

For base type114227, with spectators `2A1+A5`, the literal selected
quadrilateral is raw254, actual provenance1505, row64334. Its unique base
of this type belongs to edge row82048. Its mass is2991/100000000, twice
the current edge mass2991/200000000.

## Raw259 and raw262: a mixed capacity-one family

These two profiles share core371, whose five-root graph has one dotted
edge of normalized squared product `3+sqrt(5)` and a simple trivalent
tail. The full core automorphism group only exchanges the two tail
leaves; it fixes both profiles and both edge vertices. Each profile
therefore has one core orientation.

In the saved tree gauge the raw259 cap has scale t=1 and the raw262 cap
has scale t'=3+sqrt(5). The normalized squared products of caps on
different A1 copies are:

| Pair | Normalized two-cycle | Exclusion |
|---|---|---|
| 259 with259 | (3−sqrt(5))/8 | Not an algebraic integer |
| 259 with262 | 1/2 | Not an algebraic integer |
| 262 with262 | (3−sqrt(5))/2 | Negative product strictly between−1 and0 |

The last interval contains no norm2 Coxeter wall product: a nonzero
finite contact has magnitude at least1. In the first two rows, arithmetic
two-cycle integrality is a necessary condition in the fixed ambient field
K=Q(sqrt(5)). On the same A1 copy, all three products are strictly
positive, as checked using the exact radical sign formula above.
Thus at most one cap from the union of these two profiles can occur at
a given containing edge.

For base type48730, with spectators `2A1+A2+D5`, the two actual
provenances are1205 and1519, mapped to quadrilateral rows82719 and64876.
Each contributes exactly one base of this type. Their total mass exceeds
the current mass of edge row80010 by9/2500000000.

## Incidence and provenance accounting

The canonical full quadrilateral star determines the selected base and
its maximal containing simplex. Every source occurrence is matched using
the literal deleted normal, opposite base facet, core isomorphism and
complete finite spectator context. The source reconstructs the union of
the relevant K2,2 corners and checks that the base vertices are exactly
the two spherical deletion roots. The frozen context table gives exactly
one base of the indicated type for each selected provenance occurrence.

Summing the pointwise capacity-one statement over actual containing
edges gives

\[
\sum_Q m(Q,S)f_Q\le f_S.
\]

The saved columns retain the actual provenance map and the literal base
multiplicities. They neither select a new representative of a merged
row nor multiply a cap count by the number of spectator copies. Missing
profiles have coefficient zero; this weakens the necessary inequality.
The complete core orientation groups and every isolated A1 position are
covered by the pair arguments, so the capacities are independent of
the number of those copies.

## Validation and scope

Generation and a fresh saved replay check the exact pair arithmetic,
core gauges, automorphism and deletion censuses, all literal contact
maps, the17 occurrences and all15 columns. The saved model masses are
unchanged. The new constraints do not assert ambient realizability of
the local configurations, complete geometric coverage, an optimized
bound, or a compact arithmetic dimension theorem.

The construction agent's separate narrow source review passed: it checked
the gauge and same-coface formulas, full orientation orbits, literal
provenance maps and occurrence accounting. It did not repeat the arithmetic
or census. A positive self-product excludes two distinct walls, not the
one wall determined by an identical literal placement; that placement and
its base incidence are counted only once.

**FINAL/FROZEN:** generation and fresh saved replay both passed, as did
source syntax, all ten direct source hashes and the local report links.

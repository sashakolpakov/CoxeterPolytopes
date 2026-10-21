# Exact repair of the two scaled-core edge-cap violations

The [new scaled-core cap constraints](fresh_sqrt5_scaled_core_cap_incompatibilities_20261017.md)
reduce the certified partial-model obstruction from12.11531601432504 to
**12.11501285634343 > 12**. This is an exact repaired witness, not a new
optimum or a dimension theorem. No optimizer or repeated large joint-port
scan is used.

The [source](../scripts/fresh_sqrt5_scaled_core_cap_edge_repair_20261017.py)
and [saved certificate](fresh_sqrt5_scaled_core_cap_edge_repair_20261017.json)
bind the [preceding52-family incidence witness](fresh_sqrt5_core23_25_spectator_incidence_20261017.md),
the inherited full matrix, the complete local packing covers, and all cap
inputs. They save the full new masses, all literal packing mixtures,
and exact checks of every retained signed and incidence inequality.

## Only two maximal-edge masses change

The simultaneous minimal simplex-mass repair of the15 new capacity-one
columns changes exactly

| Simplex row | Old mass | New mass | Increase |
|---|---:|---:|---:|
|80010|150071/10000000000|150107/10000000000|9/2500000000|
|82048|2991/200000000|2991/100000000|2991/200000000|

Both rows represent maximal edges, so their curvature demand is zero.
Every one of their endpoint blocks has one port. The source checks these
facts against the literal simplex features and the rebuilt inherited
matrix, rather than inferring them from row names. All336 positive actual
quad masses and every other simplex mass remain fixed. The positive
support stays2524 faces:2188 simplex rows and336 actual quad rows.

## Complete local packing coverage at the changed masses

The certificate replays all3712 inherited packing-mixture terms, checking
their literal masks, pairwise disjointness, old required marginals, and
total masses. The new simplex masses change exactly four one-port orbit
marginals. A positive deficit is covered by adding the corresponding
literal singleton packing. This adds four terms, for3716 total terms,
and increases three vertex allowances: indices1614,1656,1875.

The resulting mixtures cover all4169 block-orbit marginals on all648
positive simplex forests. Group averaging the literal whole packings
gives the same complete nonnegative orbit-local charge guarantee as in
the frozen predecessor. No previous cover is simply reused at changed
masses: its literals and marginals are rechecked before the singleton
increments are added.

The vertex allowances start at the final incidence witness, which already
dominates both the old local-cover allowances and the old joint-port
allowances. They only increase. A direct rebuilt-matrix check verifies
every inherited positive-cost column and all94722 signed columns.

## Why the17 joint families need no new large scan

All selected joint3/4/5-port simplex coefficients depend on a block
intersecting the tested subset in at least two ports. This includes the
three-port terms2s+t, the four-port excess term, and every simplex
intersection histogram used by the15 five-port facets. A singleton block
has intersection size0 or1, so every one of these coefficients is zero.

The source checks that identity over all286 triples,715 quartets and1287
five-element subsets of a13-port set. A port relabelling covers every
singleton. Since the only changed simplex blocks are singletons, all
joint loads are exactly unchanged. The quad masses and their actual
provenances are also unchanged; the vertex allowances only increase.
Thus the predecessor's full17-family verification remains valid. This
argument does not claim that those17 families exhaust all possible joint
geometric restrictions.

## Other retained constraints and result

The exact verifier checks all2071 cap columns:2056 inherited columns plus
the15 new ones. It also checks every one of the52 forced-quadrilateral
incidence inequalities under the final masses. The common actual
provenance-to-face map is rechecked for each new cap, with no selection of
a different representative.

The objective numerator is unchanged:

\[
3842394855727/309375000000.
\]

The exact normalizer and rational ratio are saved in the JSON. Fraction
arithmetic proves the resulting ratio is strictly greater than12; its
decimal value is12.11501285634343. This remains an obstruction within the
stated partial model. Complete geometric coverage and further global
compatibility conditions are not supplied by this calculation.

**FINAL/FROZEN:** generation and fresh saved replay passed. All12 direct
source hashes, all4 local report links, and source syntax passed. An
independent agent's narrow conceptual accounting check accepted the
singleton zero-feature and monotonicity argument; it was not a source
audit or duplicate arithmetic replay. The new geometric cap family has
the separate review scope stated in its frozen report.

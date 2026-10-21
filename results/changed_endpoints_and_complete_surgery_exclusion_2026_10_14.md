# Changed endpoints and complete exclusion of the original surgery family

**14 October 2026. No new finite-volume or compact Coxeter polytope has been
found.** The main exact result is the exclusion of the entire specified
Bugaenko wedge-first, two-truncation family in compact H9. Changed H20
constructions give useful new contact formulas, but each tested chamber
still has a certified escape. This continues the
[preceding cap-layer checkpoint](balanced_layers_and_compact_carrier_closure_2026_10_14.md).

## Literature constraints remain in force

Esselmann's restriction is for **noncompact arithmetic** reflection groups:
their dimensions satisfy `n <= 19` or `n = 21`. In particular, arithmetic
H20 and noncompact arithmetic H22+ are impossible. Vinberg's bound for
all arithmetic reflection groups, and independently for all compact
Coxeter polytopes, is `n <= 29`. The dimension-21 theorem does not exclude
compact arithmetic examples in dimensions 22 through 29. See the
[dimension audit](arithmetic_dimension_bounds_and_effective_high_dimension_search.md)
and [Belolipetsky, Section 1](https://arxiv.org/pdf/1506.03111).
The two dimension-29 statements also appear explicitly in
[Vinberg's original abstract](https://www.mathnet.ru/rus/mmo444).

The sources checked still give compact dimension 8 and noncompact
dimension 21 as the records, with dimension 20 absent from the known
finite-volume dimensions. The active construction targets are compact H9
and nonarithmetic finite-volume H20. No lattice enumeration in the excluded
noncompact arithmetic dimensions is warranted.

Vinberg's algorithm enumerates mirrors in a specified arithmetic setting.
[Allcock's published edgewalking](https://web.ma.utexas.edu/users/allcock/research/edgewalking3.pdf)
also takes a Lorentzian integral lattice, rather than arbitrary incidence
data. Traversing edges is a useful geometric strategy beyond that setting,
but requires its own exact metric construction and global completeness
proof. Rational wall coordinates alone do not imply arithmeticity.

## Compact H9: the old family is completely excluded

The [full report](fresh_h9_bugaenko_remaining_20261014.md) starts with the
specified known eleven-facet Bugaenko compact H8 incidence, wedges any
facet, and truncates two proper faces. The original low-codimension census
had 73,513 safe recipe representatives: 67,527 were excluded by a minimal
nonface of order greater than five, and the remaining 5,986 represented
2,953 incidence types. Every one of those types is now excluded.

| Exact exclusion | Distinct types |
|---|---:|
| Two-Lannér-block inertia theorem | 2,365 |
| Previously checked complete searches outside that theorem | 112 |
| Newly replayed complete searches outside that theorem | 476 |
| Total | **2,953** |

The block theorem classifies 767 pairs of internal Lannér types for the
stated block orders. All 42 possible nonright cross completions force
another negative direction. Exact polynomial signs handle unbounded
finite labels; the argument does not restrict the field or dotted
distances. A separate Fraction/radical implementation and an additional
[proof and arithmetic review](fresh_h9_bugaenko_remaining_20261014_block_gate_review.md)
pass. The remaining search records have complete label exclusions,
including the formerly limited cases.

A truncation-reordering identity extends the conclusion to arbitrary
codimensions of the two proper cuts. If a first cut has codimension above
five, any compact-Coxeter-compatible second cut must remove its large
minimal nonface in the prescribed nested way. Reordering then gives two
cuts of codimension at most five, already in the closed census. The
[normalization review](fresh_h9_bugaenko_remaining_20261014_cut_normalization_review.md)
checks the argument and 6,160 literal vertex/cardinality cases. This is a
theorem about the specified surgery family, not all fourteen-facet H9
polytopes or all compact H9 polytopes.

The changed order is **truncate H8, wedge the new truncation facet, then
truncate again**. Its census gives 5,536 incidence types. Of these, 1,040
coincide with the old catalogue, 1,220 further types fail the block theorem,
and a 32-type pilot ends in complete exact label exclusions. A separate
pilot audit rebuilds its literal complexes and repeats the arithmetic.
The [frozen frontier](fresh_h9_bugaenko_remaining_20261014_checkpoint.json)
lists **3,244 still-open types**. No full necessary-label assignment or
metric candidate survived the pilot.

## H20: endpoint changes and mixed-factor contacts

The [continuous endpoint-replacement report](fresh_h20_endpoint_replacement_20261014.md)
shows why simply lengthening radial chains cannot complete the balanced
`A6 × A6 × A7` cusp. It also constructs a single-factor edge/triangle
hierarchy whose endpoint radius can be coupled to make adjacent triangle
caps meet at `pi/5`. The large numerical pair test passes, but a single
exact strict escape survives the entire continuous parameter interval,
all vertex-radial caps and all saved middle/motif caps. A
[targeted independent review](fresh_h20_endpoint_replacement_20261014_review.md)
checks the continuous argument and key rational identities. No full exact
chamber certificate is claimed for the numerical hierarchy.

The [mixed-factor construction](fresh_h20_endpoint_replacement_20261014_mixed.md)
leaves that excluded hierarchy. It shifts simultaneously along nearest
edges in both A6 factors, with radius `R = 14 - 4sqrt(7)`. Its 392-cap orbit
has exact mutual `pi/3` contacts. Replacing the endpoint radius by

```text
rho = R * (sqrt(119 + 49cos²(pi/8))/7 - cos(pi/8))
```

gives 392 exact `pi/8` contacts with the endpoint layer. Exact field
arithmetic, rational radical bounds and symmetry establish the pair
relations and a coherent 1,202-facet chamber, preserving the 392 old
finite endpoints. However, none of the mixed caps belongs to a new finite
vertex: its finite-neighbor graph has too small a clique. The explicit
strict escape near an A7-edge midpoint proves infinite volume, even with
all allowable vertex-radial caps added. The obvious matched three-factor
extension has positive pairing with the retained mixed cap and requires
replacement rather than addition. The separate
[mixed-construction review](fresh_h20_endpoint_replacement_20261014_mixed_review.md)
states the scope of its independent checks.

The [integer-gap theorem](fresh_h20_integer_inner_gap_20261014.md) is stronger
than the earlier endpoint-dependent exclusion. Every central-compatible
radius-two cap with nonnegative integer seed gaps summing to `7,7,8`
misses `3z_000/20` by power at least `1/100`. The complete exact cutoff
contains 651 recipes; only two pass the central gate, both attaining the
bound. A [separate arithmetic reconstruction](fresh_h20_integer_inner_gap_20261014_review.md)
passes. Thus changing endpoint radii cannot rescue this cap family. The
two noninteger contact scouts in that report fail only greedily and remain
numerical bounded-search results, not family exclusions.

## A changed seed and a stronger combined covering obstruction

The [weighted-box seed](fresh_h20_seed_deformation_20261014.md) changes the
nineteen cusp factors to A1, with inradii ten times 1, once `sqrt(2)` and
eight times `phi`. Its central contacts have labels 3, 4 and 5. The initial
report excludes the radius-four vertex and radius-two nonvertex families
separately. The new [combined-family theorem](fresh_h20_weighted_mixed_probe_20261014.md)
now excludes their union, and indeed permits **every seed-compatible cap
centered at any old vertex, with arbitrary radius and finite label**, as
well as every seed-compatible radius-two nonvertex cap with no label cutoff.

There are `2^19` half-vertex probes. Fixing the eleven unit/heavy signs
leaves an eight-dimensional golden-sign cube. Let `A` be its probes not
covered by a radius-four ball. Any putative complete covering must cover
each member of `A` by a radius-two cap. The exact adjacent-probe rule
assigns each induced cube edge to one endpoint's distinguished golden
coordinate; each endpoint can receive at most one edge. Thus
`e(A) <= |A|`. Degree counting also gives `e(A) >= 8(|A|-128)`, so
`|A| <= 146`. Radius-four balls must consequently cover at least
`110 * 2048 = 225280` probes. But at most 512 such caps are mutually
compatible, and each covers exactly 275 probes: their total coverage is
at most `140800`. This is a contradiction. Equivalently, radius-two caps
cover at most `146 * 2048 = 299008` probes, so the entire combined family
leaves at least **84,480** of the 524,288 probes uncovered.

The arbitrary-vertex-radius extension is short. An opposite unit seed
wall is at distance 2 from every old vertex. Its allowed contact forces
either radius 4, or radius at most `sqrt(8)`. Every half-probe is at squared
distance at least `5+2phi > 8` from every old vertex. All the other radii
therefore miss every probe and cannot repair the counting deficit.

The all-label radius-two classification uses the smaller of the two
opposite seed gaps. Since it is at most the inradius, its label is bounded
by 3, 4 or 5 in the three coordinate types. The only additional proposed
complementary gaps, `2sqrt(2)-1` and `2phi-sqrt(2)`, lie strictly between
`2cos(pi/7)` and `2cos(pi/8)` and so give no integer label. The original
coordinate catalogue is consequently complete for every finite label. A
[separate proof and arithmetic review](fresh_h20_weighted_mixed_probe_20261014_review.md)
checks the count and extensions. This rules out the combined dictionary;
nonvertex caps of other radii and a changed core remain
outside the theorem.

Useful next work must leave the closed families. The changed-order H9
catalogue gives a concrete finite list to investigate. A weighted-box
construction using nonvertex caps of other radii, or a genuinely different
H20 endpoint mechanism, remains distinct. New H20 walls should be judged by their ability to
produce complete finite vertices and full-rank cusps, followed by a global
finite-volume check. More pairwise-compatible facets or coverage of a
finite probe set alone does not establish the requested example.

## Frozen archive verification

Individual reports distinguish exact producer replays, independent
arithmetic, proof reviews and numerical scouts. The
[checkpoint auditor](../scripts/audit_changed_endpoints_checkpoint_20261014.py)
checks local links and recursively declared hashes; its
[saved result](changed_endpoints_and_complete_surgery_exclusion_2026_10_14_audit.json)
checks provenance only. It does not replace any mathematical certificate.

```sh
python3 scripts/audit_changed_endpoints_checkpoint_20261014.py --verify
```

The requested existence goal remains active. No chamber in this checkpoint
is claimed as a new finite-volume example or as a new dimension.

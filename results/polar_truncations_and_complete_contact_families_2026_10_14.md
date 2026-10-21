# Polar truncations and complete contact families

**No new finite-volume Coxeter polytope or dimension is established.** The current H20 chamber has 100 actual facets, at least 134 verified finite vertices and three preserved bounded cusps. Its strict rational escape survives, so its volume is infinite. The compact-H9 work excludes a larger family with seven unrestricted finite labels and all nine anchor-own contact patterns. It also excludes the latest saved necessary assignment; the whole incidence remains open.

This continues the [bridge-face checkpoint](bridge_face_edges_and_unbounded_label_search_2026_10_14.md), preserving its frozen mathematical artifacts. The [dimension and algorithm restrictions](arithmetic_dimension_bounds_and_effective_high_dimension_search.md) still direct the work toward nonarithmetic H20 and compact H9. Constructing more finite vertices does not by itself advance either dimension record.

## Two exact polar truncations

The [polar construction](fresh_h20_polar_pair_20261014.md) adds two rational normals to the 98-wall chamber. All 197 new pairs satisfy the Coxeter conditions: 53 right, four at pi/3, three at pi/4, twelve parallel and 125 dotted. A common strict timelike center and explicit facet points certify a coherent chamber with 100 actual facets. The three old bounded cusp links are preserved.

The [complete new-wall census](fresh_h20_polar_pair_20261014_census.json) checks 20,810 candidate twenty-wall supports with exact positive-definiteness tests. Exactly sixteen give finite vertices on the union of the two new walls: fourteen on each wall and twelve on both. The old 118 known finite vertices are strict at both additions. Hence there are at least 134 finite vertices overall; this is not a complete global census.

The [independent review](fresh_h20_polar_pair_20261014_review.md) uses a separate clique enumeration and exact component tests, reconstructs all sixteen vertices, and directly checks all 1,600 wall inequalities there. It also rechecks the 197 pairs, polar extreme-ray ranks, actual facets, cusp preservation and escape.

The unchanged rational escape has

```text
E^2 = 20181174579467302741 / 10^18 > 0.
```

Its pairings with the new normals are strictly negative, as are all 98 inherited pairings. Together with the future orientation and common timelike center, this proves infinite volume. Both new normals also project future timelike to the retained S75 face, which remains an entire H2 plane. These truncations do not close either obstruction.

The [numerical polar scout](fresh_h20_polar_frontier_20261014.md) sampled the full core and variants deleting wall 8, wall 15, or both. Its four runs used 15,288 LP objectives and found no passing escape cutter. The changed-core passes recovered deleted walls or the same two polars numerically. These are bounded diagnostics, not nonexistence results. The source coordinates require a 90-decimal-digit Lorentz chart before the floating-point stage; only the separate exact certificates accept walls.

## Global restrictions on face completion

Three exact results now constrain a retained-core completion.

1. The [S75 support gate](fresh_h20_s75_projection_gate_20261014.md) allows at most two nonright support contacts for a cap inducing a boundary geodesic on S75. Its complete necessary profile list has seventeen discrete profiles and one isolated-node family with arbitrary finite label. Every such cap must be right to ten specified support walls, including walls 8, 15 and 75.
2. The [escape dual identity](fresh_h20_s75_completion_20261014.md), checked by an [independent physical-coordinate and rank review](fresh_h20_s75_completion_20261014_gate_review.md), proves that every S75 boundary-producing Coxeter cap strictly misses the current escape. In a spanning expression for E, only the coefficients at walls 8 and 15 are negative. Every acute Coxeter escape cutter must therefore be nonright to at least one of those walls. The S75 inverse-Gram diagonal bounds force that cutter's projection to be future timelike, so it leaves the whole S75 plane strict. The statement covers arbitrary real normals and lengths satisfying the Coxeter conditions, with no label cutoff.
3. The [concurrent-sector theorem](fresh_dimension_global_routes_20261014.md) excludes every finite coherent extension of the 98-wall chamber using only normals orthogonal to S61. On S75, each effective new boundary passes through the same finite point; their feasible intersection retains an angular sector of positive width. An explicit perturbation gives an ambient strict spacelike escape. Closing S61 locally by perpendicular caps therefore cannot finish the ambient problem.

The second result continues to apply after the two polar additions, since the required old walls and the same escape are retained. The third is stated and certified for the 98-wall predecessor. No result excludes a mixed extension or a changed core. A successful retained-core construction needs distinct caps that bound S75 and cut the saved escape, with compatible mutual angles and global completeness; testing either role in isolation is insufficient. Ordinary finite doubling of an infinite-volume chamber cannot repair its volume.

## Compact H9: complete own-contact patterns and a sparse obstruction

The [complete own-contact theorem](fresh_h9_t8_asymmetric_20261014.md) covers all nine nonempty own-support patterns for anchors 11 and 17 on the last two-node block. The first internal label ranges over every finite integer at least seven; the last ranges over every finite integer at least four; all four own labels are unrestricted finite integers at least two, with a nonright contact at each anchor. Exact domain guards reduce the remaining contacts to 408 polynomial cases. All 11,736 Bernstein coefficients are strictly positive. The [independent review](fresh_h9_t8_asymmetric_20261014_review.md) checks both the literal incidence guards and separately reconstructed polynomial arithmetic.

The [middle-block extension](fresh_h9_t8_middle_20261014.md) additionally allows every finite label at least four on pair (7,8). Label four inherits the 408 cases. Label five changes 256 cases, with 4,320 strictly positive coefficients over Q(sqrt(2),sqrt(5)); 152 cases are unchanged. At labels at least six, a foreign nonright contact into this block alone costs at least the whole positive norm, before the anchor's positive own-contact cost, and is impossible. The [independent extension review](fresh_h9_t8_middle_20261014_review.md) verifies this domain split, the inherited complement and every new coefficient. These are full finite-label ranges, not numerical cutoffs.

A bounded continuation stopped at its first surviving necessary assignment after 1,450 nodes, before reaching its 100,000-node limit. That assignment changed the contacts of leaves 6 and 15. The [sparse dependency certificate](fresh_h9_t8_sparse_frontier_20261014.md) then excluded it without another polynomial bank. Two actual four-wall faces have an affine star as a relaxed Gram matrix, with a positive exact null vector. These replace the only old guards involving unrelated leaves. The resulting sufficient hypotheses for the rank obstruction are just

```text
(1,2) >= 7, (4,5) = 4, (7,8) >= 4, (12,13) >= 4;
(1,14) = 2, (2,14) = (4,14) = (5,14) = 3.
```

The own contacts remain free with nonempty support. Other labels may vary subject to the compact incidence and Coxeter conditions. The [independent sparse-domain review](fresh_h9_t8_sparse_frontier_20261014_review.md) rebuilds all 24 cross-block guards, four anchor budgets, two third-wall budgets, two affine stars and nine own-support domains from these eight hypotheses. It accounts for all 55 selected Gram entries and binds the independently checked coefficient banks. Generation and saved replay pass. The saved assignment satisfies the hypotheses and is excluded. No subsequent assignment has been searched or asserted to survive. This is a family obstruction within the case-8 T8 incidence, not exclusion of that entire incidence and not a compact H9 construction.

## Verification and remaining work

Exact producers, independent arithmetic reviews and the numerical transcripts pass generation and saved-result replay. The [checkpoint audit](polar_truncations_and_complete_contact_families_2026_10_14_audit.json) binds current sources, certificates and reports, verifies local links and preserves the preceding frozen archive. Its role is archive integrity; it does not replace the mathematical checks. Earlier README and research-plan hashes are historical snapshots.

The next H20 work must couple the two required cap roles or change the retained core; increasing the same perpendicular or polar sample banks is not an exhaustive strategy. The next compact-H9 continuation should skip the entire sparse family above and violate at least one of its eight core conditions, or use a different unresolved incidence. The new exact obstructions should be used as search gates before further metric work. Both dimension targets remain open.

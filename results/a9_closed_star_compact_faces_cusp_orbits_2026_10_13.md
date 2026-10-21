# Closed arithmetic star, compact face, and cusp-orbit repairs

**Historical checkpoint.** The [next checkpoint](a9_nonreflectivity_compact_prism_face_orbits_2026_10_13.md)
supersedes the open A9 status below: complete edge walks now prove this full
arithmetic host nonreflective. Its old 36-orbit queue is no longer a search
priority. The geometric constructions have also advanced, still without a
finite-volume example.

**No new compact or finite-volume Coxeter dimension has been found.** This continuation advances the arithmetic A9 chamber walk and constructs exact geometric repairs in H9, H20, and H22. Every saved finite wall set still has an exact strict spacelike escape. The A9 full reflection chamber remains open; the finite-prefix escape is not a nonreflectivity proof.

## Arithmetic H9: a complete central star

For the quartic host `⟨−(T−1)⟩ ⊕ A9`, `T=2cos(π/8)`, four complete nearest-mirror searches and the integral reversal resolve all eight previously unexplored directions at the central A4+A4+A1 corner. They inspect 380 finite-angle profiles and 900 exact relative-orbit interval points. Certified PARI norm/unit interfaces and an independent Fraction replay establish firstness, not merely compatibility of a chosen endpoint.

The coherent prefix has **19 walls and ten known actual finite vertices**, including eight new vertices. All nine edges at the central corner are closed. The complete outgoing queue at the ten known corners has 90 directions: 18 are the two orientations of those nine edges, and 72 remain unresolved, or **36 reversal orbits**. A bounded 360-map corner-frame census finds only identities and the known reversal among integral maps; it supplies no new global chamber symmetry.

The [A9 report](fresh_h9_a9_next_20261013.md) contains the roots, exact intervals, independent geometry replay, and persistent queue. This is the strongest current path toward determining a full new arithmetic chamber. It remains distinct from the two previously excluded F4 hosts.

## Geometric H9: an actual compact tetrahedral face

Four coupled caps extend the previous 90-wall configuration to **94 actual facets**. They cut the old escape and form an actual compact H3 tetrahedral face with diagram `3–5–3`, four finite vertices, and six finite edges. This is a lower-dimensional face, not a compact H9 example and not a claimed new tetrahedral type.

Every one of the 366 new pairs is exact. The complete 960-subset census meeting a new wall has no affine circuit. The four new vertices expose 24 additional open directions, and a new strict escape of norm approximately `0.0464919891` survives. The [construction](fresh_h9_compact_next_20261013.md), [separate symbolic and 160-digit review](fresh_h9_compact_next_20261013_review.md), and [second source review](fresh_h9_compact_next_20261013_source_review.md) distinguish the compact face from the infinite-volume whole chamber.

## Geometric H20: coupled polar repairs

The old gap polar failed only two wall-angle checks. Instead of accepting its forbidden angles, two coupled formulas combine successive exact polars with their offending old wall normals. This retains all 96 old walls and adds two compatible walls, reaching **98 actual facets**.

The first new wall lies in `Q(√6,√10)` and has norm two; the second lies in `Q(√6)` and has norm one. All 193 new pair relations, 214 signs at inherited cusps, and two successive strict target cuts are exact. A third strict escape has norm

`698631401 / 1887436800 > 0`.

The first new wall also creates a new deficient cusp: its link has quotient
rank 17 rather than the required 19, and the second wall does not repair it.
The [H20 report](fresh_h20_gap_next_20261013.md) explains the formulas,
this limitation, and the finite scouting scope. The [source review and exact replay](fresh_h20_gap_next_20261013_review.md)
pass. The construction changes the local repair mechanism, but does not
establish finite volume or close the entire edge frontier.

## Geometric H22: bounded inherited cusps, remaining global escape

The preceding H22 caps had timelike old-space projections and could not bound the added direction at any inherited cusp. Paired walls with null projections repair that defect. Completing 15 sampled symmetry orbits gives **1,065,708 bounded inherited cusp links**, with 2,131,416 new walls and **2,252,588 actual facets** in total.

The exact compatibility proof uses the negative integral product of distinct primitive null rays and orbit transfer; it does not attempt a quadratic enumeration of over two trillion pairs. One orbit uses height one quarter rather than one to avoid forbidden angles. Each old affine link of rank 20 gains an orthogonal affine A1, producing a bounded link of rank 21.

The new walls cut the preceding strict escape, but a new rational escape of square `1656895868840709 / 10000000000000000` works for every center label `m≥3`. The [construction report](fresh_h22_cusp_orbits_20261013.md) and [independent source/replay review](fresh_h22_cusp_orbits_20261013_review.md) state the limitation: these are complete orbits of sampled cusps, not an exhaustive old-cusp census or a finite-volume completion.

## A dimension-lift shortcut excluded

A [new local wedge argument](fresh_dimension_wedge_20261013.md) rules out replacing just one old facet by two sloping facets while retaining all others. Any finite nonright old label at that facet must be 4, 6, or 10, without an upper-label cutoff. A label-six facet with a retained parallel neighbor is also impossible. Exact Gram checks exclude this ansatz for all 11 facets of the stored H8 example, all 36 facets of H19, and all 210 facets of H21. General coupled lifts and changed chambers remain possible.

## Next search priorities

1. Advance the A9 arithmetic queue through its 36 unresolved reversal orbits, keeping the complete norm/unit and firstness gates. The eight new mirrors are perpendicular to the central cap, so its H8 facet may provide useful structure, but is not identified with Bugaenko's golden-field example.
2. Use the 24 open directions at the new H9 tetrahedron to seek connections to inherited vertices, with exact metric constraints and backtracking when local choices obstruct completion.
3. Investigate the new H20 and H22 gap supports. Repaired inherited cusps do not control hyperideal gaps, and increasing facet counts alone supplies no evidence for finite volume.

The [Vinberg–Allcock audit](vinberg_allcock_search_strategy_20261013.md) remains the algorithmic guide: distance-ordered mirror enumeration and edge traversal are complementary. Allcock's published algorithm is a lattice algorithm; our unrestricted diagram/metric constructions use separate proofs of compatibility and volume.

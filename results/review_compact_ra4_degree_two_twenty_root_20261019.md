# Root review: degree two forces W>=20

19 October 2026. **PASS, full written-source review.**
Source: `compact_ra4_degree_two_excess_twenty_20261019.md`.
The underlying L6 metric/incidence lemma is source-bound in the note
and has separate root and independent written reviews.

For k>=7 the inherited W>=5k-15 suffices. At k6 the fifteen selected
facets are distinct, with excess at least two, giving baseline thirty.
Each of the two disjoint lateral rows with m non-L6 neighbors contributes
at least m extra facet excess. Its remaining 6-m L6 neighbors meet the
base in actual hexagons by the unique-longest-edge and third-parent
argument. The base has excess at least 7-m, hence contributes at least
max(0,5-m) above baseline. Every row therefore adds at least five.
Adding the two distinct sets of facet excesses gives 2W>=40.

This sum uses distinct facets, not disjoint ridges; double counting
ridge excess is already part of sum d=2W. It permits unselected facets
and arbitrary additional excess. The proof is valid and conditional
on a degree-two vertex in G. The unconditional positive-excess gap
remains17. No numerical certificate, new minimum-volume theorem,
manuscript integration or Lean verification is claimed.

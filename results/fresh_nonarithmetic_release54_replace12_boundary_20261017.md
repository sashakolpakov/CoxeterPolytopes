# Exact S6 boundary: seventeen open ends become sixteen

The [exact replacement](fresh_nonarithmetic_release54_replace12_20261017.md)
closes the entire actual two-face S6+8 as a compact heptagon. Its side cycle
is `7,11,10,44,65,52,new12`. Every previously compact S6 side-face remains
compact.

The [checker](../scripts/fresh_nonarithmetic_release54_replace12_boundary_20261017.py)
and [certificate](fresh_nonarithmetic_release54_replace12_boundary_20261017.json)
replay the complete residual4 enumeration, using the exact new normal.

| Intrinsic boundary quantity | Before65 | After64 |
| --- | ---: | ---: |
| Finite vertices | 23 | 24 |
| Exterior projective rays | 9 | 8 |
| Compact geodesic edges | 28 | 30 |
| One-open edges | 13 | 12 |
| Two-open edges | 2 | 2 |
| Outside-H3 projective edges | 5 | 4 |
| Geodesic open ends | 17 | 16 |

The complete projective combinatorial graph still has18 sides,32 simple
rays and48 edges. Only three edge types change:

* `[8,12]`: one-open to compact;
* `[8,52]`: one-open to compact;
* `[12,52]`: outside H3 to one-open.

The exterior component affected by deleting34 is now the three-ray path

`[7,12,47] — [12,47,52] — [30,47,52]`.

Its two internal projective edges remain outside H3. The other five-ray
exterior component is unchanged.

Completeness uses the same certified chart as the preceding exact ledger:
the four right-to-S6 normals8,30,44,52 span the Lorentz4 complement; the sum
of positive-side covectors is a strict dual chart. All816 triples of the18
side equations give nonsingular chart systems and are either strictly
rejected or produce a simple ray. There are no ideal rays, unresolved
signs or nonsimple incidences. The complete ray list determines all edges;
their two-ray Lorentz Gram matrices classify their positive cones exactly.
The other40 nonwhite projected walls are strictly future timelike, hence
redundant on the entire future H3.

The auxiliary exterior cone omits those causal inequalities outside H3.
Its exterior rays are therefore not claimed to be the full ambient cone
slice. This is a complete intrinsic boundary calculation, with a separate
complete ambient finite census197 and strict escape supplied by the parent
construction.

Generation and fresh saved replay passed. The frozen exact boundary engine
is reused; no separate arithmetic implementation is claimed.

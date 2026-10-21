# XY-to-8Z flip: bounded triangular-facet discovery

The changed disk has no native-admissible Z in the stated finite search. This is a numerical discovery result, not an all-label obstruction. The exact 63-wall, 199-finite-vertex, 10-open-end control remains unchanged and has infinite volume.

FINAL/FROZEN. The source is [fresh_old8_triangle_flip_20261017.py](../scripts/fresh_old8_triangle_flip_20261017.py), with [saved data](fresh_old8_triangle_flip_20261017.json). Generation90528 and fresh saved replay42533 both terminated with PASS; there were no unresolved numerical records.

## Changed incidence and geometric recipe

Starting from the [old30 connector disk](fresh_old30_connector_flip_20261017.md), remove triangles8XY andXYZ and insert8XZ and8ZY. The new diagonal8Z was absent. Literal edge incidences, all vertex links, the boundary cycle7–8–44–63–30, and Euler characteristic are checked in the source. The disk has13 vertices,31 edges and19 triangles, with x=5, y=8, z=6. Its13 old–new and13 new–new edges give the conditional change2(13+13)−3·19=−5 in open ends. This identity requires an actual embedded simple patch with the stated incidences; it is not a geometric existence claim. The retained old30–44 chord lies outside the proposed disk.

Y has exactly the three neighbors8,R,Z. The recipe chooses Z with finite native contacts to8,R,C, then solves Y right toS6,8,R,Z. Since8 andR are orthogonal norm-two normals, the projected8/R/Z Gram is Lorentzian precisely when

`qZ − (theta(Z8)^2 + theta(ZR)^2)/2 < 0`.

Its perpendicular line is then spacelike and has two norm-fixed polar signs. A common strict center selects the possible chamber side. Positive-definite corner supports alone are insufficient for a recorded compact face: any surviving prefix must also pass the complete cone and independent finite-support tests. No Z reached that stage here.

## Domain and result

The27 A/C events are the frozen [C10 endpoint events](fresh_c10_endpoint_disk_20261017.md); the source rebuilds139 native X events. R44 ranges over6,7,8,9,10,12,15. Z uses all seven A2 profiles; its assigned labels to8 andR range over2,3,4,5,6,8,12, and its C label ranges over2 through6. Both norm branches are tested. This is an explicit bounded set, not a complete finite-angle domain.

There are7,280 Z equations and14,560 norm branches. All fail retained native gates:13,712 have an obtuse contact, while848 are acute but have a non-Coxeter old contact. The first recorded failing wall is10 for6,582 branches,30 for4,690,7 for3,162 and4 for126. The data save representative acute failures. Several have products with30 or63 strictly between−1 and0, rather than the dotted contacts required by the intended nonincidence.

The recipe skips5,880 positive8/R/Z frames and525 numerically zero polar limits. These are not silently accepted or promoted to an obstruction for other topologies. The source additionally checks native contact residuals, the fixed strict center and intrinsic nonincidences; assigned-label tests use full native products, including the A2 offset.

## Deletion cost and verification scope

The five deleted old walls are12,47,52,64,65. The inherited exact deletion ledger records11 lost finite-vertex certificates and six touched compact two-face certificates, with supportsS6 plus8,12,44,63,64,65. Replacements would need to rebuild these; the disk count does not cancel that debit. No old8 wall is deleted.

The abstract disk is checked combinatorially. The root calculations, angle matching and rejection values are floating-point diagnostics with explicit residual and ambiguity guards, using inherited numerical helpers. Generation/replay is reproducibility of this finite discovery run, not independent exact arithmetic. No new exact chamber, finite-volume example, all-parameter obstruction, or improvement over the exact control is claimed.

Root independently read the complete new source and passed its stated numerical domain and gate logic. Frontier separately checked the conditional disk flip and polar-triangle argument. Neither review independently reran the floating-point geometry or upgraded it to exact arithmetic.

The next changed duty is recorded separately in [the CZ-to-XW successor](fresh_old8_triangle_double_flip_20261017.md); that successor removes the finite C–Z requirement instead of enlarging this search's label list.

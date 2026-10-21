# Bounded bridge-face polygon diagnostic

The frozen numerical scout accepts 668 candidate normals, including 32 that cut the old escape. An [independent numerical diagnostic](../scripts/fresh_h20_e7_bridge_polygon_20261014_scout.py) recomputes all 120,438 pairs sharing one of the two specified bridge-face supports. It reproduces exactly 40,265 compatible pairs; every one is dotted. There are no finite-angle or parallel pairs in this pool.

Consequently this bounded pool supplies no finite or ideal corner of a Coxeter polygon on either H2 face, and cannot close either face into a finite-area polygon. The [saved diagnostic](fresh_h20_e7_bridge_polygon_20261014_scout.json) passes generation and saved-result replay and binds the frozen [source scout](fresh_h20_e7_bridge_face_20261014_scout.json).

All arithmetic in this diagnostic is floating point. It does not establish an exact pair-exclusion theorem or exclude other anchors, unsampled labels, or caps with nonright contacts to the eighteen face roots. Its purpose is to prevent a repeated search for a polygon within this already tested pool. The [exact face review](fresh_h20_bridge_faces_20261014_review.md) and the [exact 96-wall construction](fresh_h20_e7_bridge_face_20261014.md) have separate proof scopes.

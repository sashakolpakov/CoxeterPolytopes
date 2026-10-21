# A disk replacement for the whole B region

The control remains the [64-wall chamber with 198 finite vertices](fresh_nonarithmetic_replace47_close12_20261017.md), whose [complete S6 boundary](fresh_nonarithmetic_replace47_close12_boundary_20261017.md) has 15 open ends. This investigation finds a previously untested class of numerical bridges with only one old intrinsic neighbor, but no improved chamber is promoted.

The source rechecks the saved control counts, the [deletion-cost audit](fresh_nonarithmetic_b_release_costs_20261017.md), and the [three-cap incidence gate](fresh_nonarithmetic_b_release_incidence_20261017.md). Deleting52 loses five finite vertices and the four compact polygons on sides8,12,64,65; it loses no higher-dimensional compact face. The retained intrinsic boundary has 39 geodesic edges and20 finite vertices, hence18 nonfinite-end incidences.

## The conditional disk identity

The exposed B-region endpoint duties lie on the seven-edge cycle

`7 — 47 — 12 — 8 — 65 — 64 — 30 — 7`.

Consider a proposed replacement whose dual patch is a triangulated disk, with these genuine distinct old vertices and edges as its boundary, no old chords or old-only interior triangles, and f new interior vertices. Assume all old incidences are retained, the resulting boundary is simple with no ideal vertices, and the patch has no further intersections or holes. This specifies a possible topology; it does not prove the existence of compatible normals.

Let q be the number of old boundary edges; a and b the old–new and new–new edges; and x,y,z the triangles with respectively one,two,three new vertices. Every boundary edge belongs to exactly one new triangle, so x=q. Counting incidences at interior edges gives a=x+y and 2b=y+3z. Euler's formula gives z=b−f+1. Consequently

`Δ(2E−3V_f) = 2(a+b)−3(x+y+z) = −q`.

A complete filling of this seven-duty disk would therefore reduce18 ends to11, closing the B region while leaving the A region. If every new facet has at most two old intrinsic neighbors, it supplies at most one distinct old-boundary-edge triangle, so f≥7 is necessary for this whole-disk filling. This conditional statement neither excludes other patch topologies nor strengthens the earlier three-cap gate into an unrestricted finite-volume obstruction.

The script independently checks an explicit six-interior-vertex abstract disk. Its boundary cap order is `L,Q,M,R,N,P`; L,M,N handle two old endpoint duties each; P handles one. The old neighbor sets are:

| New facet | Retained old neighbors |
|---|---|
|L|30,7,47|
|Q|47|
|M|47,12,8|
|R|8|
|N|8,65,64|
|P|64,30|

Its four all-new triangles are LQP, MQR, NRP and QRP. Every interior edge has two incident triangles, the boundary is exactly the seven-cycle, every interior vertex link is a cycle, every boundary vertex link is a path, and Euler characteristic is1. The counts are `(a,b,x,y,z)=(13,9,7,6,4)`, giving the required change−7. These integer checks certify an abstract disk only. A separate agent read the disk-counting argument and its hypotheses and found no issue; no independent geometry computation is claimed.

## A bounded search using new–new framing contacts

The frozen [three-old-contact bank](fresh_nonarithmetic_b_release52_patch_20261017.json) contains three L seeds,136 M seeds and132 N seeds with the prescribed two-endpoint duties. None of its multiple-endpoint seeds have a finite mutual contact. Thus the four-cap boundary-arc layout2,2,2,1 cannot use those sampled seeds directly. This is a statement about that saved finite sample, not an all-label four-cap exclusion.

The six-cap disk avoids those direct contacts. The new search first solves Q from its contacts to **old47 and new L,M**, rather than three old walls. It retains all seven actual-S6-facet A2 profiles, including their offsets in every native old and new product. Assigned labels are2,3,4,5,6,8,10,12. Both required projected endpoint triples must be positive definite; all63 retained native-wall gates and the shared old strict interior point are checked numerically.

Of23,904 real branches with the required finite endpoints,188 distinct Q bridges survive. They all use the all-right L seed. All188 partial L–Q–M chambers have18 open ends in the complete numerical S6 triple-intersection ledger. They do not improve the frozen15-end control and are not exactified or promoted.

Before solving R, the search tests the necessary opposite-side cap P, framed by old64, old30 and L. The labels to64 and L range through every integer2–32 plus48,64,96,128; the label to30 ranges through2–5. All seven A2 profiles and both signs are retained. Of2,610 real finite-endpoint branches,807 caps pass all retained native gates. None has the required finite join to any of the132 N seeds:105,717 pairs are dotted and807 lie in non-Coxeter gaps. Moreover, all151,716 sampled P/Q products fail native Coxeter compatibility. Their numerical range is approximately `[−0.3858471,1.1982769]`; this range itself is not an exact interval certificate. The remaining R solve is therefore unnecessary for this bounded six-cap family.

The accepted one-old-neighbor bridges show that the earlier multi-endpoint seed incompatibility can be bypassed. The tested opposite-side and central contacts do not complete the disk. Missing labels, different multi-endpoint seeds and other disk triangulations remain open. No conclusion of nonexistence, ambient finite volume, or improved geometry follows from the bounded failures.

## Reproduction and scope

[Source](../scripts/fresh_nonarithmetic_b_disk_bridge_20261017.py) and [saved result](fresh_nonarithmetic_b_disk_bridge_20261017.json) preserve the exact input hashes, the abstract disk, all numerical branch records and every scored partial boundary. Generation and fresh saved replay passed. The disk incidence verification is exact integer combinatorics; normal construction, native-angle recognition and the boundary scout use floating-point arithmetic and are explicitly not geometry certificates.

```sh
python3 scripts/fresh_nonarithmetic_b_disk_bridge_20261017.py
python3 scripts/fresh_nonarithmetic_b_disk_bridge_20261017.py --verify
```

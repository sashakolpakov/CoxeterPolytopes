# Releasing wall49 and a three-pentagon hexagon patch

This bounded construction experiment found no compatible three-pentagon replacement of wall49. A separate [exact all-label endpoint gate](fresh_nonarithmetic_release49_endpoint_gate_20261017.md) subsequently excludes both placements of this particular patch. Other topologies and further replacements remain outside the result.

The [scout](../scripts/fresh_nonarithmetic_release49_hexagon_20261017.py) and [saved output](fresh_nonarithmetic_release49_hexagon_20261017.json) start from the frozen [64-wall,198-finite-vertex control](fresh_nonarithmetic_replace47_close12_20261017.md), whose [complete intrinsic S6 boundary](fresh_nonarithmetic_replace47_close12_boundary_20261017.md) has25 finite vertices,7 exterior rays and15 open ends.

## Deletion cost and proposed topology

Wall49's link is the six-cycle [7,33,30,59,60,58]. Five corners of that side are exterior; its only finite corner is [7,49,58]. Deleting49 therefore removes exactly the one ambient finite-vertex certificate S6+[7,49,58]. An exhaustive check of the old198 finite supports and all their incident edges finds no compact face of dimension at least2 touching49; no compact S6 side face is lost. These are exact statements about the frozen old incidence data.

The numerical deletion-only boundary has24 finite vertices,6 exterior rays and10 open ends. That is an enlarged region, not an improved chamber: the removed constraint and lost finite vertex must be replaced before a construction claim could follow.

For old cycle a,b,c,d,e,f, introduce L,M,N with respective old contacts(a,b,c),(c,d,e),(e,f,a). The proposed disk has triangles

`abL, bcL, cdM, deM, efN, faN, cLM, eMN, aNL, LMN`.

The exact incidence check gives9 vertices,18 edges and10 triangles, Euler characteristic1, one cyclic boundary and disk links at every vertex. Each new side would be a pentagon. This is conditional combinatorial feasibility, not geometric realization. Both alternating cyclic placements were searched.

## Bounded numerical experiment

For each of the six endpoint frames the scout uses all seven A2 attachment profiles, labels2–12,16,24,32 and both norm branches. It checks all63 retained native contacts and the common-center sign. Native root construction and products use70-digit midpoint arithmetic to avoid cancellation; preliminary filters and complete intrinsic boundary scoring use double precision. Neither method is exact acceptance evidence.

| Frame | Real finite-endpoint branches | Acute branches | Coherent candidates |
|---|---:|---:|---:|
| [7,33,30] |88|0|0|
| [30,59,60] |8|1|0|
| [60,58,7] |10|0|0|
| [33,30,59] |92|26|0|
| [59,60,58] |390|195|195|
| [58,7,33] |26|0|0|

The27 acute rejected branches fail a native Coxeter-angle gate. Since each placement has empty mandatory pools, no mutual triple or complete disk survives and no new boundary is promoted. The exact successor independently reconstructs the two small mandatory endpoint domains and removes their label cutoff.

## Verification and scope

**Final/frozen.** Generation and a fresh saved replay terminate with PASS; source syntax, direct source hashes and local report links pass. Exact old deletion/incidence and disk-topology checks are distinguished from numerical discovery. This package does not certify the195 remaining numerical roots as actual new walls, does not exclude arbitrary patches, and supplies no new finite-volume chamber or ambient finite-vertex census.

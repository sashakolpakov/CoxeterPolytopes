# Fourth-cap subdivision of the release47 interior edge

The changed-incidence search finds586 native-coherent numerical K caps and
4,128 coherent K/R additions. **Every complete target cone still has ten
open ends.** Every accepted K–R product is dotted, so the proposed four-cap
disk is not completed. No geometry is promoted; the exact63-wall,
199-finite-vertex, ten-open-end control remains unchanged and infinite-volume.

The [source](../scripts/fresh_cu_release47_interior_cap_20261017.py) and
[saved result](fresh_cu_release47_interior_cap_20261017.json) use the frozen
[quadratic moving-wall events](fresh_cu_release47_quadratic_events_20261017.md).
They remove the earlier requirement that L meet R by introducing a new
interior cap K. This is a bounded numerical test of a new incidence domain,
not an all-label or all-parameter obstruction.

## Disk and retained geometry

The patch boundary is `7–12–52–30–A–7`. P owns the two12 arcs, L the two30
arcs, and R the7/A arc. Replace the internal edge L–R and its two triangles
A/L/R and P/L/R by A/L/K, A/K/R, P/L/K and P/K/R. The other seven triangles
remain7/12/P,12/52/P,52/P/L,52/30/L,30/A/L,A/7/R and7/R/P.

The source checks the literal disk: nine vertices,19 edges and11 triangles,
Euler characteristic1, exactly the five stated boundary edges, connected
interval links at boundary vertices and cycle links at all four interior
vertices. K has link `A–L–P–R–A`; L–R is absent. These exact integer checks
certify the abstract topology only. Its realization still requires the
four finite K corners and every other prescribed incidence.

All62 control walls other than47 are retained. The exact old-incidence
debit is copied and bound from the predecessor: deleting47 loses the old
finite supports `S6+[7,12,47]` and `S6+[12,47,52]`. The old12 quadrilateral
is the only previously certified compact face of dimension at least2
touched by this deletion. In particular, deletion does not touch the
compact8 facet or its S6 heptagon. The P3 seed equals old47 and serves as a
control; P4 and P5 are genuine replacement seeds. Restoration of lost
geometry is checked numerically here, rather than inferred from the debit.

## Bounded representative domain

The preceding quadratic package has1,361 accepted moving-wall events and
4,437 R caps. This source partitions the R-bearing events by
`(A family, P7 label, L profile)`. In every nonempty group it chooses the
event with a saved R–L product closest to−2, breaking ties by event index,
and retains **all** saved R caps at that event.

There are83 such groups and83 selected events, spanning18 A families and
516 R caps. Twenty-one events use P3,28 use P4 and34 use P5. The other five
of the predecessor's23 accepted A families have no saved R cap, so they
do not enter this continuation domain. Only one event per group is tested;
the selection is not a proof that other events behave identically.

For each event, K is solved once in the frame `(A,L,P)`. Its required
finite corners are A/L/K and P/L/K. The inherited seven-A2-profile
classification reads the actual A–L and L–P bonds only: A–P is a nonedge
of the intended patch and is not silently imposed as a right bond.
For K profile \(\pi\), its residual four-vector satisfies

\[
 K\mathbin{\cdot}n_i=-2\cos(\pi/m_i)-o_{\pi,i},
 \quad i\in\{A,L,P\},\qquad K^2=q_\pi.
\]

The native offsets \(o_{\pi,i}\) and residual square \(q_\pi\) restore
the S6 projection. Both norm branches are tested. Symbolic isolated high
labels are sampled at6,8,12,24; this finite sampling does not cover their
whole continuous or integer tails.

Each branch must pass its norm and prescribed-contact residuals, all
retained old native products, C/A/L/P contacts and the shared strict
center. Every accepted K is then joined with every saved R at that event.
The assembled native products and center are rechecked. Finite K–R
contact and both A/K/R and P/K/R corner classifications are necessary for
the intended completion. Dotted native-coherent joins are also scored by
the complete cone, so an unanticipated boundary improvement would not be
discarded merely for missing the proposed incidence.

## Results and complete boundary checks

| Stage | Count | Complete S6 finite rays | Exterior rays | Open ends |
|---|---:|---:|---:|---:|
| K frame equations / real branches |1,082 /2,164|—|—|—|
| Native-coherent K additions |586|34|6|10|
| Candidate K/R pairs |5,929|—|—|—|
| Native-coherent K/R additions |4,128|36|6|10|

The remaining1,801 K/R pairs fail the native Coxeter-contact test. There
are no accepted finite K–R contacts. All4,128 accepted products are dotted,
ranging numerically from approximately−19.91850716 to−2.01310803. Thus
none supplies the required A/K/R and P/K/R finite corners, and no whole
compact K quadrilateral is obtained.

Every accepted K and K/R addition undergoes complete numerical projective
ray/edge enumeration and a separate positive-definite pair/triple incidence
calculation. Their finite support sets and open-end counts agree. No saved
cone warning or unresolved case occurs. The full8,12 and C faces have,
respectively, exactly seven, four and three finite extreme rays and no
ideal or exterior extreme ray, in every accepted record.

In4,008 assembled K/R chambers the two B-region exterior supports are
`[A,R,K]` and `[P,R,K]`; K has two finite and two exterior rays. In the
other120 they are `[A,P,R]` and `[A,P,K]`; K has two finite and one
exterior ray. Neither case supplies the required four finite K rays. The
other exterior supports are the unchanged chain `[30,33,66]`,
`[30,59,66]`, `[58,60,66]`, `[59,60,66]`. The extra finite vertices
therefore move the unresolved B-region boundary without reducing its open
ends.

The saved K profiles are0:486,1:53,2:6,4:40,5:1. Among rejected K branches,
341 have old7 as their sole native blocker and six have old30 as their
sole native blocker. These are diagnostic counts, not a proof that changing
either old wall would permit a completed chamber. The failed finite join
has moved from L–R to K–R; further endpoint stacking alone has not supplied
closure in this selected domain.

## Scope and verification

The exact disk counts and old-incidence debit are distinct from the
floating-point normal, root and cone calculations. Numerical tolerances,
near-singular frames and geometric helper implementations are not
independently interval-certified here. The recorded absence of unresolved
cases is specific to this run. Unselected moving events, other high
labels, other A families/incidences and continuous deformations remain
open. No new exact ambient normal certificate, ambient finite-vertex
census, complete ambient cone or finite-volume conclusion is asserted.

**FINAL/FROZEN.** Generation and a fresh saved replay both terminated
successfully, with Python syntax, nine direct source/data bindings and
three report links checked. A separate saved-data walk importing no
producer checked the exact disk, representative-selection rule,4,714
accepted records,39,288 saved native products, all complete-cone/PD support
agreements, open-end sums and whole8/12/C ray sets. This was not a numerical
normal, root or cone reconstruction.

Root's independent full source/data read passed within the numerical
scope: representative selection, one K bank per event, the true A/P
nonedge, seven-profile finite-corner domain, native and center gates, and
complete K-alone/KR cone and PD-support checks including dotted joins. It
did not independently reconstruct physical coordinates, solve roots or
enumerate the cones. No new exact chamber is used.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_cu_release47_interior_cap_20261017.py
OPENBLAS_NUM_THREADS=1 python scripts/fresh_cu_release47_interior_cap_20261017.py --verify
```

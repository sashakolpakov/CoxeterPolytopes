# Flip the B-region central diagonal to L–R

The changed diagonal admits22 native L–R events, but no final Q cap in the sampled domain. Of2,024 Q branches,440 pass every other retained wall and cap and fail only at R: their native Q–R product is strictly positive, between approximately+0.16842728 and+1.14262520. The obstruction has moved from a dotted closing contact to an obtuse one. Every accepted event chamber still has ten open ends. No improved exact control or finite-volume example is produced.

The [source](../scripts/fresh_b_region_lr_diagonal_20261017.py) and [saved records](fresh_b_region_lr_diagonal_20261017.json) are a new successor of the frozen [P–Q diagonal experiment](fresh_b_region_dotted7_four_cap_20261017.md). No frozen file or root entrypoint is changed. This is a bounded numerical experiment with exact old-incidence bookkeeping and an exact abstract disk check; it does not certify new exact normals.

## Changed incidence, with the same deletion debit

Only old wall47 is removed. All other62 old walls, including8 and the old66 component, are retained. The exact debit is the same two ambient finite vertices

```
S6 + 7,12,47       S6 + 12,47,52
S6 = {0,1,2,3,5,6}
```

and one touched compact two-face certificate, supported by S6+12. There is no touched higher-dimensional compact-face certificate and no lost finite support containing8. “Touched” means that a removed old vertex belonged to that certificate, not that the support face has necessarily disappeared. The source recomputes the debit from all199 exact old finite supports and compares it to the frozen predecessor. It does not erase this cost when new finite corners appear. Wall8's deletion-side compactness certificate survives; no exact preservation claim for the newly sampled walls follows from numerical angles alone.

Keep the four boundary duties P:7/12, L:12/52, Q:52/30, R:30/7. Replace the two old central triangles LPQ/PQR by LPR/LQR. The complete new list is

```
7 12 P       12 52 L      52 30 Q      30 7 R
12 P L       52 L Q       30 Q R       7 R P
L P R        L Q R
```

The boundary is7–12–52–30–7. The new diagonal is L–R; P–Q is absent and may be dotted. With labels L67,P68,Q69,R70, this remains a disk with8 vertices,17 edges and10 triangles. All interior edges occur twice, boundary edges once, and all links are connected intervals or cycles. The cap links now have sizes L5,P4,Q4,R5. A separate peer topology read checked this literal flip only, with no geometry claim.

## Enforce the new native contact before constructing Q

L(u) has full native products L·12=L·52=0 and L·7=−u, with u>2. Its intrinsic norm branches use the frozen cofactor-oriented conic formula in frame(12,52,7), including all white-profile offsets. Both signs and all seven inherited A2 profiles are sampled.

For a given L, P is reconstructed from(7,12,L), imposing only the actual corner duties7/12/P and12/L/P. R is then reconstructed from(7,30,P), imposing7/30/R and7/P/R. During this screening step, R is tested against every retained old wall and P, but its product with L is deliberately left unresolved. That is the new moving equation:

\[
 L(u)\cdot R(u)=-2\cos(\pi/m).
\]

Only targets for which L/P/R is finite under the inherited corner classifier are considered. After solving a sign-bracket event, all native products are restored, including L–R, and all three caps must satisfy the same strict-center inequality. The event receives a complete numerical S6 cone computation and an independent PD incidence check; both finite support sets and open-end counts must agree. Thus an R passing the preliminary screen is not silently accepted with a nonnative L contact.

Q is constructed last from(30,52,L), requiring30/52/Q and52/L/Q. It must be native against all62 retained walls and L,P,R. P–Q is tested as a legitimate native relation but is not required to be finite. The two remaining finite duties30/Q/R and L/Q/R, and exact whole-cap face ray sets, are checked on any surviving Q. No Q reaches those completion checks in this run.

This construction does not recycle the old168 P–Q event prefixes. It solves a different contact equation and permits a different absent diagonal.

## Sampled domain and complete event boundary

The parameter grid is the predecessor's30 values:21 evenly spaced from2.01 to3, nine from3 to4 with the duplicate3 removed, and6. Seven profiles and both signs give420 L samples;124 pass the native and strict-center gates. Low finite labels2 through5 are included, while symbolically isolated high slots and the moving L–R target are sampled at6,8,12,24 in addition to the low labels. No claim is made for the omitted high-label tail.

There are9,292 P evaluations, of which2,974 pass. They generate19,648 R screening evaluations;6,266 pass the old-wall/P/center tests without imposing L. Only recipes native at both adjacent grid endpoints can create a sign bracket. All22 retained brackets yield rechecked native L–R events. Every accepted L and R has profile0; R is right to7 and30. The moving parameter ranges from approximately2.015370045 to2.356807769.

| Native L–R label | Accepted events |
|---:|---:|
|4|1|
|5|9|
|6|3|
|8|3|
|12|3|
|24|3|

| Stage | Accepted records | S6 finite rays | Exterior rays | Open ends |
|---|---:|---:|---:|---:|
|Exact63-wall control, numerical comparison|1|26|6|10|
|Delete47|1|24|6|10|
|Add sampled L|124|25|7|11|
|Native L/P/R event|22|30|6|10|
|Add Q|0|—|—|—|

Five distinct boundary ledgers are saved. One is the inherited helper's unrelated delete58/66 comparison; it is not the chosen deletion. All22 event chambers have the same exterior supports:

```
30,52,R       52,L,R
30,33,66      30,59,66      58,60,66      59,60,66
```

The new B-region has four one-open edges, and the unchanged66 region has six. These event chambers have65 nominal ambient walls, but30 is an S6 finite-ray count, not a new ambient finite-vertex census. The full boundary still witnesses failure of the intended closure at numerical scope.

The final Q banks contain1,012 norm equations and2,024 norm branches, all rejected by native contacts. The inherited near-rejection records are retained with their event IDs. They contain440 R-only failures,27 P-and-R failures, and106 old7-and-R failures. Every one of these573 records has a negative common-center product and a strictly positive Q–R product. Thus the R-only failures are not discarded because of a wrong center sign or a near-parallel finite-label tolerance. Their native profiles0,1,4 occur301,82,57 times respectively. These are sampled, numerical sign separations, not interval proofs across the parameter family.

## Scope and verification

No numerical unresolved flag occurred, but this does not certify complete root coverage. Sign-bracket preselection can miss tangencies or thin/disconnected native windows, and the grid omits parameters outside2.01..6. Higher isolated labels remain unsampled. The conic branch and finite-corner machinery are inherited inputs; no continuous or all-label exclusion is claimed. No Q cap is accepted, no new normal is exactified, and no full ambient Gram/cone certificate or finite-volume proof is supplied. The original63-wall/199-finite/ten-end exact control remains unchanged and infinite-volume.

Initial generation3972 completed terminal exit0/PASS. Final generation97850 also completed terminal exit0/PASS after adding only serialization of already collected near-rejection values with event IDs; the search domain and results were unchanged. Fresh saved replay39376 completed terminal exit0/PASS, reproducing the complete saved object.

A separate [saved-record checker](../scripts/fresh_b_region_lr_diagonal_20261017_saved_check.py), with no producer or numerical imports, completed terminal exit0/PASS (chunk9065de). It rehashes all14 direct bindings, reconstructs the two removed supports from the199-vertex census, compares the frozen full deletion debit, checks the flipped disk and links, all five saved cone/PD support agreements,22 events,2,024 rejected Q branches and573 near-rejection records. It confirms the positive products of the440 R-only failures from saved values; it does not recompute numerical normals, native arithmetic or cones. Both Python ASTs and four local links pass.

Exacthub independently read the full producer, report, saved checker and inherited conic-root/native/bank/cone interfaces. The review passed within the stated bounded numerical scope: R omits only L during screening and restores every native product at the event; recipe caching preserves the fixed profiles and duties; every accepted cone has an independent PD support/open-end comparison; the flipped disk, wall47 debit, final Q duties and440 positive Q–R failures match the report. This was a source/mathematical review, without independent numerical generation, cone reconstruction or saved-arithmetic replay. Parameter-window and high-label coverage remain explicitly open.

The source, saved data, report and saved checker are FINAL/FROZEN. No live job or additional bank remains.

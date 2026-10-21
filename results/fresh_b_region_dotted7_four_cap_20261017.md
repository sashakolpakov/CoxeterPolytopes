# Separate the 12–52 cap from the old 7-edge duty

Deleting only wall47 and allowing the new 12–52 cap L to be dotted to7 produces168 coherent moving P–Q events and258 further native R additions. Every accepted R–Q product is dotted. The complete numerical S6 cones retain ten open ends, so this experiment gives no improvement over the exact63-wall/199-finite/ten-open-end control and no finite-volume example.

The [source](../scripts/fresh_b_region_dotted7_four_cap_20261017.py) and [saved records](fresh_b_region_dotted7_four_cap_20261017.json) leave all frozen predecessors unchanged. The old incidence debit and abstract disk are exact combinatorial checks. New normal reconstruction, contact events and complete cone calculations are numerical. No sampled normal is promoted to an exact wall.

## Changed duty and exact deletion cost

The previous [released-prism search](fresh_release55_59_61_prism_replacement_20261017.md) left the two old exterior supports(7,30,47) and(30,47,52) untouched. Here47 is removed directly, while all other62 old walls remain. In particular, the old55/58/59/61/66 region and compact wall8 are retained. The old P-type replacement forced the new cap to own both7/12 and12/52, with a finite contact to7. This new cap owns only12/52 and instead satisfies

\[
 L\cdot12=L\cdot52=0,\qquad L\cdot7=-u,\qquad u>2.
\]

These are full native products. There is no C(u), old A prefix, or inherited requirement that L meet7. The current conic branch has a genuinely variable dotted native product.

The exact debit is reconstructed from all199 certified ambient finite supports of the original control. Removing47 loses exactly

\[
S_6\cup\{7,12,47\},\qquad S_6\cup\{12,47,52\},
\quad S_6=\{0,1,2,3,5,6\}.
\]

Exactly one compact old two-face certificate is touched: the quadrilateral with support S6+12. No higher-dimensional compact face certificate and no finite support containing8 is touched. “Touched” records a removed certified vertex; it does not claim that the corresponding support face is absent in the enlarged chamber. No deletion cost is suppressed when the new caps are counted.

| Compared deletion | Lost old finite vertices | Touched compact two-/three-face certificates | Lost supports containing8 |
|---|---:|---:|---:|
|47|2|1 / 0|0|
|30|13|10 / 1|0|
|52|5|4 / 0|2|
|30,47|15|11 / 1|0|
|47,53|2|1 / 0|0|

These comparison rows are incidence debits only; no completed search or cone certificate is asserted for the unchosen deletions. Retaining all wall8 supports preserves its deletion-side compactness certificate. Since the new normals remain numerical, this statement is not an exact ambient preservation theorem for the newly sampled chambers.

The chosen deletion cone has24 finite S6 rays, six exterior rays and ten open ends, compared with26/6/10 for the original control. Its two B-region exterior supports are(7,12,52) and(7,30,52). The four other exterior supports are(30,33,66),(30,59,66),(58,60,66),(59,60,66).

## Conic reconstruction and the new disk

For each of the seven inherited A2 profiles π, let qπ be its intrinsic projected square and oπ(i) its white projection product with old wall i. In the intrinsic Lorentz4 complement, set F=(12,52,7), G=FMFᵀ and

\[
r_\pi(u)=(-o_\pi(12),-o_\pi(52),-u-o_\pi(7)),\qquad
y_0=r_\pi(u)G^{-1}F.
\]

For an orthogonal kernel vector k and κ=k², the two possible norm branches are

\[
 y=y_0\pm\sqrt{\frac{q_\pi-y_0^2}{\kappa}}\,k.
\]

The producer uses ordered cofactors of FM for a continuous kernel orientation, rather than a sign-ambiguous SVD kernel. This formula describes the search equations; its physical coefficients are evaluated numerically here. Singular, null, tangent and residual-failure thresholds are explicitly unresolved rather than accepted.

The proposed replacement disk has boundary7–12–52–30–7 and four new vertices L67,P68,Q69,R70. Its ten triangles are

```
7 12 P       12 52 L      52 30 Q      30 7 R
12 P L       52 L Q       30 Q R       7 R P
L P Q        P Q R
```

It has8 vertices,17 edges and10 triangles. Every boundary edge occurs once, every interior edge twice, and each vertex link is a connected interval or cycle as appropriate. The interior diagonal is P–Q. In a realized disk, L/R would have whole quadrilateral facets and P/Q whole pentagonal facets. The source and a separate peer topology read both check these counts and links. No disk realization follows from this topology alone.

P is solved from(7,12,L), imposing the two finite corner duties7/12/P and12/L/P. Q is solved from(30,52,L), imposing30/52/Q and52/L/Q. The outer-anchor pairs7/L and30/L are not required to be edges. The moving parameter u is then chosen by the native equation

\[
P(u)\cdot Q(u)=-2\cos(\pi/m),
\]

with the L/P/Q finite-corner classifier applied before searching. R is subsequently solved from(7,30,P), with required finite corners7/30/R and7/P/R. The remaining30/Q/R and P/Q/R corners are checked on every accepted R; they are not assumed from its frame.

## Bounded search and complete numerical boundary results

The search samples30 values of u:21 evenly spaced values from2.01 to3, nine from3 to4 with the duplicate3 removed, and6. All seven L profiles and both conic signs give420 samples, of which124 pass the retained native walls and shared strict center. All endpoint recipes use the seven inherited profiles and low labels2 through5; isolated high slots are sampled at6,8,12,24. The P–Q event labels are2,3,4,5,6,8,12,24.

There are14,696 endpoint evaluations and3,521 native endpoint records. A sign bracket is considered only when the same P and Q recipes are native at both adjacent grid endpoints. All168 retained brackets yield an accepted rechecked event; their u values lie between approximately2.014592713 and3.977209131. L profiles0,1,4 supply142,16,10 events respectively. Every event restores all profile offsets, tests all62 retained old walls and every new mutual product, satisfies the common strict-center gate, and passes a complete numerical S6 cone enumeration against an independent positive-definite incidence computation. Finite support sets and open-end counts must agree.

The R banks contain850 norm equations and1,700 branches, with258 accepted R additions. All accepted R profiles are zero and their prescribed(7,30,P) labels are(2,2,m), where m runs over the eight sampled labels. Each accepted R chamber receives the same full native, center, cone and PD-support tests. Intermediate individual P/Q samples are native screening records; the full cone claim applies to joint events and accepted R chambers, not to every endpoint sample.

| Stage | Accepted records | S6 finite rays | Exterior rays | Open ends |
|---|---:|---:|---:|---:|
|Original63-wall control|1|26|6|10|
|Delete47|1|24|6|10|
|Add L at a sampled parameter|124|25|7|11|
|Add P and Q at a native event|168|30|6|10|
|Add R|258|32|6|10|

Eight distinct complete boundary ledgers are saved. One is an inherited helper's delete58/66 control record; it is bookkeeping only and is not part of the chosen delete47 search. With L/P/Q/R present there are66 nominal ambient walls, but32 is an S6 finite-ray count, not an ambient finite-vertex count.

Every R–Q native product lies in the observed numerical range

\[
-9.037647064\ \le R\cdot Q\le\ -2.004692716.
\]

Thus this failure is separated from the near-parallel finite-label tolerance issue of older scouts. Exactly the two proposed finite corners30/Q/R and P/Q/R fail in all258 records. L has precisely its four expected finite rays and no other rays on its face; the P/Q/R whole-face checks fail. No disk is completed.

The old four-ray66 component persists. The other two exterior supports have two observed patterns:

| Records | New B-region exterior supports |
|---:|---|
|174|(30,Q,R), (P,Q,R)|
|84|(30,P,Q), (30,P,R)|

The second pattern matters: a failed prescribed corner is not automatically an actual exterior ray. The saved whole cone determines the boundary. Both patterns retain four B-region one-open edges and six one-open edges in the old66 region, for ten open ends overall. The original exterior supports involving47 have changed, but their open boundary has not closed.

## Limits and verification

There are no unresolved flags in the sampled producer output. This does not certify complete continuous root coverage. The grid and endpoint-native bracket preselection can miss tangencies, disconnected or thin admissible windows, and roots at which a recipe is not native at a bracket endpoint. Other finite labels and parameters outside2.01..6 remain open. The finite corner classifier is an inherited input, while its isolated high slots are only sampled here. No all-label or all-parameter obstruction is asserted.

This changed incidence successfully separates L from7 and makes the L quadrilateral finite, but the sampled closure still transfers the gap to R–Q. A further construction would need to change that remaining duty or parameter geometry; merely counting the compact L face is not a chamber closure. The exact63-wall/199-finite/ten-end control remains unchanged and infinite-volume. No new full ambient finite census, exact escape certificate, finite-volume certificate or dimension example is claimed.

Generation17371 and fresh saved replay34961 both completed terminal exit0/PASS. The replay recomputes the bounded numerical experiment and compares the entire saved object. A separate [saved-data checker](../scripts/fresh_b_region_dotted7_four_cap_20261017_saved_check.py), with no producer or numerical-geometry imports, completed terminal exit0/PASS (chunk025338). It rehashes all11 direct bindings, independently rebuilds the five deletion debits from199 old finite supports, checks the disk links, all eight saved cone/PD incidence agreements,168 events,258 completions, their two exterior patterns and the saved dotted R–Q range. It does not rerun normal reconstruction, native numerical arithmetic or cone enumeration.

The peer review covers the abstract ten-triangle topology only. The parent completed a full narrow source/report and saved-checker read and found no issue in the exact deletion/disk distinction, cofactor branch orientation, endpoint-native sign brackets, joint native/common-center/cone gates, or bounded coverage statements. A separate parent saved-record query (chunk aa1881) confirmed the literal174/84 exterior support sets and persistent66 component. This review did not recompute numerical normals or cones and did not prove a continuous interval exclusion. All11 direct hashes, both Python ASTs and four local links pass. Producer, saved data, durable checker and report are FINAL/FROZEN. There are no live jobs and no frozen predecessor or root entrypoint was edited.

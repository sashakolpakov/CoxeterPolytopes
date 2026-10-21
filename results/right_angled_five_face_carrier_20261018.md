# A new five-face incidence constraint from the ideal RA4 equality theorem

18 October2026. **Final universal corollary with exact identity checks and two independent mathematical reviews.** No new right-angled dimension exclusion or example is claimed. No arithmetic field restriction is imposed. This gives a concrete research constraint for the open ideal dimensions5/6 and mixed dimensions9–12, using the frozen [second-value theorem](ideal_ra4_second_volume_20261018.md) and [equality classification](ideal_ra4_second_volume_equality_20261018.md).

## A necessary expensive facet in every ideal RA5 polytope

For an entirely ideal right-angled hyperbolic4-polytope Q, write

\[
 \Delta(Q)=f_0(Q)-f_3(Q),\qquad
 \kappa(Q)=3\Delta(Q)-(f_3(Q)-24)=3f_0(Q)-4f_3(Q)+24.
\]

The24-cell has $\Delta=0$, $\kappa=0$. The reflected double has $\Delta=4$, $\kappa=-2$. The gap and equality theorems imply that any other Q has integer $\Delta\ge5$ and volume at least $3\pi^2$.

Let P be an entirely ideal right-angled hyperbolic5-polytope and put $v=f_0$, $e=f_1$, $r=f_2$, $g=f_3$, $f=f_4$. Its cusp link is a rectangular4-box. It has16 vertices,24 two-faces and8 facets, so

\[
 e=8v,\qquad \sum_{R\in\Omega_3(P)}f_0(R)=24v,
 \qquad \sum_{Q\in\Omega_4(P)}f_0(Q)=8v.
\]

Every ideal right-angled3-face R satisfies $f_2(R)=f_0(R)+2$. Every2-face of P is in exactly three3-faces, because its finite relative-interior normal diagram consists of three orthogonal walls. Hence $3r=24v+2g$. Euler's equation for the boundary4-sphere is $v-e+r-g+f=2$. These two identities give

\[
 g=3v+3f-6.
\]

Set $D=\sum_Q\Delta(Q)$ and $A=\sum_Q(f_3(Q)-24)$. Each3-face belongs to two4-facets, so

\[
 D=8v-2g=2v-6f+12,\qquad
 A=2g-24f=3D-48.
\]

Therefore the following identity holds for **every** entirely ideal RA5 polytope:

\[
 \boxed{\sum_{Q\in\Omega_4(P)}\kappa(Q)=48.} \tag{1}
\]

At least one4-facet must have positive kappa. Since the24-cell and its double have nonpositive kappa, this facet is neither. It has volume at least $3\pi^2$, and also satisfies the stricter incidence condition

\[
 f_3(Q)\le3\Delta(Q)+23,\qquad f_0(Q)\le4\Delta(Q)+23. \tag{2}
\]

For example, the first allowed delta5 branch would require at most38 facets and43 vertices. This is a necessary candidate domain, not an assertion that such a facet exists. Merely constructing ideal RA5 using24-cell and double24-cell facets is impossible, irrespective of how many are used or how they are glued: their kappa sum would be nonpositive. This strengthens the earlier all24-cell Euler obstruction.

For an entirely ideal RA6 polytope, each5-facet satisfies (1). Every4-face belongs to two5-facets, so

\[
 \sum_{Q\in\Omega_4(P)}\kappa(Q)=24 f_5(P).
\]

Each5-facet must meet at least one positive-kappa4-face; consequently there are at least $\lceil f_5(P)/2\rceil$ such4-faces. These are incidence requirements, not a contradiction in dimension6.

## A shared constraint for mixed dimensions

Let P be any finite-volume right-angled hyperbolic n-polytope, n≥5, with $v_{\rm fin}$ finite vertices. Let I5 count its entirely ideal5-faces. A5-face not entirely ideal contains at least one finite vertex. A finite vertex is simple and belongs to exactly $\binom n5$ five-faces. Thus

\[
 f_5(P)-I_5\le\binom n5 v_{\rm fin}. \tag{3}
\]

Every4-face has finite relative-interior points and an orthogonal normal diagram of rank n−4, whether or not its vertices are ideal. It belongs to exactly n−4 five-faces. Apply (1) to each ideal5-face, discard its nonpositive kappa terms, and sum. Only entirely ideal4-faces can occur inside these ideal5-faces. This gives

\[
 48I_5\le(n-4)\sum_{Q\in\Omega_4(P),\ Q\ \mathrm{ideal}}\max(\kappa(Q),0).
\]

Combining with (3),

\[
 \boxed{48\left[f_5(P)-\binom n5 v_{\rm fin}\right]
 \le(n-4)\sum_{Q\ \mathrm{ideal4}}\max(\kappa(Q),0).} \tag{4}
\]

The unweighted version is $f_5(P)\le\binom n5v_{\rm fin}+(n-4)h_4$, where h4 counts the positive-kappa entirely ideal4-faces: choose a finite vertex for every nonideal5-face and a positive-kappa4-face for every ideal5-face. The relevant incidence coefficients are:

|Ambient dimension|Five-faces at a finite vertex|Five-faces at a4-face|
|---:|---:|---:|
|9|126|5|
|10|252|6|
|11|462|7|
|12|792|8|

This couples the finite-vertex budget to a geometrically restricted family of ideal4-faces. A positive carrier must lie beyond both the24-cell and its unique second-volume double. It can be combined with weighted finite/ideal incidence inequalities, but no numerical dimension contradiction follows from (4) alone. A useful next task is a universal lower bound on $f_5-\binom n5v_{\rm fin}$ together with a geometric upper bound on the positive carrier sum, or a proof that all ideal RA4 facets have nonpositive kappa. Either would be substantially stronger than a bounded search for more examples.

## Exact and logical scope

The [producer](../scripts/fresh_right_angled_five_face_carrier_20261018.py) checks the exact linear combinations of the cusp/Euler identities, the kappa values of the two classified4-polytopes, and all incidence coefficients for dimensions5–12. It does not certify the universal geometric premises merely by sampling values. Those premises are the face/cusp arguments above plus the frozen RA4 theorems. No recent preprint vertex bound is needed as an assumption of this corollary. The2023 Alexandrov weighted inequalities and2026 Egorov lower counts remain complementary research inputs, as recorded in the [literature audit](right_angled_dimension_literature_20261018.md).

Generation and a separate-process fresh replay both completed with exit0/PASS in terminal chunk `b132c2`. Root and construction independently read the complete producer and written mathematical argument: both PASS. Their scope includes the rectangular-box incidences, Euler elimination, positive-carrier deduction and the ideal6/mixed coface inequalities. Neither reviewer executed the script or reproved the frozen second-value/equality classification. The script checks exact identity arithmetic; the universal geometric proof is the written argument, not a formal kernel verification. This packet does not strengthen the published dimension upper bounds6 or12, and it does not assert the existence of any positive-kappa ideal RA4 polytope.

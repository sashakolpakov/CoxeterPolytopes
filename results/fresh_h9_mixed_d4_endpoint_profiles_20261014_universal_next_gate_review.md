# Independent review of the 71-edge endpoint gate

**PASS.** The [extension](fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_next_gate.md) correctly carries the common-gap theorem to all **71 open directions** of the frozen 14-wall chamber. Its scope is **2,221 finite spherical profiles and 92 continuous dihedral tails**, including **459 profiles and 21 tails** on the 21 directions at the three new finite vertices. No mathematical or coverage defect was found.

The [separate checker](../scripts/fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_next_gate_review.py) and [saved result](fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_next_gate_review.json) use fresh standard-library `Fraction` arithmetic in
\[
\mathbb Q\!\left(\sqrt5,\sqrt{(16+4\sqrt5)/11}\right).
\]
They import no producer or project arithmetic. Signs are decided by exact successive quadratic comparisons, without numerical tolerances. Source hashes bind the primary PY/JSON, the complete 14-wall cone, its profile census, and the [previous interval-engine review](fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate_review.md); recursive manifests also pass.

The independent reconstruction checks every outgoing frame
\[
A=-X^2>0,\qquad Y=Aq+(X\cdot q)X,
\]
including \(Y^2>0\), \(X\cdot Y=0\), future orientation of \(X\), orthogonality to all eight retained support walls, and the strictly negative outgoing pairing with the released wall. It matches all 71 directions to their original indices in the complete 135-direction queue. All 15 finite corners have exactly nine incident walls, positive-definite normalized Grams, and entrywise nonnegative inverse. The three new corners and their 21 directions are identified from their literal incident sets.

The positive sheet is exhaustive for a valid finite-endpoint addition. In a complete old finite-corner sector, write a future strict interior point and an acute added normal as
\[
z=bU-\sum d_i\omega_i,\qquad r=hU-\sum c_i\omega_i,
\]
with \(b>0\), \(d_i>0\), \(c_i\ge0\). The corner inverse is nonnegative. If \(h\le0\), then \(r\cdot z=-hb+c^T Jd\ge0\), so this normal cannot be an inward wall of a chamber with a strict future interior point in that sector. Hence \(h>0\). Forward intersection with the old edge selects the positive outgoing coefficient, giving
\[
r=p+hU+\sqrt{h^2+\rho}\,V.
\]
The certificate includes \(h=0\) as a harmless closed enlargement. This argument does not require the originally saved center itself to remain inside the enlarged chamber.

The checker separately reconstructs all normalized support Grams and inverses, verifies every finite profile cost and positive residual norm, and checks every tail is attached to an isolated \(A_1\) with inverse diagonal \(1/2\). It verifies the exact escape norm, future sign and all 14 strict wall signs. It also proves the scalar monotonicity conditions directly: for unnormalized pairings \(a=X\cdot w\), \(c=Y\cdot w\), when \(c>0\), the comparison
\[
a<0,\qquad a^2Y^2>c^2(-X^2)
\]
certifies \(a/\sqrt{-X^2}+c/\sqrt{Y^2}<0\). When \(c\le0\), the signs suffice. Thus the extension's finite uniform height bound and endpoint evaluation use the same justified inequalities as the reviewed base gate.

The interval engine, alternative cosine audit, treatment of \(k=2\), overlapping closed children, and complete high-label relaxation are inherited unchanged. Labels 2 through 15 are isolated points; the whole ray below \(-2\cos(\pi/16)\) retains every larger label, parallel contact and dotted contact. There is no finite label cutoff. Profile-list equality, permitted outside-wall indices, binary-tree count identities and all totals are independently checked. The full primary saved replay passed **11,093 boxes / 6,703 terminal boxes**, maximum depth 9. The independent frame/coverage checker also passed generation and saved replay.

This is a separate exact frame, Gram, cost, sign and coverage reconstruction plus source/proof review. It is **not** a second implementation of every interval box decision or of the inherited complete spherical-profile classification. Those dependencies and the full shared-engine replay are stated explicitly.

The conclusion remains narrow: every finite-endpoint cap on these saved 71 directions that preserves the 14 Coxeter walls strictly misses the same rational spacelike point. A finite coherent collection of these caps therefore retains a strict gap and infinite volume. Arbitrary caps, affine endpoint profiles, later supports, vertex-free lines and changed cores are outside this theorem.

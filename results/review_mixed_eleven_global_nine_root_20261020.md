# Root review of the global nine-direction proof

20 October 2026. **PASS.** I read the complete analytic proof and its
exact-polynomial checker and independently checked every inequality.
This is a scalar theorem; geometric composition is a separate step.

Strictness correction found during TeX translation: the pinned source's
displayed clipping chain uses strict signs after multiplying a strict
constant comparison by d^4. At d=0 those signs must be weak; the
valid universal statement is (lambda-r^-2)_+<=3d^4, which is exactly
the statement subsequently used. The theorem and certificate identities
are unchanged. The final TeX uses weak signs there and states the
reciprocal coefficient is at most r^-2+4r^-1 before its strict global
constant bound. Preserve the pinned source bytes with this explicit
correction instead of silently changing its certificate binding.

Source `right_angled_mixed_eleven_nine_direction_global_proof_20261020.md`,
SHA256 `bf288985ade5acf93e13faafb25544b0642b38f13e422ad38a3f1223e9706d15`;
checker `scripts/check_mixed_eleven_nine_direction_global_20261020.py`,
SHA256 `357070546d69aaa634c01292fb15c174f6d85552e87cb42be6adbf4d2cdb0b62`.

The nine-coordinate sum gives delta in[0,12/11] and
V<=2delta/3-11delta^2/18. The centered moments V/4,
-3T3/28 and5V^2/28-3T4/28 are correct. The sharp cubic bound
follows from summing the displayed nonnegative cubic polynomials;
its sign agrees with the negative triple-sum deviation. Its deliberately
loose consequence (3/7)T3<=V/6 is valid because V<=2/11 and sqrt11>3.

The clipping correction is genuinely global. For b=10/sqrt37,
mu<=15/11 and b>18/11 give b-mu>3/11. Tangency followed by
(d-a)_+<=27d^4/(256a^3) bounds the correction by3d^4. The latter
scalar inequality also has the exact nonnegative factorization checked
by the worker. No triples above the clipping threshold are omitted.

Since mu>=1 and r>11/20, the reciprocal Taylor coefficient is at
most1280/121. Combining it with clipping and Ed^4<=V^2/6 gives
exactly1643V^2/726. All terms increase with V, so replacing V
by x(4-11x)/2, x=delta/3, is in the required direction.

I checked the final identity and the derivative of
(4-11x)^2(1+x)^5. The derivative is
(4-11x)(1+x)^4(-2-77x), nonpositive on[0,4/11]. Its maximum
is16, and the remaining polynomial terms are bounded below by233/24.
The rational margin is exactly635/968>0. This proves the global
upper average63/100-delta/18 and the stated uniform equality case.

The complete checker passed generation and a separate read-only replay
under the shared lock and59-second/1GiB limits. Ledgers:
`mixed_eleven_nine_direction_global_generation_20261020.json` and
`mixed_eleven_nine_direction_global_replay_20261020.json`. Each took
less than0.1 seconds and55MB. Its scalar polynomial construction,
moment coefficients, Taylor identity and clipping factorization were
all inspected in source. There was no search or numerical optimizer.

The earlier complete interval proof was also written and independently
replayed, and is preserved as separate supporting evidence. The global
proof is preferable for the manuscript because its universal inequality
has a short human derivation. Neither proof alone implies a dimension
bound without the other subset sizes and actual opposition transport.

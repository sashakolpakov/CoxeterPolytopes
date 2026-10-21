# D18 cubic dotted two-cycle exclusion: full outside-author review

21 October2026. **FULL MANUAL MATHEMATICAL PASS; no correction requested.** Reviewed the complete source `results/arithmetic18_whole_d18_cubic_dotted_head_two_cycle_exclusion_20261021.md`, SHA256 `f1a40f93284a52519aac2415aa2fb15e5f6b099ac2045e5a38f7c457b75f8269`. I did not author this cubic proof or the new two-cycle lemma. I contributed inherited actual-wall/prefix foundations, reviewed the D18 head-core predecessor and separately reviewed root's quartic use of the same new lemma. Root and the composition author supplied the current norm/local-field argument. I independently rederived the complete cubic chain and hand tables; no previous PASS label was used as its proof. No worker, symbolic/enumeration computation, import, Lean or TeX job, source edit or status promotion was performed.

The two explicitly named source pins match current streamed bytes:

| Input under `results/` | SHA256 | Exact used boundary |
|---|---|---|
| `arithmetic18_whole_d18_dotted_head_degree_four_and_isolated_forks_20261021.md` | `0590f2b9d8747f6edbe0cafd74f2909c54c87f62e5ebe7d52081abbdfe35fa16` | Actual cubic U, both isolated dotted forks, full conjugate strips, N(V)<=14, actual K=Q(U)=Q(V), raw time line and parity cycles |
| `arithmetic18_whole_d18_long_end_isolation_and_head_field_reduction_20261021.md` | `6ef34bb417647b846c9b26b7c35fb8018e27b999c0395396c33c6d5ba691b3ac` | Original actual D18 geometry, isolated head and actual-field normalization |

These inputs distinguish unit actual cap normals from field-valued raw roots. The new proof uses raw square2U/2V and their stated own products, not norm2 for both caps. It assumes no power basis for U and no finite global field list. Elementary algebraic-integer norm, local valuation, coprime Hensel lifting and discriminant/index facts are the explicit arithmetic foundations.

## 1. Actual cycles and all norm possibilities

The chosen head and fork walls are distinct by their sole old contacts. The actual closed old-chain cycle gives integral g=(r0,rf); its product is±g after cancellation of root rescalings. Their two-cycle gives the separately integral L=g²/(UV), with no missing scalar factor. The spatial product is UV/2, so T=UV/2−g is integral from the inherited parity conclusion and satisfies `T²=U(U−2)V(9V−4)/2`. Expanding and dividing by UV gives `L+T=19UV/4−2U−9V+4`. Thus19UV/4 is integral, and `3*(19UV/4)−14UV=UV/4` is integral. No unit condition on U is inserted. The two nonzero integral quantities UV/4 and(2−U)V/2 therefore give64|A*n and8|B*n at degree3.

I independently checked the sharper head product bound. For s=U−1 and the magnitudes r1,r2 of its other conjugates, the nonzero integral norm gives s*r1*r2>=1. Since both r_i lie strictly between0 and1, the identity for A*B and `r1²+r2²>=2r1r2` yield the displayed bound by `(s²−1)(1−1/s)²`. Both factors increase for s>1, so the strict endpoint is `298144/5625<54`. The separate bounds A<=37,B<=29 and1<=n<=14 use EVERY nonphysical embedding.

The2-adic cases are complete. Odd n contradicts A<=37; v2(n)=1 forces A>=32,B>=4, contradicting A*B<54. Thus only n4/8/12 remain. At n4/12 the sole norm pair is(16,2). At n8, A is a positive multiple of8 and B is even because p_U(2) and p_U(0) agree modulo2. The complete union is exactly(8,2),(8,4),(8,6),(16,2),(24,2). The formula `S=2T+(A−B)/2−4`,5<=T<=13, is correct, including the sign p_U(2)=−B.

For all pairs except(8,4), S is odd. Modulo2, p_U=x*(unit-constant quadratic), a COPRIME degree1/degree2 factorization. Hensel lifting therefore splits the local Q2 algebra into those degrees; on the degree2 component U and its inverse are integral. UV/4 then makes V/4 integral there, forcing its norm contribution divisible by4²=16. The other component contributes an integral norm, so16|n, impossible. This argument remains valid when Z[U] has even index or the quadratic reduction has a repeated internal factor: only the coprime outer splitting and invertibility of U are used. The remaining values are A8,B4,n8.

## 2. The complete head polynomial reduction and T11 local obstruction

Shifting X=y+1 gives `q_t=y³−ty²+y+t−6`, t=T−3 in2..10. I checked each interval elimination: for t<=5, q_t<=q5<0 on(-1,1), using the displayed negative completed square for y>=0; t6 has rational root0; for t>=9, the value at25/3 is at most−944/27, impossible when all three roots are smaller. Thus only T10/11 remain, without a hidden unexamined coefficient table.

At T11, modulo2 is x²(x+1). The coprime split gives one unit-U degree1 component and one even-U degree2 component. Their U and U−2 norm valuations on the latter are3 and2. An unramified quadratic component cannot give odd U-norm valuation. A split component would have U valuations1/2, so U−2 valuations>=2/exactly1, total>=3 rather than2. Thus it must be ramified quadratic, with v(2)=2, v(U)=3,v(U−2)=2 and v(H0²)=5.

On the degree1 unit component UV/4 forces v(V)>=2. Since total v2(n)=3, the ramified component has v(V)<=1. Therefore v(9V−4)=v(V), and v(Hf²)=2v(V)−2 is even. The time-square ratio Hf²/H0² has odd valuation there, so cannot be a square in the ACTUAL K. This excludes T11 for every possible placement of the actual V conjugates; no choice of a convenient local embedding is made. The raw fork division by2 is retained in that valuation.

For T10, substituting delta=2/(U−2) gives `delta³+5delta²+4delta−2=0`, which I re-expanded. Delta is integral and generates the SAME cubic K. Its discriminant is316=4*79. The order index can only be1 or2; index2 would give field discriminant79, impossible since a number-field discriminant is0 or1 modulo4. Thus O_K=Z[delta] and Disc(K)=316. This maximal-order argument occurs ONLY after the earlier local reasoning; it does not assume a power basis for U retroactively.

## 3. Full seven-row fork norm obstruction

The actual V generates K and has norm8. With physical v and other positive conjugates a,b<4/9, ab=8/v implies v>81/2 and hence integral trace41<=R<=75. For its polynomial X³−R X²+S_V X−8, direct substitution at4/9 gives

`m=N(9V−4)=5768+144R−324S_V>0`.

The actual raw time-norm ratio is EXACTLY m/32: N(Hf²)=n*m/8, while N(H0²)=A*B=32. This necessary rational-square condition implies m=2k²; m divisible by4 then forces k even and m=8l². I checked the complete upper estimate `(4−9a)(4−9b)<=[4−9sqrt(ab)]²`, its positive-square direction, sqrt(ab)>3sqrt2/13 and sqrt2>24/17. Together with9v−4<672 they give `m<37427712/48841<800`. Thus1<=l<=9. Modulo9 the exact equation gives l²=1, leaving ONLY l1/8 in that interval.

For l1, `9S_V=160+4R`, R=5 modulo9; for l8, `9S_V=146+4R`, R=4 modulo9. The trace interval therefore gives exactly the seven displayed rows. I hand-recomputed every discriminant from `R²*S_V²−4S_V³−32R³−1728+144R*S_V`:

| l,R,S_V | Recomputed discriminant |
|---|---:|
|1,41,36|−2704|
|1,50,40|30272|
|1,59,44|198448|
|1,68,48|617792|
|8,49,38|−250812|
|8,58,42|−256784|
|8,67,46|−72956|

Negative rows contradict total reality. Since the remaining actual V is a primitive integral generator of K, its order discriminant must be316 times an INTEGER square. The values30272 and617792 have remainders15 and12 modulo79. The last198448 equals316*628, and25²<628<26². Thus every row is impossible. The finite table has explicit domain, exact formulas, termination and exhaustion; no numerical or enumeration job is hidden in this review.

The proof excludes the ENTIRE actual cubic dotted-head case at WHOLE arithmetic D18. It uses one fork, with no required distinctness of the two fork opposites, no extra component/fork energy and no formal-Gram completion promoted to a polytope. The source itself does not settle quadratic or quartic cases, wholeD18 or an arithmetic17 ceiling; any separate quartic composition remains separate evidence. No new example, manuscript result, worker output or geometric Lean verification follows from this review.

# A restricted short-eight tail exchange is excluded by the native gap

21 October 2026. New manual supporting lemma, pending independent review. No scientific, symbolic, enumeration, worker-import, Lean or TeX job ran. Current rules, manuscripts and accepted sources remain unchanged. This excludes the complete profile stated below, not all D25 short-eight faces or a whole dimension.

**Normalization correction.** An unreviewed working draft described a norm-two product of magnitude `sqrt(120/103)` as dotted. That was false: for norm-two roots the simple3 magnitude is1, the next finite magnitude is sqrt2, parallel magnitude is2 and dotted magnitudes exceed2. The proposed star therefore violated the native gap. Its exact rejected8725-byte source, SHA `8893782b81226be5816ac2390cb5d77e53825ed608eb96e18a58b1b3423b65c2`, is preserved in `source_archive/finite57_d25_tail_exchange_rejected_norm_two_dotted_draft_20261021.md.gz` (4058B, SHA `97acb4796782c2c4143fcf8c70adbd207d6ea37404b54ba5c10124ce8ee6cc96`; roundtrip verified). It earns no Coxeter-realizability, local-obstruction or example credit. The same exact rank identity instead proves the exclusion below.

The current comparison is `finite57_d25_short8_off_face_a16_e8_financing_20261021.md`, SHA `3223a63642d0ddf3425e8f5c111b1d207e8455bfb7db3ccb5e3155a15c5ce6d5`; its different Z-orthogonal family remains unchanged. The current weight source is `manuscripts/finite_volume57.tex`, SHA `c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0`.

## 1. Precisely the profile excluded

At an actual finite vertex take a WHOLE norm-two D25 with spine alpha0,...,alpha22 and short leaves alpha23,alpha24. Put a=alpha15, b=alpha23. All other old roots form an orthogonal finite forest Z. Suppose Z contains a WHOLE A15 path with endpoint z0. Let c,d be the actual opposite side roots of a quadrilateral with opposite pairs {a,c},{b,d}. Assume the following complete products:

* c has simple3 contact with alpha16 and z0, own product `c.a=-s`, and no other old contact;
* d has own product `d.b=-t` and NO supporting-white or other old contact;
* the adjacent products `c.b=d.a=c.d=0`.

Here s,t are unrestricted actual nonnegative native magnitudes, including0, all finite labels, parallel or dotted contacts. The conclusion is that NO such actual Coxeter profile exists. No arithmeticity is used.

This is the tail-exchange branch: after deleting a, c joins the right D9 to the other A15 tail, while the original alpha0,...,alpha14 becomes the unmarked A15. If c,b is another whole D25 short-leaf-eight corner, c has no contact with the original A15, and its remaining attached tail is disjoint, finite classification forces exactly this simple A15 tail. The isolated-white hypothesis on d is an additional restriction; it has NOT been established for every universal producer.

## 2. Full rank identity, with all root norms retained

Let p be the old vertex's future unit timelike vector. The complete old finite forest spans p-perp. Let omega_i denote its D25 duals and eta the dual of z0 in the separate A15. Their exact products are

`omega15^2=16`, `omega16^2=17`, `omega15.omega16=16`,

`omega23^2=25/4`, `omega15.omega23=8`, `omega16.omega23=17/2`, `eta^2=15/16`.

The two actual norm-two roots must therefore be

`c=A*p-(s*omega15+omega16+eta)`,

`d=B*p-t*omega23`,

where A,B>0 since they do not contain the old vertex. Normalization and `c.d=0` give

`A^2=16(s+1)^2-1/16`, `B^2=(25/4)t^2-2`,

`A*B=(8s+17/2)t`.

These are consequences of the actual one-dimensional time complement, not assumptions that a formal Gram is realized. Squaring the last equation and substituting the first two yields

`t^2=(2048s^2+4096s+2040)/(2304s^2+4096s+1751)`.  (1)

The denominator is strictly positive for every s>=0. Moreover

`t^2-1=(289-256s^2)/(2304s^2+4096s+1751)`.           (2)

If t=0, B^2=-2 is impossible. Thus the native gap requires t>=1. Likewise s is0,1, or at least sqrt2. Equation(2) immediately excludes s>=sqrt2, since `289-256*2<0`. At the only two remaining native values,

`s=0: t^2=120/103`,

`s=1: t^2=248/247`.

Both lie strictly between1 and2. Thus `1<t<sqrt2`, which is not an actual Coxeter product in norm-two normalization. This contradiction covers EVERY own-contact label in the stated profile, not merely an isolated or finite-own subcase.

## 3. Why the geometric working draft was tempting but invalid

At s=0 the nonnative formal roots have `A^2=255/16`, `B^2=544/103`, `t^2=120/103`. Eliminating the white forest gives the exact projected side Gram, ordered a,b,c,d,

`[ 25/144  -2/9   -8/9       0 ]`

`[ -2/9     4/9   -2/9      -t ]`

`[ -8/9    -2/9   25/144     0 ]`

`[  0       -t      0         2 ]`.

It really is a Lorentzian quadrilateral configuration: in square X^2+Y^2-T^2 its side normals can be

`a_bar=(sqrt(17/32),0,sqrt(103/288))`,

`c_bar=(-sqrt(17/32),0,sqrt(103/288))`,

`b_bar=(0,-sqrt(60/103),sqrt(128/(9*103)))`,

`d_bar=(0,sqrt2,0)`.

Its Klein rectangle is `|X/T|<=sqrt(103/153)`, `-sqrt(32/135)<=Y/T<=0`; all four rays are feasible and finite, with Lorentz squares `-50/153` and `-206/2295` at T=1. The full formal Gram has signature(57,1,1) when Z=A15+17A1 in ambient57. Its four corner types would be two `D25+A15+17A1` and two `A24+A15+18A1`. Under the CURRENT marked rule, the first two would pay the targeted short-leaf-eight units; the other marks are unlisted A1/A24 roles and would pay0. There is no shorter adjacent-black white path.

These correct Schur/rectangle and weight calculations do NOT rescue the example: its unprojected pair b,d has the forbidden product `-sqrt(120/103)`. Intrinsic projected face angles are not the test for ambient Coxeter labels. This illustrates a genuine native obstruction to an otherwise feasible Lorentzian/edge relaxation, not a counterexample to a payment theorem.

## 4. Uniform formula and what remains open

The rank identity is not a selected-row coincidence. For D_(p+10) with two exchanged A_p tails and the same right A8 white part, put L=p+1. The identical profile gives

`A^2=L(s+1)^2-1/L`,

`t^2=8*(L^2*(s+1)^2-1)/(L^2*(9s^2+16s+7)-2L-9)`,

`t^2-1=((L+1)^2-L^2*s^2)/(L^2*(9s^2+16s+7)-2L-9)`.

For every integer p>=2 the denominator is positive, and `1+1/L<=4/3<sqrt2`. Thus s>=sqrt2 forces t<1; s0 gives `8p/(7p-2)` and s1 gives `(16p+8)/(16p+7)`, both strictly between1 and2. The same restricted native-gap exclusion holds for all these tail lengths, in particular all large D ranks relevant here, whenever the stated actual profile occurs.

The universal short-eight problem still has unclassified branches: c may contact only one original large white block with other spherical attachments, the new tail need not be a pure A_p in a different output type, and d may have nonzero contacts with W or Z. Adjacent nonzero contacts and other current corner payments also require their own exact routing. The already financed3223 family is not changed, and its A16/E8 source is not inferred from this calculation. No omitted inert root is counted as a supporting white.

**Scope:** complete exclusion of the explicitly stated native tail-exchange/white-isolated-d profile; no complete universal D25 or D_r distance-eight replacement, finite-volume existence result, coupon deletion, finite maximum or global56 claim. Existing formal budget profiles such as `2D17+D23` remain only counting obstructions, with no geometric realization credit.

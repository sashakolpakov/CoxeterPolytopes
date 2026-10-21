# Independent audit of the analytic concentration insertion

20 October2026. **PASS for the stated concentration lemma and its low8
application.** I independently derived both formulas and the convex-order
argument by hand. This is a review of these two insertions only, not a
whole-paper audit or a rerun of the terminal cusp calculation. I did not
author this analytic insertion; no scientific process, imported checker,
TeX or Lean job was run.

The actively edited source is `manuscripts/finite_volume57.tex`. Bind each
excerpt from its preceding `\begin{lemma}` through its first following
`\end{proof}`, followed by one LF, encoded UTF-8:

* `lem:concentration`:1432 bytes, SHA256
  `aaf7b5c24931fca0f780676efb0d29aa1fdde3e8bd3d96dbfd5e9cb30c102585`.
* `lem:low8`:1317 bytes, SHA256
  `913d7519017e0a67cb4387a123fa6b6bfcf7a6a24dacd33330d118fd3f939c22`.

The locations at review time were lines724–756 and758–783; the excerpt
hashes, rather than those shifting line numbers or a whole-source hash,
are the review boundary. Only source extraction/hashing was executed.

Let X=U_a+U_b and Y=U_(a-1)+U_(b+1), with1<=a<=b. Both have mean
L/2 and symmetry about L/2, where L=a+b. Equal means make the differences
of right and left stop-loss expectations identical. Symmetry then gives
Delta_t=Delta_(L-t), and the differences vanish outside[0,L].

For integer0<=t<=a the two lower-tail convolution counts are k+1 for
0<=k<t. The sum of (t-k)(k+1) is binom(t+2,3). This gives exactly

    Delta_t=binom(t+2,3)*(1/[a(b+2)]-1/[(a+1)(b+1)]),

with denominator difference b-a+1>0. For a<=t<=L/2<=b, summing over
the first uniform gives independently

    E(t-X)_+=[t²+(1-a)t+a(a-1)/3]/[2(b+1)],
    E(t-Y)_+=[t²+(2-a)t+(a-1)(a-2)/3]/[2(b+2)].

Their difference has numerator

    t(L-t)-(a-1)(a+2b+2)/3

over2(b+1)(b+2), precisely the manuscript's second formula. That numerator
is increasing on[a,L/2] and equals(a+2)(b-a+1)/3 at a. The formulas agree
at the join. They also cover a=1 and a=b; in particular a=b=1 gives
Delta_1=1/12, so the U_0 boundary has not been omitted.

Stop-loss differences are linear between consecutive integer thresholds.
On the finite integer support, every convex function is an affine function
plus sum_t [f(t+1)-2f(t)+f(t-1)](x-t)_+, with nonnegative coefficients.
The affine terms cancel by equal means, proving the asserted convex order.
Condition on any common independent summand: every translate of a convex
function is convex, so integration preserves the inequality. The summands
used in the cusp application have finite support; no integrability issue
is hidden in that application.

For the low8 transfer,1<=a<=b<=7 makes all changed positive ranks at most8,
where c(r)=r/2, and c(0)=0 for the omitted zero factor. Thus normalized
charge is unchanged and no changed rank9 value is crossed. Applying convex
order to |2S-57| decreases G/V and hence normalized curvature and normalized
margin. The link product drops by replacing(a+1)(b+1) with a(b+2).

The integer potential sum r_i² increases by2(b-a+1)>=2, total rank stays57,
and57<=sum r_i²<=57². Therefore at most(57²-57)/2=1596 steps occur.
The stopping condition is exactly at most one positive rank below8,
with every other rank at least8, matching the displayed terminal domain.

Finally old V>=terminal V' and old margin m>=terminal m'>0 imply
Vm>=V'm'. Positivity is explicitly used at the correct step. Conditional
on the separately recorded terminal proposition V'm'>=36 and the geometric
actual-charge envelope, the claimed actual cusp surplus follows. This audit
does not independently verify the2552 terminal rows, their minima or the
geometric envelope; those are correctly separate inputs in this insertion.

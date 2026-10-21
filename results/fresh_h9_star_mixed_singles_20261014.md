# Shared single contacts with unrestricted noncore leaves

The specified shared-single family of star case 0 has no compact H9 Coxeter realization. The first three internal labels may be any finite integers at least 4, and the last may be any finite integer at least 7. Leaf 14 has its two own contacts labelled 3. Leaves 11 and 17 both have one own contact labelled 3 at the same last-block root and are right to the other last-block root. **All contacts involving noncore leaves are unrestricted by this theorem.** The other orientation of the last roots is included by a checked incidence automorphism.

This covers the [saved 1,269-visit assignment](fresh_h9_star_single_right_20261014_finite_core_search_batch_009.json), including the opposite-single profiles on leaves 8 and 16. It is a family exclusion, not an exclusion of the entire incidence or a metric construction.

The [checker](../scripts/fresh_h9_star_mixed_singles_20261014.py) and [certificate](fresh_h9_star_mixed_singles_20261014.json) use the positive eight-root support

\[
S=(3,4;\ 6,7;\ 9,10;\ 12,13)
\]

and the rank core \(S\cup\{11,14,17\}\). Normalized Gram diagonals are 1. Every pair in this core is a meeting pair in the literal incidence. The proof has two parts: a continuous negative-plane argument when the first internal label is at least 6, and an inherited affine determinant calculation with new core-only vertex witnesses when it is 4 or 5.

## Core domains

Twenty-four actual four-root supports force all cross-block contacts right: a trial nonright bond gives the affine path with labels 4–3–4 and positive null vector \((1,\sqrt2,\sqrt2,1)\). Increasing any bond preserves a nonpositive value on this vector. Thus the Gram on \(S\) is a positive direct sum of four rank-two blocks.

The actual triangle \(\{11,13,17\}\) forces the anchors 11 and 17 mutually right. The actual paths

\[
(3,14,11,13,17),\qquad (3,14,17,13,11)
\]

have affine Gram at labels \((3,4,3,3)\), with positive null vector \((2,4,3\sqrt2,2\sqrt2,\sqrt2)\). Consequently each contact from 14 to an anchor has label 2 or 3. These witnesses use only the rank core; no contact on an unused leaf is assumed.

Three actual eight-root supports retain one own root for each of 11, 17 and 14, together with all six foreign block roots. The retained own contact contributes \(1/4\) to the Schur projection cost. Every foreign nonright contact contributes at least \(1/2\), because a block with internal label at least 4 has inverse diagonal at least 2 and entrywise nonnegative inverse. Each leaf therefore has at most one foreign contact, necessarily labelled 3.

If a foreign internal label is at least 6, its inverse diagonal is at least 4. One nonright contact would then cost at least 1 in addition to the retained own cost \(1/4\). Such contacts are impossible. In particular all three core leaves are right to any middle block with internal label at least 6, and both anchors are right to the first block when that block has label at least 6.

Finally, the actual triangle \(\{12,13,14\}\) forces 14 right to both last roots. An additional nonright bond would give a three-node path with determinant below \(1-3/4-1/4=0\), since the last label is at least 7. This retains the full unbounded integer-label range.

## A negative plane when the first block is high

Suppose the first internal label is at least 6. Project the vectors \(r_{14}\) and \(r_{11}+r_{17}\) off \(S\). Their own-block projection costs are respectively

\[
X=\frac{D_A(1+c_A)}2\ge2,
\qquad D_{\mathrm{last}}>4,
\]

where \(c_A\) is the first block cosine and \(D=1/(1-c^2)\). Their projections onto the other own block vanish by the domain bounds above. Write \(H,F\ge0\) for their remaining squared projection lengths in the positive middle-block space, and put \(K=X-1\ge1\), \(L=D_{\mathrm{last}}-2>2\). The two residual norms are \(-(K+H)\) and \(-(L+F)\).

The ambient pairing of \(r_{14}\) with the anchor sum has absolute value at most 1. Cauchy's inequality bounds the remaining projection product by \(\sqrt{HF}\). Hence the residual determinant is strictly larger than

\[
(1+H)(2+F)-(1+\sqrt{HF})^2
=1+(\sqrt F-\sqrt H)^2+H\ge1.
\]

Both residual diagonals are negative. The residual space would therefore contain a negative-definite plane, contrary to Lorentz index one. This argument has no upper label bound and no field assumption.

## First internal label 4 or 5

Any middle internal label at least 6 has no core contact, as proved above. Replacing its internal label by 4 leaves the residual matrix exactly unchanged. Thus the relevant residual data reduce to the three finite internal choices \(4/5\), while retaining the variable last-block quantity

\[
D=\frac1{\sin^2(\pi/m)}>4.
\]

Both anchors project to the same last root, and 14 is right to that block. The three-by-three residual determinant is consequently affine in \(D\). The [frozen shared-high certificate](fresh_h9_star_shared_high_20261014.json), with its [independent arithmetic review](fresh_h9_star_shared_high_20261014_review.md), supplies all 7,840 coefficient cases. Every slope is positive. In 7,392 cases the value at \(D=4\) is already positive, giving rank at least 11 after adjoining the positive eight-root block. Such a Gram cannot lie in the ten-dimensional ambient Lorentz space.

The remaining 448 inherited cases formerly used vertex witnesses that could involve noncore leaves. The present checker replaces **all 448** with exact LDL failures on literal actual faces entirely inside the rank core. It checks nine maximal core faces and 60 distinct matrix types; the selected failures use 39 distinct full matrices. The last internal label is relaxed to 4; any replaced high middle label is also relaxed to 4. Raising those labels or adding other acute contacts preserves failure of positive definiteness by the nonnegative-eigenvector criterion for acute Gram matrices. Thus these witnesses require no noncore contact hypothesis.

The checker also verifies the complete 7,840-case key set, the inherited slope and boundary identities, all 27 fresh positive-kernel guards, the three actual budgets, and the literal swap of roots 12 and 13 on every vertex and minimal nonface.

## Verification and limits

The primary reuses the frozen exact field helper over \(\mathbb Q(\sqrt2,\sqrt5)\) for its new LDL calculations. It inherits the separately reviewed affine determinant bank rather than recomputing that bank. The high-first proof uses the displayed square identity and Cauchy's inequality. All dependencies are hash-bound.

```sh
python scripts/fresh_h9_star_mixed_singles_20261014.py --verify
```

Generation and saved replay passed. The [independent review](fresh_h9_star_mixed_singles_20261014_review.md) also passed generation and saved replay. Its separate Fraction arithmetic reconstructs all 448 replacement matrices, checks their principal determinants and quadratic witnesses, verifies the 27 core guards and three budgets, and checks the complete inherited case-key coverage. It also reviews the continuous high-first argument and literal last-root swap. The old determinant coefficients remain a stated, independently reviewed dependency.

The previous cursor and all earlier proof files remain unchanged. Search continuation is recorded separately below; the theorem alone does not exhaust the incidence.

## Saved continuation

The [gate binder](../scripts/fresh_h9_star_mixed_singles_20261014_gates.py) adds the two shared-root orientations and their literal star-reversal images to the previous 205 predicates, giving [209 proved predicates](fresh_h9_star_mixed_singles_20261014_gates.json). It rechecks the reversal on every actual vertex and minimal nonface and binds the new independent review. Generation and saved replay pass.

The [search wrapper](../scripts/fresh_h9_star_mixed_singles_20261014_search.py) changes only the namespace and bound gate references in the frozen resumable engine. [Batch 010](fresh_h9_star_mixed_singles_20261014_search_batch_010.json) resumes the 1,269-visit cursor and stops after two new visits, at **1,271 cumulative visits with nine pending frames**. The shared-single theorem removes the prior assignment. The next necessary assignment retains internal labels \((4,4,5,\ge7)\), the opposite-single profiles \((2,4),(4,2)\) on leaves 8 and 16, and changes the last two leaves to opposite single \(\pi/3\) contacts, \((2,3),(3,2)\). Other own contacts remain \((3,3)\).

That next state is outside the present theorem. The saved high-label category still means every finite integer label at least 7, and all dotted distances remain unassigned. It is neither a metric realization nor evidence of compactness. No pending branch is discarded or restarted. The final saved batch replay passed.

```sh
python scripts/fresh_h9_star_mixed_singles_20261014_gates.py --verify
python scripts/fresh_h9_star_mixed_singles_20261014_search.py --batch 10 --verify
```

# Opposite single contacts with unrestricted noncore leaves

This certificate excludes a family of compact H9 Coxeter realizations of star case 0: the four internal labels satisfy \((m_A,m_B,m_C,m_D)\ge(4,4,4,7)\); leaf 14 has own contacts \((3,3)\); leaves 11 and 17 have opposite single own contacts \((2,3)\), \((3,2)\). Every noncore contact is free. Both last-root orientations are included. The [saved 1,271-visit assignment](fresh_h9_star_mixed_singles_20261014_search_batch_010.json) belongs to this family, but the whole incidence is not excluded.

The [primary checker](../scripts/fresh_h9_star_opposite_core_20261014.py) uses only the core
\[
S\cup\{11,14,17\},\qquad S=(3,4;6,7;9,10;12,13).
\]
It inherits the [frozen opposite-high polynomial certificate](fresh_h9_star_opposite_high_20261014.json) and its [independent arithmetic review](fresh_h9_star_opposite_high_20261014_review.md). The new work removes the noncore hypotheses from every domain guard and exceptional vertex witness; it does not repeat the large coefficient bank.

## Literal core domains

Twenty-four actual four-root supports force the blocks of \(S\) mutually orthogonal by the affine 4–3–4 path and its positive null vector \((1,\sqrt2,\sqrt2,1)\). The actual eight-root supports that retain one own root of each core leaf give retained own cost \(1/4\). A foreign nonright contact costs at least \(1/2\), with nonnegative cross terms. Thus each core leaf has at most one foreign contact, necessarily labelled 3. A foreign block with internal label at least 6 would instead cost at least 1, so all such contacts vanish.

The actual triangle \(\{12,13,14\}\) forces 14 right to both last roots because \(m_D\ge7\). The actual path \((13,11,17)\) has determinant zero at labels \((3,6)\); hence the mutual anchor label is at most 5.

Both contacts from 14 to the anchors have label at most 4. To check this without noncore leaves, use the actual paths
\[
(3,14,11,13),\qquad(3,14,17,12).
\]
At labels \((3,5,3)\), the positive vector \((1/2,1,1,1/2)\) has quadratic value \((2-\sqrt5)/2<0\). Larger trial labels or additional acute pairings cannot restore positive definiteness. The two anchor contacts cannot both be 4: the actual path \((11,14,17)\) at \((4,4)\) has positive null vector \((1,\sqrt2,1)\).

The permutation \((11\ 17)(12\ 13)\) preserves all nine maximal induced core faces. It also preserves the own-contact pattern while exchanging the two third-contact labels. Therefore, for the rank calculation, one may arrange these labels in the domain \(\{2,3\}\times\{2,3,4\}\). This is an explicitly checked symmetry of the core face list, **not a claim of a full-incidence automorphism**. Every transferred vertex witness remains on an actual face. Separately, the swap \((12\ 13)\) is checked on every full vertex and minimal nonface to cover the other own-contact orientation.

## High first block: a negative plane

Suppose \(m_A\ge6\). The anchors are right to block A, and 14 is right to block D. Project \(r_{14}\) and \(r_{11}+r_{17}\) off \(S\). Their own-block costs are
\[
X=\frac1{2(1-c_A)}\ge2+\sqrt3,
\qquad Y=\frac1{2(1-c_D)}>2+\sqrt3.
\]
The anchor sum has ambient squared norm at most 2. Write \(H,F\ge0\) for the remaining projection norms in the positive middle-block space. The residual norms are \(-(K+H)\), \(-(L+F)\), where \(K=X-1>2\) and \(L=Y-(r_{11}+r_{17})^2>1\).

Each third-to-anchor cosine is at most \(1/\sqrt2\), so their sum is at most \(\sqrt2\). Cauchy's inequality bounds the middle projection pairing by \(\sqrt{HF}\). The residual determinant is consequently strictly greater than
\[
(2+H)(1+F)-(\sqrt2+\sqrt{HF})^2
=(\sqrt{2F}-\sqrt H)^2\ge0.
\]
Strictness holds even when the displayed square vanishes, because both diagonal lower bounds are strict. The negative diagonals and positive determinant give two negative directions, impossible in a Lorentz space of index one. No label cutoff or number-field assumption enters this argument.

## Low first block: inherited rank calculation, new witnesses

If \(m_A\) is 4 or 5, any middle internal label at least 6 has no core contact and drops out of the residual. Replace it by 4 without changing that residual. The complete domain then matches the frozen 47,040-case bank: three internal choices in \(\{4,5\}\), seven foreign patterns for each anchor, five for 14, four mutual-anchor labels, and third labels in \(\{2,3\}\times\{2,3,4\}\).

The independently reviewed polynomial identity and Bernstein bounds cover every finite last label at least 7. There are 43,412 positive cases on the coarse continuous interval and 16 on the sharper certified interval. Their nonzero residual determinants would give rank at least 11 after adjoining the positive eight-root block, whereas H9 has ten-dimensional ambient Lorentz space.

The other 3,612 cases now receive fresh exact LDL failures entirely inside the rank core. The search checks 262 matrix types; its selected failures use 240 distinct full matrices. The checker reconstructs each matrix, verifies an actual containing vertex, and verifies that the local core permutation also sends its support into an actual face. These failures hold with the last internal label relaxed to 4 and with unused middle blocks similarly relaxed. Increasing the internal labels or adding acute contacts preserves failure of positive definiteness. Thus no noncore own-contact value is needed.

## Verification scope

The certificate binds the old coefficient bank and independent review recursively, checks the full case-key domain, and performs the new domain and fallback arithmetic using the frozen exact field helper over \(\mathbb Q(\sqrt2,\sqrt5)\). It is not an independent recomputation of the inherited 47,040 determinants. The continuous high-first argument and the local core symmetry are new proof steps.

```sh
python scripts/fresh_h9_star_opposite_core_20261014.py --verify
```

Generation and saved replay passed. The [independent review](fresh_h9_star_opposite_core_20261014_review.md) also passed generation and saved replay. Its separate Fraction arithmetic checks all 27 quadratic guards, the mutual-angle and projection budgets, the local core involution, the complete 47,040 inherited keys/signs, and all 3,612 new fallback matrices and placements. It reconstructs independent determinant and nonnegative quadratic witnesses for the 240 selected full matrix types. The old polynomial bank remains an explicitly inherited, separately reviewed dependency.

The [literal-core transfer package](fresh_h9_star_core_transfers_20261014.md) records the placement generalization and subsequent continuation. The frozen 1,271-visit cursor is preserved as its input.

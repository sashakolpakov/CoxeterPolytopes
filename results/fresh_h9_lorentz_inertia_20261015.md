# Which local certificates are independent of dimension?

A substantial part of the saved star analysis excludes local configurations in **every Lorentz space**, not merely in hyperbolic dimension nine. The distinction is the inertia of the residual matrix, rather than its rank alone. The isolated negative-determinant exceptions remain rank obstructions specific to dimensions at most nine.

This note changes no frozen theorem or search output. Its [small checker](../scripts/fresh_h9_lorentz_inertia_20261015.py) and [ledger](fresh_h9_lorentz_inertia_20261015.json) count exact signs using separate rational arithmetic and bind the existing primary and independent certificates. They do not recompute the original matrix censuses or actual-face witnesses.

## The reusable local lemma

Let (G\succ0) be the Gram matrix of (k) roots in a space with at most one negative direction. For three further normals, write the full Gram matrix as
\[
\begin{pmatrix}G&B\\B^T&C\end{pmatrix},\qquad
R=C-B^TG^{-1}B.
\]
A block triangular congruence gives (G\oplus R), so (R) has at most one negative eigenvalue, independently of the ambient dimension.

If (det R>0) and one diagonal entry of (R) is negative, then (R) has exactly two negative eigenvalues. Indeed, its determinant is nonzero and positive, so its number of negative eigenvalues is even; the negative diagonal excludes zero. This contradicts Lorentz signature in **every dimension**.

The analogous two-vector criterion is simpler: two negative residual diagonals and positive determinant give a negative definite plane, also impossible in every Lorentz space.

By contrast, if (det R<0) and an actual two-by-two principal residual minor is indefinite, then (R) has signature ((2,1)). For (k=8), the full Gram matrix has signature ((10,1)). Its rank excludes hyperbolic dimensions at most nine but supplies no inertia contradiction in dimension ten or higher. This is not a realization of a compact polytope: other face and global closure conditions may still fail.

## Exact classification of selected saved packages

Every assertion below keeps its original literal spherical-face and label hypotheses. In a larger-dimensional compact Coxeter polytope, the same sets of normals must still span actual spherical faces. Merely embedding the graph without those face premises is insufficient.

| Saved package or bank | Dimension-independent part | Rank-only part |
|---|---|---|
| [All-own-three, sixteen internal-4/5 banks](fresh_h9_star_zero_20261014_own_three_review.md) | 20,824 positive residual determinants | 1,128 negative residual determinants |
| [Shared-high, first-three internal-4/5 coefficient bank](fresh_h9_star_shared_high_20261014_review.md) | 7,392 positive determinant regions; 448 actual-face failures | none in this bank |
| [Opposite-high, first-three internal-4/5 coefficient bank](fresh_h9_star_opposite_high_20261014_review.md) | 43,428 positive determinant regions; 3,612 actual-face failures | none in this bank |
| [Dense third-contact theorem](fresh_h9_star_third34_20261014_review.md) | all six negative-plane regions, 66 positive coefficients | none |
| [Shared anchors, single third contact at least 5](fresh_h9_star_single_third_20261014_review.md) | all eleven determinant regions, 78 positive coefficients | none |
| [Opposite anchors, single third contact at least 5](fresh_h9_star_single_third_opposite_20261015_review.md) | fifteen positive regions | one negative region |
| [Third contact exactly 4, shared or opposite anchors](fresh_h9_star_third_four_20261015_review.md) | all no-foreign cases and eight positive foreign regions | one negative foreign region |

For the high-anchor packages, a residual diagonal is (1-D/4-F<0), with (D>4) and (F\ge0). In the newer packages (D>100/19). In the all-own-three finite bank, the own projection cost is
\[
\frac{1}{2(1-c)}>1\qquad(c\ge\cos(\pi/4)),
\]
and foreign contacts can only increase it. Thus the required negative diagonal is present, rather than inferred from the determinant sign.

The checker independently recounts all 21,952 all-own-three determinant signs and all 784 negative anchor minors from their saved exact field elements. It converts the descending coefficients at (sqrt2+\sqrt5) to the basis (1,\sqrt2,\sqrt5,\sqrt{10}), then determines signs by rational comparisons. For the polynomial packages it independently checks the relevant coefficient signs. Original independent matrix and guard reviews remain explicit dependencies.

The table does not recertify every arbitrary-high-block reduction attached to the first three rows, nor every historical family in the repository. In particular, the whole all-own-three H9 predicate cannot be promoted to an all-dimensional forbidden configuration solely from this count.

## The two negative exceptions have exactly one negative direction

The opposite-anchor exception with owner label 4 and third label 5 has last label 7, a foreign contact at the third's right owner root, and a label-3 third contact to the other anchor. The third-label-4 theorem has the same exception with the owner and third labels exchanged. Their cleared determinant polynomials are identical and strictly negative in the certified label-7 interval.

Let (F) be the foreign projection cost. It is (1/2) in the first exception and ((5+\sqrt5)/10) in the second. In both cases (1/2\le F<1), (D>100/19), and the two-anchor minor equals
\[
1-F+D\left(\frac F4-\frac7{16}\right)
<\frac12-\frac{3D}{16}
<-\frac{37}{76}<0.
\]
Thus that minor is indefinite. Together with the negative three-by-three determinant, this proves residual signature ((2,1)), and full eleven-root signature ((10,1)). The exceptional pieces really are rank-based; a positive-determinant inertia argument cannot be applied to them.

The same distinction applies to the 1,128 negative entries in the finite all-own-three bank: its separately checked negative anchor minors make the residual signature ((2,1)). This observation does not certify all other spherical faces or any global polytope for those entries.

## What is missing for a global dimension bound

These are reusable local forbidden configurations. No argument here proves that an arbitrary high-dimensional compact Coxeter polytope must contain one. Most of the strongest continuous statements require a particular eleven-node incidence pattern and labels at least 7; neither increasing dimension nor a simplicial nerve alone supplies those hypotheses.

A global deduction would need an additional unavoidable-configuration result, for example a face-counting argument forcing one of the precise labeled configurations with all its spherical-face premises. Alternatively one could enlarge the forbidden catalogue until a separate exhaustive incidence theorem applies. Neither step is supplied by this ledger. Consequently this note proves **no new upper dimension bound**, and it does not turn the bounded H9 searches into an exclusion of all compact H9 polytopes.

Generation and saved replay of the ledger pass. All referenced theorem arithmetic and independent reviews are frozen and passing.

```sh
python scripts/fresh_h9_lorentz_inertia_20261015.py --verify
```

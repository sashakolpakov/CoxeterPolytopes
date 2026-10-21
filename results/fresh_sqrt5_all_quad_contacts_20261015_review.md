# Independent complete Q(sqrt5) contact-domain review

**PASS.** The [independent checker](../scripts/fresh_sqrt5_all_quad_contacts_20261015_review.py) and [certificate](fresh_sqrt5_all_quad_contacts_20261015_review.json) reproduce every necessary one-white domain in the [8,298-core primary census](fresh_sqrt5_all_quad_contacts_20261015.md). Generation and saved replay pass. No source or proof issue was found.

The independent search reverses both variable and label order and uses the separate finite-tree classifier from the frozen independent contact review. It applies no short-path or H4 pruning. All 81,772 nonzero contact vectors agree exactly. It also checks the terminal contiguous class list, all committed SQLite records, the final cursor, and the database digest. There is no unfinished core, label cutoff, or sampled-domain substitution.

For arithmetic checks, the reviewer uses standalone Fraction arithmetic and products along a separately chosen spanning tree to rescale roots into \(\mathbb Q(\sqrt5)\). It independently inverts an old spherical corner and computes the full three-by-three remaining Schur matrix, rather than the producer's two-by-two update. All 52,159 finite field-compatible profiles receive exact physical and conjugate inertia checks. Every reported affine rejection is independently verified by its positive kernel. The finite-cycle test uses the independent squareclasses of 2, 3, and \(2+\varphi\); its treatment of dotted cycles remains explicitly symbolic.

The reproduced core split is:

| Disposition | Count |
|---|---:|
| Certified affine principal subdiagram | 319 |
| Sealed nonsingular finite core | 1,224 |
| Sealed singular finite core | 636 |
| Finite core with nonzero accepted profiles | 4,044 |
| Retained symbolic dotted core | 2,075 |

The retained profile counts are 18,964 finite and 15,532 dotted. There are 175 dotted cores with no accepted nonzero profile. These counts describe necessary domains; neither an accepted profile nor a collection of individually accepted profiles is asserted to have a simultaneous geometric realization.

The domain completeness argument is valid. Every new white node belongs to all four enlarged spherical corners. Its contact to each old node is therefore finite and belongs to the complete arithmetic alphabet 2, 3, 4, 5, 6, 10. Assigning unchosen contacts to right angles is a safe positive-definiteness relaxation: on a symmetric matrix with nonpositive off-diagonal entries, decreasing an off-diagonal entry cannot repair a nonpositive least eigenvalue, as seen by taking an absolute-value minimizing vector. The independent exhaustive recursion confirms the resulting terminal sets. The all-right vector is intentionally omitted but always allowed.

Empty accepted nonzero domains give an all-size sealing statement. Any external white node directly touching the core would supply a tested necessary profile. An external path touching the core has such a first edge, so an arbitrary collection of external whites cannot evade the test. All external whites must therefore be orthogonal to the core. For a finite nonsingular core the full quadrilateral star would then be a nonsingular block sum of that core and a positive-definite white block. A compact simple two-face star has one more normal than the ambient Lorentz dimension, giving the claimed rank contradiction. A singular finite core instead retains its entire marked spherical components. Empty symbolic dotted domains also seal uniformly, but their unspecified distance and rank prevent the nonsingular conclusion.

The compactness and arithmetic hypotheses matter: affine positive kernels exclude compact ambient polytopes; conjugate positivity and the finite alphabet use the stated arithmetic ground field. Failure to find an affine witness is not promoted to an absence theorem. The complete original union list and its marked quotient are independently reviewed inputs, not rebuilt here. No global quadrilateral payment, arbitrary compatible multi-white realization, or improved dimension bound is claimed.

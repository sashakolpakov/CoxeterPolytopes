# Independent review of the Q(sqrt5) Lannér catalogue

**PASS.** The [independent checker](../scripts/fresh_sqrt5_lanner_catalogue_20261015_review.py) and [saved certificate](fresh_sqrt5_lanner_catalogue_20261015_review.json) reproduce the complete finite catalogue in the [primary report](fresh_sqrt5_lanner_catalogue_20261015.md): 12 triangles, six order-four classes and four order-five classes. The dotted order-two family remains symbolic, with squared norm-two contact \(t\in\mathbb Z[\phi]\), \(t>4\), and \(0<t'<4\). This review does not assert a new dimension bound.

The implementation imports no producer or other project code. It uses integer arithmetic in the basis \(1,\sqrt2,\sqrt5,\sqrt{10}\), doubled Gram matrices, permutation determinants, exact nested-quadratic signs and Sylvester minors. In particular, it does not reuse the producer's finite Dynkin classifier. It independently tests all 64 rank-three and 4,096 rank-four label matrices, obtaining 25 and 182 positive patterns. It enumerates higher-order graph supports by connected edge masks, independently of the producer's Prüfer-word construction.

All 56 unordered finite-label triples and their rejection reasons agree. All 675 order-four and 13,041 order-five tree/cycle assignments agree, including every geometric class, determinant, conjugate determinant, field rejection and positive-definiteness rejection. The labelled arithmetic survivor counts are 78 and 360. Determinants are converted independently to the producer's golden-ratio basis before comparison.

The all-label and topology reductions were separately source/proof reviewed. The two-cycle field condition gives exactly labels \(2,3,4,5,6,10\). The square classes of \(2,3,\phi+2\) are independent over \(\mathbb Q(\sqrt5)\), so the stated cycle parity test is complete. A connected Lannér diagram of order four or five cannot contain a bond at least six, by a proper three-node positive-definiteness test. Any cycle must be spanning and chordless; otherwise a proper cycle already has a nonpositive all-ones quadratic form. Thus the enumerated trees and spanning cycles exhaust the relevant supports. Conjugate positivity, cyclic field membership and integrality are checked with the stated normalization. The maximum Lannér order and arithmetic Gram criterion remain the published theorem dependencies identified in the primary report; this review does not reprove those external theorems.

Generation and saved replay pass. The certificate binds the final primary source and JSON, plus the independent source; the primary prose is not hash-bound. Full ambient attachments, opposite-piece unions, quadrilateral stars and any Q(sqrt5) dimension improvement remain outside this review.

```sh
python3 scripts/fresh_sqrt5_lanner_catalogue_20261015_review.py --verify
```

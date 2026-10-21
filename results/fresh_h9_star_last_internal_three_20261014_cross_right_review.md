The twelve cross contacts between the last A2 block and the first three blocks are forced right **after the dense last-own profile reduction**. All24 oriented affine-path witnesses occur in literal actual vertices. The [independent checker](../scripts/fresh_h9_star_last_internal_three_20261014_cross_right_review.py) and [certificate](fresh_h9_star_last_internal_three_20261014_cross_right_review.json) pass generation and saved replay.

The hypotheses are literal tree case0, last internal label3, first three internal labels at least4, and the prescribed first three leaf contacts3/3. Initially the two last leaf profiles, at facets11 and17, are free. Their missing triples are11,12,13 and12,13,17; the triples11,17,12 and11,17,13 are actual faces.

If a last leaf has a single own contact, its missing triple forces that label to be at least7. The actual triangle containing its high bond forces the other leaf right to the same root. The other missing triple consequently makes that leaf an opposite single contact, also at least7. This is precisely the separately reviewed opposite-single branch; its full exclusion is inherited, not repeated here.

Otherwise both profiles are dense. Each must contain a label at least4, since the all3 triangle would be affine. At either shared root, two own labels at least4 already violate positive definiteness, even when the mutual leaf contact is relaxed to right. The two high labels therefore occur at opposite roots. Each other own label must be3. A high label at least6, together with that3, also violates the actual triangle, so the dense profiles are exactly

\[
(p,3),\ (3,q),\qquad p,q\in\{4,5\},
\]

or their root reversals. A mutual leaf label at least3 makes the all-ones quadratic value at most \(1-\sqrt2<0\) in either actual anchor triangle. Thus the two leaves are mutually right.

Choose a last root \(d\in\{12,13\}\), its high-contact leaf \(h\), a root \(b\) in a first block, and that block's partner \(b'\). The certificate checks the actual support \(\{h,d,b,b'\}\) for each of the twelve choices and both orientations. If \(d\)–\(b\) were nonright, its Gram would dominate the relaxed path with labels

\[
h\mathrel{-^{4}}d\mathrel{-^{3}}b\mathrel{-^{4}}b'.
\]

That unit-diagonal path Gram has positive null vector \((1,\sqrt2,\sqrt2,1)\). Larger prescribed labels and any omitted nonright contacts only decrease its quadratic value on this vector. It therefore cannot be a positive definite actual-vertex submatrix. This forces every crossing right.

The conclusion is a conditional domain reduction, not a rank census or complete exclusion. The order of proof matters: the dense-own reduction supplies the high bond used by each path. The old internal4–3–4 argument cannot simply be reused with the last internal bond replaced by3.

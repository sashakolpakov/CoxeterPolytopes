# Independent review of the exact H5 candidate bank

The common-escape argument is sound for the **38,209 exact formula branches from the 12 frozen source banks**. They use 121 distinct nine-normal incidence spans. This excludes finite-volume completion by any compatible subset of that bank together with the retained 125 walls. Deleting optional walls only enlarges the intersection. It does not exclude other profiles, other branches, or arbitrary completions.

The reviewed implementation is [the primary verifier](../scripts/fresh_h9_global_h5_20261013_verify.py), with its [canonical source mapping](../scripts/fresh_h9_global_h5_20261013_schema.py), [recipes](../data/enumerations/fresh_h9_global_h5_20261013_recipes.json), and [exact certificate](../data/enumerations/fresh_h9_global_h5_20261013_verify.json). A separate small [Fraction checker](../scripts/fresh_h9_global_h5_20261013_review.py) and [audit artifact](fresh_h9_global_h5_20261013_review.json) bind their final hashes.

For each incidence span, the inverse Gram matrix defines the prescribed component (P) of a new normal. The component (z) of the common interior point perpendicular to that span is nonzero. The formula

\[
r=P+\sigma\sqrt{(2-P^2)/z^2}\,z
\]

therefore has norm two and the prescribed nine pairings. The code certifies invertibility, the nonzero complement norm, positive radicand, interior orientation, and strict negative pairing with the saved escape. This works for either sign of (z^2). The sign \(\sigma\) is defined and checked by projecting the literal rational-decimal proposal into the one-dimensional complement. The rounded proposal itself is **not** asserted to be an exact norm-two normal.

The unique `(bank, row)` check and canonical incidence/profile mapping cover every row of every source bank exactly once. The three alternate anchors are identical, including branch, to bank rows 7315, 4818, and 9678. They introduce no uncensused candidate normals. The fallback cosine routine gives rigorous bounds using Machin's identity and a conservative Taylor remainder; this is not a floating-point angle decision.

The H5 projection also has the correct sign. The four fixed walls are mutually perpendicular and have norm two. If (a_i=(x,e_i)/2<0), then (x_F=x-\sum_i a_i e_i). Every bank normal satisfies \((r,x)<0\) and \((r,e_i)\leq0\), hence \((r,x_F)\leq(r,x)<0\). The projected point has positive norm and the required future orientation. Thus the four fixed walls give a common strict escape within their H5 intersection; independently, the original H9 escape survives optional deletions.

The companion checker uses only the standard library and independent Fraction arithmetic. It recomputes the exact ambient escape norm and orientation, checks the inverse Gram matrix, recomputes the H5 projection bounds, verifies the recorded future field, and checks bank/anchor provenance. Its projected norm lies near (0.01242718848669>0). The maximum certified bank pairing is below (-0.00204116858075). The companion checker uses previously certified normal enclosures and the reviewed exact orthogonality identities; it is not a second implementation of all 38,209 root evaluations.

Validation: separate checker generation and saved replay **PASS**; a fresh full primary-verifier replay of all 38,209 rows also **PASS**. The source review found no remaining mathematical issue. No pairwise Coxeter compatibility of the entire bank is claimed or needed for this obstruction.

```sh
PYTHONPATH=src:. python scripts/fresh_h9_global_h5_20261013_verify.py --verify
python scripts/fresh_h9_global_h5_20261013_review.py --verify
```

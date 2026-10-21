# Independent review of all-length single-path completion

The [independent checker](../scripts/fresh_sqrt5_single_chain_completion_20261015_review.py) and [certificate](fresh_sqrt5_single_chain_completion_20261015_review.json) reconstruct all 458 profiles in the [primary certificate](fresh_sqrt5_single_chain_completion_20261015.json). Independent generation and saved replay pass, as do syntax, source bindings and report-link checks. The review uses the frozen separate Fraction arithmetic in the basis 1,sqrt5, with independent radical gauges and exact-square-root Gram entries; it imports no producer arithmetic.

The all-length rank argument is valid. If the nonsingular core Gram is H, the first path node has contact column v and scaled norm d, and c=v^T H^(-1)v, the length-l simply-laced path block is (d/2)Cartan(A_l). Its first inverse entry is 2l/[d(l+1)]. This formula follows for every positive integer l from the tridiagonal determinant recurrence: the block determinant is (d/2)^l(l+1), and its endpoint cofactor is (d/2)^(l-1)l. The determinant lemma therefore forces

```
(l+1)d - 2lc = 0,
l = d/(2c-d).
```

The denominator-zero case has no finite solution because d is positive. Otherwise there is at most one possible length, which must be a positive rational integer in the coefficient field. This is not a finite length cutoff.

The review independently reconstructs each core matrix and its inverse, checks the inverse identities, reconstructs all first-node contact columns and recomputes all 458 values of d, c and the resulting field element l. The exact dispositions match:

| Disposition | Count |
|---|---:|
| Zero denominator; no finite path length | 37 |
| Candidate is not a positive integer | 258 |
| Unique positive integer; principal-star checks pass | 163 |

The 163 lengths are 147 of length 1, nine of length 2, five of length 3 and two of length 4. These are the values produced by the exact equation, rather than a prescribed search range.

For every surviving length, the checker rebuilds the full labeled matrix. It independently verifies physical inertia `(N-2,1,1)`, conjugate inertia `(N-1,0,1)`, and all 652 finite corner matrices in both embeddings. It directly inverts each path block and verifies its endpoint entry and zero determinant factor. The literal full-star labels, absence of a distance-two payment, and absence of the specified entire-H4 endpoint payment agree with the primary record. The primary's affine-witness search remains an inherited scoped check; this review does not turn failure to find an affine witness into a new exhaustive no-affine assertion.

The theorem covers one simply-laced external path whose first node alone touches the fixed core. Its first node may have several core contacts. Additional white components orthogonal to both the core and the entire path contribute positive-definite blocks and do not alter the rank condition. Branching attachments, multiple touching components, paths touching the core at both ends, dotted cores and arbitrary further completions remain outside scope. The 163 survivors are principal stars, not new polytopes or a global dimension bound.

```sh
python3 scripts/fresh_sqrt5_single_chain_completion_20261015_review.py --verify
```

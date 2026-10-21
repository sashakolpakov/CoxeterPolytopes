# Independent review of the two-high-block T8 obstruction

**PASS.** The literal T8 incidence has no compact H9 Coxeter realization satisfying only

```text
m(1,2)>=7,  m(4,5)>=7,  m(7,8)>=4,  m(12,13)>=4,
```

where these are finite integer labels. All other meeting labels and all divergent distances are free. The proof does not inherit the earlier fixed contacts of wall 14, the earlier own-contact patterns, or any polynomial parameter bank. This excludes a family inside T8, not the entire T8 incidence.

The [independent checker](../scripts/fresh_h9_t8_sparse_search_20261014_two_high_blocks_review.py) and [certificate](fresh_h9_t8_sparse_search_20261014_two_high_blocks_review.json) bind the [frozen primary certificate](fresh_h9_t8_sparse_search_20261014_two_high_blocks.json). They use SymPy exact symbolic arithmetic and no producer/project imports. Generation and saved replay pass. The checker independently builds all supports, checks their literal containment in saved actual vertices, reconstructs all 28 positive-null matrices, and verifies the exact scalar bounds. The original primal realization of T8 is an inherited certificate, not re-enumerated here.

## Actual supports and label dependencies

The eight roots `S={1,2,4,5,7,8,12,13}` lie in one actual vertex and therefore have positive-definite Gram matrix. Every one of the 45 pairs in `S+{11,14}` is an actual meeting pair, so all its labels are finite. The three literal missing faces `{1,2,14}`, `{4,5,14}`, `{11,12,13}` are checked directly; their three-node Gram matrices must be Lannér. In particular, each of the two contacts of wall 14 with its two high blocks is nonempty, as is wall 11's contact set with its last block.

All 24 cross-block pairs in S are forced right. Twenty are ruled out by a relaxed `(6,3)` path and four by a relaxed `(4,3,4)` path, each inside a checked actual vertex. The displayed kernel vector is strictly positive. Increasing an assumed internal label or restoring any omitted acute contact can only decrease its quadratic value. Thus these guards depend on the four displayed internal bounds alone; no outside leaf label is retained in their test matrices.

For wall 14, all four choices of one active root from each high block have an actual seven-root support consisting of these two roots, both complete low blocks, and wall 14. In the normalized Schur complement, the two active roots cost at least `1/4+1/4`; any nonright contact to either complete low block costs at least `1/2`. Positive definiteness therefore forces wall 14 right to roots 7,8,12,13.

For wall 11, both choices of an active own root have an actual eight-root support retaining that root, both high blocks, the other low block, and wall 11. A contact to either high block already costs more than one, before its positive own-root cost. Thus wall 11 is right to roots 1,2,4,5.

Finally, all four stars `{11,14,a,b}`, with `a in {1,2}` and `b in {4,5}`, occur in actual vertices. Choosing active a,b, an assumed label at least four between 11 and 14 gives the relaxed central star with cosines `sqrt(2)/2,1/2,1/2`, whose positive null vector is `(sqrt(2)/2,1,1/2,1/2)`. Consequently their mutual label is two or three. Every contact used in this argument is either one of the four hypotheses, chosen nonzero from a literal missing triple, or forced right by a preceding guard.

## Uniform own-block cost and the negative plane

For an internal label `n>=4`, write `c=cos(pi/n)`. A Lannér triple adjoining that block has projection cost strictly greater than `5/4`:

* With both own contacts nonzero, the cost is at least `1/(2(1-c)) >= 1+sqrt(2)/2 > 5/4`.
* With one nonzero contact, at n=4 its label must be at least five; the cost is at least `(3+sqrt(5))/4`.
* At n=5 it must be at least four; the cost is at least `1+sqrt(5)/5`.
* At n=6 it must be at least four, giving cost at least two. The label-three case is affine and therefore cannot be a compact missing face.
* At every n>=7, even a label-three contact costs `1/(4(1-c^2))>5/4`.

The last inequality has no finite cutoff. The identity `U6(cos(pi/7))=0` gives `f(c^2)=0` for `f(t)=64t^3-80t^2+24t-1`. The actual root lies above `3/4`; f is strictly increasing there and `f(4/5)=-29/125`, so `cos^2(pi/7)>4/5`. Monotonicity in n finishes the all-label argument. The independent checker verifies the recurrence identity and exact derivative/value bounds.

After eliminating the positive S block, wall 14 has squared residual below `1-2(5/4)=-3/2`; wall 11 has squared residual below `1-5/4=-1/4`. Their projections use disjoint blocks, so their residual mutual pairing is their ambient pairing, zero or `-1/2`. The residual determinant is therefore strictly greater than `3/8-1/4=1/8`. This negative-definite two-plane contradicts the ambient Lorentz form's negative index one. The primary uses the weaker sufficient wall-14 bound below minus one; both versions are valid.

No mathematical or dependency issue was found. The all-other-labels-free scope is justified by the literal guards themselves, rather than merely by the size of their union.

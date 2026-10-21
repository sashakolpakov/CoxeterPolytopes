# Independent complete replay of the two-norm height-100 prefix

**PASS.** A separate enumeration reproduces every candidate level and accepted root in the [height-100 transcript](fresh_h9_beta11_vinberg_20261014.md). It checks all 46 levels in both norm classes, obtaining the same twelve positive-height walls and 21 total walls. This is a complete independent bounded enumeration, not a replay using the producer's arithmetic or pruning tree.

The [checker](../scripts/fresh_h9_beta11_vinberg_20261014_review.py) and [certificate](fresh_h9_beta11_vinberg_20261014_review.json) use standalone `Fraction` field arithmetic and a separate integer-pair hot loop. They import neither `Golden`, the producer, any project Schur helper, nor symbolic algebra. Generation and saved `--verify` pass.

## Independent finite bounds

Write `q∈{2,2β}`, `δ=φβ`, and `h=δx²/(2q)≤100`. Since `q≤2β`,

`0<x²≤400/φ<400`.

At the conjugate embedding, positivity gives `x′²≤−2q′/δ′≤4φ<9`. Thus `0<x<20` and `|x′|<3`. For `x=a+bφ`, the identities `b=(x−x′)/√5` and `a=x′−bφ′` justify the independently chosen coefficient box `|b|≤12`, `|a|≤15`. Exact filters in this box reproduce ten short-norm and 36 long-norm control values. They are independently sorted by exact geometric height.

For the positive coordinates `(λ,l)`, the target is `T=2q+δx²≤2q(1+100)<2020`, while `T′≤2q′≤4β′<16`. Each coordinate has square at most the target at both embeddings. Hence `0≤z<45`, `|z′|<4`, and the independent universal coefficient box `|b|≤25`, `|a|≤29` contains every coordinate. Exact filtering reproduces every saved per-level alphabet; no producer coefficient bounds are reused.

## Different recursion and selection

The checker reconstructs the literal E8 Cartan matrix `C` from the host Gram and inverts it using separate rational Gaussian elimination. It reverses the producer's variable order, using

`(λ8,λ7,…,λ1,l)`.

The positive form is the corresponding permutation of `diag(2C⁻¹,1)`. Its inverse has diagonal one, which independently verifies the coordinate bounds above. Rather than eliminate a suffix as the producer does, the checker obtains the minimum quadratic energy at each fixed prefix by **inverting the corresponding principal block of the inverse positive form**. These rational matrices are stored in the review certificate.

Both real embeddings are used in every Schur-energy rejection. The norm-`2β` residue conditions are checked independently using `φ≡8 (mod β)`. The final `l` coordinate must satisfy the exact two-coefficient parity condition `l≡x (mod 2O)`.

The earlier-wall inequalities are accumulated only from the reviewer's own previously enumerated roots. Their omitted terms are nonnegative because all roots lie in the chosen zero sector. The producer's accepted roots are consulted only after a level is complete, to compare full sets. At every accepted leaf the physical norm and reflection-integrality condition are checked directly.

The independent recursion visits **60,587 nodes**, compared with the producer's 68,242. Its complete level sets agree, including every empty level and each of the seven nonempty levels. The different coordinate order and minimum-energy construction give a different pruning tree while preserving the exhaustive domain.

## Geometry and dependencies

After the enumeration, separate physical arithmetic checks all 441 entries of the full Gram matrix, all Cartan integrality conditions, the strict future center, and all 21 facet witnesses against all walls. No geometry failure was found.

The complete two-norm alphabet and complete height-zero arrangement remain the separately certified [norm-pool](fresh_h9_h8_cap_height_20261014_norm_pool_review.md) and [zero-sector](fresh_h9_beta11_zero_sector_20261014.md) dependencies. Their literal matrices and source hashes are checked here; their universal proofs are not presented as a new enumeration by this reviewer.

The result certifies the full positive-height reflection search **through height 100 only**. It does not claim termination, a complete global corner census, reflectivity, or finite volume. It does not constrain nonintegral geometric completions. All primary artifacts were read without modification.

```sh
python scripts/fresh_h9_beta11_vinberg_20261014_review.py --verify
```

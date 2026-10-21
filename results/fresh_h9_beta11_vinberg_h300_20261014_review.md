# Independent complete replay through height 300

**PASS.** The separate verifier exhausts all **80 candidate levels** in both primitive norm classes and reproduces every accepted root in the [height-300 transcript](fresh_h9_beta11_vinberg_h300_20261014.json). The result is 26 positive-height walls and **35 total actual walls**. All 34 levels beyond the frozen height-100 prefix are included; the verifier also re-enumerates the earlier 46 levels from scratch.

The [independent checker](../scripts/fresh_h9_beta11_vinberg_h300_20261014_review.py) and [review certificate](fresh_h9_beta11_vinberg_h300_20261014_review.json) use standalone `Fraction` field arithmetic and integer-pair search arithmetic. There are no imports of the producer, `Golden`, project Schur helpers or symbolic algebra. Generation and saved replay pass. The earlier height-100 primary and review are unchanged.

## Complete coefficient bounds at the larger height

For `q∈{2,2β}` and `h=δx²/(2q)≤300`, where `δ=φβ`,

`0<x²≤1200/φ<28²`.

Conjugate positivity still gives `x′²≤4φ<9`. Thus `0<x<28` and `|x′|<3`. The identities `b=(x−x′)/√5` and `a=x′−bφ′` place every `x=a+bφ` in the independent coefficient box

`|b|≤16`, `|a|≤19`.

Exact filtering produces 18 short-norm and 62 long-norm control values. They are sorted and compared using exact geometric height; all 80 levels are distinct.

The positive-form target satisfies `T=2q+δx²≤2q(1+300)<6020<78²`, and `T′≤2q′≤4β′<16`. Consequently every nonnegative positive-form coordinate satisfies `z<78` and `|z′|<4`. The independent universal coefficient box is therefore

`|b|≤42`, `|a|≤46`.

All inequalities used for these boxes are checked exactly. Each full box is filtered at both embeddings, reproducing the producer's coordinate alphabet at every level without importing its box construction.

## Independent exhaustive search

As in the [independent height-100 review](fresh_h9_beta11_vinberg_20261014_review.md), the coordinate order is reversed to `(λ8,…,λ1,l)`. The literal E8 Gram is inverted by separate rational Gaussian elimination. Minimum energies at a fixed prefix are obtained by inverting principal blocks of the inverse positive form, rather than by the producer's suffix-elimination routine.

The search imposes the exact norm-`2β` residue condition on all coordinates and the coefficientwise parity condition `l≡x (mod 2O)`. Both real embeddings are used for each quadratic rejection. Earlier-wall pruning uses only independently accumulated roots, with nonnegative remaining coefficients. Complete root sets are compared with the primary only after finishing each level.

The independent tree visits **792,480 nodes**, compared with the producer's 1,082,417. Every level matches, including empty levels and the fourteen additional walls above height 100. Direct arithmetic additionally verifies all **1,225 Gram entries**, Cartan integrality, the strict future center, and all 35 facet witnesses against every wall.

## Why no separate primitive-vector filter is missing

The reviewer independently checks that the literal lattice is even, that `β` has prime norm eleven, and that `β` is totally positive. If a vector in either enumerated norm class were `r=αs` with a nonunit common divisor `α`, evenness would imply

`b=r²/2=α²(s²/2)`, with `s²/2∈O`.

But `b` is either `1` or the prime element `β`, neither of which has a nonunit square divisor. Thus every integral vector of norm `2` or `2β` is already primitive. This argument uses the PID property established in the reviewed norm-alphabet proof. It is stronger than merely checking primitivity of the accepted walls and explains why a gcd test is unnecessary in the search.

Reflection integrality remains a separate condition for norm `2β`; it is enforced during enumeration and checked again on every accepted vector.

## Scope and reproduction

The complete norm alphabet and zero-height sector remain separately certified dependencies. This is a complete independent **positive-height enumeration through 300**, not a claim about greater heights, termination, global vertex closure, finite volume or reflectivity. No primary artifact was modified.

```sh
python scripts/fresh_h9_beta11_vinberg_h300_20261014_review.py --verify
```

The review certificate binds the frozen height-300 primary source and data, host basis, zero sector, norm-alphabet review, and its own checker.

# Independent review of the two-chain completion census

**PASS for the stated two-path construction.** The independent
[checker](../scripts/fresh_sqrt5_two_chain_completion_20261015_review.py) and
[certificate](fresh_sqrt5_two_chain_completion_20261015_review.json) reconstruct
all 1,455 profile pairs, all 239 exact rank equations, all 107 isolated
integer solutions, the five corner-bounded infinite rank families, and all
112 surviving principal stars. Generation and saved replay pass.

The review binds the frozen [primary source](../scripts/fresh_sqrt5_two_chain_completion_20261015.py)
and [data](fresh_sqrt5_two_chain_completion_20261015.json). It imports only
earlier independent review helpers: tuple-Fraction arithmetic in the basis
1,sqrt(5), separately reconstructed radical Gram entries, full congruence
inertia and a finite-Dynkin classifier. It imports no producer arithmetic.
Completeness of the inherited 458 one-node profiles and their original
compact/arithmetic hypotheses remains a dependency, not a new census here.

| Independent check | Count |
|---|---:|
| Unordered profile pairs, with replacement | 1,455 |
| Nonfinite two-head corners | 880 |
| Two physical negative directions | 283 |
| Bad conjugate inertia | 53 |
| Exact two-length equations | 239 |
| Isolated positive-integer solutions | 107 |
| Infinite rank families excluded above free length one | 5 |
| Surviving principal stars | 112 |
| Reconstructed finite corners of those stars | 448 |
| Direct inverses of their literal chain blocks | 224 |

All five surviving rank families have a fixed chain length of one. Four
have an affine corner at free length two; the fifth has an indefinite
corner. The independent checker rebuilds the literal matrices and verifies
the saved inertia and affine kernels. Every longer path contains this
non-positive-definite corner as a principal subdiagram, so the exclusion
is valid for all lengths at least two. It is not a finite length scan.

## All-length and normalization audit

For core Gram matrix H, contact vectors v,w and first-node norms d,e, the
two paths have Gram matrices `(d/2)Cartan(A_l)` and
`(e/2)Cartan(A_m)`. After eliminating the core their end couplings are
determined by `v^T H^-1 v`, `w^T H^-1 w` and `v^T H^-1 w`. Eliminating the
path interiors leaves diagonal corrections `d(l+1)/(2l)` and
`e(m+1)/(2m)`. Thus the cleared rank equation is exactly

`(1-a-b+f)lm + (1-a)l + (1-b)m + 1 = 0`,

with the coefficients defined in the primary report. The independent
inverse identities and all 239 coefficient reconstructions agree. Direct
inverses of the completed path blocks separately verify the first inverse
entry `2l/[d(l+1)]`, including the scale of the radical gauge. Both real
coordinates of the equation must vanish: checking only its physical
embedding would not suffice.

The new length solver derives the two rational coordinate equations. A
nonzero resultant in l gives an exact quadratic or linear integer-root
problem. When that resultant vanishes identically but some coefficient is
nonrational, the form factors as `(B l+1)(C_0 m+1)` and its nonrational
factor cannot vanish at a positive integer l. For rational coefficients,
the identity `(A l+C)(A m+B)=BC-AD` gives the complete signed-divisor list;
zero-product and linear degeneracies are treated explicitly. Every linear
case actually encountered here has a finite positivity bound or a fixed
variable. No numerical root test or imposed maximum length is used.

The source-level audit also checked the producer's more general linear
Diophantine parameterization and positivity inequalities. That unused
generic branch is not claimed as an independently exercised infinite
family in this dataset.

## Precise limitation

The result is complete only for two disjoint simply-laced white paths with
their first nodes alone touching one of the 130 specified nonsingular
cores. All full-label diagrams and paired source-profile indices are
checked literally. An orthogonal positive-definite spectator block cannot
change the singularity condition. Branching, several touching nodes in a
component, contacts between the paths, additional touching components and
dotted cores are outside the theorem.

For every one of the 112 records the review verifies all four finite
corners, physical inertia `(N-2,1,1)`, conjugate inertia `(N-1,0,1)`, and
absence of the short-path/H4 payment. These necessary principal-star tests
do not certify an actual future-oriented compact quadrilateral or a global
polytope. The producer's bounded affine-witness search is not promoted to
an exhaustive absence-of-affine theorem. No dimension bound or arbitrary
attachment closure follows from this certificate alone.

```sh
python3 scripts/fresh_sqrt5_two_chain_completion_20261015_review.py --verify
```

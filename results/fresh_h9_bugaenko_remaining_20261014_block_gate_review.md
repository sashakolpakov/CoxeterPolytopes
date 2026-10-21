# Review of the two-Lanner-block obstruction

**PASS.** The proof and its literal catalogue application are sound. Fresh
saved replays of the [primary checker](../scripts/fresh_h9_bugaenko_remaining_20261014_block_gate.py)
and the [separate Fraction checker](../scripts/fresh_h9_bugaenko_remaining_20261014_block_independent.py)
both passed. The scope is the 2,365 catalogue types selected by this block
criterion; this review does not certify the separate searches on the other
588 types.

The additional [stdlib witness checker](../scripts/fresh_h9_bugaenko_remaining_20261014_block_gate_review.py)
and [hash-bound result](fresh_h9_bugaenko_remaining_20261014_block_gate_review.json)
reconcile all 42 Schur determinant polynomials with the separately computed
Fraction results, recompute all 64 Bernstein coefficients, and prove each
strictly negative using rational enclosures obtained from integer square
roots. Eleven completions contain a parameter. This is a further witness
cross-check, not a third full category census or determinant engine.

## Why the signature argument works

The first block has order 4 or 5. Its Gram matrix `A` is invertible with
exactly one negative eigenvalue: all proper principal submatrices are
positive definite and its determinant is negative. These facts are checked
in the primary arithmetic; the separate replay also obtains the elimination
pivot signs `(+,...,+,-)`.

For the full two-block matrix the exact congruence to

\[
 A\;\oplus\;S,\qquad S=B-C^T A^{-1}C
\]

adds the inertia of the two symmetric matrices. A negative determinant of
`S` implies at least one negative eigenvalue of `S`; it is unnecessary to
infer its complete inertia from the determinant. Thus the full matrix has
at least two negative eigenvalues and cannot be a principal Gram block in
a Lorentz space with negative index one. When every cross pair is right,
the two original Lanner blocks already provide the two negatives.

## Why the category enumeration covers all labels

The hypothesis requires every union of a proper subset of each block to
be spherical. In particular every cross pair occurs in an actual
spherical face, so it cannot be dotted or parallel.

A failed trial that sets one cross label to 3 and every other cross label
to 2 excludes every nonright value for that pair. Indeed a normalized
Coxeter Gram is `I-W`, with `W` entrywise nonnegative. Increasing labels or
restoring omitted contacts increases `W`; Perron–Frobenius monotonicity
prevents a non-positive-definite trial from becoming positive definite.
This justifies testing only the listed potentially nonright pairs jointly.

Labels 2 through 6 are explicit. Category 7 stands for **every** finite
integer label at least 7. This is exact for spherical tests: such an edge
can occur only as an isolated rank-two `I2(m)` component. The classifier
accepts arbitrary rank-two labels and rejects them in larger components.
It is also exact in the rank-three Lanner table. The hyperbolic triangle
condition `1/a+1/b+1/c<1` has the same truth value when any entry at least
7 is replaced by 7. The complete order-four and order-five Lanner tables
have only labels at most 5.

All cross-category combinations are enumerated, including several
independent category-7 slots before rejection. Every surviving nonright
completion has cross labels at most 6 and at most one category-7 label in
the two internal blocks. Hence using a single variable `x=cos(pi/m)` in
its Schur determinant does not identify different unbounded parameters.
Cases with several unbounded internal labels and forced-right cross
contacts are already excluded by the orthogonal-block argument.

For the eleven parameterized completions, `cos(pi/m)` lies in `[9/10,1]`
for all `m>=7`. The stdlib check independently verifies

\[
 T_7(x)+1=(x+1)(8x^3-4x^2-4x+1)^2,
\]

the cubic value `-1/125` at `9/10`, and its positive increasing derivative
on `[9/10,1]`. Its largest root is `cos(pi/7)`. Strictly negative
Bernstein coefficients on that closed interval imply a strictly negative
Schur determinant throughout it, including the limiting endpoint 1.

## Counts and incidence application

| Block orders | Internal pairs | Forced-right pairs | Nonright completions |
|---|---:|---:|---:|
| 4,3 | 396 | 391 | 10 |
| 5,3 | 220 | 199 | 27 |
| 4,4 | 81 | 81 | 0 |
| 5,4 | 45 | 44 | 2 |
| 5,5 | 25 | 22 | 3 |
| Total | 767 | 737 | 42 |

For each selected catalogue case, the code checks two disjoint minimal
nonfaces of the required orders. More importantly, it directly verifies
that deleting one vertex from each block leaves a set contained in an
archived actual combinatorial vertex. Every smaller proper union is then
spherical in any compact Coxeter realization. The separate replay checks
all 2,365 witnesses against the literal incidence data. This does not rely
on describing the entire polytope as a product, and no unused dotted
distance can affect the forbidden principal block.

The complete classical Lanner/finite-diagram tables and the source
incidence catalogue remain shared mathematical inputs. Their bindings
are explicit. The arithmetic replay imports neither the primary
calculation nor SymPy arithmetic: it uses squarefree-radical Fractions,
Galois-conjugate inversion, permutation determinants, and a separately
written finite-tree classifier. Source review found no completeness or
signature gap within the stated hypothesis.

# Independent height-100 beta131 checks

The [independent prefix verifier](../scripts/fresh_h9_beta131_h100_20261016_review.py)
and [certificate](fresh_h9_beta131_h100_20261016_review.json) reproduce all
**114 controls and 41 walls** through height 100. Natural covector order,
separate enlarged coefficient boxes, the weaker coordinate alphabet and
project Golden arithmetic give 64,680 recursion visits, compared with the
producer's 66,878. Every accepted root and every control agrees exactly.
The first 52 controls also reproduce the frozen independent height-20
calculation literally. Generation and fresh saved replay pass.

The [separate witness checker](../scripts/fresh_h9_beta131_h100_witnesses_20261016_review.py)
and [certificate](fresh_h9_beta131_h100_witnesses_20261016_review.json)
independently establish the following:

- All 45 saved ambient finite vertices are actual, have independent
  rank-nine active supports, and satisfy every current wall inequality.
- A strict future spacelike escape still proves the current prefix has
  infinite volume.
- Wall 8 is exactly t=0. Its eleven perpendicular walls are literally
  the archived compact Bugaenko H8 roots, in the source order
  `[11,0,2,3,4,5,6,7,9,1,40]`.
- A fresh exhaustive enumeration of all **165 eight-wall supports** in
  that nine-dimensional source cone gives exactly **41 finite vertices**
  and no other rays. This uses separate Gaussian elimination rather
  than the producer's incremental double description.
- Every one of the other 29 walls projects to a strictly future timelike
  vector in t=0. Its inequality is therefore redundant on the entire
  compact H8 facet, not merely at selected vertices.
- All 66 saved integral corner maps and their integral inverses preserve
  the form and fix the same future timelike vector, of norm
  `(13-9*phi)/29`. Their generated group is finite: it lies in a compact
  point stabilizer inside the discrete integral orthogonal group.

The witness checker passes generation and fresh saved replay. A separate
source review confirms its complete H8 support argument and redundancy
test. It does **not** independently reproduce the complete ambient
41-wall cone census or the domain of the corner-map search. The producer's
45 finite / zero ideal / 1,891 hyperideal census remains the primary
complete calculation; the independent claims above have their stated
narrower scopes.

These results neither construct a finite-volume H9 chamber nor prove
nonreflectivity of the full lattice. The new lattice is distinct from the
recently excluded determinant-norm-11/55 golden hosts. Its exact compact
H8 facet and four current ambient finite vertices outside that facet give
concrete geometry for further continuation.

```sh
python3 scripts/fresh_h9_beta131_h100_20261016_review.py --verify
python3 scripts/fresh_h9_beta131_h100_witnesses_20261016_review.py --verify
```

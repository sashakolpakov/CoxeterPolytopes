# Independent review of the two-cycle quadrilateral rank obstruction

The [standalone review checker](../scripts/fresh_arithmetic_sqrt2_quad_rank_20261015_review.py) and [review data](fresh_arithmetic_sqrt2_quad_rank_20261015_review.json) independently reproduce the exact algebra and finite placement census of the [primary source](../scripts/fresh_arithmetic_sqrt2_quad_rank_20261015.py) and [certificate](fresh_arithmetic_sqrt2_quad_rank_20261015.json). The reviewer uses its own Fraction arithmetic in \(\mathbb Q(\sqrt2)\), Gaussian determinants and full-matrix LDL tests; it imports no producer or project arithmetic code.

**Scope.** The complete connected part of a proposed quadrilateral star consists of two disjoint chordless cycles of lengths 3, 4 or 5, each having exactly one bond 4 and all other bonds 3, with at least one cycle of length 4 or 5. Each contains one opposite black pair. They are joined by one white bridge, with arbitrary pendant white spherical trees attached at single white vertices; other white spherical components may be orthogonal. Shared white vertices, extra cross edges, extra black contacts, other Lannér triangles, and two bare three-cycles are not classified. This is not a dimension theorem or an exhaustion of arithmetic quadrilateral stars.

## Independent exact checks

The reviewer enumerates all **4,032** choices of the two black pairs, two white bridge endpoints, and bridge label 3 or 4. Each of the four proposed finite corners is tested by a full quadratic-field LDL computation, independently of the producer's rational leaf-elimination method. Exactly 28 placements pass all four tests:

| Ordered cycle lengths | All-four-PD placements |
|---|---:|
| 3,4 and 4,3 | 6 each |
| 3,5 and 5,3 | 8 each |
| 4,4; 4,5; 5,4; 5,5 | 0 |

The last four rows account for 3,528 tested placements. Every one of the 28 survivors has nonzero full determinant, checked both by direct Gaussian elimination and the bridge-cofactor formula. Therefore none has the rank required of the complete quadrilateral star.

The separate cycle calculation counts determinant terms by matchings, rather than by the producer's symbolic determinant engine. With unit root squares and variable rational diagonals, it gives

\[
\det L=a-b\sqrt2,\qquad
b=2^{1-r},\qquad
a(1,\ldots,1)=
\begin{cases}0&r=3,\\-1/16&r=4,5.\end{cases}
\]

The matching terms are rational. The two oriented full cycles are the only irrational terms, so their coefficient does not depend on any diagonal. Every proper principal subdiagram is independently checked positive definite, and every deleted-node cofactor is rational.

## Why the proof covers arbitrary pendant white trees

This extension is analytic, not inferred from a finite sample of branches. A pendant tree is positive definite, since it lies in the white support. Its Schur elimination changes only the diagonal of its unique attachment vertex, subtracting a nonnegative rational quantity. It creates no cross edge. Rationality follows from the matching-polynomial determinant formula for trees with bond squares \(1/4\) or \(1/2\).

The restriction to these bond squares must be justified, not assumed from the quadratic ground field alone: label 8 also has its squared pair in \(\mathbb Q(\sqrt2)\). A bridge or attached pendant bond belongs to an actual spherical corner component of rank at least three. Finite bonds 6 and 8 occur only in isolated rank-two components, so they cannot occur in this connected component. The only nonright labels left are 3 and 4. The same observation proves that the bridge square is rational.

Conjugate positivity of the arithmetic star survives Schur elimination of its positive-definite white trees. The resulting cycle principal block is positive semidefinite at the conjugate embedding. Its determinant is nonzero, because its irrational coefficient is the fixed nonzero \(b\), while its rational part is rational. Consequently this conjugate block is positive definite.

Increase its reduced diagonals one at a time to one. Positive definiteness is preserved. Each determinant derivative is a strictly positive principal cofactor, which is rational and unchanged by conjugation. The irrational coefficient remains fixed. Therefore the reduced rational part satisfies

\[
a\le0\quad(r=3),\qquad a\le-1/16\quad(r=4,5).
\]

For the two reduced cycles, write their determinants \(a-b\sqrt2\) and \(A-B\sqrt2\), with \(b,B>0\). A single bridge of magnitude \(\delta\) gives

\[
\det G=\det L\det M-
\delta^2\det(L-u)\det(M-v).
\]

The correction is rational. The coefficient of \(\sqrt2\) is thus \(-aB-Ab\), which is strictly positive because at least one of \(a,A\) is at most \(-1/16\), and the other is nonpositive. Hence the full reduced determinant cannot vanish. Positive-definite pendant and orthogonal white factors do not change this nonvanishing conclusion.

But a quadrilateral two-face has four black sides and a white support of rank \(n-2\). Its \(n+2\) wall normals lie in a Lorentz vector space of dimension \(n+1\); the complete star Gram must be singular. This gives the contradiction. It would not be legitimate to apply this rank count to an arbitrary smaller subdiagram while silently omitting other connected white nodes; that is why the complete-star topology is part of the hypothesis.

A unique connector path of length at least two is even easier: the original two disjoint Lannér cycles have no direct pairings with one another, so their indefinite spans are orthogonal. They supply two negative directions, impossible in a Lorentz space. The intermediate path does not alter those original pairings.

## Explicit limits and surviving configurations

The review independently reconstructs two singular configurations outside the proved scope:

* Two bare \((3,3,4)\) cycles joined at their white vertices by a 4-bond have determinant zero and four positive-definite F4 corners. Their rational determinant parts are both zero, so the strict irrational-coefficient argument does not apply.
* The six-node path with labels \(8,3,4,3,8\) has determinant zero, checked by the exact tridiagonal recurrence. It uses different minimal indefinite pieces.

These checks do not claim complete polytope realizations. They demonstrate why the current obstruction cannot be reported as a universal quadrilateral classification or a distance-two face-weight theorem. Shared-white configurations and multiple cross contacts require separate analysis.

## Verification

```sh
python3 scripts/fresh_arithmetic_sqrt2_quad_rank_20261015_review.py
python3 scripts/fresh_arithmetic_sqrt2_quad_rank_20261015_review.py --verify
```

Generation and saved replay pass. The primary saved replay, source hashes, syntax and local links are also checked before the final freeze. The review binds primary PY/JSON, not its report. No prior artifact is modified.

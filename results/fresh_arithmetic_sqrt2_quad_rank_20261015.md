# A restricted quadrilateral-star obstruction over Q(√2)

The [exact checker](../scripts/fresh_arithmetic_sqrt2_quad_rank_20261015.py) and [saved certificate](fresh_arithmetic_sqrt2_quad_rank_20261015.json) exclude the complete-star topology specified below. This is a local arithmetic obstruction, not a classification of all quadrilateral stars and not a new global dimension bound. The [independent review](fresh_arithmetic_sqrt2_quad_rank_20261015_review.md) reconstructs the determinant identities and every bounded placement with separate Fraction arithmetic.

**Theorem scope.** Consider an actual quadrilateral two-face of a compact arithmetic Coxeter polytope with ground field \(\mathbb Q(\sqrt2)\). Its white nodes are the full spherical support of the face, and its four black nodes are its four sides. Suppose its affected connected diagram consists of two disjoint chordless cycles of orders 3, 4 or 5, each containing one opposite black pair. Each cycle has exactly one bond 4 and all its other bonds 3. At least one cycle has order 4 or 5. Suppose the cycles are joined by a single white bridge, allowing arbitrary pendant spherical white trees, each attached at precisely one white cycle node. Other white components may be orthogonal. There are no other contacts. Then this cannot be the complete star of that quadrilateral.

A unique connector path of two or more bonds is also impossible: the two original cycle subdiagrams are orthogonal Lannér diagrams and already give two negative directions.

The direct-bridge theorem allows any number and size of the specified pendant trees; this part is proved analytically, not by sampling trees. It does **not** allow shared white nodes, extra black contacts, further cross edges, or both cycles of order 3. These exclusions from scope matter.

## Exact finite placement check

Fixing the distinguished 4-bond on each cycle still covers every labelled placement of its two black nodes and white bridge endpoint. Both possible bridge labels, 3 and 4, are tested. The census contains **4,032 labelled placements**, not 4,032 isomorphism classes.

| Ordered cycle orders | Placements | All four corner diagrams positive definite |
|---|---:|---:|
| 3,4 | 72 | 6 |
| 4,3 | 72 | 6 |
| 3,5 | 180 | 8 |
| 5,3 | 180 | 8 |
| 4,4 | 288 | 0 |
| 4,5 | 720 | 0 |
| 5,4 | 720 | 0 |
| 5,5 | 1,800 | 0 |
| Total | **4,032** | **28** |

All 3,528 placements with both cycle orders at least four already fail a finite-corner test. The other 28 survivors have nonzero full Gram determinant. The argument below also excludes all their allowed pendant-tree extensions. Some survivors have long marked paths, so the full-star rank condition supplies information beyond the four separate corner tests.

The primary uses rational tree elimination for the finite-corner checks and symbolic determinant identities. The independent review uses full-matrix quadratic-field LDL and Gaussian elimination, and derives the cycle formula from determinant matchings. Neither finite-placement census is presented as an exhaustion of arithmetic quadrilateral stars.

## Proof with arbitrary pendant trees

Use unit-square normal vectors. The finite-contact alphabet over this field is \(2,3,4,6,8\), as derived in the [frozen quadratic local package](fresh_arithmetic_sqrt2_local_20261015.md) and its [independent review](fresh_arithmetic_sqrt2_local_20261015_review.md). A bridge or pendant bond belongs to an actual spherical corner component of rank at least three. Finite labels 6 and 8 must be isolated rank-two components, so every such bond here has label 3 or 4. In particular its squared Gram entry is rational. This step is needed: a label-8 square would otherwise introduce a new irrational term in the bridge correction.

Schur-eliminate the pendant white trees. Their positive-definite Gram matrices are forest matrices with rational squared bonds. Successive leaf elimination therefore subtracts rational, nonnegative quantities from the diagonal of the unique attachment node and creates no additional edge. The resulting cycle diagonals \(p_i\) are positive rational numbers at most one. Positivity follows from the original positive-definite white support. Orthogonal white components can be removed as positive factors.

For a cycle of order \(r\in\{3,4,5\}\), its reduced determinant has the form

\[
\det L=a(p)-b_r\sqrt2,
\qquad b_r=2^{1-r}>0.
\]

The two oriented full-cycle terms are the only irrational determinant terms; all matching terms use rational bond squares. Thus the coefficient \(b_r\) is independent of the diagonals. Every proper principal minor is a forest determinant and is rational. Direct calculation gives

\[
a(1,\ldots,1)=0\quad(r=3),
\qquad a(1,\ldots,1)=-\frac1{16}\quad(r=4,5).
\]

At the other arithmetic embedding the full Gram is positive semidefinite. This property survives Schur elimination of the white trees. The conjugate cycle block is actually positive definite: its determinant is \(a(p)+b_r\sqrt2\), which cannot vanish because \(a(p)\) is rational and \(b_r\ne0\). Increasing the reduced diagonals to one preserves positive definiteness. The derivative of the determinant with respect to any diagonal is a positive proper principal cofactor, which is rational and unchanged by conjugation. Consequently

\[
a(p)\le0\quad(r=3),
\qquad a(p)\le-\frac1{16}\quad(r=4,5).
\]

Write the two reduced cycle determinants as \(a-b\sqrt2\) and \(A-B\sqrt2\). A single bridge of magnitude \(\delta\), joining nodes \(u,v\), gives

\[
\det G=(a-b\sqrt2)(A-B\sqrt2)
-\delta^2\det(L-u)\det(M-v).
\]

The second term is rational. The coefficient of \(\sqrt2\) is therefore

\[
-aB-Ab>0,
\]

because at least one cycle has order four or five. Thus the reduced determinant is nonzero, and restoring the eliminated positive-definite factors does not change that conclusion.

But the **complete** star has \(n-2\) independent white nodes and four black nodes in a Lorentz space of dimension \(n+1\). Its Gram determinant must vanish. This contradicts the nonvanishing just proved. Applying this rank count to a smaller subdiagram while omitting other connected white nodes would be invalid; completeness of the specified topology is an explicit hypothesis.

## Configurations deliberately left open

The checker verifies two singular examples outside the theorem:

- The six-node path with labels \(8,3,4,3,8\) has determinant zero. Its long marked corner is of type F4.
- Two bare \((3,3,4)\) cycles joined through their white vertices by a 4-bond have determinant zero and four positive-definite F4 corners. Here both rational determinant parts are zero, so the strict coefficient argument does not apply.

These are local diagrams, not certificates of a complete compact polytope. They suggest testing the weight rule “distance at most two, plus one extra unit on the end pair of an entire F4 component.” The component bound would remain \(2r-2\), but no exhaustive quadrilateral inequality for this rule has been proved here. Shared-white configurations, other Lannér triangles, dotted opposite pairs and additional contacts remain separate tasks.

A further checked follow-up connects the dotted-pair and label-8 cases. Take a dotted pair \(a,y\), a cycle \((b,u,v,w,d,b)\) with its only 4-bond on \(b,d\), and one extra bond \(a,v\) of label 3. All unlisted pairings are right. Give node \(a\) reduced diagonal \(p\), all other diagonals one, and write its pairing with \(y\) as \(-x\). Direct symbolic elimination gives

\[
\det G=\frac{1+\sqrt2}{64}
\bigl(4x^2-4p+1-\sqrt2\bigr).
\]

Attaching a white \(A_\ell\) tail by a 3-bond at the **black** node \(a\) gives \(p=(\ell+2)/(2(\ell+1))\). Hence singularity forces the norm-two squared contact

\[
D=4x^2=1+\frac{2}{\ell+1}+\sqrt2.
\]

For \(\ell=0\), \(D=3+\sqrt2>4\), genuinely dotted. For \(\ell=1\), \(D=2+\sqrt2<4\), exactly label 8; the minimal indefinite piece changes to a \((8,3,2)\) triangle. For \(\ell=2,3\), the forced values \(5/3+\sqrt2\) and \(3/2+\sqrt2\) are not in the allowed finite-contact alphabet and are not algebraic integers. This identity was checked separately with a symbolic seven-by-seven determinant; it is a follow-up lead, not part of the complete placement census or a global exclusion. In particular, attachment at a black node is outside the pendant-white theorem above.

## Reproduction and status

```sh
python3 scripts/fresh_arithmetic_sqrt2_quad_rank_20261015.py
python3 scripts/fresh_arithmetic_sqrt2_quad_rank_20261015.py --verify
python3 scripts/fresh_arithmetic_sqrt2_quad_rank_20261015_review.py --verify
```

Primary generation and saved replay pass. The independent generation, saved replay, source bindings and syntax checks also pass. Source manifests bind the preceding quadratic local proof and its independent review. This package is frozen; no prior frozen artifact is modified.

# Exact projected-wall no-go in three golden OG10 carrier boxes

## Result

Let

\[
F=\langle-2\phi\rangle\oplus E_8\oplus\langle2\rangle,
\qquad
H(\ell)=\begin{pmatrix}F&\ell\\ \ell^t&2\end{pmatrix}.
\]

The last basis vector \(e\) is automatically a norm-two root.  Put
\(p=F^{-1}\ell\).  If \(2p\in F\), then

\[
J|_F=1,\qquad J(e)=2p-e
\]

is an integral isometric involution fixing the old golden \(H^9\).  The
orthogonal Schur coordinate is \(e-p\), so the two required field signatures
are equivalent to

\[
0<p^2<2,\qquad (p')^2<2.
\]

This turns the geometric design problem around: search for the desired new
divisor wall first, rather than search ambient Gram matrices and hope that a
useful wall appears.

Write

\[
\ell=(A,E_8w,C),\qquad A,C\in\mathbf Z[\phi],\qquad
w\in E_8\otimes\mathbf Z[\phi].
\]

Then \(2p=(-A/\phi,2w,C)\) is integral automatically and

\[
q:=2p^2=-\frac{A^2}{\phi}+C^2+2w^2.
\]

Thus the exact signature gate is \(0<q<4\) and \(q'<4\).  A wall cuts off the
archived 61-wall escape precisely when its exact pairings with the archived
interior point and the ray \(v_{61}\) have opposite signs.

## The zero-E8 rank-three case

For \(w=0\), with all four coefficients of \(A,C\) in \([-30,30]\), exactly
66 pairs survive the two-place signature gate.  None separates the interior
from \(v_{61}\).  This case is screened separately: it is not silently omitted
from either nonzero-E8 family.

## Box A: scalar-unit E8 directions

The first box is

\[
w=\phi^n v,\qquad 1\le n\le8,\qquad
v\in E_8,\quad v^2\in\{2,4,6,8\},
\]

with all four integral coefficients of \(A=a+b\phi\) and \(C=c+d\phi\) in
\([-30,30]\).  The exact \(E_8\) shell census is

| square | vectors |
|---:|---:|
| 2 | 240 |
| 4 | 2,160 |
| 6 | 6,720 |
| 8 | 17,520 |

The enumerator generates only genuine \(E_8\) vectors and pins completeness
against the theta-series formula \(240\sigma_3(m)\) in square \(2m\).

Across the 32 power/shell choices, only 480 \((A,C)\) pairs survive the exact
two-place signature gate.  Expanding them over every vector in the relevant
shell gives 2,706,240 exact wall tests.  None separates the interior from
\(v_{61}\).

## Box B: the first genuinely mixed E8 shell

The second box is not confined to scalar multiples:

\[
w=a+\phi b,\qquad a,b\in E_8,\qquad b^2=2,\qquad 0<(w')^2<2,
\]

with the same coefficient box \([-30,30]\) for \(A,C\).  For each of the 240
choices of \(b\), all possible \(a\) are enumerated by exact rational Schur
complements of the positive-definite \(E_8\) form.  This gives exactly 13,920
eligible mixed vectors in three square classes.  The signature gate leaves 64
parameter pairs summed over those classes, hence 15,360 exact wall tests.
Again none cuts \(v_{61}\).

Here \(b\ne0\) by definition.  If \(b=0\), a nonzero \(a\in E_8\) has
\((w')^2=a^2\ge2\), contradicting the strict conjugate signature gate.  The
only remaining \(b=0\) possibility is therefore \(w=0\), already exhausted
in the preceding rank-three case.

## Meaning and scope

This is a finite-box theorem, not global non-existence.  It nevertheless
eliminates the most rigid rank-11 extensions suggested by the failed coupled
\(A_2\) carrier, including the first genuinely mixed golden \(E_8\) shell.
The no-go occurs before asking whether the rank-22 trace lattice embeds
primitively in the OG10 Beauville lattice.  Consequently discriminant-form or
Kähler-wall filters cannot recover a candidate inside these boxes.

The next meaningful search should enlarge the mixed coefficient \(b\)-shell,
or relax the one-root reflection involution ansatz; merely increasing the
scalar-unit power is unlikely to change the geometry efficiently.

The formerly exploratory run through \(b^2\le4\) has since been promoted in
a separate certificate.  Exact Schur recursion and a no-project-import
Cartesian-shell replay compare all 186,720 eligible mixed vectors and all
827,520 wall records literally and again find no cut.  That later theorem is
not folded retroactively into this certificate; see the
[mixed-shell-two no-go](og10_golden_projected_wall_mixed_shell2_no_go.md).

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/search_og10_golden_projected_wall.py \
  --coefficient-bound 30 --unit-power-bound 8 --e8-shell-bound 4

PYTHONPATH=src:. python3 scripts/search_og10_golden_projected_wall.py \
  --coefficient-bound 30 --mixed-e8-b-shell-bound 1

PYTHONPATH=src:. python3 scripts/search_og10_golden_projected_wall.py \
  --write-certificate

PYTHONPATH=src:. python3 scripts/search_og10_golden_projected_wall.py \
  --verify-certificate
```

The verifier independently recomputes all three full boxes and compares
canonical wall-record digests, the combined audit digest, the certificate
payload digest, and all pinned source-file hashes.  Frozen data are in
`data/candidates/og10_golden_projected_wall_box_no_go.json`.

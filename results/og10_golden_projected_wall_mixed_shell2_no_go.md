# Exact mixed-E8 projected-wall no-go through \(b^2\leq4\)

## The finite-box theorem

Keep the golden OG10 fixed carrier

\[
F=\langle-2\phi\rangle\oplus E_8\oplus\langle2\rangle,
\qquad
H(\ell)=\begin{pmatrix}F&\ell\\ \ell^t&2\end{pmatrix},
\]

and write

\[
\ell=(A,E_8w,C),\qquad w=a+\phi b,qquad a,b\in E_8.
\]

Let all four rational-integer coefficients of \(A,C\in\mathbf Z[\phi]\)
lie in \([-30,30]\).  In the mixed family impose

\[
b^2\in\{2,4\},\qquad 0<(w')^2<2.
\]

Then **no** such covector which passes the two-place Hodge-signature gate

\[
0<q=2p^2=-A^2/\phi+C^2+2w^2<4,
\qquad q'<4,
\]

separates the archived invariant interior point from the archived hyperideal
escape ray \(v_{61}\).  Thus this exact finite box supplies no projected wall
which repairs the 61-wall chamber.

This promotes the former exploratory \(b^2\leq4\) computation to a theorem.
It does not alter or overwrite the earlier certificate for the first shell
\(b^2=2\).

## Exact census

The complete census is

| quantity | exact value |
|---|---:|
| \(b\)-vectors with \(b^2=2\) or \(4\) | 2,400 |
| independently bounded \(a\)-vectors tested | 9,121 |
| raw independent \((a,b)\) pairs | 21,890,400 |
| eligible \(w=a+\phi b\) | 186,720 |
| eligible \(w^2\) classes | 6 |
| signature-admissible \((A,C)\) counts, summed over the six classes | 128 |
| exact wall tests | 827,520 |
| strict separators | **0** |

The eligible vectors split by \(b^2\) as

\[
13{,}920\quad(b^2=2),\qquad
172{,}800\quad(b^2=4),
\]

and by \(a^2\) as

\[
2{,}400\quad(a^2=0),\qquad
43{,}920\quad(a^2=2),\qquad
140{,}400\quad(a^2=4),\qquad
0\quad(a^2=6).
\]

For all 827,520 wall records, the exact signs on `(interior,v61)` are

\[
413{,}760\ (-,-),\qquad 413{,}760\ (+,+).
\]

There are no opposite-sign pairs and no zero evaluations.  Thus the zero-cut
conclusion is a literal sign census, not a floating-point inference.

## Why the production recursion is complete

Put \(y=a+\phi'b\).  Once the first \(k\) coordinates of \(y\) are fixed,
write the positive-definite \(E_8\) Gram matrix in the corresponding block
form

\[
G=\begin{pmatrix}G_{II}&G_{IJ}\\G_{JI}&G_{JJ}\end{pmatrix}.
\]

The minimum of \(y^tGy\) over every real suffix is exactly

\[
y_I^t\left(G_{II}-G_{IJ}G_{JJ}^{-1}G_{JI}\right)y_I.
\]

The production program computes every one of these Schur forms over
\(\mathbf Q\), clears a positive denominator, and compares the resulting
element of \(\mathbf Z[\phi]\) with 2 by exact integer arithmetic.  A prefix
is removed only when even its real-suffix minimum is at least 2.

The coordinate ranges are also necessary, not heuristic.  Positivity gives

\[
|y_i|^2<2(G^{-1})_{ii},
\]

and \(|\phi'|<1\) converts this into a conservative finite interval for the
integer \(a_i\).  Every leaf is checked again by the exact full quadratic
form.  Finally, the \(E_8\) shells for \(b\) contain respectively 240 and
2,160 vectors, pinned to the theta-series formula
\(240\sigma_3(m)\) in square \(2m\).  These facts prove necessity and
completeness of the production recursion for the declared box.

## Independent replay without translated Schur pruning

The independent program imports no project code and does not run that
translated-sphere recursion.  From

\[
a=w'-\phi'b
\]

and the parallelogram inequality it obtains

\[
a^2\leq2(w')^2+2(\phi')^2b^2
   <4+8(\phi')^2
   =20-8\phi<8.
\]

Because \(E_8\) is even, it is therefore necessary that

\[
a^2\in\{0,2,4,6\}.
\]

It constructs the square-2, square-4, and square-6 shells independently as
simple-reflection Weyl orbits.  Their sizes `240`, `2160`, and `6720` equal
the complete theta-series coefficients.  It then filters all
\(9{,}121\cdot2{,}400=21{,}890{,}400\) pairs using exact
\(\mathbf Z[\phi]\) signs.

The two implementations compare all 186,720 eligible records literally,
including \(a,b,w^2\), the full \(E_8\) covector, and both evaluations.  They
also compare the complete coefficient list and all 827,520 derived wall
records literally.  These objects are equal before any digest is consulted.
The checksums then provide compact provenance:

- eligible records: `1f139312ef48ed4d97b2a6ab3fceddf594b48ad1a884ab842718ef6704fdbf36`;
- wall records: `20184b175d5b937270f6a2da8d8336ded7e42427b785aa43f7121d532c4c8312`.

Here “checksum” means a checksum after literal equality; it is not a
substitute for semantic comparison.

## Scope and next dual move

This is only the declared mixed-E8 finite-box no-go.  It says nothing about
\(b^2>4\), coefficients outside \([-30,30]\), another fixed carrier, or a
non-involutive extension.  A later theorem separately shows that saturated
central-inversion rank-two norm-two coupling creates no new projected-wall
type; that global statement depends on its explicit hypotheses and is not a
consequence of this finite census.
It is neither a global obstruction for OG10 nor a hyperbolic Coxeter example.
Because no projected wall passes the cut gate, this shell reaches neither a
new trace-lattice discriminant calculation nor a primitive OG10 embedding or
MBM/effectivity gate.

Within this carrier, the next disjoint shell is \(b^2=6\).  The structurally
different geometric move is a square-six field root whose selected trace lift
has ambient OG10 divisibility three; merely coupling a second norm-two root is
redundant under the saturated central-inversion theorem.  See
`results/og10_golden_rank_two_reduction_no_go.md`.

## Reproduction

```bash
PYTHONPATH=src:. python3 \
  scripts/certify_og10_golden_projected_wall_mixed_shell2.py --recompute

python3 scripts/independent_og10_golden_projected_wall_mixed_shell2.py \
  --compare-certificate \
  data/candidates/og10_golden_projected_wall_mixed_shell2_no_go.json

PYTHONPATH=src:. python3 \
  scripts/certify_og10_golden_projected_wall_mixed_shell2.py \
  --verify-certificate
```

The frozen certificate is
`data/candidates/og10_golden_projected_wall_mixed_shell2_no_go.json`.

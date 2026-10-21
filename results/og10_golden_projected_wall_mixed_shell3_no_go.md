# Exact mixed-E8 projected-wall no-go through \(b^2\leq6\)

## Theorem

For the fixed golden OG10 carrier

\[
F=\langle-2\phi\rangle\oplus E_8\oplus\langle2\rangle,
\qquad
H(\ell)=\begin{pmatrix}F&\ell\\ \ell^t&2\end{pmatrix},
\]

write

\[
\ell=(A,E_8w,C),\qquad w=a+\phi b,qquad a,b\in E_8.
\]

Let every rational-integer coefficient of \(A,C\in\mathbf Z[\phi]\) lie in
\([-30,30]\), and impose

\[
b^2\in\{2,4,6\},\qquad 0<(w')^2<2.
\]

Among all such covectors passing the exact two-place signature gate

\[
0<q=2p^2=-A^2/\phi+C^2+2w^2<4,
\qquad q'<4,
\]

none separates the archived invariant interior point from the archived
hyperideal escape ray \(v_{61}\).  Thus the complete mixed-E8 shell through
\(b^2=6\) supplies no projected divisor wall repairing the 61-wall chamber.

This is a new shell-3 certificate.  It does not modify the earlier
\(b^2=2\) or \(b^2\leq4\) theorem artifacts.

## Exact census

| quantity | exact value |
|---|---:|
| \(b\)-vectors with \(b^2=2,4,6\) | 9,120 |
| independently necessary \(a\)-vectors | 26,641 |
| raw \((a,b)\) pairs before Weyl quotienting | 242,965,920 |
| representative-pair tests after exact Weyl quotienting | 79,923 |
| eligible \(w=a+\phi b\) | 569,760 |
| eligible \(w^2\) classes | 8 |
| signature-admissible \((A,C)\) counts summed over those classes | 144 |
| exact wall records | 1,472,640 |
| strict separators | **0** |

The eligible-vector distributions are

\[
\begin{array}{c|rrrrr}
a^2&0&2&4&6&8\\ \hline
\#&2{,}400&238{,}800&321{,}840&6{,}720&0
\end{array}
\]

and

\[
\begin{array}{c|rrr}
b^2&2&4&6\\ \hline
\#&13{,}920&172{,}800&383{,}040.
\end{array}
\]

Writing a golden integer as its pair of coefficients, the eight values of
\(2w^2\), their multiplicities, and their admissible coefficient counts are

| \(2w^2\) | eligible \(w\) | admissible \((A,C)\) |
|---:|---:|---:|
| `(4,4)` | 240 | 40 |
| `(8,8)` | 15,600 | 0 |
| `(8,12)` | 240 | 24 |
| `(12,16)` | 30,240 | 24 |
| `(16,20)` | 319,680 | 0 |
| `(16,24)` | 15,600 | 40 |
| `(20,28)` | 181,440 | 0 |
| `(24,36)` | 6,720 | 16 |

Every one of the 1,472,640 resulting walls is evaluated exactly.  The sign
census on `(interior,v61)` is

\[
736{,}320\ (-,-),\qquad736{,}320\ (+,+).
\]

There are no opposite signs and no zero evaluations.

## Production completeness

For each of the theta-pinned 9,120 choices of \(b\), production enumerates
every integral \(a\) with \((a+\phi'b)^2<2\).  Given a fixed coordinate
prefix \(y_I\) of \(y=a+\phi'b\), its pruning quantity is the exact minimum
over all real suffixes:

\[
y_I^t\left(G_{II}-G_{IJ}G_{JJ}^{-1}G_{JI}\right)y_I.
\]

All Schur matrices are rational, positive denominators are cleared, and every
comparison in \(\mathbf Z[\phi]\) is exact.  The coordinate boxes follow from

\[
|y_i|^2<2(G^{-1})_{ii}
\]

and \(|\phi'|<1\); they are conservative necessary bounds.  Every leaf is
then rechecked against the full quadratic form.  The \(b\)-shell sizes
`240`, `2160`, and `6720` equal the exact theta-series coefficients
\(240\sigma_3(m)\) in square \(2m\).

## Structurally independent Weyl replay

The independent replay imports no project code and performs no translated
Schur recursion.  The parallelogram inequality gives

\[
a^2\leq2(w')^2+2(\phi')^2b^2
   <4+12(\phi')^2
   =28-12\phi<10.
\]

Since \(E_8\) is even,

\[
a^2\in\{0,2,4,6,8\}.
\]

The replay constructs these complete shells as simple-reflection Weyl
orbits and pins their totals to the theta series:

| square | Weyl-orbit sizes | shell size |
|---:|---:|---:|
| 2 | 240 | 240 |
| 4 | 2,160 | 2,160 |
| 6 | 6,720 | 6,720 |
| 8 | 17,280 and 240 | 17,520 |

The second square-eight orbit consists of doubled roots; omitting it would
make the independent bound incomplete.

Each \(b\)-shell of square 2, 4, or 6 is a single Weyl orbit.  At one seed in
each shell, exhaustive comparison with all 26,641 possible \(a\)'s gives
respectively 58, 80, and 57 eligible vectors.  The replay records an explicit
integral Weyl matrix from the seed to every \(b\), transports those complete
eligible lists, and checks the exact leaf inequality again after every
transport.  Since a Weyl matrix is a lattice bijection preserving norms and
inner products, this enumerates the same family as all 242,965,920 raw pairs,
without using production's pruning.

Production and replay then compare literally:

- all 569,760 eligible records, including \(a,b,w^2\), the full covector and
  both evaluations;
- all 46,899 coefficient records; and
- all 1,472,640 derived wall records.

Only after literal equality are the compact checksums recorded:

- eligible records: `f3437fc0d251f3c0d0fb0729bb3e317e72676ba52988002b1f4c1aa759608ecf`;
- wall records: `851270ff6c0ab5ea5afefb7e7890e957c20ca31bac6fe2dcfa3d3114fe08ad73`.

## Scope

This is a finite-box theorem only.  It does not cover \(b^2>6\), coefficients
outside \([-30,30]\), other carriers, a relaxed reflection involution, or a
coupled multi-root extension.  It is neither a global OG10 obstruction nor a
new hyperbolic Coxeter example.  Because no wall survives the cut gate, no
new trace-lattice discriminant, primitive OG10 embedding, or MBM/effectivity
gate is reached.

The zero-\(E_8\) subfamily has separately been closed globally, with no
coefficient bound; see
`results/og10_golden_projected_wall_zero_e8_global_no_go.md`.  That theorem
does not make the mixed-\(E_8\) statement here global.

The next disjoint shell of this same rank-one search is \(b^2=8\), whose
independent necessity bound would allow larger \(a\)-shells.  A coupled
rank-two norm-two carrier is **not** an open alternative: the separate exact
rank-two reduction theorem proves it exceeds the OG10 Picard-rank bound and,
abstractly, adds no projected norm-two wall type under the natural
involution.  The structurally different remaining OG10 direction is instead
a primitive norm-six, divisibility-three wall outside its already excluded
finite boxes, together with the required primitive gluing; otherwise one
must change the carrier or variety.

## Reproduction

```bash
PYTHONPATH=src:. python3 \
  scripts/certify_og10_golden_projected_wall_mixed_shell3.py --recompute

python3 scripts/independent_og10_golden_projected_wall_mixed_shell3.py \
  --compare-certificate \
  data/candidates/og10_golden_projected_wall_mixed_shell3_no_go.json

PYTHONPATH=src:. python3 \
  scripts/certify_og10_golden_projected_wall_mixed_shell3.py \
  --verify-certificate
```

The frozen certificate is
`data/candidates/og10_golden_projected_wall_mixed_shell3_no_go.json`.

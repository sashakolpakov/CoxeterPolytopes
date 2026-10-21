# Independent review of the Q(sqrt2) weight accounting

The proposed weights and their **conditional dimension-20 implication pass** an
independent reconstruction. The separate global proof must still show that every
quadrilateral, including all external white components, reduces to the listed
entire-component configurations. This review does not establish that reduction.

The [standalone checker](../scripts/fresh_sqrt2_weight_bound_20261015_review.py)
imports no producer, project arithmetic, NetworkX, or graph-isomorphism routine.
It uses breadth-first distances, an intrinsic arm-length classifier, and exact
Fraction arithmetic in Q(sqrt2). The [certificate](fresh_sqrt2_weight_bound_20261015_review.json)
records all literal corner diagrams, finite Schur pivots, weights, component
totals, partition maxima, and dimension inequalities.

## Literal corner reconstruction

All seven selected singular finite-label cores are reconstructed from the raw
labelled union lists. Their four actual finite-corner graphs are tested by exact
leaf elimination: after normalizing diagonal entries to 2, only the squares of
the edge entries are needed, so the arithmetic remains in Q(sqrt2), including
label 8. This is an independent positive-definiteness check, not only a count of
graph distances.

For the three pendant survivors, the checker starts with the saved **scaled
physical Gram matrices**, computes every corner's normalized Cartan cycles, and
recovers the finite labels. It does not reuse the primary weight checker's
handwritten tail graphs. All twelve tail corners agree, with quadrilateral
totals **1, 1, 2**. Each of the seven finite-label cores likewise has total at
least 1, with all four individual weights matching the primary certificate.

The corrections are intrinsic to complete connected components. Equal-length
arms receive symmetric assignments; the E7/E8 distinguished arms have different
lengths; and the F4 endpoints are intrinsic. Thus the assignments are invariant
under all diagram automorphisms. The code also verifies a relabelling guard.
No correction is assigned merely because a proper subdiagram resembles one of
the exceptional types. This entire-component requirement is essential to the
separate classification dependency.

## All ranks and triangular faces

The base distance-two pair count is \(2r-3\) on a path of rank \(r\ge2\), and
\(2r-2\) on a D/E tree. These formulas follow directly by counting edges and
unordered neighbor pairs. A rank-one component contributes zero. After the
corrections, the only component deficits from \(2r\) below 2 are:

| Component | Total weight | Deficit |
|---|---:|---:|
| D5 | 9 | 1 |
| E6 | 11 | 1 |
| E7 | 49/4 | 7/4 |
| E8 | 57/4 | 7/4 |

F4 has weight 6 and deficit 2. Every other finite component has deficit at least
2. Consequently, at total rank \(n\ge9\), a single component has deficit at
least 2; two or more components also have combined deficit at least 2, because
each component has deficit at least 1. This proves the all-rank vertex bound
\(2n-2\). The additional exact partition calculation through rank 100 attains
that bound for every \(n\ge9\); it is a check of the proof, not a finite cutoff
standing in for it.

For a compact triangular face, a minimal indefinite diagram in its full star
must contain all three marked black nodes: omitting any one would put it in a
finite corner. The separately supplied Q(sqrt2) Lannér catalogue restricts this
diagram to a connected triangle or a cycle of length four or five. The connected
triangle supplies at least two nonright edges. In a four- or five-cycle, the
three positive cyclic gaps between black nodes sum to at most five, so at least
two gaps have length at most two. Each of these paths avoids the third black
node and is present in the relevant finite corner. These give two base weights,
regardless of any additional nodes. Hence every triangle has weight at least 2.
All cycle markings are independently checked. The catalogue's global
completeness remains a named input, rather than an inference from those markings.

## Strict dimension inequality

For a simple compact n-polytope, let \(h_0,\ldots,h_n\) be the symmetric
nonnegative h-vector of its simplicial dual. Then
\[
f_0=\sum h_i,\qquad
f_2=\sum_i\binom{n-i}{2}h_i.
\]
After symmetry, the coefficient of \(h_i\) in \(f_2\) is
\(\{\binom i2+\binom{n-i}2\}/2\). Its minimum is
\(n(n-2)/8\) for even n and \((n-1)^2/8\) for odd n. Since
\(h_0=h_n=1\) and the two endpoint coefficients are strictly larger than the
minimum, \(f_2/f_0\) is **strictly** larger than this minimum.

Conditional on quadrilateral weight at least 1, every k-gonal two-face has
weight at least \(5-k\): triangles and quadrilaterals give the nontrivial cases,
while nonnegative weights suffice for \(k\ge5\). Double-counting angles gives
\[
\frac{\text{total vertex weight}}{f_0}
\ge5\frac{f_2}{f_0}-\binom n2
>
\begin{cases}n(n-6)/8,&n\text{ even},\\
(n-1)(n-5)/8,&n\text{ odd}.
\end{cases}
\]
For \(n=21+2u\) and \(n=22+2u\), subtraction of \(2n-2\) gives respectively
\(u^2/2+5u\) and \(u^2/2+11u/2+2\). Both are nonnegative for every integer
\(u\ge0\). The first is zero at n=21, where the strict h-vector inequality
still contradicts the vertex upper bound 40. Thus the conditional conclusion
is n at most 20, with no parity or equality gap.

The finite spherical diagram classification, h-vector facts, and the separately
supplied Lannér catalogue are the mathematical inputs. Exhaustion of all
quadrilateral stars and arbitrary external attachments is outside this review.
The arithmetic and accounting do not by themselves promote a bounded local
census into an unconditional theorem.

Generation and saved replay pass against the final primary dependency hashes.
Source manifests, syntax, and local links pass. No mathematical discrepancy was found.

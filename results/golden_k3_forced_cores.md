# Forced crystallographic cores in a golden K3 9/14 candidate

## Result

Consider a compact Coxeter polytope in \(\mathbb H^9\) with 14 facets whose
walls are represented by norm-two roots of an even
\(\mathbb Z[\phi]\)-lattice \((M,h)\), with signatures \((9,1)\) and
\((10,0)\) at the two real places.  Let

\[
N=\left(M,\operatorname{Tr}_{\mathbb Q(\sqrt5)/\mathbb Q}
       \left(-\frac{\phi}{\sqrt5}h\right)\right)
\]

be its rank-20 K3 trace lattice.

Then the Coxeter diagram contains one of the following finite golden cores:

\[
H_3\quad\text{or}\quad H_4.
\]

Their full \(\mathbb Z[\phi]\)-spans have the exact trace lifts

\[
H_3\longmapsto D_6(-1),
\qquad
H_4\longmapsto E_8(-1).
\]

The \(D_6(-1)\) sublattice is primitive in the even lattice \(N\).  In the
\(H_4\) case, \(E_8(-1)\) is unimodular and hence splits orthogonally:

\[
N\cong E_8(-1)\oplus N_{12},
\qquad \operatorname{sig}(N_{12})=(1,11).
\]

Equivalently, already over \(\mathbb Z[\phi]\),

\[
M\cong H_4\perp M_6,
\]

because the determinant \(5-3\phi=\phi^{-4}\) of the norm-two \(H_4\)
Gram matrix is a unit.  Thus the missing-face combinatorics can force an
entire \(E_8(-1)\) divisor summand before a K3 surface is constructed.

This is a structural reduction, not an existence or non-existence theorem for
the full polytope.

## Proof of the forced-core dichotomy

Burcroff proves that a hypothetical compact 9-polytope with 14 facets has no
missing face of order three and no Coxeter label \(m\geq6\).  On the other
hand, Esselmann's lemma says that a 3-free compact Coxeter 9-polytope has at
least \(18\) facets.  Therefore a 14-facet candidate is not 3-free.  Since a
Lanner diagram has order at most five, the candidate has a missing face of
order four or five.

The finite subdiagram used in the standard missing-face reduction is

\[
\begin{array}{c|c}
\text{missing-face order}&\text{finite subdiagram}\
\hline
4&H_3\text{ or }B_3,\\
5&H_4\text{ or }F_4.
\end{array}
\]

With all root norms normalized to two, a label-four edge would require field
product \(-\sqrt2\), which does not belong to \(\mathbb Q(\sqrt5)\).  Hence
the \(B_3\) and \(F_4\) alternatives are unavailable: order four forces
\(H_3\), and order five forces \(H_4\).

The scope of the norm-two hypothesis matters.  This argument does not exclude
label four for an unrelated unequal-root-length construction.

## Exact trace calculation

For the path diagrams

\[
H_2=[5],\qquad H_3=[5,3],\qquad H_4=[5,3,3],
\]

use diagonal entries \(2\), a 3-edge product \(-1\), and a 5-edge product
\(-\phi\).  In the interleaved integral basis

\[
e_1,\phi e_1,\ldots,e_r,\phi e_r,
\]

entrywise restriction of scalars gives even negative-definite lattices with
determinants

| golden type | field determinant | trace rank | trace determinant | exact trace type |
|---|---:|---:|---:|---|
| \(H_2\) | \(3-\phi\) | 4 | 5 | \(A_4(-1)\) |
| \(H_3\) | \(4-2\phi\) | 6 | 4 | \(D_6(-1)\) |
| \(H_4\) | \(5-3\phi\) | 8 | 1 | \(E_8(-1)\) |

The implementation contains unimodular integral matrices \(S\) satisfying

\[
Q_{\mathrm{trace}}=S^tQ_{\mathrm{standard}}S
\]

for all three rows.  Thus the Dynkin identifications are independently
checkable integral isometries, not determinant-based guesses.

For \(D_6\), the three nonzero discriminant classes have quadratic values
represented by \(1,3/2,3/2\), up to sign.  None is isotropic modulo
\(2\mathbb Z\), so \(D_6(-1)\) has no proper even overlattice and is primitive
inside \(N\).  For \(E_8(-1)\), unimodularity proves both primitivity and the
orthogonal splitting asserted above.

## Consequences for inverse design

The first search branch is now forced by the nerve rather than guessed from a
surface:

1. If a candidate nerve has a missing 5-face, choose the forced \(H_4\) and
   solve only for the rank-six golden complement \(M_6\).  On the K3 side,
   test the rank-12 lattice \(N_{12}\) inside
   \(U^3\oplus E_8(-1)\); its discriminant group must still have length at
   most two because the transcendental complement has rank two.
2. If all non-pair missing faces have order four, start with the primitive
   \(D_6(-1)\) trace core and enumerate its finite discriminant-form gluings
   to a rank-14 complement.
3. Only after this forced-core gate should the remaining two-Gram entries,
   dotted products, rank-ten Gram equations, and Vinberg completeness be
   solved.

An order-five missing face satisfies an additional published restriction: its
\(H_4\) has at least two bad neighbours.  The next local calculation should
classify the possible projections of those neighbours to the rank-six
orthogonal complement.  For the neighbour that completes the \(H_4\) to its
order-five Lanner missing face, this calculation is already finite.  Testing
all \(3^4=81\) attachments with labels \(2,3,5\) leaves exactly

| attachment to the ordered \([5,3,3]\) path | field determinant | norm in \(M_6\) | visible curve type |
|---|---:|---:|---|
| \((2,2,2,3)\) | \(6-4\phi\) | \(-2\phi\) | \(A_5(-1)\) |
| \((2,2,2,5)\) | \(8-6\phi\) | \(-2-6\phi\) | \(A_5(-1)\) |
| \((2,2,3,2)\) | \(4-4\phi\) | \(-4-8\phi\) | \(D_5(-1)\) |

Here the complement norm is the exact Schur complement, or equivalently the
full determinant divided by \(5-3\phi\).  Consequently the rank-six golden
complement is not arbitrary: it must integrally represent at least one of
these three elements, negative at the identity place and positive at the
conjugate place.  The second bad neighbour may involve a dotted edge and needs
a separate integral-pairing classification.  This is still substantially
smaller than searching arbitrary rank-ten forms.

## Reproducibility

The certificates are implemented in
`coxeter_polytopes.geometric.golden_finite` and tested by
`tests/test_golden_finite_trace_lifts.py`.  No floating-point comparison is
used.

## Sources and scope

The combinatorial inputs are Burcroff, Lemma 10.10 and the missing-face
reduction preceding Theorem 10.7, together with Esselmann's Lemma 6.7 as quoted
there.  The statement applies to the norm-two golden K3 factory.  It neither
enumerates every 9/14 nerve nor proves that a surviving lattice has a complete
or compact eigen-nef chamber.

# Complete arithmetic Lannér input over Q(√2)

The [standalone exact checker](../scripts/fresh_sqrt2_lanner_catalogue_20261015.py) and [certificate](fresh_sqrt2_lanner_catalogue_20261015.json) give the complete Lannér catalogue compatible with a compact arithmetic hyperbolic Coxeter chamber whose ground field is \(K=\mathbb Q(\sqrt2)\). There is one parametric order-two family, nine triangle types, and one type each of orders four and five. This is classification input for a separate quadrilateral analysis; it proves no new dimension bound.

All Gram matrices here have diagonal two. A finite bond \(m\) has negative entry \(-2\cos(\pi/m)\). Dotted distances are not cut off or discretely sampled.

## Published inputs and normalization

Vinberg's arithmetic criterion requires the norm-two Gram's cyclic products to be algebraic integers in the ground field, its entry field to be totally real, and its conjugate matrices away from the physical ground-field embedding to be positive semidefinite. The criterion is recalled in [Nikulin, §2](https://arxiv.org/pdf/0708.3991). A nonsingular principal block is therefore positive definite at such an embedding. The classification of Lannér diagrams, their maximum order five, and their role as minimal indefinite subdiagrams are given in [Vinberg, §2.4, printed pages 74–75](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt). These are published inputs, not conclusions of the finite enumeration below.

The square of a finite contact is \(2+2\cos(2\pi/m)\in K\). Thus the real cyclotomic degree is at most two. The complete candidates are \(2,3,4,5,6,8,10,12\): the totient condition permits only primes 2, 3 and 5 with bounded exponents. Labels 5 and 10 generate \(\mathbb Q(\sqrt5)\), and 12 generates \(\mathbb Q(\sqrt3)\). The complete finite alphabet is consequently

\[
\boxed{2,3,4,6,8}.
\]

This is an all-label field argument, not an enumeration up to a chosen maximum label.

## Order two: a complete parametric family

Write the Gram as

\[
G=\begin{pmatrix}2&-\sqrt t\\-\sqrt t&2\end{pmatrix}.
\]

Its complete arithmetic conditions are

\[
\boxed{t=a+b\sqrt2\in\mathbb Z[\sqrt2],\qquad t>4,
\qquad 0<t'=a-b\sqrt2<4.}
\]

The first strict inequality is physical indefiniteness; integrality is the two-cycle condition. Total reality requires \(t'>0\), and conjugate positive definiteness requires \(t'<4\). The boundary values 0 and 4 cannot occur under a field embedding of a nonzero \(t\) or \(t-4\). Conversely these conditions give a totally real radical extension and the required two-cycle and conjugate-sign tests. The example \(t=3+\sqrt2\) is checked exactly. There is no finite catalogue of dotted magnitudes asserted.

## Order three: exactly nine triangle types

All proper two-node subdiagrams are finite, so the finite alphabet applies to every bond. The triangle is physically hyperbolic exactly when \(1/p+1/q+1/r<1\). If a right angle occurs, the three-cycle vanishes. Otherwise its product must belong to \(K\).

The magnitudes for labels 3 and 4 already lie in \(K\). The other two are \(\sqrt3\) and \(\sqrt{2+\sqrt2}\), whose squareclasses are independent over \(K\): 3 is not a square in \(K\), and the norms of \(2+\sqrt2\) and \(3(2+\sqrt2)\) are 2 and 18, both nonsquares in \(\mathbb Q\). Hence each of labels 6 and 8 must occur an even number of times in a nonzero three-cycle. This exact parity condition tests field membership without symbolic numerical recognition.

For magnitudes \(a,b,c\), the determinant is

\[
\Delta=8-2(a^2+b^2+c^2)-2abc.
\]

Each conjugate proper two-node block is positive definite. Therefore the remaining arithmetic sign test is exactly \(\Delta'>0\). All cyclic products that pass the field test are integral. The required radical fields are totally real because \(3>0\) and \(2\pm\sqrt2>0\).

| Unordered labels | \(\Delta\) | \(\Delta'\) |
|---|---|---|
| 2, 3, 8 | \(2-2\sqrt2\) | \(2+2\sqrt2\) |
| 2, 4, 8 | \(-2\sqrt2\) | \(2\sqrt2\) |
| 2, 6, 8 | \(-2-2\sqrt2\) | \(-2+2\sqrt2\) |
| 2, 8, 8 | \(-4\sqrt2\) | \(4\sqrt2\) |
| 3, 3, 4 | \(-2\sqrt2\) | \(2\sqrt2\) |
| 3, 8, 8 | \(-6-6\sqrt2\) | \(-6+6\sqrt2\) |
| 4, 4, 4 | \(-4-4\sqrt2\) | \(-4+4\sqrt2\) |
| 4, 6, 6 | \(-8-6\sqrt2\) | \(-8+6\sqrt2\) |
| 4, 8, 8 | \(-8-8\sqrt2\) | \(-8+8\sqrt2\) |

The certificate records all 35 unordered triples from the proved alphabet and the exact reason for every rejection.

## Orders four and five: one cycle each

An edge of label at least 6 in a connected diagram has an adjacent nonright edge. Their three-node Gram has determinant no greater than the relaxed unit-square value \(1-\cos^2(\pi/6)-\cos^2(\pi/3)=0\). Any additional bond only decreases the determinant. For Lannér order at least four this is a proper subdiagram, contradicting positive definiteness. Parallel and dotted bonds likewise contradict proper two-node positivity. With the field alphabet, only labels 2, 3 and 4 remain.

A tree with these labels cannot be arithmetic over \(K\): conjugation changes the signs precisely on the 4-bonds, and vertex sign changes on a tree undo all those changes. Its conjugate Gram is therefore congruent to its indefinite physical Gram.

If a Lannér graph contains a cycle, that cycle must use all vertices. Indeed the all-3 cycle already has a nonpositive all-ones quadratic form, and larger bonds or additional chords decrease it. A proper cycle would violate proper-subdiagram positivity. There can be no chord, because it produces a shorter cycle. The graph is consequently a chordless cycle.

An even number of 4-bonds again makes the conjugate Gram sign-switch-congruent to the physical Gram. An odd number at least three in a cycle of length four or five forces two adjacent 4-bonds; their proper three-node path is affine. Exactly one 4-bond remains. Both cycles pass all proper-subdiagram and conjugate tests:

| Order | Diagram | Norm-two determinant | Conjugate determinant |
|---|---|---|---|
| 4 | C4 with labels 4,3,3,3 | \(-1-2\sqrt2\) | \(-1+2\sqrt2\) |
| 5 | C5 with labels 4,3,3,3,3 | \(-2-2\sqrt2\) | \(-2+2\sqrt2\) |

As a separate check of this analytic classification, the script tests every labelled graph on four and five nodes after the proved 2/3/4 reduction. It checks every maximal proper principal matrix exactly, then the physical determinant and the conjugate Gram. The 729 and 59,049 graphs give respectively 12 and 60 labelled arithmetic survivors, each one isomorphism class. Six additional labelled order-four Lannér graphs fail the conjugate test. No numerical eigenvalue acceptance or project arithmetic import is used.

## What a complete no-short-corner argument must still handle

For a compact quadrilateral with full white support \(W\) and cyclic black sides \(a,b,c,d\), each opposite block \(W+a+c\) and \(W+b+d\) is indefinite. It contains a minimal Lannér diagram that includes both of its black nodes. Thus the list above is the complete input for those two chosen pieces. It does not classify their union or the complete quadrilateral star.

The frozen [disjoint-cycle bridge theorem](fresh_arithmetic_sqrt2_quad_rank_20261015.md) handles two disjoint single-4 cycles with one direct white bridge and certain pendant white trees, provided at least one cycle has order four or five. Beyond that theorem, a global weighting proof still needs the following reductions:

- One or both opposite pieces can be a dotted pair, including the entire unbounded algebraic family above. All nine triangle types must also be retained, not just the 3,3,4 cycle.
- The selected pieces can share white nodes. Their union can have further cross contacts; one cannot assume disjointness or a unique bridge before a finite-corner argument proves it.
- White vertices outside the chosen union can attach at black nodes or at several parts of the union. Eliminating such a component can change off-diagonal entries as well as diagonals. The pendant-white theorem's rational diagonal correction is not automatically valid for those attachments.
- The complete star, rather than an arbitrarily selected subdiagram, must satisfy the ambient rank defect. Full-rank local Lorentz subdiagrams can sit inside a larger singular star. A determinant calculation on a partial union alone is not a contradiction.
- Every proposed weighted exceptional corner must be the entire spherical component whenever its weight bound depends on that property. The simultaneous four actual-corner inequalities must be checked, not only the selected Lannér blocks.

For example, the already frozen singular path with labels \(8,3,4,3,8\) and the singular bridge of two 3,3,4 cycles show why the cubic-seven quadrilateral argument does not transfer unchanged. Both have relevant F4 corners. They motivate testing an entire-F4 correction, but the present catalogue does not establish that correction's quadrilateral inequality.

The separate quadrilateral support enumeration is owned by the compact-search lane. This package does not duplicate it and makes no exhaustion claim about external attachments.

## Reproduction

```sh
python3 scripts/fresh_sqrt2_lanner_catalogue_20261015.py
python3 scripts/fresh_sqrt2_lanner_catalogue_20261015.py --verify
```

Generation and saved replay pass. Primary literature was checked on 2026-10-15 for the inputs cited above. No novelty or current-best claim is made.

# The H10/16 polytope survives, but its proposed Nikulin lift does not

## Sharp theorem

Let

\[
M=U\oplus D_8(-1)\oplus A_1(-1)^2
\]

and let \(h\in O(M)\) be the second archived high involution. Its fixed
section is the genuine finite-volume noncompact Coxeter polytope
R12-D8A1x2-C2b-H10-16: 16 facets, 70 finite vertices, six cusps, and 426
elliptic rank-nine edges. Every edge has exactly two finite or ideal
endpoints, all 14 full affine zero-neighbourhood audits pass, and the
canonical Vinberg audit has SHA-256

    900df6d000ca91669cea117474cda559e4188600fb2744bafa0f69fd1541c226

The algebraic conclusion is negative but exact:

> The action \(h\) cannot be the restriction of a symplectic Nikulin
> involution on any special K3 surface obtained from a primitive embedding
> \(M\hookrightarrow\Lambda_{K3}\). More generally, no K3 involution of
> either involution type can realize \(h\) while keeping the fixed Picard
> lattice equal to the archived rank-11 lattice. Consequently, any exact
> fixed-rank-11 K3-automorphism realization of this H10/16 section must have
> order at least four.

This does **not** invalidate the hyperbolic polytope. It closes one proposed
algebro-geometric carrier route and turns the next positive search into a
higher-order problem.

## What passes before the obstruction

The archived Picard action has characteristic polynomial

\[
(x-1)^{11}(x+1).
\]

Its primitive coinvariant is \(\langle-4\rangle\) of divisibility two. The
fixed lattice has determinant 16 and Smith factors

    1^8, 2, 2, 4.

The orthogonal complement of a primitive K3-lattice embedding is the unique
indefinite even 2-elementary lattice with signature \((2,8)\) and invariants
\((10,4,1)\):

\[
T=U\oplus U(2)\oplus D_6(-1).
\]

The four-bit discriminant gate is not the problem. The certificate gives an
explicit anti-isometry

\[
P=
\begin{pmatrix}
1&0&1&1\\
0&1&1&1\\
0&0&1&1\\
0&0&1&0
\end{pmatrix}
\]

between \(A_M\) and \(A_T\), checks all 16 quadratic-form classes, and
transports the Picard discriminant action equivariantly. Thus a test based
only on discriminant groups would incorrectly leave the Nikulin route open.

The complement identification uses the classification of indefinite even
2-elementary lattices; see [Alexeev's modern formulation](https://arxiv.org/abs/2209.09110)
and [Nikulin's primitive-embedding/discriminant-form theorem](https://doi.org/10.1070/IM1980v014n01ABEH001060).

## The mod-two obstruction

On \(T/2T\), define

\[
\bar q(x)=\frac{x^2}{2}\pmod 2.
\]

Its polar form is the reduction of the bilinear form of \(T\). The exact
calculation gives:

- \(\dim_{\mathbf F_2}(T/2T)=10\);
- the polar form has rank six and radical dimension four;
- \(\bar q\) is nonzero on that radical, so its singular radical has
  dimension three;
- the nondegenerate symplectic quotient has dimension six and maximal
  isotropic dimension three.

Therefore every totally \(\bar q\)-singular subspace has dimension at most

\[
3+3=6.
\]

The certificate also displays a six-dimensional totally singular subspace,
so this bound is sharp.

Now suppose a primitive \(E_7(-2)\hookrightarrow T\) existed. Primitivity
would make

\[
E_7(-2)/2E_7(-2)\longrightarrow T/2T
\]

injective. Since the ordinary \(E_7\) root lattice is even, \(\bar q\)
vanishes identically on the seven-dimensional source. Its image would be a
seven-dimensional totally singular subspace, contradicting the sharp bound
six. Hence

\[
E_7(-2)\not\hookrightarrow_{\mathrm{primitive}}T.
\]

This also explains a tempting false positive found during the search: an
E7(-2)-shaped Gram embedding exists at index two in its saturation, but it
fails precisely the primitive-kernel gate.

## Why this excludes a Nikulin involution

For a symplectic K3 involution, the full anti-invariant cohomology lattice is
\(E_8(-2)\). The archived primitive coinvariant vector has square \(-4\), so
inside \(E_8(-2)\) it is a scaled \(E_8\) root. Its orthogonal complement is
\(E_7(-2)\). Intersecting with the primitive transcendental complement would
therefore produce exactly the primitive embedding just ruled out.

An order-two nonsymplectic escape with the *same* fixed Picard lattice is also
impossible. The invariant cohomology lattice of a nonsymplectic K3
involution is 2-elementary, while the required fixed lattice has a
\(\mathbf Z/4\) Smith factor. This second statement is deliberately scoped:
it does not exclude constructions that enlarge the fixed Picard lattice and
thereby change the automorphism-fixed hyperbolic section.

## Remaining positive completion target

The next carrier search should not repeat the E7(-2) construction. It must
find an order-four or higher cohomological action satisfying all of the
following exact gates:

1. its restriction to \(M\) is the archived C2b action;
2. its action on \(A_T\) is the transported four-bit action stored in the
   certificate;
3. the period has a non-real eigenvalue and no new invariant Picard classes
   appear, so \(NS(X)^g=M^h\) still has rank 11;
4. every added algebraic direction is nonfixed and the enlarged Picard
   coinvariant is rootless;
5. an invariant Kähler class lies in the archived strict fixed-section
   interior;
6. the already frozen complete Vinberg audit is then inherited unchanged.

This is a finite lattice-and-period criterion for a positive result, not a
global no-go for higher-order automorphisms.

## Reproducibility

The exact certificate is generated by
[certify_r12_h10_16_nikulin_obstruction.py](../scripts/certify_r12_h10_16_nikulin_obstruction.py)
and archived at
[r12_h10_16_nikulin_obstruction.json](../data/enumerations/r12_h10_16_nikulin_obstruction.json).

Run:

    PYTHONPATH=.:src python3 scripts/certify_r12_h10_16_nikulin_obstruction.py --check-archive

Novelty is restricted to the carrier obstruction and resulting order-four
criterion. No new abstract-polytope, arithmetic-class, or worldwide
priority claim is made here.

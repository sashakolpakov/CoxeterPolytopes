# A compact all-dimensional counterfamily to field-only descent

20 October 2026. Own elementary construction, then primary-source check;
no scientific execution, search worker, Lean or TeX. This is a barrier to
a proposed implication, not a new dimension bound or a novelty claim.

For every integer n>=2 put K=Q(sqrt2), O=Z[sqrt2], and

    c=(2n+sqrt2)/(2n^2),
    B_ii=1,       B_ij=-c (i!=j),       size(B)=n+1.

This is the Gram matrix of an ACTUAL compact acute-angled regular simplex
in H^n. Indeed 1/n<c<1/(n-1), so B has eigenvalues 1+c (multiplicity n)
and 1-nc<0. Every proper principal submatrix is positive definite. More
explicitly its inverse has negative off-diagonal entries and diagonal

    (1-(n-1)c)/[(1+c)(1-nc)]<0.

Thus the dual vertex rays are timelike in one common time component and
give the compact simplex; this is not an abstract corner or an unfinished
wall prefix. Its dihedral angles are theta=arccos(c), strictly acute.

At the other real embedding, 0<c'=(2n-sqrt2)/(2n^2)<1/n, so B' is
positive definite. Its normalized cyclic-Gram field is exactly K,
because a two-cycle already gives c^2 with a nonzero sqrt2 coefficient.
Even the whole form can be integral: A=2n^2 B has diagonal 2n^2 and
off-diagonal -(2n+sqrt2), all in O, with the same admissible signatures.
Hence compactness, an integral admissible form, a fixed quadratic field,
conjugate positivity and preservation of that cyclic field together
do not give any dimension bound for acute-angled polytopes.

The missing condition is reflection integrality. For i!=j the coefficient
of e_i in R_i(e_j)-e_j is

    -2A_ij/A_ii=(2n+sqrt2)/n^2,

which is not in O. More strongly, these facet reflections generate a
NONDISCRETE group. The rational trace of 4c^2 from K to Q is

    Tr(4c^2)=(8n^2+4)/n^4.

It equals9/4 at n=2 and lies strictly between0 and1 for every n>=3,
so 4c^2 is not an algebraic integer. If theta/pi were rational, then
4cos^2(theta)=2+zeta+zeta^-1 for a root of unity zeta, and would be an
algebraic integer. Thus theta/pi is irrational. A pair of adjacent wall
reflections has an infinite-order elliptic product, whose powers
accumulate at the identity. These are not quasi-arithmetic reflection
LATTICES or Coxeter polytopes.

The distinction is already built into the primary literature. Bogachev
and Kolpakov extend quasi-arithmeticity to finite-volume acute-angled
polytopes by the field/conjugate conditions in Definition2.9, p5 of
[arXiv:2002.11445v4](https://arxiv.org/pdf/2002.11445). Their Lemma3.3
and Proposition3.4, pp6-7, prove the same-ground-field property for
facets without a Coxeter hypothesis. Section3.3, pp7-8, iterates through
possibly non-Coxeter intermediate faces. Those results supply field
descent; they do not assert that the wall-reflection groups of such
intermediate faces are discrete or integral. The family above satisfies
their extended field/conjugate conditions in every dimension.

There is an important remaining constraint for the intended application.
An actual facet of a compact Coxeter polytope has the finite rank-three
projected-angle alphabet in
`coxeter_metric_facet_angle_rigidity_20261020.md`: every nonright cosine
is at least1/3. For n>=4 the present family has 0<c<1/3, so it does NOT
model such a facet. A stronger argument using that alphabet together
with the AMBIENT reflection-integrality conditions is still possible;
the counterfamily blocks only the field/conjugate/integral-form shortcut.
No universal compact or arithmetic ceiling is changed.

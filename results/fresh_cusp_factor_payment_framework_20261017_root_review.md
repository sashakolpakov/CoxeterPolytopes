# Root review of cusp-factor payments and spectator elimination

The restricted theorem and the Schur-complement lemma pass this internal
mathematical/source review. The result is dimension at most117 when every
irreducible affine cusp factor has rank at least9 and every cusp has at most
two factors. Neither hypothesis is dropped. The unrestricted163 theorem is
unchanged; external refereeing and formal verification have not occurred.

This review binds the final [proof report](fresh_cusp_factor_payment_framework_20261017.md),
[source](../scripts/fresh_cusp_factor_payment_framework_20261017.py) and
[certificate](fresh_cusp_factor_payment_framework_20261017.json):

* source SHA256: `93f554e77649cfe1c1a76696824591e59b1028756d48cd1fa9d9fb4b23969100`;
* certificate SHA256: `c55a84924abce529ce88697b611f991a2b9ae7c5c9e1e20f6948abf88239bcd9`;
* report SHA256: `c65cfc5a4409a8cc83fcf0a083a6c1a7a29bb7bc7d649c11435773877e70ccab`.

## Scope and corner multiplicity

The reviewer checked the argument from actual finite-volume incidence, not
from arbitrary positive-semidefinite Gram matrices. A two-face has n−2
supporting roots, and adding one side gives a positive-definite edge support.
At a nonsimple ideal corner the link is a product of simplices. A link edge
varies in one factor; its two omitted nodes are the marked pair in that
factor, while one node is omitted from every other factor. Thus a pair in
factor i has precisely the stated product of other factor vertex counts,
V/(r_i+1), as its multiplicity. The corner budget7sV follows from this
multiplicity and the all-rank graph bound. Counting only one copy of a pair
at a nonsimple cusp would be incorrect; that shortcut is not used here.

The positive null relation for an affine subdiagram is nonzero because it
pairs strictly negatively with a chamber center. It gives a feasible ideal
point. Containment in the complete cusp diagram then implies that a
connected affine subdiagram is an entire factor, since proper principal
subdiagrams of an irreducible affine diagram are positive definite.
This agrees with the parabolic-containment statement in
[Felikson–Tumarkin, Section1](https://arxiv.org/pdf/math/0502413).
Consequently the minimum-rank hypothesis excludes all affine diagrams of
order at most9, and the inherited Lanner/quasi-Lanner classification reduces
minimal indefinite diagrams to order at most5.

The triangle and quadrilateral short-path arguments use actual closed-side
incidence. Positive-definite opposite-side support would give a finite
intersection; singular positive-semidefinite support would give an ideal
intersection. Both are impossible for distinct opposite closed sides of
an actual quadrilateral. Two orthogonal indefinite cores would contradict
the ambient negative index1. Their overlap or joining edge therefore gives
the asserted paid adjacent pair. The finite-vertex7n and ideal7sV bounds
use the same distance-seven rule.

## Strict original-polytope inequality

The reviewer checked [Khovanskii, Theorems8–10 and Section6](https://www.math.utoronto.ca/askold/1986-Faa-1-english.pdf)
directly. The application is to the bounded projective closure, which is
simple at edges; it does not replace the original two-faces by truncation
faces. For a generic linear function, let E_w count vertex-figure edges
and E_w^0 those not cut by the associated hyperplane. The exact reduction
identities give I=sum E_w and2f2=sum E_w^0. The local ratios2E_w/E_w^0
are at most eta_n. At the global minimum vertex the hyperplane misses
the figure entirely, so this ratio is2, strictly below eta_n. All weights
E_w^0 are positive in the dimensions used. Hence I/f2<eta_n strictly.

With I=binom(n,2)v_f+(n−1)sum V_v/2, this gives the strict D5 lower bound
in the report. Combining it with the corner budget yields a strict
factor-weighted inequality. At n=56s+6 the lower side already equals s,
and the two displayed parity polynomials exclude every higher dimension.
Thus56s+5 is valid, including117 for s=2 without requiring a finite vertex.

## Spectator lemma and limitations

The upper-left affine Schur block has positive alpha and rank1. Subtracting
a positive-semidefinite spectator Gram correction only in the lower2x2
block changes the determinant by
alpha(pk²−2qhk+rh²). The terms quadratic in the correction vanish because
the affine block has determinant0. Division by h²k² gives the displayed
nonnegative quadratic form. Here h and k must be nonzero. Four residual
normals of an actual quadrilateral lie in a Lorentz three-space after all
supporting whites are eliminated, so their determinant is zero. A strictly
positive determinant obstruction therefore survives these spectators.

This proves the scalar transfer, not every step of the finer weight-three
payment argument. Changes of finite component type and newly permitted
small affine factors need their own classification. The report explicitly
retains these obligations. The rational two-factor local cusp model is a
neighborhood example, not a global finite-volume chamber in every dimension.

The review checks the written all-rank arguments, source hypotheses and
certificate interpretation. Producer generation and saved replay passed;
this review does not claim an independent rerun of every inherited graph
or a proof-assistant formalization. No mathematical gap was identified in
the stated restricted theorem or determinant lemma.

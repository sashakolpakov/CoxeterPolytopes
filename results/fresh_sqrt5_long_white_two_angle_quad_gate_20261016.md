# A long-white quadrilateral pattern is arithmetically impossible

A specific quadrilateral with exactly two discounted adjacent corners cannot occur over Q(√5). This excludes one possible source of an adjacent-discount cap of 1/2; it does not cover other quadrilaterals or certify any weight rule.

Let the white diagram be A_r, r≥4, attached by one 3-edge from an endpoint only to black root 0. Take black labels 01=03=3 and 12=23=2. Leave the two opposite Gram products −p at 02 and −q at 13 free. All roots have norm 2. Other common white components are assumed orthogonal to this entire local graph. A quadrilateral face requires this local Gram matrix to have Lorentz index one and corank at least one. We impose integral cyclic products in Z[φ], φ=(1+√5)/2, and positive semidefiniteness at the other real embedding.

Put m=r+1≥5. Eliminating A_r gives diagonal a=1+1/m at root 0 and leaves diagonal 2 elsewhere. The exact residual determinant is

\[
(2+q)\bigl((2a-p^2)(2-q)-4\bigr).
\]

Thus (a−p²/2)(2−q)=2. Acute products give q≥0. If q<2, eliminating root 2 and then the positive-definite block on roots 1,3 leaves a zero Schur complement, making the full matrix positive semidefinite. This contradicts the required Lorentz index. The value q=2 fails the determinant equation. Consequently

\[
t=q-2=\frac{4}{p^2-2-2/m}>0.
\]

The triangle 0–1–3–0 has cyclic product −q, and the two-cycle 0–2–0 has product p². Hence q,p²∈Z[φ]. Conjugate semidefiniteness of the 1,3 block gives −4≤t′<0. The strict upper inequality follows because a nonzero algebraic number cannot have zero conjugate. Multiplying the displayed relation shows

\[
mt(p^2-2)=2t+4m,\qquad m\mid 2t.
\]

For p≥2, including the parallel value, 0<|N(t)|≤8m/(m−1)≤10. Odd m≥5 forces m² to divide N(t), which is impossible. Even m≥8 forces (m/2)² to divide N(t), also impossible. At m=6 write t=3s. The norm bound forces N(s)=−1, whereas 3s(p²−2)=s+4 gives s≡−4 modulo 3 and therefore N(s)≡1 modulo 3. This is again impossible.

For finite p, the condition p²∈Q(√5) restricts the Coxeter label to 2,3,4,5,6,10. Indeed the real cyclotomic degree forces Euler's totient to be at most four; the additional possibilities 8 and 12 generate Q(√2) and Q(√3). Labels 2,3,4 have p²≤2 and contradict t>0. Label 6 has rational p²=3 and hence rational q>2, contradicting conjugate semidefiniteness.

For label 5,

\[
q=2+\frac{8m+4m^2\phi}{m^2-2m-4}.
\]

The six cases m=5,…,10 are explicitly nonintegral. For every m≥11, the rational coefficient has increment strictly between zero and one, since

\[
(m^2-2m-4)-8m=(m-11)(m+1)+7>0.
\]

For label 10, conjugate semidefiniteness requires φ≥2−2/m. This fails for m≥6 because φ<5/3. At m=5 the forced value q=(2+100φ)/31 is nonintegral. These cases exhaust every acute opposite product p.

The [standard-library source](../scripts/fresh_sqrt5_long_white_two_angle_quad_gate_20261016.py) checks the determinant polynomial, exact finite-label formulas, divisibility residues, and the endpoint inequalities used in the all-rank argument. The [certificate](fresh_sqrt5_long_white_two_angle_quad_gate_20261016.json) passed generation and a fresh saved replay. A separate agent read the normalization, integrality, conjugate signs, and exceptional m=6 argument; this was a narrow formula/scope check, not an independent census or complete certificate replay.

**Scope:** one explicit white-attachment pattern, with all opposite products free and all r≥4 covered. No dimension improvement, complete short-face classification, or general permission to discount other corners follows.

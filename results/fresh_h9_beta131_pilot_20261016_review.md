# Independent check of the determinant-norm-131 H9 pilot

The [independent checker](../scripts/fresh_h9_beta131_pilot_20261016_review.py)
and [certificate](fresh_h9_beta131_pilot_20261016_review.json) reconstruct the
complete height-20 Vinberg prefix on

\[
L=E_8\perp\begin{pmatrix}-2\phi&-1\\-1&6\end{pmatrix},
\qquad \phi=(1+\sqrt5)/2,
\]

over \(\mathbb Z[\phi]\). There are **52 height controls and 12 actual walls**:
eight zero-height E8 walls and four positive-height walls. The independent
recursion visits 1,932 nodes. This is a bounded search result, not a
finite-volume example or a proof of nonreflectivity.

The implementation imports project `Golden` arithmetic, but no producer
or producer arithmetic helper. It rebuilds the literal E8 matrix, uses
natural rather than reversed covector order, independently enlarged
coefficient boxes, and a weaker coordinate alphabet than the producer.
It recomputes every control and every accepted root, all previous-wall
inequalities, all facet witnesses and all Coxeter pair relations. Separate
rational matrix inversion supplies its prefix energy bounds. It does not
independently enumerate the cone's vertices or exterior rays.

## Norm and zero-height completeness

Put \(\delta=1+12\phi\) and \(\beta=\phi\delta=12+13\phi\).
The lattice determinant is \(-\delta\), with norm \(-131\); \(\beta\)
is totally positive with prime norm 131. The physical form has signature
(9,1), whereas its conjugate is positive definite. In coordinates
\((x,y_1,\ldots,y_8,t)\), the basepoint is \(H=(6,0^8,1)\), and
\((H,H)=-6\delta\), \((H,r)=-\delta x\).

The following standard ring facts are explicit mathematical inputs:
\(\mathbb Z[\phi]\) is a PID with unit group \(\{\pm\phi^j\}\).
The lattice is even. For a primitive reflector with norm \(2b\),
reflection integrality gives \(Gr/b\in L\). Multiplication by the
adjugate and a primitive-coordinate Bezout identity imply
\(b\mid\det G\). Total positivity and the fact that every totally
positive unit is a square leave exactly the normalized norm classes
**2 and \(2\beta\)**. The numerical pilot includes both.

At \(x=0\), a norm-\(2\beta\) reflector would have all E8 coordinates
and \(t\) divisible by \(\beta\). Division by \(\beta\) would then
give an integral lattice vector of norm \(2/\beta\), contradicting
evenness. For norm 2, the two embeddings give \(6t^2\le2\) and
\(6(t')^2\le2\), forcing the integral norm of \(t\) to vanish, hence
\(t=0\). Writing the E8 vector as \(u+\phi v\), the norm equation
gives \(u^2+v^2=2\) and \(2(u,v)+v^2=0\). E8 evenness and positive
definiteness force \(v=0\). Thus the zero-height arrangement is exactly
ordinary E8. The checker verifies 240 distinct integral norm-two vectors
against its independently rebuilt E8 matrix; completeness uses the
standard E8 root classification and this argument.

## Exhaustive bounded enumeration

Let \(C\) be the E8 Cartan matrix, \(J=C^{-1}\), \(\lambda=Cy\),
and \(l=x-6t\). In the chosen E8 sector, \(\lambda_i\ge0\), while
\(l\) has either sign. The exact equation is

\[
l^2+6\lambda^T J\lambda=6q+\delta x^2,
\qquad l\equiv x\pmod{6\mathbb Z[\phi]}.
\]

Both \(\lambda\) and \(l\) are divisible by \(q/2\). The checker
also directly verifies full reflection integrality at every accepted root.
Height is \(\delta x^2/(6q)\). Using \(q<67\), \(q'<8\),
\(\delta>20\), and \(-\delta'>6\), the finite coefficient boxes
follow from \(x^2<402T/20\), \((x')^2<8\),
\(l^2<402(1+T)\), and \((l')^2<48\).

The same bounds contain all covector coordinates: minimization of the
positive quadratic form gives \(6\lambda^TJ\lambda\ge3\lambda_i^2\)
at either embedding. The independent recursion deliberately uses the
larger factor-two coordinate alphabet and checks the producer's
factor-three alphabet separately. Prior-wall pruning starts only after
the signed \(l\) is fixed; all remaining physical coefficients are
nonnegative. No sign condition is imposed on conjugate covectors.

The determinant norm separates this form from the prior golden
determinant-norm-11 and -55 forms even up to similarity: rank 10 is even,
so scalar similarity changes the determinant only by a square, and 131
has odd valuation in the norm ratio. This does not classify other hosts
or establish literature novelty.

Generation and saved replay use:

```sh
python3 scripts/fresh_h9_beta131_pilot_20261016_review.py
python3 scripts/fresh_h9_beta131_pilot_20261016_review.py --verify
```

The proof of the norm alphabet and the bounded enumeration do not imply
that the full reflection chamber has finite volume. Higher roots and
global chamber geometry remain separate questions.

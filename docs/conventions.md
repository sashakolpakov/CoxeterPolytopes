# Mathematical and Data Conventions

These conventions are normative for code and certificates. Any imported source
using different signs or normalizations must be converted explicitly.

## Quadratic spaces

- A Lorentzian space for \(\mathbb H^d\) has signature \((d,1)\): \(d\)
  positive and one negative direction.
- Hyperbolic space is the projectivization of one component of the negative
  cone.
- A wall normal/root \(r\) is spacelike: \((r,r)>0\).
- Reflection in \(r^\perp\) is
  \[
  s_r(x)=x-\frac{2(x,r)}{(r,r)}r.
  \]
- For a lattice \(L\) over a ring \(\mathcal O\), an integral root satisfies
  \[
  \frac{2(r,L)}{(r,r)}\subseteq\mathcal O.
  \]

## Gram matrices and Coxeter relations

The wall/root Gram matrix is unnormalized:

\[
B_{ij}=(r_i,r_j).
\]

Simple roots are oriented so that \(B_{ij}\leq0\) for \(i\ne j\). For positive
norms \(a=B_{ii}\), \(c=B_{jj}\), and \(b=B_{ij}\), define

\[
q(i,j)=\frac{4b^2}{ac}.
\]

If the walls meet at angle \(\pi/m\), then

\[
q(i,j)=4\cos^2(\pi/m).
\]

In particular:

| \(m\) | \(q\) |
|---:|---:|
| 2 | \(0\) |
| 3 | \(1\) |
| 4 | \(2\) |
| 5 | \(\phi^2=\phi+1\) |
| 6 | \(3\) |
| 10 | \(\phi+2\) |

The walls are parallel when \(q=4\), and ultraparallel when \(q>4\).
Values \(0<q<4\) not recognized exactly as \(4\cos^2(\pi/m)\) are recorded as
unrecognized, not silently rounded to a label.

The normalized Gram matrix used in combinatorial Coxeter-polytope work has
diagonal 1 and intersecting entry \(-\cos(\pi/m)\). It need not lie in the
same base field as the unnormalized root Gram matrix.

## Golden field

Write

\[
\phi=\frac{1+\sqrt5}{2},\qquad \phi^2=\phi+1.
\]

Elements of \(\mathbb Q(\sqrt5)\) are stored uniquely as \(a+b\phi\), with
\(a,b\in\mathbb Q\). The two real embeddings are

\[
\sigma_0(\phi)=\frac{1+\sqrt5}{2},\qquad
\sigma_1(\phi)=\frac{1-\sqrt5}{2}.
\]

Serialized coefficients are reduced rational strings, for example
`{"a":"3/2","b":"-1"}`. Sign decisions at both embeddings are exact.

## Evidence states

Candidate computations use exactly these terminal labels:

- `compact_reflective`
- `cofinite_noncompact`
- `nonreflective_proved`
- `undecided`

Resource exhaustion, numerical solver failure, or nontermination always maps to
`undecided`.


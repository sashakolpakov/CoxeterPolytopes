# Independent exact audit of the rank-16 K3 wall fold

## Verdict

The proposed divisor-wall fold does produce a genuine hyperbolic Coxeter
polytope:

\[
P\subset \mathbb H^{14}
\]

has 17 irredundant facets, finite volume, five cusps, and is noncompact.
This conclusion passes the literal Vinberg--CoxIter endpoint and full
zero-neighbourhood tests; it is not inferred merely from a Lorentzian Gram
matrix or a list of candidate vertices.

It is nevertheless **not a new abstract example**: the edge-labelled diagram
is Vinberg's published arithmetic \(\mathbb H^{14}\) chamber, recorded by
CoxIter as `14-vinb85.coxiter`.  The fold is a useful algebro-geometric
construction of that chamber from a complete K3 divisor-wall configuration.

A separate statement must not be conflated with the Coxeter result: the small
graph involution used in the fold is not induced by an involution of a K3
surface having exactly the advertised Néron--Severi lattice.  Thus the cone is
not an actual K3 quotient-nef cone for this action.  Surface-automorphism
realizability is not needed for the abstract fixed-cone construction.

## Input and exact fixed cone

The source is the complete 19-wall nef configuration for

\[
N=U\oplus E_8(-1)\oplus D_6(-1),\qquad \operatorname{rk}N=16.
\]

The weighted wall Gram graph admits the involution

\[
(1)\cdots(15)(16\ 18)(17\ 19).
\]

Its fixed lattice is saturated of rank 15 and signature \((1,14)\).  Restricting
the complete nef inequalities gives exactly 17 orbit inequalities: singleton
walls remain unchanged, while the two equal inequalities in each two-cycle
are represented by the orbit-sum normal.  The 17 orbit normals have squares

\[
(-2)^{15},-4,-4,
\]

and their Gram matrix has rank 15 and inertia \((1,14,2)\).

An integral strict interior point, in the fixed basis, is

\[
(22,46,72,35,65,60,57,56,57,60,65,72,35,46,22).
\]

It has square 1540 and pairings

\[
(2^{15},44,44)
\]

with the orbit normals.  For every one of the 17 inequalities, the equality
extreme rays span a rank-14 hyperplane and their sum is timelike, lies on that
wall, and pairs strictly positively with every other wall.  Hence every wall
is irredundant and has strict relative interior.

The 17 hyperplane reflections were recomputed in the fixed integral basis.
Every matrix is integral, squares to the identity, and preserves the fixed
lattice Gram matrix.  The complete Coxeter relation census is

\[
119\times m=2,\quad 14\times m=3,\quad 2\times m=4,
\quad 1\times m=\infty,
\]

with no non-Coxeter or ultraparallel pair.

## Literal finite-volume test

The independent certifier enumerates all principal subdiagrams by exact
rational arithmetic.  It finds:

- 704 elliptic rank-13 edge diagrams;
- 89 elliptic rank-14 finite vertices;
- five componentwise-affine parabolic rank-13 ideal vertices;
- no hyperideal extreme ray.

Every elliptic edge has exactly two legal endpoints.  Their type distribution
is

| finite endpoints | ideal endpoints | edges |
|---:|---:|---:|
| 2 | 0 | 549 |
| 1 | 1 | 148 |
| 0 | 2 | 7 |

For an ideal incidence, the code requires the edge to be obtained by deleting
exactly one node from every affine component of the cusp diagram.  It also
performs Vinberg's other global check independently: all eight connected
affine subdiagrams have a full zero-neighbourhood which is parabolic of total
rank 13.  Those neighbourhoods give exactly the same five cusp diagrams.

Therefore the conclusion is rigorously

\[
\boxed{\text{finite volume, noncompact, with 89 finite vertices and 5 cusps}.}
\]

These numbers independently reproduce the published beginning of the
\(f\)-vector, \((94,704,\ldots)\), for Vinberg's chamber.  See
R. Guglielmetti,
[*CoxIter—Computing invariants of hyperbolic Coxeter groups*](https://doi.org/10.1112/S1461157015000273),
Appendix A, Figure A.1.

## Why this is not a quotient-nef cone

The first 16 published walls form a basis of determinant \(-4\), with Smith
diagonal

\[
1^{14},2,2;
\qquad A_N\cong(\mathbb Z/2)^2.
\]

The proposed small involution acts nontrivially on \(A_N\).  For example, if
\(y=B^{-1}e_4\in N^\vee\), then \(Jy-y\) has half-integral coordinates and so
does not lie in \(N\).  By contrast, the published nontrivial double-plane
deck action is the large reversal of the wall graph and acts trivially on
\(A_N\).

For completeness, suppose a K3 involution with full Néron--Severi lattice
\(N\) induced \(J\).  If its eigenvalue on the period is
\(\varepsilon\in\{\pm1\}\), the opposite rational eigenspace in the
transcendental lattice is orthogonal to the period and would supply an extra
integral \((1,1)\)-class.  Thus the transcendental action must be
\(\varepsilon I\).  Both \(I\) and \(-I\) act trivially on the 2-elementary
discriminant group, while unimodular K3 gluing forces the discriminant actions
on \(N\) and its orthogonal complement to agree.  The explicit witness gives
a contradiction.

This rejects only the surface-involution/quotient interpretation.  It does not
undo the exact Coxeter polytope constructed from the divisor-wall data.

The rank-16 surface and its automorphisms are described in X. Roulleau,
[*An atlas of K3 surfaces with finite automorphism group*](https://arxiv.org/abs/2003.08985),
Section 15.2.

## Reproduction

The independent certifier imports neither the production fold verifier nor
the production Vinberg-volume implementation:

```text
PYTHONPATH=src:. python3 scripts/independent_certify_k3_rho16_h14_fold.py --write-certificate
PYTHONPATH=src:. python3 scripts/independent_certify_k3_rho16_h14_fold.py --verify-certificate
PYTHONPATH=src:. pytest -q tests/test_k3_rho16_h14_fold_audit.py
```

The frozen artifacts are
`data/candidates/k3_rho16_h14_fold_audit.json` and
`tests/test_k3_rho16_h14_fold_audit.py`.

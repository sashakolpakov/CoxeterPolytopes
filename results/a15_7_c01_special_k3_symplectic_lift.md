# A genuine special-K3 lift of the 27-wall `H10` chamber

## Theorem

The fixed section `A15_7_1_C01` of Alexeev's complete 60-wall
((15,7,1)) reflection chamber is realized by a genuine automorphism:

> There is a projective K3 surface (X) of Picard rank 19 and a symplectic
> Nikulin involution (g\in\operatorname{Aut}(X)) such that the archived
> rank-15 lattice (N) embeds in (NS(X)), the restriction (g|_N) is exactly
> the archived `C01` action, and
> \[
> NS(X)^g=N^{C01}
> \]
> integrally.  The fixed section therefore remains the finite-volume,
> noncompact 27-facet Coxeter polytope in (\mathbb H^{10}).

The Picard lattice grows from rank 15 to rank 19, but all four added
algebraic directions are anti-invariant.  Hence the fixed divisor rank stays
11 and the hyperbolic dimension stays 10.

This is a positive algebro-geometric realization result, but it is **not a new abstract
hyperbolic polytope**.  Its signed normalized Gram matrix is
exactly that of the previously frozen `D16-D4-C7-H10-27` chamber.  The
dimension, arithmetic lattice, maximal reflection group, and symplectic
involution type were also known.  What is new in the project is the exact
special-K3 automorphism lift of the Alexeev `C01` carrier.

## The complement and its involution

Write (N) for the even 2-elementary lattice with invariants ((15,7,1)).
Its discriminant group has order (2^7).  In the archived rank-15 basis, a
basis of half-vector representatives is supported at the zero-based
coordinates

```text
3, 4, 7, 8, 11, 12, 13.
```

Take the orthogonal complement stock

\[
T=U(2)\oplus U(2)\oplus\langle-2\rangle^3=2T_0,
\]

where

\[
T_0=U\oplus U\oplus\langle-1\rangle^3.
\]

Inside (T_0), use the four rows

\[
K=\begin{pmatrix}
0&0&0&0&1&1&0\\
-2&-2&-2&1&-2&1&-1\\
-2&-1&-2&0&-2&1&1\\
-2&0&0&2&0&-1&-1
\end{pmatrix}.
\]

Their Gram matrix is

\[
K T_0 K^t=
\begin{pmatrix}
-2&1&1&1\\
1&-2&0&0\\
1&0&-2&0\\
1&0&0&-2
\end{pmatrix},
\]

and their Smith diagonal in (\mathbb Z^7) is (1^4), so they span a
primitive (D_4(-1)) in (T_0), equivalently a primitive (D_4(-2)) in (T).
Define

\[
B=I-2T_0K^t(KT_0K^t)^{-1}K.
\]

The exact replay proves that (B) is integral, (B^2=1), and
(B G_T B^t=G_T).  Its characteristic polynomial is

\[
(x-1)^3(x+1)^4.
\]

The fixed lattice has signature ((2,1)), while the anti-invariant lattice is
the displayed (D_4(-2)).

## The seven-bit gluing gate

Because ((15,7,1)) has (\delta=1), reducing the discriminant quadratic form
only modulo 2 would lose essential information.  The certificate instead
records

\[
q_2(x)=2q(x)\pmod 4
\]

for all 128 discriminant classes.  The source histogram for
(q_2=0,1,2,3) is ((28,36,36,28)).

In the displayed discriminant bases, the anti-isometry is

\[
P=\begin{pmatrix}
1&0&0&1&1&0&1\\
0&0&0&1&0&1&1\\
0&1&0&1&0&1&1\\
1&0&1&0&1&1&0\\
0&1&0&0&0&1&1\\
0&1&1&0&1&0&1\\
0&0&0&0&1&1&1
\end{pmatrix}.
\]

Its integral determinant is 3, hence it is invertible over
(\mathbb F_2).  If (M_N) is the archived action on (A_N), then the exact
checks are

\[
M_NP=P(B\bmod2),\qquad q_{2,T}(xP)=-q_{2,N}(x)\pmod4
\]

for every (x\in\mathbb F_2^7).  Thus the graph gluing is equivariant and is
an anti-isometry of discriminant forms, not merely an isomorphism of the
underlying groups.

Gluing (N\oplus T) along this graph gives an even unimodular lattice of
signature ((3,19)), hence the K3 lattice.  The glued integral involution has

\[
\chi_g(x)=(x-1)^{14}(x+1)^8.
\]

Its saturated anti-invariant lattice has determinant (2^8).  Dividing its
Gram matrix by two gives an even, negative-definite, unimodular rank-eight
lattice.  Consequently

\[
H^2(X,\mathbb Z)^-=E_8(-2).
\]

In particular, it contains no vector of square (-2).

## An explicit generic fixed period

A convenient saturated basis of (T^B) has Gram matrix

\[
H=\begin{pmatrix}
-6&2&14\\
2&2&-6\\
14&-6&-30
\end{pmatrix}.
\]

Set

\[
u=(0,1,0),\qquad a=(-1,1,0),\qquad b=(3,0,1),\qquad p=2^{1/4},
\]

and

\[
v=p a+\left(p+\frac1{4p}\right)b,
\qquad \omega=u+iv.
\]

An exact calculation gives

\[
u^2=v^2=2,\qquad u\cdot v=0,
\]

so

\[
\omega^2=0,\qquad \omega\bar\omega=4,
\qquad B\omega=\omega.
\]

The irrationality is deliberate.  Since

\[
v=p(a+b)+\frac{p^3}{8}b,
\]

a rational fixed vector perpendicular to (\omega) must pair trivially with
(u), (a+b), and (b) separately.  The corresponding rational constraint
matrix has determinant (-32).  Thus no nonzero rational vector in (T^B) is
orthogonal to the period.  The rational orthogonal in all of (T) is exactly
the rank-four anti-invariant (D_4(-2)).

The saturated Picard lattice recovered inside the graph glue therefore has

\[
\rho(X)=19,\qquad \operatorname{sig}NS(X)=(1,18),\qquad
\det NS(X)=32,
\]

with Smith diagonal

\[
(1^{16},2,4,4).
\]

An integral unimodular row-lattice comparison, rather than a rational-rank
comparison, proves

\[
NS(X)^g=N^{C01}
\]

exactly.  The latter has rank 11, determinant 512, and Smith diagonal
((1^4,2^5,4^2)).

## Kähler chamber and Torelli

The action permutes all 60 walls of the source reflection chamber.  The
archived fixed interior vector has square 7,360 and minimum pairing 16 with
those walls, so the invariant positive section is nonempty and strict.

The only remaining obstruction would be a K3 root hyperplane containing the
entire invariant positive cone.  Such a root would be anti-invariant, but
the anti-invariant lattice is (E_8(-2)), whose norms are divisible by four.
It contains no `(-2)` vector.  Hence every root cuts the invariant cone in a
proper hyperplane, and a generic invariant point in the strict source
section avoids all root hyperplanes.  The Weyl chamber containing that point
is invariant.

After choosing the marking so that chamber is Kähler, the glued involution is
an effective Hodge isometry fixing both the period and a Kähler class.
Strong global Torelli realizes it by a genuine K3 automorphism.  Its action
on (H^{2,0}) is the identity, and its anti-lattice is (E_8(-2)), so it is a
symplectic Nikulin involution.  For the classical lattice theory of
prime-order symplectic K3 automorphisms, see Garbagnati--Sarti,
[*Symplectic automorphisms of prime order on K3
surfaces*](https://arxiv.org/abs/math/0603742).

## Literal Vinberg finite-volume handoff

The automorphism lift does not substitute for a volume proof.  The separate
literal Vinberg certificate gives:

| datum | exact value |
|---|---:|
| dimension | 10 |
| facets | 27 |
| elliptic edges | 6,221 |
| finite vertices | 895 |
| ideal vertices/cusps | 67 |
| bad edge endpoints | 0 |
| bad affine zero-neighbourhoods | 0 |

Every elliptic edge has exactly two finite-or-ideal endpoints, and every
irreducible affine component has the required full zero-neighbourhood.
Therefore the polytope has finite volume.  Its 67 cusps prove that it is
noncompact.

## Lineage and exact scope

The lineage certificate supplies an explicit 27-facet permutation proving
signed normalized-Gram equality with `D16-D4-C7-H10-27`.  Thus this result
does not add a new abstract polytope, dimension, stable lattice, maximal
reflection group, or arithmetic commensurability class.  The 27-wall chamber
already has exact index 2,016 in the known 14-wall full reflection chamber of
the determinant-32 ((11,5,1)) over-lattice.

What this theorem adds is more rigidly algebro-geometric: the `C01` symmetry
is no longer only an isometry of a K3-admissible Picard lattice and its wall
stock.  It is the restriction of an actual symplectic automorphism of an
explicitly certified special K3 period, while the fixed divisor lattice and
the 27-wall hyperbolic section remain unchanged.

The source 60-wall object is a mixed-root reflection chamber.  It is **not
asserted to be the literal nef cone** of the special rank-19 surface, and the
four extra algebraic directions may add walls to the actual nef cone.  No
quotient-nef statement is made.

## Reproduction

The permanent files are:

- [exact K3-lift certificate](../data/enumerations/a15_7_c01_special_k3_symplectic_lift.json);
- [certificate generator](../scripts/certify_a15_7_c01_special_k3_symplectic_lift.py);
- [regression tests](../tests/test_a15_7_c01_special_k3_symplectic_lift.py);
- [source construction certificate](../data/enumerations/alexeev_remaining_high_fixed_section_construction.json);
- [literal Vinberg certificate](../data/enumerations/alexeev_remaining_high_fixed_section_vinberg.json); and
- [exact duplicate/lineage certificate](../data/enumerations/alexeev_remaining_high_fixed_section_lineage.json).

Run from the repository root:

```bash
PYTHONPATH=src:. python3 scripts/certify_a15_7_c01_special_k3_symplectic_lift.py
PYTHONPATH=src:. pytest -q tests/test_a15_7_c01_special_k3_symplectic_lift.py
```

The certificate payload SHA-256 is
`98304fee498980acae78915903155f5310778dc29cb718d842ac8e41cb2213e9`.
The frozen JSON file SHA-256 is
`dc5fad771088ac750d7e98640534b38936dc0ff46423b958b7ab6a58e06f3a87`.

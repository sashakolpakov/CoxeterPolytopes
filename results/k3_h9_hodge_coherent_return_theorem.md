# Hodge-coherent return at the quadratic K3 boundary

**Status (2026-08-09).**  This note separates two genuinely different
polytope-to-surface notions at the maximal real-quadratic K3 boundary.

- The existing **NS-only eigencarrier** asks only for the coefficient-field
  action on the divisor lattice.  Its `H9` locus is not exhausted.
- The stronger **Hodge-coherent eigencarrier** asks that the same field act
  unitally and self-adjointly on the full rational second cohomology as Hodge
  endomorphisms.  In this class the `H9` K3 branch is empty, every
  `K3^[k]` branch is empty, and the returned K3 bound is sharply `n<=8`.

The second assertion is a genuine

```text
polytope arithmetic -> complement/Hodge gate -> surface theorem -> bound on P
```

statement.  It does **not** turn the same bound into a theorem for the weaker
NS-only carrier class, much less for every arithmetic Coxeter polytope.

## 1. A surface-free coherent carrier package

Let `K` be a totally real field of degree `d`, and let `P` be a marked
arithmetic Coxeter realization whose rational full trace space `N_P,Q` is a
`K`-module of rank `n+1` and whose integral lattice is stable under the
specified order.  Thus

\[
 \operatorname{rk}_{\mathbf Q}N_P=d(n+1),\qquad
 \operatorname{sig}N_P=(1,d(n+1)-1).
\]

A **rational Hodge-coherent K3 package** for `P` is finite arithmetic data
constructed from `P`:

1. an oriented primitive embedding `N_P -> Lambda_K3` and its exact
   discriminant anti-isometry/glue;
2. the orthogonal complement `T`;
3. a faithful self-adjoint `K`-action on `T_Q` which, together with the given
   action on `N_P,Q`, gives a unital `K`-action on
   `Lambda_K3,Q=N_P,Q+T_Q`;
4. an oriented positive period plane in one real eigenspace of `T_R`, with
   no nonzero rational vector perpendicular to it; and
5. the already independent complete restricted-root, exact-Gram,
   reflection-integrality, and Vinberg predicates returning the same `P`.

At an optional stronger integral level, an order in `K` must preserve the
unimodular glue.  At an optional algebraic level, a rational codimension-two
correspondence must induce the field generator on `H^2`.

These conditions do not assume a surface.  The K3 period theorem constructs
one from items 1--4, with exact Picard lattice `N_P`; the K3 Weyl theorem and
item 5 return the same polytope.  Conversely, a Hodge-coherent surface gives
exactly this complement package.  Therefore the complement/Hodge test is an
exact `P -> surface` gate.  At `H9` it terminates before any root enumeration,
by the theorem below.

## 2. The period-multiplicity lemma

### Lemma 2.1

Let `V` be a rational Hodge structure of K3 type and let a totally real
field `K` of degree `d` act faithfully and unitaly on `V` by Hodge
endomorphisms.  Put `m=dim_K V`.  Then

\[
                         \boxed{m\geq2.}                 \tag{2.1}
\]

If the action is self-adjoint for the polarization, the real eigenspaces are
orthogonal, and the eigenspace selected by the action on `V^{2,0}` contains
the entire oriented positive period plane.

#### Proof

For `a in K`, its action on the one-dimensional complex line `V^{2,0}` is
by one embedding `sigma_0(a)`.  Since `K` is totally real, this scalar is
real.  Complex conjugation therefore gives the **same** scalar on
`V^{0,2}`.  Hence the real plane

\[
 (V^{2,0}+V^{0,2})\cap V_{\mathbf R}
\]

lies in the `sigma_0`-eigenspace.  Every real eigenspace of a `K`-module of
rank `m` has real dimension `m`, so `m>=2`.  Self-adjointness gives the
orthogonality assertion.  `square`

This elementary multiplicity argument is the part of the general
Zarhin--van Geemen theory needed here.  It is stronger than a bare Hodge-rank
count because it uses the two-real-dimensional period plane.

## 3. The returned K3 and `K3^[k]` theorems

### Theorem 3.1 -- Hodge-coherent K3 return bound

Suppose a full-trace degree-`d` arithmetic `H^n` realization has a rational
Hodge-coherent K3 package.  Then

\[
 d(n+1)+2d\leq22,
 \qquad
 \boxed{n\leq\left\lfloor\frac{22}{d}\right\rfloor-3}. \tag{3.1}
\]

In particular, for a real-quadratic coefficient field,

\[
                              \boxed{n\leq8}.            \tag{3.2}
\]

#### Proof

The divisor trace module uses `d(n+1)` rational cohomology dimensions.  Its
orthogonal complement contains the transcendental Hodge structure and is
stable because the field action is self-adjoint.  Lemma 2.1 requires at
least `2d` further dimensions.  Since a K3 surface has `b2=22`, (3.1)
follows.  Extra algebraic classes only make the transcendental complement
smaller and cannot weaken the inequality.  `square`

At the quadratic `H9` boundary the contradiction is especially transparent:
the trace lattice has rank twenty and the transcendental lattice rank two.
It would be one-dimensional over `K`, but one real eigenspace cannot contain
the period plane.

### Theorem 3.2 -- full-`H^2` `K3^[k]`, `k>=2`, obstruction

No manifold of `K3^[k]` type with `k>=2` admits a unital action of a
real-quadratic field on all of `H^2(X,Q)`.

Indeed `b2(K3^[k])=23`, whereas the rational dimension of every vector space
over a quadratic field is even.  Consequently the entire `K3^[k]` column
disappears from the **full-`H^2` coherent** compiler.  This says nothing
about an action confined to the divisor lattice or to a proper primitive
cohomology summand.

### Transfer back to polytopes

Let `C` be a surface-free class of marked real-quadratic Coxeter
realizations.  If an exhaustive compiler proves that every `P in C` has a
full-`H^2` Hodge-coherent carrier on either a K3 or `K3^[k]` (`k>=2`), and the complete
wall test returns the same `P`, then Theorems 3.1--3.2 give

\[
                              \dim P\leq8.               \tag{3.3}
\]

This is precisely the requested `P -> variety -> theorem on P` implication.
The unproved premise is essential surjectivity of this stronger compiler.
An NS-only carrier does not satisfy it automatically.

## 4. TOOT: Bugaenko makes the coherent bound sharp

The `H8` bound is not a vacuous strengthening.  Start with the already
certified Bugaenko K3 package

\[
 N=\langle2\rangle+\langle-2\rangle+E_8(-1)^2,
 \qquad
 T=N^\perp=\langle-2\rangle+\langle2\rangle+U.
\]

In the displayed basis of `T`, put

\[
G_T=
\begin{pmatrix}
-2&0&0&0\\0&2&0&0\\0&0&0&1\\0&0&1&0
\end{pmatrix},\qquad
B=
\begin{pmatrix}
-3&2&0&0\\-2&3&0&0\\0&0&0&5\\0&0&1&0
\end{pmatrix}.                                         \tag{4.1}
\]

Exact multiplication gives

\[
B^2=5I,\qquad B^tG_T=G_TB,
\]

and

\[
G_TB=
\begin{pmatrix}
6&-4&0&0\\-4&6&0&0\\0&0&1&0\\0&0&0&5
\end{pmatrix}>0,                                       \tag{4.2}
\]

with leading minors `6,20,20,100`.  Therefore the `+sqrt(5)` eigenspace of
`B` is positive definite of dimension two and the `-sqrt(5)` eigenspace is
negative definite.  Put `tau_T=(I+B)/2`; then

\[
\tau_T^2=\tau_T+I.
\]

The period can be made completely explicit as a CM joint eigenline.  The
integral matrix

\[
J=
\begin{pmatrix}
0&0&1&-5\\0&0&-1&-5\\-10&10&0&0\\2&2&0&0
\end{pmatrix}                                          \tag{4.3}
\]

satisfies

\[
J^2=-20I,\qquad JB=BJ,\qquad J^tG_T=-G_TJ.              \tag{4.4}
\]

For a completely executable period, write `s=sqrt(5)` and put

\[
u=(2,3+s,0,0),\qquad
v=(0,0,-5-s,-1-s),\qquad \omega=u+iv.                 \tag{4.5}
\]

Exact arithmetic in `Q(s)` gives

\[
 Bu=su,\quad Bv=sv,\quad Ju=-2sv,\quad Jv=2su,         \tag{4.6}
\]

so `omega` spans the joint complex eigenline

\[
             B=+\sqrt5,\qquad J=+i\sqrt{20}.            \tag{4.7}
\]

Moreover

\[
u^2=v^2=20+12s,\qquad u\mathbin{\cdot}v=0,
\]

and therefore

\[
\omega^2=0,\qquad
\omega\mathbin{\cdot}\bar\omega=40+24\sqrt5>0.        \tag{4.8}
\]

All identities in (4.6)--(4.8) are recomputed by the certificate, rather
than stored as declared eigenvalues.  No nonzero rational vector is
orthogonal to this line: its real orthogonal plane is the irrational
`-sqrt(5)` eigenspace of the rational operator `B`.  Period surjectivity
therefore produces a projective K3 with **exact** Picard lattice `N`.

The four matrices `I,B,J,BJ` are rationally independent, so they generate a
quartic copy of

\[
 E=\mathbf Q(\sqrt5,\sqrt{-20})=\mathbf Q(\sqrt5,i).     \tag{4.9}
\]

Because the period has exact Picard lattice, `T` is the transcendental K3
Hodge structure and is simple.  Zarhin's theorem makes `End_Hdg(T)` a number
field.  It contains the degree-four field `E`, while `dim_Q(T)=4`; hence

\[
                        End_{\rm Hdg}(T)=E.             \tag{4.10}
\]

Thus this is not an unjustified assertion that a random line in the positive
eigenspace is CM; the commuting CM generator is (4.3).  Nor does it conflict
with van Geemen's stronger small-multiplicity obstruction when the **entire**
endomorphism field is totally real: here `Q(sqrt(5))` is a proper totally
real subfield of the quartic CM field (4.10).

On `N`, use the original Bugaenko action `B_N=2tau_N-I`.  Both `A_N` and
`A_T` are `(Z/2)^2`.  Exact action on the two half-basis generators is the
identity on each discriminant group, so the actions agree across the
anti-isometry defining the unimodular K3 glue.  The executable check is
stronger still: it forms the explicit index-four basis matrix
`S=[N|T]` inside `U^3+E8(-1)^2`, verifies

\[
 S^tG_{K3}S=G_N\oplus G_T,
\]

and conjugates `B_N+ B_T` by `S`.  The resulting archived sparse `22 x 22`
matrix is integral, squares to `5I`, and is self-adjoint for `G_K3`.
Therefore the suborder `Z[sqrt(5)]` acts integrally and self-adjointly on
full `H^2(X,Z)`; this is a computed glue statement, not an assertion inferred
only from the group orders.

The maximal order does **not** act integrally on full cohomology in this
certificate.  On the `U` block, `(I+B)/2` has half-integral entries.  The
full `Q(sqrt(5))` action is rational; its restriction to `NS` is the original
integral `Z[phi]` action.

Finally, the algebraicity statement splits cleanly into its transcendental
and algebraic pieces.  Ramón-Marí proves the Hodge conjecture for arbitrary
self-products of a K3 whose transcendental Hodge endomorphism field is CM.
Applied to (4.10), the class of `B_T` in the `T tensor T` Künneth summand of
`H^4(X x X,Q)` is algebraic.  The prescribed operator `B_N` on `NS_Q` is a
rational linear combination of divisor-product correspondences.  Adding
these two cycles gives a rational codimension-two correspondence inducing the
full operator `B`, and therefore `(I+B)/2`, on `H^2`.

Buskin supplies a complementary modern route.  Put
`i_T=BJ/10`, `t=1+sqrt(5)`, and

\[
 z=(t-i_T)(t+i_T)^{-1}\in E.
\]

The Rosati involution fixes `sqrt(5)` and sends `i_T` to `-i_T`, so `z` and
`z^{-1}` are rational Hodge isometries.  Directly in `E`,

\[
 \sqrt5={29\over8}(z+z^{-1})-{15\over4}.               \tag{4.11}
\]

Extend `z,z^{-1}` by the identity on `NS_Q`; Buskin's theorem makes their
classes algebraic, and divisor products correct the resulting action on
`NS_Q` to `B_N`.  This independently gives the same rational
correspondence conclusion.  Neither route makes the whole quartic CM field
a unital action on integral `H^2`, and neither produces a surface
automorphism or algebraic self-map.  The certificate explicitly computes the
rational matrix of `z`, verifies that `z` and `z^{-1}` preserve `G_T`, and
checks (4.11) as a matrix identity.

Changing the period did not change `N`, its roots, or the field action on
`NS`.  The previously certified complete-root/Weyl framing therefore returns
the same compact eleven-facet Bugaenko `H8` chamber.  This proves sharpness of
(3.2) inside the Hodge-coherent algebraic-correspondence class.

## 5. What the original NS-only `H9` branch actually contains

The stronger theorem must not be retrofitted onto the weaker ledger.  The
current rank-twenty packages are:

| NS-only package | binary K3 complement | complete geometric-wall status |
|---|---|---|
| silver `Q(sqrt2)` carrier | `<2>+<2>` | the complete **intrinsic field-norm-two** reflection chamber has infinite volume; saturation adds a `(-2)` divisor of field norm `3-sqrt2`, and the complete restricted K3 `(-2)` arrangement is not closed |
| golden discriminant-11 carrier | `[[2,1],[1,6]]` | exact carrier and nonempty Weyl-conjugate eigen-nef marking; 97-wall norm-two prefix still has a hyperideal escape; complete root arrangement and volume remain open |
| golden discriminant-19 carrier | `[[2,1],[1,10]]` | exact singular-K3 carrier; 22-wall bounded prefix has a spacelike escape; complete root arrangement and volume remain open |

For any **fixed** rank-twenty trace lattice `N`, positive even binary forms of
determinant `|det N|` and discriminant anti-isometries give a finite exhaustive
K3 complement/glue list.  There is no finite global list until the input
polytope class supplies determinant/dotted-product bounds.  Likewise, the
`K3^[k]` NS-only branch ranges over unbounded deformation parameter `k`; for
fixed `(N,k)` it is a ternary complement/glue and complete-MBM problem, but
the union over all `k` is not presently a finite exhaustion.

Therefore the honest weak-branch verdict is:

```text
no new compact H9 survivor;
no exhaustive NS-only H9 obstruction;
silver norm-two closed negatively, full K3 (-2) refinement open;
golden discriminants 11 and 19 still open at complete-wall/volume.
```

## 6. Reproduction and sources

Run

```sh
python3 scripts/certify_bugaenko_hodge_coherent_carrier.py
python3 -m pytest -q tests/test_bugaenko_hodge_coherent_carrier.py
```

The frozen certificate is
[`data/candidates/bugaenko_hodge_coherent_carrier.json`](../data/candidates/bugaenko_hodge_coherent_carrier.json),
with payload digest

```text
69bb7458d1fc58972a08a8366f58c4c775ac2e77b7beb3715a0dc207bfe7a34d
```

The Hodge endomorphism-field input is Zarhin's theorem and van Geemen's
explicit treatment of real multiplication; in particular van Geemen's
Lemma 3.2 isolates the small-multiplicity phenomenon:
[van Geemen, *Real multiplication on K3 surfaces and Kuga--Satake
varieties*](https://arxiv.org/abs/math/0609839).
The algebraicity step uses the stated self-product theorem in
[Ramón-Marí, *On the Hodge conjecture for products of certain
surfaces*, Collect. Math. 59 (2008), 1--26](https://doi.org/10.1007/BF03191179),
with the [arXiv source](https://arxiv.org/abs/math/0505357), and the
isometry-spanning cross-check uses
[Buskin, *Every rational Hodge isometry between two K3 surfaces is
algebraic*](https://arxiv.org/abs/1510.02852).

The same-polytope root and volume inputs remain the previously frozen
[`bugaenko_k3_eigen_nef.md`](bugaenko_k3_eigen_nef.md).  The weak silver
scope is recorded literally in
[`silver_k3_carrier.md`](silver_k3_carrier.md) and the separate norm-two
result in
[`silver_k3_norm2_infinite_stabilizer.md`](silver_k3_norm2_infinite_stabilizer.md).

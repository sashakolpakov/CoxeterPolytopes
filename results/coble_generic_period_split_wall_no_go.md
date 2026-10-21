# Generic 2-elementary periods do not realize one-boundary Coxeter walls

**Status (2026-08-09): proved no-go theorem and exact finite replacement.**

This note closes one tempting but false route to essential surjectivity for
the surface-free denominator `D_1_split`.  The invariant lattice

\[
              M=N(2),\qquad N\simeq I_{1,10},
\]

does produce a one-boundary strong Coble quotient at a very general
2-elementary K3 period.  It does **not** produce any of the twelve required
`(-2)` walls, and therefore cannot return the input polytope.  Exact return
requires a special Noether--Lefschetz period carrying a finite system of
anti-invariant algebraic classes.

The positive replacement is the already proved finite compiler
`C_1_split(P,h)`: its nonemptiness is equivalent to existence of an exact
one-boundary Coble carrier and each passing branch returns literally the
same `P`.  What remains open is class-wide nonemptiness on `D_1_split`.

## 1. Generic-period no-go theorem

Let `(P,h)` be any object of `D_1_split`.  Thus `N=I_(1,10)` has a unique
characteristic norm-`-1` boundary direction `b`, and `b^perp=E10` contains
twelve norm-`-2` wall classes spanning it over `Q`.  Put `M=N(2)`.
Its 2-elementary invariants are

\[
                         (r,a,\delta)=(11,11,1).
\]

> **Theorem 1 (generic invariant-lattice obstruction).**  Let `(X,tau)` be
> a very general complex 2-elementary K3 surface with invariant lattice
> `M`.  Then `NS(X)=M`, the quotient `S=X/<tau>` is a strong Coble surface
> with exactly one boundary component, and `S` contains no `(-2)`-curve.
> Consequently none of the twelve root facets of `P` is realized and
> `P(Nef(S))` is not `P`.

### Proof

The anti-invariant period domain has dimension `20-r=9`.  Outside the
countable union of rational Noether--Lefschetz hyperplanes, its only
algebraic classes are the prescribed invariant classes, so `NS(X)=M` and
`tau` acts as `+1` on all of `NS(X)`.

Nikulin's fixed-locus formula gives

\[
 g={22-r-a\over2}=0,\qquad k={r-a\over2}=0.
\]

Hence the fixed locus is one rational curve.  The quotient theorem for a
tame anti-symplectic K3 involution with rational fixed locus makes `S` a
strong Coble surface with one boundary curve `B`.

Suppose that `C` were a `(-2)`-curve on `S`.  Adjunction gives
`K_S.C=0`, hence `B.C=(-2K_S).C=0`; the two distinct effective curves are
disjoint.  The K3 double cover is therefore etale over `C=P1`.  Since `P1`
has no connected etale double cover over `C`,

\[
       pi^{-1}(C)=C^+\sqcup C^-,\qquad
       (C^+)^2=(C^-)^2=-2,\qquad tau(C^+)=C^-.
\]

The nonzero class

\[
                         delta=[C^+]-[C^-]
\]

is algebraic and anti-invariant.  It is nonzero because the two curves are
disjoint while each has square `-2`.  This contradicts `NS(X)=M`, on which
`tau` is the identity.  Thus no such `C` exists.  In particular the twelve
specified root walls cannot occur.  QED.

The theorem is uniform in the raw polytope: choosing a generic period fails
for **every** nontrivial `D_1_split` input, not merely for an exceptional
diagram.  An exact carrier must lie on a Noether--Lefschetz locus.

## 2. The first obstruction is entirely polytope-side

Let `u_1,...,u_12` be the twelve downstairs roots and put

\[
                         d_{ij}=u_i.u_j\quad(i\ne j).
\]

If an exact carrier exists, its K3 cover splits each wall as

\[
 pi^*u_i=c_i^++c_i^-,\qquad
 (c_i^\pm)^2=-2,\quad c_i^+.c_i^-=0,\quad
 tau(c_i^+)=c_i^-.
\]

Choose one sheet over each wall and define `delta_i=c_i^+-c_i^-`.  If

\[
 a_{ij}=c_i^+.c_j^+,
\]

then nonnegativity of intersections between distinct irreducible curves and
the projection formula force

\[
 0\le a_{ij}\le d_{ij},\qquad
 c_i^+.c_j^-=d_{ij}-a_{ij}.
\]

Thus the anti-difference Gram matrix is selected from the finite alphabet

\[
 A_{ii}=-4,\qquad
 A_{ij}=4a_{ij}-2d_{ij},\qquad
 a_{ij}\in\{0,1,\ldots,d_{ij}\}.                 \tag{2.1}
\]

Exchanging the two sheets above wall `i` changes `delta_i` to `-delta_i`,
so assignments are quotiented exactly by diagonal sign switching.

There are immediate, surface-free rejection tests.

1. `A` must be negative semidefinite.  Its nonzero quotient is contained in
   the negative-definite anti-invariant algebraic space orthogonal to the
   invariant ample class.
2. `rank(A)<=9`, since the invariant lattice has rank eleven and a complex
   K3 surface has Picard rank at most twenty.
3. Write the `u_i` in an integral `E10` basis as the rows of `W`.  Every
   congruence

   \[
           zW=0\pmod2
   \]

   forces `(sum z_i delta_i)/2` into the saturated anti-invariant lattice:
   it is the difference between `sum z_i c_i^+` and the integral invariant
   half-sum.  Exact HNF saturation of all these forced half-sums must remain
   integral and negative definite.
4. The saturated anti-lattice may contain no vector of square `-2`.  For
   such a vector `v`, K3 Riemann--Roch makes one of `v,-v` effective, while
   the invariant ample class has intersection zero with it, a contradiction.

Failure at any of these four tests is a rigorous obstruction expressed only
in the marked Gram data of `P`.  Passing them is not sufficient.

## 3. The complete finite compiler and its exact logical strength

For clarity, the full arithmetic compiler `C_1_split(P,h)` is:

1. enumerate every matrix (2.1), modulo every sheet switch;
2. apply exact inertia, radical compatibility, rank, forced-half-sum HNF
   saturation, and anti-rootlessness;
3. enumerate every compatible even equivariant finite-index closure through
   the finite discriminant groups;
4. enumerate primitive K3 complements, discriminant anti-isometries and
   graph glues, and retain only branches whose full saturated fixed lattice
   is the primitive marked `M=N(2)` of type `(11,11,1)`;
5. use the supplied strict interior vector `h` and a complete
   positive-definite short-vector enumeration to test root regularity and
   simultaneous irreducibility of all twenty-four split roots; and
6. invoke the K3 period map and global Torelli on each passing lattice
   branch, then the one-boundary Coble and negative-curve/cusp theorems to
   return the exact original nef chamber.

All discrete choices are finite.  At anti-rank nine the transcendental
complement is a positive even binary form of fixed determinant, so reduced
binary-form enumeration is finite.  At smaller anti-rank the standard
finite discriminant-form/genus and primitive-embedding algorithms apply.
The last root test is finite because the ample class turns it into a
positive-definite bounded-vector problem.

The exact theorem already proved in
[`coble_one_boundary_raw_k3_promotion.md`](coble_one_boundary_raw_k3_promotion.md)
is

\[
 \boxed{\quad
 (P,h)\text{ has an exact one-boundary strong-Coble carrier}
 \iff C_{1\_split}(P,h)\ne\varnothing .
 \quad}                                                   \tag{3.1}
\]

Every passing branch constructs the surface and returns the same `P`; no
genericity shortcut is used.  For the two calibrated inputs the exhaustive
anti-Gram layers collapse as follows.

| raw chamber | assignments / switch orbits | unique rootless survivor | K3 Picard rank |
|---|---:|---|---:|
| Type I, 22 full walls | `27,648 / 18` | `E8(-2)+<-4>` | 20 |
| Type II, 30 full walls | `32,768 / 16` | `D9(-2)` | 20 |

The two censuses are independently replayed in
[`enriques_type_i_downstairs_anti_gram_census.md`](enriques_type_i_downstairs_anti_gram_census.md)
and
[`enriques_type_ii_downstairs_anti_gram_census.md`](enriques_type_ii_downstairs_anti_gram_census.md).

## 4. Honest score of the inverse theorem

With the strengthened five-part contract, the status is:

```text
C          = independent surface-free D_1_split denominator;
Car        = finite exhaustive split-root/equivariant-K3 compiler;
Return     = exact same-P return on every passing branch;
AG         = Kondo classifies a returned complex finite-Aut one-boundary Coble;
Pushback   = the compiler-positive locus has exactly the Type-I and Type-II P's.
```

The missing implication is

\[
          P\in D_{1\_split}\quad\Longrightarrow\quad
          C_{1\_split}(P,h)\ne\varnothing.              \tag{4.1}
\]

Neither the generic period nor the two positive calibrations prove (4.1).
Theorem 1 proves that the generic-period argument cannot prove it.  Therefore
the current result is a **uniform non-singleton L3 exact-return iff and AG
classification of the carrier-positive locus**, not strict L4 essential
surjectivity on all of `D_1_split`.

## 5. Decisive next finite theorem test

Use the independently stated subdenominator

```text
D_1_split^nodot:
    P in D_1_split and every pair of its twelve E10 roots has
    intersection 0, 1, or 2.
```

This contains both calibrated chambers and has a finite pair alphabet.  A
certificate must enumerate its marked twelve-root Grams up to permutation,
then apply, in this order:

1. exact `E10` integral saturation and signature;
2. exact finite-volume/no-hyperideal and complete cusp-lift tests;
3. the complete anti-Gram and equivariant-K3 compiler above.

There are only two logically meaningful outcomes.

* If every raw survivor has a passing branch, (4.1) holds on this
  non-singleton subclass; Kondo then pushes back and proves that its members
  are exactly Types I and II.  This is strict L4 for the scoped class.
* If one certified finite-volume survivor has no branch, its failed earliest
  gate is a new polytope-only obstruction and disproves automaticity even in
  the nodot class.

A preliminary nauty lane with connected twelve-node support, underlying
degree at most three, exactly thirteen nonzero pairs, and at most three
parallel pairs found `100` exact Lorentzian corank-two weighted graphs; `38`
have saturated radical-quotient determinant `4`.  These counts only size the
denominator.  They have **not** passed the `E10`, face, finite-volume, cusp,
or lift gates and are not examples.  The scout is
[`scout_coble_nodot_weighted_graphs.py`](../scripts/scout_coble_nodot_weighted_graphs.py).

## Sources

- V. Alexeev and V. Nikulin,
  [*Classification of log del Pezzo surfaces of index at most 2*](https://arxiv.org/abs/math/0406536),
  for 2-elementary K3 involutions and their quotient correspondence.
- S. Kondo,
  [*Coble surfaces with finite automorphism group*](https://arxiv.org/abs/2102.00591),
  especially the nine-class theorem and the one-boundary Types I and II.
- The project
  [`K3-cover bridge`](coble_k3_cover_polytope_to_surface_bridge.md), for the
  fixed-rational-locus quotient theorem and the orbit-root taxonomy.
- The project
  [`raw one-boundary promotion theorem`](coble_one_boundary_raw_k3_promotion.md),
  for the finite compiler, primitive K3 glues, simultaneous irreducibility,
  and branchwise exact return.

# Hostile audit of the coherent full-`H^2` bound and Bugaenko return

**Status (2026-08-09): pass, with sharply limited scope.**

The coherent branch contains two correct results:

1. a unital totally-real action on the whole weight-two Hodge structure
   forces the uniform slot inequality

   \[
                         d(n+3)\le b_2;
   \]

2. Bugaenko's classical compact `H8` chamber has an exact K3 realization at
   equality, for which the quadratic action is rationally Hodge-coherent and
   is induced by a rational algebraic correspondence.

The second statement is a genuine `polytope -> K3 -> same polytope` return.
It is not a new abstract hyperbolic polytope.  The first statement is a bound
only for the Hodge-coherent carrier category: no essential-surjectivity
theorem currently puts every compact arithmetic Coxeter polytope, or every
NS-only eigencarrier, in that category.

## 1. Audit verdict

| gate | verdict | exact reason |
|---|---|---|
| period multiplicity | pass | the real period plane lies in one real eigenspace of the transcendental `K`-module, so that eigenspace has dimension at least two |
| K3/OG10 numerical bounds | pass | `d(n+1)` divisor slots plus at least `2d` transcendental slots give `d(n+3)<=b2` |
| Bugaenko complement | pass | `T=<-2>+<2>+U` has the displayed self-adjoint square-root-of-five action with a positive two-plane |
| discriminant gluing | pass after literal ambient replay | `A_N=A_T=(Z/2)^2`; both square-root actions are the identity on the discriminant modules, and the resulting operator is an integral self-adjoint `22x22` operator on `Lambda_K3` |
| exact Picard period | pass | the period plane has irrational orthogonal complement, hence no nonzero rational `(1,1)` class occurs in `T` |
| transcendental endomorphism field | pass | exact Picard makes `T_Q` irreducible; the four independent operators `I,B,J,BJ` form a quartic CM field on a four-dimensional rational space, so Zarhin's endomorphism field is exactly that field |
| algebraic correspondence | pass over `Q` | Ramón-Marí proves the Hodge conjecture for all powers of a K3 when `End_Hdg(T)` is CM; the NS summand is already algebraic by products of divisor classes |
| same-polytope return | pass by independent certificate | the special period changes neither `NS`, its complete restricted root arrangement, nor the eleven-facet compact Vinberg chamber; Weyl re-framing conjugates the full Hodge action and returns the same marked metric polytope |
| universality for arbitrary compact arithmetic polytopes | **not proved** | the Hodge-coherent carrier hypothesis is additional and is known to be stronger than an NS-only divisor action |

The earlier executable recorded the two discriminant actions and the glue
conclusion as prose/Boolean outputs.  The strengthened certificate now
computes the index-four `N+T` inclusion in the K3 lattice, conjugates the
block operator through it, and checks integrality, self-adjointness, and the
square relation on all 22 coordinates.  This closes the only material
certificate gap found by this audit.

## 2. The period-multiplicity bound

Let a totally real field `K` of degree `d` act unitally on
`H^2(X,Q)` by Hodge endomorphisms, self-adjoint for the polarization.  The
action preserves `NS(X)_Q`; self-adjointness then preserves
`T(X)_Q=NS(X)_Q^perp`.  Both nonzero spaces are faithful `K`-modules: the
kernel of either restricted unital representation would be an ideal in the
field `K`.

The action on the complex line `H^(2,0)` is a field embedding
`sigma_0:K->C`.  Since `K` is totally real, complex conjugation has the same
character on `H^(0,2)`.  Therefore

\[
 \langle\operatorname{Re}\omega,\operatorname{Im}\omega\rangle_{\mathbf R}
       \subset T_{\sigma_0}\subset T_{\mathbf R}.
\]

If `m=dim_K T_Q`, every real eigenspace has real dimension `m`; hence
`m>=2`.  A faithful rank-`n+1` divisor block and this transcendental block
give

\[
   b_2=d\dim_K H^2_Q\ge d(n+1)+2d=d(n+3).
\]

The consequences used in the ledger are exact:

- K3, `b2=22`, `d=2`: `n<=8`;
- OG10, `b2=24`, `d=2`: `n<=9`;
- `K3^[k]`, `k>=2`, `b2=23`: no unital quadratic action on all of `H2_Q`.

At a hypothetical coherent quadratic K3 `H9` boundary, `NS` already has
rank 20 and `T` has rank two, so `T` would have `K`-rank one.  That cannot
contain the period plane.

There is an important quantifier distinction.  If the **entire** Zarhin
field `End_Hdg(T)` is totally real, van Geemen's stronger result gives
`dim_K T>=3`.  The present theorem assumes only an action by a totally-real
subfield.  Bugaenko attains `dim_K T=2` precisely because this quadratic
field sits properly inside a quartic CM endomorphism field.  There is no
contradiction.

## 3. Literal discriminant and K3-glue check

The already certified Bugaenko trace lattice and its complement are

\[
 N=\langle2\rangle+\langle-2\rangle+E_8(-1)^2,
 \qquad
 T=\langle-2\rangle+\langle2\rangle+U.
\]

The canonical primitive embedding used by the project puts the two rank-one
summands into the first two hyperbolic planes and the two `E8(-1)` summands
identically into `Lambda_K3`.  Appending the displayed basis of `T` gives a
`22x22` matrix of determinant four.  Thus `N+T` has index four, exactly the
size required by

\[
                         A_N\simeq A_T\simeq(\mathbf Z/2)^2.
\]

On `N`, let `tau_N` be multiplication by `phi` and put
`B_N=2tau_N-I`.  Since `tau_N` preserves `N` and its dual, for every
`x in N^vee`

\[
 B_Nx\equiv-x\equiv x\pmod N,
\]

where the last equality uses the exponent-two discriminant group.  On `T`,

\[
B_T=
\begin{pmatrix}
-3&2&0&0\\-2&3&0&0\\0&0&0&5\\0&0&1&0
\end{pmatrix}
\]

sends each of the two half-generators of `A_T` to itself modulo `T`.
Therefore `B_N` and `B_T` agree across every discriminant anti-isometry.

The hostile replay does not stop at that abstract argument.  If `S` is the
index-four basis matrix with columns the embedded bases of `N` and `T`, it
forms

\[
 B_{K3}=S\,(B_N\oplus B_T)\,S^{-1}
\]

and checks directly that

\[
 B_{K3}\in M_{22}(\mathbf Z),\qquad B_{K3}^2=5I,
 \qquad B_{K3}^{t}G_{K3}=G_{K3}B_{K3}.
\]

Thus `Z[sqrt(5)]` acts integrally on the full K3 lattice.  The maximal-order
generator `(I+B_K3)/2` is rational and self-adjoint but not integral.  The
claimed coherent action is accordingly a `Q(sqrt(5))`-action on
`H^2_Q`, with precisely the stated smaller integral order; no maximal-order
integrality is being smuggled into the theorem.

## 4. The special period and exact Picard lattice

The complement calculation uses

\[
 G_TB_T=
 \begin{pmatrix}
 6&-4&0&0\\-4&6&0&0\\0&0&1&0\\0&0&0&5
 \end{pmatrix},
\]

whose leading minors are `6,20,20,100`.  Since `B_T` is self-adjoint and
`B_T^2=5I`, its `+sqrt(5)` eigenspace is a positive-definite real two-plane
and its `-sqrt(5)` eigenspace is negative definite.

The commuting operator

\[
J=
\begin{pmatrix}
0&0&1&-5\\0&0&-1&-5\\-10&10&0&0\\2&2&0&0
\end{pmatrix}
\]

satisfies `J^2=-20I` and is skew-adjoint.  Writing `s=sqrt(5)`, the executable
certificate uses

\[
 u=(2,3+s,0,0),\qquad v=(0,0,-5-s,-1-s),\qquad
 \omega=u+iv.
\]

It verifies exactly

\[
 Bu=su,\quad Bv=sv,\quad Ju=-2sv,\quad Jv=2su,
\]

and

\[
 u^2=v^2=20+12s,\quad u\cdot v=0,
 \quad \omega^2=0,
 \quad \omega\cdot\bar\omega=40+24s>0.
\]

Thus the joint line

```text
B_T = +sqrt(5),       J = +i*sqrt(20)
```

therefore obeys the K3 period equations.  Its real plane is the entire
positive eigenspace.  Its orthogonal complement is the negative irrational
eigenspace of the rational matrix `B_T`; a nonzero rational vector cannot be
an eigenvector with irrational eigenvalue.  Hence there is no extra rational
`(1,1)` class and `NS(X)=N` exactly.

This exact-Picard condition also makes the transcendental Hodge structure
irreducible.  The matrices `I,B_T,J,B_TJ` are rationally independent and
generate

\[
 E=\mathbf Q(\sqrt5,\sqrt{-20})=\mathbf Q(\sqrt5,i).
\]

Zarhin's endomorphism algebra is a field; since it contains this degree-four
field and acts on the four-dimensional rational space `T_Q`, it equals `E`.
This is the missing justification behind the otherwise too-quick inference
from four independent matrices to the phrase “exactly a quartic CM field.”

## 5. What “algebraic correspondence” means here

Ramón-Marí's published theorem states that if
`End_Hdg(T(X))` is a CM field, then the Hodge conjecture holds for every
self-power of the K3 surface.  In particular it holds for `X x X`.  Under
Poincare duality, the rational Hodge endomorphism of `T_Q` defined by `B_T`
corresponds to a rational `(2,2)` Hodge class on `X x X`, hence to a rational
codimension-two algebraic cycle.  The prescribed action on `NS_Q` is
independently algebraic: products `D_i x D_j` of divisor classes, combined
using the inverse intersection matrix, realize arbitrary rational
endomorphisms of `NS_Q`.  Adding the two pieces gives the required action on
all of `H^2_Q`; rational linear combination also gives `(I+B)/2`.

Buskin's theorem gives an independent algebraicity cross-check that does not
mistake `B` itself for an isometry.  Put `i_T=BJ/10`, `t=1+sqrt(5)`, and

\[
 z=(t-i_T)(t+i_T)^{-1}.
\]

Then `z` and `z^{-1}` are rational Hodge isometries and

\[
 \sqrt5={29\over8}(z+z^{-1})-{15\over4}.
\]

Extend these isometries by the identity on `NS_Q`; Buskin algebraicizes their
classes, and divisor-product correspondences correct the `NS_Q` action to
`B_N`.  This recovers the same rational correspondence while keeping the
transcendental CM field separate from the unital full-`H^2` quadratic action.

This proves an algebraic correspondence in the ordinary cohomological sense.
It does **not** prove:

- an integral Chow correspondence;
- a canonical or effective cycle;
- a Hodge isometry (`B^tGB=5G`, not `G`);
- an algebraic self-map or surface automorphism; or
- a correspondence preserving the ample cone as a map.

None of those stronger statements is needed for the Hodge-coherent carrier
definition.

The primary result is J. J. Ramón-Marí, *On the Hodge conjecture for products
of certain surfaces*, Collectanea Mathematica 59 (2008), 1--26,
[DOI](https://doi.org/10.1007/BF03191179), Theorem 4.3 in the published
numbering (the theorem labelled `k3` in the arXiv source).  The arXiv source
used for this audit is [math/0505357](https://arxiv.org/abs/math/0505357);
its downloaded gzip and decompressed TeX SHA-256 hashes are respectively

```text
6a0dc93a25d63a5c049c9bcc3a28cc37a4f0a36532ba05dd6aa85408ceaedb1b
eeb6a89c3d4c1b07c438f436a7bcedcf160d7a36881291f0040d2d5a7ddd1070
```

The isometry cross-check is N. Buskin, *Every rational Hodge isometry between
two K3 surfaces is algebraic*, J. Reine Angew. Math. 755 (2019), 127--150,
[arXiv:1510.02852](https://arxiv.org/abs/1510.02852),
[DOI](https://doi.org/10.1515/crelle-2017-0027).

## 6. Return of the same Bugaenko chamber

The CM period lies in `T`; it does not alter the integral Picard lattice,
the quadratic action on it, or its roots.  The independent Bugaenko
eigen-nef certificate already proves that the complete restricted K3
`(-2)` arrangement equals the field norm-two arrangement and that the
selected chamber is the compact eleven-facet `H8` polytope with 41 finite
vertices.  A K3 Weyl re-framing sends the ambient root chamber to the nef
chamber.  Root reflections act trivially on `A_N`, so that re-framing extends
over the K3 glue; conjugating the full Hodge action leaves the period fixed
and returns the same metric Coxeter polytope.

The logical dependency is therefore

```text
Bugaenko marked H8 datum
  -> exact rank-18 Picard embedding and complete restricted roots
  -> explicit rank-4 coherent CM period and full-lattice glue
  -> projective exact-Picard K3 + rational algebraic correspondence
  -> Weyl-framed nef eigensection equal to the original P.
```

This is a positive inverse-carrier theorem and makes the coherent quadratic
K3 bound sharp.  It neither discovers a new hyperbolic polytope nor supplies
the missing universality implication that would turn the bound into a
worldwide `n<=8` theorem.

## 7. Reproduction

The exact coherent-complement/glue replay and frozen payload are

- [`certify_bugaenko_hodge_coherent_carrier.py`](../scripts/certify_bugaenko_hodge_coherent_carrier.py);
- [`bugaenko_hodge_coherent_carrier.json`](../data/candidates/bugaenko_hodge_coherent_carrier.json); and
- [`test_bugaenko_hodge_coherent_carrier.py`](../tests/test_bugaenko_hodge_coherent_carrier.py).

Run

```bash
python3 scripts/certify_bugaenko_hodge_coherent_carrier.py
python3 -m pytest -q tests/test_bugaenko_hodge_coherent_carrier.py
```

The independent same-polytope input is
[`bugaenko_k3_eigen_nef.md`](bugaenko_k3_eigen_nef.md).  The general bound is
[`coherent_full_h2_real_multiplication_bound.md`](coherent_full_h2_real_multiplication_bound.md),
and the sharp return theorem is
[`k3_h9_hodge_coherent_return_theorem.md`](k3_h9_hodge_coherent_return_theorem.md).

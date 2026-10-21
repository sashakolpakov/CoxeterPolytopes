# Polytope-to-eigencarrier Hodge-rank return theorem

**Status (2026-08-09).**  This note supplies the missing converse rank layer
for the existing totally-real eigencarrier criterion.  It applies to a
faithful field module contained in the divisor lattice, not only when the
trace lattice is the whole Picard lattice.  It also states the exact
no-hidden-wall condition needed to return the **same** polytope.

This is a carrier-class theorem.  It becomes a theorem about every polytope
in a surface-free class only after proving that the class has an exhaustive
compiler into one of these carriers.

The hypotheses and exact-cone conventions have been independently checked in
the
[`hostile audit`](polytope_surface_inverse_theorems_hostile_audit.md), with a
qualified-pass verdict after explicit repairs.

## 1. Hodge-rank return

Let `X` be either a smooth complex projective surface, with the intersection
form on `NS(X)/tors`, or a projective irreducible holomorphic symplectic
manifold, with the Beauville--Bogomolov form.  Write this lattice as `L` and
its rank as `rho`.  In both cases

\[
                    \operatorname{sig}L=(1,\rho-1).       \tag{1.1}
\]

Let `K` be a totally real field of degree `d`, and let an order `A` in `K`
act faithfully and self-adjointly on a primitive sublattice
\(N\subset L\).  Thus `N` is `A`-stable, and the action extends to the
field `K` on `N_Q`.
Assume that

\[
 N_{\mathbf Q}\quad\hbox{is a free }K\hbox{-module of rank }m=n+1,
 \tag{1.2}
\]

and that the eigenspace `E_0` for one real embedding `sigma_0` has signature
`(1,n)`.  No action by automorphisms of `X` is assumed.

### Theorem 1.1 -- faithful eigencarrier rank theorem

Under these hypotheses:

1. `rank_Z N=d(n+1)`, and therefore

   \[
   \boxed{d(n+1)\leq\rho(X).}                            \tag{1.3}
   \]

2. There is an orthogonal real decomposition

   \[
   N_{\mathbf R}=\bigoplus_{\sigma:K\hookrightarrow\mathbf R}E_\sigma,
   \qquad \dim E_\sigma=n+1.                             \tag{1.4}
   \]

   Every `E_sigma` other than `E_0` is negative definite.  Thus `N` has
   signature `(1,d(n+1)-1)`, and (N^\perp\cap L) is negative definite.

3. The rational intersection form uniquely recovers a symmetric
   `K`-bilinear form `h` by

   \[
   (ax,y)=\operatorname{Tr}_{K/\mathbf Q}(a h(x,y))
   \quad(a\in K).                                       \tag{1.5}
   \]

   Up to the harmless global sign convention, `h` is Lorentzian at
   `sigma_0` and definite at every conjugate place.  Hence the Hodge index
   theorem forces precisely the admissibility condition used in arithmetic
   hyperbolic reflection theory; it is not merely a convenient input to the
   surface construction.

4. If `N=NS(X)/tors`, then `rho=d(n+1)`.  For a K3 or a fixed
   hyperkaehler deformation type with Beauville rank `b2`, a primitive
   embedding of `N` has a marked lattice-polarized period locus of dimension

   \[
                     b_2-d(n+1)-2.                      \tag{1.6}
   \]

   Away from the countable union of rational Noether--Lefschetz
   hyperplanes its Picard lattice is exactly `N`; the positive class in `N`
   makes that generic carrier projective.  In particular, if
   `d(n+1)=b2-2`, the orthogonal complement is the positive period plane,
   the period locus is zero-dimensional, and a primitive full-Picard
   realization has no divisor classes outside `N` capable of adding hidden
   walls.

#### Proof

The first assertion is restriction of scalars.  For two distinct embeddings
`sigma,tau`, choose `a in K` with `sigma(a) != tau(a)`.  Self-adjointness
gives

\[
(\sigma(a)-\tau(a))(E_\sigma,E_\tau)=0,
\]

which proves (1.4).  Since `E_0` already contains the unique positive
direction allowed by (1.1), every other summand and the orthogonal complement
are negative definite.  This proves (1.3) and the signature assertion.

For fixed `x,y`, the map (a\mapsto (ax,y)) is a rational linear functional
on `K`.  Nondegeneracy of the field trace gives a unique `h(x,y)` satisfying
(1.5).  Self-adjointness and symmetry make `h` symmetric and `K`-bilinear.
The real-place signatures just proved give arithmetic admissibility.
Finally, the period line lies in `N^perp`, of signature
`(2,b2-d(n+1)-2)`; the usual period-surjectivity, generic
Noether--Lefschetz avoidance, and projectivity arguments give (1.6) and the
exact-Picard assertion.  `square`

### A discriminant sharpening at the K3 boundary

If `N` embeds primitively in the unimodular K3 lattice, with
`r=d(n+1)`, its complement has rank `22-r` and opposite discriminant form.
Consequently

\[
               \ell(A_N)\leq22-d(n+1).                 \tag{1.7}
\]

For a real-quadratic `H9` carrier this is the sharp binary-complement gate
`ell(A_N)<=2`; exact existence is equivalently the existence of a positive
rank-two even lattice with discriminant form `-q_N`.  This is stronger than
the rank inequality.  For a maximal real-quadratic OG10 `H10` carrier, the
orthogonal complement is again positive rank two, but its glue must leave
the order-three discriminant form of the OG10 Beauville lattice; the full
rank-two complement/glue enumeration is the exact gate.

## 2. Deformation-type ceilings

For a projective irreducible symplectic manifold, `rho<=b2-2`.  The theorem
therefore gives

\[
 \boxed{n\leq\left\lfloor\frac{b_2-2}{d}\right\rfloor-1.} \tag{2.1}
\]

For the classical smooth complex carrier types this specializes as follows.

| carrier type | divisor-rank ceiling `R` | general ceiling | real-quadratic ceiling |
|---|---:|---:|---:|
| K3 surface | 20 | `floor(20/d)-1` | `H9` |
| Enriques numerical lattice | 10 | `floor(10/d)-1` | `H4` |
| `K3^[k]`, `k>=2` | 21 | `floor(21/d)-1` | `H9` |
| generalized Kummer, dimension at least four | 5 | `floor(5/d)-1` | `H1` |
| OG6 | 6 | `floor(6/d)-1` | `H2` |
| OG10 | 22 | `floor(22/d)-1` | `H10` |

Here `R=10` for Enriques means `Num(Y)=U+E8(-1)`; for the hyperkaehler rows
`R=b2-2`, using respectively `b2=23,7,8,24`.  Thus:

- among this entire list, a real-quadratic faithful eigencarrier has
  `n<=10`;
- OG10 is the only listed type capable of reaching `H10`;
- after removing OG10, the ceiling is `H9`;
- generalized Kummer type cannot carry a real-quadratic hyperbolic section
  of dimension at least two.

These are sharp **Hodge-slot bounds**, not claims that the boundary chamber
exists or has finite volume.  The current boundary calibration is exact:

| datum | trace rank / period dimension | return status |
|---|---:|---|
| Bugaenko golden `H8` on a K3 | `18 / 2` | same compact 11-facet eigen-nef chamber |
| silver `H9` on a K3 | `20 / 0` | weak/full-Picard carrier; its complete **intrinsic field-norm-two** chamber has infinite volume, while the full restricted K3 `(-2)` arrangement remains open |
| golden `H10` on OG10 | `22 / 0` | weak/full-Picard carrier, but its complete chamber has infinite volume |

Thus the strong positive real-quadratic frontier remains Bugaenko `H8`.
The OG10 ceiling object fails its audited complete reflection-volume gate.
For the silver K3 object only the intrinsic field-norm-two subarrangement is
closed negatively: saturation adds a `(-2)` divisor of field norm
`3-sqrt(2)`, and the complete geometric restricted-root arrangement is not
known.  Its tested Alexeev local-star sections are non-Coxeter.  The twelve
compact golden Lanner simplices in
`H2`--`H4` give further exact same-polytope K3 returns, but are not
high-dimensional boundary examples.

### Coherent full-`H2` sharpening

The table above permits an action only on the divisor lattice.  If the
degree-`d` totally-real action is instead required to extend as a **unital
self-adjoint Hodge action on all of `H2_Q`**, the transcendental lattice is a
`K`-module and its distinguished real eigenspace must contain the positive
period two-plane.  Therefore `dim_K T>=2`, and the stronger inequality is

\[
                         d(n+3)\leq b_2.                \tag{2.4}
\]

For a real quadratic field this lowers the coherent K3 ceiling to `H8` and
the coherent OG10 ceiling to `H9`; in particular a maximal OG10 `H10`
NS-action cannot extend coherently across its rank-two transcendental
lattice.  Odd-`b2` types cannot carry any unital quadratic action on all of
`H2_Q`.  The known NS-only `H10` OG10 carriers are not contradicted.  See the
[coherent full-cohomology theorem](coherent_full_h2_real_multiplication_bound.md).
On K3 the quadratic bound is sharp: the Bugaenko `H8` package has an explicit
exact-Picard CM period, a computed integral `22 x 22` full-cohomology
`Z[sqrt(5)]` action, and a rational algebraic correspondence, while its
complete restricted-root return remains the same compact eleven-facet
polytope.  See the
[sharp coherent K3 return](k3_h9_hodge_coherent_return_theorem.md) and its
[hostile audit](coherent_full_h2_bugaenko_hostile_audit.md).

## 3. Exact return of the same polytope

Rank and primitive embedding produce only a weak divisor carrier.  The
following criterion isolates the exact missing condition.

Let `D_X` be the **complete** set of oriented geometric wall classes defining
the closure `A_X` of the chosen nef, movable, or Kähler-type chamber, so that
`A_X` is the intersection of their oriented halfspaces.  Fix a point
`p in E_0` in the relative interior of `A_X intersect E_0` and assume it is
regular: no wall in `D_X` contains `p` (in particular, no complete wall
contains all of `E_0`).  Orient every wall to be positive at `p`.  For
`delta in D_X`, let

\[
        \lambda_\delta=(\delta,\,\cdot)|_{E_0}\in E_0^*. \tag{3.1}
\]

Write the oriented cone over the proposed Coxeter polytope as

\[
             P^+=\{x\in \overline{E_0^+}: \ell_i(x)\geq0,
                          \ 1\leq i\leq f\}.             \tag{3.2}
\]

Here `overline(E_0^+)` is the closed chosen positive-cone sheet, including
its lightlike boundary and the origin; the same convention is used for the
closed geometric chamber `A_X`.

### Theorem 3.1 -- same-chamber return criterion

The equality

\[
                    A_X\cap E_0=P^+                    \tag{3.3}
\]

holds if and only if both conditions below hold.

1. **Facet realization.**  For every irredundant `ell_i`, some
   `delta_i in D_X` has `lambda_delta_i` a positive multiple of `ell_i`.
2. **No hidden wall.**  For every `delta in D_X`, including classes outside
   `N` when `N` is a proper Picard sublattice,

   \[
                    \lambda_\delta(x)\geq0
                    \quad\hbox{for every }x\in P^+.     \tag{3.4}
   \]

Equivalently, every complete geometric wall restriction is nonnegative on
`P^+`, while every irredundant input facet occurs among those restrictions.
The proof is just the two inclusions in (3.3), but the criterion is useful
because it separates the exact AG issue from Gram matching.

To return the same **Coxeter** datum one must additionally check that the
restricted facet Gram is the oriented input Gram and that each projected
reflection preserves the intended marked arithmetic module (the
`A`-lattice before choosing a real embedding).  An individual irrational
real eigenspace usually has no discrete full-rank `Z`-lattice of its own.
Neither the arithmetic preservation nor the Gram matching follows from
(3.3).

For a full-Picard K3 carrier, `D_X` is the complete `(-2)`-root arrangement.
Provided `E_0` is not contained in a root mirror and `p` is regular as above,
K3 Weyl transitivity lets one frame the ambient root chamber containing `p`
as the nef chamber.  Hence a component of the complete restricted root
arrangement satisfies (3.3) automatically; only root completeness and the
Coxeter/arithmetic-preservation gates remain.  For hyperkaehler carriers,
`D_X` must be
the complete MBM set (Kähler) or prime-exceptional set (movable), together
with the appropriate monodromy/Kähler-type chamber theorem.  Bounded MBM
squares give finitely many numerical wall types, not automatically a finite
or complete wall list.

The boundary saturation in Theorem 1.1 is real leverage: for K3 `H9` and
OG10 `H10` real-quadratic full-rank carriers, there are no outside Picard
classes, so the no-hidden test is entirely internal to the trace lattice.

## 4. Finite folding is a separate compiler

If a finite lattice group `G` preserves a geometric chamber and `E=L_R^G`
contains a positive class, then for `x in E`

\[
 (x,\delta)=(x,g\delta),\qquad
 (x,\delta)=\frac1{|G\delta|}
       \left(x,\sum_{\eta\in G\delta}\eta\right).       \tag{4.1}
\]

Thus the complete folded cone is obtained by restricting the primitive orbit
sums of **all** geometric wall orbits.  The same criterion of Theorem 3.1 is
necessary and sufficient.  A proposed fold is Coxeter only after separate
angle, irredundancy, and reflection-integrality tests on those projected
normals.  A plain fixed-space fold yields only `n+1<=rho`; the multiplier
`d(n+1)` appears precisely when the rational representation carries the
faithful self-adjoint degree-`d` field structure above.

This paragraph does not require the field action to be a surface
automorphism.  It merely records how a genuinely geometric finite fold, when
one is supplied, differs from an abstract divisor-lattice eigenspace.

## 5. The theorem returned to polytopes, with its quantifier exposed

Let `C(K,T)` be a surface-free class of marked arithmetic Coxeter polytopes
over `K`.  Suppose one proves an exhaustive compiler which sends every
`P in C(K,T)` to at least one faithful carrier of a type in `T`, and suppose
Theorem 3.1 plus the Coxeter/integrality tests return exactly `P`.  Then

\[
       \dim P\leq\max_{T\in\mathcal T}
          \left(\left\lfloor\frac{R_T}{[K:\mathbf Q]}\right\rfloor-1\right).
\tag{5.1}
\]

This is the desired `P -> variety -> theorem on P` implication.  For the six
carrier types in the table and a real-quadratic field it gives `dim P<=10`.
If an exhaustive argument eliminates the OG10 branch, it immediately lowers
the returned bound to nine.

The exact terminal route is therefore:

1. classify all admissible integral normal scalings and saturated trace
   lattices of the surface-free polytope class;
2. use (1.3), (1.7), and the exact Beauville complement/glue problem to
   exhaust its possible carriers;
3. at the saturated boundary, classify the internal root/MBM arrangements
   and apply Theorem 3.1;
4. either find a compact boundary survivor or eliminate every carrier
   package; and
5. only after compiler exhaustivity, return the result as a universal
   polytope theorem.

At present Step 1 is proved only for restricted classes such as the connected
golden normal-scaling class.  Therefore the `H10`/`H9` bounds above are exact
carrier-class boundaries, not universal bounds for compact arithmetic
Coxeter polytopes.

## References and relation to existing reports

The full-Picard primitive-embedding, finite visible-shell, K3 Weyl-framing,
hyperkaehler MBM, and arithmetic stabilizer theorems are in the
[`totally-real eigencarrier criterion`](totally_real_eigencarrier_criterion.md).
The present note adds the arbitrary-sublattice Hodge-rank converse, the
deformation-type table, the exact no-hidden-wall iff, and the conditional
bound returned to a surface-free polytope class.

The period inputs are the standard marked K3 and hyperkaehler period and
projectivity theorems; see [Dolgachev](https://arxiv.org/abs/alg-geom/9502005),
[Huybrechts](https://arxiv.org/abs/alg-geom/9705025), and
[Markman](https://arxiv.org/abs/1101.4606).  The OG6 value `b2=8` is due to
[O'Grady](https://arxiv.org/abs/math/0010187), and the OG10 Beauville lattice
and `b2=24` to
[Rapagnetta](https://arxiv.org/abs/math/0606409).

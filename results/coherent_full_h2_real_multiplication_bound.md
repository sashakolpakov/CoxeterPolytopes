# Coherent full-`H^2` real multiplication removes the OG10 `H10` branch

## Result

There is a sharper Hodge-slot obstruction when the totally-real action is
required to be a genuine action on the **whole weight-two Hodge structure**,
not merely a self-adjoint endomorphism of the Neron--Severi lattice.

Let `X` be a projective irreducible holomorphic symplectic manifold, let

\[
K\hookrightarrow \operatorname{End}_{\mathrm{Hdg}}
       (H^2(X,\mathbf Q))
\]

be a unital action of a totally real field of degree `d`, and suppose the
action is self-adjoint for the Beauville form.  If `NS(X)_Q` contains a
faithful `K`-submodule of rank `n+1` whose distinguished real eigenspace has
signature `(1,n)`, then

\[
             \boxed{d(n+3)\leq b_2(X).}                 \tag{1}
\]

For real quadratic multiplication this gives

| deformation type | `b2` | coherent quadratic ceiling |
|---|---:|---:|
| K3 | 22 | `H8` |
| `K3^[m]`, `m>=2` | 23 | no unital quadratic action on all `H2_Q` |
| generalized Kummer | 7 | no unital quadratic action on all `H2_Q` |
| OG6 | 8 | `H1` |
| OG10 | 24 | `H9` |

In particular:

> A quadratic `H10` eigencarrier on an OG10 variety has Picard rank `22`
> and transcendental rank `2`.  Its real-quadratic action on `NS_Q` cannot
> extend to a unital self-adjoint Hodge action on `H2_Q`.

This uniformly eliminates the maximal OG10 branch from the **coherent
full-cohomology** carrier category.  It lowers the carrier ceiling from
`H10` to `H9` in that category, independently of the Coxeter diagram, normal
scaling, discriminant gluing, or wall enumeration.

It does **not** eliminate Neron--Severi-only eigencarriers.  The existing
golden `H10`, A2-coupled `H10`, and norm-six `H10` OG10 packages are therefore
not contradicted: their displayed real-quadratic action is only an action on
the divisor lattice.  A divisor-supported correspondence that realizes one
operator on `NS_Q` but does not make all of `H2_Q` a unital `K`-module also
lies outside the theorem.

## Proof

Write

\[
T(X)_\mathbf Q=NS(X)_\mathbf Q^\perp\subset H^2(X,\mathbf Q).
\]

A rational Hodge endomorphism preserves `NS(X)_Q`.  Self-adjointness then
makes `T(X)_Q` stable as well.  Because the action of the field is unital,
both stable rational spaces are vector spaces over `K`.  In particular,
their rational dimensions are multiples of `d`.

Let `sigma` span `H^{2,0}(X)`.  Every `a in K` acts on this complex line by a
scalar.  Since `K` is totally real, that scalar is `rho_0(a)` for one real
embedding `rho_0:K->R`.  Hence the real period plane

\[
\Pi=\langle\operatorname{Re}\sigma,
          \operatorname{Im}\sigma\rangle_\mathbf R
\]

lies in the `rho_0`-eigenspace of `T(X)_R`.  If

\[
r_T=\dim_K T(X)_\mathbf Q,
\]

then every real eigenspace of `T(X)_R` has dimension `r_T`.  The containment
of the two-dimensional period plane forces

\[
                         r_T\geq2.                       \tag{2}
\]

The divisor eigencarrier consumes at least `n+1` dimensions over `K`, while
(2) consumes at least two more.  Therefore

\[
b_2(X)=d\dim_KH^2(X,\mathbf Q)
       \geq d((n+1)+2)=d(n+3),
\]

which proves (1).

For maximal quadratic OG10 `H10`, the divisor block already has rational
rank `2(10+1)=22`.  Its orthogonal complement has rational rank two and
would consequently have `K`-rank one.  But its distinguished real
eigenspace would then be one-dimensional and could not contain the period
plane.  This is the promised contradiction.

## Polytope-to-variety meaning

This is a uniform branch obstruction of the required direction:

```text
quadratic H10 Coxeter/trace datum
  -> proposed coherent OG10 carrier
  -> period plane needs two K-slots in T
  -> only one K-slot is available
  -> coherent carrier moduli locus is empty.
```

Thus any exhaustive `P -> variety` compiler whose top-dimensional allowed
species require coherent full-`H2` real multiplication returns the bound
`n<=9`.  At present, however, the project also allows the strictly weaker
and demonstrably nonempty **NS-only** OG10 species.  Until one proves that a
strong carrier must be coherent, or eliminates every NS-only package by its
complete MBM chamber, this is not a universal nonexistence theorem for
compact arithmetic Coxeter polytopes in `H10`.

The argument is elementary Hodge linear algebra.  Its standard context is
the theory of real multiplication on K3-type Hodge structures; see van
Geemen, [Real multiplication on K3 surfaces and Kuga--Satake
varieties](https://arxiv.org/abs/math/0609839), and Zarhin's classification
of endomorphism fields of K3-type Hodge structures.

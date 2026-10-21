# Quadratic-CM and algebraic-correspondence upgrade of the golden Lanner returns

**Status (2026-08-09): TOOT — four exact positive upgrades, with a complete
quadratic-extension audit of all twelve rows.**

The full-`H^2` coherent return theorem has a genuinely algebraic refinement
for four of the twelve golden Lanner chambers:

```text
P -> projective exact-Picard K3 with End_Hdg(T)=E (a quartic CM field)
  -> a rational algebraic correspondence inducing the golden action
  -> the same P.
```

The four survivors are `(2,5,5)`, `[5,3,3,3]`, `[5,3,3,5]`, and
`D5(terminal 5)`.  This is new project-level structure on classical Lanner
simplices, not a claim that the abstract polytopes themselves are new.  The
other eight rows are negative only for a **quadratic CM extension of the
already certified `Q(sqrt(5))` action**; higher-degree CM fields and different
full-cohomology actions are not excluded.

## 1. Finite twelve-row ledger

Put `K=Q(phi)`, `phi^2=phi+1`.  If `m=dim_K T`, the exact outcomes are:

| chamber | dim `P` | `m` | quadratic-CM verdict | exact reason / `c` |
|---|---:|---:|---|---|
| `(2,5,5)` | 2 | 8 | **pass; algebraic correspondence** | `c=8-4phi=(2phi-2)^2` |
| `(3,3,5)` | 2 | 8 | impossible for this `K`-form | `det C=4+2phi` is not a square in `K` |
| `(3,5,5)` | 2 | 8 | impossible for this `K`-form | `det C=6phi^2` is not a square in `K` |
| `(5,5,5)` | 2 | 8 | impossible for this `K`-form | `det C=10phi^2` is not a square in `K` |
| `[3,5,3]` | 3 | 7 | impossible | `2` does not divide `dim_K T=7` |
| `[5,3,5]` | 3 | 7 | impossible | same degree-divisibility obstruction |
| `D4(5,3,3)` | 3 | 7 | impossible | same degree-divisibility obstruction |
| `C4(one 5)` | 3 | 7 | impossible | same degree-divisibility obstruction |
| `C4(opposite 5s)` | 3 | 7 | impossible | same degree-divisibility obstruction |
| `[5,3,3,3]` | 4 | 6 | **pass; algebraic correspondence** | `c=10-6phi=2(2-phi)^2` |
| `[5,3,3,5]` | 4 | 6 | **pass; algebraic correspondence** | `c=14-8phi` |
| `D5(terminal 5)` | 4 | 6 | **pass; algebraic correspondence** | `c=8-4phi=(2phi-2)^2` |

There are no undecided rows in this quadratic-extension audit.

## 2. Exact positive certificate

For each positive row, let `C` be the diagonal `K`-form used by the coherent
return theorem.  The executable constructs a `K`-linear matrix `J` and checks
exactly

\[
 J^2=-cI,\qquad J^tC+CJ=0,\qquad c>0
\]

at both real places of `K`.  It independently expands `J` to a rational
matrix and checks skew-adjointness for the rational trace form.  It also
checks that `J` commutes with multiplication by `phi`, that the determinant
has the required CM square class, and that the distinguished positive
two-plane is `J`-stable.  Hence

\[
 E=K(\sqrt{-c})
\]

is a quartic CM field acting on `T_Q`, with complex conjugation represented
by the adjoint involution.

The construction is small rather than search-dependent.  On the first
two-dimensional diagonal block use

\[
 J_0=\begin{pmatrix}0&-c\\1&0\end{pmatrix}.
\]

For `(2,5,5)`, the remaining six coordinates split into three standard
skew blocks using `sqrt(c)=2phi-2`.  In each rank-six case, the remaining
`<1,1,1,1>` block uses left multiplication by a pure quaternion of norm
`c`.  The exact triples are

\[
 (2-\phi,0,\phi-2),\quad (1,0,2\phi-3),\quad
 (2\phi-2,0,0)
\]

for `[5,3,3,3]`, `[5,3,3,5]`, and `D5(terminal 5)`, respectively.  The full
matrices are emitted by the certifier.

## 3. CM period, exact Picard lattice, and algebraicity

Write `m=2s`.  At the distinguished real place, the `E`-Hermitian form has
signature `(1,s-1)`: the first block above is the positive real two-plane,
and its orthogonal complement is negative.  The positive `+i sqrt(c)`
eigenlines of `J` therefore form a complex ball of dimension

\[
 s-1=3\quad\hbox{for }(2,5,5),\qquad
 s-1=2\quad\hbox{for the three }H^4\hbox{ rows}.
\]

A generic line in this ball avoids two countable collections of proper
subvarieties:

1. the rational orthogonality loci that would add a `(1,1)` class outside
   the frozen Picard lattice; and
2. the eigenspace loci of rational Hodge endomorphisms not in `E`.

Thus one can choose a period with `NS(X)` equal to the certified Picard
lattice and `End_Hdg(T(X))=E` exactly.  The first equality makes `X`
projective and leaves the complete Picard root arrangement unchanged.  The
second is the standard generic-endomorphism statement on an `E`-linear K3
period domain: an extra rational operator imposes a proper locus, while an
operator scalar on the whole domain descends to an element of `E`.  For the
background facts that the transcendental Hodge structure is irreducible and
its Hodge endomorphism algebra is a number field, see Huybrechts,
[*Lectures on K3 Surfaces*, Section 3.3](https://www.math.uni-bonn.de/people/huybrech/K3.html).

Ramón-Marí proves that if the transcendental Hodge endomorphism field of a
K3 surface is CM, then the Hodge conjecture holds for every self-power of
the surface.  In particular, every rational Hodge endomorphism of `T_Q`
gives a rational codimension-two algebraic cycle on `X x X`.  This
algebraicizes multiplication by `phi` on `T_Q`.  On `NS_Q`, arbitrary
rational endomorphisms are realized by rational combinations of products of
divisor classes, using the inverse intersection matrix.  Adding the two
pieces yields a rational algebraic correspondence inducing the full golden
action on `H^2_Q`.

The primary source is J. J. Ramón-Marí, *On the Hodge conjecture for
products of certain surfaces*, Collectanea Mathematica 59 (2008), 1--26,
[arXiv:math/0505357](https://arxiv.org/abs/math/0505357),
[DOI](https://doi.org/10.1007/BF03191179).

This does not produce an algebraic self-map, an automorphism, a canonical or
effective cycle, or an integral Chow correspondence.  The correspondence is
rational and cohomological, which is exactly the level needed for the
full-`H^2` Hodge action.

## 4. Why the eight negative gates are proofs

Suppose a nondegenerate `K`-quadratic space `V` of dimension `2s` admits a
skew-adjoint `J` with `J^2=-c`, where `c` is totally positive.  Choose an
anisotropic `v`.  The plane spanned by `v,Jv` is orthogonal with Gram

\[
 \langle a,ca\rangle,\qquad a=(v,v),
\]

and its orthogonal complement is `J`-stable.  Iteration gives

\[
 \det V\equiv c^s\pmod {K^{*2}}. \tag{4.1}
\]

For the three failed rank-eight triangles, `s=4`, so (4.1) forces the
determinant to be a square.  It is not: `N(4+2phi)=20` is not a rational
square, while squareness of `6phi^2` or `10phi^2` would put `sqrt(6)` or
`sqrt(10)` in `Q(sqrt(5))`.

For every rank-seven row, an action of the quadratic field
`K(sqrt(-c))` would make `T` a vector space over a degree-two extension of
`K`; hence `2` would divide `7`, an impossibility.

These arguments exclude exactly the attempted quadratic CM enhancement
commuting with the frozen golden action.  They do not exclude:

- a CM endomorphism field of larger degree;
- a different rational `K`-action on the same K3 lattice; or
- an algebraicity proof not passing through CM.

## 5. Same-polytope return and reproduction

The CM specialization changes neither the integral Picard lattice nor its
complete restricted `(-2)`-root arrangement.  The upstream certificate
therefore still gives compactness by Vinberg's finite-vertex criterion, no
interior root cuts, and the exact same metric simplex after K3 Weyl framing.
No volume conclusion is inferred from a partial wall list.

Files:

- base coherent theorem: `results/golden_lanner_full_h2_coherent_returns.md`;
- exact CM certifier: `scripts/certify_golden_lanner_cm_upgrade.py`;
- regression tests: `tests/test_golden_lanner_cm_upgrade.py`.

Run

```bash
PYTHONPATH=src:. python3 scripts/certify_golden_lanner_cm_upgrade.py
PYTHONPATH=src:. python3 -m pytest -q tests/test_golden_lanner_cm_upgrade.py
```

The output contains four exact positive matrices, five degree-divisibility
failures, and three determinant-square-class failures.

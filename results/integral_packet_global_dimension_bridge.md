# Integral packets: the exact global-dimension bridge and its local obstruction

**Status (2026-08-10): unconditional carrier theorem, exact local-to-global
no-go, and one finite falsifiable completion lemma.  Universality for all
reflection-maximal arithmetic finite-volume chambers remains open.**

This report uses the project-native interval polynomial theorem, elementary
linear algebra, and explicit algebraic packet bounds.  It does not use a
classification or finiteness theorem for arithmetic reflection groups, an
arithmetic volume formula, or rigidity.

The answer is sharp in logical form:

```text
the interval probe makes the coefficient-field packet finite;
a faithful global cohomological packet makes its rank at least d(n+1);
an independent AG rank/degree ceiling then bounds n;
four-node packet types alone do not supply that global ceiling.
```

## 1. The faithful integral packet

Let `P` be an arithmetic finite-volume Coxeter chamber in `H^n`, with
ground field `k` of degree `d`, and let `L_k` be the `k`-span of its facet
normals.  The normals span the ambient Lorentzian space, so

\[
                         \dim_k L_k=n+1.                 \tag{1}
\]

A **faithful linear cohomological packet** for `P` consists of a finite free
abelian group `H_Z`, an integral pairing and marked integral endomorphisms on
it, together with an injection

\[
 \iota:\operatorname {Res}_{k/\mathbf Q}L_k
             \hookrightarrow H_\mathbf Q                 \tag{2}
\]

which intertwines the marked `k`-action, sends the facet normals to their
marked cohomological classes, and returns their complete Gram pairing at the
distinguished real eigenspace.  It may be a submodule of `NS(S)`, a subquotient
of integral cohomology, or a monodromy local system.  What matters is the
linear injection (2), not the name of the carrier.

> **Theorem 1 (faithful-rank inequality).**  If `R=rank_Z H_Z`, then
>
> \[
>                    \boxed{d(n+1)\le R}.                \tag{3}
> \]
>
> If the packet comes with an independently certified cellular or
> cohomological ceiling `R<=B`, then
>
> \[
>              \boxed{n\le\left\lfloor{B\over d}\right\rfloor-1
>                         \le B-1}.                      \tag{4}
> \]

**Proof.**  Restriction of scalars in (2) has rational dimension
`d(n+1)`.  An injection into `H_Q` cannot increase dimension, giving (3),
and (4) follows.  `square`

There is an equally elementary degree bound for facet markings.  Let a
projective carrier have an integral ample class `A`, one distinct irreducible
curve `D_i` per facet, and

\[
                       \sum_i A\mathbin\cdot D_i\le D.    \tag{5}
\]

Each summand is a positive integer.  Thus

\[
                f\le D,\qquad n+1\le f,qquad n\le D-1. \tag{6}
\]

Equations (3) and (6) are genuinely different.  The first uses faithful
global linear return; the second uses total marked-cycle degree.  A bound on
the degree of each curve separately does not imply (5).

## 2. What the interval probe adds

Assume the packet marks integral totally real spectral coordinates
`alpha_1,...,alpha_t` which generate `k`.  For each coordinate assume one
distinguished conjugate lies in `[2,14]` and every other conjugate lies in
`[-2,1]`.  The exact project-native theorem

\[
 \|[x(x+1)]^4[(x-1)(x+2)]^3[x^2+x-1]^3\|_{[-2,1]}
 ={884736\over9765625}                                  \tag{7}
\]

and its integer-norm calculation give

\[
 [\mathbf Q(\alpha_j):\mathbf Q]\le23,
 \qquad [k:\mathbf Q]\le23^t.                          \tag{8}
\]

If one coordinate is primitive for `k`, then `d<=23`.  Combining with (3)
gives the exact packet inequality

\[
 d(n+1)\le R\le B,
 \qquad d\le23^t.                                      \tag{9}
\]

The upper bound on `d` does **not** by itself bound `n`; the independent
ceiling `B` or `D` is indispensable.  The role of (8) is effectivity.  Every
minimal polynomial has degree at most `23`, or at most `23^t` after a
primitive-element elimination, and its coefficients are bounded from the
conjugate envelope.  Thus the possible integral spectral coordinates form a
finite, directly enumerable set.  See the
[exact interval theorem](fekete_interval_minus2_1_probe_theorem.md) and its
[hostile audit](fekete_interval_minus2_1_probe_hostile_audit.md).

## 3. A project-native finite-packet theorem

The following is the strongest noncircular conditional replacement presently
available.

> **Theorem 2 (bounded faithful packet finiteness).**  Fix integers
> `B,D,C,t`.  Consider chambers admitting faithful packets for which:
>
> 1. `rank H_Z<=B`;
> 2. there are at most `D` facet vectors;
> 3. in a certified integral basis, the pairing matrix, the at-most-`t`
>    marked endomorphism matrices, and all facet-vector coordinates have
>    absolute value at most `C`;
> 4. the distinguished real eigenspace, time sheet, halfspaces, finite/ideal
>    flags, and exact-return rule belong to one fixed finite grammar.
>
> Then only finitely many complete marked Grams, fields, and chambers occur.
> They are exhaustively obtainable by finite integer enumeration followed by
> exact Gram-inertia and finite-volume tests.  Moreover every survivor obeys
> `n<=min(B-1,D-1)`.

**Proof.**  For each `R<=B`, there are at most

\[
 (2C+1)^{R(R+1)/2+tR^2+DR}                              \tag{10}
\]

choices before the finite flags, if unused facet slots are allowed.  Summing
over `R` is finite.  Exact linear algebra selects the packets satisfying the
marked algebra and return identities.  The returned Gram and its cyclic
products determine the field inside this finite list.  The spanning,
irredundancy, signature, endpoint and cusp gates select another finite subset.
Concretely, for each candidate one:

1. isolates the finitely many real eigenvalues of the integral marked
   operators and constructs every distinguished projector allowed by the
   finite grammar;
2. checks Gram rank and inertia, the acute/Coxeter entry conditions, a strict
   time-oriented interior point, and irredundancy of every listed halfspace;
3. enumerates all relevant facet subsets and applies the full finite-volume
   endpoint and affine-neighbourhood test, rejecting hyperideal rays and
   incomplete parabolic links.

All three are terminating exact algebraic calculations because `R` and the
facet count are bounded.  A mere count of cusps is not substituted for the
last test.  Here “complete” means that the packet's exact-return rule supplies
the whole ordered facet list defining the returned intersection, not a wall
prefix.  Reflection-maximality is neither assumed by these filters nor needed
for the finiteness conclusion; a proposed universal maximal-chamber compiler
must prove that its exact return really preserves maximality.  Equations (3)
and (6) give the dimension bound.  `square`

This theorem is not disguised arithmetic-volume finiteness: its finite input
is an explicit collection of bounded integer matrices.  In AG language,
`B,D,C` can be certified by a bounded cellular presentation, a bounded Hilbert
packet plus Bezout bounds, or a fixed bounded cohomological lattice and total
cycle degree.  Merely saying “surface of dimension two” or “bounded degree per
curve” does not supply them.

## 4. Four-node packets cannot force the global ceiling

The obstruction is exact.  For `m>=1`, let

\[
 G_m=-2I+\operatorname {Adj}(K_{1,m}),                   \tag{11}
\]

where `K_{1,m}` is the star with one centre and `m` leaves.  Its entries are
only `-2,0,1`, the rational Coxeter alphabet for right angles and
`pi/3` angles.  Direct determinant expansion gives

\[
 \chi_{G_m}(lambda)
  =(lambda+2)^{m-1}\bigl((lambda+2)^2-m\bigr).           \tag{12}
\]

Consequently, for every `m>=5`, `G_m` is nonsingular Lorentzian of signature
`(1,m)` and rank `m+1`.

Every principal submatrix on exactly four nodes is nevertheless one of only
two matrices:

- `-2I_4`, if the centre is absent;
- `-2I+Adj(K_{1,3})`, if the centre is present.

The latter has eigenvalues `-2+sqrt(3),-2-sqrt(3),-2,-2`, so both types are
negative definite.  Even the collection of all at-most-four-node isomorphism
types is the fixed seven-element list

```text
edgeless s K_1,  1<=s<=4;    K_{1,j},  1<=j<=3,
```

independently of `m`.  The one-node edgeless graph is not counted again as
`K_{1,0}`.

> **Theorem 3 (four-packet no-go).**  A finite alphabet of exact four-node
> Gram packet types, bounded label field, bounded local Belyi degree, and
> ellipticity of every four-node principal packet do not imply a bound on the
> rank of the global Lorentzian Gram or on the rank of a faithful global
> cohomological packet.

**Proof.**  The family (11) has all the stated local data fixed, while its
global rank is `m+1`.  Theorem 1 forces any faithful global packet for `G_m`
to have rank at least `m+1`.  `square`

This does not manufacture an unbounded family of finite-volume polytopes.
It pinpoints the missing global condition.  Deleting the centre of `G_m`
leaves a negative-definite matrix.  Deleting a leaf leaves `G_{m-1}`: it is
parabolic for `m=5`, but indefinite for every `m>=6`.  Thus the first row is
the sole row not killed by the simplex finite/ideal-vertex deletion gate;
the unbounded tail fails exactly at a global face condition invisible to
four-node packets.  We use `m=5` only as a passed **deletion gate**, not as a
new-polytope or maximality claim.  A noncompact proof must likewise retain
full affine cusp links, not just their four-node shadows.

There is also a literal surface version.  On

\[
 S_r=\operatorname {Bl}_{p_1,\ldots,p_r}\mathbf P^2,
\]

the exceptional curves have intersection matrix `-I_r`, and an ample class
can be chosen with degree one on every exceptional curve.  Every at-most-four
curve intersection packet is therefore a principal submatrix of `-I_4`, but
`rho(S_r)=r+1` is unbounded.  The required polarizations have unbounded global
Hilbert data.  Hence fixed carrier dimension and bounded **individual** marked
degree do not give `B` or `D`; total degree or global Hilbert/cohomology data
must be bounded.

## 5. The next finite falsifiable lemma

The missing statement should not say merely that local packet types are
finite.  The star family proves that insufficient.  The minimal useful AG
statement is the following.

> **Fixed-Hilbert faithful descent lemma (open).**  There exist absolute
> constants `B,D,C,t` and a finite algebraic construction which, from every
> reflection-maximal arithmetic finite-volume Coxeter chamber `P`, produces:
>
> 1. one global integral cohomological/monodromy packet `H_P` and the
>    injection (2), not separate unrelated four-node realizations;
> 2. a certified chain presentation with `rank H_P<=B`;
> 3. distinct facet cycles of total ample degree at most `D`;
> 4. pairing, correspondence and facet coordinates bounded by `C`;
> 5. at most `t` integral interval generators satisfying the hypotheses of
>    the exact `[-2,1]` probe; and
> 6. complete finite-vertex and affine-cusp descent, proving that the returned
>    halfspaces give the same finite-volume chamber.

This lemma is finite and falsifiable for any proposed compiler.  Its checks
are:

```text
overlap cocycle -> one global lattice;
spanning minor -> injectivity of the full trace root module;
chain/Hilbert certificate -> the rank ceiling B;
total intersection -> the facet ceiling D;
integer matrices -> the coordinate ceiling C;
minimal polynomials and isolating intervals -> the 23^t field gate;
full vertex/cusp links -> finite-volume coverage.
```

If the lemma is proved, Theorem 2 gives a global, effective finiteness theorem
and

\[
                       \boxed{n\le\min(B-1,D-1)}.        \tag{13}
\]

If it fails along a sequence of chambers, at least one geometrically
meaningful quantity is forced to grow: global cohomology/cell rank, total
facet degree, matrix height, number of spectral generators, or cusp/overlap
decoration.  That is the rigorous surface/dessin dichotomy.  It is not yet a
global Coxeter theorem because no universal construction is currently known
to satisfy all six clauses.

## 6. Exact replay

- symbolic verifier:
  [`scripts/certify_integral_packet_dimension_bridge.py`](../scripts/certify_integral_packet_dimension_bridge.py);
- frozen archive:
  [`data/enumerations/integral_packet_dimension_bridge.json`](../data/enumerations/integral_packet_dimension_bridge.json);
- focused tests:
  [`tests/test_integral_packet_dimension_bridge.py`](../tests/test_integral_packet_dimension_bridge.py).

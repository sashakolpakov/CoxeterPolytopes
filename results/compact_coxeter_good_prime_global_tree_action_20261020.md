# A conditional global tree action from good white two-faces

20 October 2026. Written research theorem proposed for independent review.
No scientific job, census, optimizer, Lean or build was run. This gives a
global structural restriction, not a new dimension bound.

## Statement

Let P be a compact hyperbolic Coxeter polytope of dimension n>=3, with
reflection group Gamma and a realization over a nonarchimedean local
field F of odd residue characteristic. Work in its actual Lorentz
representation, base-changed to F; no arithmeticity is assumed.

For EVERY two-face of P, assume its n-2 white root normals satisfy:

1. their squared norms are units;
2. all Gram entries are integral and every nonzero one is a unit;
3. their Gram determinant is a unit.

Then Gamma admits an action on a simplicial tree T and a Gamma-equivariant
map from T to its orthogonal Bruhat--Tits building X. If its action on X
is unbounded, its action on T has no global fixed point. In particular,
if Gamma has the fixed-point property for tree actions, its F-representation
is bounded. A nonintegral matrix trace is enough to give the unbounded
case. No fixed-point property is asserted here for general Coxeter groups.

The conclusion is deliberately conditional on all white faces being good
at the SAME place. A good face chosen separately for each prime or each
vertex is not sufficient.

## Inputs and the graph-subcomplex point

Use the self-dual norm / almost-self-dual lattice-chain model of X, as in
Cornut, *On p-adic norms and quadratic extensions, II*, Section3, printed
pages6--8:
https://webusers.imj-prg.fr/~christophe.cornut/papers/npeq2.pdf .
Use also the standard facts that X is a complete CAT(0) space, a bounded
group orbit has a fixed point, and the Coxeter tessellation of H^n is a
contractible polyhedral complex with pointwise cell stabilizers. These
are the external geometric/topological inputs.

The root-projector proof in
`general_coxeter_good_prime_two_face_trees_20261020.md`, Sections1--2,
shows that the fixed set Y of a two-face's white finite group is exactly

    alpha=alpha_W direct-sum alpha_U,

where alpha_W is the standard norm on its unimodular white root lattice
L_W and U has dimension three. Thus Y is a tree or a point.

Here it is essential to check more than abstract tree topology: Y is a
subcomplex of the SAME lattice-chain complex X. In multiplicative norm
notation, the white closed ball for every radius q^t, 0<=t<1/2, is
exactly L_W. Therefore the almost-self-dual lattices associated to a
norm in Y are exactly

    L_W direct-sum L_U,

as L_U runs through its complement's almost-self-dual chain. A vertex
of Y is a vertex of X; an edge of Y is an edge of X. This is a chain
identity, not a guess that an arbitrary geodesic tree in a higher-rank
building lies in its one-skeleton. It works because the white norm is
standard and unimodular at the selected place.

The orthogonal group acts without inversions on this chain complex. The
type of a vertex lattice L includes dim(L*/L); these dimensions strictly
change along a strictly increasing almost-self-dual chain, while
isometries preserve them. Consequently an isometry preserving a simplex
setwise fixes each of its vertices. A finite subgroup has a fixed point
by CAT(0), and hence a fixed vertex, by taking the vertices of the unique
open simplex containing that point.

## Equivariant map of the actual tessellation

Let K be the Coxeter tessellation of H^n by copies of P. Its two-skeleton
is simply connected: attaching cells of dimension at least three does
not change the fundamental group, and K is contractible.

For every Gamma-orbit of vertices of K choose one vertex of X fixed by
its finite stabilizer, and extend the choice equivariantly. This is
well-defined. Map each edge to the unique CAT(0) geodesic between its
endpoint images. Its finite edge stabilizer fixes that path pointwise,
so these choices too are equivariant.

For an actual two-face F, both endpoint stabilizers of each boundary
edge contain its white group H_F. Thus every boundary vertex image is
in Y_F=Fix(H_F). The unique geodesic on every boundary edge is also in
Y_F, by convexity. Because Y_F is the actual tree subcomplex just
identified, this is a finite edge path in X's one-skeleton.

Every two-face boundary therefore extends continuously across its disk
into Y_F, with no image outside X's one-skeleton. Choose these extensions
on orbit representatives and extend equivariantly. This is legitimate
because cell stabilizers act pointwise on the source face and fix Y_F
pointwise. We obtain a continuous Gamma-equivariant map

    f: K^(2) --> G subset X^(1),

where G is the connected invariant graph formed by all the image edge
paths (isolated image points are included). Each disk extension can be
taken in the finite subtree spanned by its boundary, so its image lies
in that same graph G. If necessary subdivide source edges along their
finite image edge paths. No black reflection has been restricted to a
white orthogonal complement in this construction.

## Lift and the nontrivial-action conclusion

Let pi:T-->G be the universal cover of the connected graph G. It is a
tree. Since K^(2) is simply connected, f has a lift f_tilde. Fix a
basepoint x0. For each gamma in Gamma, choose the unique lift gamma_tilde
of its graph automorphism for which

    gamma_tilde(f_tilde(x0))=f_tilde(gamma*x0).

Uniqueness of lifts on the connected domain K^(2) gives

    gamma_tilde composed with f_tilde = f_tilde composed with gamma.

Applying this identity twice proves the group law for the chosen lifts:
two lifts of the same graph automorphism agreeing at one point are equal.
Thus Gamma genuinely acts on T; this is not merely an unspecified
extension of Gamma by the deck group. The covering map pi, followed by
G subset X, is equivariant. A global fixed point of T would map to a
global fixed point of X. Therefore an unbounded F-image gives a tree
action without a global fixed point. Barycentric subdivision removes
any possible tree-edge inversions.

## What this does and does not settle

This upgrades the local chordal restrictions to an actual global action
under uniform good-white hypotheses. It is a possible route to studying
reflection-group splittings or decomposition into pieces; such a
decomposition and its dimension consequences have NOT been proved here.
General compact Coxeter groups are not simply assumed to have the
fixed-point property for trees. The explicit bad Lambert quadrilateral
in the companion note remains a warning that short-face counting alone
does not remove nonintegrality.

Primes above2, bad white determinants, nonunit bonds, arbitrary ideal
vertices, and other real embeddings remain separate issues. In
particular this theorem does not prove arithmeticity or a uniform
compact/finite-volume dimension ceiling. Those would require a further
argument, using an actual obstruction to the resulting group action or
a controlled geometric decomposition.

# A movable Coxeter cone and complete Kähler wall shells on the golden OG10 carrier

> **Subsequent decision.**  The complete norm-two chamber is now proved
> infinite-volume and noncompact by an infinite stabilizer of its complete H9
> fixed slice.  Finite-sidedness remains undecided.  See
> [`og10_golden_h9_infinite_stabilizer.md`](og10_golden_h9_infinite_stabilizer.md).

## Result

The golden OG10 construction now reaches the geometric chamber level, not
only the abstract divisor-lattice level.  Let

\[
H=\begin{pmatrix}-2\phi&1\\1&2\end{pmatrix}
  \oplus E_8\oplus\langle2\rangle
\]

over `O=Z[phi]`, and let `N` be its rank-22 trace lattice with pairing

\[
B(x,y)=\operatorname{Tr}_{\mathbf Q(\sqrt5)/\mathbf Q}
\left(-\frac{\phi}{\sqrt5}H(x,y)\right).
\]

The carrier theorem realizes `N` as the full Neron--Severi lattice of a
projective OG10-type manifold and realizes the `phi`-eigenspace as a
hyperbolic `10`-space `V_phi`.

Fix the positive-cone component `C_N^+` compatible with the marked-moduli
component and write

\[
V_\phi^+=V_\phi\cap C_N^+.
\]

Below, a restricted chamber always means a connected component of the
mirror complement in `V_phi^+`.

The new theorem is:

> **Complete-wall and chamber-selection theorem.**  In the explicit
> primitive embedding
> \(N\hookrightarrow\Lambda_{OG10}\), no primitive class of `N` has ambient
> divisibility `3`.  Consequently the complete numerical OG10 wall
> arrangement in `N` consists exactly of the classes of square `-2` and
> `-4`, all of ambient divisibility one.  Its restriction to `V_phi^+` is
> exactly the union of the mirrors of field vectors having one of the three
> norms
> \[
> H(r,r)\in\{2,\ 4,\ 6-2\phi\}.
> \]
> Every connected component of the complement of this restricted arrangement
> is the intersection with
> `V_phi^+` of the Kähler cone of a same-period projective marked OG10 model.

There is a sharper Coxeter statement at the movable level:

> The open movable cone is cut only by the square `-2` stably
> prime-exceptional
> divisors.  After a change of marking, its intersection with `V_phi^+` can be
> made equal to any prescribed chamber of the norm-two field reflection
> arrangement.  It is therefore an exact, possibly infinite-sided,
> hyperbolic Coxeter chamber.  Every one of its generating reflections is
> the restriction of a product of two commuting divisorial OG10 monodromy
> reflections.

The last assertion permits a Hodge-monodromy re-marking.  It does not say
that every chamber is obtained by the canonical graph-pullback marking from
one fixed model.

This is a genuine algebro-geometric Coxeter-chamber construction theorem in
hyperbolic dimension `10`.  It is not yet a finite-volume or compact
Coxeter-polytope theorem: the norm-two chamber may have infinitely many
facets or infinite volume.  Its Kähler subdivisions require two additional
norm orbits and need not themselves have only Coxeter dihedral angles.

## Why the divisibility-three wall types vanish

The explicit carrier gluing uses the ambient integral basis consisting of
the rank-18 unimodular core followed by the six columns of the glued
rank-six basis.  In this same basis, the primitive embedding is literally

\[
N=\langle e_1,\ldots,e_{22}\rangle
  \subset\Lambda_{OG10}=\langle e_1,\ldots,e_{24}\rangle.
\]

The ambient Gram matrix has determinant `-3` and Smith nonunit `[3]`.
For a primitive ambient vector, its divisibility is the order of the
corresponding class in the discriminant group, so it is therefore either
`1` or `3`.
Exact row reduction modulo `3` gives

\[
\operatorname{Rad}(\Lambda_{OG10}/3\Lambda_{OG10})
=\mathbf F_3(0^{19},2,1,0,0,1).
\]

The last coordinate of this generator is nonzero, whereas both final
coordinates of every class coming from `N/3N` vanish.  Thus

\[
\operatorname{Rad}(\Lambda_{OG10}/3\Lambda_{OG10})
\cap (N/3N)=0.
\]

Equivalently, the `24 x 22` pairing map from `N` to the dual of the ambient
lattice has full column rank `22` modulo `3`.

If `D in N` has ambient divisibility three, all its ambient pairings vanish
modulo three.  The preceding injectivity forces `D=3D_0` for some `D_0 in N`.
In particular `D` is not primitive.  There is also a square-level
contradiction: `N` is even, so

\[
D^2=9D_0^2\in18\mathbf Z.
\]

It can therefore have neither square `-6` nor square `-24`.

The same radical gives a useful scope theorem.  If `tau` extended to an
integral self-adjoint golden endomorphism of the full Beauville lattice,
self-adjointness would preserve its one-dimensional radical modulo `3`.
But

\[
t^2-t-1
\]

has no root in `F_3` (its values at `0,1,2` are `2,2,1`), so it cannot act on
a nonzero one-dimensional invariant space while satisfying the golden
polynomial.  Therefore **no integral self-adjoint extension of `tau` to the
full OG10 Beauville lattice exists**.  The enhancement is intrinsically an
endomorphism of the divisor lattice.  This does not exclude an arbitrary
non-self-adjoint integral extension, which would have no compatible
Beauville-form interpretation.

Mongardi--Onorati's Proposition 5.4 gives the full list of primitive OG10
wall-divisor types

\[
(-2,1),\quad(-4,1),\quad(-6,3),\quad(-24,3),
\]

where the second entry is divisibility in the full Beauville lattice, not in
`N`.  The calculation above removes the last two types for this particular
Neron--Severi lattice.  Conversely, every class of square `-2` or `-4` is
primitive because `N` is even, and its ambient divisibility cannot be three.
It is therefore a wall divisor.  This proves completeness of the numerical
wall list on the carrier.

## The finite field-norm shell theorem

Write

\[
q=H(D,D)=a+b\phi\in\mathbf Z[\phi].
\]

The trace scalar has

\[
\operatorname{Tr}\left(-\frac{\phi}{\sqrt5}\right)=-1,
\qquad
\operatorname{Tr}\left(-\frac{\phi^2}{\sqrt5}\right)=-1,
\]

so

\[
B(D,D)=-(a+b).
\]

The field lattice `H` is even: its diagonal coefficients lie in `2O` and its
off-diagonal coefficients lie in `O`.  Hence both `a` and `b` are even.
The conjugate form is positive definite, so `q'` is positive for every
nonzero vector.  The restricted hyperplane meets the positive cone of
`V_phi` exactly when `q` is positive at the identity embedding.  Thus the
possible wall norms are precisely the totally positive even elements with
`a+b=2` or `a+b=4`.

For `s=a+b`, total positivity is equivalent to

\[
-s\phi<b<s/\phi.
\]

There are consequently only seven raw values:

| Beauville square | Field value `q` | Field norm | Unit-square orbit |
|---:|---:|---:|---:|
| `-2` | `4-2phi` | `4` | `2` |
| `-2` | `2` | `4` | `2` |
| `-4` | `10-6phi` | `4` | `2` |
| `-4` | `8-4phi` | `16` | `4` |
| `-4` | `6-2phi` | `20` | `6-2phi` |
| `-4` | `4` | `16` | `4` |
| `-4` | `2+2phi` | `4` | `2` |

Multiplication of a vector by `phi^k` multiplies its field norm by
`phi^(2k)` and does not change its hyperplane in `V_phi`.  The seven values
therefore form exactly three mirror orbits:

\[
\begin{aligned}
\mathcal O_2&=\{10-6\phi,\ 4-2\phi,\ 2,\ 2+2\phi\},\\
\mathcal O_4&=\{8-4\phi,\ 4\},\\
\mathcal O_{6-2\phi}&=\{6-2\phi\}.
\end{aligned}
\]

This proves that a complete restricted-wall search requires exactly three
field norm equations, represented by

\[
2,\qquad4,\qquad6-2\phi.
\]

All three families are nonempty.  In the field basis ordered as the coupled
rank-two block, eight `E8` coordinates, and the final positive line, explicit
vectors are

\[
\begin{array}{c|c}
r&(H(r,r))\\ \hline
(0,0,0^8,1)&2\\
(0,0,1,0^7,1)&4\\
(1,1,1,0^7,0)&6-2\phi.
\end{array}
\]

The latter two mirrors cannot secretly be norm-two mirrors.  Equality of two
field hyperplanes would make their normals field-proportional.  Relative to
a norm-two normal this would require `2` or `3-phi` to be a square in
`Q(sqrt(5))`.  The first is not a square in that quadratic field; the second
has field norm `5`, which cannot be the norm of a square.  The ratio
`(3-phi)/2` between the last two families has norm `5/4`, so those families
are distinct from each other as well.

The two extra Kähler-wall families are genuinely nonreflective over `O`.
Pairing either displayed vector with an adjacent simple `E8` root gives
`-1`.  The corresponding field-reflection coefficients include

\[
\frac{2}{4}=\frac12,
\qquad
\frac{2}{6-2\phi}=\frac{2+\phi}{5},
\]

which are not in `Z[phi]`.  This is why the full Kähler subdivision is not
automatically Coxeter, while the norm-two movable arrangement is.

Thus the existing norm-two Vinberg prefix is provably a strict subarrangement
of the complete restricted **Kähler** wall arrangement.

## Every golden reflection is a divisorial monodromy product

The norm-two arrangement has a stronger geometric interpretation than mere
equality of hyperplanes.  Let `r` be a field root with `H(r,r)=2`, let
`D=r` be its underlying divisor, and put

\[
E=\tau^{-1}D=(\tau-1)D.
\]

The trace calculations

\[
D^2=-2,\qquad (\tau D)^2=-4,\qquad(D,\tau D)=-2
\]

give

\[
E^2=-2,\qquad(D,E)=0.
\]

Both `D` and `E` are therefore divisorial, stably prime-exceptional OG10
classes.  Let `R_D` and `R_E` be their Markman reflections.  They commute,
their product is `-1` on the rank-two lattice

\[
\mathbf ZD\oplus\mathbf ZE=\mathcal O r,
\]

and it is the identity on the common orthogonal complement.  The latter is
exactly the field orthogonal `H(-,r)=0`: if `H(x,r)=a+b phi`, pairing with
`D` gives `-(a+b)` and pairing with `E` gives `-b`, so both pairings vanish
if and only if the field pairing vanishes.

It follows on `N` that

\[
(R_D R_E)|_N=s_r
\]

where `s_r` is the integral field reflection in `r`.  In particular

\[
s_r\in \operatorname{res}_N(W_{Exc})
       \cap\operatorname{Cent}_{O(N)}(\tau),
\]

and its restriction to `V_phi` is the desired hyperbolic Coxeter reflection.
Thus every golden facet reflection is the restriction of a product of two
commuting divisorial monodromy reflections of the algebraic variety.  No
extension of `tau` from `N` to the full Beauville lattice is asserted.

## The movable cone is the canonical Coxeter object

The distinction between divisorial and small-contraction walls isolates the
right algebraic cone.  Mongardi--Onorati Proposition 3.1 says that the
stably prime-exceptional OG10 classes have types

\[
(-2,1)\quad\hbox{and}\quad(-6,3).
\]

Their Proposition 5.4 adds the two small-contraction Kähler-wall types
`(-4,1)` and `(-24,3)`.  Since the carrier contains no primitive
divisibility-three classes, only square `-2` classes bound its fundamental
exceptional, hence movable, cone.  The square `-4` classes subdivide that
cone into Kähler chambers but do not bound the movable cone.

More precisely, Markman's Definition 5.2, Propositions 5.6 and 6.10, and
Theorems 6.17--6.18 identify the fundamental exceptional chamber `FE` as a
component of the positive cone cut by all stably prime-exceptional
hyperplanes and as a fundamental chamber for `W_Exc`.  Definition 6.21 and
Lemma 6.22 identify

\[
\operatorname{Mov}^0(X)=FE_X\cap NS(X)_{\mathbf R}
\]

and make the closed movable cone a `W_Exc` fundamental domain.  Here the
Picard rank is maximal, `rho=22=h^(1,1)`, so the open movable cone is exactly
the relevant fundamental exceptional chamber.

On `V_phi`, the visible square `-2` shells are `2` and `4-2phi`; multiplying
the second by `phi^2` gives the first without changing the mirror.  Hence

\[
FE_X\cap V_\phi^+
\]

is, whenever nonempty, exactly a connected component of the complement of
the complete norm-two field reflection arrangement.

Nonemptiness and chamber selection are automatic after marking.  Given any
norm-two field chamber `P`, choose a regular point `x in P` and let `C` be
the ambient exceptional chamber containing it.  Hyperplane signs and
convexity give `C\cap V_\phi=P`.  The exceptional reflection group acts
simply transitively on exceptional chambers, so `C=w(FE_X)` for some
`w in W_Exc`.  Postcomposing the marking by `w` yields

\[
\eta'(\operatorname{Mov}^0(X))\cap V_\phi^+=P.
\]

This re-marking conjugates the golden endomorphism on the underlying divisor
lattice; `w` need not commute with the previously marked `tau`.  The target
carrier datum and its fixed eigenspace remain the prescribed ones.

Thus a projective OG10 variety with the designed divisor lattice
automatically produces a hyperbolic Coxeter chamber: its marked eigen-movable
cone.  The field reflections are integral because `H(r,r)=2`, and the
previous section shows that they are restrictions of genuine divisorial
monodromy products.

This statement concerns the movable/fundamental-exceptional cone, not the raw
birational Kähler union.  Markman's inclusions

\[
\mathcal{BK}_X\subset FE_X\subset\overline{\mathcal{BK}}_X
\]

allow the square `-4` flop hyperplanes to lie inside `FE_X`; they disappear
only after passing to the movable closure/interior.

## Why an actual Kähler slice exists

Fix the carrier period and let `C_N^+` be its compatible positive-cone
component.  Let `Delta_N` be the set of all primitive classes in `N` having
one of the four numerical OG10 wall types.  For this carrier, the preceding
modulo-three theorem reduces `Delta_N` to the square `-2` and `-4` classes.

Markman's Kähler-type chamber theory and the OG10 wall classification imply:

> Every component
> \[
> Q\subset C_N^+\setminus\bigcup_{\delta\in\Delta_N}\delta^\perp
> \]
> is the marked Kähler cone of a marked OG10 manifold in the fixed period
> fiber.

One way to assemble the published inputs is to choose a very general point
of `Q`.  Markman's Lemma 5.11 places it in a Kähler-type chamber.  By
Mongardi--Onorati Proposition 5.4 and Theorem 5.5, that Kähler cone is a
component of exactly the same numerical wall arrangement; because it meets
`Q`, the two components coincide.  Here parallel transport preserves the
Beauville square, ambient divisibility, and hence the four numerical OG10
wall types.  Markman's Theorem 5.16 realizes the
Kähler-type chamber by a marked pair in the fixed period fiber.  Its marked
Neron--Severi lattice is again `N`, and the positive integral class in `N`
makes the model projective.

Now let `U` be a connected component of the complement of the arrangement
restricted to `V_phi^+`.  A
point of `U` lies in a unique ambient chamber `Q`.  Ambient chambers are
convex, and their defining signs restrict to the signs defining `U`, so

\[
U=Q\cap V_\phi^+.
\]

No nonzero rational divisor hyperplane contains all of `V_phi`: otherwise
its normal would be a rational vector in the conjugate irrational
eigenspace, which is impossible for the integral golden action.  The marked
model realizing `Q` therefore satisfies

\[
\eta(\mathcal K_Y)\cap V_\phi^+=U.
\]

This closes the earlier chamber-selection gap.  An arbitrary chamber can
require a Hodge-monodromy re-marking; canonical pullbacks from birational
models of one fixed marked variety remain inside Markman's fundamental
exceptional chamber.

The Vinberg control vector used in the bounded search is itself on the ten
height-zero walls, so it is not an interior Kähler class.  Local finiteness
allows any adjacent full chamber to be selected; the control vector then
lies in its nef closure.

## What remains open for the finer Kähler subdivision

For the complete norm-two movable Coxeter chamber, finite volume and
compactness are subsequently decided negatively by the infinite-stabilizer
theorem.  If one separately wants the finer Kähler subdivision, the remaining
finite-alphabet arithmetic tasks are:

1. enumerate simple mirrors for the three norm equations `2`, `4`, and
   `6-2phi`;
2. choose a component and determine its actual facets;
3. verify all mixed-shell dihedral angles rather than assuming they are
   Coxeter angles;
4. determine the volume of the selected finer chamber; and
5. certify compactness independently before making any compact-polytope claim.

A numerical wall class is a facet only when its hyperplane is adjacent to
the selected chamber with codimension-one intersection.  The theorem does
not make every wall divisor a facet, nor does it imply finite generation or
finite volume for the finer Kähler subdivision.

## Sources and reproducibility

The geometric inputs are Eyal Markman's
[Definition 5.2, Propositions 5.6 and 6.10, Lemma 5.11, Theorems 5.16,
6.17, and 6.18, and Definition 6.21/Lemma 6.22](https://arxiv.org/abs/1101.4606)
on Kähler-type, exceptional, and movable chambers, and Giovanni
Mongardi--Claudio Onorati's
[Propositions 3.1 and 5.4 and Theorem 5.5](https://arxiv.org/abs/2010.12511)
on OG10 prime-exceptional classes, wall divisors, and Kähler classes.
Onorati's computation of the OG10
monodromy group is available at
[arXiv:2002.04129](https://arxiv.org/abs/2002.04129).  The published
correction to the OG10 birational-geometry paper does not alter the wall
classification used here.

Run the exact certificate and its frozen-artifact test from the repository
root:

```sh
PYTHONPATH=.:src python3 scripts/certify_og10_golden_wall_shells.py
PYTHONPATH=.:src pytest -q tests/test_og10_golden_wall_shells.py
```

The archived output is
`data/candidates/og10_golden_complete_wall_shells.json`.  It certifies the
ambient modulo-three radical, the full-rank pairing map from `N`, the trace
formula, the exhaustive seven-value shell list, the three unit-square
orbits, and explicit representatives of all three families.  The
Markman/Mongardi--Onorati chamber theorem is an external mathematical input,
not replaced by computation.

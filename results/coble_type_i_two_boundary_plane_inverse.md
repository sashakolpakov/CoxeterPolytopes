# TOOT, qualified: raw `H11/26` forces a line--quintic Coble model

**Status (2026-08-09): exact `POLYTOPE -> SURFACE -> same marked
POLYTOPE` theorem for the two-boundary Type-I chamber, after an explicit
hostile-audit repair.**  The original root-seed handoff constructed only 16
of 26 negative curves and was insufficient as written.  The repaired primary
certificate now constructs and absolutely factors the ten missing cusp-lift
curves before invoking hidden-wall exhaustion.  The abstract
finite-volume noncompact `H11/26` chamber and the Coble surface isomorphism
class were already frozen.  The new project result is a raw-wall derivation
of a horizontal plane model, including its complete determinantal existence
proof and literal same-chamber return.  Worldwide priority for this particular
line--quintic birational model is open.

The executable
[`certify_coble_type_i_two_boundary_plane_inverse.py`](../scripts/certify_coble_type_i_two_boundary_plane_inverse.py),
the frozen
[`JSON`](../data/enumerations/coble_type_i_two_boundary_plane_inverse.json),
and the
[`tests`](../tests/test_coble_type_i_two_boundary_plane_inverse.py)
certify every count and identity below.  The full replay peaks below `0.12
GiB`.  It invokes the helper-free
[`raw hostile replay`](coble_type_i_two_boundary_plane_inverse_raw_hostile_audit.md),
which performs fresh exact CM and full-cone double descriptions and peaks at
`0.086 GiB`.

## 1. Exact inverse theorem

> **Theorem.**  Start with only the odd-unimodular lattice Gram and the 26
> primitive oriented walls of the marked `H11/26` chamber.  Its exhaustive
> numerical contraction compiler has 52 sequential `P2` endpoints and nine
> endpoint-degree signatures.  Exactly four endpoints have the intrinsic
> profile
>
> ```text
> boundary degrees:              1,5
> actual roots after contraction: 7 of degree 0, 3 of degree 1
> contracted walls:              7 square -2, 4 square -1.
> ```
>
> All four recover the same proximity/root program.  Its compiler-generated
> quintic jet matrix has rank-drop locus `a^2-b^2=0`.  The component
> `a=-b` has nonreduced boundary; the component `a=b!=0` constructs over
> `Q` a smooth rational line boundary and the resolution of an absolutely
> irreducible rational quintic boundary.  Blowing up the derived cluster gives
> a strong two-boundary Coble surface.  Its projectivized full nef cone, with
> the induced integral marking and orientation, is literally the input
> 26-wall finite-volume noncompact polytope in `H11`.

No surface name, Kondō type, equation, known contraction, or stored output
Gram participates in selecting the four endpoints.  Surface theory enters
only after the raw program has constructed all twelve Coble--Mukai root
seeds.

## 2. What the raw walls reconstruct

For each selected branch the compiler obtains a determinant-one basis

\[
                   [H,E_1,\ldots,E_{11}],\qquad
                   \operatorname{Gram}=\operatorname{diag}(1,-1^{11}),
\]

and

\[
                         K=-3H+\sum_{i=1}^{11}E_i.
\]

The degree-zero roots are seven proximity differences.  Their forest is the
disjoint union of four chains of lengths

\[
                              3,\ 2,\ 3,\ 3.             \tag{2.1}
\]

The two first chains lie on the degree-one boundary.  The other two resolve
the singularities of the quintic.  The three positive-degree roots are
lines through exactly three assigned proper or infinitely-near centers.

After a projective normalization, take

```text
A=(0:0:1),  B=(1:0:1),  C=(0:1:1),  D=(-1:2:1).
```

The boundary line is `y=0`.  The root lines are

\[
                 AC:x=0,\qquad AD:2x+y=0,
                 \qquad BCD:x+y-z=0.                    \tag{2.2}
\]

The line-boundary chains have lengths three at `A` and two at `B`.
At `C` and `D`, the first tangent is respectively `AC` and `AD`; the
remaining second-near slopes are the only parameters `a,b`.

This reconstruction is numerical and invariant, although the four raw
endpoint index sets and their exceptional orderings differ.  In every branch
the two deepest line-contact exceptionals have boundary profile `(1,1)` and
give the two crossing `(-1)` curves.  The deepest exceptionals over `C,D`
have profile `(0,2)` and are two already visible cusp lifts.

This assertion is a coefficientwise compiler check, not a comparison to a
hard-coded plane template.  For example, one raw branch first yields the
four proper-to-deep chains

```text
[1,0], [4,3,2], [7,6,5], [10,9,8].
```

The degree-one boundary support intrinsically singles out the length-three
and length-two chains as `A,B`; the quintic multiplicity-two support singles
out the other two as `C,D`.  The positive-root supports are then

```text
[1,7,10], [4,6,7], [4,9,10],
```

which canonically distinguish `BCD` and the two `A`-tangent lines.  The same
algorithm runs independently on each of the four endpoints.  After its
derived relabeling, the two boundary classes, all ten actual-root classes,
the two crossing exceptionals, and the two visible cusp exceptionals agree
coefficient by coefficient.  The four resulting canonical profiles are
identical.

## 3. The complete symbolic existence calculation

A plane quintic has 21 coefficients.  The proximity multiplicities impose

```text
A line-contact chain:       1+1+1 = 3 rows
B line-contact chain:       1+1   = 2 rows
C double chain:             3+3+3 = 9 rows
D double chain:             3+3+3 = 9 rows
                                      -------
total                                  23 rows.
```

Let `M(a,b)` be this canonical `23 x 21` jet matrix.  Deleting rows 11 and
20 gives the exact maximal minor

\[
                 \det M_{\widehat{11},\widehat{20}}
                    =-64(a-b)(a+b).                    \tag{3.1}
\]

Direct symbolic row reduction gives

\[
 \operatorname{rank}M(a,a)=20,\qquad
 \operatorname{rank}M(a,-a)=20,
\]

while, for example, `rank M(1,2)=21`.  Hence (3.1) is not just one
necessary determinant:

\[
               \boxed{\operatorname{rank}M(a,b)\le20
                      \iff a^2-b^2=0.}                 \tag{3.2}
\]

On `a=1,b=-1` the one-dimensional kernel is

\[
 y\bigl(x^2-xy+2xz-y^2+2yz-z^2\bigr)^2.               \tag{3.3}
\]

It contains the line boundary and is nonreduced, so the raw boundary gates
canonically discard this component.  No chosen surface coordinate witness
is being smuggled into membership: the compiler generates the jet matrix,
the factorization (3.2), both kernel vectors, and the reducedness test from
the raw divisor/proximity program.

The other factor is horizontal.  At `a=b=1`, its homogeneous quintic is

\[
\begin{aligned}
Q={}&16x^5+39x^4y-32x^4z+34x^3y^2-52x^3yz+16x^3z^2\\
 &+9x^2y^3-16x^2y^2z+6x^2yz^2-2xy^4+8xy^3z\\
 &-10xy^2z^2+4xyz^3-y^5+4y^4z-6y^3z^2+4y^2z^3-yz^4.
                                                               \tag{3.4}
\end{aligned}
\]

In particular this branch has a component dominating an open of
`Spec Z[1/2]`; it is not characteristic-three or characteristic-five
geometry.

## 4. Geometry of the quintic and boundary

Exact `Singular` calculations over `Q` prove that `Q` is absolutely
irreducible and that its saturated projective Jacobian has exactly the two
points `C,D` in its radical.  Its singular scheme has length ten, five at
each point.

The local blowup calculation is more informative.  At both `C,D` the
multiplicity sequence is `(2,2,2)`.  At `C`, the three tangent cones are

\[
                 -v^2,\qquad -(u-w)^2,
                 \qquad7u^2+2uq-q^2,
\]

whose terminal discriminant is `32`.  At `D` they are

\[
                 -2v^2,\qquad -2(u-w)^2,
                 \qquad-2(3u^2+2uq+q^2),
\]

with nonzero terminal discriminant.  Thus both singularities are of type
`A5` in characteristic zero and are resolved by the three declared blowups.
Each has delta invariant three.  Since a quintic has arithmetic genus six,
its strict transform is smooth, geometrically integral and rational.

The restriction to the line boundary is

\[
                         Q(x,0,1)=16x^3(x-1)^2,          \tag{4.1}
\]

which proves the required contact orders three and two.  After the five
declared blowups, the two strict boundary transforms are disjoint and smooth.
Their classes are

\[
\begin{aligned}
B_5&=5H-E_1-\cdots-E_5-2E_6-\cdots-2E_{11},\\
B_1&=H-E_1-\cdots-E_5.
\end{aligned}
\]

Consequently

\[
 B_5^2=B_1^2=-4,\qquad B_5.B_1=0,
 \qquad B_5+B_1=6H-2\sum E_i=-2K.                      \tag{4.2}
\]

The usual negative-intersection argument is now sufficient and avoids a
fragile expected-dimension assertion.  Every divisor in `|-2K|` must contain
both `(-4)` boundaries, whose sum already has the full class; hence

\[
                   |-2K|=\{B_5+B_1\}.
\]

If `D` lay in `|-K|`, then `2D=B_5+B_1`, contradicting reducedness.  Thus

\[
                       |-K|=\varnothing,
                       \qquad|-2K|=\{B_5+B_1\},          \tag{4.3}
\]

and the constructed rational surface is strong Coble.

## 5. Root seed and exact return

The seven proximity exceptionals and the three lines (2.2) are ten smooth
actual `(-2)` curves.  The two deepest line-contact exceptionals are smooth
`(-1)` curves meeting the two boundary components once each; they produce
the other two effective irreducible Coble--Mukai roots.  The two deepest
quintic exceptionals are also explicit `(-1)` curves, each meeting one
boundary twice.  Thus the plane model explicitly constructs four `(-1)`
walls.  The determinant-one marking identifies these classes and all ten
actual-root classes coefficientwise with the raw wall stock on every one of
the four endpoints.  This is a 16-wall root seed, not yet a same-nef-cone
proof.

The hostile audit found the missing implication: enumerating an integral
class `f-B_i/2` does not itself make that class effective.  At `a=b=1`, the
ten remaining positive-degree `(-1)` classes have degree multiset

```text
2,2,2,2,3,3,4,4,5,5.
```

For each class the repaired compiler builds its four exact cluster-jet
blocks.  Every linear system has projective dimension zero.  Literal chart
replay verifies the exact proper and infinitely-near multiplicities, and
`Singular` gives one absolute component of exponent one for all ten curves.
Since every strict-transform class has `E^2=K.E=-1`, its arithmetic genus is
zero.  Thus all ten strict transforms are smooth geometrically integral
rational `(-1)` curves.  The plane model now constructs all 26 raw wall
curves before hidden-wall exhaustion is invoked.

The exact-return gates are now:

1. From the raw vectors, `C` for each of the ten `(-2)` walls and
   `2E+B_1/2+B_2/2` for each crossing reconstruct a rank-ten, twelve-root CM
   diagram.  Its Gram is literally the pinned Type-I diagram Gram.  The
   literal Vinberg audit has `0` bad finite-volume edges and `0` invalid
   affine neighborhoods, with `14` finite and `9` ideal vertices, so its
   reflection group is finite index.
2. The twelve single-boundary walls group by the raw formulas
   `f=E+B_i/2` into exactly six distinct isotropic classes, namely CM cusps
   `2,4,5,6,8,9`; every class has both boundary partners and each wall is
   recovered exactly as `E=f-B_i/2`.
3. The independent interpolation constructs the other ten lifts
   coefficientwise.  This supplies the formerly missing
   `listed class => actual curve` implication.
4. Kondō--Vinberg root exhaustion excludes hidden CM roots.  Any hidden
   same-boundary `(-1)` curve gives a primitive nef isotropic CM ray; the
   complete cusp calculation identifies its class with one of the twelve
   already constructed lifts.  The negative-curve trichotomy leaves no
   other nef wall.

Therefore

\[
                         \boxed{\mathbf P(\operatorname{Nef}S)=P} \tag{5.1}
\]

with the original marking and orientation.  The independent exact cone and
literal Vinberg certificate remains

```text
dimension / facets:           H11 / 26
finite / ideal / hyperideal:  508 / 42 / 0
volume:                       finite
compactness:                  noncompact.
```

Finite volume is not inferred from the finite wall list or from the mere
presence of ideal rays.  In particular, `0` hyperideal rays is used only in
the separate full-cone volume gate; it is not used as a surrogate for the
same-`P` transfer.  The latter requires the twelve-root finite-index gate,
the six complete cusp pairs, the explicit `16+10` geometric construction,
the hidden-wall theorems, and literal equality of all 26 marked wall vectors.
The executable no longer accepts the theorem by searching for phrases in an
upstream JSON.

## 6. Significance for the inverse programme

This supplies the second genuinely horizontal low-row calibration demanded
by the corrected inverse target:

| raw input | generated incidence equation | surviving support |
|---|---|---|
| `H10/30` | `c^2-1=0`, with the collision branch removed | horizontal rational point `c=-1` |
| `H11/26` | `a^2-b^2=0`, with the nonreduced branch removed | horizontal component `a=b!=0` |
| `H14/60` | polar/chord equations have integer content `3` after saturation | vertical at characteristic `3` |
| `H17/48` | complete-projectivity equations have integer content `3` after saturation | vertical at characteristic `3` |

This is strong evidence for a P-derived lattice-determinantal grammar, but it
does **not** by itself prove the arbitrary-input realization theorem for that
grammar.  The separate candidate-class report records precisely what remains.

The result is wholly cusped: a Coble full-nef chamber has primitive nef
isotropic rays.  It supplies no compact polytope and no compact
nonexistence theorem.

## 7. Lineage and replay

Kondō's [classification paper](https://arxiv.org/abs/2102.00591), Example
4.1 and Proposition 4.2, already constructs and classifies the two-boundary
Type-I Coble surface through a `P1 x P1` double-cover presentation; Section 7
recovers those divisors in the uniqueness proof.  Hence neither the surface
nor its isomorphism class is new.  The project-new claim is deliberately
narrower: the raw Coxeter chamber forces the plane cluster, the line--quintic
rank locus, an explicit horizontal root seed, and the exact return (5.1).

Replay with the memory guard:

```bash
PYTHONPATH=src:. python3 scripts/run_with_rss_limit.py --memory-gib 1 --timeout 120 -- \
  python3 scripts/certify_coble_type_i_two_boundary_plane_inverse.py \
  --check-output data/enumerations/coble_type_i_two_boundary_plane_inverse.json

PYTHONPATH=src:. python3 scripts/run_with_rss_limit.py --memory-gib 1 --timeout 120 -- \
  python3 -m pytest -q tests/test_coble_type_i_two_boundary_plane_inverse.py
```

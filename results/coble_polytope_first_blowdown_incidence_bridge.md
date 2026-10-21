# TOOT: raw Coble polytopes compile plane and grid incidence problems

**Status (2026-08-09): two exact single-chamber `POLYTOPE -> SURFACE`
blowdown compilers, including a new raw-MII compiler; no non-singleton `L4`
essential-surjectivity theorem.**

The principal new result is that the Type-MII `H17/48` chamber can be fed to
an inverse algorithm as only an even lattice and forty-eight mirrors.  The
algorithm discovers its eight Coble boundaries uniquely, reconstructs
`Pic(S)=I_(1,17)`, contracts sixteen recovered exceptional curves to
`P1 x P1`, and turns the other twenty-four walls into all twenty-four
permutation graphs on a `4 x 4` grid.  The incidence scheme is nonempty if
and only if the characteristic is three.  The resulting grid blow-up is
proved directly to be a strong Coble surface; its existence is not imported
from a K3 cover or from Kondō's nine-row catalogue.

The raw Type-VI `H14/60` chamber has an equally concrete plane output.  It
has exactly five boundary-preserving blowdowns to `P2`, one for each choice
of boundary as a conic.  The other four boundaries become tangent lines.
Two of the six residual line-root incidences already give the incompatible
characteristic-zero equations

```text
lambda = -1,    lambda = 2,
```

so they meet precisely in characteristic three.  This explains the
characteristic from the polytope-compiled plane geometry rather than merely
reading it from the classification table.

The exact calculation is in
[`certify_coble_polytope_first_blowdown_incidence_bridge.py`](../scripts/certify_coble_polytope_first_blowdown_incidence_bridge.py),
with frozen
[`JSON`](../data/enumerations/coble_polytope_first_blowdown_incidence_bridge.json)
and
[`tests`](../tests/test_coble_polytope_first_blowdown_incidence_bridge.py).

## 1. Logic ledger: what is new and what is still imported

The three logical layers must not be merged.

| layer | Type MII | Type VI |
|---|---|---|
| `P`-only compiler | unique eight-frame, Picard de-scaling, `K4,4`, grid, 24 matchings | independent unique five-frame recognizer, then exactly five plane blowdowns and six line incidences |
| direct geometric existence | grid blow-up; all 24 graphs iff `p=3`; strong-Coble linear systems | conic plus four tangents; six residual lines iff `p=3`; strong-Coble linear systems |
| separate exact-nef gate | Kondō root exhaustion, all `202*8` cusp lifts, marked finite-volume return | remaining 35 exceptional curves, Kondō root exhaustion, all `57*5` cusp lifts, marked finite-volume return |

Thus the incidence construction supplies simultaneous effectivity of every
listed MII wall without a template-existence appeal.  It does **not** prove
by itself that no further negative curve exists.  The latter statement is a
separate theorem, already certified in the full-NS MII archive.  For Type VI
the cross-ratio construction supplies the boundary skeleton and six
residual line roots, while the other thirty-five exceptional walls still
use the published Type-VI geometry and the exhaustive cusp lift.

This separation is the reason the result is rated `L3`, not `L4`.

## 2. Raw Type MII: `735471 -> 21211 -> 1`

Start with only the even arithmetic lattice `M` of signature `(1,17)` and
its forty-eight primitive oriented mirrors.  No boundary, exceptional,
grid, or surface labels are used.  There are twenty-four short mirrors of
square `-2`.  Exhaust all eight-subsets:

```text
C(24,8) candidate frames                         735,471
pairwise-orthogonal frames                        21,211
frames making q_M/2 integral on the Picard module      1
```

Here is the exact integrality criterion.  For a proposed frame
`H={h_1,...,h_8}`, form

\[
 K_H=-\sum h_i,\qquad
 N_H=\langle 2h_i,\;v_j\ (j\notin H),\;K_H\rangle\subset M.
\]

Modulo two, the selected generators `2h_i` vanish.  Therefore `q_M/2` is
integral on `N_H` exactly when

1. `H` meets every odd-pairing edge among the forty-eight raw walls; and
2. `K_H` pairs evenly with every surviving wall generator.

The certifier verifies this equivalence twice for **all 21,211** orthogonal
frames.  One pass uses the odd-pair graph.  The hostile pass independently
forms the complete generator Gram modulo two—the surviving wall block, the
`K_H` column, and `K_H^2`—and demands that every entry be even.  Both tests
have the same unique survivor.

For that survivor,

\[
 [M:N_H]=2^7,\qquad \det(q_M/2|_{N_H})=-1,
 \qquad \operatorname{sig}N_H=(1,17).
\]

The reconstructed lattice is the odd unimodular lattice `I_(1,17)`.
Moreover `K_H` is characteristic, `K_H^2=-8`, and the raw mirrors de-scale
intrinsically into

```text
8 boundary curves of square -4
16 curves of square -1
24 curves of square -2.
```

This is a genuine raw-polytope boundary discovery, not the surface ordering
being read backwards.

## 3. The raw mirrors force the `4 x 4` grid

The sixteen recovered `(-1)` curves are pairwise disjoint.  Each meets
exactly two boundaries.  Their boundary-incidence graph is `K4,4`, whose
bipartition is unique up to switching its two sides.  Contract all sixteen.
The four boundaries on either side acquire one common isotropic image; if
these two images are `F_1,F_2`, then

\[
 F_1^2=F_2^2=0,\qquad F_1F_2=1.
\]

The contraction endpoint is therefore `P1 x P1`.  The exceptional curves
are the sixteen grid points `(a_i,b_j)`.

Every recovered `(-2)` wall meets four exceptional curves, exactly one in
each row and one in each column.  After contraction its class is `F_1+F_2`.
The twenty-four incidence profiles are all the perfect matchings of `K4,4`,
once each.  In other words, the raw long mirrors become

\[
 C_\sigma\in |\mathcal O_{P^1\times P^1}(1,1)|,
 \qquad C_\sigma\supset\{(a_i,b_{\sigma(i)}):1\le i\le4\},
 \quad \sigma\in S_4.
\]

No Kondō diagram is used to obtain this grid or the permutation list.

## 4. MII incidence is nonempty exactly in characteristic three

Assume the base field is algebraically closed and has characteristic
different from two.  A smooth `(1,1)` curve is the graph of a projectivity
between the two projective lines.  If all twenty-four `C_sigma` exist, choose
one to identify the two ordered four-point sets.  Composing with the other
twenty-three shows that **every** permutation of one four-point set is
induced by `PGL2`.

Normalize its cross-ratio to `lambda != 0,1`.  Full `S4` symmetry contains
the anharmonic transformations forcing

\[
 \lambda=-1,\qquad 2\lambda=1.
\]

Their integral ideal contains three:

\[
 2(\lambda+1)-(2\lambda-1)=3.
\]

Hence the incidence fiber is empty unless `p=3`.  In characteristic three
the solution is `lambda=2=-1=1/2`; the quadruple is projectively
`P1(F3)`, and

\[
 PGL_2(F_3)\cong S_4

\]

realizes every permutation.  Conversely, taking both quadruples to be
`P1(F3)` constructs all twenty-four graph curves.  Thus the statement is an
actual if and only if:

> The raw-MII incidence stack over an algebraically closed field of
> characteristic not two is nonempty exactly in characteristic three, and
> modulo `PGL2 x PGL2` it has one geometric point.

### Direct strong-Coble proof

Let `S` be the blow-up of `P1 x P1` at the sixteen grid points.  Let `V_i`
and `H_j` be the strict transforms of the four fibers of the two rulings.
Then they are eight disjoint smooth rational `(-4)` curves and

\[
 \sum_iV_i+\sum_jH_j
 =4F_1+4F_2-2\sum E_{ij}=-2K_S.
\]

An anticanonical divisor would push down to a `(2,2)` curve through all
sixteen grid points.  Restriction to any one of the four fibers has degree
two but four distinct zeros, so each of four parallel fibers would have to
be a component.  This is impossible in bidegree `(2,2)`.  Hence
`|-K_S|` is empty.

An anti-bicanonical divisor pushes down to a `(4,4)` curve with multiplicity
at least two at every grid point.  Restriction forces the four fibers of one
ruling to be components.  After subtracting them, the residual divisor has
the class of four fibers of the other ruling and must still pass through
each of their four coordinate values.  It is therefore exactly their sum.
Thus

```text
|-K_S| = empty,
|-2K_S| = {V_1+...+V_4+H_1+...+H_4}.
```

This proves strong-Coble existence directly from the `P`-compiled grid.

The twenty-four-wall boundary skeleton is not the desired polytope.  Exact
double description gives

```text
2674 rays = 776 finite + 434 ideal + 1464 hyperideal.
```

The twenty-four graph walls are essential: together with the skeleton they
give the full forty-eight-wall chamber.  Absence of still more walls uses
Kondō's effective-root exhaustion and the independent test of all
`202*8=1616` cusp/boundary lifts; none is integral.  The marked Alexeev
return then proves finite volume.  This terminal step is not being smuggled
into the cross-ratio lemma.

## 5. Type VI: the five plane blowdowns are visible in `P`

The independent Type-VI inverse recognizer has already found the unique
five-boundary frame without surface labels.  Starting from that reconstructed
Picard marking, enumerate every extreme ray of the sixty-wall chamber and
ask for a boundary-preserving plane class

\[
 H^2=1,\qquad K.H=-3,\qquad H.B_i>0\quad(1\le i\le5).
\]

Since `sum B_i=-2K`, the five positive boundary degrees sum to six.  They
must be `(2,1,1,1,1)`.  Among all `28,872` exact cone rays there are exactly
five such classes, one for each position of the `2`.

For the blowdown in which `B_1` has degree two, the fourteen contracted
classes are

```text
all ten E_pair_ij,
C_12, C_13, C_14, C_15 after the E_pair contractions.
```

In the original Picard lattice their total exceptional classes are

```text
E_pair_ij                                      (10),
C_1j + E_pair_1j,  j=2,3,4,5                  (4).
```

Together with `H` they form a determinant-one basis with Gram
`diag(1,-1^14)`.  Reversing the contraction gives six ordinary blow-ups at
the intersections of four lines and two consecutive blow-ups at each of
their four tangencies with a conic.  The five boundary images are therefore
a smooth conic and four distinct tangent lines.  The other four choices of
conic give the other four blowdowns.

The nineteen automatic skeleton walls have exact cone census

```text
198 rays = 9 finite + 28 ideal + 161 hyperideal.
```

Again the automatic surface skeleton is very far from finite volume.

### The six line roots force `p=3`

Parameterize the four points of tangency on the conic by

```text
t_2=infinity, t_3=0, t_4=1, t_5=lambda.
```

The raw divisor classes say, for example,

```text
C_23 is the line through L_2 intersect L_3 and t_4,t_5;
C_24 is the line through L_2 intersect L_4 and t_3,t_5.
```

For the Veronese conic `xz-y^2=0`, the tangent at a finite parameter `t` is
`x-2ty+t^2z=0`.  The first collinearity gives `lambda+1=0`; the second gives
`lambda-2=0`.  The six complementary-pair conditions are

```text
lambda+1,
lambda-2,
2*lambda-1,
2*lambda-1,
lambda*(lambda-2),
(lambda-1)*(lambda+1).
```

With `lambda != 0,1`, all six vanish exactly when

```text
characteristic = 3,  lambda = 2 = -1.
```

Thus one raw plane blowdown already recovers the Type-VI characteristic.
The other four blowdown rays coexist in the same raw chamber.  On the exact
surface they are nef square-one classes and hence semiample by the general
nef-semiampleness theorem for rational surfaces with `|-K|` empty but a
positive anticanonical multiple effective.  A separate coordinate-level
construction of all five plane maps is not needed for the lattice census,
but the present cross-ratio calculation should not be advertised as proving
the remaining thirty-five exceptional curves.

### Direct strong-Coble skeleton

Blow up the six line-line intersections and twice at each conic-line
tangency.  The five strict transforms are disjoint `(-4)` curves and sum to
`-2K`.  A cubic satisfying the fourteen simple cluster conditions would
meet each tangent line in at least its three line intersections plus the two
successive tangency conditions, more than degree three; it would contain all
four lines, impossible.  Hence `|-K|` is empty.

For a sextic with double cluster conditions, the same infinitely-near
Bezout count forces all four lines as components.  After subtracting them,
the residual conic has the prescribed tangent direction at all four points.
Its intersection with the original conic is then at least eight, exceeding
`2*2`, unless it is that conic.  Therefore the boundary is the unique member
of `|-2K|`.

This gives a one-dimensional automatic family before the residual-root
conditions and a single characteristic-three incidence point afterward.
Exact return of all sixty walls still invokes the independent thirty-five
cusp lifts and no-hidden-root theorem.

## 6. The general finite blowdown/incidence compiler

The two calculations instantiate a general construction that is rigorous
but has one important terminal gap.

Let `P` pass the surface-free Coble numerical gates C0--C2, with a finite
candidate list of boundary, `(-1)`, and `(-2)` wall classes.  For a
`(-1)`-curve `E`, blowing it down changes the numerical data by

\[
 D'.F'=D.F+(D.E)(F.E),\qquad K'=K-E
\]

after identifying the new Picard lattice with `E^perp`.  Enumerate every
wall that becomes a smooth numerical `(-1)` candidate, update, and repeat.
Every branch has length at most the Picard rank, so this is a finite tree.

Why is it exhaustive in the direction that matters?  If `P` really is the
**complete** nef chamber of a strong Coble surface, every curve appearing in
a minimal-model contraction sequence pulls back to an irreducible negative
curve on the original surface.  Completeness then puts it among the facets
of `P`.  Hence some enumerated branch reaches `P2` or a Hirzebruch surface.
This argument is exhaustive among exact carriers; it does not say that a
partial wall list contains a blowdown of every weak realization.

For a fixed branch, reversing the contractions gives a finite proximity or
Enriques diagram.  The space of strict clusters with a fixed diagram is a
smooth finite-type scheme.  Over its universal iterated blow-up, relative
Hilbert schemes or projectivized relative linear systems parameterize the
finitely many pushed-down wall curves.  After finite cohomology/base-change
stratification, the following are open or closed algebraic conditions:

```text
specified multiplicities and infinitely-near incidences;
smoothness, geometric irreducibility, and distinctness;
simultaneous effectivity of all listed walls;
h0(-K)=0 and h0(-2K)=1;
the required boundary divisor is the unique anti-bicanonical member.
```

Thus `P` compiles a finite union of finite-type incidence stacks.  Over a
fixed algebraically closed computable field, nonemptiness of every explicit
branch is a saturated polynomial-ideal question.  Over `Spec Z`, its image
also records the generic and exceptional characteristics.  The MII and VI
equations above are small exact examples of that characteristic feedback.

The cluster-space assertion is supported by Kleiman--Piene--Tyomkin's
[representability theorem for fixed Enriques diagrams](https://arxiv.org/abs/0905.2169).
Their Hilbert map is an embedding in characteristic zero but can be purely
inseparable in positive characteristic.  Consequently the compiler must use
the cluster functor itself in characteristics three and five; it may not
silently transplant a characteristic-zero Hilbert-stratum argument.

## 7. The no-hidden-wall obstruction is genuinely separate

The finite incidence stack decides whether all **listed** divisor classes
can be effective simultaneously.  It does not automatically decide that
they are all negative curves on the surface.

A hidden `(-1)` or `(-2)` class ranges over the entire Picard lattice.  In
an indefinite rational-surface lattice there is no P-derived bound currently
known on its degree relative to an interior polarization.  Therefore
“there exists a hidden curve” is a countable union of relative Hilbert
images, not automatically one finite-type condition.  Replacing it by
“the guessed cone is the nef cone” merely restates the desired conclusion.

There are three honest ways to close this gate:

1. prove a height bound from `P`, then enumerate every possible negative
   class;
2. provide and verify a finite Cox-ring presentation whose GIT chamber is
   exactly `P`; or
3. use an independent exhaustive surface theorem, such as Kondō's
   effective-root classification plus the complete cusp-lift calculation.

Computational Cox-ring algorithms can verify a proposed presentation by
saturation and graded primality, but their termination is not automatic for
an arbitrary blow-up; see Hausen--Keicher--Laface,
[*Computing Cox rings*](https://arxiv.org/abs/1305.4343).  This is precisely
why a finite blowdown stack is not yet the flagship theorem.

## 8. What the positive-dimensional families do—and do not—prove

There are now two natural, non-singleton automatic surface families:

```text
conic + four tangent lines on P2:       dimension 1;
four fibers x four fibers on P1xP1:    dimension 2.
```

Every member of either family is a strong Coble surface by the direct linear
system arguments above.  This is genuine `POLYTOPE -> SURFACE` leverage at
the boundary-skeleton level.  But their skeleton Coxeter cones have 161 and
1,464 hyperideal rays.  They are infinite-volume and are not the full nef
cones.  The special cross-ratio equations create the extra Type-VI/MII wall
configurations; Kondō's theorem and the cusp audits then exclude all hidden
walls.

Accordingly the exact status is

```text
raw Type-MII -> unique surface incidence -> exact same P:  proved, one P;
raw Type-VI  -> plane incidence -> exact same P:           proved, one P,
                                                           with 35-wall
                                                           effectivity imported;
automatic positive-dimensional strong-Coble skeletons:    proved;
non-singleton raw finite-volume class with automatic
all-wall effectivity and no-hidden-wall return:             open.
```

Calling the union of the two successful raw chambers a “class” would be a
two-template reformulation, not `L4`.  The next real completion step is to
find an intrinsic blowdown condition that forces the additional wall stock
and supplies a uniform no-hidden-curve theorem on a denominator chosen
before these two answers are known.

## 9. Replay

```bash
PYTHONPATH=src:. python3 scripts/certify_coble_polytope_first_blowdown_incidence_bridge.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_coble_polytope_first_blowdown_incidence_bridge.py
```

Primary surface input used only at the separately marked completeness gate:

- Shigeyuki Kondō,
  [*Coble surfaces with finite automorphism group*](https://arxiv.org/abs/2102.00591).
- Antonio Laface and Damiano Testa,
  [*Nef and semiample divisors on rational surfaces*](https://arxiv.org/abs/1104.4270).


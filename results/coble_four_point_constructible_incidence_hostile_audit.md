# Hostile audit of the four-point constructible-incidence theorem

**Status (2026-08-09): the two explicit four-point inverse calibrations
survive; the advertised extension-closed class theorem does not.  Strict
rating: `L3` two-member inverse calibration, not `L3+` and not `L4`.**

The production theorem contains a real and useful result.  Starting from two
raw Coxeter wall systems, it recovers two different presentations of the same
four-point finite geometry, realizes both over `P1(F3)`, constructs strong
Coble root seeds, and hands them to the general `D_CobVin` exact-return
theorem.  I independently replayed that calculation without importing the
production helpers or reading its precompiled blowdown bridge.  The two cores
are sound.

The original claimed closure under arbitrary “standard-free extensions,”
however, used a false Bertini implication.  Base-point-free plus a generically
separable associated map does not imply that the general member is
geometrically irreducible.  Both certified rows have zero extension steps, so
this gap does not damage either explicit surface; it does invalidate the
advertised theorem on every member of the extension-closed class.

The production artifacts have now been corrected: `D_CobVin^CI4` is
explicitly core-only, its rating is `L3`, and the repaired extension verifier
is marked prospective.

The exact independent executable is
[`audit_coble_four_point_constructible_incidence_hostile.py`](../scripts/audit_coble_four_point_constructible_incidence_hostile.py),
with frozen
[`JSON`](../data/enumerations/coble_four_point_constructible_incidence_hostile_audit.json)
and
[`tests`](../tests/test_coble_four_point_constructible_incidence_hostile_audit.py).

## 1. Bottom-line five-gate verdict

| gate | hostile verdict | reason |
|---|---|---|
| independent `C` | **partial pass** | the proposed certificate language is phrased in P-derived data, but the executable reads exactly two raw archives plus a precompiled bridge; it has no arbitrary-input interface, generic dependency canonicalizer, or finite extension-certificate verifier |
| `Car` nonempty for every member | **fail for the advertised extension class; pass for the two cores** | the extension Bertini step is false, and both certified rows have zero extensions |
| exact `Return` | **conditional pass** | once a genuine strong-Coble root seed exists and the `D_CobVin` gates hold, Vinberg/Kondō/Halphen exhaustion returns literally the marked input `P` |
| independent `AG` | **conditional pass** | Kondō is invoked after construction and finite-index reflection theory, not in class membership |
| new `Pushback` on every member | **strict fail; formal qualification only** | on the only two certified members the complete core already fixes characteristic three, `b=5` or `8`, and hence `H14` or `H17` |

In tuple form this is

```text
(C, Car, Return, AG, Pushback)
= (partial, fail/pass-two-cores, conditional-pass,
   conditional-pass, fail).
```

The defensible theorem level is therefore

```text
L3 two-member inverse calibration;
not an extension-class theorem;
not L4.
```

This is not a claim that the construction merely renamed two matrices.  The
raw reconstructions and common finite-geometry theorem are substantive.  The
more limited point is operational: the current accepted denominator contains
only the two known zero-extension configurations, and its attempted route to
further members is not yet proved.

## 2. Independent symbolic replay

Normalize four points of `P1` as

\[
                         (\infty,0,1,\lambda),
                         \qquad \lambda\ne0,1.
\]

### 2.1 The complete graph design

For every permutation of the four points I independently formed the `4 x 4`
evaluation matrix of `H0(O(1,1))`.  The determinant histogram is exactly

| determinant | multiplicity |
|---|---:|
| `0` | 4 |
| `-(lambda-1)(lambda+1)` | 4 |
| `2lambda-1` | 4 |
| `lambda(lambda-2)` | 4 |
| `lambda^2-lambda+1` | 4 |
| `-(lambda^2-lambda+1)` | 4 |

On the stated domain, simultaneous vanishing includes
`lambda+1=0` and `2lambda-1=0`; the integer identity

\[
                 2(\lambda+1)-(2\lambda-1)=3
\]

forces characteristic three.  Conversely, at `lambda=2` in `F3`, all 24
matrices have rank three.

The hostile replay adds a geometric check omitted by the production
certificate.  It computed each one-dimensional modular kernel, normalized it
projectively, and obtained 24 distinct bilinear sections.  For a section

\[
 a x_0y_0+b x_0y_1+c x_1y_0+d x_1y_1,
\]

all 24 coefficient determinants `ad-bc` are nonzero in `F3`.  Hence these are
not merely nonzero sections: they are 24 distinct smooth irreducible `(1,1)`
graphs.

### 2.2 The complete polar/chord design

For the conic `xz-y^2=0`, the six polar/chord equations independently give

```text
(lambda-1)(lambda+1)  twice;
lambda(lambda-2)      twice;
1-2lambda             twice.
```

The domain and

\[
                    (\lambda+1)-(\lambda-2)=3
\]

again force characteristic three and `lambda=2`; all six evaluation matrices
then have rank two.  The six kernel lines were also computed explicitly.
They are pairwise distinct and none coincides with one of the four tangent
boundary lines.  Thus the residual-line geometry, not just its rank count, is
certified.

This proves the common four-point theorem exactly as advertised:

\[
 \text{either complete design is realizable}
 \iff \operatorname{char}k=3,
 \quad Q\simeq\mathbf P^1(\mathbf F_3).
\]

## 3. Raw dependency extraction without the bridge

### 3.1 Rank 18, 48 walls

Starting only with the raw even lattice Gram and 48 mirror vectors, the audit
recomputed

```text
short mirrors                                  24
all eight-subsets                         735,471
pairwise-orthogonal eight-frames           21,211
odd pairing edges                              32
orthogonal frames covering every odd edge       1
integral de-scaling frames                       1
```

The unique frame is raw indices `0,...,7`.  A fresh Hermite-normal-form
de-scaling gives an integral Gram of determinant `-1` and exact inertia
`(1,17,0)`, hence `I_(1,17)`.  Its physical walls split as

```text
8 disjoint (-4) boundaries;
16 disjoint (-1) exceptional curves;
24 (-2) roots.
```

Boundary incidence of the 16 exceptional classes is exactly `K4,4`.
Contracting them sends each side to one ruling, whose Gram is

\[
                         \begin{pmatrix}0&1\\1&0\end{pmatrix}.
\]

Every root meets a perfect matching of the `4 x 4` centres, pushes down to
class `(1,1)`, and the 24 profiles are all elements of `S4`, once each.  No
surface label or frozen bridge output was used.

### 3.2 Rank 15, 60 walls

Starting only with the full raw `I_(1,14)` Gram and 60 wall vectors, square
and boundary-intersection tests reproduce

```text
5 disjoint (-4) boundaries;
10 actual (-2) roots;
10 pair-boundary (-1) crossing curves;
35 single-boundary (-1) cusp lifts.
```

The ten crossings index all two-subsets of the five boundaries, and each has
a unique associated root.  For each possible conic boundary, I formed the 14
exceptional classes—ten crossings and the four total tangency exceptionals—
and found the primitive generator of their common orthogonal directly.  In
all five cases it has

\[
 H^2=1,\qquad K.H=-3,\qquad
 \{H.B_i\}=\{2,1,1,1,1\}.
\]

The plane class plus the 14 exceptionals is a determinant-one basis with Gram
`diag(1,-1^14)`.  Four roots contract in the tangency chains; each of the six
remaining roots is a plane line supported on one ordinary pair centre and
the two complementary tangency centres.  Every two-subset/complement occurs
once, for every one of the five plane branches.

Again this was reconstructed without reading the bridge on which the
production CI4 script relies.

## 4. The two strong-Coble constructions survive

For the grid core, blow up the sixteen intersections of four fibres in each
ruling of `P1 x P1`.  For the polar core, blow up the six tangent-line
intersections and twice along each of four conic tangencies.  The preceding
computations give the claimed odd-unimodular exceptional bases and every root
seed as a smooth irreducible curve.

In either case the boundary `B=sum B_i` is reduced, smooth componentwise,
disjoint, and satisfies

\[
                         B=-2K,\qquad B_i^2=-4.
\]

If `D` lies in `|-2K|`, then `D.B_i=-4`; negativity forces every `B_i` to be
a component of `D`, and the divisor class is exhausted.  Hence `D=B`.  If
`A` lay in `|-K|`, then `2A` would lie in `|-2K|`, so `2A=B`, impossible
because `B` is reduced.  Therefore

\[
                         |-K|=\varnothing,
                         \qquad|-2K|=\{B\}.
\]

This argument also applies to any future extension that is actually
constructed and retains the exact reduced, disjoint, `B=-2K` hypotheses.  It
does not by itself prove that such an extension exists.

## 5. The standard-free extension lemma is false as stated

The pre-audit production report made the following implication twice:

> Off the prescribed cluster, the residual system is base-point-free and
> defines a generically separable map; therefore Bertini gives a nonempty open
> of smooth geometrically irreducible members.

Here is an exact counterexample in the very characteristic under discussion.
Let

\[
               X=\mathbf P^1\times\mathbf P^1
\]

over an algebraically closed field of characteristic three and take

\[
                         L=\mathcal O_X(2,0)=2F_1.
\]

The complete linear system is base-point-free.  Its associated morphism is
the first projection followed by the quadratic Veronese embedding,

\[
 X\xrightarrow{\operatorname{pr}_1}\mathbf P^1
  \xrightarrow{\nu_2}\mathbf P^2.
\]

The projection is smooth and separable and the Veronese map is a closed
immersion, so the associated map is generically separable.  A general section
is a square-free binary quadratic.  Over the algebraic closure it factors
into two distinct linear factors, and its divisor is two disjoint ruling
fibres.  It is smooth but disconnected and geometrically reducible.

Thus separability controls generic smoothness, not irreducibility when the
linear system is composed with a pencil.

At minimum the extension definition must add:

1. image dimension at least two, or a direct certificate that the generic
   member is geometrically integral;
2. for a pencil, a geometrically integral generic fibre, not merely a
   separable map;
3. relative constant rank and local freeness of the jet-kernel family;
4. irreducibility/dominance of the sequential parameter tower, with every
   nonempty open lying on the same component; and
5. explicit, replayable smoothness, factorization/irreducibility,
   transversality, and distinctness witnesses.

The JSON wording that a residual system “has a separable Bertini open set”
is ambiguous.  If this phrase literally assumes a nonempty open of
geometrically integral members, the logical implication is repaired, but
membership then contains much of the desired existence conclusion and still
lacks a finite verifier schema.  The formal report instead states only
base-point-freeness plus a generically separable map, which is insufficient.

Both archived signatures report

```text
standard_free_extension_steps_after_core = 0.
```

Consequently this failure does not touch the two core surfaces.  It removes
the claimed uniform theorem for the larger extension-closed class.  The
corrected production theorem no longer claims that closure.

## 6. Is membership genuinely P-only and input-computable?

There are two different answers.

Mathematically, the proposed language can be made P-only: enumerate numerical
contractions, derive incidence hypergraphs from the intersection table,
recognize a finite projective design, and verify explicit algebraic
certificates.  No Kondō label or output Gram is intrinsically necessary.

Operationally, that algorithm is not the current artifact.  Static inspection
of the pinned production executable finds only the options

```text
--check-archive
--write
```

and four hardwired path constants

```text
BRIDGE, TYPE_VI, TYPE_MII, OUTPUT.
```

It has no arbitrary wall-matrix input, all-contraction enumerator, generic
dependency canonicalizer, or verifier for a standard-free extension
certificate.  It also reads the precompiled bridge rather than extracting
either dependency graph itself.

The correct verdict is therefore:

```text
conceptually P-only certificate language;
not yet an implemented input-computable membership test.
```

## 7. Does the D_CobVin handoff return the exact same P?

Yes, conditionally on a genuine constructed root seed.  This part is not the
source of the downgrade.

The exact logic is:

1. the constructed curves realize every raw Coble--Mukai root seed with the
   original integral marking;
2. the raw Vinberg gate gives finite index of their reflection group;
3. finite index gives finite automorphism group, and Kondō's root-exhaustion
   proposition excludes hidden effective irreducible CM roots;
4. Kondō's CM-nef criterion makes each primitive ideal ray of the complete
   root chamber nef on the surface;
5. primitive-isotropic and Halphen theory realizes every integral
   single-boundary lift;
6. adjunction sends any hidden `(-1)` curve either to a hidden two-boundary
   root or to one of those primitive cusp lifts; and
7. all negative curves, hence every nef inequality, are exhausted with the
   original orientation and marking.

Therefore the conclusion is literal equality of marked chambers, not merely
an abstract isomorphism.

The pinned raw audits independently keep finite volume separate from this
argument:

| core | CM roots | CM cusps | integral lifts | full-cone rays |
|---|---:|---:|---:|---:|
| `H14/60` | 20 | 57 | 35 | `27,845` finite + `1,027` ideal + `0` hyperideal |
| `H17/48` | 40 | 202 | 0 | `190,808` finite + `5,244` ideal + `0` hyperideal |

For `H14`, the literal Vinberg audit also has zero bad edges and zero bad
parabolic neighborhoods.  For `H17`, the reflection index in `O+` is
`1,152`.  No finite-volume claim is inferred merely from finite sidedness or
from the presence of a cusp.

The qualification is artifact-level: the CI4 executable states this handoff
in prose rather than replaying it.  The pinned general `D_CobVin` theorem and
the two raw cone archives supply the actual conditional certificate.

## 8. What does Kondō add on every current member?

The logical order is clean.  After construction and exact return, finite
automorphism follows from reflection finite index, and only then is Kondō's
classification applied.  It is not hidden in P-membership.

Formally, characteristic three leaves seven rows of the nine-class theorem
and dimensions

\[
                           \{10,11,14,17\}.
\]

For the current CI4 denominator, however, the polar core already has five
boundary components and rank 15, while the grid core already has eight
boundary components and rank 18.  Thus they already encode `H14` and `H17`
before Kondō is invoked.  Classification identifies those known rows but
does not materially narrow either individual input.

This gate would become nontrivial if a third P-derived design or a repaired
genuine extension entered the class without already fixing its classified
endpoint.  Until then, the pushback is logically independent but
mathematically thin.

## 9. A separate enlarged finite-geometry core class

The rank-12, 206-wall Type-MI datum supplies a genuinely different third
core.  This does not repair the false CI4 extension lemma; it suggests a
better replacement: a **core-only** finite-geometry class with no extension
closure until the Bertini certificates are fixed.

Ignoring every stored surface label and stored block list, the hostile replay
reads only the physical Gram and 206 wall vectors.  It intrinsically recovers

```text
2 disjoint (-4) boundaries;
30 actual (-2) roots;
10 pairwise-disjoint (1,1)-profile (-1) curves;
164 single-boundary cusp-lift walls.
```

The common orthogonal of the ten crossing curves is a saturated copy of `U`;
with those ten exceptionals it forms a determinant-one basis.  Contraction
sends the boundaries to bidegrees `(1,3)` and `(3,1)` and every root to
bidegree `(1,1)`.  Each root passes through four crossing centres.  The thirty
four-subsets satisfy

```text
every triple of the ten centres belongs to exactly one block,
```

so the raw wall datum itself reconstructs the Steiner system
`S(3,4,10)`.

### 9.1 Independent finite-field realization

I built `F9=F3[r]/(r^2-2)` from scratch, enumerated all `720` elements of
`PGL2(F9)`, and obtained its thirty images of `P1(F3)`.  A direct design
isomorphism transports the raw thirty blocks to exactly these thirty
`F3`-sublines in `P1(F9)`.

On `P1 x P1`, take the two Frobenius graphs

\[
                         y=x^3,\qquad x=y^3.
\]

They meet transversely at `x=x^9`, the ten reduced points of `P1(F9)`.  For
each `F3`-subline, there is exactly one Möbius transformation in
`PGL2(F9)` agreeing with `t -> t^3` on its four points.  The hostile replay
found thirty distinct such transformations; their graphs are thirty distinct
smooth irreducible `(1,1)` curves.  Blowing up the ten centres gives the
complete two-boundary strong-Coble root seed.  The existing `D_CobVin`
theorem then generates all 164 integral cusp lifts and returns the literal
`H11/206` chamber, whose independent cone census is

```text
327,566 finite + 18,354 ideal + 0 hyperideal.
```

The finite-field Type-MI geometry is published in Kondō and
Katsura--Kondō.  The defensible project-level gain is the label-free
`P -> contraction -> Steiner design` extraction and its uniform
`P -> surface -> same P` use, not discovery of the surface, the finite
geometry, or a new abstract polytope.

### 9.2 The parametric Möbius factory and its exact limit

The construction has a clean prime-power parameter.  For a prime power `q`,
the graphs

\[
                         y=x^q,\qquad x=y^q
\]

meet at the `q^2+1` points of `P1(F_(q^2))`.  The images of `P1(F_q)` form

\[
              S(3,q+1,q^2+1),
              \qquad q(q^2+1)\text{ blocks}.
\]

After blowing up all centres, direct divisor arithmetic gives

\[
\begin{aligned}
 B^2&=-(q-1)^2,\\
 C^2&=1-q,\\
 K.C&=q-3,\\
 B+B'&=-2K+(q-3)(F_1+F_2).
\end{aligned}
\]

Consequently the simultaneous Coble/root identities

```text
B^2=-4,  C^2=-2,  K.C=0,  B+B'=-2K
```

hold exactly at `q=3`.  This is a useful necessary-and-sufficient criterion
inside a parametric algebro-geometric factory.

For `q!=3`, the construction has simply left the strong-Coble `(-4)/(-2)`
grammar.  It proves **no nonexistence theorem for an abstract Coxeter
polytope** in those parameters or dimensions.  Any ledger statement must
keep that distinction explicit.

### 9.3 Provisional enlarged-class rating

Define `D_CobVin^FGcore` structurally: a raw member of `D_CobVin` whose
P-derived contraction hypergraph is either one of the complete CI4 cores or
a complete Möbius-plane `S(3,q+1,q^2+1)` core satisfying the Coble numerical
alphabet.  Do not include the broken standard-free extension clause.

If the production artifact freezes a genuinely arbitrary-input recognizer
for this definition, the formal gates become

| gate | enlarged core-only verdict |
|---|---|
| `C` | scoped pass: a parametric finite-geometry grammar, not a list of Grams |
| `Car` | pass: explicit `P1(F3)` and `P1(F9)` constructions for every accepted core |
| `Return` | pass through the marked `D_CobVin` theorem |
| `AG` | pass after carrier construction |
| `Pushback` | formal but thin pass: Kondō applies class-wide, while each rigid core already fixes much of the rank/boundary data |

That would warrant the deliberately qualified label

```text
scoped formal L4 core class;
not broad or flagship L4;
original CI4 extension theorem still false.
```

This provisional promotion depends on the production recognizer being
structural and input-computable.  **The currently frozen remaining-seven
executable does not meet that condition:** it is hardwired to six archives
containing the seven classified rows, has no arbitrary-input option, and
selects the Möbius recognizer by its internal `two_boundary_206` output slot.
Its extraction after selection is genuinely label-free, and its theorem
definition does not compare an output Gram, but it is not yet a normalized
raw-input membership implementation.

Therefore the actual frozen rating remains

```text
L3 structural P-first calibration.
```

The `scoped formal L4` label above is a precise conditional target, not a
current claim.  It is earned only when the same normalized recognizer is run
on an arbitrary supplied `D_CobVin` wall datum and the core-only theorem is
frozen independently of the seven known archives.

## 10. Exact repair path

There are two clean ways forward.

1. Restrict the theorem to the two complete four-point cores.  This yields a
   correct `L3` inverse-calibration theorem immediately.
2. Repair the extension certificate as in Section 5, implement a genuine
   arbitrary-P dependency extractor/verifier, and exhibit at least one
   nonzero extension or a third structural core.  Only then reassess a scoped
   formal `L4` claim.

A finite union of named output Grams would not solve the latter problem.  A
new finite-geometry grammar can: it must be recognizable from raw `P`, have
an explicit uniform nonemptiness theorem, construct all root seeds, and enter
the same marked `D_CobVin` return argument.  Whether such an enlarged class is
formally `L4` then depends on whether its recognizer accepts an independently
defined structural locus rather than simply spelling out the known answers.

## 11. Replay

```bash
python3 scripts/audit_coble_four_point_constructible_incidence_hostile.py --check-archive
python3 -m pytest -q tests/test_coble_four_point_constructible_incidence_hostile_audit.py
```

The hostile archive freezes the exact modular kernels, both raw lattice
reconstructions, all five plane contractions, the strong-Coble divisor
argument, the extension counterexample, the conditional marked return, and
the strict five-gate rating.

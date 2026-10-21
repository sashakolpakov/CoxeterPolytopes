# A two-input raw-polytope-to-Coble promotion theorem

**Status (2026-08-09): proved exact-return criterion with two executable
calibrations; `L3`, not strict `L4` under the strengthened five-part
contract.**

This is the first genuinely non-singleton instance in this project in which
the logical route is

```text
arbitrary raw (P,h) in an independently defined surface-free denominator
  -> exhaustive finite K3 lift compiler
  -> passing branch if and only if an exact carrier exists
  -> complex surface existence
  -> algebraic-geometric classification feedback
  -> the full nef cone is exactly the original P.
```

The carrier-admitting locus has exactly two underlying polytope isometry
types.  They are the
22- and 30-facet noncompact chambers in `H10` which, **after** the construction,
are identified with the one-boundary Coble Types I and II.  Thus this is a new
polytope-first construction and criterion for two known positive chambers, not
a new abstract Coxeter polytope.

The exact archive is
[`coble_one_boundary_raw_k3_promotion.json`](../data/enumerations/coble_one_boundary_raw_k3_promotion.json),
rebuilt by
[`certify_coble_one_boundary_raw_k3_promotion.py`](../scripts/certify_coble_one_boundary_raw_k3_promotion.py)
and tested by
[`test_coble_one_boundary_raw_k3_promotion.py`](../tests/test_coble_one_boundary_raw_k3_promotion.py).

## 1. The surface-free denominator, without surface templates

Fix a marked oriented integral Coxeter datum, with its integral strict
interior orientation vector,

\[
                 P=(N;w_1,\ldots,w_m;h),
                 \qquad N\simeq I_{1,10}.
\]

The denominator `D_1_split` consists of the data satisfying the following
finite conditions.

1. The `w_i` are primitive, their reflections are integral, and their
   irredundant cone is a finite-volume Coxeter polytope.  Exact double
   description has no hyperideal extreme ray.
2. There is a unique facet normal `b` which is characteristic and has square
   `-1`.
3. Exactly twelve normals `u_i` have square `-2`, are orthogonal to `b`, and
   span `b^perp` over `Q`.  Their cone is a finite-volume Coxeter chamber in
   the even unimodular lattice `b^perp ~= E10`.
4. The twelve-root chamber has at least one ideal ray.  Every remaining facet
   is a norm-`-1` vector `E` with `b.E=1`, and the primitive vectors `f=E+b`
   are **exactly all** its ideal extreme rays.
Nothing in `D_1_split` asks that a lift exist.  On this independent
denominator define the finite compiler `C_1_split(P,h)` below.  A branch passes
only if it has a primitive equivariant K3 graph glue, its full fixed lattice is
exactly the primitive marked `N(2)` with Nikulin invariants `(11,11,1)`, and
`h` passes the complete root-regular and bounded irreducibility test for all
twenty-four proposed split roots.  The accepted locus

```text
R_1_split = {(P,h) in D_1_split : C_1_split(P,h) is nonempty}
```

is compiler output, not a hypothesis disguised as a class definition.  All
these tests are arithmetic: neither `D_1_split` nor `C_1_split` assumes a K3
surface, a Coble surface, effectivity, a ramification curve, absence of hidden
surface curves, or equality with a nef cone.  The period map, Torelli,
ramification, and exact return are conclusions below.

Including `h` is not surface input: an oriented marked hyperbolic chamber is
routinely supplied with a timelike interior vector.  It also removes an
otherwise artificial existential search.  For a fixed decorated input every
branch is decided by a single finite positive-definite enumeration.  If the
same chamber is an exact nef cone, its strict interior `h` is ample; hence the
actual K3 cover necessarily supplies a passing branch for this same `h`.

### Quantifiers

The theorem is not a replay statement about two named matrices:

> **Theorem A (polytope-to-surface iff).** For every marked, `h`-decorated
> `(P,h) in D_1_split`, there exists an exact complex finite-automorphism
> strong one-boundary Coble carrier taking `h` to an ample class if and only
> if `C_1_split(P,h)` has a passing branch.  Every passing branch constructs
> a surface `S(P)` and an integral marking `NS(S(P)) ~= N` for which
> \[
>             \mathbf P(\operatorname{Nef}S(P))=P.
> \]

Necessity comes from the actual K3 cover: ampleness makes `h`
root-regular and the actual lifted curves are simultaneously irreducible.
Sufficiency is the construction below.  Classification then proves that, after forgetting
the harmless `h` decoration, the accepted polytopes have exactly two isometry
types.

The executable freezes the exhaustive replay for the two resulting positive
types.  The **generic** compiler and its termination proof are mathematical;
the program does not pretend to enumerate every as-yet-unsupplied raw wall
matrix in the a priori class.

## 2. Why the compiler is finite and exhaustive

Each of the twelve `(-2)` walls must split on the K3 cover into two disjoint
exchanged roots.  Choose one root over each wall and put

\[
 c_i^++c_i^-=\pi^*u_i,\qquad
 \delta_i=c_i^+-c_i^-.
\]

If `u_i.u_j=d`, nonnegative same-sheet and cross-sheet intersections give

\[
 c_i^+.c_j^+=a,\quad c_i^+.c_j^-=d-a,quad
 \delta_i.\delta_j=4a-2d,qquad 0\le a\le d.       \tag{2.1}
\]

Thus every pair contributes only `d+1` possibilities.  A finite input has a
finite product of such possibilities.  Independent exchanges of the two
sheets act by diagonal sign switching, so the compiler quotients every
labelled assignment without losing a branch.

For each difference Gram it then performs, in this order:

1. exact inertia and radical-parity tests;
2. the HNF saturation forced by all integral half-sums;
3. the anti-invariant rootlessness and rank budget;
4. all even equivariant finite-index closures, equivalently isotropic
   subgroups of finite discriminant groups;
5. all primitive K3 complements and discriminant anti-isometries, followed
   by every graph glue, the extended involution, the saturated fixed lattice,
   and the marked-`N(2)` primitivity test; in rank 20 the complement is a
   positive even binary form, whose reduced forms of fixed determinant are
   finite;
6. the supplied integral strict interior `h`; and
7. one complete positive-definite short-vector enumeration for each branch.

This is a terminating membership algorithm.  It enumerates all foldings and
K3/Coble carriers compatible with the stated one-boundary split-root
dictionary.  It does **not** enumerate carriers by unrelated algebraic
species; no conclusion below obstructs such a carrier.

The rank bookkeeping is part of the general algorithm, not an unmentioned
rank-20 assumption.  If the saturated anti-invariant difference lattice has
rank ten, then `rank NS=21`, so the branch violates the complex K3 Picard
bound and is rejected.  At anti-rank nine the orthogonal complement is the
positive binary lattice explicitly treated below.  At smaller anti-rank the
complement has higher rank; standard finite discriminant-form, genus, and
Nikulin primitive-embedding algorithms enumerate its isometry classes,
graph glues, and involution extensions.  The same positive-definite `Q_h`
test then decides root regularity and bounded irreducibility.  Thus the
mathematical compiler covers every possible anti-rank, although the supplied
executable only calibrates the rank-nine survivors of the two archived
inputs.

For the two positive raw inputs the complete split denominators are:

| raw input | labelled Grams | sheet-switch orbits | budget survivors | final admissible orbit | labelled assignments in it |
|---|---:|---:|---:|---:|---:|
| 22 facets | `27,648` | `18` | one rank-9 plus one rank-10 residue | **1** | `1,024` |
| 30 facets | `32,768` | `16` | two rank-9 residues | **1** | `2,048` |

For 22 facets, the rank-nine orbit saturates to
`E8(-2)+<-4>`; the rank-ten residue violates the complex K3 rank budget and
also creates a saturated anti-invariant `(-2)`-root.  For 30 facets, the two
rank-nine residues saturate to `D9(-2)` and `<-2>^9`; exact short-vector
enumeration rejects the second because it has eighteen signed
anti-invariant `(-2)`-roots.  Hence the construction consumes **all**
survivors, not a preferred familiar one.

## 3. The two rank-20 equivariant lattices

Put `M=N(2)`.  In each surviving difference lattice `A`, form

\[
 L=\left\langle M\mathbin\perp A,
       \frac{u_i+\delta_i}{2}:1\le i\le12\right\rangle .       \tag{3.1}
\]

The new boundary direction is the orthogonal norm-`-2` vector corresponding
to `b`; equivalently these are the raw boundary extensions of the two
rank-19 Enriques lift packages.  HNF gives index `512` over the invariant
plus anti-invariant direct sum in both cases.

| raw input | anti lattice | `rank L` | `det L` | Smith diagonal |
|---|---|---:|---:|---|
| 22 facets | `E8(-2)+<-4>` | 20 | `-8` | `1^18,2,4` |
| 30 facets | `D9(-2)` | 20 | `-16` | `1^18,2,8` |

The involution is `+1` on `M`, `-1` on `A`, and exchanges the two roots over
each `u_i`.  The certificate checks it entry by entry as an integral order-two
isometry and verifies that its fixed rank is eleven.

The generated `L` is not silently chosen among possible Picard
overclosures.  The discriminant-form census enumerates every isotropic
subgroup.  The 22-facet `L` has no nontrivial isotropic subgroup.  The
30-facet `L` has one proper even index-two overclosure, of determinant `-4`
and Smith diagonal `1^18,2,2`.  The hostile audit found, rather than hid, a
genuine extra lattice branch: its complement `diag(2,2)` has two
discriminant anti-isometries, and both graph glues are even unimodular K3
lattices on which the involution extends with fixed rank eleven.

That branch nevertheless fails the carrier compiler twice, independently.
The marked `M=N(2)` has primitive Smith data `1^10,2`, so the saturated fixed
lattice is larger; its computed Nikulin invariants are `(11,9,1)`, not
`(11,11,1)`.  Moreover the complete `Q_h` enumeration through bound `118`
returns `96` signed vectors and finds exactly `18` degree-zero `(-2)`-roots,
all anti-invariant.  An invariant ample class would be orthogonal to each of
them, which is impossible.  The same enumeration records `39` positive
degree-one roots and twelve negative-pairing witness incidences obstructing
split-root indices `6,7,14,15,22,23`.  Thus the generated `L` is the only
even closure passing the **fixed-lattice and period/ample gates**, not the only
one admitting a primitive K3 embedding.

### Primitive K3 embeddings, not a determinant guess

The certifier enumerates every Minkowski-reduced positive even binary form of
the required determinant and every discriminant anti-isometry.  The results
are:

| raw input | `A_L` | reduced complements tested | unique compatible `T` | graph-glue index |
|---|---|---|---|---:|
| 22 facets | `Z/2 + Z/4` | `diag(2,4)` | `diag(2,4)` | 8 |
| 30 facets | `Z/2 + Z/8` | `diag(2,8)`, `diag(4,4)` | `diag(2,8)` | 16 |

There are two discriminant anti-isometries in each compatible row.  The
archive records explicit fractional generators.  The involution acts as
`-1` on every element of `A_L`; extending it by `-1` on `T` therefore
preserves the graph glue.  The glued rank-22 lattice is checked to be even,
unimodular, and of signature `(3,19)`, hence is the K3 lattice.  Exact Smith
forms prove both

```text
L -> Lambda_K3 is primitive;
M=N(2) -> Lambda_K3 is primitive.
```

The extended involution has fixed rank eleven.  Since the displayed `M` is
primitive and already has that rank, it is the entire fixed lattice, not a
finite-index sublattice.  Its invariants are `(r,a,delta)=(11,11,1)`.

## 4. The period and simultaneous irreducibility gate

Let `G` be the Gram matrix of `L`, let `h^2=H>0`, and define

\[
       Q_h(x)=2(h.x)^2-Hx^2,
       \qquad Q_h=-HG+2(Gh)(Gh)^t.                      \tag{4.1}
\]

The Hodge index theorem makes `Q_h` positive definite.  For every K3 root
`x^2=-2` of degree `d=h.x`,

\[
                       Q_h(x)=2H+2d^2.                 \tag{4.2}
\]

Consequently:

- enumeration through `2H` sees **every** root orthogonal to `h`, invariant,
  anti-invariant, or mixed; and
- if the largest proposed degree is `D`, enumeration through
  `2H+2(D-1)^2` sees every positive root of smaller degree.

The exact results are:

| raw input | `h^2` | proposed degrees | PARI bound | signed vectors | degree-zero roots | irreducibility obstructions |
|---|---:|---|---:|---:|---:|---:|
| 22 facets | 138 | `20*1,4*8` | 374 | 474 | 0 | 0 |
| 30 facets | 58 | `18*1,6*2` | 118 | 36 | 0 | 0 |

Every returned vector is replayed with exact integer arithmetic; in these
two runs all returned vectors are `(-2)`-roots.  The full degree histograms
and vector-set hashes are archived.

Why does this prove simultaneous irreducibility?  Surjectivity of the refined
(polarized) K3 period map---equivalently, choosing the marking and Kähler
chamber---realizes the prescribed root-regular chamber containing `h` as the
ample cone.  K3 Riemann--Roch then makes
every proposed positive-degree root effective.  If an effective divisor
`D` with `D^2=-2` were reducible, one of its irreducible `(-2)` components
`C` would satisfy

\[
                   0<h.C<h.D,\qquad C.D<0.
\]

The complete bounded enumeration finds no such `C` for any of the twenty-four
roots.  All twenty-four are therefore irreducible in the **same** ample
chamber.  Their pairwise intersections are nonnegative, and each prescribed
pair is exchanged and disjoint.

The absence of any degree-zero root proves that `h` is regular against the
entire root system of `L`, not only the proposed curves.  This is exactly what
is needed for the ample-chamber and Torelli step; positive-degree hidden roots
may exist, but they neither obstruct the period nor get assumed away.  Their
downstairs effect is eliminated later by algebraic geometry.

Choose the oriented period in the positive complement `T`.  Surjectivity of
the refined complex K3 period map gives `NS(X)=L` with the prescribed
root-regular chamber containing `h` as the ample cone.  The involution acts by
`-1` on the period and fixes `h`; hence it preserves that same chamber, so
global Torelli realizes it as an anti-symplectic automorphism of `X` without
conjugating the selected split-root data.

## 5. Ramification is forced; it is not a class hypothesis

Nikulin's fixed-locus formula for `(11,11,1)` gives exactly one rational
fixed curve `R`.  The quotient theorem therefore makes

\[
                  \pi:X\longrightarrow S=X/\langle\tau\rangle
\]

a smooth rational surface with one irreducible branch component `B=pi(R)`.
Here is the short argument that the word **strong Coble** is earned rather
than assumed.  The ramification formula gives

\[
                         R=\pi^*(-K_S),\qquad \pi^*B=2R.
\]

Thus `B ~ -2K_S`, and comparison of squares under the degree-two pullback
gives `B^2=-4`.

Since `R` is an irreducible `(-2)`-curve, negative intersection and induction
give `|mR|={mR}` for every `m>=1`.  If `D in |-K_S|`, then
`pi^*D=R`; pushing forward gives `2D=B`, impossible because `B` is reduced
and irreducible.  If `B' in |-2K_S|`, then
`pi^*B'=2R=pi^*B`, whence `B'=B`.  Therefore

\[
                         |-K_S|=\varnothing,
                 \qquad |-2K_S|=\{B\},
\]

so `S` is a strong Coble surface with exactly one boundary component.

There is a small integral point which must not be skipped.  Pullback embeds
`NS(S)(2)` into the fixed lattice `M`.  Both lattices have rank eleven and
absolute determinant `2^11`: `NS(S)` is unimodular because `S` is rational,
while the certificate gives `M=N(2)`.  Hence the inclusion has index one and
identifies `NS(S)` with the raw `N`.

The twenty-four irreducible roots descend in twelve exchanged pairs to
twelve irreducible `(-2)`-curves `C_i` with classes `u_i`.  Put
`b_act=-K_S=B/2`.  Adjunction gives

\[
                         b_{act}.C_i=0
                         \quad(1\le i\le12).
\]

But the `u_i` span `b^perp`.  Therefore `b_act` lies in `Qb`; since both have
square `-1`,

\[
                         b_{act}=\pm b.
\]

The fixed curve has positive intersection with the invariant ample class,
and the raw orientation has `h.b>0`, so the sign is `+`.  Thus

\[
                         -K_S=b,qquad \pi^*b=R.         \tag{5.1}
\]

This replaces the formerly separate “unique ramification frame” hypothesis.
Notice that `b=-K_S` itself is not effective: `|-K_S|` is empty.  The
effective branch divisor is `B=2b`.

## 6. Algebraic geometry returns the exact same polytope

Let `Gamma` be the twelve descended curves.  Their exact root cone is a
finite-volume Coxeter chamber.  Hence the reflection group `W(Gamma)` has
finite covolume in `O(E10)` and therefore finite index.  Since
`W(Gamma)` is contained in the group `W(S)` generated by all effective
irreducible roots, `W(S)` also has finite index.  Dolgachev--Kondō's finite
index theorem makes `Aut(S)` finite.

Now the feedback is genuinely algebro-geometric:

1. Kondō's classification says that over the complex numbers a
   finite-automorphism one-boundary Coble surface is one of exactly two
   classes.
2. Kondō's Proposition `Vinberg2` applies to the already constructed
   `Gamma`: finite automorphisms plus finite index of `W(Gamma)` imply that
   these twelve are **all** effective irreducible Coble--Mukai roots.
3. For one boundary, every `(-1)`-curve `E` gives the primitive nef isotropic
   class `f=E+b`.  Conversely every primitive ideal ray `f` of the CM chamber
   is a genus-one fibration.  The one-boundary Halphen description forces the
   index-two `I1/II` case and produces its unique exceptional curve
   `E=f-b`.  Thus

   ```text
   all (-1)-curves  <---->  all primitive CM cusp rays.
   ```

4. Item 4 in the raw definition listed exactly these `f-b` facets.  Kondō
   also proves that the only negative curves on a Coble surface are the
   boundary, `(-2)`-curves, and `(-1)`-curves.  Curves of nonnegative square
   lie on the closed positive sheet and add no hyperbolic wall.

Therefore the complete negative-curve list is exactly the raw wall list, and

\[
                     \boxed{\mathbf P(\operatorname{Nef}S)=P}. \tag{6.1}
\]

No projected-hidden-root equality was assumed in `D_1_split` or used as a
compiler gate; (6.1) is the classification pushback.

## 7. Classification of the carrier-admitting locus

Apply Theorem A to an arbitrary `(P,h) in R_1_split`, equivalently to an
object of `D_1_split` on which `C_1_split` passes.  The resulting surface is
complex, one-boundary, and finite-Aut, so Kondō puts it into one of the two
classes above.  Root exhaustion and the cusp bijection then determine its
entire wall Gram.  Thus there are at most two underlying raw-polytope
isometry types.

The exhaustive certificates above show that both types actually pass, and
their facet counts differ.  Hence:

> **Theorem B (AG classification pushed back to raw polytopes).** Up to
> integral marked isometry, the compiler-positive locus `R_1_split` inside
> `D_1_split` has exactly two underlying Coxeter
> polytopes.  Both are finite-volume noncompact polytopes in `H10`, with 22
> and 30 facets respectively, and both have exact complex strong-Coble
> carriers constructed from their raw data.

After this theorem is proved, the two rows may be named:

| raw row | returned surface | status of abstract polytope | new content here |
|---|---|---|---|
| `H10/22` | one-boundary Type I | exact project chamber; worldwide literal priority still open | surface-free-at-input `P -> K3 -> Coble -> same P` certificate; new primitive rank-20 equivariant K3 certificate |
| `H10/30` | one-boundary Type II | known `D16-D5-H10-30` chamber | surface-free-at-input exact return certificate; new primitive rank-20 equivariant K3 certificate |

This is a directionally new certification, not a claim of independent
discovery.  The live program imports specialized arithmetic split-census
archives for these two calibrations, and those archives have geometric
Type-I/II lineage.  The generic theorem derives the finite pair alphabet
directly from the intersections in `P` via (2.1); the executable does not yet
accept an arbitrary wall matrix and synthesize that census end to end.

This is positive progress even though it creates no new abstract polytope:
it gives two new directionally correct `P -> surface -> same P` realizations
and an intrinsic iff criterion on an independently defined denominator with a
non-singleton carrier-admitting locus.

## 8. Hostile audit and exact scope

| required gate | result |
|---|---|
| surface-free denominator does not name the two answers | **pass**; the two types appear only after AG feedback |
| all split configurations, not selected ones | **pass**; `27,648/18` and `32,768/16` exhaustive censuses |
| finite compiler exhaustivity | **pass mathematically** for the `h`-decorated datum; finite pair alphabets, discriminant groups, binary reduction at anti-rank nine, finite genus/primitive-embedding algorithms below it, positive-definite enumeration, and every compatible branch retained |
| all even closures and primitive equivariant K3 embeddings | **pass**; the extra det-4 closure is carried through both graph glues rather than discarded |
| exact fixed lattice | **pass**; the det-4 branch is rejected by `M` Smith `1^10,2` and actual invariants `(11,9,1)` |
| all 24 roots simultaneously effective and irreducible | **pass**; one root-regular ample chamber and complete bounded-root theorem |
| actual branch is the raw `b` | **pass**; adjunction plus the rank-ten span forces it |
| no hidden downstairs negative wall | **pass by AG feedback**; `Vinberg2`, classification, and cusp bijection |
| branchwise exact same-`P` return | **pass**; every passing branch constructs the carrier and returns the original full nef chamber |
| intrinsic carrier iff on `D_1_split` | **pass**; a carrier exists exactly when the compiler has a passing branch |
| realized branch for every `P in D_1_split` | **not proved**; the two positive calibrations and the classification of `R_1_split` do not show `D_1_split=R_1_split` |
| AG pushback on the whole a priori class | **not proved**; the exact two-isometry-class conclusion applies to the compiler-positive locus `R_1_split`, not every object of `D_1_split` |

The rating is therefore **uniform non-singleton `L3` exact-return criterion
and classification of the compiler-positive locus; not strict `L4` under the
2026-08-09 strengthened five-part contract**.  In particular, finite
exhaustivity of possible carrier branches is not class-wide essential
surjectivity: no theorem here says that every `P in D_1_split` has a realized
branch.  The scope is also limited as follows:

- every member is **noncompact**, because primitive isotropic cusp rays and
  their lifts are part of the definition;
- the dimension is fixed at `H10`;
- the compiler exhausts the one-boundary complex K3/Coble quotient species,
  not every imaginable algebraic variety attached to the same polytope;
- this gives no new higher-dimensional bound and does not solve the
  characteristic-three `H14/H17` crystalline transition; and
- it proves no abstract-polytope nonexistence outside `R_1_split`.

Thus it closes a real two-input theorem requested by the programme, but it
does not close the programme's global “all carrier species / higher
dimensions / compact examples” goal.

## 9. Sources and replay

The algebraic-geometric inputs are Shigeyuki Kondō,
[*Coble surfaces with finite automorphism group*](https://arxiv.org/abs/2102.00591),
especially the main classification theorem, Propositions `CME10`,
`NefCone`, `FiniteIndex`, `Vinberg2`, `Halphen`, and Lemma
`fibration-isotropic`; Nikulin's classification of non-symplectic K3
involutions; the complex K3 period map and global Torelli theorem; and
Vinberg's finite-index criterion for hyperbolic reflection groups.

Replay:

```bash
PYTHONPATH=src:. python3 scripts/certify_coble_one_boundary_raw_k3_promotion.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_coble_one_boundary_raw_k3_promotion.py
```

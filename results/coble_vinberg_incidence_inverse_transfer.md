# A uniform Coble--Vinberg incidence inverse theorem, and the exact remaining L4 gate

**Status (2026-08-09): a non-singleton, surface-free, finite
`POLYTOPE -> incidence stack -> strong Coble -> exact same P` criterion is
proved.  It is a uniform L3 exact-return theorem, not yet the requested
strict L4 theorem, because the raw numerical hypotheses do not force the
incidence stack to be nonempty for every input.**

The theorem applies simultaneously to the characteristic-three Type-VI
`H14/60` and Type-MII `H17/48` inputs, without putting either diagram, either
characteristic, or any surface in the definition of the raw class.  It also
isolates the remaining step with no euphemism:

```text
proved:
raw P -> finite exhaustive algebraic incidence compiler
      -> [a branch exists iff a compatible strong-Coble carrier exists]
      -> every existing branch returns exactly the same P
      -> Kondo classifies the compiler-positive locus;

not proved:
the P-only numerical conditions force a branch to exist for every P in the
a-priori denominator.
```

Thus this closes the formerly open *hidden-wall* part of the general
blowdown compiler.  It does not close essential surjectivity on the whole raw
denominator.

The two exact calibrations are frozen by
[`certify_coble_vinberg_incidence_inverse_transfer.py`](../scripts/certify_coble_vinberg_incidence_inverse_transfer.py),
with
[`JSON`](../data/enumerations/coble_vinberg_incidence_inverse_transfer.json)
and
[`tests`](../tests/test_coble_vinberg_incidence_inverse_transfer.py).

## 1. The independent surface-free denominator

An input is a marked oriented integral Coxeter datum with a strict integral
interior vector `h`.  Begin with the raw arithmetic lattice and mirrors; do
not supply a surface, boundary marking, or divisor labels.  Exhaustively
enumerate every admissible boundary de-scaling.  A branch consists of an odd
unimodular lattice

\[
             N\simeq I_{1,9+b},\qquad b\geq1,
\]

and pairwise orthogonal primitive classes `B_1,...,B_b` satisfying

\[
 B_i^2=-4,\qquad
 K=-\frac12\sum_iB_i\in N,\qquad
 K^2=-b,
\]

with `K` characteristic.  Every raw wall must de-scale to one of the purely
numerical possibilities

```text
B_i^2=-4;
C^2=-2, K.C=0, C.B_i=0 for every i;
E^2=-1, K.E=-1, with boundary profile (1,1) or (2).
```

Put `h_i=B_i/2` in the extended Picard module and form the rank-ten
orthogonal lattice

\[
 \operatorname{CM}_{\rm raw}=\{x:x.B_i=0\text{ for every }i\}.
\]

The candidate root associated with a numerical `(-2)` wall is `C`; the
candidate root associated with a numerical `(-1)` wall meeting `B_i,B_j`
once is

\[
                 \alpha_E=2E+h_i+h_j.                 \tag{1.1}
\]

Let `Gamma_raw` be this finite set.  The branch is in the surface-free
denominator `D_CobVin` when all of the following exact arithmetic tests pass.

1. The classes in `Gamma_raw` are distinct roots of square `-2`, span
   `CM_raw` rationally, and define a nondegenerate acute chamber.
2. Their graph passes Vinberg's finite-index criterion: every connected
   parabolic subdiagram is contained in a maximal parabolic subdiagram of
   the required rank.  Equivalently, an exact cone/Vinberg certificate proves
   that `W(Gamma_raw)` has finite index in `O(CM_raw)`.
3. Every primitive ideal extreme ray `f` of this root chamber and every
   boundary half `h_i` is tested for `E=f-h_i in N`.  The single-boundary
   numerical `(-1)` walls are exactly all integral outputs of these tests.
4. Reassembling the boundaries, root walls, and cusp lifts gives literally
   the original oriented wall set, and the full raw cone passes a direct
   finite-volume test.

These conditions never mention Kondō Types I--VII, MI, or MII.  They do not
mention an effective curve, a K3 surface, a crystalline period, a
characteristic, or a nef cone.  They also do not impose
`b in {1,2,5,8}`.  The class is non-singleton before any geometric
construction: the archived `H14/60` and `H17/48` inputs are distinct members.

## 2. The finite all-contraction incidence compiler

For every numerical boundary branch, enumerate sequences of raw wall
classes using the update

\[
 D'.F'=D.F+(D.E)(F.E),\qquad K'=K-E,                  \tag{2.1}
\]

whenever the current image of `E` has square `-1`.  A branch stops at
`P2` or at a minimal Hirzebruch surface `F_e` (including every `e` produced
by the finite numerical branches).  The enumeration is finite because there
are finitely many raw classes and every contraction drops the Picard rank.

It is also exhaustive in the direction needed for an inverse theorem.
Suppose `S=S_0` is an exact strong-Coble carrier and choose any minimal-model
contraction

\[
 S_0\longrightarrow S_1\longrightarrow\cdots\longrightarrow S_t.
\]

If `E_j` is contracted on `S_j`, its strict transform `C_j` on `S_0` is a
smooth irreducible rational curve and

\[
 C_j^2=-1-\sum_{i<j}m_i^2<0.                           \tag{2.2}
\]

On a strong Coble surface every irreducible negative curve is a boundary,
a `(-2)`-curve, or a `(-1)`-curve.  Exact equality with the input nef cone
therefore puts `C_j` among the raw facets.  Explicitly, for an ample `A` and
an irreducible negative curve `C`,

\[
        A_C=A+\frac{A.C}{-C^2}C
\]

is timelike, has `A_C.C=0`, and pairs strictly positively with every other
irreducible curve because distinct curves have nonnegative intersection.
Thus the inequality of `C` is genuinely a facet rather than a redundant
negative class.  Hence the compiler enumerates at
least one contraction sequence of every exact carrier.  This avoids an
unbounded search through arbitrary exceptional bases.

Reverse a fixed sequence as a universal iterated point blow-up of `P2` or
`F_e`.  Relative Hilbert schemes parameterize one curve in every pushed-down
raw wall class.  Take the locally closed locus on which

```text
all listed curves are smooth, geometrically irreducible, and distinct;
their prescribed incidences and infinitely-near conditions hold;
h is ample;
the B_i are disjoint reduced smooth rational curves with sum B=-2K;
h0(-K)=0 and h0(-2K)=1.
```

There are finitely many sequences and finitely many wall classes, so this is
a finite union `C_CobVin(P,h)` of finite-type locally closed stacks over
`Spec Z[1/2]`.  Cohomology/base-change stratification turns the two `h0`
conditions into finite determinantal conditions.  Standard Hilbert and
elimination algorithms decide geometric nonemptiness over a fixed
characteristic and compute its constructible image in `Spec Z[1/2]`.

This is an actual algebraic object compiled from `P`, not a search through
an unbounded collection of surfaces.

## 3. Every geometric branch is a strong Coble carrier

A geometric point of the compiler is an iterated blow-up of a rational
minimal surface, hence a smooth rational surface `S`.  The imposed boundary
conditions give

\[
 |-K_S|=\varnothing,\qquad
 |-2K_S|=\{B_1+\cdots+B_b\}.                           \tag{3.1}
\]

Thus `S` is a strong Coble surface.  Its marked Picard lattice is the raw
`N`, and every numerical wall is represented by the simultaneously chosen
smooth irreducible curve.

There is a small irreducibility point for the roots (1.1).  It cannot merely
be cited from the one-way classification of effective irreducible roots.
Let `E` meet `B_i,B_j` once.  The support matrix of
`E,B_i,B_j` is

\[
 \begin{pmatrix}
 -1&1&1\\1&-4&0\\1&0&-4
 \end{pmatrix},
\]

which is negative definite.  The effective representative of
`alpha_E=2E+B_i/2+B_j/2` is unique: successively intersecting an alternative
integral double with `E,B_i,B_j` forces the same components and
multiplicities.  If an effective root `beta` could be subtracted, its support
would be contained in these three curves.  Orthogonality to `B_i,B_j` would
write

\[
       \beta=xE+\frac{x}{4}B_i+\frac{x}{4}B_j
             =\frac{x}{2}\alpha_E.
\]

The equation `beta^2=-2` and effectivity force `x=2`, so
`beta=alpha_E`.  Therefore every compiler-produced candidate in
`Gamma_raw` is an effective irreducible Coble--Mukai root.  A listed
`(-2)`-curve is irreducible for the simpler support reason.

## 4. Vinberg closes all hidden root walls uniformly

Let `Gamma` be the now-effective root set on `S`.  The raw Vinberg gate gives

\[
             [O({\rm CM}(S)):W(\Gamma)]<\infty.         \tag{4.1}
\]

The reflection group `W(S)` generated by *all* effective irreducible roots
contains `W(Gamma)`, so it too has finite index.  Dolgachev--Kondō's
finite-index theorem makes `Aut(S)` finite.  Kondō's Proposition 3.7
(`Vinberg2` in the source TeX) then applies exactly: a finite-index set of
mutually distinct effective irreducible roots on a finite-automorphism Coble
surface is the entire effective irreducible root set.  Hence

\[
              \Gamma=\{\text{all effective irreducible CM roots}\}. \tag{4.2}
\]

This single argument works in every characteristic different from two.  It
is not a Type-VI or Type-MII template lookup.

Equation (4.2) excludes every hidden `(-2)`-curve and every hidden
`(-1)`-curve meeting two different boundaries: the latter would yield the
effective irreducible root (1.1).

## 5. The `b>1` cusp argument, including primitivity

Let `E` be any remaining hidden `(-1)`-curve.  Adjunction and
`B=-2K` give

\[
                         B.E=2.                         \tag{5.1}
\]

Distinct irreducible curves have nonnegative intersection, so after the
`(1,1)` case eliminated above, `E` meets one boundary `B_i` with total
multiplicity two and is disjoint from every other boundary.  Put

\[
                         f=E+h_i.                       \tag{5.2}
\]

Then `f` is an integral vector of the extended CM lattice,

\[
 f^2=0,\qquad f.B_j=0,qquad 2f=2E+B_i.                \tag{5.3}
\]

It is nef: `2f` has intersection zero with `E` and `B_i`, and nonnegative
intersection with every other irreducible curve.

The word *primitive* needs proof.  Suppose `f=mg` with primitive nef
isotropic `g` and `m>=1`.  The primitive-isotropic/fibration theorem gives a
genus-one fibration for `g`; the integral fiber class is `epsilon g` with
`epsilon in {1,2}`, the two Halphen possibilities.  The connected vertical
divisor

\[
                         2E+B_i
\]

has intersection matrix with radical generated by the multiplicity vector
`(2,1)`, whose coefficients have gcd one.  Hence it is one primitive fiber,
not a nontrivial multiple of a fiber.  But its class is
`2mg=(2m/epsilon)(epsilon g)`.  Therefore `2m/epsilon=1`; since
`epsilon<=2`, this forces `m=1` and `epsilon=2`.

Thus `f` is a primitive nef isotropic CM vector.  A null vector is an extreme
ray of any cone in one Lorentz time sheet, so `f` is one of the finitely many
primitive ideal rays already enumerated in `D_CobVin`.  Since
`E=f-h_i in N`, the complete cusp-lift gate already lists `E`.  This
contradiction excludes every hidden `(-1)`-curve, including the case of two
distinct intersection points versus tangency on the same boundary.

Together with the unique anti-bicanonical boundary and (4.2), all negative
curves have now been exhausted.  Curves of nonnegative square lie on the
closed future sheet and add no further hyperbolic inequality.  Consequently

\[
                  \boxed{\mathbf P(\operatorname{Nef}S)=P}       \tag{5.4}
\]

under the original integral orientation and marking.

## 6. The exact inverse theorem

Combining the preceding sections gives the following statement.

> **Theorem (Coble--Vinberg incidence inverse criterion).** Let `k` be an
> algebraically closed field of characteristic different from two.  For
> every marked `(P,h) in D_CobVin`, there is an exact strong-Coble full-nef
> carrier over `k`, taking `h` to an ample class, if and only if one of the
> finite incidence stacks `C_CobVin(P,h)` has a geometric point over `k`.
> Every such point produces a carrier whose projectivized nef cone is
> literally the input `P`.

Necessity is the exhaustive contraction argument.  Sufficiency is the
strong-Coble construction, Vinberg/Kondō root exhaustion, and the complete
`b>1` cusp argument.  This is stronger than a one-polytope foldback and
stronger than a finite list of guessed surface templates: it is a terminating
carrier criterion on an independently defined, non-singleton raw class.

It is nevertheless not automatic existence for every raw input.  Define

\[
 R_{\rm CobVin}=\{(P,h)\in D_{\rm CobVin}:
                   C_{\rm CobVin}(P,h)\ne\varnothing\}. \tag{6.1}
\]

The theorem constructs and returns every member of `R_CobVin`.  It does not
prove `R_CobVin=D_CobVin` or give a new P-only subdenominator on which this
equality holds.

## 7. The K3/crystalline gate is bypassed, not assumed

For every geometric branch, the canonical double cover of `S` branched
along the reduced divisor `B in |-2K_S|` is a smooth K3 surface with a tame
anti-symplectic deck involution.  This is a conclusion of the rational
surface construction.  In characteristic three, the resulting K3 therefore
comes with its actual crystalline period and a period-preserving involution.

This matters for the high-dimensional endpoints.  An abstract rank-22
lattice isometry would still have to preserve an Ogus characteristic
subspace, but the incidence route never guesses such a subspace.  It first
constructs `S`; the K3 and its crystalline data then exist geometrically.
For Type MII the cover is the supersingular `Kum(E x E)` of Artin invariant
one.  For Type VI the cover is the supersingular Hessian K3 in
characteristic three.

Accordingly:

```text
standalone lattice-only crystalline compiler:       still not constructed;
crystalline existence for every accepted branch:    automatic from the cover;
crystalline Torelli as a barrier to exact return:    bypassed.
```

## 8. Algebraic-geometry pushback and the honest L-level

Every compiler-positive surface has finite automorphism group by Section 4.
Kondō's exhaustive classification therefore gives

\[
 b\in\{1,2,5,8\},\qquad
 \dim P=9+b\in\{10,11,14,17\}.                         \tag{8.1}
\]

In particular the compiler-positive locus has no members in `H12`, `H13`,
`H15`, `H16`, or above `H17`; its `H14` and `H17` members exist only in
characteristic three.  These restrictions were not placed in
`D_CobVin`.

But (8.1) is a theorem about the carrier/stack-positive locus
`R_CobVin`, not abstract nonexistence of raw members of `D_CobVin`.  This is
the precise strict-goal score:

Every geometric branch is **noncompact**.  If its rational minimal model is
a Hirzebruch surface, pull back a ruling; if it is `P2`, use the ruling on
the penultimate `F1` model.  This gives a nonzero primitive nef isotropic
class and hence an ideal ray of the full nef cone.  Finite volume is wholly
separate: it comes from the raw exact Vinberg/cone gate in Section 1, not
from finite sidedness and not merely from finding an isotropic ray.

| stage | result |
|---|---|
| `C` independent surface-free non-singleton class | **pass** |
| finite enumeration inside `Car` | **pass**; every compatible carrier occurs on an enumerated branch |
| `Car` under the strict contract | **fail**; no realized branch is forced for every `P in D_CobVin` |
| `Return` exact same `P` from every geometric branch | **pass** |
| `AG` finite Aut and Kondō classification | **pass after branch existence** |
| `Pushback` theorem on every a-priori raw input | **fail**; only `R_CobVin` is classified |

The rating is therefore **carrier-relative L3 exact-return criterion; strict
programme L2--L3**, not L4.  The exact missing L4 theorem is:

> Find a surface-free, non-singleton P-only subdenominator, not defined by
> the known answers or by stack nonemptiness, whose numerical/combinatorial
> conditions force at least one compiled incidence stack to have a geometric
> point.  Then (5.4) and (8.1) immediately upgrade it to L4.

## 9. The two high-dimensional positive branches

The theorem has two distinct certified positive inputs.

| raw input | compiled incidence | root/cusp exhaustion | exact return |
|---|---|---|---|
| `H14/60` | five plane branches; conic plus four tangent lines; residual equations force `p=3` | 20 CM roots; 57 cusps; exactly 35 integral single-boundary lifts | `27,845` finite + `1,027` ideal + `0` hyperideal rays |
| `H17/48` | unique eight-frame; `4x4` grid on `P1xP1`; all 24 permutation graphs iff `p=3` | 40 CM roots; 202 cusps; no integral lift among `202*8` tests | `190,808` finite + `5,244` ideal + `0` hyperideal rays |

These are known abstract Coxeter polytopes and known Coble surfaces.  The new
project theorem is the common surface-free finite carrier criterion and the
uniform exact-return proof; no worldwide discovery claim is made for either
polytope.

## 10. Sources and replay

The published ingredients are:

1. Vinberg's finite-index criterion for hyperbolic reflection groups.
2. Dolgachev--Kondō's theorem that finite index of `W(S)` implies finite
   `Aut(S)`.
3. Kondō, Proposition 3.7 (`Vinberg2` in the source), exhausting effective
   irreducible roots from a finite-index root stock on a finite-Aut Coble
   surface.
4. Kondō's primitive-isotropic/genus-one-fibration lemma and the Halphen
   index-one/index-two theorem.
5. Kondō's nine-class theorem for finite-automorphism Coble surfaces in
   characteristic different from two:
   [*Coble surfaces with finite automorphism group*](https://arxiv.org/abs/2102.00591).

Replay the frozen high-dimensional calibration:

```bash
PYTHONPATH=src:. python3 scripts/certify_coble_vinberg_incidence_inverse_transfer.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_coble_vinberg_incidence_inverse_transfer.py
```

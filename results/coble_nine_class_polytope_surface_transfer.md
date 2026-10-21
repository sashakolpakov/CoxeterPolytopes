# The nine-class Coble return program: from a raw Coxeter chamber to a classified surface and back

**Status (2026-08-09): surface-side bound proved; all nine exact full-NS
templates are complete, including the sharp `H17` endpoint; the broad raw
Coble compiler remains `L2--L3`, while the separate nine-chamber Alexeev
denominator has a carrier-relative `L3` exact carrier iff.**  The Type-VI
and Type-MII computations supply complete positive chambers in `H14` and
`H17`; Type MI supplies a complete `H11/206` chamber exactly equal to the
independently constructed Fermat `C2-H11/206` chamber; the
characteristic-five Types VI and VII additionally supply
surface-first `H10/78` and `H10/72` full-nef chambers.  They do not by themselves
prove that every member of an independently defined raw polytope class has a
Coble carrier.  The nine-template recognizer is therefore kept separate from
the missing carrier-universality theorem.

Type VI now has a fully independent raw-polytope recognizer: from only its
even `(15,7,1)` lattice and sixty primitive mirrors it exhausts all `3,003`
candidate five-boundary frames, finds a unique integral Picard de-scaling,
reconstructs the complete divisor/CM/cusp package, and returns the identical
sixty oriented mirrors.  This proves a Type-VI carrier iff at L3; the final
`Gamma_VI`/Kondō gate still makes it a calibration rather than the
template-free L4 transfer theorem for a broader class.

Type MII now has a raw compiler of the same directional strength.  It
searches all `735,471` eight-subsets of the twenty-four short mirrors,
retains `21,211` orthogonal frames and one integral Picard de-scaling, then
forces the `4 x 4` grid on `P1 x P1` and all twenty-four permutation graphs.
Their incidence stack is nonempty exactly in characteristic three.  The grid
blow-up is directly strong Coble; the already separate root/cusp exhaustion
returns the identical `H17/48` chamber.  See the
[`blowdown/incidence compiler`](coble_polytope_first_blowdown_incidence_bridge.md).

There is nevertheless one noncircular finite-class theorem which must not be
lost in that broader caveat.  Alexeev's nine complete `r+a=22` chambers form
an independently defined surface-free class.  With the canonical Coble
long-wall decoration, exact strong-Coble full-nef return holds precisely for
the `(15,7,1)` and `(18,4,0)` rows, both in characteristic three.  Kondō's
boundary-count classification excludes ranks `13,14,16,20`, while the raw
split-root Gram theorem excludes the actual complete rank-11 and rank-12
chambers; the two full-NS constructions close the converse.  Under the
strengthened five-part contract this is carrier-relative `L3`, not `L4`:
the seven carrier-negative inputs have no realized Coble branch and no
jointly exhaustive alternative carrier.  It is not an abstract-polytope no-go
or essential surjectivity for the larger numerical screening class.  See the
[`exact nine-chamber iff`](alexeev_g0_exact_strong_coble_carrier_iff.md).

A second, structurally stronger but lower-dimensional theorem is now frozen
for the complex one-boundary species.  The independent surface-free
denominator `D_1_split` is defined by raw `I_(1,10)` wall, characteristic
facet, `E10` root-chamber, cusp-lift, and interior-vector data.  Its finite
mathematical compiler exhausts split-root assignments, saturations, every
even overclosure, primitive equivariant K3 graph glues, fixed lattices,
period chambers, and bounded irreducibility.  Exact same-polytope
one-boundary strong-Coble return is equivalent to compiler nonemptiness, and
Kondō then classifies the positive locus as exactly the `H10/22` and
`H10/30` chambers.  The extra determinant-`4` overclosure in the latter row
really is K3-embeddable, but has fixed type `(11,9,1)` and eighteen
anti-invariant degree-zero roots, so it fails at the correct geometric gates.
This is a carrier-relative `L3` exact-return criterion, with a mathematical
generic compiler and two executable calibrations.  Class-wide branch
nonemptiness is unproved, so it is not strict `L4` and not the still-missing
global characteristic-three or compact all-carrier theorem.  See the
[`one-boundary raw promotion theorem`](coble_one_boundary_raw_k3_promotion.md).

## 1. The algebraic classification already gives a sharp carrier bound

Let `S` be a strong Coble surface over an algebraically closed field of
characteristic different from two:

```text
|-K_S| is empty,
|-2K_S| = {B_1+...+B_b},
B_i are pairwise disjoint smooth rational curves.
```

Adjunction gives `B_i^2=-4` and

\[
 K_S^2={1\over4}\Big(\sum_iB_i\Big)^2=-b.
\]

Since `S` is rational, Noether's formula gives

\[
 \rho(S)=10-K_S^2=10+b,\qquad
 \dim \mathbf P(\operatorname{Nef}S)=9+b.                 \tag{1.1}
\]

Kondō classifies **all** Coble surfaces with finite automorphism group in
characteristic different from two into exactly nine isomorphism classes.
Their boundary counts are `1,2,1,2,1,5,1,2,8`.  Consequently every full-nef
finite-automorphism Coble carrier satisfies

\[
 \boxed{\dim P\le17},                                    \tag{1.2}
\]

and the unique top boundary-count class is Type MII in characteristic three.
This is AG-native information: neither `b<=8` in the finite-automorphism
locus nor the nine-row list is a consequence of the input Coxeter matrix.

The finite-automorphism hypothesis is forced once the proposed carrier is an
exact finite-sided finite-volume full-nef Coxeter chamber.  Every surface
automorphism preserves the nef cone.  Its image in the integral orthogonal
group therefore lies in the setwise stabilizer of a finite-volume
finite-sided hyperbolic polytope; that stabilizer is finite because it acts
faithfully on the finite spanning facet set.  Kondō's natural map
`Aut(S)->O(CM(S))` has finite kernel, so `Aut(S)` is finite.  Hence (1.2)
applies to **every** strong Coble surface whose complete nef cone has the
polytope property under study, without adding finite automorphism as an
independent surface hypothesis.

All these full-nef chambers are necessarily noncompact.  A smooth rational
surface of Picard rank at least two has a pulled-back nef isotropic ruling
class, which gives an ideal boundary ray.  The Coble lane is therefore kept
separate from the compact determinant-line/eigencarrier lane.

The classified carrier rows are:

| type | characteristic | boundary count `b` | CM roots | full-NS hyperbolic dimension |
|---|---:|---:|---:|---:|
| I | any `!=2` | 1 | 12 | 10 |
| I | any `!=2` | 2 | 12 | 11 |
| II | any `!=2` | 1 | 12 | 10 |
| V | 3 | 2 | 20 | 11 |
| VI | 5 | 1 | 20 | 10 |
| VI | 3 | 5 | 20 | 14 |
| VII | 5 | 1 | 20 | 10 |
| MI | 3 | 2 | 40 | 11 |
| MII | 3 | 8 | 40 | 17 |

Here “CM roots” means Kondō's complete effective irreducible root system in
the rank-ten Coble--Mukai lattice.  It is not yet the full Picard wall count.

## 2. The intrinsic raw-polytope compiler

The input is a marked rational Coxeter wall system `(N,{v_a})` with an
integral Picard marking; no surface is supplied.  Here `N` is the candidate
Neron--Severi lattice.  Reflections in boundary vectors need not preserve
`N` itself.  Arithmetic reflection integrality is checked on the canonical
scaled extended lattice `2*<N,B_i/2>` constructed below.  The compiler
enumerates every subset of walls that can serve as a Coble boundary frame
and applies the following finite tests.

### C0. Rational-surface lattice and boundary frame

For a candidate frame `B_1,...,B_b`, require

```text
N is the odd unimodular lattice I_(1,9+b),
B_i^2=-4 and B_i.B_j=0 for i!=j,
K=-(B_1+...+B_b)/2 belongs to N,
K is characteristic and K^2=-b.
```

These conditions, including every possible half-integral saturation, are
decidable by Smith normal form.  They reconstruct the canonical class from
`P`; it is not supplied as surface metadata.

### C1. Adjunction alphabet

Every non-boundary facet intended as a negative curve must have square `-1`
or `-2` and satisfy

\[
                       K.v=-2-v^2.                       \tag{2.1}
\]

A square-minus-two wall is disjoint from all boundaries.  A square-minus-one
wall has total boundary intersection two, hence has exactly one of the two
patterns

```text
(1,1): it meets two different boundaries once;
(2):   it meets one boundary twice.
```

Any other pattern rejects the package.

### C2. Fold to the Coble--Mukai lattice

Form the finite saturation

\[
 \widetilde N=\langle N,B_1/2,\ldots,B_b/2\rangle,
 \qquad CM=B_1^\perp\cap\cdots\cap B_b^\perp\subset\widetilde N. \tag{2.2}
\]

It has rank ten in every surviving package.  Translate the full-NS walls to
CM data by the intrinsic formulas

\[
 \begin{array}{rcl}
 C^2=-2&\longmapsto&C,\\
 E^2=-1,\ E.B_i=E.B_j=1
   &\longmapsto&2E+B_i/2+B_j/2,\\
 E^2=-1,\ E.B_i=2
   &\longmapsto&f=E+B_i/2.
 \end{array}                                               \tag{2.3}
\]

The first two outputs have square `-2`; the last is primitive isotropic.
Compute the complete finite-volume chamber cut out by the root outputs and
enumerate all of its primitive ideal rays.  The isotropic outputs in (2.3)
must be exactly the integral boundary lifts permitted by those rays.  Thus
the compiler checks roots **and** all cusp-to-exceptional-curve lifts, rather
than recognizing only Kondō's codimension-`b` CM section.

### C3. Nine finite AG packages

For each surviving CM chamber, compare its labelled root Gram, the
`(-1)`/`(-2)` root marking, discriminant glue, and Kondō `R`-invariant with
the nine classified packages.  This is a finite exact comparison; graph
shape alone is insufficient.  The package fixes the allowable
characteristic and, by Kondō's uniqueness theorem, its surface moduli is
either empty or one isomorphism class.

### C4. Exact return

For every nonempty package, reconstruct in `N`:

1. all boundary `(-4)` curves;
2. all actual `(-2)` curves represented by CM roots;
3. every `(-1)` curve meeting two boundaries, obtained from a `(-1)` root;
4. every `(-1)` curve meeting one boundary twice, obtained from a primitive
   CM cusp and the integral glue.

Kondō proves that these are the only possible negative curves on a Coble
surface.  Hence this list cuts out the entire nef cone.  The normalized wall
Gram, oriented cone, Coxeter labels, extreme rays, and literal Vinberg
finite-volume criterion are then replayed independently.  A package is
accepted only when this recovered full nef chamber is exactly the input
`P`, not merely its CM face or its unlabelled nerve.

## 3. Target return theorem and the non-tautology gap

Let `C_Cob^raw` be the surface-free screening class of marked wall systems
that pass C0--C2, including literal Coxeter and finite-volume gates but
without a surface supplied at input.  Define `C_Cob` intrinsically as the
subclass whose folded labelled Gram, glue, cusp lifts, and `R`-invariant pass
one of the nine finite C3 templates.  This definition is a finite calculation
on `P`; it neither names nor assumes the existence of a surface.  The
completed compiler is to prove

```text
for every P in C_Cob^raw:

  P has a finite-Aut strong Coble full-nef carrier in char != 2
    <=>
  P belongs to C_Cob and C4 returns P.
```

The forward implication is exhaustive **among actual Coble carriers** because C0--C2 are forced by
adjunction and the definition of the Coble--Mukai lattice, while Kondō's
classification exhausts all finite-automorphism Coble surfaces.  The reverse
implication uses Kondō's existence/uniqueness plus the complete negative-curve
and cusp-lift reconstruction.  When all nine templates are certified, the
surface-free accepted sublocus is therefore finite, has dimension at most
`17` by (1.2), and contains the nonempty rows in the table.  Equivalently,
every member of the independently checkable, non-singleton class `C_Cob`
has an exact carrier and the nine AG strata exhaust that class.  This is the
exact recognition form

```text
raw P -> finite exhaustive Coble packages -> exact same-P return
      -> AG classification -> n<=17 and positive boundary strata.
```

There is a crucial logical limitation.  `C_Cob` is currently defined by
passing one of the nine classified templates.  Consequently the resulting
nine-type statement is a necessary-and-sufficient **carrier recognizer**, but
does not yet give flagship leverage on the larger, independently defined
screening class `C_Cob^raw`.  Defining the accepted class by the answer and
then observing that it has nine answers would be circular as an `L4`
classification theorem.

The missing bridge must instead prove one of the following without importing
the nine templates:

1. every `P in C_Cob^raw` has a nonempty exhaustive blowdown/incidence moduli
   package whose resulting surface is strong Coble and whose full nef cone is
   exactly `P`; or
2. a general K3 double-cover compiler, derived from `P`, exhausts all
   invariant and anti-invariant lift lattices, periods (including the
   supersingular cases), effectivity branches, and hidden-root alternatives,
   and its intrinsic nonemptiness conditions force a strong Coble quotient.

Only after such a bridge does Kondō's classification transfer `b<=8`,
`dim(P)<=17`, and the nine rows to **every** polytope in a class defined
before the classification.  A failed Coble package remains only a carrier
obstruction until then; it is not an abstract-polytope nonexistence proof.

## 4. Current exact status

The following full-NS templates are already proved and replayable:

| class | complete walls | full chamber | lineage status |
|---|---:|---|---|
| Type I, one boundary | 22 | finite-volume noncompact `H10` | literal abstract priority open |
| Type II, one boundary | 30 | finite-volume noncompact `H10` | known `D16-D5-H10/30` |
| Type I, two boundaries | 26 | finite-volume noncompact `H11` | literal abstract priority open |
| Type V, `p=3`, two boundaries | 36 | finite-volume noncompact `H11` | literal abstract priority open |
| Type VI, `p=5`, one boundary | 78 | finite-volume noncompact `H10` | literal abstract priority open |
| Type VII, `p=5`, one boundary | 72 | finite-volume noncompact `H10` | literal abstract priority open |
| Type MI, `p=3`, two boundaries | 206 | finite-volume noncompact `H11` | known project `F3-C2-H11/206` |
| Type VI, `p=3`, five boundaries | 60 | finite-volume noncompact `H14` | known Alexeev `(15,7,1)` |
| Type MII, `p=3`, eight boundaries | 48 | finite-volume noncompact `H17` | known Alexeev/Keum--Kondō `(18,4,0)` |

For characteristic-five Type VI, the one-boundary cusp bijection gives
`1*(-4)+20*(-2)+57*(-1)` and a 78-wall chamber with
`21,661+1,797+0` rays and 166,320 Vinberg-certified edges.  Its labelled
roots span `E10` integrally.  Type VII similarly gives
`1*(-4)+20*(-2)+51*(-1)`, a 72-wall chamber with
`18,761+1,731+0` rays and 143,565 certified edges.  Both are project-new
surface-first complete realizations, `N/A/L0` on the inverse track, with
known surfaces and known `20E`/`20C` CM seeds; worldwide literal priority is
open.

For characteristic-three Type MI, the exact curve split is
`2*(-4)+30*(-2)+174*(-1)`.  Its forty-root Coble--Mukai chamber has
`6,210+247+0` rays; exactly 82 cusps lift over both boundaries, producing
164 single-boundary exceptional curves in addition to the ten ordinary
exceptional curves.  The complete `H11/206` cone has `327,566+18,354+0`
rays and 2,688,957 literal Vinberg edges, with zero edge or affine-neighborhood
failures.  A 206-facet permutation and a determinant-one lattice isometry
identify it entrywise with the known project chamber `F3-C2-H11/206` on the
Fermat quartic.  This is a project-new complete Type-MI full-nef realization
of a known polytope; worldwide priority for the surface interpretation is
open.

For characteristic-three Type V, the four `(-1)`-root vertices are recovered
intrinsically as the unique four roots pairing evenly with the other nineteen.
The visible divisor lattice has determinant four and a unique integral
index-two overclosure `NS=I_(1,11)`.  Exactly seven of its 29 Coble--Mukai
cusps lift over both boundaries, giving the complete split
`2*(-4)+16*(-2)+18*(-1)=36`.  Exact double description gives
`6,356+432+0` rays.  Literal Vinberg audits 47,929 elliptic edges and 1,154
connected affine components, with zero failures.  This is a project-new
surface-first complete `H11/36` realization; worldwide literal priority is
open and its inverse score is `N/A/L0`.

For characteristic-three Type VI, the complete curve split is
`5*(-4)+10*(-2)+45*(-1)`.  Exact double description gives 27,845
finite vertices, 1,027 ideal vertices, and no hyperideal rays; literal
Vinberg audits 244,330 elliptic edges and 2,704 affine components with zero
failures.  For Type MII, the complete split is
`8*(-4)+16*(-1)+24*(-2)`.  All 202 primitive Coble--Mukai cusps and all 1,616
single-boundary lifts are exhausted, with no additional exceptional curve;
an integral marked isometry identifies the resulting chamber with the known
finite-volume `(18,4,0)` system.  Independently, a raw full-48-wall cdd audit
gives `190,808` finite rays, `5,244` ideal rays, and zero hyperideal rays on
one causal sheet, directly proving finite volume.  These two constructions establish the
positive surface-side endpoint.  The Type-VI row additionally has the
independent polytope-first exact-return theorem described above; neither row
yet supplies template-free essential surjectivity for the broader Coble
screening class.  The uniform `D_CobVin` incidence compiler now closes exact
return for every geometric branch, including the crystalline `H14/H17`
cases, but neither it nor the complex one-boundary `D_1_split` compiler forces
a realized branch for every a-priori raw input.  See the
[`uniform incidence inverse theorem`](coble_vinberg_incidence_inverse_transfer.md).

Completion requires:

1. extend the now-proved raw Type-VI boundary-frame enumeration to every
   remaining template and prove it enumerates every marking without importing
   a surface-first wall ordering;
2. freeze the nine template isometry/`R`-invariant comparison and exact-return
   archive; and
3. lineage every resulting full-NS chamber independently of Kondō's known
   rank-ten CM chamber; and
4. replace template-defined acceptance by the general blowdown/incidence or
   K3-cover universality bridge described above.  This fourth item, not
   completion of the now-closed wall tables, is the `L4` gate.

The characteristic-five one-boundary templates are documented in
[`coble_type_vi_char5_one_boundary_full_ns_h10.md`](coble_type_vi_char5_one_boundary_full_ns_h10.md)
and
[`coble_type_vii_char5_one_boundary_full_ns_h10.md`](coble_type_vii_char5_one_boundary_full_ns_h10.md).
The Type-MI completion and its exact Fermat same-polytope bridge are
documented in
[`coble_type_mi_full_ns_h11.md`](coble_type_mi_full_ns_h11.md).
The final Type-V two-boundary completion is documented in
[`coble_type_v_two_boundary_full_ns_h11.md`](coble_type_v_two_boundary_full_ns_h11.md).
The two high-dimensional positive templates are documented in
[`coble_type_vi_char3_five_boundary_full_ns_h14.md`](coble_type_vi_char3_five_boundary_full_ns_h14.md)
and its
[`raw inverse recognizer`](polytope_first_coble_type_vi_h14_recognizer.md),
and
[`coble_type_mii_full_ns_h17.md`](coble_type_mii_full_ns_h17.md).
The latter's independent volume/curve/isometry replay is the
[`Type-MII hostile audit`](coble_type_mii_h17_hostile_audit.md).

## 5. Primary source

Shigeyuki Kondō,
[*Coble surfaces with finite automorphism group*](https://arxiv.org/abs/2102.00591),
especially Theorem 1.1/Table 2, Proposition `CME10`, Lemmas `roots` and
`roots2`, Propositions `NefCone`, `Vinberg2`, and `Halphen`, Lemma
`fibration-isotropic`, Sections 5--7, and the uniqueness argument.  The main
theorem states that there are exactly nine isomorphism classes in
characteristic different from two.

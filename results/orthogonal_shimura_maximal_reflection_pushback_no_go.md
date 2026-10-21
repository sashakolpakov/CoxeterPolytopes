# Orthogonal Shimura carrier for Coxeter chambers: exact return and the maximality pushback obstruction

**Status (2026-08-12): exact universal real-framed marked carrier; negative
terminal result for the direct Shimura-boundedness route.**  This note
addresses the actual proposed implication

\[
 \mathsf {Cox}^{\max }_{\rm arith}
 \longrightarrow \mathsf {AGCat}
 \longrightarrow \text{finiteness of the original chambers}
 \tag{0.1}
\]

without importing arithmetic reflection-group finiteness as the second
arrow.

There is a genuine first arrow to a **real-framed arithmetic AG category**.
A marked admissible Lorentz lattice has a canonical determinant-line
extension of type `(2,n)`.  Its orthogonal Shimura stack, special divisors,
real descent datum, and a framed lift of a totally real fixed component
recover the original Coxeter chamber exactly.  The construction is
functorial for marked lattice isometries.  If `AGCat` is required to mean
bare complex schemes/stacks with no real-component or uniformizing-sheet
marking, the forgetful first arrow still exists but is not faithful and has
no exact inverse.

There is no known, or currently applicable, second arrow.  The obstruction
is structural rather than terminological:

1. the original reflection group fixes the added positive line and lies in
   `Res SO(1,n)`, a proper non-Hermitian subgroup of `Res SO(2,n)`;
2. it has infinite index in every full arithmetic group defining the
   Shimura variety;
3. its hyperbolic space is a totally real fixed locus, not an algebraic
   Shimura subvariety;
4. reflection maximality is a statement inside `SO(1,n)` and forces neither
   a pure reflective modular divisor nor a bounded automorphic slope on the
   ambient Shimura variety; and
5. in field degree `d>1`, restriction of scalars turns one `F`-reflection
   into a rational isometry with a `d`-dimensional `(-1)`-eigenspace, so the
   classical rational reflective-form finiteness theorems do not apply.

The sharp control is not merely absence of a theorem.  The lattice

\[
 S=U\oplus E_8(-1)\oplus E_7(-1)
\]

has a finite-volume reflection-maximal chamber in `H^16`, but no
automorphic correction in Wang's sense.  Thus the most direct proposed
bridge

```text
maximal real arithmetic reflection group
  => pure Borcherds/reflective automorphic correction
```

is false.

This is a no-go for the **present direct orthogonal-Shimura strategy**, not
a theorem that no algebro-geometric proof can ever exist.

## 1. The exact functor

### 1.1 Input category

An object is marked admissible Coxeter data

\[
 (F,\sigma _0,L,b;r_1,\ldots ,r_s;x_0),                       \tag{1.1}
\]

where:

- `F` is totally real and `sigma_0` is a distinguished real embedding;
- `L` is an integral projective `O_F`-lattice of rank `n+1`, of signature
  `(1,n)` at `sigma_0` and negative definite at every other real place;
- the primitive negative vectors `r_i` define integral reflections;
- the oriented inequalities `(x,r_i)>=0`, oriented by the interior ray
  `x_0`, cut out a finite-volume Coxeter chamber `P`; and
- morphisms are field and lattice isometries preserving all markings.

The marking is essential.  An unmarked real polytope does not canonically
choose an integral lattice, and a complex divisor arrangement does not
remember an oriented real chamber.

### 1.2 Determinant-line extension

Put `K_L=det_{O_F}L` and give it the bilinear form

\[
 b_K(v_1\wedge\cdots\wedge v_{n+1},
     w_1\wedge\cdots\wedge w_{n+1})
 =2(-1)^n\det(b(v_i,w_j)).                                    \tag{1.2}
\]

Then

\[
 \widehat L=L\perp K_L                                      \tag{1.3}
\]

has signature `(2,n)` at `sigma_0` and is negative definite at
every other real place.  Every wall reflection acts on `widehat L` by its
given action on `L` and trivially on `K_L`.  Both (1.2) and (1.3) commute
with marked lattice isometries; no auxiliary square class is chosen.

Let

\[
 G=\operatorname {Res}_{F/\mathbf Q}\operatorname {GSpin}(\widehat V),
 \qquad \widehat V=\widehat L\otimes F.                       \tag{1.4}
\]

Its distinguished real factor has the type-IV domain `D^+` of oriented
positive two-planes, of complex dimension `n`; all other factors are
compact.  With, for example, the full lattice stabilizer `Gamma`, the
quotient stack

\[
 \mathscr X_{\widehat L}=[\mathcal D^+/\Gamma]                \tag{1.5}
\]

is an arithmetic Deligne--Mumford stack.  Its coarse quotient is
quasi-projective and its Baily--Borel compactification is a normal projective
variety.  If `[F:Q]>1`, `widehat V` is anisotropic: an `F`-isotropic vector
would remain isotropic at a definite real place.  The quotient is then
projective without a cusp boundary.

Toroidal compactifications may be added after choosing admissible cone data;
that choice is not needed for the minimal functor and should not be called
canonical.

### 1.3 Algebraic divisors and the totally real return

Choose a positive generator `e` of `K_L` over `F`.  Normalize every positive
ray `[x]` in `L_{sigma_0}` by `x^2=e^2`.  In the projective model of the
type-IV domain define

\[
 j_e([x])=[e+i x].                                            \tag{1.6}
\]

The line in (1.6) is isotropic and positive, and this gives a totally real
embedding

\[
 j_e:\mathbb H(L_{\sigma _0})\hookrightarrow\mathcal D^+.
 \tag{1.7}
\]

Let `iota` act as `+1` on `L` and `-1` on `K_L`.  The map
`[z] -> [iota(bar z)]` is an antiholomorphic involution; (1.7) is a connected
component of its fixed locus.  This packages the real slice as algebraic
descent data plus a marked component, rather than pretending it is a complex
subvariety.

For a wall root `r`, let

\[
 \mathcal D_r=\{[z]:(z,r)=0\}.                                \tag{1.8}
\]

Its arithmetic orbit defines an algebraic special divisor (or, to retain a
particular sheet, a marked special-cycle morphism) on (1.5).  Directly from
`(e,r)=0`,

\[
 j_e^{-1}(\mathcal D_r)=r^\perp\cap\mathbb H(L_{\sigma _0}).  \tag{1.9}
\]

Let `AGCat_R^fr` denote real algebraic arithmetic stacks enhanced by their
analytic uniformization, a framed fixed-component lift, and marked
special-cycle sheets.  This is a standard kind of arithmetic-geometric
object with extra real/analytic marking; it is deliberately not identified
with the category of bare complex varieties.  The following object is
faithful in `AGCat_R^fr`:

\[
 \mathcal A(P)=
 (\mathscr X_{\widehat L},\lambda,\tau,
  \widetilde {j_e(\mathbb H^n)},\{\widetilde{\mathcal Z(r_i)}\},
  x_0,\text{sheet labels}).
 \tag{1.10}
\]

The tildes in (1.10) mean that the chosen sheets live in the framed
uniformizing domain, before the full arithmetic quotient identifies their
translates.  Intersect the marked special sheets with the marked real
component, orient them using `x_0`, and take the common halfspaces.  Equation
(1.9) returns `P`.  Pairwise root Gram entries return every finite Coxeter
label; on the reflection stack the same label appears as the dihedral
inertia at a generic two-wall stratum.  Thus

\[
 \mathcal R(\mathcal A(P))\simeq P                            \tag{1.11}
\]

in the marked groupoid.

The unenhanced complex variety `X_widehat L` alone is not faithful.  Even
the full arithmetic quotient can identify distinct wall sheets.  The real
structure, its framed component lift, the sheets, and the oriented
halfspaces are part of the exact carrier.  This is the precise categorical
cost of exact return.

## 2. Why the Coxeter quotient is not the Shimura quotient

Let

\[
 H=\operatorname {Res}_{F/\mathbf Q}SO(V)
 \subset \operatorname {Res}_{F/\mathbf Q}SO(\widehat V)      \tag{2.1}
\]

be the stabilizer of the positive line `F e`.  The extended Coxeter group
`W(P)` is contained in `H`.

### Theorem 2.1 -- non-Hermitian stabilizer and infinite index

For `n>=3`:

1. the unique noncompact real factor of `H` is `SO(1,n)`, which is not of
   Hermitian type;
2. `H` is a proper algebraic subgroup of the ambient type-IV group, with

   \[
   \dim G_{\rm ad}-\dim H=d(n+1),\qquad d=[F:\mathbf Q];       \tag{2.2}
   \]

3. `W(P)` is Zariski dense in `H` but has infinite index in every full
   arithmetic subgroup of the ambient type-IV group; and
4. the orbit (1.7) is totally real: at every point

   \[
   Tj_e(\mathbb H^n)\cap J Tj_e(\mathbb H^n)=0.               \tag{2.3}
   \]

Hence `W(P)\backslash D^+` is not the Baily--Borel arithmetic quotient in
(1.5), and the hyperbolic cycle is not an algebraic Shimura subvariety.

#### Proof

The classification of Hermitian symmetric spaces excludes `SO(1,n)` for
`n>=3`.  The two orthogonal groups in (2.1) have relative dimensions
`n(n+1)/2` and `(n+1)(n+2)/2`, giving (2.2) after restriction of scalars.
Since `W(P)` is a finite-covolume lattice in `SO(1,n)`, Borel density gives
its Zariski closure `H`.  A finite-index subgroup of an ambient arithmetic
group is Zariski dense in the ambient group, so properness of `H` rules out
finite index.  Finally, (1.7) is a fixed component of an antiholomorphic
involution, which proves (2.3).  `square`

For `n=2`, `SO(1,2)` is itself of Hermitian type and the hyperbolic quotient
has the separate modular-curve route.  The fixed-positive-line embedding
(1.7) is nevertheless still a totally real embedding in the larger type-IV
domain.  It does not repair the all-dimensional argument.

This theorem pinpoints the categorical compromise.  Quotient by the group
which actually has `P` as a chamber, and Baily--Borel algebraicity is lost.
Quotient by a full ambient arithmetic group, and one gets an algebraic
variety, but real reflection maximality is no longer a maximality statement
about that group.

## 3. Restriction of scalars does not recover the classical rational theory

One might try to replace the native Hilbert-orthogonal variety by the
rational trace quadratic space

\[
 T=\operatorname {Res}_{F/\mathbf Q}\widehat V
\]

with trace form.  Put `m=n+2` and `d=[F:Q]`.  Then

\[
 \operatorname {sig}(T)=(2,N),\qquad N=d(n+2)-2.              \tag{3.1}
\]

The full rational type-IV domain has complex dimension `N`, while the native
domain has dimension `n`, so its codimension is

\[
 N-n=(d-1)(n+2).                                              \tag{3.2}
\]

An `F`-reflection `s_r` has `(-1)`-eigenspace `F r`.  After forgetting `F`,
this eigenspace has rational dimension `d`.  It is negative definite for the
trace form, and therefore the fixed type-IV locus has codimension `d`, not
one.  Thus

\[
 s_r\text{ is a rational pseudoreflection}\iff d=1.           \tag{3.3}
\]

There is also an ambient-group index obstruction:

\[
 \dim SO(T)-\dim\operatorname {Res}_{F/\mathbf Q}SO(\widehat V)
 ={d(d-1)(n+2)^2\over2}>0                                    \tag{3.4}
\]

for `d>1`.  The native arithmetic group is consequently infinite index in
every full trace-orthogonal arithmetic group containing it.

Equations (3.2)--(3.4) rule out the tempting shortcut from a native facet
branch divisor to Ma's or Wang's classical rational reflective-divisor
classification.  A coincidental rational root may exist for a particular
input, but it is an additional root-by-root fact, not a functorial transfer
of the original reflection.

## 4. Reflection maximality does not force automorphic purity

A reflective modular form on (1.5) is an automorphic section whose divisor
is supported on global ramification mirrors of the **ambient** arithmetic
group.  Hyperbolic reflection maximality says that `W(P)` has no proper
discrete reflection overgroup inside `Isom(H^n)`, equivalently inside the
real subgroup `H` of (2.1).  These are statements about different groups.

The `W(P)`-orbit of its walls is not invariant under the ambient arithmetic
group.  Taking the ambient orbit makes an algebraic special divisor but can
add sheets not controlled by maximality.  It can also change the divisor
class.  Nothing in real reflection maximality proves

\[
 \sum_i a_i\mathcal Z(r_i)\sim w\lambda,                      \tag{4.1}
\]

with no other support, let alone bounds `a_i/w`, `w`, or
`w lambda^n`.

### Theorem 4.1 -- a maximal negative control

For

\[
 S=U\oplus E_8(-1)\oplus E_7(-1),                            \tag{4.2}
\]

the complete norm-minus-two reflection chamber is the finite-volume
`Gamma_16` chamber in `H^16`.  It is maximal among all discrete hyperbolic
reflection groups.  Nevertheless, Wang's Theorem 7.2 says that `S` has no
automorphic correction: there is no positive integer `m` for which
`U(m)\oplus S` carries the required `2`-reflective modular form.

Therefore

\[
 \boxed{\text{reflection maximality}
 \not\Rightarrow\text{automorphic correction}.}              \tag{4.3}
\]

The exact maximality and same-lattice identification are recorded in
[`wang_h16_maximal_reflection_no_correction.md`](wang_h16_maximal_reflection_no_correction.md).
Statement (4.3) concerns Wang's standard automorphic-correction extension;
it does not assert nonexistence of every modular form on every imaginable
extension.  It is exactly enough to invalidate the proposed universal
maximality-to-correction lemma.

## 5. What the available AG theorems actually assume

| theorem | valid input | output | missing implication from `P` |
|---|---|---|---|
| Baily--Borel | a full arithmetic quotient of a Hermitian symmetric domain | canonical projective algebraic compactification | `W(P)` is infinite index in the full type-IV group |
| Bruinier's Hilbert-orthogonal theta lift | a native `(2,n)` datum over a totally real field and principal-part data satisfying the cusp obstruction | meromorphic automorphic forms and special divisors | no automatic pure support, bounded slope, or finiteness of carrier lattices |
| Gritsenko--Hulek | a strongly reflective form already supplied | canonical-class and uniruledness consequences at specified weight | no existence or weight bound from reflection maximality |
| Ma | integral **rational** lattices of signature `(2,N)`; reflective form of bounded vanishing/slope | finiteness up to scaling in the stated high-dimensional ranges | fails (3.3) for `d>1`; purity and bounded slope are not forced even for `d=1` |
| Wang | classical rational orthogonal lattices, often with split `U`/`2U` hypotheses | nonexistence, finiteness, and explicit reflective classifications | the anisotropic native Hilbert carrier need not split; (4.3) kills automatic correction |
| Maeda--Ohara (2025) | unitary lattices over imaginary quadratic fields of odd discriminant; a reflective form with slope at most fixed `r` | finitely many such lattices up to scaling | not the general CM-over-totally-real carrier, and again assumes the bounded reflective form |

Thus no theorem in this list has the premise

```text
maximal arithmetic real-hyperbolic reflection group
```

and the conclusion

```text
bounded pure reflective divisor, bounded carrier, or finite carrier lattice.
```

General-type theorems do not fill this gap.  A neat compact determinant-line
quotient already has

\[
 K_X\equiv n\lambda                                             \tag{5.1}
\]

with `lambda` ample, independently of whether the original Lorentz lattice
is reflective.  General type is therefore an output shared by the positive
and negative inputs, not a contradiction to a finite-volume totally real
cycle.

## 6. Why compensated divisors, surfaces, and BMY do not push back

Suppose a selected target divisor `T` contains the facet sheets.  On a
projective carrier, ampleness alone gives, for all sufficiently large `k`,

\[
 T+C_k\sim k\lambda,\qquad C_k\ge0.                            \tag{6.1}
\]

Hilbert-orthogonal Borcherds theory can often make (6.1) much more explicit
with special-divisor compensation chosen to miss the marked real slice.
That preserves the exact return (1.11), but it does not create boundedness:
the level, support, matrix denominators, weight, top degree, and arithmetic
volume may all grow with `P`.  Unbounded compensation is the automorphic
version of Serre's tautology (6.1).

On an exact sheet-separating neat cover, an effective divisor
`D equiv w lambda` gives only

\[
 \#\operatorname {Facets}(P)
 \le {w\lambda^n\over q_{\min}},                              \tag{6.2}
\]

where `q_min` is the least target multiplicity.  Equation (6.2) becomes a
polytope theorem only after a noncircular, class-uniform upper bound for its
right side.

Cutting by a canonical complete-intersection surface does not improve
(6.2).  Every surface intersection contributes the common factor
`k^(n-2)`, which cancels.  Ordinary BMY on that surface is a strict positive
proportionality identity involving the ambient volume but no facets; with
the divisor arrangement it reproduces (6.2).  On the quotient stack the
reflection inertia is already absorbed into the orbifold Chern classes, so
it supplies no separately signed penalty per Coxeter label.  The exact
calculation is in
[`compact_shimura_ci_bmy_feedback_no_go.md`](compact_shimura_ci_bmy_feedback_no_go.md).

Nor can general boundedness of pairs be invoked without its hypotheses.
KSBA/HMX boundedness fixes dimension and volume (or a finite coefficient set
plus the relevant numerical data); BAB fixes dimension for log Fano pairs.
Here a dimension bound is one of the desired conclusions, `lambda^n` is an
unbounded arithmetic covolume, and the natural compact neat carrier is of
general type rather than log Fano.

## 7. A unitary-ball variant removes one defect, not the pushback defect

There is an alternative universal AG carrier that does not add a positive
line.  Let `E=F(i)` and put

\[
 V_E=V\otimes_F E,\qquad
 h(v\otimes a,w\otimes b)=a\bar b\,b(v,w).                    \tag{7.1}
\]

At `sigma_0`, `h` has signature `(1,n)`; at the other real places it is
definite.  The unitary Shimura domain is a complex `n`-ball.  Complex
conjugation has `H^n` as a totally real fixed component, and every real wall
reflection extends to a genuine order-two complex pseudoreflection whose
fixed set is a ball divisor.  Thus (7.1), with the real component and marked
divisor sheets retained, is another exact functor to AG data.

This avoids the determinant line and the rational trace codimension-`d`
phenomenon **on the native unitary variety**.  It does not provide the second
arrow in (0.1).  The real orthogonal group is a proper subgroup of the
unitary group, and `W(P)` is again infinite index in the full arithmetic
unitary group.  Real reflection maximality forces neither unitary reflection
generation nor a reflective unitary modular form.  Maeda--Ohara's theorem
then cannot be applied: it assumes such a form with bounded slope and, in
its stated finiteness theorem, works over imaginary quadratic fields rather
than arbitrary CM extensions `E/F`.

The unitary variant is therefore a cleaner carrier option, not hidden AG
leverage.

## 8. Minimal honest repair

The direct Shimura program would become a proof of Coxeter finiteness only
after two new, independent results.

### A. A maximality-to-bounded-divisor theorem

For every reflection-maximal input, produce either:

1. a pure reflective automorphic divisor on a canonical native carrier,
   with slope and top degree bounded independently of the number of facets,
   and with exact no-hidden-wall return; or
2. an exhaustive alternative AG carrier branch with the same bounded data.

Theorem 4.1 means this cannot simply be Wang's automorphic correction for
all inputs.  Arbitrary slice-disjoint compensation is insufficient by
Section 6.

### B. A native Hilbert-reflective finiteness theorem

Prove, for

\[
 \operatorname {Res}_{F/\mathbf Q}SO(2,n),                    \tag{8.1}
\]

that the bounded divisors from A occur on only finitely many fields,
dimensions, and lattice similarity classes.  The theorem must be native
over arbitrary totally real `F`; equations (3.3)--(3.4) prevent automatic
reduction to the classical rational list.

After A and B one must still enumerate the finite determinant/gluing data
and run the exact real-slice hidden-mirror and finite-volume return tests.

The final ledger is therefore

```text
marked P -> algebraic Shimura carrier:          PASS
exact real-slice/divisor return:                PASS
purely complex AG object alone faithful:       FAIL
Coxeter group itself defines Shimura quotient: FAIL for n>=3
reflection maximality -> pure correction:      FALSE in the standard route
trace reduction to rational reflections:       FAIL for [F:Q]>1
compensated effective carrier:                  PASS, but unbounded
general type / ordinary BMY pushback:           FAIL
native Hilbert bounded-slope finiteness:        MISSING
new finiteness theorem for original P:          NOT OBTAINED by this route
```

## 9. Primary sources and exact project controls

1. J. H. Bruinier,
   [*Regularized theta lifts for orthogonal groups over totally real fields*](https://arxiv.org/abs/0908.3076):
   native Hilbert-orthogonal Shimura varieties, special divisors, and the
   Whittaker/Borcherds obstruction.
2. S. Ma,
   [*On the Kodaira dimension of orthogonal modular varieties*](https://arxiv.org/abs/1701.03225):
   classical rational orthogonal scope and finiteness for reflective forms
   of bounded vanishing order.
3. H. Wang,
   [*The classification of 2-reflective modular forms*](https://arxiv.org/abs/1906.10459),
   especially Theorem 7.2: the hyperbolically reflective lattices without
   automorphic correction, including (4.2).
4. V. Gritsenko and K. Hulek,
   [*Uniruledness of orthogonal modular varieties*](https://arxiv.org/abs/1202.3361):
   consequences after a strongly reflective form is supplied.
5. Y. Maeda and K. Ohara,
   [*Finiteness of Free Algebras of Modular Forms on Unitary Groups*](https://arxiv.org/abs/2505.13698),
   Theorem 2.3.1: unitary bounded-slope reflective finiteness in its stated
   imaginary-quadratic scope.
6. The complete determinant-line construction and exact compensation scope
   are in
   [`polytope_to_automorphic_divisor_transfer_theorem.md`](polytope_to_automorphic_divisor_transfer_theorem.md).
7. The trace and native-Hilbert obstruction is independently developed in
   [`automorphic_pushback_l4_hostile_audit.md`](automorphic_pushback_l4_hostile_audit.md).

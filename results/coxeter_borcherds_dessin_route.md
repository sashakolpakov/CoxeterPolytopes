# Coxeter--Borcherds dessins: universal carrier and rigid classification locus

**Status (2026-08-10): proved synthesis of the existing determinant-line,
Whittaker-obstruction, safe-compensation, field-split, and cusp-transfer
theorems, now complemented by the literal Gram--Belyi surface theorem.**  It
gives general arithmetic-polytope-to-algebraic-object constructions and exact
marked return.  It also identifies the extra absolute-degree rigidity needed
for unconditional classification pushback.

## 1. Three levels of dessin

Let \(P\subset\mathbb H^n\) be a marked compact or finite-volume arithmetic
Coxeter polytope, supplied by an even admissible integral realization
\((L,q,\{r_i\})\) over its Vinberg field \(k\).
Vinberg's reflection lemma is the reason this quadratic/simplest-type input is
the relevant general arithmetic lane: an arithmetic hyperbolic lattice which
contains real hyperplane reflections must be of simplest type.  The marking
and even integral normalization are still explicit input choices; the theorem
does not claim that an unmarked real polytope canonically chooses one of them.

1. The **arrangement dessin** is the determinant-line orthogonal Shimura
   stack, its marked Heegner mirrors, the integral root/reflection local
   system, the chosen real chamber, and its rational-isotropic cusp marking.
2. A **compensated Borcherds dessin** also carries a Whittaker/Borcherds
   divisor
   
   \[
       D_{\rm fac}+D_{\rm safe},
   \]
   
   where \(D_{\rm fac}\) contains the separated facet sheets and every
   component of \(D_{\rm safe}\) is disjoint from the complete canonical real
   slice.
3. A **rigid Borcherds dessin** has no compensation, or has compensation with
   a facet-count-independent bound on its normalized automorphic slope
   \(w\lambda^n/q_{\min}\).  This is the analogue of the three-branch-value
   rigidity in classical dessins.  Unbounded compensation is analogous to
   allowing arbitrary branch data and is not classifying.

The distinction is essential: the first two levels are general; the third is
selective and is the only presently credible route to a polytope bound.

## 2. Universal compensated-dessin theorem

> **Theorem.**  For every marked \(P\) as above, the direct determinant-line
> extension \(\widehat L=L\perp K_L\) gives a type-IV arithmetic stack and a
> marked arrangement dessin whose canonical real slice contains exactly the
> original hyperbolic space.  After a computable congruence level, its facet
> sheets are separated.  The marked roots, reflections, dihedral inertia,
> oriented chamber, and rational-isotropic marking recover the normalized
> Gram datum, finite Coxeter labels, ordinary face poset, and ideal
> vertices/cusps of \(P\).
>
> At a direct level \(N\widehat L\), Fourier columns with nonzero
> auxiliary-line discriminant coordinate span the full dual cusp obstruction
> space.  Therefore every finite facet principal part admits a finite complex
> Whittaker correction supported on special divisors disjoint from the whole
> real slice.  If the Weil representation has a real \(q\)-structure, the
> safe columns positively span.  With an integral \(q\)-structure, clearing
> denominators produces an effective integral compensated Borcherds divisor
> on the open orthogonal stack after multiplying the target multiplicity.

The integral statement is unconditional over \(k=\mathbb Q\) by McGraw's
integral-basis theorem.  For a general totally real field it currently uses
the Hilbert--Weil integral-basis adaptation stated, but not written out in
full, by Bruinier; the complex Whittaker/Green statement does not need that
upgrade.
For \(k=\mathbb Q\), orders on the toroidal/Baily--Borel boundary are a
separate compactification gate and may include poles; “effective” here refers
to the interior special-divisor coefficients.  For the anisotropic
degree-greater-than-one lane there is no cusp boundary.

For \(k\ne\mathbb Q\), the definite conjugate place forces anisotropy, so the
arithmetic quotient and a finite-covolume reflection chamber are compact.
Finite-volume noncompact inputs therefore lie over \(\mathbb Q\).  At their
cusps, the corrected Witt-coordinate theorem splits candidate mirrors into:

- nonvertical mirrors, whose penetration height is explicitly bounded; and
- vertical mirrors, which form finite positive-definite shells modulo the
  parabolic translation lattice.

The literal Vinberg full-parabolic-neighbourhood test remains part of each
input's exact return.  Only ideal vertices/cusps lie on the sphere at
infinity; ordinary faces may merely have ideal vertices in their closures.

## 3. Exact obstruction criterion

For a finite principal part \(c\), let \(\operatorname{ob}(c)\) be its pairing
with the dual vector-valued cusp space.  Bruinier's theorem gives

\[
  c\text{ is a weakly holomorphic Whittaker/Borcherds principal part}
  \quad\Longleftrightarrow\quad
  \operatorname{ob}(c)=0.
\]

This is a genuine necessary-and-sufficient algebraic criterion for the
chosen carrier species.  It is not a converse for every meromorphic
automorphic form in the compact anisotropic Hilbert setting.  It also does
not say that the pure facet part passes: safe compensation may be necessary.

Thus the proposed classification invariant attached to \(P\) is not merely
the existence of a Shimura variety.  It is the finite package

\[
 (k,L,A_L,N_{\rm sep},c_{\rm fac},\operatorname{ob},
   c_{\rm safe}^{\min},w,\lambda^n,\partial_P),
\]

including the least safe-support/weight problem and exact cusp return.

## 4. Why the rigid locus is selective

Two rational norm-two controls prove that exact determinant-line
reflectivity is neither empty nor automatic.

- For \(S=II_{1,9}=U+E_8(-1)\), the determinant extension is isometric to
  \(2U+E_7(-1)\), which is completely \(2\)-reflective.  Its weight-165 form
  gives an exact complete-arrangement return of the known finite-volume
  \(H^9\) chamber.
- For \(S=U+E_8(-1)+D_4(-1)\), the finite-volume \(H^{13}\) norm-two chamber
  exists, but its determinant extension \(S+\langle8\rangle\) has signature
  \((2,13)\).  Wang's classification contains no \(2\)-reflective lattice in
  that signature.  Hence this determinant dessin has no globally pure
  norm-two reflective form.

For compact Bugaenko \(H^8/11\), the pure standard principal part and the
complete first compensator shell are obstructed.  Some higher safe
compensation exists at the stated arithmetic level, but its least support
and slope are not yet known.

These examples give the first useful trichotomy:

```text
pure rigid dessin              -> reflective classification can act;
bounded compensated dessin     -> a facet/degree theorem may act;
unbounded compensated dessin   -> universal but AG-nonrestrictive.
```

Indeed, on any projective carrier an arbitrary effective target divisor can
be completed to a high multiple of an ample divisor.  Therefore mere
existence of uncontrolled compensation cannot restrict \(P\).

## 5. Finite reconstruction from a height bound

The finite permutation analogy becomes exact once the arithmetic geometry
controls coefficient height.

> **Height-to-level lemma.**  Let \(k\) have degree \(d\), and let
> \[
>   \mathcal O_k(B)=
>   \{a\in\mathcal O_k:|\sigma(a)|\le B
>       \text{ for every embedding }\sigma:k\hookrightarrow\mathbb C\}.
> \]
> If a nonzero ideal \(\mathfrak a\) satisfies
> \(N\mathfrak a>(2B)^d\), then reduction
> \[
>       \mathcal O_k(B)\longrightarrow\mathcal O_k/\mathfrak a
> \]
> is injective.

Indeed, if distinct \(a,b\) have the same residue, then
\((a-b)\subseteq\mathfrak a\), whence
\[
 N\mathfrak a\le |N_{k/\mathbb Q}(a-b)|\le(2B)^d,
\]
a contradiction.  For a fixed denominator \(D\), apply the lemma to the
integral entries \(Dg_{ij}\).

Consequently, fix the field/order, root norms and denominators, and suppose a
P-only class has a uniform all-embedding bound \(B\) for its integral Gram
entries.  At one good level of norm \(>(2B)^d\), the marked finite quadratic
module, the reduced root lines/reflections, and their pairings determine the
exact Gram matrix uniquely.  The oriented halfspaces and literal Vinberg
test then recover the polytope.  Exact finite orders of the finitely many
elliptic products can simultaneously be preserved by excluding the finitely
many prime ideals dividing the nonzero matrices
\((s_i s_j)^e-I\) for proper divisors \(e\mid m_{ij}\).

The finite-order assertion and its stack interpretation are now proved in
the [congruence-dessin theorem](coxeter_congruence_dessin_theorem.md).  It
also proves that a fixed prime-power tower separates any prescribed finite
list of nonidentity matrix words.  In the \(E_{10}\) calibration, level five
preserves all ten marked involutions, all 45 finite pair orders, and recovers
the exact bounded Gram matrix.  This is finite marked monodromy in the
literal dessins sense, not an analogy imposed after the fact.

No one finite quotient proves that an infinite-order product is infinite.
That information remains in the exact lifted Gram, or is separated through
the congruence tower.  Thus the missing bridge to a finite dessin
classification is now sharply located:

\[
\boxed{\text{bound the Gram height from the Borcherds/Hilbert complexity}.}
\]

Once that bound is available, the continuous-looking number-field problem
reduces to a finite marked quadratic/reflection enumeration at a computable
level.

For a fixed root-norm/divisibility alphabet, this is principally a **dotted
height** problem.  Intersecting Coxeter pairs have bounded normalized
pairing, parallel pairs have the fixed limiting value, and the definite
conjugate forms bound their pairings by Cauchy--Schwarz.  The unbounded
entries at the distinguished real place are the ultraparallel/dotted
pairings \(\cosh d(H_i,H_j)\).  A successful Borcherds-degree theorem need
therefore control those distances (and the finite denominator/root-norm
alphabet), rather than all Gram data independently.

There is now one exact way to supply that bound.  Choose a maximal arithmetic
reflection overgroup \(W(P)\leq W_0=W(Q)\) and set

\[
 \deg_{\rm Cox}(P)=min_{W_0}[W_0:W(P)].
\]

The [bounded-degree Schreier--Coxeter dessin theorem](bounded_degree_schreier_coxeter_dessin_theorem.md)
proves that the finite coset action, together with a reflection certificate
and the marked base carrier, recovers \(P\).  For a fixed integral base and
\(I=[W_0:W(P)]\), every facet Gram entry has all-embedding height at most
\(C_QM_Q^{2I-2}\).  Hence bounded Coxeter degree gives a computable faithful
congruence level and a finite exact census.  Combined with the finiteness of
maximal arithmetic reflection groups, bounded absolute covolume gives global
finiteness, although the required maximal-base catalogue is still
existential rather than compiled.

The degree parameter is indispensable.  An exact rational ideal-polygon
family keeps \(n=2\), \(k=\mathbb Q\), the lattice \(I_{2,1}\), and the
primitive norm/divisibility alphabet fixed, but has index and area \(N\pi\)
and dotted entry \(1-2(N-1)^2\).  Its wall roots also collide modulo every
fixed principal level.  Thus dimension, field, discriminant, compactness
type, and root alphabet alone cannot support the requested classification.
The remaining automorphic bridge is the concrete implication

\[
\boxed{\text{bounded marked Borcherds complexity}
        \Longrightarrow \deg_{\rm Cox}(P)\text{ bounded}.}
\]

There is now a direct route which bypasses both the maximal carrier and a
Gram-height lift.  The real-marked
[Gram--Belyi surface](gram_belyi_coxeter_dessin_theorem.md) stores each
scaling-invariant squared Gram entry as a marked point on a minimal normalized
Belyi-polynomial leaf, joins all leaves to one Belyi backbone, and takes the
product with `P1`.  Its natural polarization satisfies

\[
A_P^2=2\left(\binom f2+\beta_{\rm GB}(P)\right),\qquad
\beta_{\rm GB}(P)\ge f(f-1).
\]

The real-place, acute-sign and time-sheet markings return the exact geometric
unit-normal Gram; the construction does not pretend to recover primitive
integral root lengths.  Bounded absolute \(A_P^2\) gives a finite effective
passport, chamber and Vinberg-field list without ABSW, Prasad, or Mostow.
The precise missing bridge is now equally concrete:

\[
\boxed{\text{maximality or intrinsic bounded geometry}
       \Longrightarrow A_P^2\text{ bounded}.}
\]

The prime-index `H10` redoubling calculation proves that no such bound can
hold for unrestricted finite-index descendants with fixed ambient arithmetic
data.  The [hostile audit](gram_belyi_coxeter_dessin_theorem_hostile_audit.md)
checks the gluing, intersection formula, real return and explicit field bound.

The exact low-label calibration gives packet degrees
\(5,3,\ge29,7\) for \(m=3,4,5,6\), and \(2\) for a parallel pair.
Thus the local arithmetic cost can be much larger than the Coxeter label:
the \(m=5\) lower bound comes from the prime \(29\) in
\(\lambda=(7+2\sqrt5)/29\).  This weighted edge statistic is now the first
finite experimental input to the maximality-versus-\(A_P^2\) problem.

The initial low-vertex degree-\(29\) sharpness reduction is exact but
incomplete.  After complementing so \(B(\lambda)=0\), it proved at that stage
that a witness needs at least **four** vertices in the lambda-colored fibre.
The frozen marked passport denominator is `42,610,678`: `3,977` records fell
in the first low-vertex pass, and the
follow-up eliminates all `1,372` three-black/`WW` records after collapsing
them to `784` coefficient systems.  Thus `5,349` marked records are excluded
and exactly `42,605,329` \(k\ge4\) records remain.  All `210`
quadratic-factor ansatz trials are also excluded.  This advances the finite
search boundary but leaves sharpness **OPEN**; see the
[degree-29 reduction](gram_belyi_m5_degree29_reduction.md) and the
[three-black theorem](gram_belyi_m5_degree29_three_black_ww.md).

The four-black pass is likewise quotient-first and exact.  Its `42,994`
marked passports reduce to `12,752` necessary master systems.  The first
univariate pass removes all two-black-endpoint and all
one-black/critical-white cases.  A sharp degree-`756` critical resultant and
`10,696` exact tower-field endpoint-gcd checks then eliminate every `WW`
status, another `2,560` systems / `8,960` records.  The final `BW+WB`
calculation tests both fixed white-endpoint values on all `28,159` critical
tower orbits and eliminates the last `3,458` systems / `10,374` records.
Thus four black vertices are impossible: a witness needs at least five, and
the global degree-`29` remainder is `42,562,335`.  Sharpness remains **OPEN**;
see the [four-black reduction](gram_belyi_m5_degree29_four_black_reduction.md),
the [four-black `WW` theorem](gram_belyi_m5_degree29_four_black_ww_simple.md),
the [four-black completion theorem](gram_belyi_m5_degree29_four_black_one_black_simple.md),
and its [hostile audit](gram_belyi_m5_degree29_four_black_one_black_simple_hostile_audit.md).

At five black vertices, the exact marked-passport quotient has `190,047`
raw records and `33,079` binary-status master systems.  The complete `BB`
block (`11,011` systems / `55,055` records) is now impossible by exact
projective finite-field certificates: the residual critical ideal is unit
on the affine chart and has no common leading-form point at infinity.
Properness supplies the characteristic-zero implication.  The exact
white-endpoint equations are status-free: `G=DQ+R` gives `Q=0` at every
white endpoint, simple or critical, collapsing `BW/WB/WW` to
`4,225/4,225/1,292` distinct necessary systems.  The first exact partition pass eliminates
`(5,1^24)` completely by the factorization `D=29(x-t)^4`, removing `3,876`
marked records.  The five-black residual is `9,742` status-free necessary
systems (`22,068` in the binary-status ledger) / `131,116` records, and the
global degree-`29` remainder is `42,503,404`.  For `(3^2,1^23)`, a uniform
cluster-count argument excludes projective infinity, and the exact `a>=21`
sublocus eliminates another `12` systems / `48` records.  The resulting raw
remainders are `5,120` in that partition, `131,068` in five-black, and
`42,503,356` globally.  The fixed `CC` factor `K=x(x-1)` gives exact
dimension-`24` base towers for all `1,292` black records; adjoining the value
equation always gives the unit ideal.  All `CC` passports are therefore
impossible, leaving `3,840` `SS/SC/CS` records in the partition, `129,788`
five-black records, and `42,502,076` globally.  The separate one-critical
towers then eliminate all `1,280` unresolved `CS` rows plus the exact
`41`-row `SC` tranche, leaving `2,519` partition records (`1,280 SS` plus
`1,239 SC`), `128,467` five-black records, and `42,500,755` globally; sharpness stays
**OPEN** and the broader census is stopped.  On the theory track,
`e(S_P)=f(f-1)+4` gives a genuine surface-native dimension bridge, while an
explicit fixed-volume stable-pair family with unbounded fields shows that
KSBA boundedness alone cannot give the field list.  See the
[five-black quotient](gram_belyi_m5_degree29_five_black_quotient.md), the
[five-black `BB` theorem](gram_belyi_m5_degree29_five_black_bb.md), and the
[five-black `WW` partition theorem](gram_belyi_m5_degree29_five_black_ww_partition5.md),
with its [hostile audit](gram_belyi_m5_degree29_five_black_ww_partition5_hostile_audit.md),
the [`(3^2,1^23)` high-degree theorem](gram_belyi_m5_degree29_five_black_ww_partition33_high_a.md),
its [hostile audit](gram_belyi_m5_degree29_five_black_ww_partition33_high_a_hostile_audit.md),
the [complete `CC` tower theorem](gram_belyi_m5_degree29_five_black_ww_partition33_cc_tower.md),
its [hostile audit](gram_belyi_m5_degree29_five_black_ww_partition33_cc_tower_hostile_audit.md),
the [`CS` theorem / `SC` tranche](gram_belyi_m5_degree29_five_black_ww_partition33_cs_full_sc_tranche.md),
its [hostile audit](gram_belyi_m5_degree29_five_black_ww_partition33_cs_full_sc_tranche_hostile_audit.md),
and the [surface-bound obstruction theorem](gram_belyi_carrier_dimension_field_bound_obstruction.md).

There is a complementary field-effectivity theorem which does not use
Prasad's arithmetic-volume formula or ABSW. The separate Mostow--Prasad
strong-rigidity theorem enters only as a uniqueness result after fixing a
Coxeter matrix. For a bounded universal coefficient-order and
Hilbert/Borcherds grammar \(\mathcal D\), form the exact-return stack
\(\mathscr H^{\rm rig}_{n,\mathcal D}/\mathbb Q\) on a scalar
splitting/eigencharacter cover.  The fields are outputs, not a finite list
hidden in \(\mathcal D\), and the Lorentzian embedding is a real enhancement.
For polytope fields it suffices that the scheme-theoretic image of the
outward unit-normal geometric Gram map (equivalently Cartan data modulo
diagonal root rescaling) be zero-dimensional: its finite Galois set
enumerates unit-normal geometric Gram fields and cyclic-product Vinberg
subfields even when surface fibres
move.  Zero-dimensionality and trivial stabilizers of the whole carrier
stack are needed only to enumerate carrier fields and descend the complete
surface object.  Bounded Hilbert polynomials provide only finite type.

The [field-effectivity report](dessin_native_field_effectivity_theorem.md)
freezes the exact hypotheses and an executable H11/26 calibration.  Its
rank locus \(a^2-b^2=0\) leaves, after rejecting the nonreduced branch, one
reduced projective point \([1:1]\) with coordinate and Vinberg field
\(\mathbb Q\).  This is one chart point, not a proof that the variable-class
Gram image is zero-dimensional.  Bugaenko's natural K3 carrier has a
two-dimensional period family, so it fails the stronger full-carrier gate;
it may nevertheless have constant Gram image and is not a counterexample to
polytope-field effectivity.

This is not advertised as the first abstract finiteness mechanism.  Agol,
Belolipetsky, Storm, and Whyte proved finiteness of maximal arithmetic
hyperbolic reflection groups, and Nikulin produced explicit finite ground-field
sets and an effective-classification route in principle.  The purpose of the
Coxeter dessin is sharper: turn those finite arithmetic possibilities into a
computable algebraic/automorphic invariant, a finite marked monodromy object,
and necessary-and-sufficient carrier equations which can produce explicit
survivors and exclusions.

## 6. Completion path

The Borcherds-dessin expedition has five exact terminal steps.

1. Compile the obstruction matrix and the least positive safe correction for
   each calibration row: compact Bugaenko \(H^8\), toric \(H^6\), raw Coble
   \(H^{10}\) and \(H^{11}\), and characteristic-three \(H^{14}/H^{17}\).
2. Express separation level, safe detection height, obstruction-matrix
   denominators, weight, and arithmetic volume in terms of
   \(n,[k:\mathbb Q],\operatorname{disc}(k),|A_L|\), level, and wall data.
3. Define a non-singleton P-only class on which the normalized slope is
   uniformly bounded, or prove that every member is pure reflective, and
   prove that this bounds \(\deg_{\rm Cox}\).
   Independently test whether the same class bounds the absolute
   Gram--Belyi polarization degree; if it does, enumerate the finite passport
   list directly and skip the maximal-base/Mostow route.
4. Enumerate the resulting bounded-degree coset dessins and apply an
   independent reflective/Hilbert-orthogonal classification to their
   algebraic carriers.
5. Replay every survivor through the complete Vinberg finite-volume test and
   return either a new chamber, a sharper dimension/facet bound, or a genuine
   nonexistence theorem.

Until Steps 2--5 close, the result is a universal algebraic carrier and a
sharp carrier criterion, not the requested final classification theorem.

## 7. Dependencies

- [Arithmetic Coxeter-dessin programme](arithmetic_coxeter_dessin_program.md)
- [Polytope-to-automorphic-divisor transfer](polytope_to_automorphic_divisor_transfer_theorem.md)
- [Compact inverse-leverage audit](borcherds_compact_inverse_leverage_audit.md)
- [Finite-volume cusp transfer](finite_volume_coxeter_dessin_cusp_transfer.md)
- [Independent cusp hostile audit](finite_volume_coxeter_dessin_cusp_transfer_hostile_audit.md)
- [Finite congruence-monodromy theorem](coxeter_congruence_dessin_theorem.md)
- [Bounded-degree Schreier--Coxeter dessin theorem](bounded_degree_schreier_coxeter_dessin_theorem.md)
- [Exact six-row calibration ledger](arithmetic_coxeter_dessin_calibration_ledger.md)
- [Dessin-native field-effectivity theorem](dessin_native_field_effectivity_theorem.md)
- [Universal Gram--Belyi surface theorem](gram_belyi_coxeter_dessin_theorem.md)
- [Gram--Belyi hostile audit](gram_belyi_coxeter_dessin_theorem_hostile_audit.md)
- [Gram--Belyi standard-packet calibration](gram_belyi_standard_packet_calibration.md)

Primary classification/obstruction sources are Vinberg's *Hyperbolic
reflection groups*, Bruinier's
[*Regularized theta lifts for orthogonal groups over totally real fields*](https://arxiv.org/abs/0908.3076),
Wang's
[*2-reflective lattices of signature `(n,2)` with `n>=8`*](https://arxiv.org/abs/2301.11536),
and Metzler's
[*A new converse theorem for Borcherds products*](https://arxiv.org/abs/2512.20701).
Metzler assumes a hyperbolic split and therefore does not close the compact
anisotropic Hilbert converse.
For the pre-existing arithmetic finiteness boundary see
[Agol--Belolipetsky--Storm--Whyte](https://arxiv.org/abs/math/0612132) and
[Nikulin's ground-field classification programme](https://arxiv.org/abs/0710.2340).

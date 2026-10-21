# Nikulin as a hostile benchmark: an AG translation ledger

**Status (2026-08-11): project-native research specification with one
explicit downstream classification import.**  The AG carrier/rank/field
theorems below do not use a Nikulin--Vinberg reflective-lattice
classification or dimension theorem.  The projected K3 construction uses the
standard primitive-embedding criterion.  Separately, the new passage from
finitely many native pairs to finitely many abstract maximal chambers
explicitly imports Nikulin's 1982 fixed-\((k,m)\) finiteness theorem for
reflective lattices.  That last step is not AG-only.  Published arguments are
split into exact finite statements that we can replay, scope failures that we
must not inherit, and analogies that become new algebraic-geometric proof
obligations.

The central separation is now

```text
raw Coxeter packet
        -/-> exact AG trace carrier;

fixed rational arithmetic P
        -> finite RLogComp(P)
        -> exact weak-log carrier stack or branch-specific EMPTY;

pinned exact compact nonrational realization P
        -> finite K3Comp(P);

K3Comp(P) != empty
        <-> primitive full-trace, full-Picard K3 eigen-nef return
            with an ample interior point;

K3Comp(P) = empty
        -> exact PHYS, EMB, or ROOT failure in this K3 branch,
           not nonexistence of P or of every alternative AG carrier;

maximal compact (2,5,6)
        -> K3Comp(P) = empty already at same-vector PHYS;

absolutely maximal compact Delta(2,4,6)
        -> direct rational = weak-log rational
                         = exact-return minimalizable trace = empty
        -> the current exact three-branch Property-4 disjunction is false;

the same Delta(2,4,6)
        -> rank-six K3 Picard lattice + rational self-adjoint projector
        -> exact projected eigen-nef return with three physical root facets
           and no hidden restricted cuts
        -> positive one-input projected calibration with n<=19;

rank-20 X3, as one input to the all-finite-volume exact-module proposal
        -> projector complement vanishes
        -> norm-(-2) and divisibility-three norm-(-6) lines force
           a rational square root of 3
        -> the universal all-finite-volume complex projected-root K3
           assertion is false, while the maximal-only assertion is unsettled
           because X3 is not proved absolutely maximal;
        -> mixed-divisor K3 return remains positive;

absolutely maximal Borcherds H21/210
        -> direct rational = RLogComp = complex full-trace/projected K3
                           = empty
        -> characteristic-2 mixed curve/Cremona K3 carrier is positive
        -> the next category must include mixed divisor/action data or
           another carrier species;

pinned declared-host P17, P18, P19
        -> the same three branches are empty in that relative-host scope,
           without a claim of absolute arithmetic maximality;

exact admissible finite-volume full restriction-of-scalars
trace/order-module complex-K3 carrier in n>=2
        -> rank/dimension bound and finite field set.

finite-sided finite-volume literal full-Picard/full-nef K3 return
        -> rho physical (-2) facets spanning NS_Q with intersections <=13
        -> finite Picard hosts
        -> one complete-root Weyl orbit per oriented host, without maximality;

compact primitive-trace eigenslice return in a larger Picard lattice
        + physical in-trace (-2) star
        + exact native return
        -> finite native hosts;

finite native hosts + full same-host reflection maximality
        -> finite unmarked normalized chamber classes.

finite native pair set (k,m) + reflection maximality
        -> Nikulin 1982 fixed-data reflective-lattice finiteness
        -> finite unmarked normalized chamber classes,
           without finite projected ambient hosts but with a
           classification import.
```

The old proposal for the missing first arrow—essential surjectivity and exact
return into exactly the direct rational, weak-log rational, and minimalizable
full-trace branches—is now disproved by \(\Delta(2,4,6)\).  Its projected
rational K3 repair is one nonempty calibration of a broader projected K3
eigenslice category: a rank-six Picard lattice and rational self-adjoint
projector return the exact triangle with three physical projected roots and
no hidden cuts.  The broader K3 category inherits \(n\le19\); nonrational
members must retain the full restriction-of-scalars trace/order-module and
admissibility data needed by the field theorem.  The rational calibration
lies outside the old literal rational
full-nef Property-2 contract, and it cannot be universal.  Rank-$20$ **X3**
is one exact-module square-class counterexample to the universal assertion
over all finite-volume inputs while retaining a mixed-divisor K3 realization.
It is not proved absolutely maximal, so maximal-only universality is not
settled by this example.
The absolutely maximal rank-$22$ Borcherds chamber then empties every current
complex rational/weak-log/projected-root K3 branch, while a
characteristic-$2$ mixed curve/Cremona carrier remains positive.  The live
category must therefore contain mixed divisor/action data or another carrier
species.
Once an exact full restriction-of-scalars trace carrier is present, a
collection of unrelated local diagrams is unnecessary for the K3
dimension/degree bound.  If the carrier is additionally admissible,
finite-volume, equipped with the finite projective trace/order module, and in
dimension at least two, compact narrow-facet geometry or the noncompact cusp
fixed line gives the field set.  In the literal finite-sided finite-volume full-Picard/full-nef
K3 lane, the same narrow-facet geometry now supplies a bounded spanning
physical-root basis even at cusps, hence finite Picard hosts and one
complete-root Weyl orbit per oriented host without maximality.  A single
cusp link can leave isotropic divisibility unbounded, but the two
extreme-point packets in one global narrow star close that escape.
Primitive trace eigenslices in a larger Picard lattice remain separate:
projected physical walls do not suffice, and their compact host theorem
still requires that the bounded star consist of the same trace-module
vectors as physical square-\((-2)\) curves together with exact native
return.

For abstract maximal-chamber counting, finite projected ambient hosts are
unnecessary.  A finite set of native pairs \((k,n+1)\), together with
Nikulin's 1982 fixed-data reflective-lattice theorem and reflection
maximality, gives finitely many normalized chambers.  This conditional
downstream theorem imports classification and neither constructs the carrier
nor proves global ABSW from AG.

Within the compact nonrational full-Picard K3 subbranch, a pinned exact
realization has become a finite nonemptiness question.  It is not
automatically nonempty: the absolutely maximal \((2,5,6)\) triangle fails
already at `PHYS`.  The finite PHYS/EMB/ROOT compiler accepts exactly when
\[
 \operatorname {NS}(S)=M_{\mathbf Z},\qquad
 p\in\operatorname {Amp}(S)\cap E,\qquad
 P=\mathbf P(\operatorname {Nef}(S)\cap E),
\]
and otherwise returns a branch-specific `EMPTY` certificate.  This is full
Picard but generally a proper eigen-nef section, not the literal full-nef
object governed by the FV13 spanning-root theorem.

This logical bridge is now a theorem.  A faithful injection of the full
restriction-of-scalars root module into a rank-`B` integral packet gives
`{[k:Q]}(n+1)<=B`; a total ample facet degree `D` independently gives
`n<=D-1`.  Fixed `B,D,C,t` plus one finite return grammar gives a terminating
exact enumeration.  What remains open is deriving the carrier from every raw
maximal chamber after enlarging the now-disproved three-way carrier
disjunction, not the linear-algebraic implication itself.  The
[integral-packet theorem](integral_packet_global_dimension_bridge.md) and its
[hostile audit](integral_packet_global_dimension_bridge_hostile_audit.md)
freeze the exact boundary.

**Current AG replacement status (2026-08-11).** The global faithful packet
is no longer merely formal on the verified carrier branches. Direct rational
full-nef geometry gives `n<=11`; the certified weak-log/orbifold branch gives
`n<=25`; and a primitive full-trace complex-K3 carrier gives
`{[k:Q]}(n+1)<=20`. These are independent AG rank theorems, not translations
of the Nikulin--Vinberg reflective-lattice classification. Essential
surjectivity into the present three exact-return branches is false; essential
surjectivity into a revised broader AG carrier stack remains open.  A
rational member positively calibrates the broader projected K3 candidate on
the decisive counterexample and has \(n\le19\), but **X3** disproves
exact-module projected-root universality and the absolutely maximal rank-$22$
Borcherds chamber forces a mixed divisor/action or other branch.
Nonrational K3 members retain the full restriction-of-scalars trace/order-module
and admissibility hypotheses whenever the field theorem is invoked.

The field-height problem is now closed after exact admissible finite-volume
full restriction-of-scalars trace/order-module complex-K3 return in every
dimension `n>=2`. Noncompactness gives a rational null line fixed by
the full cusp packet and forces `k=Q`. Compactness gives a minimal-aperture
star with pairings `<14` spanning the Lorentz space; its connected Lorentz
basis block generates the full coefficient field, so `dm<=20` yields
`|Disc(k)|<=(32*14^10)^30`. The older nondotted-core, simplex and polygon
theorems remain sharper subcases. The unresolved packet is geometric:
simultaneous reduced irreducible effectivity, the actual low-degree section
algebra/divisor primality, contraction/Picard compatibility, and exact
chamber return. An NS embedding and reflection maximality do not supply it.

The literal full-Picard host problem now has a uniform compact/cusped
solution after return.  For \(3\le\rho\le20\), one minimal-aperture closed
star contains \(\rho\) actual square-\((-2)\) facet curves spanning
\(\operatorname {NS}_{\mathbf Q}\) and satisfying
\(0\le E_i.E_j\le13\).  Hence
\(R\subseteq\operatorname {NS}\subseteq R^*\) leaves finitely many Picard
hosts, and completeness of the full root arrangement gives one Weyl orbit
per oriented host without maximality.

The rational weak-log search is now uniform as well, but only internally to
that species.  `RLogComp` finitely enumerates physical squares, rational
Picard hosts, all facet-supported MMP programs, simultaneous prime-divisor
loci and three noncircular nef gates.  In the big case Fujino gives the fixed
necessary system `|192D|`; in the null case a positive effective affine-cusp
divisor proves nefness without semiampleness.  The result is a terminating
constructible carrier stack or `EMPTY`, with exact full-nef return after
adapted-host maximality.  Its generic elimination is not implemented, and an
`EMPTY` answer has no force on other AG species.

The compact **finite-native-host** problem for a primitive-trace eigenslice
has a different, precise solution.  To invoke that host-finiteness theorem,
every wall normal in a spanning minimal-aperture star must be represented by
that same in-trace physical square-\((-2)\) vector; the resulting placewise
trace-energy identity bounds the field Gram, coefficient order and all
integral same-span modules.  Exact native return, with
same-host maximality when only a selected subarrangement was initially
given, supplies the complete native mirror statement.  A fixed golden
\((2,5,5)\) scaling/glue family proves that roots merely projecting to the
same walls leave infinitely many hosts.  Same-vector normalization is not
necessary for exact eigenslice **existence**: the projected
\(\Delta(2,4,6)\) K3 package is the explicit counterexample, but it does not
inherit finite-host control.  The affine
\(U(d)\oplus\langle-2\rangle\) family remains a local warning, not a
finite-volume counterexample to the full-Picard theorem.  More strongly, one
fixed Jacobian K3 has projective torsors with
\(\operatorname {NS}=U(q)\oplus E_8(-1)^2\) for every \(q\ge2\), while
the two \(II^*\) packets, Euler number, Picard rank, Mordell--Weil rank and
discriminant length stay fixed.  This selected-fibration family is still
only cusp-local.  If finite volume is added, the global basis theorem gives
\(q^2\le2877^9\); the sharper \(q\le390\) needs a cusp-adapted
transverse-wall hypothesis.

## 1. Hostile classification of the inherited steps

| Published-style step | Independent status here | What may be retained |
|---|---|---|
| A facet of minimal visual aperture gives pairings `<14` among walls meeting it | Elementary cone calculation survives; normalization and closure-at-infinity must remain literal | A canonical bounded local Gram window |
| The complete closed star of that facet spans the Lorentz space | Valid for every finite-sided finite-volume chamber in dimension at least two: two distinct extreme points of the facet have complete incident packets spanning their orthogonal hyperplanes; bounded horospherical links handle nonsimple ideal vertices | A bounded wall-normal spanning packet, including all cross-pairings between the two endpoint packets |
| The Lorentz basis block inherits the whole coefficient field | Proved without a minimal-field premise: a proper cycle subfield gives a second real embedding and Cartan diagonal similarity, contradicting admissible definiteness | Integral bounded generators on at most `n+1` nodes |
| A finite-volume literal full-Picard/full-nef K3 return gives finitely many Picard hosts | Proved for `3<=rho<=20`: extract `rho` physical `(-2)` facets spanning `NS_Q` with intersections in `{0,...,13}`, then use `R subset NS subset R*` | Finite Picard hosts and one complete-root Weyl orbit per oriented host, without maximality |
| A spanning star of physical K3 `(-2)` roots gives finitely many compact primitive-trace hosts | Proved only when those are the same vectors in the trace module: the fixed square gives a placewise energy budget, the codifferent bounds the field Gram, and `F subset M subset F#` bounds saturations | A finite complete native-host list for the eigenslice scope; exact native return remains literal |
| Existence of a compact full-Picard K3 eigen-nef return from a fixed exact raw realization requires an unbounded search | False in the pinned nonrational scope: `PHYS`, exact discriminant-complement `EMB`, and compact hidden-root `ROOT` are finite and exhaustive; acceptance includes `[p] in P°` with `p in Amp(S) cap E` | A finite iff/`EMPTY` decision for that K3 branch; `EMPTY` remains compatible with the polytope or another AG branch, and a proper eigenslice is not FV13 literal full-nef scope |
| Raw reflection maximality forces the compact K3 compiler to be nonempty | False: the absolutely maximal \((2,5,6)\) triangle has a complete five-element codifferent trace-\((-2)\) shell, but inert-prime valuation parity makes its same-vector `PHYS` square-class intersection empty | Maximality can remove selected hidden native walls only after physical normalization, embedding, zero-projection exclusion and same-host preservation are supplied; projected larger-Picard returns remain open |
| The weak-log rational species requires an unbounded surface search | False: `RLogComp` gives finite `HOST/MMP/PRIME/NEF` constructible loci and exact maximal return; the \(D^2=0\) affine relation and the \(D^2>0\) fixed \(|192D|\) gate close nefness | A finite iff/`EMPTY` compiler for this species, not a proof that its output or any broader Property-4 stack is nonempty |
| The current direct/weak-log/minimal-trace exact-return disjunction might be universally nonempty after more implementation | False: the absolutely maximal compact \(\Delta(2,4,6)\) chamber empties all three branches by rational anisotropy/compactness and the \(K=0\) \((-2)\)-angle obstruction | Replace the disjunction with mixed divisor/action data or another AG category; the projected \(\Delta(2,4,6)\) packet is only a scoped calibration, **X3** kills universal projected-root return, and the bare carrier-relative third-row implication is unchanged |
| The natural projected rational K3 repair might also be empty on \(\Delta(2,4,6)\) | False: the even rank-six lattice with Gram \(\left(\begin{smallmatrix}15G&G\\G&-2I_3\end{smallmatrix}\right)\), a rational self-adjoint projector, and an exact finite root-window census returns precisely that triangle | A positive one-input calibration with \(n\le19\), not literal rational full-nef Property 2, an integral raw-group action on all of the Picard lattice, a class-level P4 result, or a universal projected-root branch |
| That projected-root repair can be promoted to every finite-volume chamber in the complex K3 rank range | False: at rank \(20\), **X3** saturates Picard rank, so the projector is the identity; its primitive norm-\((-2)\) and divisibility-three norm-\((-6)\) lines would force \(3\) to be a rational square under one common exact-module similarity | One counterexample refuting the universal assertion over the all-finite-volume class, not a theorem that every input fails and not a maximal-only no-go; **X3** has a weaker mixed-divisor K3 realization |
| Every absolutely maximal chamber is routed by the current complex rational/weak-log/projected-root K3 stack | False: the rank-$22$ Borcherds $210$-wall chamber is absolutely reflection-maximal, has no complex K3 injection by rank, and exhausts the direct rational and weak-log physical-square branches to emptiness | A characteristic-$2$ supersingular K3 mixed curve/Cremona carrier remains positive; the next universal species must include mixed divisor/action data or another carrier |
| Finite projected ambient hosts are needed before abstract maximal chambers can be counted | False after finite native pairs are known: Nikulin 1982 gives finitely many reflective \(\mathcal O_k\)-lattice similarity classes for fixed \((k,m)\), maximality gives the full native reflection group, and the complete arrangement has one chamber orbit | A classification-importing conditional finite-chamber theorem; not AG-only, not a carrier compiler, and not global ABSW |
| The standard `P17`--`P19` chambers show the same failure absolutely | Too strong: they empty all three branches for the pinned `I_(1,n)` hosts, but only declared-host reflection completeness is certified | A broad relative-host obstruction and a regression corpus; no claim of absolute/maximal-arithmetic reflection maximality, all normalizations, or all carrier species |
| Physical roots merely projecting to the star give the same finite-native-host conclusion | False: the golden `(2,5,5)` family keeps the chamber and native cycles fixed while trace scale, glue index and ambient discriminant diverge | Require literal in-trace physical normalization only for `K3-HOST-PHYS/COMP` finiteness, not for exact projected eigenslice existence; projected ambient-host finiteness is unnecessary for abstract maximal-chamber counting once finite native pairs and reflection maximality are known |
| A full finite-volume cusp packet forces the field to be rational | Proved: bounded horosphere section gives a one-dimensional common fixed null line, impossible after definite conjugate base change | Every exact noncompact admissible carrier has `k=Q` |
| The affine physical cusp packet therefore gives finitely many K3 hosts | False locally: `U(d)+<-2>` has a fixed affine packet and unbounded primitive isotropic divisibility `d` | Global finite volume plus literal full-Picard/full-nef return selects a bounded spanning physical-root basis; the local family is not a counterexample |
| Fixing the Jacobian, all fibre schemes, full local parabolic rank, Euler `24`, Picard rank `18`, MW rank zero and discriminant length two bounds the cusp index | False: projective K3 torsors of one `2 II*+4 I1` Jacobian have `NS=U(q)+E8(-1)^2` and fibre divisibility `q` for every `q>=2` | These are not finite-volume chambers; with that extra global hypothesis one gets `q^2<=2877^9`, while `q<=390` additionally requires a cusp-adapted transverse wall |
| The narrow facet is simple and its two-face average is bounded by its `h`-vector | Valid for compact Coxeter chambers | A compact incidence inequality |
| The same simplicity argument works for every finite-volume chamber | False: an ideal vertex need not be simple; an ideal right-angled octahedron already has four walls at a cusp in `H^3` | Replace by a separate affine-link/cusp theorem |
| A minimal hyperbolic subdiagram in the `u<2` branch is Lanner | **False as printed:** the marked arithmetic orthoscheme `[3,3,6]` has `u=1`, is minimal indefinite, and has an affine deletion | Split compact Lanner and cusped quasi-Lanner branches; the latter has field `Q` by the isotropic-vector argument |
| A non-negative-definite principal packet has the full ground field | Independently proved by cyclic-product sign switching; includes hyperbolic, parabolic, and quasi-Lanner packets | A primitive coefficient subpacket, provided its full cyclic-product algebra is retained |
| Exactly five marked four-node families remain | Correct schematically: exact generation gives seven literal support strata, with two label-`2` edge-deletion strata inside the five drawings | A complete local finite grammar, not a global rank bound |
| The auxiliary polynomial printed in arXiv `0910.5217v4`, equation (24), has norm `0.0905969664` | False as printed: its true norm is about `0.1598126957` | Nothing from that printed row |
| Replacing its last exponent `2` by `3` restores the degree-23 branch | Exact; the final journal equation makes this repair silently | Our independent degree-20 interval certificate, not the citation |
| A finite high-dimensional field list gives a global dimension bound | Invalid implication | No translation; a new rank theorem is mandatory |
| Reflection maximality forces a Borcherds correction or rigid carrier | False: the maximal reflection chamber `Gamma16` has neither implication | Maximality must be translated into a different bounded AG invariant |

The compact/cusped rows above already force a split.  The last three
prevent us from mistaking a field funnel or a modular form for the desired
ABSW-type theorem.

## 2. The useful AG dictionary

The surviving ideas have the following nonformal translations.

| Hyperbolic datum | AG/cohomological datum | Exact gate |
|---|---|---|
| Integral root and its Gram products | Integral divisor class or integral correspondence action | Effectivity, integral action, and bounded positive/negative cycle degrees |
| Distinguished Lorentzian embedding | Distinguished real eigenspace of a trace-form action | Signature `(1,n)` and orientation are finite selectors |
| Negative-definite conjugate Grams | Hodge-index-negative conjugate eigenspaces | Every nonidentity embedding must pass an exact inertia test |
| Small hyperbolic subgraph generating the ground field | Primitive bounded coefficient subalgebra of the cohomological packet | Prove generation; do not identify `Q(alpha)` with a larger ambient field |
| Auxiliary integer polynomial on `[-2,1]` | Characteristic-polynomial/norm obstruction for an integral cohomological eigenvalue | Exact resultants or integer norm, never decimal extrema |
| Narrow face | Canonically selected bounded local intersection packet | Selection must be functorial from raw `P` |
| Ideal vertex and parabolic link | Labeled exceptional boundary cluster plus affine monodromy | Full parabolic rank and compact horosphere-cell test |
| Face-incidence average | Boundary/cluster incidence inequality | A surface SNC support has dimension at most one; higher faces require bounded extra decoration |

This dictionary is deliberately asymmetric.  Hodge theory can certify the
conjugate sign condition once the integral packet exists, but it does not
construct a bounded packet from an arbitrary chamber.

## 3. The proved local engine

The project-native polynomial

\[
 P_{20}(x)=[x(x+1)]^4[(x-1)(x+2)]^3[x^2+x-1]^3
\]

has exact norm

\[
              \lVert P_{20}\rVert_{[-2,1]}
                   ={884736\over9765625}.
\]

If an integral spectral coordinate `alpha` is primitive for its field, its
distinguished embedding lies in `[2,14]`, and all other embeddings lie in
`[-2,1]`, the integer-norm inequality proves

\[
                       [\mathbb Q(\alpha):\mathbb Q]\le23.
\]

For `t` such coordinates which jointly generate `k`, the current crude bound
is `[k:Q]<=23^t`.  This is executable and finite.  It does **not** bound `n`
unless the packet containing those coordinates is also globally faithful and
has bounded cohomological rank.

## 4. The exact global theorem target

The next publishable theorem is not another diagram census or another field
height lemma.  The dimension implication is proved once the exact full-trace
carrier exists; the field implication additionally requires the exact
admissible finite-volume full restriction-of-scalars trace/order module in
dimension at least two.  The target is the missing first arrow:

> **Exact carrier compiler target.**  From every reflection-maximal
> arithmetic finite-volume Coxeter chamber `P`, construct functorially a
> finite exhaustive list of integral AG/cohomological packets `H(P)` such
> that:
>
> 1. each nonrational packet contains the entire restriction-of-scalars trace
>    module, so `rank H(P) >= {[k:Q]}(n+1)`;
> 2. each old rational packet returns the literal full nef cone; a projected
>    K3 packet, over either a rational or nonrational input, instead records
>    its self-adjoint projector and returns the exact eigen-nef section;
>    every nonrational packet also retains the full trace and admissibility
>    data from item 1;
> 3. compact faces and cusped affine links, the canonical/log class, and the
>    complete hidden-wall arrangement are recovered exactly;
> 4. actual section kernels, multiplication maps and divisor primality certify
>    simultaneous reduced irreducible effectivity and every required
>    incidence;
> 5. each finite-sided finite-volume literal full-Picard/full-nef K3 packet
>    supplies its automatic bounded spanning physical square-`(-2)` star;
>    a nonrational primitive-trace eigenslice seeking the existing
>    finite-native-host theorem instead supplies a literal physical in-trace
>    star and exact native return; a projected K3 packet records
>    actual ambient root lifts, projector/complement/glue data, and exact
>    zero-projection/hidden-cut tests.  Exact return requires those data, not
>    class-wide finiteness of projected ambient hosts; and
> 6. literal complete full-Picard root return needs no maximality for the
>    one-Weyl-orbit conclusion; for a weaker selected native subarrangement,
>    maximality is measured in the complete marked native `k`-reflection
>    host, not only among ordinary trace-lattice roots or
>    complement/period-extendable isometries; and
> 7. for abstract maximal-chamber counting, the output native pairs
>    \((k,n+1)\) lie in a finite set.  Nikulin's fixed-data theorem then
>    replaces any need to finite-ize projected ambient Picard hosts, but is a
>    separate classification input.

The compact nonrational full-Picard K3 component of this target is now a
finite inner compiler for each pinned exact realization.  Its three terminal
gates are `PHYS`, `EMB`, and `ROOT`; a survivor gives the exact ample
eigen-nef return, while failure gives a rigorous branch-specific `EMPTY`
certificate.  The rational weak-log component likewise has the finite
`RLogComp` iff/`EMPTY` compiler.  The outer quantifier is no longer merely
unproved in its advertised three-way form: \(\Delta(2,4,6)\) disproves it.
A rational member of the broader projected K3 branch now returns that
counterexample on a rank-six Picard lattice, and K3 rank gives \(n\le19\) in
the whole category.  Nonrational members must carry the full
restriction-of-scalars trace/order-module and admissibility data needed by the
field theorem.  The revised target cannot be
universal projected-root K3: rank-$20$ **X3** is one exact-module
counterexample to the universal assertion over all finite-volume inputs, but
does not settle maximal-only universality because it is not proved absolutely
maximal.  The absolutely maximal rank-$22$ Borcherds
chamber lies outside the whole current complex stack.  It must instead add
mixed divisor/action data or another category and prove that revised
disjunction exhaustive.  Nothing here excludes all possible surfaces or
varieties.

If `H(P)` lies in the narrowed log-Fano carrier class, boundedness gives

```text
f <= C_epsilon/delta,
{[k:Q]}(n+1) <= rho(S) <= R_epsilon.
```

The weak-del-Pezzo special case gives the explicit
`{[k:Q]}(n+1)<=9`.  These are genuine dimension inequalities after the packet
exists.  Essential surjectivity and the scheme-theoretic item 4 remain open.
The bounded spanning-root step is no longer open after literal finite-volume
full-Picard/full-nef K3 return, but it must not be applied to a projected
eigenslice or renamed a consequence of maximality.

One nonempty alternative lane is now a theorem.  If `H(P)` is the numerical
lattice of a minimal complex surface and every simple reflection extends to
a pairing isometry fixing the canonical class, then canonical orthogonality,
Hodge index, and minimal-surface classification give
`{[k:Q]}(n+1)<=20`; the characteristic-`p` ceiling is `22`.  A Wehler K3
realizes the rational `H2` ideal triangle by actual automorphisms.  The open
step is essential surjectivity from raw maximal chambers, not the dimension
inequality.  Effective nef-wall reflections belong only to the broader
canonical-isometry lane, since an actual automorphism cannot send an
effective root to its negative.

That lane is now sharpened in two directions.  First, finite numerical orbit
of one geometric exceptional lattice is equivalent to orthogonal/isometric
descent of the trace block, so complement rank is not the missing invariant.
Second, operators are unnecessary when geometry places the full trace block
in `K^perp`: a minimal carrier gives the same `20/22` ceiling directly.  On
a rational carrier, the complete Mori alphabet can itself force this
canonical geometry.  A `(-1),(-2)` alphabet makes `-K` nef and gives
`rho<=10`; the disjoint `(-4)` Coble alphabet makes every Mori generator
orthogonal to `K+B/2`, hence forces `B=-2K` and a K3 double cover.  These are
independent AG rank theorems, not translations of a Nikulin dimension proof.
Their global raw-P existence gate must now be reformulated through that
broader outer disjunction; see the
[equivariant descent theorem](equivariant_minimalization_trace_descent_theorem.md),
[canonical-cover theorem](canonical_orthogonal_k3_cover_dimension_theorem.md),
and [rational Mori theorem](rational_mori_boundary_k3_dimension_dichotomy.md).

The rational reflection alphabet is now reconstructed without a diagram
classification.  Unimodularity, integral reflection, adjunction and
characteristic parity give only `(-1),(-2)` and the twice-primitive `(-4)`
case.  Coherent preservation of `K+B/2`, not reflection maximality, is the
exact Coble incidence gate; a maximal `(2,3,infinity)` full-nef surface is an
explicit counterexample to the latter implication.  A separate adapted-cover
calculation proves that every disjoint standard-coefficient `(-2m)` attempt
with `m>2` violates Noether.  These are project-native AG deductions, not a
translation of a Nikulin/Vinberg dimension theorem; see the
[reflection-integrality gate](rational_nef_reflection_integrality_gate.md)
and its
[hostile audit](rational_nef_reflection_integrality_gate_hostile_audit.md).

## 5. Finite experiments that decide the next lemmas

1. **No-import four-node generator.**  Generate marked symmetric Gram graphs
   from the endpoint-deletion and inertia predicates alone.  Compare the
   output with, but do not seed it from, the alleged five families.  Complete:
   `32` masks give `17` admissible masks and seven marked support orbits; the
   five schematic families cover all seven, with no omission.
2. **Exact field inheritance.**  Given a hyperbolic principal subpacket,
   enumerate embeddings over its coefficient subfield and verify the
   diagonal-similarity/signature contradiction.  Include the parabolic and
   quasi-Lanner branches.  Complete for every non-negative-definite principal
   packet by a general spanning-tree sign-switch proof.
3. **Compact face replay.**  Certify the exact `h`-vector average inequality
   and the component-size-to-dimension inequality.  Archive every scope
   hypothesis.  This replay is complete: the exact `<14` cone window and the
   strict `h`-vector average give `n<=9` under the explicit hypothesis that
   every full finite vertex graph along the narrow facet has components of
   size at most two.  Deriving that hypothesis remains open.
4. **Cusp replacement.**  Complete for exact finite-volume admissible
   carriers: boundedness of the horosphere cell kills its recession cone, so
   the incident fixed hyperplanes meet in exactly one null line; definite
   conjugate base change forces `k=Q`.  This is independent of cusp
   simplicity and does not reuse the compact proof.
5. **Six-carrier rank table.**  On the `H8`, toric `H6`, Coble `H10/H11`, and
   characteristic-three `H14/H17` controls, record the smallest integral
   packet which contains the full trace module, its cycle/monodromy degree,
   and whether the Gram is intrinsic or external.  This is complete: every
   row has an intrinsic trace injection, with ranks `18,7,11,12,15,18`; only
   toric `H6` and MII `H17` have certified total ample facet degrees, `12`
   and `168`.  No row supplies a maximality-derived uniform bound.
6. **Local-to-global counterexample search.**  Look for a family with the
   same bounded local interval packets and unbounded faithful trace rank.  A
   positive family disproves the proposed packet theorem and prevents wasted
   surface work.  This experiment is now complete at the matrix/AG-obstruction
   level: `-2I+Adj(K_{1,m})` has a fixed seven-type at-most-four-node packet
   but unbounded Lorentzian rank, and exceptional curves on `Bl_r(P^2)` have
   fixed local intersection packets and individual degree one but unbounded
   Picard rank.  The star tail fails the global finite-volume deletion gate,
   so this is a local-to-global obstruction, not a new polytope family.
7. **Maximality pushback.**  Enumerate the actual AG extension operations
   (finite marked quotient, log-MMP contraction, lattice/correspondence
   saturation) and test whether absence of all of them bounds any numerical
   packet invariant.  Primitivity alone is already excluded by prime-degree
   Belyi maps.

The stopping rule is strict: a local field list, a timeout, a positive
Borcherds form, or a bounded census is not a dimension theorem.  Promotion
requires either the bounded faithful packet implication or an exact
counterfamily identifying the next obstruction.

## 6. Linked exact artifacts

- [project-native interval theorem](fekete_interval_minus2_1_probe_theorem.md);
- [independent interval hostile audit](fekete_interval_minus2_1_probe_hostile_audit.md);
- [transition-polynomial hostile audit](nikulin_transition_polynomial_hostile_audit.md);
- [`[3,3,6]` conceptual counterexample and dependency audit](nikulin_conceptual_dependency_hostile_audit.md);
- [exact four-node family and field-inheritance audit](nikulin_four_node_transition_families_exact_audit.md);
- [global dependency hostile audit](ag_global_bound_and_transition_constant_hostile_audit.md);
- [narrowed log-Fano carrier theorem](log_fano_cohomological_coxeter_carrier_theorem.md);
- [log-Fano hostile audit](log_fano_cohomological_coxeter_carrier_hostile_audit.md);
- [global integral-packet dimension bridge](integral_packet_global_dimension_bridge.md);
- [integral-packet hostile audit](integral_packet_global_dimension_bridge_hostile_audit.md);
- [six-carrier packet calibration](six_carrier_integral_packet_rank_calibration.md);
- [six-carrier hostile audit](six_carrier_integral_packet_rank_calibration_hostile_audit.md);
- [compact narrow-face/cusped-scope audit](narrow_face_minimality14_combinatorial_reduction_hostile_audit.md);
- [minimal canonical-carrier dimension theorem](project_native_minimal_surface_automorphism_dimension_bound.md);
- [minimal-carrier hostile audit](project_native_minimal_surface_automorphism_dimension_bound_hostile_audit.md);
- [equivariant-minimalization trace-descent theorem and maximal-complement counterfamily](equivariant_minimalization_trace_descent_theorem.md);
- [equivariant-minimalization hostile audit](equivariant_minimalization_trace_descent_hostile_audit.md);
- [canonical-orthogonal/K3-cover dimension theorem](canonical_orthogonal_k3_cover_dimension_theorem.md);
- [rational Mori-boundary/K3 dimension dichotomy](rational_mori_boundary_k3_dimension_dichotomy.md);
- [rational Mori-boundary hostile audit](rational_mori_boundary_k3_dimension_dichotomy_hostile_audit.md);
- [rational-nef reflection-integrality and adapted-cover gate](rational_nef_reflection_integrality_gate.md);
- [reflection-integrality hostile audit](rational_nef_reflection_integrality_gate_hostile_audit.md);
- [direct full-NS AG dimension theorem](rational_direct_full_ns_n11_dimension_theorem.md);
- [weak-log/orbifold AG dimension theorem](rational_wrong_log_orbifold_rank_gate.md);
- [four-property completion audit](canonical_chamber_four_property_completion_audit.md);
- [`RLogComp` finite weak-log iff/`EMPTY` compiler](rational_weak_log_finite_iff_compiler.md),
  SHA-256 `48cf57f67dea15ed7a307047506e1d23ffd6841af6dcb29798d99e9eb5e34ffe`,
  and its
  [independent hostile audit](rational_weak_log_finite_iff_compiler_independent_hostile_audit.md),
  SHA-256 `648108238da96ea27bc3f16a3e4c69d1eb18ceef3539b84857ade85948b4c54e`;
- [maximal compact \((2,5,6)\) K3 `PHYS`-empty theorem](k3_compact_maximal_256_phys_empty.md),
  SHA-256 `9fa8afdb6f8e3f5225c494d48e3e9130a9b86cfe2aa8f58b38b177b9f972b83e`,
  and its
  [independent hostile audit](k3_compact_maximal_256_phys_empty_independent_hostile_audit.md),
  SHA-256 `6e249060fcef82242eeeb9bb7e6ae87328ad702f933ec28f184a275cfa84f29f`;
- [absolutely maximal \(\Delta(2,4,6)\) exact three-branch counterexample](compact_246_raw_maximal_property4_counterexample.md),
  SHA-256 `860d82c50ecfe46bd12c0f2fba851e85c6d08f6bf5e6b9b6ef575c75f31006b0`,
  and its
  [independent hostile audit](compact_246_raw_maximal_property4_counterexample_independent_hostile_audit.md),
  SHA-256 `688d5f81eb3c941228058f1e2a795d07e46581412584f5248e473e017e7a274d`;
- [projected rational K3 eigenslice repair for \(\Delta(2,4,6)\)](compact_246_projected_k3_eigenslice_carrier.md),
  SHA-256 `f4c4a038610e7c47600af971154ca056d7ae143546812582619e83ea7286cd1b`,
  and its
  [independent hostile audit](compact_246_projected_k3_eigenslice_carrier_independent_hostile_audit.md),
  SHA-256 `ca67ba5e6ce92e706506a3a8d870e03dfb21bd7c490ced65b1e809e8a4783d94`;
- [declared-host `P17`--`P19` three-branch exhaustion](p17_p19_three_branch_property4_escape.md),
  SHA-256 `6bc82285d0ca42c4e7eeaed40d2b6302182083557e88db9c57e20e842d67b4b0`,
  and its [independent hostile audit](p17_p19_three_branch_property4_escape_independent_hostile_audit.md),
  SHA-256 `587390d61f7a91de9579cc6bba19c6cd005da437bdf6830d4666f1aef82b3ebd`;
- [finite-native-pair maximal chamber theorem](projected_k3_maximal_abstract_chamber_finiteness.md),
  SHA-256 `59dd458acf14d18021f6c5583096da13f39c83115fc136d6ba91bcc850ce3603`,
  and its
  [independent hostile audit](projected_k3_maximal_abstract_chamber_finiteness_independent_hostile_audit.md),
  SHA-256 `edc663c0e7f35880123149d81d832812d5f28837c837df16eaa0f1daf725916a`;
- [**X3** rank-saturation exact-module no-go](projected_k3_universality_rank_saturation_no_go.md),
  SHA-256 `7ed5176c0c0a36e036c0e12dfa6acf1ee15d90d5ec7143ca5df59ebc3b7143d2`,
  and its
  [independent hostile audit](projected_k3_universality_rank_saturation_no_go_independent_hostile_audit.md),
  SHA-256 `89296f04bac6b0b4ed68d9b86936a08ede515ba5bf9c847abeaadbc628926071`;
- [Borcherds rank-$22$ outer-routing counterexample](borcherds_h21_outer_routing_counterexample.md),
  SHA-256 `fa899977f418d6fa93f7d5270e067c28a648c5e55f2636cb2b362e381534821d`,
  and its
  [independent hostile audit](borcherds_h21_outer_routing_counterexample_independent_hostile_audit.md),
  SHA-256 `8d705a619e9c1ffff3ad69765e90bb1a1a2a5d8713015ecca79382e8a1ccbc86`;
- [fixed-data closure and outer-routing synthesis](ag_to_absw_fixed_data_closure_and_outer_routing_obstruction.md);
- [exact noncompact ground-field collapse](exact_trace_noncompact_ground_field_collapse.md);
- [compact narrow-facet field-finiteness theorem](k3_compact_narrow_facet_field_finiteness.md);
- [compact theorem independent hostile audit](k3_compact_narrow_facet_field_finiteness_independent_hostile_audit.md);
- [compact physical-star trace-host theorem](k3_compact_physical_star_trace_host_finiteness.md);
- [physical-star independent hostile audit](k3_compact_physical_star_trace_host_finiteness_independent_hostile_audit.md);
- [same-native-host chamber composition](compact_trace_host_chamber_finiteness_composition_hostile_audit.md);
- [projected-root scaling/glue no-go](k3_compact_projected_root_scaling_no_go.md);
- [noncompact cusp-divisibility host gate](noncompact_k3_cusp_divisibility_host_gate.md);
- [noncompact cusp independent hostile audit](noncompact_k3_cusp_divisibility_host_gate_independent_hostile_audit.md);
- [fixed-Jacobian index-blindness theorem](noncompact_k3_fixed_jacobian_index_blindness.md);
- [fixed-Jacobian independent hostile audit](noncompact_k3_fixed_jacobian_index_blindness_independent_hostile_audit.md);
- [compact full-Picard K3 iff/EMPTY compiler](k3_compact_full_picard_raw_compiler_criterion.md);
- [compact compiler independent hostile audit](k3_compact_full_picard_raw_compiler_criterion_independent_hostile_audit.md);
- [finite-volume bounded spanning-root theorem](noncompact_k3_bounded_spanning_root_lemma.md);
- [bounded spanning-root independent hostile audit](noncompact_k3_bounded_spanning_root_lemma_independent_hostile_audit.md);
- [bi-anticanonical section-algebra property-4 gate](rational_pfirst_bianticanonical_section_algebra_gate.md);
- [nondotted field-generating-core theorem](k3_nondotted_field_generating_core_bound.md);
- [finite-volume simplex field theorem](k3_finite_volume_simplex_field_bound.md);
- [general compact-polygon field theorem](k3_n2_compact_polygon_field_bound.md);
- [prime-binomial raw-recognizer failure and `PB-Atlas` target](coble_prime_binomial_k3_rank_l4_pilot.md);
- [maximal reflection/no-correction counterexample](wang_h16_maximal_reflection_no_correction.md).

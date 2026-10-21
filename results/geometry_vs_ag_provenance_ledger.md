# Geometry versus AG: implication and provenance ledger

**Date:** 2026-08-12  
**Scope:** the composite theorem in
`agcat_global_maximal_arithmetic_coxeter_finiteness.md` and its frozen
dependencies.  This is a provenance audit, not a new proof and not an edit of
the owner theorem.

## Executive finding

The result is accurately described as an **AG-mediated arithmetic--geometric
hybrid**, but not as an AG-only proof and not as a proof independent of
Vinberg--Nikulin techniques.

The cleanest non-numerical accounting is by indispensable proof gates:

1. **Exact carrier:** algebraic-geometry packaging (projective quadric,
   polar hyperplanes, invariant quotient) plus finite real framing.  It is
   exact, but by itself proves no finiteness.
2. **Finite ground fields:** the only global finiteness gate substantially
   mediated by the AG carrier.  It is itself hybrid: hyperbolic geometry
   forces small wall packets; Cartan-cycle invariant functions make them
   integral points; arithmetic Fekete/norm arguments make those points
   finite.  The compact `H3` general-fork coordinate is the most clearly new
   part of this gate.
3. **Finite dimensions:** imported hyperbolic/polytope geometry
   (Vinberg and Khovanskii--Prokhorov), not a consequence of the AG carrier.
4. **Finite native hosts for fixed field and rank:** arithmetic lattice and
   order theory following the architecture of Nikulin's fixed-data method.
   The conclusion is reproved rather than cited; the conductor/projective-
   module implementation is a reconstruction and simplification.
5. **Maximal group to one chamber:** elementary reflection-group and
   hyperbolic gallery geometry.

Thus one of the four finiteness bottlenecks (fields, dimensions, hosts,
chamber orbit) is AG-arrangement-mediated; one is pure geometry, one is
arithmetic order theory, and one is elementary reflection geometry.  This is
more meaningful than assigning percentages, since deleting any one of the
four breaks the global theorem.

There is also an important methodological qualification.  The arrangement
coordinates are genuine regular functions on an algebraic quotient, but the
field proof can be rewritten entirely in the traditional language of Gram
matrices and Vinberg cyclic products.  No BAB, KSBA, K-moduli, Hilbert-scheme
boundedness, or higher-dimensional AG compactness theorem supplies the
finiteness.  In that sense the proof uses **algebraic invariant geometry and
arithmetic geometry**, not deep birational/moduli AG.  The dimension-two
orbifold-curve branch is the one place where an actual algebraic DM curve and
its genus are used; even there, terminal finiteness is the imported
Long--Maclachlan--Reid arithmetic/spectral theorem.

## Implication-by-implication ledger

The category labels below are:

- **AG/ArG:** algebraic geometry, invariant theory, or arithmetic geometry;
- **Hyp/Poly:** hyperbolic, Coxeter, or polytope geometry;
- **Lat/Ord:** arithmetic lattices, number fields, orders, and modules;
- **Elem:** elementary linear algebra or group theory.

| # | implication actually used | primary category | exact provenance | status in this project |
|---:|---|---|---|---|
| 1 | arithmetic reflection input `=>` admissible quadratic/lattice-form model | Lat/Ord | Vinberg Lemma 7, quoted as ABSW Lemma 5.2 | **imported structural lemma**; no finiteness content |
| 2 | chamber `=> (P(V),Q,H_i,sigma_0,C)` | AG/ArG + Hyp/Poly | projective quadric polarity and real half-space data | **new categorical packaging** |
| 3 | framed arrangement `=>` identical chamber | AG/ArG + Hyp/Poly | polarity recovers root lines; framing recovers inward rays/component | **new exact-return formulation**; elementary proof |
| 4 | diagonal root rescaling quotient `=>` closed Cartan-cycle coordinates | AG/ArG | classical torus invariant monomials; Vinberg cyclic-product language | **repackaged/reconstructed**, not a new invariant in general |
| 5 | stable projective lattice `=> b_C=tr(prod(1-s_i)) in O_k` | AG/ArG + Lat/Ord | Vinberg integrality plus projective-module trace | **modernized reconstruction** avoiding free-root assumptions |
| 6 | hyperbolic/parabolic principal packet `=>` packet sees the full field | Lat/Ord + linear algebra | Nikulin 1982, Theorem 2.1.3(2), with sign-cocycle proof replayed | **classical Nikulin principle, independently reproved and repaired in scope** |
| 7 | noncompact arithmetic chamber `=> k=Q` | Lat/Ord + Hyp/Poly | standard cusp-isotropy versus nonselected definiteness argument in the Vinberg model | **classical structural argument** |
| 8 | compact chamber `=>` simple; minimal-aperture facet `=>` pairings `<14` | Hyp/Poly | Nikulin narrow-face/minimality-14 architecture; elementary cone-plane inequality | **Nikulin technique, replayed** |
| 9 | finite edge packet `=>` small/Lanner/quasi-Lanner or one of five transition families | Hyp/Poly | Nikulin's finite-edge grammar; finite spherical/Lanner calculations | **independently enumerated**, including seven literal support strata and a repaired quasi-Lanner omission |
| 10 | transition tail/fork chart `=>` finite integral locus | AG/ArG + Lat/Ord | determinant regular functions, Vinberg integrality, classical relative Fekete and cyclotomic norms | **new exact norm implementation on a Nikulin-style chart**; does not import Nikulin's published `N(14)<=25` |
| 11 | compact `n>=6` `=>` a small/transition/polygon packet exists | Hyp/Poly, with optional Hodge interpretation | simple-polytope `h`-vector/Dehn--Sommerville face averages and Coxeter incidence; Khovanskii/Timorin-style Hodge package can explain the inequalities | **reconstructed global forcing**; not a general AG boundedness theorem |
| 12 | compact `n=4,5` `=>` exact lower-dimensional face or transition packet | Hyp/Poly + Lat/Ord | Nikulin-style narrow-face descent; literal orthogonal restriction over `k` | **independently reconstructed with rational-eigenspace repair** |
| 13 | compact `H3` arbitrary-label fork `=> k=Q(beta)` and finite labels/fields | AG/ArG + Lat/Ord + Hyp/Poly | new determinant coordinate `beta`, sign separation, norm/capacity estimates, Fekete | **strongest new step in the reconstruction**; uniform cutoff `<10^6` is certified but crude |
| 14 | Coxeter polygon `=>` exact genus-zero real DM orbifold curve and return | AG/ArG + Hyp/Poly | orbifold Riemann existence, GAGA/real descent, uniformization | **new exact-carrier formulation** |
| 15 | maximal arithmetic genus-zero Fuchsian groups `=>` finite polygon field/group list | AG/ArG + Lat/Ord/spectral | Long--Maclachlan--Reid Theorem 1.1; congruence spectral gap, Zograf inequality, bounded-area arithmetic finiteness | **imported terminal theorem**; the stack `M_0,r` alone would not suffice |
| 16 | all finite-volume Coxeter polytopes `=> n<=995` | Hyp/Poly | Vinberg compact nonexistence (`n>=30`) and Khovanskii--Prokhorov finite-volume noncompact theorem (`n>995`) | **imported wholesale**; independent of arithmetic and AG carrier |
| 17 | arithmetic chamber `=>` full `W`-stable projective lattice and full reflection host `R(L)` | Lat/Ord | normal core plus finite translate sum; arithmetic-lattice commensurability | **standard/reconstructed** |
| 18 | full host `=>` connected rank-`m` frame with two-cycles `<62^2` | Hyp/Poly | Nikulin 1981 Appendix Theorem 1 / 1982 Theorem 5.2.1 opening | **explicitly imported Nikulin geometric lemma** |
| 19 | bounded frame `=>` finitely many cycle packets | Lat/Ord | integral traces, selected/nonselected bounds, Minkowski discreteness | **Nikulin fixed-data architecture, replayed in invariant-cycle language** |
| 20 | cycle packet `=>` finitely many reflection representations | Elem + Lat/Ord | spanning-tree gauge, absolute irreducibility, Burnside | **reconstructed/simplified implementation** |
| 21 | fixed representation `=>` finitely many stable projective lattices | Lat/Ord | full order, conductor `cB subset Lambda`, ideal class group, finite `M/cM` submodules | **new conductor-sandwich implementation of a Nikulin-style fixed-data gate** |
| 22 | reflection maximality `=> W=R(L)` | Elem + Hyp/Poly | finite-index arithmetic overgroup argument | **elementary bridge** |
| 23 | full reflection host `=>` one unmarked chamber orbit | Hyp/Poly + Elem | local finiteness and finite-gallery reflection argument | **standard Coxeter geometry** |

The table deliberately permits mixed labels.  For example, row 13 would not
exist without hyperbolic incidence, but its finiteness is supplied by an
integral regular function and arithmetic capacity.  Conversely row 16 has a
Chow/Hodge interpretation in polytope algebra, but the composite theorem
imports the geometric dimension theorem rather than deriving it from the
exact arrangement carrier.

## Named-source firewall

### ABSW

The Agol--Belolipetsky--Storm--Whyte **finiteness conclusion is not used**.
Neither is the characteristic ABSW mechanism: conformal volume, eigenvalue
bounds, arithmetic volume/discriminant estimates, and bounded covolume
followed by arithmetic-lattice finiteness do not appear in the composite
proof.

One statement numbered in ABSW is cited: ABSW Lemma 5.2.  The primary text
itself labels it Vinberg's result and cites Vinberg Lemma 7.  Its role is only
to place arithmetic reflection groups in the quadratic-form model.  The
precise wording should therefore be:

> No ABSW-specific finiteness theorem or proof mechanism is used.  A
> structural Vinberg lemma is cited through ABSW Lemma 5.2.

Saying simply “ABSW is not used” is defensible at the level of conclusions,
but imprecise at the level of citations.

### Nikulin

Nikulin's **fixed-field/fixed-rank finiteness conclusion**, global ground-
field theorem, and published transition bound `N(14)<=25` are not invoked as
black boxes.  Nevertheless Nikulin's techniques are substantial and
indispensable to the present proof architecture:

- narrow facets and minimality `14`;
- full-field inheritance from a small hyperbolic packet;
- finite-edge/four-node transition grammar;
- codimension descent;
- arithmetic Fekete on conjugate Gram intervals;
- the connected frame with constant `62`; and
- the fixed-data strategy “bounded frame coordinates, then finitely many
  integral lattices.”

The project independently reconstructs these steps, catches/repairs a
quasi-Lanner gap and support-stratum ambiguity, replaces published numerical
transition extrema with exact norm arguments, and gives a cleaner
projective-order/conductor treatment.  That is a noncircular proof, but it is
not a technique-independent proof.  The honest phrase is:

> Nikulin's conclusions are not imported; major parts of his method are
> reconstructed, corrected, and extended.

### Vinberg

Vinberg is foundational, not absent.  Used are:

- the quadratic-form structural reduction for reflection arithmetic;
- the cycle-field/arithmeticity criterion and definiteness at nonselected
  embeddings;
- standard reflection/Gram geometry; and
- the compact dimension theorem (`n<30`).

The proof does not use a separate Vinberg arithmetic ground-field or
reflective-lattice finiteness conclusion.  Without Vinberg's arithmetic
criterion, the AG arrangement would have no integral-point condition to
exploit.

### Long--Maclachlan--Reid

LMR Theorem 1.1 is imported and essential for the dimension-two and polygon
branches.  It is not merely an algebraization theorem: its proof uses the
congruence spectral gap, Zograf's inequality, and finiteness of arithmetic
Fuchsian groups of bounded area.  Therefore the polygon branch is genuinely
curve-mediated but not a new AG finiteness proof.

### Khovanskii--Prokhorov

The universal finite-volume ceiling `n<=995` is imported.  The result is
polytope/Coxeter geometry and is logically independent of the ground-field
argument.  Its face-average ingredients admit a polytope-Chow/Hodge
interpretation, but for irrational polytopes that is a combinatorial Hodge
object, not automatically the Chow ring of an actual projective toric
variety.  Calling this step “AG” without that qualification would overstate
the proof.

## What is genuinely new here, with the safe novelty claim

Subject to ordinary literature/peer-review caution, the following are new
**in this reconstruction**:

1. the exact real-framed AGCat functor and literal return on morphisms;
2. projective-trace realization of all integral cycle coordinates without
   choosing free root modules;
3. exact independent enumeration/correction of the four-node transition
   supports and the quasi-Lanner branch;
4. exact norm/capacity replacements for the published transition numerics;
5. the compact-`H3` general-fork invariant
   `beta=t12+t13+t23-z-4`, its sign-primitive identity `k=Q(beta)`, and the
   uniform arbitrary-label cutoff;
6. the resulting all-dimensional field theorem assembled without ABSW or a
   pre-existing global field theorem; and
7. the projective full-order/conductor sandwich closing native hosts.

Items 2, 3, 4, 6, and 7 are best called reconstructions, corrections, or new
implementations of classical architecture unless a systematic literature
search establishes priority.  Item 5 is the clearest candidate for genuinely
new mathematics.  The global theorem's achievement is therefore not
“finiteness from AG with no Vinberg/Nikulin,” but rather:

> an exact algebraic-arrangement factorization, a new compact-`H3` field
> contraction, and an independently audited reconstruction of the remaining
> classical geometric/arithmetic bridges, avoiding the ABSW and Nikulin
> finiteness conclusions.

## Computational consequences and present limits

The factorization is more computational than ABSW's qualitative volume
proof, but the current artifacts are certificates, not yet a practical
classification engine.

### Already finite/combinatorial

- The four-node support search is exact: `32` raw masks give `7` marked
  support orbits grouped into Nikulin's `5` schematic families; the label
  grammar replay has `73` categorical marked label orbits.
- Tail labels satisfy the explicit but crude cutoff `d<20000`.
- General-fork labels satisfy the explicit but very crude cutoff
  `m,n<10^6`.
- Once a field, rank, frame representation, and conductor are fixed, ideal
  classes and submodules of `M/cM` give an actual finite enumeration problem.

### Not yet effective enough for classification

1. The LMR polygon input is used qualitatively; the finite list of maximal
   genus-zero arithmetic Fuchsian groups has not been materialized in this
   workspace.
2. Relative Fekete gives degree/point finiteness, but the complete field and
   algebraic-integer lists have not been enumerated.
3. The `62`-frame and `995`-dimension constants make naive global search
   astronomical.
4. The conductor/class-group proof is algorithmic in principle but no
   uniform practical conductor bounds or host enumerator are implemented.
5. The existing tests verify determinant identities, inequalities, support
   coverage, and module examples; they do not enumerate all chambers.

### Realistic computational program

1. **Field stage:** enumerate cyclotomic bases and contracting algebraic
   integers, using exact interval arithmetic, norm pruning, and LLL; replace
   the million-label rectangle by analytic and congruence pruning before any
   brute force.
2. **Polygon stage:** make the LMR maximal genus-zero list explicit, or use
   explicit triangle/quadrilateral classifications where appropriate.
3. **Frame stage:** for each surviving `(k,m)`, enumerate connected frame
   graphs and cycle packets subject to inertia, Coxeter-angle, determinant,
   cusp, and nonselected-definiteness constraints.  SAT/SMT can enumerate
   support; exact number-field linear algebra certifies candidates.
4. **Host stage:** reconstruct Cartan matrices by spanning-tree gauge,
   enumerate `Lambda`-stable projective lattices via ideal classes and
   `M/cM`, then run a certified Vinberg algorithm to test reflectivity and
   recover the maximal chamber.
5. **Deduplication:** quotient by field automorphisms, lattice similarity,
   diagram automorphisms, and chamber isometry; certify finite volume from
   vertex links.

This route can classify bounded slices and plausibly find examples long
before the global theoretical constants become practical.  New examples are
most realistically sought at small degree/rank by solving exact Gram and
order constraints, not by scanning all ranks up to `996`.

### Improving the dimension ceiling

The global `995` is a geometric import, so better field enumeration alone
does not lower it.  Two observations focus the problem:

- compact chambers already satisfy Vinberg's much stronger `n<30`;
- every noncompact arithmetic input in this quadratic-form model has
  `k=Q`.

Therefore the only range relevant to improving the global upper bound is the
rational, noncompact branch.  A computational attack should combine affine
cusp-link grammar, forbidden elliptic/quasi-Lanner subdiagrams, face-incidence
inequalities, and exact Gram-rank/signature feasibility.  A proof-quality
lower ceiling would require turning the exhaustive infeasibility search into
finite certificates.  Independently, sharpening the connected-frame
constant `62` and the fork/tail label bounds would make classification much
more practical, but would not by itself improve `n<=995`.

## Recommended one-paragraph wording

> The proof is not ABSW in disguise: it uses neither the ABSW finiteness
> theorem nor its conformal-volume/arithmetic-volume mechanism.  But it is
> also not independent of classical reflection theory.  Vinberg supplies the
> quadratic-form and cyclic-integrality framework; Nikulin's narrow-face,
> minimality-14, transition, and bounded-frame architecture is reconstructed
> rather than invoked as a finiteness black box; Long--Maclachlan--Reid closes
> the curve branch; and Khovanskii--Prokhorov supplies the universal dimension
> ceiling.  The genuinely new core is the exact real-framed arrangement
> factorization and, especially, the compact-`H3` determinant/Fekete field
> contraction.  Thus the result is an AG-arrangement-mediated reconstruction
> with new arithmetic geometry, not a purely AG replacement of
> Vinberg--Nikulin.


# Arithmetic dimension bounds and an effective high-dimensional search

**Date:** 2026-08-12  
**Status:** provenance and reduction theorem; no new numerical nonexistence
claim is made.

**Search-direction recheck, 2026-10-13.** Esselmann's original publisher
abstract gives the maximal reflective integral Lorentzian lattice rank as
22, meaning hyperbolic dimension 21. The compact arithmetic case is a
different category and remains subject to Vinberg's bound 29. The precise
noncompact range, including the exclusion of dimension 20, is restated in
[Belolipetsky, Section 1](https://arxiv.org/pdf/1506.03111), with Borcherds'
21-dimensional example described in Section 7. These are search restrictions,
not new results of this repository. The current H22 geometric branch is
nonarithmetic if it ever achieves finite volume; its present volume is
infinite. No rational-lattice search in the excluded dimensions is warranted.

**Current-frontier cross-check.** The introduction of Guo–Ma–Zang–Zheng's
[July 2026 preprint](https://arxiv.org/pdf/2607.14715), *Several families of
incommensurable noncompact hyperbolic Coxeter polytopes*, still records known
compact examples only through dimension 8 and noncompact examples through
dimension 21. Its new constructions are families in dimensions 4, 5, 6, 7
and 9; they do not raise either dimension record. This supplies a recent
primary-source check rather than relying only on the 2015 survey for the
state of known examples.

**Dimension 20 is also a target.** The known dimensions are not simply every
dimension up to 21. [Felikson's research overview](https://www.maths.dur.ac.uk/users/anna.felikson/Polytopes/polytopes.html)
lists finite-volume examples in dimensions at most 19 and in dimension 21,
leaving 20 absent. Thus a nonarithmetic finite-volume H20 example would also
fill a new dimension in the sources checked, although it would not raise the
maximum-dimension record. Esselmann's arithmetic exclusion in dimension 20
does not prove unrestricted nonexistence there. The existing geometric
H19-to-H20 lane should be judged by its specific construction obstructions,
not discarded under the arithmetic theorem.

## 1. The four dimension statements must be kept separate

| scope | strongest bound used or immediately available | provenance | sharpness / open range |
|---|---:|---|---|
| arbitrary finite-volume hyperbolic Coxeter polytopes | `n <= 995` | Khovanskii--Prokhorov for the noncompact case, together with Vinberg for the compact case | geometry-only; intentionally ignores arithmetic |
| arbitrary compact hyperbolic Coxeter polytopes | `n <= 29` | Vinberg's 1981/1984 nonexistence theorem | examples are known through `n=8`; for each `n=9,...,29`, the unrestricted existence problem remains open |
| arithmetic finite-volume reflection groups, compact or noncompact | `n <= 29` | Vinberg's arithmetic nonexistence theorem | no maximality hypothesis is needed |
| noncompact arithmetic reflection groups in the quadratic/lattice model | `n <= 19` or `n=21` | cusp isotropy gives `k=Q`; Esselmann's theorem then bounds reflective integral Lorentzian lattices | the exceptional dimension `21` is realized by Borcherds, so the numerical ceiling `n<=21` is sharp; `n=20` is impossible |

Consequently the displayed `n<=995` in the current AGCat composite theorem is
not the best arithmetic bound.  It is the bound one gets while importing only a
theorem about **all** finite-volume Coxeter polytopes.  If Vinberg's arithmetic
dimension theorem is admitted, `995` can immediately be replaced by `29`.  If
the chamber is noncompact, Esselmann improves this further to the disconnected
range `n<=19` or `n=21`.

The compact and noncompact arithmetic cases behave differently.  For a
noncompact arithmetic chamber a cusp supplies a nonzero isotropic vector over
the ground field.  At every nonidentity real embedding the admissible form is
definite, so the ground field must be `Q`.  After scaling, one obtains a
reflective integral Lorentzian lattice, and Esselmann applies.  This argument
does not apply to compact chambers: in dimensions at least four they cannot be
defined over `Q`, since the rational Lorentzian form would have at least five
variables and hence be isotropic by Meyer.  For a genuine totally real ground
field, the nonidentity definite embeddings defeat the rational-isotropy
argument.

The strongest safe statement about the unknown compact range is therefore:

> No compact hyperbolic Coxeter polytope is presently known in dimension at
> least nine, arithmetic or otherwise, and for each dimension `9,...,29` the
> unrestricted existence/nonexistence problem remains unresolved.  This is not
> a nonexistence theorem.

There are many conditional or bounded-facet nonexistence results inside this
range, but they do not settle arbitrary compact Coxeter polytopes.

## 2. What the AG/arrangement theorem does and does not contribute

The exact carrier and field funnel make the classification problem more
discrete, but they do not presently prove a rank bound.  Their contribution is

```text
quadric--hyperplane arrangement
      -> integral closed Cartan-cycle coordinates
      -> finite ground-field set
      -> finite bounded-frame representations for fixed rank
      -> finite projective native hosts.
```

The missing arrow for a new dimension theorem is still

```text
integral arrangement data -> a rank-independent contradiction above N.
```

Thus a claim that AG has already lowered `29` would be false.  The present
dimension bounds are geometric/arithmetic imports.  The arrangement theorem
does, however, turn a proposed improvement into an exact finite computation in
each of the finitely many remaining ranks.

## 3. A rigorous finite reduction for dimensions 22 through 29

Vinberg's refined high-dimensional theorem makes the first target much smaller
than the global field theorem alone suggests.  For an arithmetic reflection
group in dimension `n>=22`, the complete high-dimensional ground-field list is

\[
 \mathbf Q,\qquad \mathbf Q(\sqrt2),\qquad
 \mathbf Q(\sqrt5),\qquad \mathbf Q(\cos(2\pi/7)).       \tag{3.1}
\]

The rational field cannot occur in the **compact** case here.  Arithmetic
reflection groups are of simplest type; in dimension at least four their
defining rational Lorentzian form has at least five variables and is isotropic
by Meyer, so the corresponding arithmetic lattice is noncocompact.  Thus a
compact arithmetic reflection group with `n>=22` has ground field in precisely
the three-field candidate list

\[
 \mathbf Q(\sqrt2),\qquad
 \mathbf Q(\sqrt5),\qquad
 \mathbf Q(\cos(2\pi/7)).                              \tag{3.2}
\]

Fix one of these fields and `22<=n<=29`, and put `m=n+1`.  Pass from any
finite-covolume arithmetic reflection subgroup to its full integral reflection
host.  Compactness is preserved under this finite-index enlargement.  The
bounded connected-frame lemma supplies `m` chamber walls spanning the
Lorentzian space and satisfying, for every frame edge,

\[
 c_{ij}<62^2                                                   \tag{3.3}
\]

at the distinguished embedding.  For every simple directed cycle `C` of
length `r`, including length-two cycles,

\[
 b_C\in\mathcal O_k,\qquad
 |b_C|<62^r,\qquad
 |\sigma(b_C)|\le 2^r\quad(\sigma\ne\sigma_0).          \tag{3.4}
\]

The Minkowski embedding therefore gives an **explicit finite box** of possible
values for every `b_C`.  There are finitely many connected graphs on `m`
vertices.  Equality of all simple cycle products determines the ordered
reflection representation up to diagonal gauge.

There is a small effectivity point here which must not be hidden behind an
assumed invariant lattice.  It can be resolved directly from the same cycle
data.  Write `T_i=1-s_i`.  The trace of a product
`T_{i_1}\cdots T_{i_r}` is the product along its closed directed walk.  If it
is nonzero, that walk decomposes into simple directed cycles (and loops, for
which `a_ii=2`), so its trace is a product of integral `b_C` and powers of
`2`.  Expanding words in the `s_i` therefore shows

\[
             \operatorname {tr}(w)\in\mathcal O_k
             \quad\hbox{for every word }w\hbox{ in the }s_i.       \tag{3.5}
\]

Burnside's theorem says that the word matrices span
`A=End_k(V)`.  Enumerate words until `m^2` of them,
`e_1,\ldots,e_{m^2}`, form a `k`-basis, and compute the trace-dual lattice

\[
 D=\{x\in A:\operatorname {tr}(xe_j)\in\mathcal O_k
                    \text{ for every }j\}.                         \tag{3.6}
\]

The trace pairing on `A` is nondegenerate, so `D` is a computable full
`O_k`-lattice.  Equation (3.5) gives

\[
 \sum_j\mathcal O_k e_j
       \ \subseteq\ \Lambda:=\mathcal O_k[s_1,\ldots,s_m]
       \ \subseteq\ D.                                             \tag{3.7}
\]

In particular, `D/(sum_j O_k e_j)` is finite.  Starting with the span of
`1,s_1,\ldots,s_m,e_1,\ldots,e_{m^2}` and repeatedly adjoining all left
products by the `s_i` computes `Lambda`: the ascending chain stays inside
`D` and must stabilize in that finite quotient.  Choose any full seed lattice
`L_b` in `V` and set `L_0=Lambda L_b`; this is now a **computed** full
projective `Lambda`-stable lattice, not an oracle input.

Put `B=End_Ok(L_0)`.  The annihilator of the finite module `B/Lambda`
computes a nonzero conductor ideal `c` with `cB subset Lambda`.  Standard
exact number-field algorithms compute ideal-class representatives, and the
remaining lattices are obtained by enumerating the submodules of the finite
modules `M/cM` that are stable under `Lambda`.  Thus the full-order conductor
sandwich really gives a finite algorithm for the stable projective lattices,
rather than only a qualitative finiteness statement.

Finally, the effective termination theorem of Belolipetsky--Kapovich turns
Vinberg's semi-algorithm into a decision procedure for each resulting
arithmetic lattice.  Exact cone and vertex-link calculations decide finite
volume and compactness.

This proves the following reduction, but not its hoped-for verdict.

> **Finite high-dimensional decision reduction.**  The existence and
> nonexistence of compact arithmetic reflection groups in dimensions
> `22,...,29`, and the classification of the resulting full integral reflection
> hosts (which suffices for the reflection-maximal chambers in this reduction),
> reduce to a finite exact enumeration over the three fields (3.2).  A complete empty
> enumeration would improve the dimension theorem for all arithmetic
> reflection groups from `n<=29` to the sharp numerical ceiling `n<=21`.

This does **not** assert a finite classification of all nonmaximal reflection
subgroups inside a surviving host.  Such subgroups are not what the host
enumeration or the Belolipetsky--Kapovich termination theorem classifies, and
there can be infinite reflection-subgroup constructions.  They are irrelevant
to the proposed nonexistence decision: every arithmetic reflection subgroup
has a full host, and compactness survives the finite-index enlargement.

This reduction is stronger and more concrete than merely saying that the
global set is finite.  It identifies eight ranks and three fields and gives a
certificate format.  It is still likely enormous: at `m=30` there are
`2^435` labelled support graphs before structural pruning, and the bounds in
(3.4) are deliberately crude.

For dimensions `14,...,21`, Vinberg similarly restricts the compact ground
field to a short explicit list of quadratic, biquadratic, and real cyclotomic
fields.  Dimensions `9,...,13` can still be made finite by the global field
funnel, but its finite field list must first be materialized rather than left
as a Fekete finiteness assertion.

## 4. Proof-quality computation and certificates

A credible implementation should enumerate constraints rather than raw graphs.

1. **Cycle/frame stage.**  Incrementally enumerate connected Cartan supports
   modulo graph automorphisms and spanning-tree gauge.  Enforce exact Coxeter
   values for intersecting walls, symmetrizability cycle identities, selected
   signature `(n,1)`, positive definiteness at the other embeddings, and the
   bounds (3.3)--(3.4).  SAT/SMT may propose supports, but every accepted or
   rejected branch must be checked in exact number-field arithmetic.
2. **Arithmetic-host stage.**  Construct the reflection order, enumerate its
   stable projective lattices through ideal classes and conductor submodules,
   and quotient by field, lattice, and diagram automorphisms.
3. **Vinberg stage.**  Run the effective Vinberg algorithm for each host.
   Certify root ordering, completeness, exact Gram inertia, all finite and ideal
   vertex links, and absence of hyperideal extreme rays.  Test compactness and
   deduplicate the resulting full mirror chambers.
4. **Coverage certificate.**  Freeze a manifest for every branch with one exact
   rejection reason: non-Coxeter entry, failed symmetrization, wrong inertia,
   failed nonselected definiteness, no stable lattice, nonreflectivity,
   noncompactness, or duplication.  A small independent verifier should replay
   the finite-box coverage and every terminal witness without floating-point
   tolerances.

For finding examples, run the same pipeline in the opposite order: prioritize
small discriminant, small two-cycle values, sparse connected supports, and
small conductor.  A survivor is then completed by Vinberg's algorithm.  This
is much more plausible than scanning facet lists directly.

## 5. Where small ridges help, and where they do not

Bogachev proves that a compact Coxeter polyhedron in `H^3` has a ridge of width
less than `5.75`; the same bound holds in higher dimension when the polytope
has a three-dimensional face which is itself Coxeter.  This improves Nikulin's
general width bound `7` dramatically in the affected classifications (the paper
reports, in its rank-four examples, reductions from roughly `280` to at most
`85` candidates over `Q`, and from roughly `160` to `83` over
`Z[sqrt(2)]`).

This is excellent pruning for low-rank classification and for any
higher-dimensional branch in which a Coxeter three-face has been certified.
It is not currently a universal dimension argument: a face of a Coxeter
polytope need not itself be a Coxeter polytope, and no theorem used here forces
the required three-face in every compact high-dimensional chamber.

A promising way to obtain a **new numerical** bound is instead to revisit
Vinberg's weighted two-face proof.  Its combinatorial proposition assigns
weights to planar angles and requires **both** a vertex inequality (the total
weight at every vertex is at most `cn`) and, for every `k`-gonal two-face, a
face inequality (the total weight around it is at least `5-k`).  Under those
global hypotheses its conclusion has the form

\[
                         n<8c+6.                         \tag{5.1}
\]

Vinberg constructs a globally valid assignment with `c=3`, giving `n<30`.
Improving this is not merely the calculation of a smaller local maximum.  In
the high-dimensional arithmetic range one may enumerate candidate local
configurations over the restricted fields using the integral cycle
coordinates, but one must then exhibit and verify an assignment satisfying
the vertex and face inequalities simultaneously for every possible chamber.
This can be organized as a finite exact linear-programming/covering
certificate only after exhaustive coverage of the configurations has itself
been proved.  For example, such a global certificate with `c<=23/8` gives
`n<29` and therefore rules out dimension `29`.  Each further decrease gives a
rigorous incremental improvement.  This target is potentially far smaller
than classifying every full chamber while still retaining an auditable finite
certificate.

The obstruction is precise: field finiteness and small ridge bounds alone do
not decrease `c`; one must prove exhaustive coverage of all configurations
that can occur and a globally valid weight assignment satisfying both families
of inequalities.

## 6. Primary references checked

- E. B. Vinberg, *Absence of crystallographic groups of reflections in
  Lobachevskii spaces of large dimension*, Functional Analysis and Its
  Applications **15** (1981), 128--130,
  <https://doi.org/10.1007/BF01082285>.  This is the source for nonexistence of
  compact Coxeter and arithmetic reflection groups in dimensions `n>=30`.
- E. B. Vinberg, *Absence of crystallographic groups of reflections in
  Lobachevskii spaces of large dimension*, Transactions of the Moscow
  Mathematical Society **47** (1984), 75--112.  The expanded paper contains
  the weighted-face argument and the high-dimensional ground-field restrictions.
- F. Esselmann, *Ueber die maximale Dimension von Lorentz--Gittern mit
  coendlicher Spiegelungsgruppe*, Journal of Number Theory **61** (1996),
  103--144, <https://doi.org/10.1006/jnth.1996.0141>.
- N. V. Bogachev, *From geometry to arithmetic of compact hyperbolic Coxeter
  polytopes*, European Journal of Mathematics **8** (2022), 1271--1293,
  <https://arxiv.org/abs/2003.11944>.
- M. Belolipetsky and M. Kapovich, *Effective bounds for Vinberg's algorithm
  for arithmetic hyperbolic lattices*, Sao Paulo Journal of Mathematical
  Sciences **18** (2024), 558--568,
  <https://arxiv.org/abs/2206.01282>.
- N. Bogachev, S. Douba, and J. Raimbault, *Infinitely many commensurability
  classes of compact Coxeter polyhedra in `H^4` and `H^5`*, Advances in
  Mathematics **452** (2024), 109826.  Its introduction records that dimension
  eight remains the highest dimension in which a compact Coxeter polytope has
  been exhibited.

## 7. Bottom line

The present AGCat proof should be described as computationally enabling, not as
already computationally classifying.  It makes fields, frame representations,
orders, and lattices finite and exact.  Vinberg supplies the current compact
ceiling and the best high-dimensional field restriction; Esselmann completely
settles the noncompact dimension range.  The most defensible next theorem target
is an exact certified empty search in dimensions `22,...,29`, or an improved
local-weight constant in (5.1).  Either would be a real new dimension result;
merely enumerating more fields would not.

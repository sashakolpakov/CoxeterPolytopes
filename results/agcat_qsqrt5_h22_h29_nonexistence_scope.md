# What AGCat would—and would not—prove for `Q(sqrt(5))` in dimensions 22–29

## Verdict

The exact AGCat functor makes a rigorous nonexistence proof possible, but it
does not by itself prove nonexistence.

Let `C_5,22:29` be the full subgroupoid of compact arithmetic Coxeter
chambers in dimensions 22 through 29 with minimal Vinberg field
`Q(sqrt(5))`.  The project has proved a faithful functor

```
F_arr : Cox_arith^mark -> AGCat_arr^R
```

and a return functor on its image with

```
R_arr F_arr(P) = P.
```

Consequently

```
C_5,22:29 is empty  <=>  F_arr(C_5,22:29) is empty.
```

This is the useful transfer statement.  It means that an **exhaustive**
empty computation on the corresponding integral arrangement image is a
theorem about the original polytopes, not a heuristic lattice search.

For bare nonexistence, faithfulness and exact return are stronger than
necessary: a total functor from a nonempty source category to an empty target
category cannot exist.  The real issue is therefore not faithfulness.  It is
whether the relevant target slice is independently empty and whether the
functor is defined on every hypothetical source chamber.  The universal
arrangement carrier is total but its target slice is not independently empty;
the full-`H^2` K3 carrier has the desired rank obstruction but is not known to
be total on raw chambers.

## Why the existing theorem is not yet that computation

The global AGCat theorem proves that the reflection-maximal arithmetic
chambers form a finite set.  Its proof produces existentially finite field,
cycle, order, and lattice stages.  It does not yet materialize all of them as
finite files.  The project’s computational architecture lists five remaining
enumeration layers:

1. the complete field manifest;
2. all Minkowski cycle-packet boxes for every field/rank;
3. the general-field ideal-valued Cartan gauge;
4. all full orders, conductors, ideal classes, and stable quotient modules;
5. complete, independently replayed Vinberg transcripts for every retained
   native host.

For the fixed field `Q(sqrt(5))` and ranks 23 through 30, the first layer is
already fixed rather than existential.  The remaining four layers are still
finite, but they have not all been enumerated.  Existing `D_n`, `D24+`,
Leech, and other host searches cover individual native lattices; their empty
prefixes do not exhaust the conductor-lattice stage.

One genuine part of the second layer is now complete.  Restoring all dotted
edges, the full local nonorthogonal-triangle census starts with
54,158,527,500 unordered edge multisets.  Exact element-squareclass joining
leaves 4,768,141 integral golden cycle packets; 2,633,885 pass conjugate
positive-definiteness, and all of those are Lorentzian at the selected
embedding.  An independent implementation reproduces every count, all 17
strata, and both packet-stream digests.  This is strong local pruning, not a
rank-23--30 obstruction: a connected frame need not contain a complete
nonorthogonal triangle, and overlapping triangles and longer cycles have not
yet been made globally compatible.

The first overlap stage is now complete too.  Modulo exact four-vertex
symmetry, every connected finite/orthogonal packet and every connected packet
with exactly one dotted edge has been checked simultaneously on all triangles
and Hamilton four-cycles and at both real embeddings.  These two tranches
retain 119 and 180,402 packets.  This is genuine four-node compatibility,
but it does not exhaust cyclic packets with at least two dotted edges or
assemble every complete connected frame of ranks 23 through 30.  A no-import
hostile replay independently reproduces the orbit counts, all 48 sign strata,
and all four compatible/retained stream hashes.

A further bounded tranche closes every connected four-node **tree** with
exactly two dotted edges.  Of 471,865,955 canonical assignments, 97,172,716
pass the two-place Gram gate (25,023,281 stars, 81 disjoint-dotted paths, and
72,149,354 adjacent-dotted paths).
This large surviving set reinforces rather than resolves the global issue:
cyclic packets with two or more dotted edges and all full-rank joins remain.

The census nevertheless supplies an exact full-rank **family** in every
target dimension.  Take a path in native rank \(r=23,\ldots,30\), put one
terminal dotted product \(c\) on it, and label every later edge 3.  Its
leading determinants are

\[
D_j=2j-(j-1)c.
\]

The selected form has index one for every \(c>4\), while the conjugate form
is positive definite precisely when \(c'<2r/(r-1)\).  The complete dotted
alphabet contains respectively
3,590, 3,584, 3,577, 3,571, 3,566, 3,561, 3,556, and 3,552 such products;
\(c=2+2\varphi\) works uniformly.  The integral Cartan gauge
\((-1,-c)\) on the first edge and \((-1,-1)\) thereafter gives a genuine
arithmetic Lorentzian lattice.  Conjugate definiteness makes the ambient
arithmetic orthogonal lattice cocompact, and the dotted cyclic product makes
its minimal Vinberg field exactly \(\mathbf Q(\sqrt5)\).

These are new full-rank cocompact arithmetic **host/frame seeds**, not new
compact reflection chambers: the displayed path cone is not finite-volume,
and no complete Vinberg computation has shown its reflection subgroup to
have finite index.  Their abundance proves that frame compatibility alone
cannot be the golden nonexistence obstruction.

## The categorical distinction

Faithfulness and exact return on the image imply that no source chamber is
lost or identified with a different chamber.  They do **not** imply that:

- every abstract real-framed quadric arrangement is in the Coxeter image;
- the image is empty in high dimension;
- the target category is bounded by a K3 Picard-rank constraint; or
- a finite set supplied existentially has already been enumerated.

In category-theoretic language, the proved functor is an exact embedding of
the source into its image.  It is not an essential-surjectivity theorem onto
a separately classified empty category.

The stronger full-`H^2` K3 carrier does give the conditional quadratic bound
`n<=8`, but its carrier is not universal: essential surjectivity for raw
arithmetic chambers is false/open in the audited scope.  That conditional
K3 rank bound therefore cannot be applied to all hypothetical golden-field
chambers in dimensions 22–29.

## Correct nonexistence theorem to target

The exact implication now worth implementing is

```
complete Q(sqrt(5)), ranks 23–30 cycle packets
  -> complete ideal gauges and reflection orders
  -> complete conductor/stable-lattice list
  -> complete Vinberg hosts
  -> no compact survivor
  => no compact arithmetic Coxeter polytope over Q(sqrt(5)) in H22–H29.
```

Passing to a full native host does not lose compact existence.  For the stable
lattice supplied by the fixed-data theorem,
`W(P) <= R(L) <= O^+(L)` has finite indices.  If `W(P)` is cocompact, then its
finite-index reflection overgroup `R(L)` is cocompact and has the same
commensurability field and dimension.  Thus an empty enumeration of **all**
full native hosts returned by the conductor-lattice stage excludes nonmaximal
compact examples too.  This formulation does not require assuming that an
abstract maximal reflection overgroup exists or is unique.

The present status is therefore: **theorem architecture complete; finite
golden slice not yet enumerated empty**.

The first new direct golden host launched alongside this audit is
`<-a> + S23` in `H23`, with `S23` the determinant-three odd-Leech section.
For `a=phi` and `a=2phi`, the complete first packets contain 553,104 and
53,682 walls, respectively.  Exact strict hyperideal witnesses survive both
packets.  This expands the searched golden host set but, because later
Vinberg packets remain open, it is neither a compact survivor nor a field-wide
nonexistence certificate.

## Authoritative project records

- `results/agcat_global_maximal_arithmetic_coxeter_finiteness.md`
- `results/agcat_global_maximal_arithmetic_coxeter_finiteness_independent_hostile_audit.md`
- `results/computational_agcat_coxeter_classification_architecture.md`
- `results/canonical_chamber_four_property_completion_global_agcat_addendum.md`
- `results/compact_h22_h29_field_map_20260812.md`
- `results/golden_agcat_h22_h29_decision_manifest.md`
- `results/golden_agcat_full_triangle_packets.md`
- `results/golden_agcat_full_triangle_packets_independent_hostile_audit.md`
- `results/golden_agcat_four_node_at_most_one_dotted.md`
- `results/golden_agcat_four_node_at_most_one_dotted_hostile_audit.md`
- `results/golden_agcat_two_dotted_four_node_trees.md`
- `results/golden_agcat_two_dotted_four_node_trees_independent_hostile_audit.md`
- `results/golden_odd_leech_section_h23_first_packets.md`
- `results/golden_odd_leech_section_h23_first_packets_independent_hostile_audit.md`

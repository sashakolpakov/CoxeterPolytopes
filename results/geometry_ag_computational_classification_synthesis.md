# Geometry, AG, and a computational Coxeter-classification program

**Date:** 2026-08-12  
**Status:** synthesis of independently audited provenance, executable pilot
censuses, certified constant improvements, and dimension-search reductions.

## 1. Short verdict

The global proof is neither a purely algebraic-geometric proof nor ABSW in a
new notation.  It is an **AG-mediated reconstruction of a
Vinberg--Nikulin-style arithmetic-geometric argument**.

Its indispensable gates are:

| gate | dominant mathematics | what it contributes |
|---|---|---|
| exact return | projective quadric arrangements + real hyperbolic framing | a faithful functor from marked chambers; no finiteness by itself |
| finite fields | narrow-face geometry + invariant cycle coordinates + integrality, norms, and relative Fekete | the genuinely AG/arithmetic-geometry-mediated finiteness gate |
| finite ranks | Coxeter/polytope geometry | the dimension ceiling |
| finite hosts | bounded root frames + arithmetic orders, conductors, and projective modules | finitely many native lattices for fixed field and rank |
| maximal chamber | elementary reflection geometry | one full host and one unmarked chamber orbit |

ABSW's finiteness theorem and its conformal-volume/spectral-volume mechanism
are not used.  A structural result cited as ABSW Lemma 5.2 is Vinberg's
quadratic-form lemma and has no finiteness content.

Nikulin's final fixed-data/global-field finiteness theorems and the published
bound `N(14)<=25` are not imported.  However, the proof substantially uses or
reconstructs Nikulin's methods: narrow faces, minimality 14, transition
packets, field inheritance, Fekete contraction, the connected frame constant
62, and the fixed-data strategy.  The honest novelty claim is therefore
reconstruction and extension, not technique independence.

Vinberg is foundational: the quadratic/lattice-form model, cyclic
arithmeticity and nonselected definiteness, reflection/Gram geometry, and the
compact dimension theorem all enter.  The clearest new mathematical step is
the compact-`H3` determinant coordinate and its sign-primitive/Fekete field
contraction.

## 2. What is computational now

The implemented module
[`cartan_packets.py`](../src/coxeter_polytopes/classification/cartan_packets.py)
does exact mathematics rather than checking theorem prose.  In a finite
rational edge-product box it:

1. enumerates all connected labelled supports;
2. constructs every closed Cartan-cycle packet;
3. rejects nonrational or nonintegral cycles;
4. reconstructs the Cartan matrix up to diagonal root gauge;
5. constructs an integral gauge using prime-valuation difference constraints;
6. quotients every facet permutation;
7. checks exact Lorentzian inertia and all spherical/affine simplex links; and
8. emits a chamber certificate which a separate geometry engine replays.

For the rational alphabet

\[
             b_{iji}\in\{0,1,2,3,4\},
\]

corresponding to Coxeter labels `2,3,4,6,infinity`, the frozen censuses are:

| stage | rank 3 / hyperbolic plane | rank 4 / hyperbolic 3-space |
|---|---:|---:|
| labelled assignments | 125 | 15,625 |
| connected assignments | 112 | 15,104 |
| integral rational labelled packets | 68 | 2,976 |
| unlabelled packets | 18 | 207 |
| Lorentzian packets | 13 | 200 |
| finite-volume simplices | **13** | **18** |
| compact / noncompact | 4 / 9 | 1 / 17 |
| independent chamber replays | 13 / 13 | 18 / 18 |

These are exhaustive results in the displayed box.  They are not yet claims
of new examples: maximality and comparison with a literature database require
the native-host and certified Vinberg-completion stages.

The executable evidence and its independent clean-room audit are described in
[`computational_agcat_coxeter_classification_architecture.md`](computational_agcat_coxeter_classification_architecture.md)
and
[`computational_cartan_classifier_independent_hostile_audit.md`](computational_cartan_classifier_independent_hostile_audit.md).

## 3. New effective cutoffs

The norm arguments can be made much smaller before any Fekete enumeration.

* The compact-`H3` general-fork label box improves from
  `m,n<1,000,000` to
  \[
                         \max(m,n)\le 420.
  \]
  Exactly 525 unordered pairs pass the necessary norm inequality; the only
  pair at the endpoint is `(3,420)`.
* The transition-tail label improves from `d<20,000` to
  \[
                              d\le126.
  \]
  Exactly 66 labels pass the necessary norm inequality.

These are exact integer sieves and were independently reimplemented.  A
survivor is only a candidate: it need not support an integral determinant
point, admissible field, or Coxeter chamber.  See
[`coxeter_finiteness_constant_audit_and_h3_cutoff.md`](coxeter_finiteness_constant_audit_and_h3_cutoff.md)
and its
[`independent audit`](coxeter_finiteness_constant_improvements_independent_hostile_audit.md).

## 4. Dimension bounds and a finite high-rank experiment

Four different statements must not be conflated.

1. All finite-volume Coxeter polytopes satisfy the geometry-only ceiling
   `n<=995` used in the AGCat composite theorem.
2. Importing Vinberg's stronger arithmetic theorem immediately gives
   `n<=29` for every arithmetic reflection group, without maximality.
3. In the noncompact arithmetic branch, cusp isotropy forces `k=Q`; Esselmann
   then gives `n<=19` or `n=21`.  The dimension-21 Borcherds example shows
   that the numerical ceiling 21 is sharp, while dimension 20 is excluded.
4. Compact Coxeter examples are known only through dimension 8; unrestricted
   compact existence remains open in each dimension 9 through 29.

The first proof-quality high-rank target is dimensions 22 through 29.  Vinberg
restricts all arithmetic ground fields there to

\[
 \mathbf Q,\quad\mathbf Q(\sqrt2),\quad\mathbf Q(\sqrt5),\quad
 \mathbf Q(\cos(2\pi/7)).
\]

Meyer excludes `Q` for compact forms in these dimensions.  Hence compact
existence reduces to eight ranks and three fields.  Bounded connected cycle
packets, a computable trace-dual reflection order, conductor/class-group
lattice enumeration, and effective Vinberg termination make this a finite
exact decision problem for full integral hosts/reflection-maximal chambers.
Emptying all eight ranks would improve the global arithmetic ceiling from 29
to the sharp value 21.  A survivor would instead raise the known lower bound
for compact examples beyond 8.

The proof and its effectivity qualifications are in
[`arithmetic_dimension_bounds_and_effective_high_dimension_search.md`](arithmetic_dimension_bounds_and_effective_high_dimension_search.md).

## 5. Practical next implementation

The shortest path to a genuine classification/new-example engine is:

1. materialize the finite field manifest, beginning with the 525 fork pairs,
   66 tail labels, and the explicit polygon branch;
2. add a `Q(sqrt(5))` backend for Minkowski boxes, real-embedding signs, and
   prime-ideal/ideal-class gauges;
3. enumerate connected cycle-basis packets before full graphs, with exact
   inertia and nonselected-definiteness pruning;
4. compute reflection orders and conductor-stable projective lattices;
5. run certified Vinberg completion and compare canonical full-host chambers;
6. for a dimension theorem, either empty ranks 22--29 or improve Vinberg's
   global weighted-two-face constant with an independently replayable finite
   linear-programming certificate.

Bogachev's small-ridge computations demonstrate that this style can cut real
classification searches sharply: in rank-four examples, approximately
`280 -> 85` candidates over `Q` and `160 -> 83` over `Z[sqrt(2)]`.  That
bound is universal in `H3` and applies in higher dimensions only when a
Coxeter three-face has been certified; it is pruning, not by itself a new
global dimension theorem.

## 6. Reproduction

```bash
PYTHONPATH=.:src python3 scripts/classify_bounded_cartan_packets.py --rank 3 --summary-only
PYTHONPATH=.:src python3 scripts/classify_bounded_cartan_packets.py --rank 4 --summary-only
PYTHONPATH=.:src python3 -m pytest -q tests/test_computational_cartan_classifier.py
python3 -m scripts.certify_coxeter_finiteness_constant_improvements
python3 -m pytest -q tests/test_coxeter_finiteness_constant_improvements.py
```

The detailed provenance ledger is
[`geometry_vs_ag_provenance_ledger.md`](geometry_vs_ag_provenance_ledger.md).

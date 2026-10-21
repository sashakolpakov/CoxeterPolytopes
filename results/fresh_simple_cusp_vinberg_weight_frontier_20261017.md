# Vinberg's compact scope and a 3n affine weight interface

The checked primary literature does **not** supply a mixed-simple
finite-volume dimension bound29 by direct citation. There is a concrete
extension candidate with the correct vertex budget: retain Vinberg's
finite-vertex weights and use the rule below at simple ideal vertices.
Its total vertex cost is at most3n, and the triangle payment transfers.
**Quadrilateral payment remains unproved.** Thus this package gives no
new ambient dimension bound; the frozen all-simple bound61 is unchanged.

The [checker](../scripts/fresh_simple_cusp_vinberg_weight_frontier_20261017.py)
and [data](fresh_simple_cusp_vinberg_weight_frontier_20261017.json) verify
the graph counts exactly. Their finite ranges are confirmations of the
all-rank formulas below, not cutoff-based proofs.

## What the primary theorem says

Vinberg's original Theorem1 on printed p.69 explicitly concerns compact
type. Its §6.1 begins with a bounded fundamental polytope, and §6.2
calculates costs for elliptic vertex diagrams. The original §4 is the
coarse distance7 argument; §§5--6 provide the much finer weights giving
the compact bound29. These are distinct proof steps.
[Original paper](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt).

One modern wording that could suggest a broader scope occurs at the
opening of Burcroff's §9. But its following Proposition9.2 assumes a
compact Coxeter polytope, and Theorem10.2 again states the bound for the
compact class. That introductory wording does not provide a proof for
mixed finite/simple-ideal vertices.
[Primary preprint, §§9--10](https://arxiv.org/pdf/2201.03437).

The primary Felikson--Tumarkin introduction likewise attributes29 to the
compact theorem and separately proves8 for the simple all-ideal class.
[Primary paper](https://arxiv.org/pdf/math/0502413).
This is a bounded source check, not a claim to have established the best
published bound for the mixed-simple class.

## The proposed ideal-vertex rule has cost exactly3n

Assume n>=14. A simple ideal vertex has a connected affine diagram on n
nodes. Exceptional affine types have at most9 nodes, so only the four
classical underlying shapes occur: a cycle, a path, a single fork, or a
double fork. Normalize every ideal corner weight to0 or1 as follows:

* weight1 if the two marked roots have graph distance at most3;
* also weight1 if they belong to an end segment of order at most6;
* otherwise weight0.

On a fork an end segment must contain both interchangeable fork leaves
or neither. The maximal allowed fork-end segment contains the two leaves
and the first four vertices of the central path. A cycle has no end
segments. Distances ignore bond multiplicity, just as in the preceding
[distance proofs](fresh_simple_cusp_affine_exclusion_20261017.md).

The exact costs are:

| Underlying affine graph | Pairs at distance<=3 | Additional end pairs | Total |
|---|---:|---:|---:|
| cycle |3n|0|3n|
| path |3n-6|6|3n|
| single fork |3n-5|5|3n|
| double fork |3n-4|4|3n|

For a path, each end contributes two extra distance4 pairs and one
distance5 pair. A fork end contributes two distance4 pairs and no
distance5 pair. These extra pairs are disjoint for n>=14. The short-pair
counts follow directly from the all-rank path/cycle/fork formulas already
proved in the frozen distance package. This proves the table in every
rank. The checker reconstructs all348 graphs of orders14 through100 and
counts their literal pair sets independently.

At finite vertices use Vinberg's published weights, including the four
E8 exceptions. Their total is at most3n by the
[frozen finite-component check](fresh_compact_dimension_bound_20261015.md)
and its inherited all-rank proof. The new affine rule only concerns ideal
vertices; it does not reassign weights to finite forests.

## Triangle payment does transfer

By the [global affine exclusion](fresh_simple_cusp_affine_exclusion_20261017.md),
each triangle star in the all-simple class, n>=14, contains a Lanner
diagram L of order at most5 containing its three black roots. Follow the
literal cases in Vinberg §6.4, rather than importing the whole compact
face theorem.

If L has order3 or is a4-cycle, there are two adjacent black pairs, so
two corners receive weight1. If L contains a bond5, the proof chooses a
corner whose connected diagram contains that bond; such a corner cannot
be a classical affine n-node diagram. It is finite of type H3 or H4 and
receives weight2 under the unchanged compact rule.

The remaining case is a5-cycle. Vinberg's two displayed marked patterns
give two corner contributions of weight1 by adjacency and his local
endpoint/special-edge lemmas. If either contributing corner is ideal,
the required marked pair is connected within L minus the omitted black
root, on at most4 vertices, so its distance is at most3. Our ideal rule
therefore gives at least the same unit payment. The finite corners retain
the original local lemma. Hence every triangle receives at least2.

The classification and local finite-corner lemmas here are explicit
published inputs. The graph checker does not claim to reclassify Lanner
diagrams or independently replay their angle arithmetic.

## The precise remaining quadrilateral obligation

The absence of affine diagrams of order below n does not make the
whole corner diagram elliptic. A quadrilateral can have a connected
affine n-node corner. In particular an affine A cycle has two marked-to-
marked arcs. Vinberg's §6.5 description of W_ij as a tree with the two
marked roots at its ends does not apply verbatim to that corner. His
later classification of each T_ij as a finite classical diagram or E8
also needs new cases. Thus the compact quadrilateral theorem has not
been transferred by the cost calculation above.

The exact target is now clear: prove that every actual quadrilateral
in the all-simple class has at least one corner paid by this mixed rule,
or give a literal realizable counterexample. Any proof must retain the
full original star, all finite/affine vertex conditions, negative index1,
corank1, and the global exclusion of affine subdiagrams of order below n.
No label cutoff or bounded white-tail scan suffices.

Two tempting shortcuts have specific defects:

1. Taking the maximum compact weight over all single-root deletions of
   an affine A cycle marks exactly every pair at cyclic distance at most5.
   Its cost is5n for n>=14, not3n. The checker verifies all27 sizes14..40;
   the proof is that a pair lies in a6-node end segment of a deleted cycle
   precisely when one of its cyclic connecting arcs has length at most5.
2. Deleting a supporting white root does not preserve a quadrilateral
   star's corank1. For a symmetric N-by-N Gram matrix of rank N-1 with
   kernel vector k, its adjugate is a nonzero scalar multiple of kk^t.
   Deleting a node w with k_w!=0 therefore leaves a nonsingular principal
   matrix. Vinberg §6.10 explicitly uses det H=0; a spherical deletion
   cannot silently be substituted for a full quadrilateral star.

These are obstructions to particular transfer arguments, not
counterexamples to the proposed mixed rule or to a possible bound29.

## Status

Generation and a fresh saved replay pass, covering348 affine graph checks
and27 cycle-deletion envelopes. Syntax, all nine direct source hashes,
and report links pass. An internal mathematical/source read independently
confirms the triangle transfer, including the original5-cycle case.
The certificate status explicitly leaves quadrilateral payment open and assigns no new
ambient dimension bound. The older61/127/163 theorems and all frozen
predecessors are untouched.

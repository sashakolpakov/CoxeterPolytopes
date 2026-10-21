# Independent review of the cross-coface cap capacities

**PASS.** The geometric capacity and the selected literal mappings in the
[primary report](fresh_sqrt5_cross_coface_cap_capacity_20261017.md) are sound.
The [independent checker](../scripts/fresh_sqrt5_cross_coface_cap_capacity_20261017_review.py)
and [saved certificate](fresh_sqrt5_cross_coface_cap_capacity_20261017_review.json)
reconstruct all 241 selected occurrences and capacity columns. They recover
70 positive loads and 53 violations of the unchanged saved masses.
No optimizer is invoked and no dimension conclusion follows from this review.

## Why the capacity is one

Fix an actual compact maximal simplex face S in the compact Coxeter polytope
and normalize all its outward facet normals to square 2. Its full star has
n+1 independent normals: a nonsingular connected Lorentz core and an
orthogonal positive definite spectator forest. Thus the star spans the
entire ambient Lorentz space. This full-star statement, and the canonical
prism/base incidence interpretation, are inherited from the frozen inputs.

In each selected raw prism star, the only roots outside that core are the
deleted cap b and one root q. After deleting b, q is an isolated A1. All
further roots in the full dimension-13 realization are appended orthogonal
spectators of the raw star; in particular they are orthogonal to b, not merely
to the core. This last fact is essential: without it, a cap could have
additional spectator components and the argument would not apply.

The complete ordered contacts with the core determine its projection v
uniquely. With the specified squared norm-2 product tau between b and q,
acuteness fixes the negative sign. Consequently

\[
b=v-\frac{\sqrt\tau}{2}q,
\qquad v^2=2-\frac\tau2\in(-1,0).
\]

The core has trivial full weighted-graph automorphism group, so its
identification supplies a single contact profile, not several independently
chosen placements. If the same profile occurs at distinct isolated A1
copies q and q', the copies and the core are mutually orthogonal. The two
caps therefore have product v², strictly between -1 and 0. This is forbidden
for distinct norm-2 facets of a Coxeter polytope: a nonzero negative finite
product has magnitude at least 1, while parallel or dotted products have
magnitude at least 2.

For the same q, the spanning full star and prescribed contacts determine
the normalized cap uniquely. The fixed canonical base-vertex subset B
determines one face of S. At any vertex of B the inherited vertex/coface
capacity permits only one selected maximal prism base in that immediate
coface. Equivalently, adjoining this unique b to the full S-star specifies
the same canonical n+2-root prism aggregate. It cannot be counted a second
time as another occurrence. This addresses same-copy repetitions as well
as different-copy incompatibility.

Hence at most one selected incidence occurs per actual S. Summing gives
the signed inequality with positive prism-incidence coefficients and
negative simplex coefficient 1. Omitting other positive types only weakens
it. Different profile columns are separate valid inequalities; their
right sides are not asserted to be disjoint resources.

## Independent mapping and arithmetic checks

The checker imports no producer or project arithmetic. Its only arithmetic
is standard-library rational arithmetic in Q(sqrt(5)). It does the following:

* Reconstructs each raw complete weighted graph and checks the saved core
  map edge by edge. Removing b leaves exactly the mapped core and the sole
  isolated q, with no additional internal component.
* Reconstructs each base's vertices from the union of the canonical K2,2
  partitions. A vertex labelled p is the corner omitting b and p; it is
  **not** a retained black-root label. The canonical image agrees with every
  saved `canonical_vertices` list, without duplicates.
* Checks the full spectator rank, every isolated-A1 multiplicity, and the
  provenance-to-source-core and incidence-to-simplex mappings.
* Exhausts weighted adjacency-preserving automorphisms separately for all
  13 used core types. Every group is trivial. It independently computes
  every positive definite one-root deletion, as well as a negative
  determinant; together these confirm the required nonsingular Lorentz
  signature.
* Computes every selected projection norm directly by inverse matrices.
  Cores with only 3/5 bonds use their literal norm-2 matrices. The remaining
  cores are trees: along a parent-child edge with squared product w, set
  s_child=s_parent*w, diagonal 2*s_i, and edge -s_child. This is a positive
  diagonal congruence of the physical Gram. These cases have one pi/3 cap
  contact at i, so its projection norm is s_i times the corresponding
  inverse diagonal. This avoids adjoining square roots and independently
  recovers 2-tau/2 in every occurrence, including the three cases with two
  pi/3 contacts, which use literal matrices.
* Independently regroups the selected occurrences, counts incidence
  multiplicities once, and recomputes all exact column loads, slacks and
  required simplex mass lower bounds.

The particular nine-copy witness also has a separate
[independent reconstruction](review_sqrt5_cross_coface_full_base_witness_20261017.md).
This family review does not rely on that checker's arithmetic.

## Verification and limits

Generation and fresh saved replay pass, as do syntax and direct source-hash
checks. The review binds the primary Python/JSON and all five literal input
JSONs; it does not bind the primary prose, so adding a review link there
cannot change this certificate.

This is an independent check of every selected capacity, its relevant
core geometry and saved diagnostic loads. It does not rerun the original
arithmetic/local-realization census, establish simultaneous realization of
the saved face masses, enumerate omitted cap types, or certify a repaired
dual or primal solution. No mathematical or mapping correction was found.

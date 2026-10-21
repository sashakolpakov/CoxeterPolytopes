# Two more redundant pairs in the order14 affine B factor

In affine B13, whose diagram has14 nodes, delete the two pairs(10,12)
and(10,13). Its cost falls from86 to84. No finite weight or coupon changes,
and no other affine type changes. The order14 maximum is now84, shared
by B and C, so its normalized rank13 coefficient becomes6.

The [source](../scripts/fresh_order14_affine_b_fork_projection_transfer_20261017.py)
and [certificate](fresh_order14_affine_b_fork_projection_transfer_20261017.json)
check all23 affected selected contexts and all52 affected induced maps.
This is a complete local payment refinement, with no ambient dimension
claim. Further ideal and finite obstructions remain in ambient66.

## The intrinsic pair rule

Number the spine0--1--...--11 with doubled bond10--11 and attach the
two short leaves12,13 to0. The changed pairs have distance11 and form
one orbit under the diagram's automorphisms. Start from the completed
[singleton deletion](fresh_b14_b15_singleton_completion_transfer_20261017.md),
which already gives cost86. Keep its five previously absent pairs absent;
remove only the additional two pairs above.

All82 basic ideal5/feature5 pairs remain at unit weight. In particular,
every elementary triangle payment and every marked path of length at
most8 remains a unit. No fractional replacement or new finite allowance
is being used.

## Complete support and necessary contact domains

The complete distance8 small-factor certificate enumerates every induced
map of every eligible selected connected marked component. Exactly23
contexts contain one of the two changed marked images:

    624,629,637,674,713,719,736,778,842,
    1555,1559,1567,1572,1584,1641,1647,1667,
    1848,1855,1934,2241,2248,2327.

For each context the new checker reconstructs the native labelled pattern
from its saved full core Gram, reruns the exhaustive induced-map routine,
and compares its complete output with the frozen map certificate. There
are52 maps with a changed marked image. Unaffected maps in the same
context retain their previous proof.

Every affected core H is nonsingular Lorentzian, all old roots have
square2, and all its nonzero off-diagonal entries are-1. The imported
constraint checker is the full one from the preceding singleton proof:
it recomputes the exact inverse J, the actual preferred corner and its
complete necessary contact domain, its partition into two end sets,
and the entrywise-positive inverse blocks on those ends. It also checks
that a root meeting both ends yields a retained marked path of length
at most8. These are the complete contact-domain statements of the
frozen overlap proof, not bounds inferred only from the selected pattern.

Thus, unless a retained short path already pays the face, any additional
supporting white has all its old-core contacts on one end. If a is its
nonnegative vector of native contact magnitudes, then

    projection_H(e)^2 = a^T J a <= e^2 = 2.

The complement of the Lorentzian span of H is positive definite. Within
one end, positivity of every entry of J means that a forced contact q
at node i implies

    a^T J a >= q^2 J_ii,

regardless of any further contacts on that end. Contacts on the other
end are already resolved by the short-path alternative. No unspecified
contact is assigned zero, and no bound on spectator rank, multiplicity
or length of indirect chains is imposed.

## Every actual factor forces an impossible new doubled terminal

In every one of the52 maps, target node10 belongs to the old selected
component and is marked, while target node11 is missing. To be the actual
whole affine B13 factor, the corner must supply a root at node11 joined
by label4 to node10.

That root is genuinely outside H. Old roots omitted from the preferred
corner are the two other black roots and cannot belong to that corner's
factor. Every old corner root outside the selected connected component
is orthogonal to it, so none can supply the required nonzero bond to
the old node10. This argument remains valid if other target nodes lie
in old components later joined by indirect paths: it uses only the
forced first terminal contact, not a completed Gram matrix.

With roots normalized to square2, a label4 bond has squared contact
magnitude2. The complete exact obstruction table is:

|Old inverse diagonal J_ii|Forced projection2*J_ii|Maps|
|---:|---:|---:|
|2|4|12|
|7/2|7|12|
|8|16|20|
|26|52|8|

Every projection lower bound exceeds2. Therefore no such actual affine
completion exists in an otherwise unpaid face. The checker verifies the
forced terminal, old marked neighbor, native label and exact inequality
for each individual map. This excludes the affected ideal occurrence;
it does not reinterpret its proper finite subdiagram as a whole finite
component.

## Coverage beyond the selected maps

The original complete theorem first pays triangles and marked paths of
length at most8; these units are retained. Its low-core domain has whole
marked-component bound at most12. Its disjoint high-core branch is
excluded geometrically. Protected/basic ideal outcomes have distance at
most9. Finite whole-component, adaptive and protected-leaf payments do
not change. The remaining affected selected cases are precisely the
23 contexts checked above.

The universal tree and D13/D14 immediate affine outputs have order at
most13, so they do not introduce an order14 B obligation. The later
cycle12 proof returns finite payments or ideal paths of length at most4;
the D11 proof returns finite payments or paths of length at most5. The
order12 B half-transfer and C deletion have their stated order12 outputs
and finite alternatives. They never call for either changed B13 pair.
The finite rank17 deletions used in the completed ceiling66 change no
ideal weight. The older same-order partial-weight replacements use
affine D factors, not these B pairs.

Consequently every triangle still receives at least2 and every
quadrilateral at least1. The only change to cusp costs is B13:86 to84.
The order14 A/B/C/D costs are70,84,84,82. Actual cusp-fiber multiplicities
remain those of the preceding theorem.

## Scope and verification

Generationdca80c and separate fresh replay68c6b2 both terminated
exit0/PASS on all52 maps,23 necessary contact domains, exact projection
bounds, retained elementary/basic pairs and automorphism invariance.
A separate full source/report read, including the inherited constraint
and map-witness helpers, passed the norm-two Lorentz projection semantics,
complete end-contact partition, positive inverse blocks, retained
cross-end short path, genuinely new doubled-terminal witness, arbitrary
additional same-end contacts, complete saved-map comparison and later
smaller-factor/finite interfaces. That reviewer did not duplicate the
52-map arithmetic or the upstream geometric census.

Root also completed a full source/report read, including the inherited
constraint and map-witness helpers. That review passed the norm-two
new-terminal argument, exhaustive affected-pattern map reconstruction,
end positivity with arbitrary additional contacts, and all preserved
interfaces. Root did not repeat the map arithmetic or geometric census.
All31 direct source bindings, Python AST and three local report links
pass. Source, certificate and report are final/frozen.

The complete all-label core census, no-extra-cross-edge proof and
arbitrary-support contact-domain geometry remain inherited frozen inputs;
their original searches are not rerun here.

The previous duplicate-rank12 concentration is invalid after this cost
change. A separate complete scalar diagnostic is required; the improved
five-rank13 profile alone does not give an ambient bound. No new
polytope, external refereeing or formalization is asserted.

    python3 scripts/fresh_order14_affine_b_fork_projection_transfer_20261017.py
    python3 scripts/fresh_order14_affine_b_fork_projection_transfer_20261017.py --verify

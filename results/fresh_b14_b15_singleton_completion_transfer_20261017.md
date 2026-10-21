# Two redundant singleton payments in affine B factors

Remove the canonical pair `(0,10)` in affine B13 of diagram order14,
and `(1,11)` in affine B14 of diagram order15. Their costs fall from87
to86 and from97 to96, respectively. No finite coupon or replacement half
is needed: every selected-core occurrence of either pair would require
a norm-two supporting root whose forced projection has square greater
than2. The required actual affine factor is therefore impossible in an
otherwise unpaid face.

The [source](../scripts/fresh_b14_b15_singleton_completion_transfer_20261017.py)
and [certificate](fresh_b14_b15_singleton_completion_transfer_20261017.json)
check the complete saved support, all64 affected induced maps, and their
exact necessary projection inequalities. This package makes no ambient
dimension assertion. It changes only these two singleton weights; other
possible reductions remain outside its scope.

## The affected rules and complete map support

For order14 the B diagram has spine0–1–…–11, a doubled terminal bond10–11,
and two short leaves12,13 at0. For order15 it has spine0–1–…–12, a doubled
terminal bond11–12, and short leaves13,14 at0. The proposed singleton in
each diagram has distance10 and is fixed by every diagram automorphism.

The order14 rule starts with its87 retained pairs in the frozen
[complete small-factor map certificate](fresh_distance8_ideal5_small_factor_20261017.md).
The order15 rule starts with the97-pair
[previous B14 reduction](fresh_order15_affine_b_pair_reduction_20261017.md),
so its three previously deleted pairs keep their existing proofs.
The new source checks all elementary triangle pairs, every pair at distance
at most8, and every basic ideal5 pair remain paid.

|Target diagram|Deleted pair|Affected selected contexts|Affected maps|New cost|
|---|---|---:|---:|---:|
|Affine B13, order14|0,10|23|46|86|
|Affine B14, order15|1,11|11|18|96|

Five contexts occur for both targets, so these are29 distinct old-core
contexts and34 target/context occurrences. The complete context lists and
every marked map are saved. The map counts include only maps with the
changed marked image; maps of the same unmarked pattern to other pairs
are unchanged. This avoids treating every map in an affected context as
an affected payment.

## A necessary projection bound with arbitrary additional contacts

The [complete overlap theorem](fresh_distance8_high_overlap_20261017.md)
applies after the retained short8 alternatives have been handled. Its
cross-core exclusions fix the literal old core H. In all contexts used
here, every old core root has square2, every nonzero off-diagonal entry
is−1, and H has signature `(order(H)-1,1,0)`. The checker reconstructs
its exact inverse J and its actual preferred corner. In particular it
checks the necessary core-contact domain for any additional supporting
white root, using the corner's finite inverse diagonals or affine kernels.

The saved contact domain splits into two end sets. Both diagonal inverse
blocks are entrywise strictly positive. A further white root meeting both
ends would give an adjacent marked path of length at most8. The source
recomputes the distance bound; its actual maximum in these cases is6.
Such a path already resolves the face, so in an otherwise unpaid face
all contacts of a new white e with H lie in one end set.

Let a be the nonnegative vector of contact magnitudes of e with H. Since
H is nonsingular and Lorentzian, its orthogonal complement in the ambient
Lorentz space is positive definite. Thus

`projection_H(e)^2 = a^T J a <= e^2 = 2`.

If the actual target factor forces a contact of magnitude q at old node i,
entrywise positivity on that end gives

`a^T J a >= q^2 J_ii`.

This remains a necessary inequality when e has further contacts anywhere
on the same end. Contacts on the opposite end are already resolved by the
short-path alternative. No unspecified contact is set to zero, and there
is no restriction on the number, labels or indirect chains of additional
supporting whites.

All projection comparisons use the norm-two normalization before any
diagonal congruence: a label4 bond has magnitude sqrt(2), hence square2.
The new source does not introduce a rational gauge or change the squared
norm of e. Consequently the obstruction for such a bond is `2 J_ii>2`.

## Every map forces the offending doubled terminal

Let K be the displayed connected marked component of the preferred core
corner, and suppose it lies in an actual affine B factor F of the indicated
order. Its embedding into the fixed target is one of the saved induced
maps. The target's doubled terminal is absent from K, while its adjacent
node is in K and is one of the two marked nodes. This is checked for
every one of the64 maps. The required terminal is therefore a supporting
white e with a forced label4 contact to that old marked node i.

This e must be a genuinely new root outside H. An old root absent from
the preferred corner is one of its two omitted black roots and cannot
belong to F. Every other old corner root outside K belongs to a different
connected component and is orthogonal to K. None can be the terminal,
which must have a nonzero bond to i. The source verifies these old corner
components and their zeros against K. This elementary adjacency argument
identifies the witness even if other nodes of F were to lie in old
components reached by indirect paths; no no-bridge lemma or completion
of the whole factor is needed.

For every map the exact old full-core inverse diagonal satisfies J_ii>1.
The complete obstruction values are:

|Target order|J_ii|Forced projection lower bound2 J_ii|Number of maps|
|---:|---:|---:|---:|
|14|5/4|5/2|2|
|14|2|4|8|
|14|7/2|7|8|
|14|8|16|16|
|14|43/2|43|8|
|14|26|52|4|
|15|5/4|5/2|2|
|15|2|4|6|
|15|7/2|7|8|
|15|8|16|2|

Every lower bound is strictly greater than the root's square2. Thus the
extra terminal cannot exist unless a retained short payment already
resolved the face. This excludes the selected affine-factor occurrence
itself; it is not an assertion that a finite displayed subdiagram has
become a whole finite component.

## Coverage beyond the selected maps

The [complete finite8/ideal5 theorem](fresh_finite8_ideal5_combined_payments_20261017.md)
and its successors first cover triangles and marked paths of length at
most8. All their required elementary B pairs retain unit weight here.
Universal cycle rules are unchanged. The universal tree and D-transfer
immediate whole-affine alternatives have order at most13, so they impose
no new order14 or15 obligation. Their existing successors also remain
unchanged.

After the short-path reduction, low-core configurations have whole marked
component order at most12, and disjoint high cores are excluded by their
all-label theorem. Finite whole-component, protected-leaf and adaptive
finite payments keep their rules. Protected or basic ideal outcomes have
marked distance at most9 and retain unit weight. The remaining high-core
overlap cases are precisely those covered by the complete small-factor
map certificate. It lists the23 and11 affected contexts here, and every
affected marked map has the obstruction just proved.

Every triangle therefore still receives at least2 and every quadrilateral
at least1. The finite budget is unchanged, including all existing halves
and units. The intrinsic ideal pair costs become86 and96. Other types
are unchanged, so the order14 maximum is86, and the order15 maximum is96
(shared with its unchanged D rule). A global dimension comparison is a
separate task.

## Verification scope

Primary generation86321 and the separate fresh saved replay (terminal
command7d085b) both completed exit0/PASS. Syntax, all16 direct source bindings
and all6 report links passed. A separate exact SymPy calculation, without
importing producer arithmetic, independently reconstructed the inverses of
all13 distinct old core Grams, checked their Lorentz signatures using
characteristic-polynomial Sturm counts, and checked all64 forced projection
inequalities across the34 target/context occurrences.
The checker validates every direct input binding, reconstructs each
literal Lorentzian core and preferred corner, verifies the exact inverse
end blocks and contact-domain coverage, and checks every affected induced
map and terminal projection inequality. The complete overlap and map
censuses remain frozen inputs; they are not independently re-enumerated.
No external refereeing, formalization or new dimension claim is made.

A separate full internal source/report mathematical read passed. It checked
the genuinely new terminal witness, the norm-two doubled-contact bound,
entrywise-positive inverse blocks with arbitrary additional contacts, the
retained short-path alternative and the complete local-domain composition.
It did not repeat the frozen overlap or induced-map censuses. Source, data
and report are final/frozen.

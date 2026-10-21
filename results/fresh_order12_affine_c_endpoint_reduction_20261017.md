# Removing the endpoint pair in the12-node affine C factor

The finite8/ideal5 face-payment rule of the completed
[dimension69 proof](fresh_ideal5_weighted_dimension69_20261017.md) remains
valid after removing just the two-terminal pair of an affine C11 factor
of order12. In the canonical numbering this is `(0,11)`, at distance11.
Every other one of its65 pairs retains weight1. Triangles still receive
at least2 and quadrilaterals at least1. The change introduces no finite
coupon and leaves the existing finite budget unchanged.

This is a local composition theorem. It does not by itself lower the
maximum over all order12 affine families or the ambient dimension bound.
In particular, the A/B order12 obligations and any coupons introduced
by their separate transfers require their own proofs and joint budget.

The [source](../scripts/fresh_order12_affine_c_endpoint_reduction_20261017.py)
and [certificate](fresh_order12_affine_c_endpoint_reduction_20261017.json)
make the complete provenance check below. They use the frozen
[order12 support diagnostic](fresh_order12_payment_support_20261017.md)
and reconstruct its C-target embeddings with a separate standard-library
path algorithm. No producer arithmetic module is imported.

## The invariant pair and the elementary payments

The canonical diagram is the path0--1--...--11 with label4 at0--1 and
10--11, and label3 elsewhere. Its two automorphisms are identity and
reversal, so the unordered endpoint pair is intrinsically defined.
Scaling the norm2 Gram by the positive vector with endpoint entries
1/sqrt2 and all other entries1 gives the path Laplacian. This verifies
the stated affine type, including its one-dimensional positive kernel.

All65 pairs other than the endpoints are retained. The basic ideal5/R5
rule requires63 pairs; its only additional retained pairs are
`(0,10)` and `(1,11)`, both at distance10. Every pair of distance at most8
is basic and therefore remains paid.

The [complete triangle theorem](fresh_ideal5_finite8_triangle_payments_20261017.md)
uses ideal distance at most5, or distance at most8 with a mark within5
of a feature. Its whole exceptional factors have order at most9.
Thus it never needs the removed distance11 pair. All other factor rules
and finite corner payments are unchanged, so its two triangle payments
remain valid, including at nonsimple cusps.

## The short8-first quadrilateral reduction

Use the complete
[finite8/ideal5 composition](fresh_finite8_ideal5_combined_payments_20261017.md).
It first treats any adjacent marked path of length at most8. A finite
corner pays by the unchanged finite rule. An ideal C11 corner of order12
pays directly, since all its short8 pairs remain. The other affine
families retain their old short-path rules and universal transfers.

Those universal transfers also remain valid when their *conclusions*
are examined. The cycle replacements at orders13,14 and at least15 use
finite whole-component, finite short-path or finite coupon payments;
their ideal alternatives are an adjacent black bond or a shared
spectator path of length2. Their classical-tree alternatives have marked
distance at most9 with a persistent nearby feature. None can call for
the C12 endpoint pair of distance11. This includes the protected
finite-leaf9 argument: its whole-component sealing is a finite outcome,
not an implicit all-small-affine payment.

It is therefore legitimate to pass to a star having no marked path of
length at most8, just as in the complete predecessor. The choice of
minimal indefinite opposite cores, closed-opposite-edge argument and
actual-edge positive definiteness are inherited unchanged. There is no
new simplicity, compactness or arithmetic hypothesis.

## Complete low and high core provenance

The complete low-core domain has66 chosen-corner closures and three
optional sealed simple-root extensions, giving69 actual component
variants. These are actual closed components, rather than finite
patterns allowed to sit in an arbitrary larger affine diagram. The
frozen exact signatures and positive kernels give12 affine variants
of order12: two B and ten D. There is no C variant. All finite outcomes
remain paid by the whole-component rule. Consequently no low-core
branch loses a payment.

The complete high-domain predecessor excludes disjoint cores and treats
all2,614 surviving overlaps, with every additional cross-core edge either
a short path or excluded by its saved quadratic certificate. For a
possible C factor of order12, every relevant selected marked component
must have size at most12 and a complete component bound at least12.
The diagnostic includes **all758** such patterns. A smaller bound rules
out the proposed whole factor; a larger selected component cannot embed
in it. Finite, protected and basic outcomes are unchanged.

The new source reconstructs all induced weighted embeddings of these
758 connected patterns into the canonical C path. A graph embedded in
a path must itself be a path, and its image is a contiguous segment;
the algorithm tries every segment and both orientations and checks every
edge and nonedge. Hence this is a complete embedding search, not a
sampling of the saved successful maps.

Exactly15 contexts have embeddings, with30 maps in total:

    148,159,164,575,602,620,1081,1328,1352,1450,1452,1549,1630,1695,2573.

Their marked images are precisely

    (0,9), (0,10), (1,10), (1,11), (2,11).

All are retained. Only164 and1452 need the nonbasic distance10 orbit.
The source compares every reconstructed map and marked image with the
frozen diagnostic, including the absence of a map in every other
pattern. Thus the removal is supported by the complete actual-star
classification and closure arguments, rather than mere absence from a
partial map list.

## Immediate alternatives and every later transfer

The current rule also reaches some small factors indirectly. These must
be checked separately from the selected high/low patterns.

The [D13 long-pair transfer](fresh_affine_d13_long_pair_transfer_20261017.md)
has six marked immediate alternatives on whole order12 factors. Its
`inner4` alternative is B and its `inner_next3` alternative is D. All
six induced maps and their positive-kernel provenance are bound by the
diagnostic; none is C. The corresponding
[D14 alternatives](fresh_affine_d14_long_pair_transfer_20261017.md)
have whole factor order13. The positive affine kernel makes these whole
factors, so they cannot conceal a C factor of order12.

The later reductions used in the dimension69 rule have these outputs:

|Transfer|Payments other than its original selected ideal pair|
|---|---|
|D13/D14 universal long pairs|The immediate B/D factors just checked; short paths; sealed finite units or halves|
|D12 opposite-leaf transfer|Finite units/full A1--D12 coupons/finite halves, or ideal paths of length at most4|
|B12 singleton and its later short-leaf upgrade|Finite A1--D12 or A2--B10 coupons/units; short paths; the old immediate B13 branch is excluded by a strict Schur gap|
|B14 three-pair reduction|Whole finite D14, A12 or A13 components|
|D16 singleton removal|Finite A4--D13 neighbor halves, or a shared-neighbor path of length2|
|Order17 half-weight reduction|Finite halves plus the retained half at the same order17 D factor|
|Order16 half-weight reduction|An impossible completion or finite halves plus the retained half at the same order16 D factor|
|Latest B13/B14 singleton completions|Projection contradictions, with no new payment|

Here finite type names use rank and affine orders count nodes. The source
binds each precise local theorem and checks the literal terminal payment
interfaces:144 D13/D14 equality configurations,96 D12 configurations,
six B12 singleton/conditional-long cases, seven B14 cases, the D16
singleton's four corners, all saved order16/order17 half-weight completions,
and the latest64 impossible B completion maps. Positive finite payments
are checked against finite corner signatures; positive ideal halves
remain on the stated D factor. The analytic short-path branches and
whole-component sealing are inherited from the cited proofs. The older
18 finite-only D replacements and their all-label selected support are
already checked by the complete predecessor and are unchanged here.

This list covers the successive changes from the complete finite8/ideal5
rule through the frozen dimension69 rule. No branch uses the removed
C12 endpoint as an alternate unit or half-payment. In particular, an
old theorem phrased with a whole-small-factor cutoff is not imported
without inspecting its actual output factors or paths.

Two newer independent local interfaces are also checked for compatibility:
the [cycle12 transfer](fresh_cycle12_ideal5_transfer_20261017.md) has only
finite coupons or ideal paths of length at most4, and the
[D11 leaf transfer](fresh_affine_d11_leaf_pair_transfer_20261017.md) has
only finite payments or ideal paths of length at most5. Their conclusions
are safe for this C deletion. Their new finite coupons are **not** charged
to the unchanged budget claimed here; composing either newer rule still
requires its separate finite-budget proof. No unfinished B12 reduction
is assumed by this package.

The separate [finite A17 endpoint deletion](fresh_finite_a17_endpoint_payment_support_20261017.md)
is compatible as well. It introduces no alternative ideal pair. Its only
debited old literal case662 retains two whole finite D17 units, checked
again here; every later typed terminal unit avoids that finite A17 pair.
Its other selected, short and low routes are the same routes audited
above. Thus combining these two deletions does not compensate one removed
pair by the other. The exact revised finite maximum remains the separate
finite-forest budget's conclusion.

## Conclusion and verification scope

The modified rule pays every triangle and quadrilateral with the same
minimum as the completed predecessor. Its C11 factor cost is65, so that
factor contributes `65V/12` at a cusp with V link vertices, by the
unchanged actual factor-fiber multiplicity. Other affine costs, all
finite payments and their established finite budget are unchanged.
There is no new ambient bound or construction claim.

Final generation and a fresh separate-process saved replay both terminated
exit0/PASS. All61 direct bindings,11 report links and Python syntax passed.
Root completed a full source and mathematical report read, passing the
contiguous-path/reversal search including all rejected patterns, the69
actual low variants, whole B/D immediate factors, all later finite,
short-path and same-factor-half interfaces, elementary payment retention
and explicit A17 case662 compatibility. That review did not independently
repeat the758-pattern embedding arithmetic or the inherited classifications.
Source, data and report are final/frozen.

The new computation reconstructs the C-target embedding search and checks
saved interface arithmetic. It does not repeat the
all-label minimal-core classifications, the2,614 supporting-white
closures, low-core exact inertia calculations or the local transfers'
continuous spectator proofs. Those frozen proofs and certificates are
explicit inputs. No external refereeing or proof-assistant formalization
is claimed.

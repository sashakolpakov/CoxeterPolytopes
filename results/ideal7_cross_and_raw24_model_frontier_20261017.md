# Finite cross-payments, dimension 91, and the enlarged raw24 model

The unrestricted ceiling is now **91** for finite-sided, finite-volume
real hyperbolic Coxeter polytopes. This is a complete internally checked
computer-assisted proof, improving the completed 97 argument. Primary
generation and fresh replay pass; a separate exact implementation agrees
on all 477 curvature profiles; separate internal mathematical reviews find
no gap. External refereeing and proof-assistant formalization have not
been established.

| Class | Current ceiling | Change at this checkpoint |
|---|---:|---|
| Arbitrary finite-volume Coxeter polytopes | 91 | Improved from 97 |
| Compact arithmetic Coxeter polytopes | 21 | Unchanged |
| Arbitrary compact Coxeter polytopes | 29 | Unchanged |
| Finite-volume with every affine cusp factor of rank at least 9 | 29 | Unchanged |

The compact-facet ambient ceiling remains 79. The compact arithmetic
claim is distinct from Esselmann's noncompact arithmetic result and
Borcherds' noncompact example. This checkpoint makes no new construction
or publication-priority claim.

## The dimension proof

The [local payment theorem](fresh_ideal7_finite_cross_payments_20261017.md)
retains finite unit payments for connected marked pairs at distance at
most 9 or in components of order at most 17. It adds half-payments between
whole A1 and D16 components, or between A2 and D14 components, at specified
fork marks. The total finite cost remains at most 9n for n at most 108.
The added payments permit ideal cost at most
`min(15m/2, m(m-1)/2)` for every factor of order m.

The proof preserves payments of at least 2 on every triangle and 1 on
every quadrilateral. Complete induced labelled embeddings of 80 inherited
preferred components isolate nine exceptional stars. Five already have
a finite unit payment; four are covered by two finite half-payments.
An affine kernel and a Lorentzian radical relation show why outside
support cannot invalidate these component types. The separate
[rational half-payment certificate](fresh_h17_half_coupon_transfer_20261017_review.md)
reconstructs all nine stars and 27 finite corners and gives a separate
proof of this support argument. It does not claim to duplicate the
80-pattern embedding census.

The [dimension transfer](fresh_ideal7_cross_dimension91_20261017.md)
uses actual factor multiplicities at nonsimple cusps, then cap-15
convex-order concentration. All 477 reduced profiles for ambient
dimensions 92 through 97 have positive curvature minus cost. The smallest
margin is exactly `10121/4194304`, at ranks `(1,15,15,15,15,15,15)`.
The completed 97 theorem supplies the tail. The
[independent checker](fresh_ideal7_cross_dimension91_20261017_review.md)
uses different shape enumeration, sliding-window convolution and a
cumulative-tail identity, with no project arithmetic imports.

The six-rank-15 profile at ambient dimension 91 has negative margin
`-428955/2097152`. Thus this particular scalar comparison cannot lower
the ceiling to 90. This formal profile is not an existence result.
Further progress requires another payment improvement or global
incidence information; rerunning the same comparison cannot help.

## Context-dependent cusp payments

The [complete H17 allocation calculation](fresh_h17_contextual_cusp_payments_20261017.md)
checks all 3,456 faces of the known pyramid while retaining the earlier
finite payments. When ideal payments can depend on actual corners, its
exact budget region is

    B_D16 >= 109,
    B_E8xE8 >= 612,
    B_D16 + B_E8xE8 >= 723.

All three integral boundary allocations are certified. The smallest
common normalized coefficient for this chamber is `109/17`. This
explains how contextual payments can evade the preceding invariant
136-pair obstruction. It is not a universal coefficient theorem.

The report also proves the general capacitated allocation criterion by
a finite network-flow argument. To turn that criterion into another
dimension reduction, one must bound its face-deficit sums for every
cusp subset of an arbitrary polytope. That geometric bound is open.

## The new arithmetic optimum

One solve of the [enlarged finite model](fresh_mixed_raw24_full_model_20261017.md)
preserves all earlier obligations and adds every exported local, joint,
same-family and mixed raw24 column. It has 435,793 variables, 82,806 face
obligations and 6,558 vertex budgets. All 96,176 zero-budget signed columns
and 52 signed columns with vertex cost retain their coefficients.

Exact primal and dual certificates give

    5755387047017860/462000155393947 <= optimum
        <= 3114386878673/250000000000,

or approximately `[12.4575435307, 12.4575475147]`. The target is 12, so
removing the preceding fixed witness did not solve the enlarged model.
No arithmetic dimension improvement follows.

The [new support diagnostic](fresh_mixed_raw24_full_model_20261017_support.md)
finds 4,647 positive face masses: 2,747 simplex/edge rows and 1,900
quadrilateral rows. Reusing old fixed-mass covers would be invalid.
Retaining all simplex/edge rows but only the 138 quadrilateral rows
already in the mixed export gives an exact restricted ratio below 12.
Retaining instead all 1,504 quadrilateral rows with explicit or unique
actual provenance still gives `12.3950149133 > 12`. Thus ambiguity of
merged rows is not necessary for the present obstruction.

The next concrete actual-corner cohort comprises raw families 178, 296
and 297: 47 positive classes and 376 corner occurrences, with core orders
9 or 10. Their dotted contact roots lie in Q(sqrt(5)), but compatibility
still requires full contact vectors, height tests and actual multiplicity.
This cohort has been identified, not yet enumerated. Complete implicit
local and joint families also need reconstruction at the changed masses.

## Construction results and exact scope

The [old-30 connector scout](fresh_old30_connector_flip_20261017.md)
subdivides the previously blocked B-C edge, then flips the independently
blocked Y-B edge to R-Z. Both new abstract disks have 13 vertices, 31
edges and 19 triangles. Under the stated embedding hypotheses their
conditional open-end change would be minus 5.

The bounded scout has 1,940 coherent six-cap prefixes, but no surviving Z.
The [exact native gate](fresh_old30_connector_native_gate_20261017.md)
then excludes all 295 low-label recipes and both entire high-label
families on the 20 specified algebraic prefixes. It tests both norm
branches and proves 594 strict native-contact rejections. It does not
exclude continuously changed prefixes or other incidences.

The new scout corrects one sentence in the frozen preceding report:
the 20 prefixes do not all have zero A2 profiles. The correct profile
multiplicities are `(0,0):14`, `(0,1):2`, `(0,4):2`, `(1,0):2`.
The predecessor's source and data already use those profiles; their
bytes and conclusions are preserved.

The certified control remains **63 walls, 199 finite vertices and 10
open target-face ends**, with infinite volume. No new finite-volume
example has been found.

## Archive and verification scope

The [archive auditor](../scripts/audit_ideal7_cross_and_raw24_model_frontier_20261017.py)
and its [saved result](ideal7_cross_and_raw24_model_frontier_20261017_audit.json)
extend the [preceding checkpoint](ideal7_and_raw24_frontier_20261017.md)
without changing its frozen files. They bind all nine new source/data/report
packages and preserve the previous recursive dependency graph.

The small independent checks reconstruct the four intrinsic order-17
paid sets, compare all nine exceptional stars against the separate rational
certificate, check 19,684 finite-budget endpoints and every affine linear
branch, and compare the two 477-profile certificates. They also recompute
the known chamber's face deficits and three allocations, check exact
arithmetic-certificate ratios and restrictions, and verify abstract disk
counts and all 594 saved strict rejection inequalities.

These archive checks do not replace the producers' complete embedding
census, full arithmetic matrix verification or outward interval
reconstruction. Those have generation and fresh replay checks and the
scoped mathematical/source reviews recorded in their reports. The
arithmetic model also has a separate narrow review of its append,
scaling and recovery wiring; that review does not rerun the matrix or
solver. All reviews described here are internal.

Final archive generation passes: 1,382 predecessor files preserved, nine
new manifests, 515 cumulative manifests, 1,402 bound sources and all 64
local links across ten reports. The separate fresh replay must reproduce
the saved archive exactly; its terminal result is recorded by the run.

# Source proposal: complete-contact affine-B terminal tail

25 October2026. PROPOSAL ONLY, requiring different-author source review and
root acceptance before integration. No frozen helper, final reader, input,
receipt or manuscript is changed. No import, test, native enumeration,
numerical matrix check, scientific job or build ran for this proposal.
No global bound or executed coverage is claimed.

## The observed failure and correct existing refusal

The preserved current failure is
results/finite57_current_rule_selected_04_coalesced_owned_execution_20261025.stdout.
Its exact full key is affine B/order16, high1668, core
(25, 25, ((5, 9), (6, 8), (7, 7), (8, 6), (9, 5))),
base_nodes=[0,1,4,5,6,7,8,9,10,11,12,14],
selected_marks=[1,12],
induced_map=[8,9,7,6,5,4,3,2,0,14,15,1], marked_image=[9,15].

The missing native nodes are10,11,12,13. The mandatory first node10 has
one simple contact to old target9, original root1. The saved inverse entry
is5/4, giving residual square3/4. The existing two-root certificate is
correct to refuse: node10--11 is simple, whose squared contact1 does not
exceed2*(3/4). The terminal-orthogonality certificate is inapplicable too:
the old actual corner includes another old root outside the mapped base,
and the doubled terminal's immediate predecessor has no direct old entry.

## Exact reusable source hypotheses

The complete existing source contact controls are independent of whether
the hypothesized final factor is B15, B16, B17 or a D target. Use ONLY an
exact original HIGH row/core/base/ordered marks and an explicitly saved
control from either:

1. fresh_d16_partial_weight_transfer_20261017.json, cases[i].constraints;
2. fresh_order15_bd_terminal_coupon_transfer_20261018.json,
   B_cases[i].proof.maps[j].residual_contact_control.

The source rows must name this same high index, core key, method, complete
base and physical Gram/inverse. The source's own HIGH pin must equal the
current frozen HIGH pin. If multiple matching controls exist, they must
agree completely. No old B15 map is asserted to be the new B16/B17 map.

Source pins:
D16 data6d3bbabeba06ea257f4d1e783ae359b1cbf9448287ec8b058969152c827f3ff8;
B15 data dda5e90c98c99b99bdf4b667421adfb2dc1637862e7e70a2b906bc1dd5a178ec;
D16 contact producer4f60793bce5c1c5a91d59eda946bbc71cb3956bbb35fa6bc68b61bca28cc55ee.
The complete different-author argument is e740b5ccf7ed66546ceca59507a3c3fdedd6179069732d12f73ca6dbb37c3555,
Sections2.1--2.4 of audit_finite57_selected_affine_completion_chain_fresh_20261024.md,
accepted at exact finite predicates by root6fe04d21f048a55997b583bc7bc54689d0f79ac994f008526795fb981cbf719f.
The residual/no-bridge report is4a074a5c70634d64581c90c1884c142bd6b560c05c1760c416db712f03062677.

For high1668 the B15 saved control is one_relevant_root_with_optional_exterior_root:
allowed=[0,1,2,3,13], ends=[[0,1,2,3],[13]], relevant node1,
exterior node3, residual matrix[[3/4,-1/2],[-1/2,1]],
joint_tail_bound=1 and joint_component_order_bound=14.
The exterior-root whole-component bound14 is below the ACTUAL target order16.
Its old source field actual_factor_order=15 is not copied as an order16 fact.
The target-independent complete contact theorem supplies the needed bound.

The original physical core is norm2, nondegenerate Lorentzian; its complete
necessary contact set, entrywise-positive inverse ends and actual adjacent
white-interior S8 exits remain mandatory. All two-contact and ordered
residual/Cauchy tests include repeated coordinates. The exact source inverse
identity may reuse the existing hash-keyed GJ cache; saved signature and
actual-space/contact completeness remain declared accepted inputs.

## Narrow given-map domain

This proposal is restricted to a current intrinsic-zero affine B factor of
order N=16 or17. The complete supplied native map must have image exactly

    {0,...,t-1} union {N-2,N-1},

and missing nodes exactly t,t+1,...,N-3, with at least TWO missing nodes.
These are a single terminal spine chain. Native labels are simple except
the final bond N-4--N-3, which is doubled. The fork leaves are already mapped.
Every old native edge and nonedge, full ordered map and physical source
association must be verified. A marked pair or index alone is insufficient.

The first missing root e at t touches the mapped node t-1 by a simple bond.
Let i be its original physical root. Require1<J_ii<=2 and require i to be
a permitted candidate; in the relevant branch it is the unique relevant node.
The saved complete ordinary uniqueness theorem must apply, or the complete
relevant/exterior/no-bridge theorem must have bound STRICTLY below N.

## Newness and old-core orthogonality BEFORE any zeros

The first missing root is new: omitted original blacks do not occur in this
actual corner; every other old corner component is orthogonal to the mapped
base and cannot supply the required first edge.

Under NO-S8 the complete contact theorem restricts e to its one unit contact
at i. Ordinary complete ordered residual tests allow no second old-core
contact root. In the relevant branch only its single exterior root could
coexist. If present, the FULL no-bridge result bounds the marked whole
component by the recorded source bound<N, contrary to the whole B_N target.
This uses exterior multiplicity and all indirect-tail control, not immediate
neighbor counting. Hence that root is absent.

Every later required tail root is now genuinely new and old-core-orthogonal.
If a first later root were an old white in another old corner block, its
preceding new tail root would have another old-core contact, impossible.
The omitted blacks are absent and mapped roots have their distinct fixed
images. Induction therefore excludes every old-root shortcut along the tail.
This is exactly the accepted first-contact argument in e740 Section2.4.

Consequently the old projection of e has square J_ii and its positive-space
residual has square rho=2-J_ii<1. Every later tail root lies in the positive
orthogonal complement of the old Lorentz core. No unspecified product was
set to zero before proving this conclusion.

## Terminal-tail vector: a uniform physical calculation

Write the remaining tail after e as f1,...,fk, where k>=1.
Every root has PHYSICAL norm2. For k>=2, consecutive bonds are simple except
the last f{k-1}--fk bond, which is doubled; all nonconsecutive products vanish.
The e--f1 bond is simple and e is orthogonal to all other tail roots.
For k=1, e--f1 itself is the doubled terminal bond.

Set

    v=f1+...+f{k-1}+fk/sqrt(2).

For k>=2 the diagonal contribution is2(k-1)+1 and consecutive cross
contribution is-2(k-1). Thus v^2=1. For k=1, v^2=2/2=1 as well.
The target incidences give e.v=-1 in BOTH cases.
This uses the coefficient1/sqrt2 on the ACTUAL norm2 terminal, not a
fictitious physical norm1 terminal.

All f_j are old-core-orthogonal, so v is in the positive complement and
e_perp.v=e.v=-1. Positive-space Cauchy gives

    1=(e_perp.v)^2 <= e_perp^2 * v^2 = rho <1,

a contradiction. Arbitrary additional supporting roots cannot alter this
necessary principal-vector calculation. No inertia computation, bounded
tail search, new numerical source receipt or geometric realization search
is needed. The vector is used only after complete old-contact/no-bridge
control; a two-root Cauchy inequality was not assumed to handle a longer tail.

## Scope and prospective application evidence

The same source conditions may apply to the later observed B16/B17 cores,
including the mirrored relevant case1935 and the ordinary first-entry cases
with source residual2/3. They are examples, not a high-index allowlist or
a claim of actual new key coverage. Each runtime application must find its
own matching complete saved contact control, verify its given terminal
map/path and retain every exact source predicate and actual target order.

A prospective output must identify the original full key/physical context,
the saved CONTACT-PROFILE selector (not a false current-target map selector),
the verified current target/native map, first root and coefficient list for v,
source exterior bound when relevant, complete source/kernel/inverse premises,
and independent current same-face S8. A requiresS8 flag alone proves none
of this disjunction. The conclusion is an exact completion exclusion under
NO-S8; no intrinsic unit or invented actual paying receipt is added.

The separate source fragment is definitions only and requires peer/root
review before a coalesced successor. Frozen92b9 helper, F3b984 and all
receipts/failures remain unchanged.

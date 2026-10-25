# Finite57 unbounded tail: source-bound algorithm audit

Date: 24 October 2026. Author: consolidation agent; independent root review pending.

This report assesses the unbounded-tail argument in Appendix “The unbounded
tail: actual framing walls and two budgets” of
[finite_volume57.tex](../manuscripts/finite_volume57.tex:2522), with its
specified mathematical reports and certificate producers/reviewers. It
does not assess the remaining bounded finite-volume proof, replace the
root geometry audit, change a dimension claim, or certify a new execution.
No worker was imported, no enumeration or certificate arithmetic was run,
and no scientific, Lean or TeX process was launched for this assessment.
The canonical manuscript and all compact17 draft/frozen sources were
left unchanged.

The mathematical reduction is sound at the explicit computed boundaries
described below. In particular, the two original cost budgets are two,
and the smaller residual budgets are justified only after selecting an
exceptional failing factor. The local finite algorithms have a complete
support domain, and the all-rank tails are analytic arguments rather than
extrapolations from boundary ranks. The singular-active-system sentence
in the manuscript admits the complete justification in Section 4 below.

Two qualifications require precise treatment. The four-thirds A-cycle
paragraph needs its trivial-threshold restriction spelled out before
asserting nonnegative margins. Also, the latest same-row E7 verifier
checks four valid supplied certificates without reconstructing their four
required keys. The fixed JSON does contain the required four distinct
keys; this is a verifier coverage gap, not an observed missing certificate.
Some auxiliary receipt arrays use presence lookups rather than exact
key-set comparison. Section 7 records the distinction in detail.

## 1. Assessed, inherited and still assumed

| Item | Assessment in this report | Numerical or external boundary |
|---|---|---|
| Original budgets, scale cancellation, framing-wall applicability | Interface checked against manuscript and root geometry audit | Actual geometry and Nikulin input were separately audited by root |
| Complete finite support construction | Source proof checked | The saved support/KKT propositions have not been freshly executed here |
| Convex KKT certificate and singular-active skipping | Complete mathematical argument checked below | Rational supplied witnesses remain numerical inputs |
| Full-budget classical dispatch | Contact-table proof and all-rank A reduction checked | 580 small-A witness inequalities remain an explicit computed proposition |
| Exceptional first-failing classification | Complete enumeration algorithm and its use checked | Exceptional census and same-row inequalities remain computed propositions |
| Residual domination and ordered packing | Domination direction, recursion and independent reachability checked | 1395 dominations and three exact maxima remain computed propositions |
| Analytic tails at every rank | Literal coordinate/cycle arguments checked | Saved symbolic coefficient outputs were read, not recomputed |
| Four-thirds far range | Proof checked with the threshold clarification below | Finite exceptional quadratic outputs remain computed propositions |
| Combined finite57 theorem | Outside this report | No full-manuscript runtime PASS or end-to-end formalization is claimed |

The root interface report is
[audit_finite57_geometry_core_and_tail_interfaces_root_20261024.md](audit_finite57_geometry_core_and_tail_interfaces_root_20261024.md:166),
Sections 5–7. Its direct reading of Nikulin supplies the unrestricted
narrow-framing statement. This report does not turn that statement into
an arithmetic hypothesis and does not independently claim another
literature verification. The compact-facet alternative, bounded by 143
in this tail argument, likewise remains at that audited interface.
No assertion that an arbitrary facet is itself Coxeter is needed.

The main source proof is
[fresh_nikulin_two_budget_frontier_20261015.md](fresh_nikulin_two_budget_frontier_20261015.md).
Its independent review and the unit-charge, defect-frontier and
height-charge reviews were read together with the actual script functions.
Counts and inherited status strings are evidence of earlier records,
not new runtime results of this audit.

## 2. The dependency order and the two genuine budgets

Let the selected narrow facet contain a cusp. Its affine factors have
ranks summing to N=n−1. Normalize old and new roots to square two and
each positive affine kernel δ to minimum entry one. Choose both omissions
in each factor before choosing the outgoing edges. Only the factor
containing the selected facet normal has an avoid-node requirement.

The two terminating walls give nonnegative contact rows p and q, with
each nonzero entry at least one. At the actual other endpoints, the
positive-semidefinite principal Grams give separate Schur budgets

    sum_a C_0a ≤ 2,              sum_a C_1a ≤ 2.

No original endpoint budget is 3/2. The inverse of each retained
spherical Gram is entrywise nonnegative. Actual affine heights
σ_a=δ_a·p and τ_a=δ_a·q are positive. The factor kernel vectors are
positive multiples of the same cusp vector, so σ_a/τ_a=q_0 is common
across all factors. Put M=max(q_0,q_0^−1) and
W_a=min(σ_a,τ_a)^2 E_a. The Lorentz identity is

    −(x,y) = (1/2) sum_a σ_a τ_a E_a − q_0 − q_0^−1
           ≥ M((1/2) sum_a W_a − 2).

The targets are T_a=2r_a/9 ordinarily and 2(r_a−1)/9 in the sole
avoid-node factor. Their sum is (2N−2)/9. The proof has the following
acyclic order:

1. At the ORIGINAL budgets two, prove W+C_0+C_1≥T for every classical
   factor in the near range M≤21/20.
2. At the same original budgets, enumerate every exceptional support
   relevant to T−W>C_0+C_1. Remove the two same-row E7 support occurrences
   using the four required full-support certificates.
3. If a factor fails, select ANY one of the remaining exceptional
   occurrences. Its positive cost floors leave residual budgets at most
   3/2 on each side.
4. Apply the complete residual continuous envelope and ordered packing
   to all OTHER factors. Any further failures are included there.
5. If no factor fails, summing the two original budgets already bounds
   total loss by four.

An ordinary E8 failing factor has floors (8/9,2/3), target 16/9 and
leaves (10/9,4/3). The other possible failures, ordinary E7 or modified
E8, have floors (1/2,1) or (1,1/2), target 14/9 and leave
(3/2,1) or (1,3/2). The selected factor's own loss is at most its target
because W≥0. Residual losses 20/9 and 22/9 therefore give total loss
at most four in both cases. Allowing one modified residual factor even
when the selected factor is itself modified only enlarges the bound.

For the far range the separate local coefficient-four-thirds estimate
gives total loss at most (4/3)(2+2)=16/3. The two resulting lower bounds
on −(x,y) are M(N−37)/9 and M(N−43)/9. At N=163+t they are at least
14+t/9 and 14+7t/60 respectively. The framing bound is STRICTLY less
than 14, so the equality boundary t=0 is also excluded.

This tail excludes n≥164. It is only the all-rank closure of the finite57
proof; its numerical ceiling must not be confused with the bounded
dimension-58 contradiction elsewhere in the manuscript.

## 3. Finite support domains: why arbitrary magnitudes are covered

Write R(i,j)=||v_i−v_j||^2 for the displayed point potentials. For an
omitted node o, the inverse entry on retained nodes k,l is

    J^o_kl = (δ_k δ_l/2)(R(o,k)+R(o,l)−R(k,l)).

The independent reviewer reconstructs the physical norm-two Gram,
its positive kernel, every retained inverse and the resistance
identity. It checks positive definiteness of retained Grams and
entrywise nonnegativity of their inverses. Thus a support S has
minimum possible cost at unit nonzero contacts

    c_o(S) = sum_{k in S} δ_k^2 R(o,k)
           + sum_{k<l in S} δ_k δ_l
               (R(o,k)+R(o,l)−R(k,l)).

This is a LOWER bound for every row with the same support and arbitrary
contact magnitudes at least one. There is no finite angle alphabet in
this reduction.

For a cost cap c, every retained contacted node must satisfy
δ_k^2 R(o,k)≤c. Enumerate every subset of those eligible nodes, then
retain exactly the subsets whose complete same-row cost is at most c.
In particular, it is not enough to enumerate singleton or double
supports unless the complete subset calculation has first proved that
larger ones exceed the threshold.

For each pair of retained supports, independently include or omit the
omitted-node contact on each side. Included contacts also have magnitude
at least one. Require each full support U,V to be nonempty. Closed
convex hulls allow zero barycentric weights even when an actual support
has positive weights; this enlarges the domain and gives a valid lower
bound. If both retained supports are empty, the only possible full
rows are the two omitted-node contacts; the separately checked
point bound W≥T applies.

For B=21/20, full supports have affine heights at least
s=sum_U δ_k and t=sum_V δ_k. Their common ratio gives

    min(σ,τ)^2 ≥ H_B(s,t)
               = max(min(s,t)^2, max(s,t)^2/B^2).

This function is nondecreasing in each positive argument. Hence certified
rational LOWER approximations to the square roots defining δ give valid
lower bounds for H_B. The square-root routine certifies its interval by
squaring. In support rejection, a LOWER approximation to a positive
cross term can reject a support only when that lower cost already
exceeds the cap; no feasible support is lost. Every surviving producer
cost is then checked exactly. The reviewer independently computes the
physical exact cost and rejects any surviving nonrational cost rather
than treating a rounded quantity as exact.

The residual finite model list covers:

| Models | Domain |
|---|---|
| A | ranks 1 through 16 |
| Ordinary B,C | ranks 2 through 7 |
| D | ranks 4 through 17 |
| Modified B | ranks 3 through 7 |
| Modified C | ranks 2 through 7 |
| Exceptional E6,E7,E8,F4,G2 | all ordinary and prescribed avoid-node cases |
| Modified B2 | separate literal Gram, not an overlapping prefix-index convention |

There are 91 model descriptors: 16+12+14+5+6+37+1. Larger ranks are
covered by the analytic tails in Section 8. Rotations for A and the
end/fork automorphisms for D choose omissions avoiding a prescribed
node while retaining the stronger ordinary estimate. Modified A1 has
target zero, even though its allowed omission is repeated. It needs
no positive-loss record.

The algorithm terminates because there are finitely many descriptors,
finitely many eligible nodes and their subsets, four omitted-node
states per support pair, and finitely many rational active sets per
convex-distance calculation. A failed convex certificate raises an
error; it is never silently discarded as an impossible geometric row.

The exact finite output proposition used in the manuscript is 1116
retained support pairs, 3630 nonempty continuous domains, 1395 positive
losses and a 28-record Pareto frontier. These numbers have not been
recomputed here. Exhaustive domain construction, rather than the
numbers themselves, is the coverage argument.

## 4. KKT certificate, including singular active systems

For full supports U,V, let K be the signed point Gram and A the two
row-sum matrix. The signs are positive for U and negative for V.
Because K is a Gram of actual Euclidean point vectors, K is positive
semidefinite even if points coincide or the Gram is singular. A
supplied certificate consists of z≥0, Az=(1,1), multipliers λ, and

    s = Kz+A^t λ ≥ 0,       z_i s_i=0,
    e = z^t Kz = −λ_1−λ_2.

For every other feasible z′,

    (z′)^t Kz′ − e
       = (z′−z)^t K(z′−z) + 2(z′)^t s ≥ 0.

This follows by expanding, using A(z′−z)=0 and z^t s=0. Thus the
certificate proves the global convex-hull distance minimum. No
invertibility hypothesis occurs in verification.

The producer searches all active subsets containing at least one
index on each side and solves the bordered KKT system. A singular
system is skipped. Here is a complete reason this cannot lose every
minimizer.

Choose a global minimizer with the fewest positive coordinates; one
exists because the product of two finite simplices is compact. On its
positive support I, both row constraints are nonempty and independent.
First-order variation between two positive coordinates on the same
side shows that Kz is constant there, while a variation into a zero
coordinate gives the required nonnegative slack. Thus a KKT multiplier
exists.

Suppose the bordered matrix on I were singular. A null pair (h,ν)
satisfies A_I h=0 and K_I h+A_I^t ν=0. Multiplying by h gives
h^t K_I h=0. Positive semidefiniteness implies K_I h=0, and independence
of the two rows of A_I then gives ν=0. Consequently h is nonzero.
Since A_I h=0, its entries have both signs within a nonzero row block.
Move the minimizer along h until a positive coordinate first reaches
zero. The constraints and nonnegativity remain valid; the objective
is unchanged because K_I h=0. This produces a minimizer with smaller
positive support, a contradiction.

Therefore some minimizing active set has a nonsingular bordered
matrix. Its rational system has a rational solution. The finite
active-subset search will find a valid witness; singular attempts
need not be retained. Conversely, if all attempts fail, the producer
raises its explicit “no convex certificate” assertion. A support is
never removed on the ground that one active system was singular.

The reviewer need not trust this optimization search: it verifies the
simple rational identities of the supplied witness against its own
physical Gram. Search completeness is needed to explain producer
termination with output; certificate soundness follows from the
displayed expansion alone.

## 5. The full-budget coefficient-one bridge

The residual 3/2 census cannot justify selection of the first failing
factor. This bridge is supplied separately by the manuscript's
full-budget classical lemma and its scripts.

For omissions i,j, set

    c_ik = δ_k^2 R(i,k),
    ell_ik = R(i,j)+R(i,k)−R(k,j).

For normalized full-row means v_i+a and v_j+b, polarization gives

    E = R(i,j) − sum_k w_k ell_ik − sum_l z_l ell_jl
                    + ||a−b||^2.

Every height is at least one, so W≥E. Unit diagonal costs bound the
actual retained costs below, including all nonnegative same-row
cross terms. Convex weights are at most one. For every B,C,D contact
in the complete 24-class table, either ell_ik≤c_ik or c_ik>T.
The latter branch is already trivial because C_0+C_1>T.
The omission distance R(i,j) is at least T. Thus the former branch
gives W+C_0+C_1≥T. The separate modified B2 Gram has retained diagonal
at least 1/2>T=2/9 and omitted-only distance one. The whole B,C,D
coefficient-one argument is independent of M.

For A, the ORIGINAL per-side cap is two. At each 2≤m≤17 the scripts
reconstruct all retained subsets of cost at most two, keep pairs whose
total unit cost is below T, and form all nonempty full supports after
both omitted-node choices. A rational KKT witness plus the height
factor H_21/20 proves the coefficient-one inequality. The 580-case
proposition is explicit; it is not obtained by applying the residual
3/2 census.

The older height-charge checks used B=31/30. The defect-frontier
review reconstructs the same support keys and recomputes the height
factor at B=21/20 from the saved exact distance witnesses. Merely
citing the older 31/30 scalar inequalities would not justify the
larger near range.

For m≥18 even or m≥19 odd, cap two permits distances at most two
from the omission. Put ε=m mod 2 and h=(m−ε)/2. The linear-loss to
unit-cost ratio is at most L=(m+ε)/(m−2). Hence

    W+C_0+C_1−T ≥ R(h)−T−4(L−1).

The even numerator is m^2+6m−304 over 36(m−2), and after m=18+2u
it is 128+84u+4u^2. The odd numerator is
m^3+6m^2−457m+18 over 36m(m−2), and after m=19+2u it is
360+1708u+252u^2+8u^3. Denominators are positive and all coefficients
are nonnegative. Thus this is an all-rank coefficient-one theorem,
not a finite-rank observation. Cycle rotation gives the avoid-node
version at the weaker modified target; repeated avoid-A1 has target
zero.

## 6. Exceptional selection, domination and exact ordered packing

The exceptional ORIGINAL-budget census uses the literal A1,E6,E7,E8,
F4,G2 matrices and all ordinary/avoid-node cases: 40 descriptors,
including the three A1 descriptors. These are not the 37 exceptional
descriptors in the residual model list.

For each descriptor the reviewer computes every retained inverse,
collects all contacts with diagonal below T, and enumerates ALL subsets
of those contacts across the two ordered sides. It includes every
same-side cross term and keeps only total unit costs below T. A
failure must be in this list because actual magnitudes are at least
one and inverse entries are nonnegative.

The fixed output has 71 singleton occurrences and 16 multiple-support
occurrences. Each singleton is treated by its exact one-variable
convex quadratic, with separate omitted-contact states and the widened
height factor. The special two-sided E7 support has energy
(1+u^2+v^2)/2 and is already closed. The remaining 15 candidate failure
occurrences are ordinary E7/E8 or modified E8. Two same-row ordinary
E7 occurrences at nodes 4,7 are closed by the four exact full-support
certificates in Section 7. The remaining 13 have precisely the
ordered floors and targets stated in Section 2. Counts and
inequalities are inherited computed propositions here; the complete
subset algorithm, not a count of 15, proves the reduction's domain.

For the residual product, each computed positive loss has a dominating
record of the same modified flag, no larger LEFT and RIGHT cost floors
and no smaller loss bound. The reviewer independently derives the exact
Pareto frontier from all 1395 positive cases, checks the 28 records are
distinct and equal that frontier, and checks the supplied dominator
for every positive case. It also checks that every record has total
cost at least 1/2.

Dropping nonpositive-loss factors can only enlarge the sum to be
bounded. Replacing the others by their dominating records uses no
geometric realization assumption. A finite multiset of records can
be sorted by record index without interchanging the two endpoints.
The producer recursion, which permits repeating its current index,
therefore enumerates every such multiset, including repeated factors.

Its recursion terminates: each added record spends at least 1/2 of
combined budget, so depth is at most floor(2(P+Q)). All costs are
rational with a common denominator, leaving only finitely many
budget states. Induction on remaining possible depth proves that
the maximum of zero and every admissible next-record continuation
is the exact packing maximum.

The reviewer uses a different forward reachability formulation on
(used-left-cost, used-right-cost, modified-count). Transitions
strictly increase total used cost, and the state space is finite.
Keeping the largest achieved loss at a state is valid because all
future admissibility depends only on that state. Processing every
improvement reaches the exact same multiset maximum. A supplied
equality witness is checked as well as the maximum.

The three required exact maxima are

    budgets (10/9,4/3): 20/9,
    budgets (1,3/2):   22/9,
    budgets (3/2,1):   22/9,

with at most one modified record. The last two ordered budgets are
separate entries. Sorting does not authorize silently swapping a
record's left and right floors. These maxima were not recomputed
during this source audit.

## 7. Exact key coverage: what existing verifiers do and do not reject

| Table or support domain | Existing domain reconstruction | Exact rejection of missing/duplicate/out-of-domain supplied keys |
|---|---|---|
| 91 residual model descriptors | Independent literal descriptor list, positional metadata, omissions, kernel squares and targets | Yes for the model list |
| Retained domains in each residual model | All eligible-node subsets; exact physical inverse costs; saved list equality | Yes, including exact saved order |
| 3630 residual full-support cases | Rebuild both omitted-node states; construct expected (U,V) map; consume by pop; require empty remainder | Yes |
| 1395 positive losses and 28 frontier records | Derive positive list and independent exact Pareto set; check every dominator | Yes for frontier keys and all positive source cases |
| Small-A original 580 distance witnesses | Rebuild supports at cap two for each m; consume (U,V) keys; compare count and require empty remainder | Yes INSIDE each selected old cycle row |
| Sixteen old small-cycle wrapper rows | Select required m by lookup in 2,…,17 | Presence checked; extra or duplicate wrapper rows are not uniformly rejected |
| Widened 580 scalar cycle rows at B=21/20 | Each independently reconstructed witness key is looked up and its scalar compared | Required-key presence/value checked; extra or duplicate scalar rows are not uniformly rejected |
| Full-budget exceptional cheap supports | All literal descriptors and all cheap-contact subsets reconstructed | Mathematical support domain complete; some saved receipt arrays use lookups and length checks rather than a uniform exact key-set guard |
| Four same-row E7 full-support certificates | Supplied supports individually checked with physical KKT/height and length four | NO required four-key reconstruction in the latest reviewer |
| Thirteen remaining failure floors | Rebuilt from independently reconstructed earlier 15 occurrences; remove precisely two same-row occurrences; compare complete ordered output | Yes for this derived list, conditional on the earlier census |
| Three packing receipt rows | All three maxima independently computed; corresponding saved row looked up | Required presence/value checked; extra/duplicate receipt rows not uniformly rejected |
| Five B/C/D tail receipt rows | Five analytic classes independently reconstructed; matching receipt looked up | Required presence/value checked; extra/duplicate receipt rows not uniformly rejected |
| Original 24 classical coordinate receipt rows | Literal geometric contact partition proved in manuscript/root review; script validates supplied 24 entries | The length check alone is not a universal exact 24-key rejection guard |

The strongest full-domain loop is
[two_budget_frontier_review.py:finite_census](../scripts/fresh_nikulin_two_budget_frontier_20261015_review.py:136):
its expected dictionary contains every admissible full-support tuple;
each saved case removes that tuple exactly once and the dictionary must
end empty. Its independent physical-Gram reconstruction does not trust
producer-side support rejection.

The small-A keyed witness loops are
[height_charge_review.py:cycles](../scripts/fresh_nikulin_height_charge_20261015_review.py:169)
and
[defect_frontier_review.py:cycle_replay](../scripts/fresh_nikulin_defect_frontier_20261015_review.py:140).
They reconstruct the actual supports; the distinction between their
strict inner witness maps and their lookup-selected outer receipt rows
is material. A statement that EVERY array in the whole chain rejects
duplicate keys would overstate the present scripts.

### The four same-row E7 keys

The required omission tuple is (3,6). Support tuples use the producer's
order: retained nodes first, then the optional omitted node.

| side | own_omitted | left tuple U | right tuple V |
|---:|:---:|---|---|
| 0 | false | (4,7) | (6) |
| 0 | true | (4,7,3) | (6) |
| 1 | false | (3) | (4,7) |
| 1 | true | (3) | (4,7,6) |

The producer explicitly generates exactly these keys. Direct inspection
of the fixed primary JSON confirms that it contains all four, with these
supports. This observation checks the actual domain entries, not their
numerical energies or margins afresh.

The independent review's
[bad_factors](../scripts/fresh_nikulin_two_budget_frontier_20261015_review.py:277)
instead iterates over the four supplied rows, verifies each KKT witness
and positive margin, and checks length four. It does not bind
side/own_omitted to the required U,V or compare the complete key set.
For example, four copies of one valid row could satisfy that particular
local loop. This does not describe the actual fixed JSON, whose keys
are correct. A small separate successor guard should enforce exactly
the four rows above before the inherited arithmetic verification is
credited as a generic complete certificate verifier.

For other arrays using a lookup, missing REQUIRED keys still fail,
and required values are checked in a freshly reconstructed loop.
Permitting redundant extra receipt rows usually does not invalidate
the mathematical inequality. Nevertheless, an advertised strict
certificate contract should reject them explicitly. The same-row
E7 case is more substantial because its required domain keys are
not reconstructed by that local loop at all.

A supplement should preserve every pinned legacy producer/reviewer
byte, compare key sets and cardinalities before invoking the existing
numerical logic, and clearly identify which expected keys it constructs
independently and which it cross-binds to an already strictly verified
inner witness table. It must not call a cross-binding operation a fresh
geometric enumeration.

## 8. Analytic tails under the RESIDUAL cap 3/2

The larger ranks are not missing from the 91-model finite census.

### B,C,D coordinate tails

For r≥8 in B,C and r≥18 in D, the complete eligible retained-node
lists follow directly from the literal prefix potentials:

| Type | Omissions | Eligible retained nodes on the two sides | Lower E | Lower min-height squared |
|---|---|---|---:|---:|
| Ordinary B | (0,r) | {1}; {r−1} | (r−2)/4 | 2/B^2 |
| Modified B | (r−1,0) | {r−2,r}; {1} | (r−3)/4 | 2/B^2 |
| Ordinary C | (0,r) | {1}; {r−1} | (r−2)/2 | 1 |
| Modified C | (1,r) | {0,2}; {r−1} | (r−3)/2 | 1 |
| D | (0,r) | {1}; {r−1} | (r−2)/4 | 1 |

The listed sets are eligible NODE lists; a row may use any subset
consistent with its full inverse cost and may additionally contact
its omitted node. Enlarging to the convex hull of all eligible nodes
and the omission is harmless for the lower energy bound.

Every other node is excluded by a literal inverse-diagonal class:
r/4, (r−1)/4, r/2, (r−1)/2 or an integer at least two, as applicable.
At the stated starting ranks these exceed 3/2. The full convex hulls
differ on fixed middle coordinates: magnitude 1/2 in B,D and one in C,
with the half metric in C. Counting those coordinates gives the
displayed E. For B, all possible full nodes on one side have δ^2≥2;
the near-range ratio therefore gives min-height squared at least
2/B^2. Other types need only height at least one.

Subtract T=2r/9, or 2(r−1)/9 when modified. In each row the resulting
margin is affine in the excess rank with a nonnegative constant and
a positive slope. The D18 starting constant is zero, which suffices
because the required inequality is nonstrict.

Both scripts also inspect three boundary ranks per class. Those
checks illustrate the literal formulas. They are not the argument
at arbitrary rank: that argument is the exhaustive coordinate-class
partition, fixed middle coordinates and symbolic margin.

### A-cycle residual tail

For even m≥18 or odd m≥19, a distance-two retained contact costs
2(m−2)/m>3/2. By cycle concavity, every more distant contact also
exceeds the cap. Two adjacent retained contacts already have diagonal
cost sum at least 2(m−1)/m>3/2, even before nonnegative cross terms.
Thus each side has at most one adjacent retained node and its optional
omitted-node contact.

Put h=floor(m/2), ε=m mod 2. If both full supports are singletons,
their cycle distance is at least h−2. The margin W−T is bounded below by
R(h−2)−T, with numerator

    m^2+8m−9(4+ε)^2

over 36m. If either full support has two nodes, that row's height is at
least two, and the near ratio gives min-height squared at least 4/B^2.
Moving either normalized mean by at most one adjacent edge loses at
most (m+ε)/m in the polarization's linear term. Hence

    E ≥ R(h)−2(m+ε)/m.

The resulting W−T numerator is

    (9−2B^2)m^2+(2B^2−72)m−81ε

over 9B^2m. After substituting m=18+2u or 19+2u, both numerators have
nonnegative coefficients and positive constant term. Denominators
are positive. Thus every residual A factor outside the finite list
has no positive loss. These are two separate cases; the proof does
not use an old coefficient-one margin as a substitute for W≥T.

## 9. The four-thirds A paragraph needs its exact threshold

The desired all-ratio inequality is

    W+(4/3)(C_0+C_1) ≥ T,       T=2(m−1)/9.

Write C=C_0+C_1. If C≥3T/4=(m−1)/6, W≥0 proves it. Only the STRICT
nontrivial branch C<(m−1)/6 needs the contact/polarization argument.

For 3≤m≤6, even the smallest retained unit diagonal (m−1)/m is at
least (m−1)/6. Consequently there is no retained contact in the
nontrivial branch. The m=2 ordinary case is covered by the separate
literal A1 argument; its modified target is zero.

For 7≤m≤10, distance two costs 2(m−2)/m≥(m−1)/6. Thus a nontrivial
row can use only distance one. A distance-two contact can first occur
in the nontrivial branch at m≥11. For m≥10 the original per-side
budget two excludes every distance at least three, since
3(m−3)/m>2 and the cyclic resistance increases up to the antipode.

For cyclic distance d=1,2 the linear loss is at most d(m+ε)/m,
while its unit diagonal cost is d(m−d)/m. Therefore

    (4/3)c−ell ≥ d(m−4d−3ε)/(3m).

This margin is nonnegative on the preceding NONTRIVIAL domain:
distance one first needs m=7 in odd parity, with zero margin there;
distance two first needs m=11 in odd parity, also with zero margin.
Smaller odd cases where the displayed expression would be negative
were already removed by C≥3T/4. Even m=10 permits a nonnegative
distance-two margin analytically, although such a contact is already
trivial at that order. Enlarging that even analytic subdomain is safe.

The unit-charge reviewer implements this split: m=3,…,6 has no
nontrivial retained contact; m=7,8,9 uses distance one; even m=10+2u
and odd m=11+2u use the nonnegative d=1,2 formulas. The manuscript
currently compresses all of this to “below the trivial threshold”
and “remaining margin … ≥0.” The threshold value and domain exclusion
should be inserted explicitly, so readers do not infer the false
unrestricted claim for small odd m.

For B,C,D the coefficient-one contact argument already proves a
stronger estimate independent of M. The exceptional four-thirds
census uses the cheaper threshold 3T/4, all cheap subsets and both
omitted-node states. Its fixed output has 40 descriptor cases,
46 single-contact occurrences and one nontrivial two-contact support.
For one contact, endpoint and stationary minima of the exact convex
quadratic are checked. The sole multiple support is the two-sided E7
case with energy (1+u^2+v^2)/2. No common-height near-range hypothesis
is used here. Those finite inequalities remain explicit computed
inputs in this source assessment.

## 10. Exact computed propositions left at the boundary

The following propositions were READ as inherited output statements,
not newly run or arithmetically revalidated in this audit.

1. Original small-A supports for m=2,…,17: 580 exact distance witnesses,
   reconstructed at per-side cap two, satisfy their KKT conditions and
   the coefficient-one inequality with B=21/20.
2. Original exceptional census: 40 descriptors, 71 singleton and 16
   multiple-support occurrences; after singleton and special-double
   inequalities there are 15 candidate failure occurrences.
3. The four required same-row E7 support certificates prove
   H_B e+3/2>14/9; removing their two support occurrences leaves exactly
   the 13 exceptional ordered floor records.
4. Residual finite census: all 91 prescribed descriptors, 1116 support
   pairs and 3630 full continuous domains have the recorded exact
   costs, rational KKT witnesses and certified height lower bounds.
5. The resulting 1395 positive losses are dominated by the complete
   28-record ordered Pareto table.
6. Ordered residual packing has the three exact maxima in Section 6.
   Producer memoized-state counts are checksums, not coverage proofs.
7. Finite exceptional four-thirds outputs and the saved exact
   coefficient signs in the literal all-rank formulas are as stated.

Items 1–6 are distinct domains and must not be conflated. In particular,
580 is not a residual support count; 3630 is not an original-budget
exceptional census; 40 includes A1 descriptors that the residual
37 exceptional list does not include. The all-rank polynomial formulas
are supplied above so their mathematical role is visible without
treating a status string as a proof.

Any conditional Lean kernel should name these interfaces, the actual
geometry, narrow-framing statement and compact-facet alternative as
explicit parameters. This report supplies no new axiom audit.
Printing Lean axioms alone would not enumerate theorem parameters or
discharge an enumeration premise.

## 11. Necessary followups, in order

1. Root independently reviews this report, especially the first-failure
   split, the minimal-support proof for singular KKT systems and the
   small-A threshold paragraph. Until that review, this is a proposed
   source audit, not a new accepted global status.
2. Prepare a SMALL separate source-only key-coverage supplement.
   Require the four literal E7 keys, exact 16 small-cycle wrapper keys,
   exact widened scalar keys cross-bound to the strictly verified old
   580 witnesses, all 40 literal exceptional descriptor keys, three
   ordered packing receipt keys and five B/C/D tail keys. Any additional
   inner support checks must say explicitly which source table supplies
   the expected set. Preserve every pinned producer/reviewer byte.
3. Have another author review that supplement before root schedules its
   minimal bounded execution. It should parse/cross-bind existing
   certificates, not regenerate optimizers or an all-rank enumeration.
   No launch is authorized by this source report; admission and the
   exclusive runner remain root responsibilities.
4. In a reviewed same-family manuscript revision, add the short
   singular-active justification and replace the compressed four-thirds
   A paragraph with the threshold/domain argument in Section 9.
   The canonical TeX/PDF must remain unchanged before that revision is
   reviewed and intentionally adopted.
5. Keep this tail's external statements, internal logic and computed
   propositions separate in the reproducibility and Lean logic maps.
   Existing successful modules should not be rerun merely to update
   prose. No post57-to56 research lemma belongs in the finite57 paper
   simply because that lemma is accepted elsewhere.

The compact17 frozen batch remains at
[source_archive/compact17_foundation_transcription_20261024](source_archive/compact17_foundation_transcription_20261024),
with its previously delivered ledger and freeze manifest. Its incomplete
working appendices were not integrated or modified during this audit.

## 12. Source binding for this assessment

SHA-256 values below bind the bytes inspected. They are provenance
metadata, not runtime receipts.

| Path | SHA-256 |
|---|---|
| manuscripts/finite_volume57.tex | c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0 |
| results/fresh_nikulin_two_budget_frontier_20261015.md | e2a54cd23b9467b5ee6bd3ce04f2db4f1cc2f2a82521102f42c44568c67c07fa |
| results/fresh_nikulin_two_budget_frontier_20261015_review.md | 27ad96416ddd7e531845560c0d43e7b7850c6316706264af31beb008daa034ec |
| results/fresh_nikulin_unit_charge_20261015_review.md | 999a3694e1889ecdb7932c97125e7f4be8b2b55b45a14a215a9b88b188e8d4fc |
| results/fresh_nikulin_defect_frontier_20261015_review.md | b644fdc5d9ed26a8c6eb3917dc211957e980f56d09fc003022b5cd3c8d7116c8 |
| results/fresh_nikulin_height_charge_20261015_review.md | 05d9dd9e050ac0bb6f5896bac1c47c00d796cd35ea3b82ff16ca9e0c9ba590a3 |
| scripts/fresh_nikulin_two_budget_frontier_20261015.py | 87187a3fb1202bc566de688defdcc57d1d0aa231b12f6c180ab5aaff0ffd9075 |
| scripts/fresh_nikulin_two_budget_frontier_20261015_review.py | 355eff73dbce58b079471efe09f0930544d1807c75621557be71b609384e8c8a |
| scripts/fresh_nikulin_unit_charge_20261015_review.py | 362dd81bb17dc5505485c1f71c8711c35ccde3472bc8b2cfbd68bfe6b4d8a930 |
| scripts/fresh_nikulin_defect_frontier_20261015_review.py | da34db7aef35f73a7bb516078aea5387fc1a1fab3435836bdb79109744e0525c |
| scripts/fresh_nikulin_height_charge_20261015_review.py | 6afcde9bd55bb88e62221ce0fa949a300880fde15dc9fe08359de8260b5b1f74 |
| results/fresh_nikulin_two_budget_frontier_20261015.json | 1a32a2aa37b44b73ab03a322c136e846cea5c37028718e5241cfdad779edfdc9 |
| results/fresh_nikulin_two_budget_frontier_20261015_review.json | 4df5e7a7bd604e9a279a075aa02886cc05402930c4627b49e3dc21d93e9413d0 |
| results/fresh_nikulin_height_charge_20261015.json | 6427dadcaa1577d0416b31f48da7b67425eb647d4661ab464c2e7f55f7098b6b |
| results/fresh_nikulin_height_charge_20261015_review.json | bebad068dac2e18fc70830adbb3fc9613487414ad2f1b3e0061dd399e9bfd449 |
| results/fresh_nikulin_defect_frontier_20261015.json | 3c8897e52d781b97c9b5de222e276234c84b8a95770f7039d2d2b99880614d04 |
| results/fresh_nikulin_defect_frontier_20261015_review.json | a45b3764e1ba9fdc71e3c28a03d36368e1b4df0ded8f26f4dd585a26e9c63da2 |
| results/fresh_nikulin_unit_charge_20261015.json | c6b85ef3c5abbcd53c9f90347f7e129ecd52285fac4672cea0b67563a67bd888 |
| results/fresh_nikulin_unit_charge_20261015_review.json | 7b1d281d52b2d2e171b2ffe1b4155ca55708e333889c4b2525af579de5f88c14 |


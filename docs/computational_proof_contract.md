# How a computation enters a mathematical proof

20 October 2026. This document gives a common mathematical interface for
the current papers. It is not a claim that every existing checker already
establishes every obligation below. An execution report, a proof of an
algorithm and a geometric coverage theorem establish different things.

## The statement before the script

For each proof-relevant computation, the paper must specify:

1. **Objects and bounds.** Define the admissible objects, parameters,
   equivalence relation and exact arithmetic domain. Prove that every
   geometric object needed in the argument maps into this domain. If a
   parameter is initially unbounded, prove its reduction to finitely many
   cases or supply a uniform symbolic argument.
2. **Algorithm.** Give mathematical pseudocode, including initial states,
   transitions, stopping tests, symmetry reductions and output. A filename
   is an implementation reference, not the definition of an algorithm.
3. **Termination.** Give a decreasing integer, a finite state space or an
   explicitly bounded loop. A root-count cap chosen for convenience is
   not a proof of termination with complete output.
4. **Soundness and completeness.** Prove that each accepted object has the
   stated property and every admissible case reaches a retained branch
   or a branch whose conclusion has already been proved uniformly.
   Every pruning rule and symmetry quotient needs this justification.
5. **Computed proposition.** State the actual exact output needed by the
   theorem: a maximum, a list of exceptional orbits, or a collection of
   rational inequalities. Name the checked input, output certificate,
   checker and execution record in the reproducibility supplement.

The main proof cites the computed proposition. The supplement explains
how to regenerate and independently check it, with explicit software
dependencies and resource limits. It distinguishes regeneration from
checking a saved certificate and from merely checking its hash.

## A basic forest algorithm, with proof

Let A be a finite alphabet of irreducible component types. Each type T
has a positive integral rank r(T) and exact rational cost w(T). For a
fixed nonnegative integer N, let F(s) be the multisets of A of total
rank s, for 0 <= s <= N. Repetition is allowed. Empty maxima below are
minus infinity; in an implementation they are an explicit unreachable
state, not a large negative numerical sentinel.

**Algorithm.** Set U(0)=0. In increasing order of s from 1 to N set

    U(s) = max { U(s-r(T)) + w(T) : T in A, r(T)<=s,
                s-r(T) is reachable }.

Store one predecessor for each finite maximum if an attaining multiset
is required.

**Proposition.** The algorithm halts and U(s) is exactly the largest cost
of a multiset in F(s), or is unreachable if F(s) is empty.

**Proof.** There are N stages and finitely many transitions per stage.
Each transition uses a smaller rank, since r(T)>0. Appending T to a
multiset of rank s-r(T) gives a multiset of rank s and the displayed
cost, so every candidate is attainable. Conversely, remove any one
component T from a nonempty multiset of rank s. Its remaining cost is
at most U(s-r(T)) by induction on s. Thus every admissible multiset is
bounded by a candidate considered by the algorithm. The same induction
proves reachability and validates predecessor reconstruction. ∎

The compact arithmetic paper uses this last-component argument after
fixing the hub counts. Its rank limit is 21. The alphabet, excluded hub
types and priority adjustments are part of the statement; they cannot
be inferred from a sample of maximizing forests.

## Checking an upper certificate without trusting a maximizer

For the same additive problem, a table B on all ranks can certify an
upper bound without recomputing a maximizing search. On reachable ranks
require B(0)>=0 and

    B(s) >= B(s-r(T)) + w(T)

for every admissible final component T whose predecessor is reachable.
Reachability itself must be reconstructed by the exact Boolean recurrence
with the same alphabet. Removing the last component proves by induction
that B bounds every multiset. A valid predecessor witness attaining B(s)
proves equality, if equality is claimed. The witness must be an actual
multiset of that exact cost; a tight predecessor chain suffices when it
starts at B(0)=0. An upper bound alone needs no attainment witness.

For context-dependent costs the state must retain all information needed
to compute transitions: for example fixed hub multiplicities, presence
guards and reserved exceptional components. A cost involving two component
types is not automatically additive. The paper must either prove that
the retained state determines every increment or fix those types before
applying the additive recurrence. A bonus for the presence of a component
is paid once, not once per repeated component.

A finite directed state graph has the same certificate proof when each
transition strictly increases used rank. The certificate must dominate
**every** allowed transition, including transitions outside the maximizing
path. Initial costs and terminal adjustments must be included explicitly.
Proving that a geometric forest reaches an allowed terminal state is a
separate coverage obligation. A checker that only replays one winner does
not establish a universal maximum.

## The terminal cusp enumeration used for dimension 57

Here is a concrete instance beyond the generic forest recurrence. The
terminal domain consists of partitions of 57 with every part at least 8,
except for at most one part in {1,...,7}. This is a domain of affine
**ranks**; the geometric reduction and the envelope bound are separate
lemmas of the finite-volume paper.

Define H(n,l) recursively to be the following list of nondecreasing
sequences, for integers n>=0 and l>=8:

    H(0,l) = [the empty sequence]
    H(n,l) = concatenate over r=l,...,n of
             [ (r) followed by p : p in H(n-r,r) ],  if n>0.

For each b in {0,...,7}, take all sequences in H(57-b,8), prefixing b
when b>0. To evaluate a resulting profile (r_1,...,r_k), start h=[1]
and multiply successively by 1+t+...+t^(r_i). Coefficients are computed
by the exact convolution

    new_h(j) = sum { h(j-a) : 0<=a<=r_i },

with out-of-range coefficients zero. Compute V=product(r_i+1),
G=sum_j h(j)*min(j,57-j), and the rational envelope margin using
the stated table c(r).

**Proposition.** This algorithm halts and lists each terminal profile
exactly once. It computes its h-polynomial and rational margin exactly.

**Proof.** Every recursive call with n>0 lowers n by at least 8. Each
call has finitely many children. The first part r of a nonempty sorted
partition is unique, and deleting it gives a unique element of
H(n-r,r); induction proves that H has neither omissions nor duplicates.
The unique part below 8, if present, determines b, so the outer loop
also has neither omissions nor duplicates. The convolution is the
coefficient formula for polynomial multiplication, proving the
evaluation claim by induction on the number of factors. All sums,
products and divisions in the margin are exact integers or rationals. ∎

The existing endpoint implementation reconstructs this domain, compares
it with the saved rows and recomputes every polynomial. Its recorded
output is **2552 profiles**, minimum normalized margin **87/10240**,
and minimum absolute surplus **36**. Separately, the concentration
argument reaches this domain by transfers (a,b)->(a-1,b+1) with
1<=a<=b<=7. Each transfer increases the squared-rank sum by
2(b-a+1)>0, while that sum is at most 57^2. This proves termination
of the reduction, rather than just termination of the finite checker.
The stop-loss inequalities and positive-margin argument prove that
the terminal bound also bounds the original cusp. The enumeration
alone does not establish those geometric or comparison claims.

## Direct finite enumerations and geometric queues

For a bounded integer-vector domain, nested loops are a valid complete
algorithm once all coordinate bounds are proved. Evaluate the defining
predicate exactly on every vector. Quotienting by symmetry is optional:
if used, prove that the group preserves the predicate and that each orbit
has a representative in the selected domain. An output count alone proves
neither claim.

For a queue of partial root/contact configurations, the mathematical
description must additionally prove that the initial states cover every
admissible object and that every valid extension is produced or discharged
by a uniform proof of the desired conclusion. Explain why each pruned
branch is either impossible or already satisfies the required payment or
other property for all its extensions. Give a proved bound on depth or a
genuinely decreasing measure.
Bounded depth with an unbounded label alphabet is still insufficient;
each extension domain must be finite or reduced by a uniform argument.

In particular, every singular stopping state needs an interpretation.
One must show that an actual geometric configuration either has already
been covered or follows a retained branch. Silently dropping singular
Gram matrices can omit precisely the boundary cases at issue.

## Discovery, exact acceptance and completeness

A numerical optimizer may propose a rational angle rule. Exact rational
evaluation of every necessary inequality can certify that proposal, so
the numerical optimizer need not be trusted. This works only after the
necessary domain has been proved complete.

If floating-point tests prune candidates, exact checks of the survivors
exclude false positives but do not exclude false negatives. Such a filter
requires a rigorous error/separation bound, or it must be replaced in the
proof computation by an exhaustive exact filter. Random checks and extra
parameter values provide diagnostics, not completeness.

## Current computed propositions and their limits

These are statements from the existing proof sources and execution
records, not fresh runs performed for this document.

| Paper | Exact finite statement used | What the endpoint replay does not establish |
| --- | --- | --- |
| [Finite-volume 57](../manuscripts/finite_volume57.tex) | Finite charge <=377; 2552 terminal ideal profiles have normalized margin at least 87/10240 and absolute surplus at least 36; the final strict correction is 4205. All 28 stop-loss transfers are checked. | The endpoint replay binds inherited tables; it does not re-enumerate the earlier finite Bellman domain or prove the geometric routing, factor-fiber or index interpretation. |
| [Compact arithmetic 20](../manuscripts/compact_arithmetic_dimension_bound.tex) | The rational vector pays at least 60007431727/60000000000 on 1979 native faces and 102917 common spectator cases. The complete rank-21 forest maximum is 40 across 68 types and 729 exact hub tuples. | Exact evaluation does not prove that every actual unpaid quadrilateral reaches the native iterator, or that the geometric incidence transfers and field reductions hold. |
| [Compact 23](../manuscripts/compact_farthest_vertex_bounds.tex) | The new proof is written geometric reasoning and displayed exact inequalities, not a new computer enumeration. | Its source-integrity replay proves byte identity only. A separately compiled conditional Lean kernel verifies the final accounting logic; the actual component energies, cap bounds and geometric obstruction remain inputs. |

For the first row, see the [source](../results/fresh_d13_context_dimension57_20261018.md)
and [registered commands](../manuscripts/catalogue.json), check IDs `finite57`,
`finite57_frontier` and `finite57_low8`. For the second, see the
[independent exact-check record](../results/fresh_compact21_self_hub_exact40_independent_check_20261018.json)
and the [compact-20 proof map](compact_arithmetic20_proof_map.md).
The [compact-23 proof map](compact_farthest_vertex_bounds_proof_map_20261020.md)
states its written-proof boundary.

The right-angled volume scripts check finite models and identities. They
do not enumerate all hyperbolic polytopes or replace the universal cutting,
recognition and equality arguments. Their computed propositions should
therefore say exactly which finite model or identity was checked.

## The conditional Lean interface

Keep three sorts of premises visible in the theorem signature and the
[assumption register](manuscript_assumptions.md): published external facts,
precisely stated script outcomes, and remaining internal geometric lemmas.
The user permits script/enumerator correctness to remain an explicit
assumption in Lean. The written paper must still explain the algorithm's
mathematical correctness and its connection to the geometric problem.

Lean should check how these facts imply the intermediate claims and final
contradiction. A kernel theorem conditional on a script result is useful;
it is not a kernel verification of that script. A theorem assuming an
unproved internal geometric lemma must retain that qualification. Standard
axiom lists alone do not expose assumptions passed as theorem parameters.

## Reproduction is serial and bounded

The [verification guide](manuscript_verification.md) supplies the existing
single entry point. Run one check or Lean module at a time under the shared
lock, record its declared mathematical scope, and remove temporary build
products after the batch. Start only with at least 40 GiB free disk and
6 GiB available memory; the abort thresholds are 8 GiB and 3 GiB. Below
20 GiB free disk, the coordinator must start a cleanup agent before further
proof jobs. Editorial metadata checks do
not relax these gates or count as scientific replay.

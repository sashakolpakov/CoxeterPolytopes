# Independent terminal finite-attachment review

**PASS: independent generation and fresh saved replay.** The
[independent checker](../scripts/fresh_sqrt5_nonsingular_rank29_fast_20261015_review.py)
and its [certificate](fresh_sqrt5_nonsingular_rank29_fast_20261015_review.json)
reconstruct the complete frozen queue: 5,172 states, 51,599 contact profiles,
10,085 literal seed/child maps, and 8,280 corner patterns. All 2,070 singular
states remain explicit weight obligations. The source and certificate are
final and frozen. This review does not itself prove a dimension bound.

The arithmetic implementation imports neither the fast producer nor its
tuple-field, minimal-gauge or white-leaf elimination routines. It uses the
older separately reviewed Fraction implementation in the basis
`1,sqrt(5)`, its tree-product Gram gauge, and the original forward finite
Dynkin domain enumeration. This contrasts with the fast producer's reverse
finite-leaf domains and minimal-gauge sparse elimination. Previously
reviewed graph typing, marked-map checking and H4-payment helpers are shared
and identified as such; their old classifications are not re-proved here.

## Independent matrix identities

For an `N`-node core choose the finite corner consisting of all white roots
and black roots0,1. In block notation write

`M = [[A,V],[V^t,B]]`, with `A` of order `N−2`.

The checker performs direct Gauss–Jordan elimination on `A`. Every pivot is
strictly positive at both embeddings, and the final inverse is checked by
literal multiplication. Thus both full inertias are those of
`S = B−V^t A^-1 V`, plus `N−2` positive directions. Only this2×2 matrix uses
the inherited general inertia routine. The usual block inverse gives a
cached inverse of `M` when nonsingular. If `S` has corank one, its kernel
lifts by `(-A^-1 V y,y)`; it is checked against every row of `M`. Deleting
one suitable outside node gives a nonsingular Lorentz ground and a cached
grounded inverse. No white-leaf elimination is used.

The tree-product gauge may change when a new root is attached. The old
block `B0` of the enlarged Gram is therefore checked coefficientwise to be

`B0 = kappa R M R`,

where `R` is diagonal with physically positive entries and `kappa` is
positive at both embeddings. Connectivity of every retained core is
checked explicitly; the diagonal factors are propagated along its actual
edges and then verified on every matrix entry. After this congruence and a
positive scalar division, the new block has old part `M`, contact vector
`u = (kappa R)^-1 v`, and new square `d/kappa`.

For a singular old core, its kernel must annihilate `u`. The reviewer
checks this directly, independently of the producer's encoded relation
products. Otherwise the profile is excluded. If the relation is respected,
the new scalar Schur complement is

`s = d/kappa − u_ground^t M_ground^-1 u_ground`.

For a nonsingular core the same formula uses the full inverse. Its signs
at both embeddings give the new inertias exactly. When a new kernel is
created, `(-M^-1 u,1)` is transported back through the diagonal congruence
and checked against the literal enlarged Gram. Its physical signs test the
nonnegative-root-dependence obstruction. A transported inherited kernel
has the same physical signs because all diagonal factors are positive.

The checker matched nine original-engine states ranging through20 nodes,
including a28-profile case, and separately reconstructed a31-node member
of the frozen subdivision families with inertia `(29,1,1)` and conjugate
inertia `(30,0,1)`. These are development checks, not a substitute for the
terminal replay.

Generation can start from the immutable completed3,473-core snapshot while
the producer continues. Those arithmetic results are retained only in the
same process, keyed by the complete literal record hash. Final use requires
exact equality with the terminal record; no external cache or assumed
producer decision is trusted. The saved replay disallows this prefetch
option and reconstructs every arithmetic case afresh.

## Coverage and database contract

The terminal run is required to have an empty queue, every result present,
and a cursor equal to its total state count. The independent checker
reconstructs all311 literal seed maps, every complete forward contact
domain, every rejection or accepted child, all literal child maps and
reachability. It directly compares the preserved original-engine
[2,007-core snapshot](../data/enumerations/fresh_sqrt5_nonsingular_rank29_exhaustion_20261015_batch2007.json)
with the corresponding final rows and decisions; gauge-specific relation
coordinates and traversal counters are not mistaken for invariant data.

Every singular state, including one with accepted children or at the
31-node boundary, must appear as a zero-further-contact obligation. All
four full corner forests, marked component orbits, existingA1 counts and
dimension-specific spectator rank budgets are recomputed. No older weight
payment is inherited as a reason to drop a case. The final database is
opened read-only and checked for integrity, literal rows/results, cursor,
seed records and its bound hash. All counters and obligation indices are
reconstructed from the records.

The terminal status counts are 2,760 expanded nonsingular states,
132 sealed nonsingular rank exclusions, 210 nonsingular states excluded
at the dimension boundary, 1,871 sealed singular obligations, 43 expanded
singular obligations, and 156 singular boundary obligations. The queue has
no pending state. There are 9,774 accepted child profiles; the other
41,825 profiles are excluded by exact field, inertia or inherited-relation
tests, or paid by the already reviewed entire-H4 rule. The 311 seed maps
plus 9,774 child maps account for the full map total.

The final primary source, JSON and database are immutable. The database is
cryptographically bound through the primary JSON's explicit SHA256 and
checked directly by this reviewer. The primary Markdown report is not
bound, so it can receive this review link without changing any mathematical
input. The producer's terminal structural audit is distinct from this
independent arithmetic replay; a full saved replay of the fast producer is
not claimed here.

The dimension boundary is a logical use of the inherited compact arithmetic
upper bound29. A quadrilateral star in dimension `d` has `d+2` nodes and
Gram rank `d+1`. Thus a nonsingular31-node core is impossible for `d≤29`,
while a singular31-node core remains a weight obligation. If there is no
further touching white root, all outside white components are positive
definite and orthogonal to the displayed core; they cannot rescue a
nonsingular full-star rank, and they must be included in a singular
state's possible forest-sensitive payment. Complete one-white continuation
is valid because any touching outside component has a first root adjacent
to the current core.

The stored spectator budgets are for dimensions 28 and 29. The same queue
is a necessary superset for a smaller ambient dimension: retain only states
with `N≤d+2` and recompute the spectator budget `d+2−N`. Such a use still
requires its own weight-payment checks; the recorded dimension-28/29
budgets are not substitutes for that calculation.

This review does not establish geometric realization of the necessary
principal stars, pay their weights, cover the other finite/dotted seed
families, or independently prove a global dimension bound. Those are
separate inputs to the final coverage certificate. Independent generation
and the fresh saved replay both pass on the frozen terminal data. The saved
replay rebuilt all arithmetic without the generation process's prefetch
cache. Source syntax, recursive manifests and local report links also pass.

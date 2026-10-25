# Finite57: context exhaustion and bounded-tail mathematics

24 October 2026. This is a source-only coordinating audit of the mathematical
algorithms and formulas in the finite-context and bounded-tail appendices of
`manuscripts/finite_volume57.tex`, SHA256
`c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0`.
It supplements the earlier common-geometry and independent low-eight audits.
It does not claim fresh execution, full current-table/code binding, or a
complete semantic review of every originating worker.

I read the complete finite-context appendix and complete bounded-tail
appendix, reconstructed the context removal/reservation arguments, derived
the window counts below directly from the canonical graphs, and checked
the concentration domains and finite-cap comparisons. The recorded
maxima, minima and preferred-image outputs remain named exact computed
propositions at their stated domains. No worker was imported or run.

## 1. The finite calculation is an upper relaxation followed by actual fibers

The rank-58 spherical alphabet is finite after using the proved common-cost
dihedral class. Representatives are an upper-cost construction, not chamber
realizations. The A12/B12 representative count is distinct from the actual
A12 count in the contextual polynomial. A17 is fixed separately, while
B17/D17 remain residual. Those distinctions are present in the text and
are needed for the subsequent guards.

Every actual multiplicity vector fixes c9 through c17, the A12/D12 split
and the rank-five count inside the printed finite rank inequalities. The
remaining components have ranks in the residual set. Finite coefficient
domination must hold for EVERY actual type and pair, after including the
specified g and h allowances; a list of maximizing forests would not prove
that property. The explicit extra 22/64 per A17 remains outside its coarse
base coefficient. No other allowance can silently absorb it.

The coarse residual weights are additive because all remaining pair
interactions are zero in that upper model. Positive ranks justify both
the Bellman maximum and strong-induction verification by all transitions
and tight predecessors. The finite coefficient inequalities, their
connection to literal corner weights, and zero-interaction assertions are
still concrete finite-input obligations. The abstract induction alone
does not establish them.

The G0/G1 reservations cover the two remaining Boolean families. If the
A2/A2 context occurs, two actual A1 components can be reserved and s=2
is legal. If only A2/D5 occurs, reserve one A1 and take s=1. Choosing j
as the COMPLETE rank-two representative count dominates the actual A2
count, so j(j-1) and j*k5 dominate the two contextual multiplicities.
The remaining rank strictly decreases by those reserved components.
Extra A1 components stay available in the no-rank-two table; they do not
multiply a Boolean witness. Forests outside both contexts use the ordinary
residual maximum. The separately bounded A1/D5 and A1/A2/A4/A10 families
cover their overlap as well as their individual domains.

## 2. Actual-fiber refinement fixes every remaining interaction

Expanding each of the 39 costly representative classes into sorted actual
type choices is exhaustive for repeated components, including mixed
assignments. At each such fixed multiset, all A1 and A2 counts are then
enumerated over their full remaining-rank intervals. The A1/A2/A4/A10
family is excluded from this complement precisely when its separate
calculation applies.

The remaining A3-presence guard is split into no A3 and one reserved A3.
The latter may have further A3 components but charges its D14 presence
bonus once. This is a Boolean condition, rather than a separate charge
per witness. The residual A4/A12 context is added linearly exactly when
A1 is present. All other contextual terms are fixed.

There is an additional zero-interaction premise worth making explicit in
the final binding checker: residual types must have ordinary interaction
zero with the separately fixed A1 and A2, not just with each other. The
printed ordinary roles support this: their large hosts lie in the fixed
rank9--16 set or fixed A17, or use fixed D5. The contextual A4/A12 term
and the separately excluded A1/A2/A4/A10 family are handled explicitly.
Nevertheless this should be checked as its own universal finite-alphabet
condition; residual/residual zero is not logically the same statement.

With those stated premises, removing the fixed components, all A1/A2,
and a required single A3 gives the displayed residual rank and cost.
Conversely every reachable residual multiset reconstructs a forest in
that branch. This proves the refinement's exactness and completeness.
All count loops and rank states are bounded by 58.

The two separate exceptional families are upper relaxations. Permitting
extra A1 or A4 in their residual table only enlarges the domain. The
reported 759/539/77 conditional outputs and the 1549 residual tables
are computed premises, not execution performed in this audit. Their
correct application requires the same table binding as the complement.

## 3. The comparison holds at every bounded-tail dimension

For n=2m, the simplex figure has G=m(m-1); for n=2m+1 it has G=m^2.
Substitution into (5G-(n-1)n)/2 gives exactly

    beta_n=n(n-6)/8                    when n is even,
    beta_n=(n-1)(n-5)/8               when n is odd.

The generic-function argument supplies two distinct extrema. Each gains
at least n-1 in the uncut-edge count, since its endpoint h coefficient
is one. Thus the weaker tail correction Delta_n>=5(n-1)>0 is justified
without assuming a finite extremal vertex. Finite cap U_n<=beta_n and
positive ideal surplus therefore suffice; no dimension decrement or
arbitrary Coxeter-face inheritance is used.

For a concentration move, convex order decreases normalized curvature.
The envelope condition is separately c(a-1)+c(b+1)>=c(a)+c(b). It makes
the normalized margin nonincreasing. The rank-square potential strictly
increases while remaining bounded, so the process halts. Absolute
surplus is lifted only after terminal positivity is established.

## 4. Direct graph derivation of the all-rank window counts

For the ideal-five rule at order m>=18, the cycle contributes 5m pairs.
The C spine contributes 5m-15 short pairs. At each additional distance
k=6..9, the first and last permitted windows contain seven starts each;
their union has min(m-k,14) positions.

In affine B, the spine has m-2 nodes. Its short pairs plus two fork
leaves and their mutual pair total 5m-14. The eight additional long
leaf--spine pairs give 5m-6. The remaining spine windows have lengths
six and seven. At original distance k=6..9 they therefore contribute
min(m-(k+2),13), giving the displayed sum over indices 8..11.

In affine D, the spine has m-4 nodes. Short spine, leaf--spine and
same-fork leaf pairs total 5m-13. Sixteen additional long leaf--spine
pairs give 5m+3. The two six-start spine windows contribute
min(m-(k+4),12), giving the displayed indices 10..13. Opposite-fork
leaf pairs are too distant in this order range. No pair is counted
twice when windows overlap: the union size is the minimum of the
available interval and the sum of the window lengths.

These derivations give all four p_X(m), not only their totals at a
few initial ranks. Their branch changes occur at finitely many window
saturation thresholds. Direct substitution gives the maximum 163/23
at m=23. At m>=25 the D count is 5m+51, hence the common envelope
(5r+56)/(r+1), r=m-1. This is a graph-count argument, not a claim
that graph counts alone prove actual face coverage.

For dimensions98--108 the same calculation with short radius eight
and distance-nine feature windows gives

    C: 8m-36+min(m-9,32)=min(9m-45,8m-4),
    B: 8m-33+min(m-11,31)=min(9m-44,8m-2),
    D: 8m-30+min(m-13,30)=min(9m-43,8m).

The cycle count is 8m, so c(r)<=min(8,r/2). The distance-only
families with radius k=9..13 similarly have at most k*m pairs;
small orders have complete-pair bound at most 8m, and larger
classical spines/forks give the stated distance count. The local
replacement theorem at the relevant threshold remains an independent
geometric/finite-certificate input.

## 5. Complete terminal domains, including exceptional repeated ranks

For a linear/cap envelope, low/low moves leave at most one part below
the cap threshold. High/high moves leave at most one part above it.
All other parts can equal the threshold. For low-only concentration,
only the first reduction is used and arbitrary high parts remain.

For the piecewise envelopes in dimensions66--91, low/low moves stop
at the printed linear endpoint; high/high moves stop at rank17.
Each additional nonrepeated exceptional rank j is reduced in pairs
by (j,j)->(j-1,j+1). The printed exact gains are positive. For
example the rank13 gain at dimension66 is 11/65, and the rank12
gain at dimensions67--69 is 23/1092. The rank12 move at66 would
decrease the envelope and is correctly excluded.

Two such exceptional moves may create another exceptional rank, but
this does not cause a termination gap: every move still raises the
same bounded rank-square potential. When no move remains, every
exceptional rank has multiplicity at most one and the low/high
restrictions are exactly those printed in the terminal table.
Thus the sorted-partition recursion covers all stopping profiles.
Its strictly decreasing remaining rank gives termination, and unique
first part and tail give completeness without duplicate profiles.

The coefficient convolution and rational margin calculations are
then the same finite mathematical algorithm as the audited rank57
endpoint. Their numerical output for each interval is a separate
computed proposition. This report does not regenerate those outputs.

## 6. Small tables and layer replacements

The preferred-image algorithm enumerates every induced labelled
injection, checking nonedges as well as edges. Target size at most17
gives finite depth and branching. Every injection has one branch,
so the domain and marked-image union are complete. Preferred removals
need their own actual replacement theorem; a paid-set total does not
supply that theorem.

The dimensions59--65 use ordered SET operations, not accumulating
debits. Full pair arrays, retained-pair lists and patches have different
semantics, and zero is the stated missing-entry default. The D17
patches are accumulated without overwriting an earlier full array.
The D5 long-end/short-leaf distinctions, the two Boolean contexts at59,
and the absence of final58 rank15 guards are explicit. The supplement
selectors identify literal mathematical arrays, rather than requiring
an executing producer merely to define a corner.

I checked the final finite-cap arithmetic. In particular equality
occurs at n=59 with U=783/2, n=61 with420, n=62 with434, n=64
with464, n=66 with495, n=69 with544, n=70 with560 and n=78
with702. The positive global correction still gives a contradiction
at equality. All other printed interval comparisons are nonnegative
with the stated beta_n. These scalar comparisons are also the scope
of an existing compiled conditional tail kernel; it was not rerun.

## 7. What this does and does not close

The context-exhaustion logic, terminal-domain mathematics, window
formulas, count scaling and tail comparison have no identified defect
at their explicitly stated coefficient, local-coverage and numerical
output premises. The required order of operations is preserved.

A full fresh finite57 acceptance still requires precise binding of
the literal simultaneous rule to the aggregate actual/coarse functions,
soundness and completeness of its local receipt interfaces, and the
exceptional finite propositions in the unbounded tail. The independent
local audit is addressing that binding. Existing numerical receipts
remain evidence with their original hashes; no new runtime PASS or
whole transitive regeneration is claimed here. All source/PDF bytes
and prior failures remain unchanged.

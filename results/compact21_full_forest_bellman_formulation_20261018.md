# A complete finite-forest LP without repeated forest cuts

18 October2026. Root mathematical construction for the rank21 golden-field
model. This is an exact extended formulation, not a solved dimension bound.
It is intended to replace support-only forest constraints in the next
bounded research solve. Source implementation and numerical outcome are
separate obligations.

## Additivity after fixing the four hubs

Fix the entire A1 count a and the capped I2(5), H3 and A4 counts
(min(b,2),min(h,1),min(d,2)). Denote this context by c. The completed
four-hub rule supplies, for each actual finite irreducible component T,
an affine weight w_c(T) in the nonnegative bonus variables. The baseline
is its constant term. The normalization of each cross family assigns
its complete vertex contribution to the destination component. Thus
all weights depend only on c, not on the uncapped b,h,d.

In particular the A4-to-other source-root contributions have the factor
1/2 per rooted destination. For d>=2 the A4-to-A4 contribution is linear
in d: each A4 receives one eighth of the ordered root-pair bonus sum.
The denominator4(d-1) cancels the number of other A4 components. These
identities, already used by the exact forest separator, are essential;
arbitrary pair-dependent costs would not admit this formulation.

Let H={A1,I2(5),H3,A4} and let I be the list of actual allowed components
outside H. For an exact hub tuple (a,b,h,d), the fixed contribution is

    L = a*w_c(A1)+b*w_c(I2(5))+h*w_c(H3)+d*w_c(A4),
    r = 21-a-2b-3h-4d.

All remaining components lie in I and have total rank r. There are
components of ranks2 and3 in I, so the reachable residual ranks are0
and every integer>=2; rank1 is omitted.

## Exact Bellman epigraph

For each c introduce U_c(r) for reachable r through the largest residual
rank allowed by that context, with U_c(0)=0. Impose

    U_c(r-rank(T))+w_c(T) <= U_c(r)

for each T in I and every pair of reachable ranks in the range. These
are affine sparse LP rows. Rank strictly increases along each step.
Induction shows that any feasible U bounds the total weight of every
remaining forest. Conversely the exact maximum forest weight at each
rank satisfies all rows. Therefore these potentials describe precisely
the epigraph of all unbounded-knapsack forest maxima, with no missing
multiset or multiplicity assumption.

For every exact hub tuple impose L+U_c(r)<=z. Different exact tuples with
the same c share potentials; their weights w_c(T) are identical. No
rank1 budget is imposed because no corresponding forest exists.

## Removing exactly the two decreasing physical forests

The proved D5 and H4 transfers decrease exactly two physical forests:

    A1+H4+2E8,     D5+2E8.

Their unmodified upper rows must not remain in the model. Only two exact
hub tuples are affected: (1,0,0,0),r20 and (0,0,0,0),r21 respectively.
Each has a unique capped context, so removing its ordinary final budget
does not remove a constraint for another exact hub tuple. Intermediate
Bellman inequalities are retained, including those at these ranks.

For the first exception, every unchanged remaining forest either
contains a component T outside {H4,E8}, or uses only H4 and E8. For each
such off-target T impose

    L+w_c(T)+U_c(20-rank(T)) <= z

when the residual rank is reachable. These rows are exact: they range
over all forests containing that reserved T, and none is the removed
target. The only H4/E8-only rank20 multisets are5H4,3H4+E8,H4+2E8.
Add explicit ordinary rows for the first two and the proved old-cost
minus t row for the target. Here0<=t<=1.

For the second exception reserve each T outside {D5,E8} in the same way
at rank21. The only D5/E8-only rank21 multiset is D5+2E8. Add its proved
active and inactive budget rows: old cost minus1, and old cost minus
twice its extreme E8 endpoint bonus. Their maximum is exactly the
conservative uniform upper used by the reviewed active/inactive rule.

These replacements neither omit an unchanged forest nor impose an
obsolete ordinary bound on a decreased forest. They also avoid an
automaton or enumeration of every exceptional-context multiset.

## Four increasing source forests

Retain all ordinary budgets for the four increased physical forests,
then add their stronger old-cost plus t rows, using the independently
proved A5 capacity one. The forests are A1+A4+2E8,A5+2E8,
A1+B4+2E8,B5+2E8. The extra ordinary budgets are redundant because t>=0.

This gives a finite LP equivalent to all old finite-forest budgets plus
the reviewed conservative transfer allowances. It can be combined with
the51443 exact native face rows. The numerical solver only proposes
bonuses and t; the existing exact complete four-hub forest separator
must re-evaluate any rounded candidate independently of LP potentials.

There is no optimizer result, memory estimate asserted as a measured
fact, Lean verification or dimension improvement in this note. The
additivity premise must be checked against literal forest rows before
the implementation is used. This construction is separate from every
frozen support-only or failed candidate packet.

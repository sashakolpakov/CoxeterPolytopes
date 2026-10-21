# Complete forest formulation with normalized self-hub pools

18 October2026. Mathematical design for a new model, prior to optimization.
The admissible angle-pool lemma, exact old-dual escape and whole-contact-star
coverage are in the frozen `fresh_compact21_self_hub_pool_dual_probe_20261018.md`.
The new exact source is `scripts/compact21_self_hub_pool_model_20261018.py`.
No optimizer result or dimension improvement is asserted in this note.

## The pool charge is a fixed hub-tuple cost

For an exact tuple q=(a,b,h,d), let c=(a,min(b,2),min(h,1),min(d,2)).
The old component weights still use c. The NEW pools use the finer exact
tuple q. Write

  rho(q) = [a>=2]u_A1(q)+[b>=2]u_I25(q)+[h>=2]u_H3(q).

The exact pair-counting lemma says this is the entire new vertex charge.
In particular this is not a cost per hub component. All old component
weights w_c(T), including their rooted cross normalizations, are unchanged.
For the old hub set {A1,I2(5),H3,A4}, let L be the exact sum of old hub
component costs and let r=21-a-2b-3h-4d. The new ordinary budget is

  L+rho(q)+U_c(r)<=z.

The Bellman inequalities defining U_c are precisely the old inequalities:
only residual component weights occur. A pool charge must never be added
while appending each residual item. Although different exact h counts may
share capped c, each has its own final budget and its actual eligibility
indicator in rho, so no count is conflated.

This selector refinement is valid by the same exact counting lemma:
the lemma permits any vertex invariant to select a nonnegative pool.
Setting every fine u_H(q) equal to the corresponding coarse u_H(c)
recovers every coarse-pool rule, so the refined model is at least as
strong. There are at most3*729=2187 eligible new columns. A valid extension
of the old dual to every fine column would also extend to each coarse
column, since a coarse column is the sum of its fine columns. The frozen
coarse-pool forced-load certificate therefore also rules out extension
to all these fine columns, without a new numerical calculation.

For every old special target/source forest rho is zero. Those forests
have at most one A1 and no I2(5) or H3. Thus the two negative physical
exceptions, four positive sources, pure H4/E8 budgets and all their exact
exclusion arguments remain unchanged. The off-target budgets have the
same exceptional exact hub tuples and consequently also have rho zero;
the implementation nonetheless computes rho literally rather than relying
on a hidden special case. Saved old forest support rows receive rho
according to their actual whole hub counts.

## Equivalence to all full finite forests

Fix nonnegative bonus/pool variables and t in[0,1]. Every feasible Bellman
potential bounds each residual multiset, so every ordinary final budget
bounds every corresponding full forest after adding its fixed rho.
The previous reserved-component/pure-target argument covers exactly the
two exceptional physical forest types, with the proved state allowances.
Hence a feasible extended vector gives all stated forest upper bounds.

Conversely, if every such full-forest budget is at most z, choose U_c(r)
as the exact maximum old residual cost at rank r. These maxima satisfy
all Bellman rows by removing/appending a last component. Then all final
budgets hold because their only new term is the fixed rho for that exact
tuple. The exceptional replacement argument is identical. Thus the new
LP is exactly the epigraph of the enlarged conservative forest model.
No independently optimized vertex states, truncation of repeated types,
or finite list of selected forests is substituted.

The exact postcheck computes minimal old Bellman states and adds rho at
each final forest offer. Numerical auxiliary values are discarded. It
also evaluates the six physical exceptions directly with the old D5/H4
rules, checking that none has an eligible new pool.

## Native rows and implementation invariants

Each actual corner uses its full native component diagram plus the exact
common spectator counts. If its marked roots lie in two distinct same-H
components, append coefficient -1/(r_H^2 binom(k_H,2)) to that case's face
row. Here k_H and q are actual counts; only the preexisting bonuses use
the capped selector c.
For a contact row with N roots the corner rank is N-2, so the common
spectator rank in dimension21 is exactly23-N. Removing the four tracked
hub families from those spectators leaves rank zero or at least two:
there is no untracked rank-one type. The enumerated spectator tuples
therefore include every actual case; no choice of an untracked component's
identity affects any new denominator.
The source keeps every one of the102917 cases over1979 native faces and
only then deduplicates exactly equal rational rows. It asserts that
projecting all new face rows to the old bonuses gives exactly51443 old
classes. Old deduplication before the new coefficients would be unsound.

All510 saved forest rows are reconstructed and augmented literally. The
source retains the established68 actual finite types,729 exact hub tuples
and1659 auxiliary
rank states. New pools are appended after old global column94705; z,t and
all old global coordinates retain their meanings. A context has a pool
variable only if an eligible exact tuple exists. Unused new columns can
be eliminated by the same monotone proof as unused old bonus columns:
face coefficients are nonpositive and every other coefficient is
nonnegative. The worker checks these signs during matrix assembly.

The guarded proposed solve uses one thread,45s maximum numerical time,
59s total,1GiB live process-tree RSS, resource start/abort gates, bounded
artifacts and signal-safe owned process-group cleanup. It is not allowed
to run before the source/whole-star coverage review and pending fresh
probe replay are complete. These execution controls are distinct from
the mathematical equivalence proved above.

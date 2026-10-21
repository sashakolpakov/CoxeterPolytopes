# Changing the rigid 12–52 bond after releasing 47

This bounded numerical construction test changes the native 12–52 angle
from π/6 to π/3, π/4 or π/5. It finds more endpoint families, including
nonright A2 profiles, but none of its 1,138 coherent three-cap collections
reduces the intrinsic S6 boundary below the control's 15 open ends.
No replacement chamber is promoted and no nonexistence or finite-volume
claim follows.

The control is the frozen [64-wall chamber](fresh_nonarithmetic_replace47_close12_20261017.md),
with 198 finite vertices and S6 = [0,1,2,3,5,6]. The preceding
[four-duty test](fresh_nonarithmetic_release47_cycle_patch_20261017.md)
kept the rigid π/6 bond. The present test changes that geometry; it does
not apply the preceding finite-profile obstruction to a changed wall.

## Exact deletion bookkeeping

Delete 47 and replace 12, initially keeping 12 right to S6, 7 and 8.
For each of the three new labels at 52, the norm equation has two
branches. One branch has the old strict common center and has exactly
one additional numerical native-angle blocker, old 53. The other branch
fails the center/old-wall screen. The accepted branch's product with 53
is approximately −1.38445547, −1.77415875 or −1.96595132 for labels 3,
4 or 5 respectively. These floating-point observations are not exact
all-label gap certificates.

Releasing 53 permits the three candidates numerically. Its deletion has
no intrinsic S6 cost: the frozen exact projection certificate has
negative square (approximately −10.38177107) and negative product with
the projected common center, so 53 is future timelike and redundant on
the entire future H3. The complete old finite census contains no vertex
incident to 53.

The exact incidence ledger for the changed/deleted set {12,47,53}
removes the following four old finite-vertex certificates:

- S6 + [7,8,12];
- S6 + [7,12,47];
- S6 + [8,12,52];
- S6 + [12,47,52].

The only previously certified compact faces touched are the heptagon
S6 + [8] and quadrilateral S6 + [12]. No higher-dimensional compact
face certificate is touched. This records certificate loss; it does not
assert that an altered face is still compact before checking its new
geometry. In particular, redundancy on S6 alone is never used to infer
preservation of another ambient face.

## Bounded changed-geometry search

The new 12 has norm 2 and prescribed contacts [2,2,m] to [7,8,52],
with m in {3,4,5}. Both norm branches are tested. After deleting 47 and
53 there are 61 unchanged walls and the replacement 12. All unchanged
native inequalities and Coxeter/dotted relations are checked numerically.

The next-stage caps use all seven previously proved necessary finite
S6-facet A2 profiles. Their offsets against every retained native wall,
including their common-center offsets, are included. The finite-anchor
alphabet is {2,…,12,16,24,32}; both real norm branches are tested. No
claim covers omitted labels or contact patterns.

For every rotation of the boundary cycle [7,12,52,30], a seed L is
required to terminate the first two consecutive old edges. Two further
caps P and Q terminate the complementary old edges and meet L along
the respective adjacent finite corner. All old contacts, prescribed
new contacts and the mutual P–Q relation are checked. Every coherent
collection is then passed to the complete numerical S6 triple-intersection
and geodesic-edge ledger, irrespective of whether it realizes the
prescribed disk topology.

| New 12–52 label | Seeds | Coherent P–Q collections | S6 finite/exterior/open ends |
| --- | ---: | ---: | --- |
| 3 | 22 | 469 | 29 / 7 / 15 |
| 4 | 9 | 364 | 29 / 7 / 15 |
| 5 | 7 | 305 | 29 / 7 / 15 |

Every accepted P–Q relation is dotted. There are no finite P–Q joins
in this domain. In particular, six new rotated seeds at label 3 do
terminate [30,7] and [7,12], but their complementary P pool is empty.
The nonzero-A2 seed cases have 27 coherent dotted P–Q frames and remain
a distinct possible fourth-cap target. That successor test is outside
this package.

These computations show that lowering the 12–52 angle opens additional
native-coherent families at the price of deleting 53, but the tested
three-cap collections merely move the open ends. The old exact chamber
remains the authoritative construction; numerical finite counts here
are not an ambient finite-vertex census or an exact promotion.

## Reproducibility and scope

The [source](../scripts/fresh_nonarithmetic_release47_replace12_frontier_20261017.py)
and [saved data](fresh_nonarithmetic_release47_replace12_frontier_20261017.json)
bind the frozen coordinate/projection inputs and deletion census. Exact
assertions concern that inherited projection and integer incidence data.
All new wall equations, native-angle classifications and boundary ledgers
are floating-point discovery calculations. A successful future candidate
would still require exact normal reconstruction, all native gates,
common interior, actual facets, the entire intrinsic boundary and an
explicit ambient volume status.

**FINAL/FROZEN.** Generation and a fresh saved-data replay both terminated
successfully (sessions 95823 and 14286). Source bindings, Python syntax and
local report links pass. This is a reproducible numerical search result
with exact inherited deletion/projection bookkeeping, not an exact
geometric exclusion or a new finite-volume example.

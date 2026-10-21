# Seven cap profiles sharing a simplex coface

This bounded exact pilot proves a new mixed capacity for seven complete
source families, raw23/24/25/33/36/38/39. It records1,645 unambiguous
typed source classes over235 simplex types. Of the62 target simplex rows
already in the current model, three new inequalities are violated by
the frozen12.3313814866 dual. All62 are satisfied by the separate
[raw23-only12.1248285767 witness](fresh_raw23_surviving_coface_witness_20261018.md).
No column has been installed and no optimizer has been run.

The [producer](../scripts/fresh_core1_mixed_coface_capacity_pilot_20261018.py)
and [data](fresh_core1_mixed_coface_capacity_pilot_20261018.json) preserve
the geometry and complete chosen provenance. This is the seven-profile
sole-A1 subfamily of core1, not a census of all245 coface profiles.

## Actual faces and intrinsic roles

The frozen [prism/base incidence proof](fresh_sqrt5_prism_simplex_base_incidence_20261017.md)
supplies an actual compact maximal4-simplex S for the relevant base of
each source aggregate. For the smaller prisms, the original base can
be a proper face of this maximal simplex. The unique connected Lorentz
core of S is core1: five norm2 walls with native bonds
0–1–2–3 of label3 and0–4 of label5. Its native automorphism group is
trivial. All remaining normal components are orthogonal to that core.

In all seven selected source templates, deleting wall b=2 leaves this
core and a single additional whole A1 wall q=0, plus the full spectator
forest. Thus q is an actual facet containing S. Each actual source
aggregate assigns one such actual maximal S and one q. The inherited
full-star classification determines that assignment; orientation maps
are not separately counted. The complete current provenance check
rejects merged source rows and verifies one occurrence at this base
for each of the1,645 selected source classes.

For a fixed actual S and a chosen whole A1 normal q, all the core roles
are intrinsic. With q²=2, the omitted wall has the form

`b = v − sqrt(tau) q/2`.

Its complete contacts with the five core walls determine v in their
nondegenerate Lorentz span. It is orthogonal to every other normal
component of S. Together the core and normal components span the
ambient Lorentz space, so b is unique. Its complete star then recovers
the source aggregate. There is no extra multiplier for the opposite
base endpoint, a native orientation, or the number of quadrilaterals
inside a prism. Those latter demand numbers remain10,6,6,3,3,1,1.

## Exact same-coface exclusions

The producer independently constructs the five-root Gram matrix G over
Q(sqrt5), inverts it, checks the native maps and reconstructs the seven
contact vectors from their literal raw records. Their squared dotted
products tau are

|Raw family|Nonzero core contacts of magnitude1|tau|
|---:|---|---|
|23|3|(7+sqrt5)/2|
|24|1|(23+9sqrt5)/2|
|25|2|6+2sqrt5|
|33|0|20+8sqrt5|
|36|4|14+5sqrt5|
|38|2,4|30+13sqrt5|
|39|3,4|(41+17sqrt5)/2|

For every profile, v²+tau/2=2. For all21 distinct profile pairs at
the same q, the exact checks give

`v_i·v_j <0` and `tau_i tau_j −4(v_i·v_j)^2 >0`.

Hence `b_i·b_j = v_i·v_j +sqrt(tau_i tau_j)/2 >0`.
Distinct outward normals of an acute Coxeter polytope cannot have this
product. The same profile at the same q is the same recovered wall and
aggregate, rather than two occurrences. Thus at most one of the seven
source aggregates can use any given actual q.

For distinct whole A1 copies, q_i·q_j=0 and the product is v_i·v_j.
The28 unordered tests exclude only raw23 against itself, whose product
lies strictly between−1 and0. That is the inherited
[raw23 distinct-copy obstruction](fresh_sqrt5_cross_coface_cap_capacity_20261017.md).
Every other product is at most−1 and is retained. This pilot makes no
additional physical compatibility claim from a retained outcome.

## Complete typed coefficients and scope of the new rows

Let a be the number of whole A1 components in S's actual normal forest.
The same-coface argument gives the necessary count inequality

`sum_{r in {23,24,25,33,36,38,39}} N_r(S-type) <= a N_S`.

There are exactly seven source contexts for each of235 target types.
Each contributes coefficient1, not its orientation count or its prism's
quadrilateral demand. The current provenance map is checked against
the entire context list, so each selected positive source row is
unambiguous even after the historical explicit-provenance splits.

The62 represented target rows yield signed zero-vertex-cost columns.
The other173 targets are saved as absent and uninstalled; four have
positive current source demand. An absent row is not silently treated
as an existing zero-mass simplex variable. The raw23-only witness has
zero source demand at all173 absent targets.

The three represented violations are:

|Simplex face|A1 capacity|Current load minus capacity|
|---:|---:|---:|
|80159|5|15473/400000000|
|80171|3|72969/1250000000|
|2891|2|449/2500000000|

These are genuine necessary inequalities excluding the old dual vector,
but not evidence of an improved optimum. The raw23-only witness
satisfies all62. Independently, its50 source masses are bounded by both
actual simplex-base masses at capacity1, giving100 exact comparisons.
This identifies why a larger scan of local cap alternatives is not
automatically useful: these particular one-simplex restrictions leave
an above12 witness. A stronger relation across actual bases, vertices
or another recoverable feature remains a different possible route.

## Verification

Generation40015/fc53df and fresh replay24334/11b8aa both terminated
exit0/PASS; the independent saved checker5abc17 passed. No old geometry, vertex graph,
unknown outcome or multiple-owner charge is changed. The independent
checker verifies saved exact Gram/inverse and sign arithmetic, all
saved seven-source type blocks, coefficient/load accounting and the
raw23 singleton comparisons. It inherits the historical classification,
actual-face/full-star theorem and completeness of their typed catalogue.
Root's full producer/report/checker mathematical read passed, with the
complete inherited prism/base and cross-coface reports reread. It
checked the actual-face/recovery argument,21 exact same-q exclusions,
separate different-q test, typed multiplicities,62/173 target split and
qualified witness interpretation, without execution or historical
census reconstruction. Construction's independent full producer/report/
checker mathematical review also passed, including both inherited
incidence/capacity reports and generic full-budget interfaces. It
repeated no execution or historical incidence/census. Its requested
explicit saved-list lengths and duplicate-free two-base-per-source
guards were added and the hardened checker957448 passed before freezing.

    python3 scripts/fresh_core1_mixed_coface_capacity_pilot_20261018.py --verify
    python3 scripts/check_raw23_coface_witness_saved_20261018.py

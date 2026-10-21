# A common cap cannot close the two edges30/44 and30/52

For the retained55-wall H9 chamber, **no Coxeter wall can simultaneously
give finite endpoints to the two edges with supports S6+[30,44] and
S6+[30,52]**. Here S6=[0,1,2,3,5,6]. This remains true in every chamber
retaining these55 walls, including the current60-wall construction.

This strengthens the earlier single-wall obstruction for the whole
[8,30,44,52] patch: only two of its four boundary edges are needed. In
particular it rules out the proposed two-cap bridge across the exterior
44/52 edge whose second cap has these two endpoint duties. It does not
exclude other incidence patterns with additional walls.

The [source](../scripts/fresh_nonarithmetic_s6_edges3044_3052_gate_20261017.py)
and [certificate](fresh_nonarithmetic_s6_edges3044_3052_gate_20261017.json)
use the frozen exact physical normals and rational outward intervals. The
[complete endpoint domains](fresh_nonarithmetic_h3_endpoint_domains_20261017.md)
give38 low-label contact profiles and two high-label tails after intersection
of the two necessary spherical conditions. In both tails the contact to44
has arbitrary integer label m≥6; all contacts to S6 and30 are right, and
the contact to52 is either right or label3.

The nine prescribed normals S6+[30,44,52] are independent. If p is the
unique particular vector in their span and k is the unit spacelike vector
in their orthogonal line chosen with k·h<0, every norm-two solution is

    r = p ± sqrt(2-p²) k.

Both branches are tested. Among the76 low-profile branches,45 have a
strictly obtuse old contact and31 have an old contact strictly between
consecutive Coxeter values. These gaps exclude every possible finite label;
the diagnostic angle loop does not impose a hypothesis on the unknown wall.

For each infinite tail put c=2cos(π/m) and p=p0+c p1. The residual norm is
a quadratic polynomial ρ(c), and each old product is

    g_j(c) = α_j c + β_j ± γ_j sqrt(ρ(c)).

Exact interval inequalities over the whole continuous interval
sqrt3≤c≤2 prove ρ>0, strict obtuseness to old4 on the negative branch,
and strict decrease of the old8 product on the positive branch. One
interval suffices for each of the two tails. The endpoint evaluations give:

| Contact to52 | Labels to44 | Strict old8 angle gap |
|---|---|---|
| right | every m≥6 | between labels3 and4 |
| label3 | 6≤m≤10 | between labels4 and5 |
| label3 | every m≥11 | between labels5 and6 |

The second split uses the discrete fact that every integer m is at most10
or at least11; no claim is made about excluding every continuous value
between those two label parameters. Both unbounded integer tails are
covered without an angle cutoff.

All fields are allowed. The conclusion requires finite endpoints: the
retained S6 support contains H4, so its face must be compact in a
finite-volume completion by the frozen
[cusp obstruction](fresh_h9_fixed_core_cusp_obstruction_20261013.md).
No new finite-volume polytope or dimension bound is claimed.

Generation and fresh saved replay passed. A separate agent's read-only
mathematical/source review also passed, checking profile completeness,
both norm branches, and the continuous-tail gap argument. This review
does not constitute an independent arithmetic implementation.

```sh
python3 scripts/fresh_nonarithmetic_s6_edges3044_3052_gate_20261017.py --verify
```

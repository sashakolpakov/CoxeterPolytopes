# Record 31: the entire compact `a=3` schedule is excluded

Status: **certified compact-only no-go**. This is a theorem-track advance, not
a new hyperbolic Coxeter polytope.

Record 31 has missing-edge graph
`M?????????O?O?o??`, namely the path
`11--0--13--1--12` together with nine isolated vertices.  In the compact
Dehn--Sommerville schedule `a=3`, the required face vector is

```text
(14, 87, 319, 768, 1269, 1449, 1110, 522, 116).
```

The already certified positive-A4 archive contains 143 full graph-action
orbits for this schedule.  The new exact completion search covers every one
of them and gives the following chain.

1. There are 378 literal completions having exactly 116 order-five
   transversals. Five have full blocker size 118 and are invalid. The other
   373 have a pure blocker `(5,116)`.
2. The 373 literal targets form 33 orbits under the complete graph action
   `C2 x S9`.
3. Thirty of the 33 orbits fail the exact ridge condition. Exactly three have
   ridge-degree histogram `{2: 522}`.
4. All three ridge survivors pass the complete global selection audit: the
   scheduled face vector, connected dual graph, every-link Euler checks, and
   exact mod-two homology-sphere checks.
5. All three nevertheless have an uncapped empty sparse Coxeter-label join
   on their actual rank-nine vertices. By monotonicity of the nonnegative
   cosine adjacency, adding any omitted finite non-right labels cannot repair
   a failed elliptic principal restriction. Thus each is an unrestricted
   compact Coxeter-label no-go.

The three especially attractive but non-Coxeter minimal-nonface families,
written as vertex sets, are

```text
{{2,3,4,5}, {6,7,8,9,10}, {2,6,11,12}, {3,4,5,13}, {7,8,9,10,13}}
{{2,3,4,5}, {6,7,8,9,10}, {6,7,11,12}, {8,9,10,13}}
{{2,3,4,5}, {6,7,8,9,10}, {5,9,10,11,12}, {2,3,4,13}, {6,7,8,13}}.
```

The independent replay uses a least-destroyer partition rather than the
production chosen-set memoization. It reproduces all 378 scheduled states,
all 373 pure families, and the five impure blockers. Its downstream audit
then checks all 373 literal targets without any symmetry quotient, recovers
the same three ridge survivors, rebuilds the classical 108 order-four and
420 order-five labelled Lanner placements, changes the support-join order,
and again obtains terminal state count zero for all three.

Scope is strict: this excludes only the compact `a=3` schedule of Record 31.
The compact schedules `a=0,1` remain open; a later exact certificate closes
`a=2`. No finite-volume noncompact,
parabolic-cusp, or abstract-polytopality theorem is inferred. In particular,
there is still no new hyperbolic example.

Frozen certificate:
`data/enumerations/h9_14_record31_a3_compact_no_go.json`

Certificate SHA-256:
`ed6fc237d8869e17c135fddae1e61a69712bf93e6e11eae0dcd2e5e9f7227464`.

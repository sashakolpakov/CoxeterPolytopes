# Joint N/R motion closes J and X but reopens M

The four-prefix pilot admits **80 native numerical R walls**, all with
clean target-cone/PD incidence agreement. Every result closes the whole
J,L,X faces, but the M/R/S corner fails because R/S is dotted. Whole M
reopens and S remains open. All80 target cones have24 finite rays and14
open ends. There is no completed disk or finite-volume example.

The [producer](../scripts/fresh_oblique_joint_n_r_pilot_20261018.py),
[saved result](fresh_oblique_joint_n_r_pilot_20261018.json) and
[standalone saved checker](../scripts/check_oblique_joint_n_r_saved_20261018.py)
bind the frozen [forward N packet](fresh_oblique_left_boundary_connector_20261018.md).
They use a separate stem and do not modify its negative fixed-R result.

## Joint motion and bounded selection

Use J67,L68,X69,R70,M71,S72,N73. The forward packet found1,100 N
branches that were native against every prior wall except R. Here R is
removed, one such N is retained, and a new R is solved from X/M/N. Its
genuine prescribed duties are X/M/R and X/N/R. No M/N edge is imposed.
The right-side M/R/S duty is deferred; every retained S product must
still pass the native gate before any cone is scored.

This jointly changes N and R instead of increasing the N label cutoff.
The same disk has boundary8–10–N–R–S–52–8 and V/E/F=10/21/12. Its
only interior vertices are J,L,X,M. The topology and all links are checked
again; no10/X,10/R,J/R,J/L or52/M edge is introduced.

Choose one forward branch closest to a right N/R join for each N/X label
3,4,5, then the next closest branch with a different frozen prefix.
Ties are broken by prefix, recipe and branch index. This gives:

|Forward prefix|N labels in10/J/X|Former N/R product|
|---|---|---:|
|140|3,2,3|0.9014174915|
|133|3,2,4|0.0549621119|
|166|3,2,5|0.0124131367|
|348|3,2,5|0.0124131367|

All four N profiles are0. The selection is deterministic and covers
three distinct N/X geometries; it is not a complete test of the1,100
otherwise-native branches. Their complete eight-step instruction chains
are reconstructed from the frozen inputs, including both earlier R
normals used as anchors. The last prior R is removed after N is built.
Fresh R70 is not an original wall ID, so receives its actual new profile.

The seven-profile label2–5 recipe list is independently compared with
the entire duty-filtered cube at each of the four prefixes. There are
140 recipes:132 have two real norm branches and eight have no real root
in this numerical calculation. All264 real branches are checked.
The four native counts are0,36,22,22. Symbolic high families are saved
but not searched or excluded. The forward packet's high-N argument is
not applied to this changed R geometry.

## Mandatory retained-wall and face checks

The original scaffold has57 projected normals plus six fixed whites.
After deletion30/44/63/64/65,52 projected original normals and six fixed
original whites remain. Alongside the six pre-R new normals J,L,X,M,S,N,
the candidate has58 explicit product gates and fixed-white profile gates:
64 prior facets altogether. Native spacings, the norm and frame residuals,
all new/new offsets and the full common-center product are checked.
Every pre-existing new wall is retested against each admitted R.

Only accepted normals reach the target-cone/PD comparison, all12 prescribed
corner tests, and full ray sets for every recorded face. The first
removed-R prefix has23 finite rays and13 open ends; the other three have
22 and14. The original63-wall control remains26 finite target rays and
ten open ends. None of these target-face counts is the original199
ambient finite-support count.

No additional original facet is removed. The exact debit remains19
finite-support certificates,12 compact two-face certificates and one
compact three-face certificate, including three finite supports on8.
The saved checker independently reconstructs it. This debit is not a
claim that every other old ambient vertex survives the new numerical cuts.

All80 final cones have24 finite rays and14 open ends. Whole compact faces
among the recorded walls are exactly8,11,J,L,X. The complete J,L,X ray sets
equal their intended disk links. M,R,S,N remain open. Exactly one of the
12 prescribed triangles fails in every case:

    M/R/S = (71,70,72).

R/S ranges from−19.6182788525 to−3.0385180469 and is always dotted.
R/M and R/N are finite by the prescribed duties. R/10,R/J and R/52 are
dotted too, but are deliberately absent from the disk.

There are two complete exterior-support patterns. Their ten common
supports are

    (7,47,55), (10,33,N), (33,55,59), (33,55,N), (33,59,66),
    (47,52,N), (47,55,N), (52,S,N), (58,60,66), (59,60,66).

The46-case pattern additionally contains(M,S,N) and(M,N,R). The34-case
pattern instead contains(M,S,R) and(S,N,R). The stored one-open and
two-open edges give14 ends in both cases. Five deduplicated boundary
ledgers include the original control, pre-R prefixes and final cones.

The result replaces the previous two failed J/X triangles by one failed
M/R/S triangle, while adding a boundary wall and preserving the total
14 open ends. It is a genuine numerical joint-normal construction but
does not improve whole-chamber closure. No larger prefix bank or further
connector is executed in this packet.

## Verification scope

Generation21704/chunkadcaad, separate fresh replay84555/chunk69bca0 and
standalone saved checker chunk497e15 all terminate exit0/PASS.

The ideal-lane peer completed a full mathematical producer/report and
standalone-helper source read, including the inherited forward topology
and recipe reconstruction. It checked the deterministic1,100-to-four
selection, both prior R normals in the eight-step chain, original/fixed-white
counts, fresh70 profile offsets, complete declared low cube, all native
and center gates, cone/PD and full-face accounting, and the numerical-only
scope. The review did not execute the programs or repeat norm solving,
cones, intervals, the historical scaffold or spherical classification.
It found no correction. This is an internal agent review.

Final metadata checks pass42 direct bindings, four local links and both
Python sources. The packet is final and frozen.

The standalone checker uses NumPy and the standard library, with no
project imports. It checks all1,100 source-branch candidates and the exact
selection, eight-step saved provenance, old support debit, disk/links,
all140 low recipes and264 branch values, the eight saved no-real statuses,
80 native caps, strict centers, complete recorded face sets and both
exterior patterns. It recomputes stored corner PD tests numerically and
checks saved cone/PD agreement. It does not solve the norm quadratics,
independently prove the no-real statuses, reconstruct geometric cones,
repeat the historical scaffold/classifier or enclose normals by intervals.
No exact new-normal or ambient H9 volume certificate is claimed.

    python3 scripts/fresh_oblique_joint_n_r_pilot_20261018.py --verify
    python3 scripts/check_oblique_joint_n_r_saved_20261018.py

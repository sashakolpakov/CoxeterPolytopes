# Three-cap replacement of the open wall-66 region

Both tested disk rotations fail to close the region. The bounded numerical
search finds 73 caps individually compatible with all 62 retained walls;
every one still leaves ten intrinsic S6 open ends. One or more required
cap banks are empty in each rotation, and all 315 available mutual pairs
are dotted. There is no completed three-cap disk or promoted chamber.
The exact 63-wall, 199-finite-vertex, 10-open-end infinite-volume control
remains unchanged.

The [source](../scripts/fresh_release66_three_cap_disk_20261017.py) and
[saved data](fresh_release66_three_cap_disk_20261017.json) replace wall66
itself, rather than retaining it or returning to the previous
[A/P and R/J closure problem](fresh_cu_native_kr_split7_20261017.md).
They use the frozen [exact control](fresh_nonarithmetic_release49_62_replacement_20261017.md)
and its [complete S6 boundary](fresh_nonarithmetic_release49_62_boundary_20261017.md).

## Different incidence and exact deletion debit

The old wall66 has cyclic intrinsic neighbors

    7–33–30–59–60–58–7.

Its old projective face has two finite vertices, at(7,33,66) and(7,58,66),
and four exterior vertices, at(30,33,66),(30,59,66),(58,60,66),(59,60,66).
The separate exterior component(7,30,47)–(30,47,52) is also retained.

The new caps X=67,Y=68,Z=69 each own two consecutive boundary edges.
There are two rotations, up to cyclic renaming of the three caps:

|Rotation|X frame|Y frame|Z frame|
|---|---|---|---|
|0|(7,33,30)|(30,59,60)|(60,58,7)|
|1|(33,30,59)|(59,60,58)|(58,7,33)|

Each frame gives two required finite boundary corners. Three further
finite corners join adjacent caps at their common old boundary wall;
XYZ must also be finite. The ten listed triangles form a disk with
9 vertices and18 edges. Its boundary is the six-cycle above, its boundary
links are intervals, and its three interior links are pentagons. Both
literal disks are checked in the source. An independent combinatorial
review passed these counts and links, without a geometry claim.

If this disk were realized with precisely its ten finite corners and no
other new incidences, the S6 count would be34 finite vertices and53
geodesic edges, hence four remaining open ends. Those four belong to the
separate B component. Thus the proposed disk targets all six ends of the
wall66 region, rather than merely adding finite vertices. These are
conditional incidence counts, not a realized candidate.

Deleting66 loses exactly these two old ambient finite-support certificates:

    [0,1,2,3,5,6,7,33,66]
    [0,1,2,3,5,6,7,58,66].

The complete inherited 199-vertex incidence ledger finds no certified
compact face of dimension at least two touched by this deletion, and no
compact S6 side is touched. This conclusion uses the full face/edge
incidence test; the fact that66 has only two finite vertices would not
alone exclude compact faces adjacent to66. Retained compact-face
preservation under any replacement still requires exact acuteness and
common-center certificates. No new ambient vertex census is claimed here.

## Numerical domain and acceptance

The inherited finite-corner classifier tests all seven A2 profiles.
It supplies all low labels2 through5 and isolated symbolic high slots;
only6,8,12,24 are substituted in each high slot. The two outer frame
anchors are a nonedge of the proposed disk and are not artificially made
right. In each prescribed frame the native products are converted to
intrinsic products using the actual profile offsets, then both real norm
branches are considered. A Lorentz boost sends the inherited strict
center to its rest frame to reduce cancellation in the large old59/60
coordinates.

Every accepted cap passes the six S6 contacts encoded by its profile,
all56 remaining old native-product tests, its projected norm/contact
residual tests, and the same strict-center test. Unassigned finite
products use the smaller of2e-7 and one hundredth of the neighboring
Coxeter-cosine spacing. Near-parallel products are not rounded to large
finite labels. These are numerical discovery gates, not exact angle
certificates.

For the control, deletion alone, and every accepted single cap, the
complete numerical S6 projective cone is compared with a separate
positive-definite-support census. Open-end counts, finite counts and
the entire finite-support sets must agree without warnings. A completed
triple would additionally require the ten target finite corners and the
whole finite pentagonal ray set of each new facet. No triple reaches
that stage.

|Stage|Count|S6 finite rays|Exterior rays|Open ends|
|---|---:|---:|---:|---:|
|Frozen control|1|26|6|10|
|Delete66|1|24|6|10|
|Accepted single caps|73|26|6|10|
|Completed three-cap disks|0|—|—|—|

The search has153 norm equations,302 real branches and229 rejected
branches. The six symbolic recipe-bank sizes are26,4,5 and28,26,13;
their accepted cap-bank sizes are5,0,0 and0,63,5. All accepted caps have
profile0, hence are right to S6. All315 mutual tests in rotation1 are
dotted, with products approximately between−131.408605 and−5.107358; no finite
junction pair exists. There are zero recorded unresolved numerical
cases and zero improvements. An empty unresolved list does not supply
interval coverage of unsampled high labels or other incidence types.

## Full boundary and the observed blockers

Five complete combinatorial ledgers are saved: control, deletion, and
three accepted single-cap types. The exterior supports after deletion are

    [7,30,47], [7,33,58], [30,33,59],
    [30,47,52], [33,58,60], [33,59,60].

The accepted singles change these supports but keep six exterior rays
and ten open ends. Their finite-ray increase only restores the two
finite S6 corners lost on deletion; it is not boundary closure.

Rotation0's frame(30,59,60) has only four finite-profile recipes and no
high slot. One outward branch, with profile0 and labels(3,2,2), fails only
at old58, with native product approximately−0.8093460084, strictly between
the right and label3 values at numerical precision. The other branches
fail further retained native gates. Rotation1's missing frame(33,30,59)
has a recurring group of angle-gap blockers at old50,51,58; several saved
near failures isolate exactly those three walls. This records concrete
obstacles for this incidence. It is not an exact no-go theorem or a
reason to delete those walls without computing their costs.

Changing the frame duties or the retained boundary geometry remains a
different possible route. Increasing only the mutual-cap label cutoff
cannot repair the empty banks in this sampled domain. No all-label,
all-parameter, other-topology or ambient finite-volume conclusion is
asserted. Even a closed S6 face would not alone certify ambient finite
volume.

Primary generation93470 and separate fresh saved replay26361 both
terminated exit0/PASS. A separate checker without producer imports
verified both literal disks, all recorded counts, the73 saved native
product lists, all315 saved-coordinate pair products, the five complete
cone/PD support ledgers, nine direct bindings, five local links and Python
syntax. It did not reconstruct the physical normals, solve their norm
equations or enumerate cones again.

A separate full source and finite-corner-interface read passed the two
boundary duties per frame, the absence of a fictitious outer-anchor
condition, all three junction mappings and central label order, full
native/center checks, complete cone/PD comparison and whole-facet
completion criterion. This review did not rerun numerical arithmetic.
Source, data and report are final/frozen. The exact control and all
predecessor files are unchanged. No owned process remains live.

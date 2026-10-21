# Native K–R contacts after replacing the old 7–R duty

The changed incidence admits native K–R contacts, but none of the sampled
assemblies improves the complete S6 boundary. There are 1,998 coherent R
caps and 8,125 coherent R+J assemblies; all retain ten open ends. Every
accepted J–R pair is dotted. No new chamber is promoted. The exact
63-wall, 199-finite-vertex, 10-open-end infinite-volume control is unchanged.

The [source](../scripts/fresh_cu_native_kr_split7_20261017.py) and
[saved numerical data](fresh_cu_native_kr_split7_20261017.json) use frozen
[pentagonal prefixes](fresh_cu_release47_pentagonal_cap_20261017.md),
but do not reuse their R normals. In particular the six
[exactly rejected near-parallel products](fresh_cu_release47_pentagonal_native_check_20261017.md)
are not accepted candidates here.

## Changed disk and prescribed native closure

Write A=68,L=69,P=70,R=71,K=72,J=73. Subdivide the old interior edge7–R:
replace7AR and7PR by7AJ,AJR,7PJ,PJR. The boundary remains
7–12–52–30–A–7. The new disk has 10 vertices, 22 edges and 13 triangles.
Every boundary link is an interval, every interior link a cycle, and each
interior edge has multiplicity two. In particular J has link
7–A–R–P–7, and7–R is absent. A separate combinatorial review passed
this stellar subdivision; it made no geometric claim.

For each fixed C/A/L/P/K prefix, R is solved in the frame(A,P,K).
The native labels RA,RP,RK are prescribed first. The finite-corner
classifier requires both ARK and PRK, with their actual A2 profiles.
A/P is a nonedge and is never artificially declared right. The projected
normal has the profile's prescribed square; both real numerical norm
branches are considered. This releases the former constraint R·7=0
and directly enforces a native K–R product instead of fitting a label
to an accidental near-parallel product.

J is then solved once per prefix in the frame(7,A,P), requiring finite
7AJ and7PJ. Every J is tested against every accepted R at that prefix.
A finite J–R relation would also have to make ARJ and PRJ finite.
The full cone must contain all 13 target finite corners, and the complete
L/P/R/K/J facet ray sets must be precisely the prescribed finite polygons.
No such completion occurs.

## Bounded numerical domain and checks

The input has 134 nonempty groups indexed by(A branch,P–7 label,L profile,
K profile). One prefix is selected per group by maximizing its previous
K–R product, used only as a ranking. These are 134 distinct events across
24 A branches, with P–7 labels 3 or 4. The frozen input has no qualifying
P–7 label 5 prefix. This is a representative subset of the earlier event
domain, not coverage of every moving C/A parameter.

The complete inherited finite-corner classifier supplies all seven A2
profiles, all low labels 2 through 5, and symbolic isolated high-label
slots. Only 6, 8, 12, 24 are sampled in each high slot. Numerical singular
or tangent norm cases are recorded as unresolved rather than accepted.
The recorded unresolved list is empty; this does not prove that no
unrepresented tangent or other continuous branch exists.

Every accepted new wall passes all 62 retained old-wall and all other
new-wall native-product checks, projected norm/contact residual checks
and the common strict-center test. An unprescribed finite product must
be within the smaller of 2e-7 and one hundredth of the adjacent cosine
spacing. This rejects the previous false 1094 matches. It is still a
floating discovery criterion, not an exact angle certificate.

For every coherent R addition and every coherent R+J assembly, the
complete numerical S6 cone is enumerated and compared with an independent
positive-definite-support census. The open-end count, finite count and
complete finite-support sets must agree, with no cone warning. Dotted
J–R assemblies are included in this full boundary check.

|Stage|Equations or pairs|Real norm branches|Accepted records|S6 finite rays|Exterior rays|Open ends|
|---|---:|---:|---:|---:|---:|---:|
|Native-first R|4,157|8,314|1,998|36|6|10|
|J bank at prefixes with R|1,411|2,822|869|—|—|—|
|R/J pairing|14,499|—|8,125|38|6|10|

All 1,998 R normals are dotted from old7, as the changed incidence permits.
All 8,125 accepted R/J products are dotted, ranging numerically from
about−132.6770 to−2.004028854818. The remaining 6,374 pairs fail the native
product gate. There are zero finite J–R joins, zero completed disks and
zero open-end improvements.

## Complete boundary diagnosis

There are four saved combinatorial ledgers. The four exterior supports
on old 66 remain in every ledger:

    [30,33,66], [30,59,66], [58,60,66], [59,60,66].

Before J is added, the other two exterior supports are either
[7,A,R] and[7,P,R], or[7,A,P] and[A,P,R]. After J is added, they are
either[A,R,J] and[P,R,J], or[A,P,R] and[A,P,J]. Thus the new caps move
the exposed region while preserving its open-end contribution. The extra
finite rays are not evidence of closure. In this sampled domain, enforcing
native K–R transfers the closing failure to the R/J separation.

## Deletion costs and remaining scope

The exact inherited debit from deleting old 47 is two old finite-support
certificates and the compact quadrilateral supported by S6+12. P with
label 3 is exactly old 47, so those control prefixes restore it; P with
label 4 is a changed wall and does not inherit restoration automatically.
All other retained-wall and compact-face preservation arguments remain
conditional on exact norm, acute-native and strict-center certificates.
This package does not count the new ambient finite vertices or certify
ambient closure or volume. The numerical S6 census alone would not prove
finite volume even if it closed.

Other prefixes, unsampled high labels, continuous deformations of A/P/K,
different duties for J, and other disk incidences remain open. The observed
dotted R/J separation is not claimed as an all-label or all-parameter
obstruction. A further search should change the geometry that constrains
both caps relative to A/P, rather than reuse the rejected 1094 contact.

Primary generation 55634 and separate fresh saved replay 9072 both
terminated exit0/PASS. A separate saved-data check, without producer
imports, verified the literal disk, all recorded counts, the four complete
saved cone/PD support ledgers, the retained deletion debit, all 11 direct
bindings, four local links and Python syntax. That check did not repeat
physical projections, norm-root construction or cone enumeration.

A separate complete source and mathematical report read passed the
changed incidence, representative selection, native-first frame with the
A/P nonedge, both finite-corner domains, assembled native/center gates,
whole-face/cone checks, spacing-relative product tolerance, and the precise
deletion and exposed-region scope. It did not independently rerun numerical
roots, cones or the census. Source, certificate and report are final/frozen.
No exact chamber certificate, finite-volume construction or all-parameter
nonexistence conclusion is asserted.

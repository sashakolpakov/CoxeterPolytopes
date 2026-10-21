# A quadrilateral separator for the dotted K/R gap

The bounded numerical search finds no separator wall in850,978 permitted
contact recipes over all32,019 saved dotted-K/R prefixes. The four prescribed
contacts determine a unique projected normal, and none passes its required
norm gate. This is numerical discovery evidence, not an exact or all-label
nonexistence theorem. The63-wall/199-finite-vertex/ten-open-end control
remains unchanged and infinite-volume.

The [source](../scripts/fresh_pentagonal_kr_separator_20261017.py) and
[result](fresh_pentagonal_kr_separator_20261017.json) extend the frozen
[pentagonal scout](fresh_cu_release47_pentagonal_cap_20261017.md). Its six
false finite-K/R records are excluded from the input: only the32,019 records
already classified dotted are used. No predecessor file is changed.

## A different incidence problem

Delete the two proposed triangles P/R/K and R/A/K, and insert

    P/R/S, P/S/K, A/R/S, A/S/K.

The new cap S has quadrilateral link P–R–A–K–P; K and R remain disjoint.
This removes the impossible finite K/R duty rather than trying to round
its old product to a nearby Coxeter angle. The complete abstract disk has
ten vertices,22 edges,13 triangles, Euler characteristic1, the original
five boundary edges, and connected path/cycle links at every boundary/
interior vertex. These exact combinatorial checks do not assert realization.

Unlike the parallel native-first search that changes R and subdivides7/R,
this experiment keeps each entire old prefix fixed and asks whether one
extra S wall closes its K/R gap.

## Complete bounded contact domain and linear norm test

Use the seven inherited native A2 profiles for S and the finite-label
alphabet

    2,3,4,5,6,8,12,24.

For each prefix enumerate every four-label tuple for S against A,K,P,R
whose four adjacent corners A/K/S,K/P/S,P/R/S,R/A/S have positive-definite
native Coxeter diagrams. The inherited classifier checks the full A2 plus
three-root diagram with each profile's native contacts. The two diagonal
bonds A/P and K/R are not declared finite or orthogonal. There are1,435
distinct profile/boundary-bond domains across the32,019 prefixes.

Let F have rows A,K,P,R in the residual Lorentz coordinates, with metric M.
For an S profile, subtract its native A2 offsets from the four prescribed
products to obtain b. A full-rank frame forces

    (F M) S = b.

The source solves this linear system and separately tests whether
S^T M S equals the profile's required projected square. There is no
remaining free coordinate or norm branch in this four-dimensional span.
No rounded parameter scan replaces the four-contact equation.

All850,978 finite-corner recipes were evaluated. No absolute norm residual
was at most2e-7. The nearest saved row is old R record24712/event2938,
profile0, labels(4,3,2,2) against(A,K,P,R). It has computed square
1.9999506362281152 against required square2, a residual approximately
4.93637718848e-5. Related earlier records give the same near miss within
floating error.

No recipe reached the native/common-center or full cone stage. Had one
done so, every old and new contact would be checked, with finite-angle
tolerance capped by one thirty-second of adjacent cosine spacing. The
code also requires the shared strict center, complete numerical ray/edge
enumeration, independent positive-definite support agreement, all13 disk
triangles and the whole intended cap faces. These later gates were not
exercised by a candidate in this run.

## Numerical limits and verification

Generation34739 and fresh separate-process saved replay85926 both terminated
exit0/PASS, reproducing all inputs, domains and residual decisions. The
replay log is `/private/tmp/kr_separator_verify.log`.

The maximum four-frame condition number is approximately3.024194053e8;
the nearest row's is approximately1.291545348e4. No frame exceeded the
source's explicit1e10 unresolved cutoff. Neither that cutoff nor a small
linear residual certifies exact norm separation. Input coordinates,
finite-label values and the linear norm test are floating point, with
no outward error bounds. Thus zero saved unresolved cases is a property
of this numerical algorithm, not a proof of exact exhaustion. Other
finite labels, other prefixes, moving R or additional caps remain open.

A separate complete source/mathematical read passed the cyclic finite-corner
domain, native A2 offsets, unique four-frame solution and separate norm
gate, spacing-aware native checks and prospective complete cone/PD/whole-face
gates. That review did not repeat the850,978 evaluations or treat a numerical
norm miss as an exact exclusion. Source, saved data and report are final/frozen;
direct binding, Python syntax and local-link checks are recorded in the
checkpoint audit.
No exact wall certificate, ambient finite-vertex census, volume theorem
or new example dimension follows from this search.

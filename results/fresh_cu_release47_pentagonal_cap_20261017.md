# Pentagonal K patch with the A–L and L–R incidences removed

**FINAL/FROZEN after generation and fresh numerical replay. No construction
improvement survives the exact native check.** The exact63-wall/199-finite-
vertex/ten-open-end infinite-volume control is unchanged. Six provisional
nine-open-end records are all rejected by the separate
[exact K/R check](fresh_cu_release47_pentagonal_native_check_20261017.md).
No new chamber or finite-volume example is asserted.

The [source](../scripts/fresh_cu_release47_pentagonal_cap_20261017.py)
implements a new incidence problem suggested by the frozen
[A/P-circle diagnosis](fresh_cu_release47_ap_circle_20261017.md).
The earlier four-cap patch removed L–R but still fixed A–L; reaching a
finite K–R angle then left K–L in a native angle gap. Here both A–L and
L–R are removed from the required patch. The moving parameter is selected
by A–K contact, rather than by the old A–L event equation.

## Literal disk

The boundary remains `7–12–52–30–A–7`. Starting from the nine-triangle
three-cap disk, remove30/A/L, A/L/R and P/L/R. Insert

\[
30/A/K,\quad30/L/K,\quad L/P/K,\quad P/R/K,\quad R/A/K.
\]

The six retained triangles are7/12/P,12/52/P,52/P/L,52/30/L,A/7/R
and7/R/P. The source checks nine vertices,19 edges,11 triangles, Euler
characteristic1, exactly the five boundary edges, edge multiplicities and
all connected vertex links. A–L and L–R are absent. K's link is the cycle
`30–L–P–R–A–30`, so its intended face is a pentagon. P,L,R have intended
face sizes6,4,4. These are exact abstract disk checks, not geometric
realization claims.

## Fixed K and moving A equations

The input is the frozen [quadratic release47 package](fresh_cu_release47_quadratic_events_20261017.md):
three P seeds and98 native-coherent L seeds. P is right to12 and52 with
7-label3,4 or5. L was solved in30/52/P using finite30/52/L and52/P/L.
The P3 seed equals old47 and is an unchanged-wall control in this family.

For each fixed P,L, solve K once in `(30,L,P)`. The mandatory finite
corners are30/L/K and L/P/K. The inherited seven-A2-profile classifier
uses the actual30–L and L–P bonds;30–P is not declared right or required
to be a finite edge. Both norm branches are tested with the native S6
offsets restored. All retained old products, L/P products and the common
strict center must pass. This produces466 native-coherent K seeds in the
completed bank stage.

The116 saved moving-A branches have

\[
A(t)=a_0+t a_++a_-/t,\qquad
u=(t-t^{-1})/(2\gamma),\quad t\ge t_0\ (u\ge3).
\]

The finite30/A/K corner classification supplies allowed target A–K
labels. For fixed K and target \(-\theta_m\), the native A–K equation is

\[
(a_+\cdot K)t^2+
(a_0\cdot K+o_{AK}+\theta_m)t+
a_-\cdot K=0.
\]

All returned positive roots in the stated parameter domain are tested;
there is no parameter grid or upper u cutoff. An identically zero equation
is explicitly unresolved. This algebraic reduction does not certify the
floating-point quadratic thresholds, discriminants or near tangencies.

At an event, C,A,L,P,K must satisfy every retained old and mutual native
contact, the strict shared center and norm/contact residual checks. A–L
has no assigned angle, but still must be a permitted native contact.
R is then solved in7/A/P with finite7/A/R and7/P/R, again using the full
profile classification and both norm branches. R must pass all old and
new native gates. The last required contacts are finite A/K/R and P/K/R;
there is no L–R prescription.

High finite-label families at every sampled stage use only6,8,12,24.
Thus the symbolic classifier derives the finite domain, while this search
evaluates a bounded subset of its high-label values. Other high labels
remain open.

## Complete geometry checks and debit

Every event and R addition passing the floating native gate is scored by complete numerical S6
projective-ray/edge enumeration and an independent positive-definite
pair/triple incidence calculation. Finite support sets, finite counts
and open-end counts must agree without warnings. Duplicate saved boundary
ledgers are interned to keep the output compact; this does not skip any
candidate's cone computation.

The intended completion additionally requires all11 literal disk
triangles to be actual finite supports, and the entire L,P,R,K faces to
have exactly their prescribed4,6,4,5 finite rays with no ideal or exterior
rays. A compact K face alone is insufficient. Native-coherent R additions
which miss this topology are still scored, so an unexpected reduction in
the complete open boundary can be detected.

The only removed control wall is47. Its exact old-incidence debit is two
finite supports, `S6+[7,12,47]` and `S6+[12,47,52]`; the old12
quadrilateral is the only previously certified compact face of dimension
at least2 touched by deleting it. All other62 control walls are retained.
P3 restores the old47 normal, whereas P4/P5 are replacement seeds.
Restoration or change of faces in a numerical candidate must be read from
its full cone, not assumed from this old-certificate debit.

The inherited floating native classifier uses a fixed absolute tolerance.
Near the parallel threshold this can confuse adjacent very large finite
labels. In this run it provisionally accepted six unprescribed K/R
contacts as label1094. They are all false positives: the exact successor
places each product strictly in the gap between labels1094 and1095.
Consequently, passing the floating native gate is not an exact Coxeter
contact certificate, even when the numerical cone and PD checks agree.

No numerical result is promoted without an exact native-normal and
common-center certificate. An improved intrinsic S6 boundary alone would
not establish finite volume for the ambient H9 chamber. No ambient
finite-vertex census, complete ambient exterior cone or new volume claim
is part of this working search.

## Results and status

The [saved result](fresh_cu_release47_pentagonal_cap_20261017.json) records
2,429 K equations and4,858 norm branches, producing466 K seeds (417 with
profile0,27 with profile1 and22 with profile4). All116 moving-A branches
give32,674 A–K quadratics,10,226 positive parameter roots and5,852 events
passing the floating gates. Of these,5,626 use P3 and226 use P4; none use
P5. All event cones have34 finite and6 exterior rays, with10 open ends.

The R continuation evaluates78,822 equations and157,644 norm branches,
with32,025 additions passing the floating gates. The32,019 dotted-K/R
records have36 finite rays,6 exterior rays and10 open ends. The remaining
six records have37 finite rays,5 exterior rays and9 open ends, but their
unprescribed finite K/R contact fails the exact check. No record realizes
the intended disk: only the L quadrilateral is whole and compact in every
record, while P,R,K fail their prescribed complete-face tests. No unresolved
cases were saved within the numerical algorithm's stated thresholds.
This does not certify discarded tangencies or unsampled high labels.

The provisional six are R records18084,18938,21864,22718,25644,26498.
They use P3=old47; L has profile1 or4 and labels(3,2,2) in30/52/P;
K has profile0 and labels(4,2,4) in30/L/P. A has profile0, A7=3,
A30=AK=2, and R is right to7/A/P. The three moving-A families25,27,29
give u approximately3.07185859,3.54231841,3.74852790. Their common local
K/R product is independent of this choice of C(u), so changing u alone
cannot repair these particular branches.

Generation70781 and fresh numerical replay32092 both terminated with
exit0/PASS, reproducing all116 families,5,852 events and32,025 R additions.
The separate exact native reconstruction generation93838 and fresh98435
also terminated with exit0/PASS, rejecting all six provisional improvements.
A no-producer saved-data walk checked the literal disk and connected links,
all five saved boundary ledgers,37,877 event/R records, support/edge/count
identities, whole-face decisions, improvement references and the six strict
rational angle gaps. It did not recompute physical intervals, numerical
roots, cones or PD incidences. Direct bindings, Python AST and local links
passed. The source and numerical output remained unchanged during replay.

Root read the complete source and report: the disk/link domain, finite K
corners without imposing an artificial30/P bond, changed A–K quadratic,
native profile offsets, strict center, complete cone/PD comparison,
all11 corners and whole-face tests, and the exact old47 deletion debit.
The review also checked that the six false positives are explicitly
rejected and zero saved numerical unresolved cases is not treated as
all-parameter coverage. It did not rerun the native equations, interval
engine or numerical cone. No predecessor file was modified.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_cu_release47_pentagonal_cap_20261017.py
OPENBLAS_NUM_THREADS=1 python scripts/fresh_cu_release47_pentagonal_cap_20261017.py --verify
```

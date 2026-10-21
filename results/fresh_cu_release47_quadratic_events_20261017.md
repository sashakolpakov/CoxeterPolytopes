# Release47 bridges from quadratic moving-wall events

The changed-A search finds1,361 native-coherent numerical events and4,437
coherent third-cap additions. **Every full target cone still has ten open
ends.** None is promoted. The exact63-wall/199-finite-vertex/ten-open-end
control remains unchanged and infinite-volume.

The [source](../scripts/fresh_cu_release47_quadratic_events_20261017.py) and
[saved result](fresh_cu_release47_quadratic_events_20261017.json) extend the
[rotated release47 disk](fresh_cu_release47_rotated_bridge_20261017.md) by
changing the moving A recipe, including its A2 attachment profile. They do
not enlarge the preceding parameter grid: the required contact is reduced
to an explicit quadratic equation, with no upper parameter cutoff.

The exact deletion debit is unchanged: old47 participates in precisely two
old finite supports, `S6+[7,12,47]` and `S6+[12,47,52]`. The old12
quadrilateral is the only previously certified compact face of dimension
at least2 touched by deleting it. All62 other control walls are retained.
The three numerical P seeds are right to12 and52 and have7-label3,4 or5;
the label3 seed is old47 itself. Thus it remains a control within the
replacement domain. The new P4/P5 geometry can rebuild the12-face.

## Fixed-L reduction

The boundary disk is `7–12–52–30–A–7`, with P owning the two12 arcs, L
the two30 arcs and R the7/A arc. Its nine triangles and all links are
checked in the predecessor. A complete realization would require the
two finite central corners A/L/R and P/L/R in addition to the endpoint
corners.

For each P, this source first solves L in the **fixed** frame `(30,52,P)`.
The required finite corners are30/52/L and52/P/L. The inherited seven-A2
profile classifier derives the finite label recipes, and both norm
branches are evaluated. High labels are represented only by6,8,12,24.
The138 frame equations produce276 real branches and98 native-coherent
fixed L seeds:37 for P3,32 for P4 and29 for P5. These checks include all62
retained old native products, P and the common strict center.

The [moving companion parameterization](fresh_cu_companion_contact_events_20261017.md)
has

\[
 A(t)=a_0+t a_++a_-/t,\qquad
 u=\frac{t-t^{-1}}{2\gamma},\qquad t\ge t_0\ (u\ge3).
\]

The A2 profile supplies a constant native offset \(o_{AL}\). Once L is
fixed, the target native product \(A\cdot L=-\theta_m\) is equivalent to

\[
 (a_+\cdot L)t^2+
 (a_0\cdot L+o_{AL}+\theta_m)t+
 a_-\cdot L=0.
\]

Both moving A branches for each saved recipe are retained:116 branches in
total. The finite30/A/L classification supplies the permitted target
labels, with the same bounded high-label representatives. This gives16,616
quadratic equations and3,242 returned positive roots in the parameter
domain. An identically zero equation is explicitly treated as unresolved,
not silently discarded; none occurs in this saved run.

At each candidate root the source checks the norm/contact residuals and
every native product among all62 old walls and C/A/L/P, with the shared
strict center. The1,361 accepted events belong to23 moving-A branches:
143 use P3,555 use P4 and663 use P5. Their saved parameter values range
from approximately3.00009553 to9.49565808. That observed range is not a
proven bound on admissible parameters. Each accepted complete S6 cone has
32 finite rays, six exterior rays and ten open ends.

R is then solved in `(7,A,P)`, using the derived finite7/A/R and7/P/R
recipes, both norm branches, full native offsets and all old/new gates.
The6,879 equations give13,758 branches and4,437 coherent caps. Every R–L
product is dotted in the saved numerical data; the nearest is about
−2.00147514. Hence neither prescribed three-cap completion nor a lower
open-end score appears. Every resulting complete cone has34 finite rays,
six exterior rays and ten open ends.

## Whole-cone checks and limits

Every accepted event and third-cap addition is scored by complete
numerical projective-ray/edge enumeration and independently by
positive-definite pair/triple incidence. The finite support sets and
open-end counts agree, with no saved warnings or unresolved records.
The whole8-face has exactly seven finite extreme rays, the whole12-face
four, and the C-face three, with no ideal or exterior ray on these faces.
Thus the numerical replacement restores the12 quadrilateral while
preserving the compact8 geometry; its B-region still escapes.

The exact old-incidence debit and the algebraic form of the parameter
equation are separate from the floating-point root and cone computations.
Quadratic degree thresholds, discriminant cancellation, tangencies,
near-singular norm branches and the shared geometric helpers are not
independently interval-certified here. The source has no parameter grid or
upper cutoff, but this does **not** make it an exact all-parameter
obstruction. Other high labels and moving-A recipes beyond the saved
finite domain also remain open. No new exact ambient normal certificate,
finite-vertex census, complete ambient cone or finite-volume claim is made.

The concrete remaining failure is the R–L adjacency. In the stated disk,
subdividing that internal edge by a fourth cap would replace A/L/R and
P/L/R by four triangles and remove the demand that R and L meet. That is
a different incidence problem, not a conclusion or a tested extension of
this package. The saved changed P/A/L/R geometry makes that future domain
available without reusing the original old47 frame.

**FINAL/FROZEN.** Generation and a fresh saved replay both passed, with
Python syntax, ten direct source/data bindings and report links. Root's
independent complete source/report read passed within the numerical scope:
fixed-L recipes, quadratic coefficients, parameter and identity handling,
native/center gates, cone/PD agreement and explicit coverage limits. It did
not rerun the16,616 equations or4,437 third-cap cones. A separate saved-data
walk, importing no producer, checked all5,798 accepted records,587,416 saved
native products, finite supports, edge end counts and whole8/12/C ray sets;
it did not reconstruct normals, solve roots or enumerate cones. No
optimizer or new exact chamber is used.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_cu_release47_quadratic_events_20261017.py
OPENBLAS_NUM_THREADS=1 python scripts/fresh_cu_release47_quadratic_events_20261017.py --verify
```

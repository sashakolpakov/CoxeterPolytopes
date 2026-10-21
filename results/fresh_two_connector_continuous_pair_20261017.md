# Two continuous zero-arc connectors and a flipped cavity disk

Status: completed bounded **numerical diagnostics**, with generation and fresh saved replay for all six producers. These files do not certify a new Coxeter chamber or an all-parameter obstruction. The authoritative exact control remains **63 walls, 199 finite vertices, ten intrinsic S6 open ends, and infinite volume**.

The new incidence is the chain `A–X–Y–R` along old side8. Both X and Y own no original boundary edge. It escapes the frozen single-X obstruction: no connector has both endpoint duties A and R. Nine numerical four-cap prefixes satisfy all retained native Coxeter/dotted gates and a common strict center, but each has thirteen open ends after deletion. Attempts to complete the original disk and a changed diagonal both fail in the stated finite domains.

## Authoritative input and deletion cost

The input is the frozen [release49/62 replacement](fresh_nonarithmetic_release49_62_replacement_20261017.md) and its [complete intrinsic boundary](fresh_nonarithmetic_release49_62_boundary_20261017.md). S6 is the positive H4+A2 support `[0,1,2,3,5,6]`. The proposed replacement deletes `[12,47,52,64,65]`, leaving the boundary cycle `[7,8,44,63,30]`. The old chord30–44 and its adjacent retained triangle are outside the proposed disk.

The inherited exact deletion ledger removes eleven old finite-support certificates and touches six previously certified compact two-faces, on S6 plus sides8,12,44,63,64,65. This is a loss of old certificates, not a claim that later additions restore them. In particular the compact old8 polygon is not presumed preserved. Each JSON repeats the literal ledger and binds the old exact source/data. No new ambient finite-vertex census is claimed here.

The original six-cap triangulation has five boundary vertices and six new vertices A,R,B,C,X,Y. Its fifteen triangles consist of five old-edge triangles, six old/new/new triangles and four new-only triangles. It has V=11, E=25, F=15, hence disk Euler characteristic1. Geometric existence and embedding are separate requirements.

## Continuous pair equations

The first three producers use zero A2 profile for A,X,Y,R, with X and Y right to old8. For fixed A,R, a connector normal has the form

`X(s)=pX+aX*s+bX/s`, `Y(z)=pY+aY*z+bY/z`, with s,z positive.

For endpoint contact `X·A=−theta`, the fixed projection is `pX=−theta*A/2`, since old8 and A have Gram2I. The residual square is `2−theta²/2>0`. A future unit timelike kernel vector and an orthogonal unit spacelike kernel vector give both conic sheets. All other required native dotted products and the common-center inequality become quadratic inequalities after multiplication by the positive parameter. Their complete numerical interval intersections are retained; zero/degenerate unresolved cases are not accepted as exclusions.

For zero endpoint contacts, the mutual product is

`a*s*z+b*s/z+c*z/s+d/(s*z)`, where a,d<0 and b,c>0.

Any interior stationary point is a saddle in logarithmic coordinates. Boundary-ray extrema are endpoints or the stationary point of `A*t+B/t`. This gives the numerical rectangle maximum used in the first producer.

For nonzero endpoint contacts, the mutual product adds `e*s+f/s+g*z+h/z+i`. Put `P=a*s²+g*s+c`, `Q=b*s²+h*s+d`. At an interior maximum in z, P,Q<0 and `z=sqrt(Q/P)`. Stationarity in s implies the degree-eight equation

`(e*s²−f)²*P*Q − (s*(P*Q)'−2*P*Q)² = 0`.

The producer checks real numerical roots lying in the parameter interval and evaluates the actual feasible points, along with all boundary candidates. Squaring can introduce candidates but does not make an infeasible point acceptable. Identically zero or degenerate polynomials are unresolved. For unbounded rectangles, the strict lower-end coefficient tests and a<0,b,c>0 establish coercivity: neither a one-variable escape nor simultaneous escape can hide a larger maximum. This analytic reduction was independently read by frontier_literature; the floating root isolation and optimizer were not independently implemented.

All outer optimization is bounded discovery. It is not exhaustive interval coverage. Continuous endpoint angles, and the third producer's A53 relaxation, are not themselves Coxeter angles.

| Producer | Explicit domain | Numerical outcome |
|---|---|---|
| [continuous_pair.py](../scripts/fresh_two_connector_continuous_pair_20261017.py), [JSON](fresh_two_connector_continuous_pair_20261017.json) | Twelve low R seeds, right8/profile0; AX=YR=right; 3,438 outer A evaluations and conic rectangles | Best XY −29.2688565; no finite pair; no unresolved rectangle |
| [angle_deformation.py](../scripts/fresh_two_connector_angle_deformation_20261017.py), [JSON](fresh_two_connector_angle_deformation_20261017.json) | Same R seeds; AX/YR products relaxed continuously through `[0,1.999]`; 4,104 outer evaluations | Best XY −3.9417380; no finite pair; no unresolved rectangle |
| [native53_deformation.py](../scripts/fresh_two_connector_native53_deformation_20261017.py), [JSON](fresh_two_connector_native53_deformation_20261017.json) | Additionally allow A53≤−1, including the finite-contact segment; 4,968 outer evaluations | Best XY −2.5753885; no finite pair; no unresolved rectangle |

The outer A domains, fixed optimizer seeds, bounded iteration counts and every best vector are saved. These negative findings apply only to the evaluated domains. They do not close the high-R two-connector family.

## Nine coherent numerical prefixes

The [prefix producer](../scripts/fresh_two_connector_prefix_witnesses_20261017.py) and [JSON](fresh_two_connector_prefix_witnesses_20261017.json) move to high R–44 labels. They set A right to7,8 and A53=π/3, R right to8,63, X/Y right to8, and XY=π/3. The accepted `(R44,AX=YR)` recipes are `(6,5)`, `(7,3)`, `(8,3)` and `(15,2)`, with nine distinct interior parameter choices in total. The tested `(10,2)` recipe produces no saved prefix at its three chosen parameter values.

For example, `(6,5)` has AR≈−71.9167, AY≈−16.1675, XR≈−16.8592, while AX=YR=−phi and XY=−1. Thus the two-connector incidence has a nonempty numerical domain despite the complete single-X exclusion. The output checks every old/new native product, all six new/new products, the norm equations, the intended frame equations and a common strict center. Each output also records the complete numerical S6 cone and a separate positive-definite-support incidence count.

All nine give **21 finite and nine exterior intrinsic rays, thirteen open ends**, with no numerical cone warnings. The second incidence calculation agrees on21 finite vertices,38 geodesic edges and13 open ends. They are not improvements on the frozen ten-end control, and B,C are not yet installed. No prefix has been promoted to exact geometry or a new control.

## Coupled disk and diagonal flip

The [coupled-disk producer](../scripts/fresh_two_connector_coupled_disk_20261017.py), [JSON](fresh_two_connector_coupled_disk_20261017.json), solves B first from its contacts to30,63,R. It tests R44 in `{6,7,8,9,10,12,15}`, all seven B A2 profiles, B30 labels2..6 and the other B labels2..5. It obtains700 norm equations,1,400 branches and422 native-coherent B normals. Solving the finite AB events then gives3,376 A candidates. Each future candidate fails old9; the other branch has the wrong center. Consequently the downstream X/Y/C search is not reached. This is a finite-event diagnostic, not an all-label AB theorem.

The [flipped-disk producer](../scripts/fresh_two_connector_flipped_disk_20261017.py), [JSON](fresh_two_connector_flipped_disk_20261017.json), changes the incidence: it removes AB and replaces triangles AXB and ABC by AXC and XBC. The new diagonal is CX, with AB now dotted. This is a disk edge flip, preserving the boundary and V/E/F counts; it removes the preceding mandatory finite AB condition.

Starting from the same422 B normals, the flipped source solves Y from8,R,B, X from8,B,Y, A from7,8,X, then C from7,30,A if reached. It checks the actual finite corner involving each assigned pair before retaining a branch. It finds615 native Y branches and888 native XYB branches. All3,552 tested future A events again fail old9; their opposite branches fail the center. No five-cap prefix or complete disk survives. There are no unresolved numerical norm frames. New assigned finite labels are2..5, so isolated high-label tails and other A2 profiles for A,X,Y remain open. This failure does not establish a global two-connector obstruction.

## Verification and scope

All six producer generation runs and fresh `--verify` runs reproduce their saved JSON exactly. Source hashes are embedded in each JSON. A final metadata check validates direct hashes, Python syntax and local links. These checks establish reproducibility of the bounded numerical investigations, not exact arithmetic certification of their geometric output.

The meaningful new finding is a feasible two-connector numerical incidence, followed by explicit native-angle failures when the missing interior is coupled in two different disk triangulations. There is no new finite-volume polytope, dimension claim, or change to the exact63/199/10 infinite-volume control.

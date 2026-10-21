# A missing long mirror changes the discriminant-19 H9 seed

19 October 2026. **New exact written construction and scoped obstruction;
independent review pending.** No search, census, optimizer, checker, Lean
process or TeX build was run. Frozen inputs are unchanged. No compact or
finite-volume example is claimed.

The concrete change is to use a previously omitted integral mirror of norm
2(4+phi), replace three incompatible norm-two walls, and add the missing
Bugaenko H8 wall. The resulting 21-facet Coxeter seed has an actual compact
Bugaenko facet and at least 42 certified finite vertices. Its known spacelike
escape survives. We prove that no number of E8-orthogonal additions retaining
the initial vertex can complete this seed to finite volume. A useful new long
mirror must have E8 contact and height greater than 124(4+phi)/19.

## 1. Why this host, and what the existing search actually covers

Here H19 denotes the **discriminant-19 rank-ten form**, not hyperbolic
dimension nineteen. In coordinates (x,y,z), z in E8, put

    G = [[-2phi,phi],[phi,2]] + E8,
    phi^2=phi+1,  b=4+phi,  D=1+5phi=phi*b,
    p=(2,-phi,0^8),  p^2=-2D,  Gp=(-D,0^9).

It has physical signature (9,1) and positive-definite conjugate. The
archived 22-wall search is explicitly **norm-two only**, complete through
positive control x<=5. Its status does not decide the full integral
reflection group. In contrast, the discriminant-11 and discriminant-55
integral hosts have later full-norm infinite-order chamber symmetries and
are closed construction routes; their old open-status notes are superseded.
The norm131 pilot remains bounded and open. The natural binary unimodular
glue was already identified with an earlier lattice, so it is not a new host.

For this even H19 form the complete primitive norm alphabet, up to unit
rescaling, is {2,2b}. Indeed a primitive reflector has norm 2d and
Gr in d O^10. The adjugate and a primitive-coordinate Bezout identity give
d|det G. The determinant ideal (D) is prime, of rational norm 19. Since
O=Z[phi] is a PID and totally positive units are squares, the two possible
totally positive norm classes are 2 and 2b. Here b is totally positive
and D=phi*b. This is the same determinant argument as in the existing
complete prime-norm host proofs, specialized to the present form.

## 2. The changed root, exact debit, and full compatibility

Let

    r=(1,2,0^8).

Then r^2=2b and Gr=(0,b,0^8), so r is primitive and its reflection is

    S(x,y,z)=(x-y,-y,z).

This matrix and its inverse are integral. Its halfspace is exactly y<=0.
Its height, in the convention h=-(p,r)^2/(p^2*r^2), is phi/4. It is not
one of the norm-two mirrors enumerated in the old archive.

Use the literal old wall indices 0,...,21. Delete walls 10,20,21 and add r.
These deletions are necessary: their y coordinates are 1,phi,phi, so their
products with r are positive. The exact relations to retained walls are

    old10 = old0 + r,
    old20 = old12 + phi*r,
    old21 = old13 + phi*r.

Consequently the resulting 20-wall cone is exactly the old 22-wall cone
intersected with y<=0; the dropped inequalities are redundant after the
cut. This is a replacement, not an invalid append-only continuation.

Every retained wall has y in {0,-1,-phi,-phi^2}. Its product with r is b*y.
For y=0 the new angle is right. Otherwise its normalized negative product
is at least sqrt(b)/2>1, so the mirrors are strictly ultraparallel. Every
other pair retains its archived Coxeter relation. Thus no finite angle,
parallel junction or cusp condition is inferred from a numerical fit.

The archived strict center

    c=p+(0,1/10,(E8^(-1)1)/10)

also satisfies the new wall strictly, since its y coordinate is negative.
Acute orthogonal projection of c onto any wall is timelike and strict on
every other inequality. All 20 inequalities are therefore actual facets.

The initial finite vertex p, with its whole E8+A1 incident diagram,
survives. There is a new finite vertex

    v_new=(1,0,0^8),  v_new^2=-2phi,

incident precisely to the eight E8 walls and r. The old A1 wall is strict
there; every retained positive-height wall has product
phi*(y_i-2x_i)<0. Both vertex claims thus include all retained inequalities,
not merely the positive definiteness of their incident diagrams.

This does not preserve every old finite vertex. Exactly those old vertices
on y<=0 survive as points; any proposed complete incidence ledger must
record the debit and new intersection faces rather than reuse the old one.

## 3. One further wall supplies an actual compact Bugaenko facet

The facet on r has y=0 and intrinsic form <-2phi>+E8. Its ten existing
bounding walls are the eight E8 simples and old9,old11. Append the missing
eleventh wall of the archived complete Bugaenko H8 chamber:

    k=(2+3phi,0;
       6+10phi,12+20phi,18+30phi,15+24phi,
       12+18phi,8+12phi,4+6phi,9+16phi).

Here the spatial semicolon only separates coordinates; all vectors use
the same H19 basis. The E8 part z_k satisfies

    E8*z_k = e_4 + 2phi*e_7.

Thus k^2=2, and its root reflection preserves the H19 lattice. Its complete
pair table against the 20 retained seed walls is:

|Wall|Product with k|
|---|---:|
|old0|-3-5phi|
|old5|-1|
|old8|-2phi|
|old1,2,3,4,6,7,9,11|0|
|old12|-1-3phi|
|old13|-3-4phi|
|old14|-6-9phi|
|old15|-4-7phi|
|old16|-4-6phi|
|old17|-4-7phi|
|old18|-5-7phi|
|old19|-2-2phi|
|r|0|

Every nonzero entry except -1 is strictly less than -2, so all pairs are
right, pi/3, or strictly ultraparallel. The old strict center has
(c,k)=-(83+139phi)/10<0, making all 21 walls actual by acute projection.

Every other wall with y<0 is ultraparallel to r. It cannot bound any
part of r's hyperplane. Consequently the full r-facet is exactly the
11-wall compact H8 chamber, rather than a selected finite subregion of
it. To match its authoritative coordinates, the Bugaenko spatial indices
0,...,7 map to H19 indices (0,7,1,2,3,4,5,6); this takes its three
positive-height roots to old9,old11,k. The archived 41 finite vertices
therefore give 41 actual finite vertices of the 21-wall H9 chamber.
Together with the preserved p off that facet there are at least 42.
This is not asserted to be the complete ambient vertex census.

The old escape still has (v,k)=-5/2-3phi<0. Thus even the entire new
compact facet does not make the ambient chamber finite-volume.

## 4. A whole-family obstruction to binary-only completion

The archive's exact spacelike extreme ray is

    v=(1,-phi/2,-3+4phi,-6+8phi,-9+12phi,
       -15/2+10phi,-6+8phi,-9/2+6phi,2+phi/2,-9/2+6phi),
    v^2=(115-71phi)/2>0.

It remains feasible after the cut: its y is negative. Its nine old active
walls (0,1,2,3,4,5,8,9,13) were all retained, so it remains an extreme ray.
Mixing it slightly with the strict center gives a strict spacelike escape.
The new finite vertices do not establish finite volume.

More generally the binary part of v is exactly half that of p. For every
new wall w orthogonal to the E8 subspace,

    (v,w)=(p,w)/2.

If such a wall retains p, it retains v. Hence **every finite full-dimensional
extension by E8-orthogonal walls retaining p still has infinite volume**:
v remains feasible and can be perturbed toward a strict interior point.
There is no norm, label, integrality or height restriction in this geometric
obstruction. A new compact completion retaining this seed must use a wall
with genuine E8 contact or replace additional inherited geometry.

## 5. A lower height for every potentially useful long mirror

For a long reflector w=(x,y,z_actual) of canonical norm 2b, unimodularity
of E8 and Gr in b O^10 force z_actual=b*z with z in E8 over O. The binary
reflection congruences give y=2x modulo b. Put

    L=(phi*x+2y)/b.

Then L is integral, L=x modulo 2O, and the norm equation is

    b*(L^2+2q(z))=4+phi*x^2,  q(z)=z^T E8 z,
    h=phi*x^2/4.

At the conjugate embedding it gives

    0<=q(z)'<2/b'<2(phi')^2

for x>0; the last inequality follows from b'=5-phi>phi^2.

If z is nonzero, write z=a+phi*d with ordinary E8 vectors, and set
A=a^2, J=a.d, E=d^2. These have A,E nonnegative even integers and J
integral. Then

    T=Tr q(z)=2A+2J+3E
     =2||a+d/2||^2+(5/2)E

is a positive even integer. Suppose T<=12. Then E is 0,2 or4.
For E=0, q(z)'=A>=2. For E=4, T<=12 gives A+J<=0, so
q(z)'=A+2phi'J+4(phi')^2>=4(phi')^2.
For E=2, A+J<=3. If A=0, q(z)'=2(phi')^2. If A=2, J<=1 and
q(z)'>=2+2phi'+2(phi')^2=4(phi')^2. If A>=4, J<0 and
q(z)'>=A>=4. All cases contradict the displayed conjugate bound.
Therefore T>=14 whenever the long mirror has genuine E8 contact.

Tracing the norm equation gives the strict inequality

    14 <= T < 18/19 + 2h/b,

because Tr(2/b)=18/19 and phi'<0. Thus every such long mirror has

    h > 124b/19, approximately 36.7.

In particular a short bounded continuation of only low-height long mirrors
cannot repair the archived escape. This is an all-vector inequality, not
a finite alphabet census. At x=0 the long reflection congruences would
force the whole vector into b times the lattice, with impossible even norm
2/b after division; no extra long height-zero wall was omitted.

For orientation, one explicit E8-contacting long root is obtained from the
archived highest root theta=(2,4,6,5,4,3,2,3):

    w_*=(2+5phi,-2-phi,(5+6phi)*theta),
    w_*^2=2b,  h=(45+74phi)/4.

It is **not** an admissible cap for the retained seed. It pairs with old9,
the root (1,0,phi*theta), by +D. In fact

    w_* = s_old9 s_old0(r),

so it is an inherited reflected mirror across existing walls. This records
a concrete first failed nonorthogonal proposal; it is not relabeled a
compactness improvement. The next actual task is an E8-contacting mirror
which passes every retained native-product inequality, using both norm
classes, or a proved necessity to replace old9 as well.

## 6. Reproducibility and precise status

The new construction and bounds are paper algebra only. The old exact
prefix and its saved escape remain dependencies; no fresh replay or
independent review of the new equations was run. Source pins are:

* `data/candidates/h19_vinberg_prefix.json`:
  `91c16b581e5763d6a283b6818f0e087e941fbb987b6b05a7e283febf1e32d0f1`.
* `src/coxeter_polytopes/search/h19_vinberg.py`:
  `991eb7a6078e5880c8d33778a2f425e5d3d8d29023f61bec772ebcc36d46be38`.
* `results/coupled_discriminant19_candidate.md`:
  `5242d0b8333a357164ae671525882762a58834bebb5430d1a6447f67a66d3b16`.
* `data/known_examples/bugaenko_compact_h8.json`:
  `c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c`.

This is a concrete changed mixed-norm arithmetic seed and a rigorous
obstruction to its entire binary-only completion family. It does not
establish global nonreflectivity, a compact H9 polytope, K3 nef equality,
or a new dimension record. No complete mixed-norm first-height or final
vertex/cusp census is asserted.

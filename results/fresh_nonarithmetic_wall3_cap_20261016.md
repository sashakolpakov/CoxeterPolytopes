# A new escape-cutting wall on the retained geometric H9 chamber

An exact new norm-two cap cuts the saved escape of the retained33-wall octic bridge chamber. It is right to

\[
T\cup\{3\}=\{0,1,2,3,5,6,10,30\},
\]

has angle \(\pi/3\) with wall7, and is strictly dotted to every other old wall. All34 facets are actual. The complete finite-corner update gives **95 finite vertices and257 finite edges**, up from92 and253. One old open finite edge closes, and the three new vertices form a finite triangle.

A new strict spacelike escape survives, so the34-wall chamber still has **infinite volume**. This is constructive geometric progress, not a new finite-volume example or a new dimension. The complete exterior-ray cone was not recomputed.

The [producer](../scripts/fresh_nonarithmetic_wall3_cap_20261016.py) and [explicit certificate](fresh_nonarithmetic_wall3_cap_20261016.json) pass generation and saved replay. They reuse the frozen exact octic arithmetic from the [independent old bridge review](fresh_h9_joint_corner_20261015_bridge_review.md). This producer is not itself another independent implementation of that arithmetic.

## Why this changes the previous endpoint recipe

The frozen [target](fresh_nonarithmetic_escape_target_20261016.md) uses the spherical support

\[
T=[0,1,2,5,6,10,30],\qquad \text{type }H_3+A_3+A_1.
\]

The previous [128-recipe obstruction](fresh_nonarithmetic_two_caps_20261016.md) required the first cap to be right to wall31. Here the new wall is right to wall3 instead. It meets the other escape-support anchor7 at the first nonright Coxeter angle, \(\pi/3\). Thus it is outside the excluded family.

The choice has a geometric reason. The projection of old wall3 to \(T^\perp\) has square \((3-\sqrt5)/4\). A norm-two wall perpendicular to T can form a finite corner on that boundary only if its contact with wall3 is right: every nonright ambient Coxeter contact has squared pairing at least1, whereas \(2(3-\sqrt5)/4<1\). The new cap realizes this right endpoint and cuts the actual saved escape.

## Exact cap recipe and arithmetic

Use the original physical form Q and old normals \(n_i\). Set

\[
I=T+[3,7]=[0,1,2,5,6,10,30,3,7].
\]

The nine corresponding covectors have rank9. Let k be their kernel vector in the normalization stored in the certificate. Solve uniquely for p by

\[
(p,n_i)=0\ (i\in T\cup\{3\}),\qquad
(p,n_7)=-1,\qquad(p,k)=0.
\]

Define

\[
w=\frac{337156+149731\sqrt5}{4691}\,k,
\qquad D=1+\sqrt5,
\qquad \boxed{r=p-\sqrt D\,w}.
\]

All square roots use the positive real embedding. The exact identities are

\[
p^2=-8-4\sqrt5,\quad
w^2=\frac{5+3\sqrt5}{2},\quad(p,w)=0,
\quad D w^2=10+4\sqrt5,
\]

so \(r^2=2\). Both quadratic branches are checked; the negative branch has all old pairings nonpositive, retains the old strict timelike center H, and satisfies \((r,E)>0\) for the old strict escape E. The opposite branch is obtuse to retained walls.

The new scalar representation is a formal pair \([a,b]\), meaning \(a+b\sqrt D\) over the frozen octic field. The coefficients use its saved primitive-element encoding. No claim that this compositum has minimal degree16 is needed. Products use \((\sqrt D)^2=D\); all divisions used by the certificate are by nonzero old-field scalars. A sign of opposite-sign terms is decided by comparing \(a^2\) and \(Db^2\), using the inherited outward rational sign bounds in the base field. This establishes exact signs without expanding another large primitive element.

There are33 new wall pairs: eight right, one \(\pi/3\), and24 strictly dotted. Every unspecified contact is proved dotted; no guessed finite-angle alphabet is used for acceptance. The common center is strict for all34 walls. Its projection to each individual wall is checked directly to be future timelike, zero on that wall and strict at the other33, proving every facet actual.

## Complete finite-vertex update

The new cap's complete finite-neighbor set is

\[
\{0,1,2,3,5,6,7,10,30\}.
\]

A finite vertex on the cap must use eight old finite neighbors. All nine possible eight-element subsets are tested by exact positive-definiteness. Exactly three pass. Their full incident sets, with the new cap indexed33, are:

| New vertex | Full incident set |
|---:|---|
| 92 | 0,1,2,3,5,6,7,10,33 |
| 93 | 0,1,2,3,5,6,10,30,33 |
| 94 | 0,1,2,3,5,7,10,30,33 |

Each physical vertex is reconstructed by projecting H off its nine active roots. Its norm is negative, it is future oriented, and all34 pairings have exactly the claimed zero set and otherwise strict negative signs. The old92 vertices are strictly preserved by the new cap. Together with the independently certified old complete census, this proves the new total95; it is not merely a list of known vertices.

The four new finite edges are:

- old vertex2 to new vertex92, on support \([0,1,2,3,5,6,7,10]\);
- vertices92–93, on support \([0,1,2,3,5,6,10,33]\);
- vertices92–94, on support \([0,1,2,3,5,7,10,33]\);
- vertices93–94, on support \([0,1,2,3,5,10,30,33]\).

Thus one old open finite-corner direction becomes a bounded edge, and the new triangle attaches to the old finite graph. There are257 finite edges in total. Since each full spherical corner is simple of rank9, the count of open finite-corner directions is

\[
9\cdot95-2\cdot257=341.
\]

The connected-component count stays seven. These are finite-vertex, finite-edge and open-germ counts only; they do not assert a complete cusp or hyperideal-ray census for the34-wall chamber.

## A renewed escape and the limitation

Take the inherited32-wall hyperideal ray with active support

\[
[0,1,2,3,5,6,7,12,30].
\]

The checker reconstructs it physically, normalizes \((V,H)=-1\), and verifies its positive norm and all34 nonpositive wall pairings, including both later caps. Then

\[
E_{\mathrm{new}}=V+H/8192
\]

has positive norm, negative pairing with H, and **strictly negative pairings with all34 walls**. Its complete coefficients, exact norm and pairings are stored in the certificate. The segment from H toward this spacelike point crosses an open portion of the future ideal boundary inside every halfspace, proving infinite volume.

The added wall therefore removes a verified old escape and changes finite geometry, but leaves another unbounded region. Further caps or a changed construction are required for finite volume. The inherited Cartan two-cycle \(4\sqrt5\) remains, so the target retains the previously proved obstruction to a quasi-arithmetic finite-covolume completion. No integral-host search or arithmeticity assumption is imposed on r.

```sh
python3 scripts/fresh_nonarithmetic_wall3_cap_20261016.py
python3 scripts/fresh_nonarithmetic_wall3_cap_20261016.py --verify
```

The [independent exact review](fresh_nonarithmetic_wall3_cap_20261016_review.md) passes generation and saved replay. It separately verifies all33 new contacts, all34 facet witnesses, strict preservation of the92 old finite vertices, all9 new-facet support tests and their3 vertices, the257 finite edges, and the renewed strict escape. It uses its own formal extension arithmetic over the inherited frozen octic sign engine and does not import this producer. No frozen predecessor was modified.

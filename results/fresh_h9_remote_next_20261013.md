# Metric endpoint and positive-plane bridge in H9

The previous [85-wall chamber](fresh_h9_remote_caps_20261013.md) has been extended to **88 actual Coxeter facets**, retaining every old wall. A metric-coupled cap creates an actual \(H_4+5A_1\) endpoint; two further caps create two \(9A_1\) vertices joined by a finite edge in a different region. Each insertion cuts the immediately preceding exact escape. **The 88-wall chamber still has infinite volume**, certified by a new strict future spacelike point. These are constructive extensions, not a new finite-volume polytope.

All normals have square two. Let \(\phi=(1+\sqrt5)/2\) and
\[
S=\{\mathrm{seed}1,\mathrm{seed}2,\mathrm{seed}3,
       \mathrm{seed}5,\mathrm{seed}6,\mathrm{seed}7\}.
\]
These six seed normals are mutually orthogonal. Exact algebraic recipes establish equalities; outward-rounded integer intervals at scale \(10^{80}\) establish strict signs. The checkers import no numerical search code.

## The 86th wall closes an actual edge

The preceding construction supplies an orthogonal frame \((P,R,C)\) of squares \((-1,1,1)\), with old caps \(R+C\) and \(R-C\). The frame is orthogonal to \(S\) and to \(U_0=\mathrm{seed}0\), whose square is two. Define
\[
L=\sqrt{2+\phi},\quad b=-\frac{1+L}{2},\quad c=\frac{L-1}{2},
\]
\[
a=-\frac{b(R,w_{17})+c(C,w_{17})}{(P,w_{17})},\qquad
u=-\sqrt{\frac{2+a^2-b^2-c^2}{2}},\qquad
r=aP+bR+cC+uU_0.
\]
The denominator is strictly negative and the radicand strictly positive. The formula gives \(r^2=2\), \((r,w_{17})=0\), \((r,R+C)=-1\), and \((r,R-C)=-L\): the latter two angles are exactly \(\pi/3\) and \(\pi/10\). Its 85 old-wall relations comprise seven right angles, these two finite angles, and 76 strictly dotted pairs.

The new wall, named `nextH4`, cuts the old escape and closes the edge with common walls \(S+w_{17}+(R+C)\). Its new endpoint has diagram \(H_4+5A_1\), with chain
\[
\mathrm{seed}3\overset5{-}w_{17}\overset3{-}(R+C)
\overset3{-}r.
\]
The point is timelike, future oriented, and strict on every nonincident inequality. The [86-wall recipe](../data/enumerations/fresh_h9_remote_next_20261013_caps.json) and [replay certificate](../data/enumerations/fresh_h9_remote_next_20261013_verify.json) also contain its remaining escape.

The [exact frontier check](../data/enumerations/fresh_h9_remote_next_20261013_frontier.json) gives one finite edge and eight open rays at this new vertex. A separate rank-eight \(I_2(10)+6A_1\) face has no ninth-wall finite or ideal endpoint. Every positive-semidefinite subset involving `nextH4` lies in one of two specified positive-definite blocks, so this insertion creates exactly one finite vertex and no positive affine circuit.

A [bounded numerical endpoint search](../data/enumerations/fresh_h9_remote_next_20261013_scout.json) released each of the nine incident walls, tested 766 sparse spherical profiles and 314,268 quadratic branches, and found no fully compatible target cutter. This is not a complete endpoint-extension theorem; its exact finite alphabet and anchor list are saved.

## Two caps in a positive normal plane

The successful next construction uses a Lorentzian span instead of a spherical edge. The eight normals \(S+w_{22}+r\) span signature \((7,1)\); their orthogonal complement \(E\) is a positive two-plane. Let
\[
u=\operatorname{proj}_E(w_{38}),\quad d=u^2,
\qquad
v=\operatorname{proj}_E(o)-\frac{(\operatorname{proj}_E(o),u)}d u,
\quad e=v^2,
\]
where \(o\) is the saved common strict timelike center. The certificate proves \(d,e>0\) and \(q=2-1/d>0\). Set
\[
r_0=-\frac ud-\sqrt{q/e}\,v,
\qquad
r_1=-\sqrt{q/d}\,u+\frac{v}{\sqrt{de}}.
\]
Orthogonality of \(u,v\) gives \(r_0^2=r_1^2=2\), \((r_0,r_1)=0\), and \((r_0,w_{38})=-1\). Both caps are right-angled to all eight defining walls. Their other 155 old-wall pairings are strictly dotted, and both cut the 86-wall escape. Thus they retain all 86 walls while adding a genuine finite-angle constraint.

The new actual vertices are
\[
S+w_{22}+r_0+r_1,\qquad S+r+r_0+r_1.
\]
Both have Gram matrix \(2I_9\), and every nonincident inequality is strict. Their shared eight walls \(S+r_0+r_1\) define the finite edge between them. The [88-wall recipe](../data/enumerations/fresh_h9_remote_next_20261013_bridge_caps.json) and [certificate](../data/enumerations/fresh_h9_remote_next_20261013_bridge_verify.json) retain the projection formulas and all strict inequalities.

The [bridge frontier certificate](../data/enumerations/fresh_h9_remote_next_20261013_bridge_frontier.json) checks all 1,536 relevant subsets containing a new wall. Its maximal positive-definite blocks are the two rank-nine blocks above and six rank-three \(H_3\) blocks involving \(r_0,w_{38}\). Hence there are exactly two new finite vertices and no new affine circuits. Their component has one finite edge and sixteen open rays. This is a local frontier census, not a finite-volume argument.

## Exact remaining limit

All 3,828 wall pairs in the final chamber have certified Coxeter or dotted relations. A common strict timelike center proves every inequality is a true facet: projection onto one norm-two wall remains timelike and strictly satisfies every other inequality, because all cross pairings are nonpositive. The inherited affine census and the two new subset checks find no positive affine circuits. Their absence does not prove compactness or finite volume.

The final saved rational dual-coordinate escape has
\[
q(y)=\frac{4953206992184125711621881
-2553473693175672985717741\phi}{16000000000000000000000000}>0,
\qquad
(y,h)_{G^{-1}}=-\frac{565127893357}{250000000000}<0.
\]
All 88 inequalities are strict. This exact point proves infinite volume; the numerical vertex census is used only to propose it. Its numerical proposal has support \([0,3,5,7,9]\).

One further pair-insertion attempt has a stronger, exact local obstruction. Consider the nine named walls
\[
\{\mathrm{seed}1,\mathrm{seed}2,\mathrm{seed}4,
\mathrm{seed}6,\mathrm{seed}8,w_{16},w_{17},w_{35},r\}.
\]
Every eight-element subset spans signature \((7,1)\). In each positive two-dimensional complement, the [projection certificate](../data/enumerations/fresh_h9_remote_next_20261013_bridge_projection_obstruction.json) exhibits two projected retained normals with positive-definite Gram matrix \(H\) and negative off-diagonal entry. Thus every entry of \(H^{-1}\) is strictly positive. Two nonzero proposed cap normals satisfying those two retained inequalities have strictly positive mutual pairing. **No mutually acute-compatible cap pair can lie in any one of these nine planes, regardless of angle labels.** This does not exclude a single cap, caps from different planes, lower-rank orthogonality, or replacement of retained walls. The named sets came from a numerical incident-set proposal; the theorem concerns their exact spans and does not assert incidence of the deliberately strict escape.

## Replay and review

All five producer-free replays pass:

```sh
PYTHONPATH=. python scripts/fresh_h9_remote_next_20261013_verify.py --verify
PYTHONPATH=. python scripts/fresh_h9_remote_next_20261013_frontier.py --verify
PYTHONPATH=. python scripts/fresh_h9_remote_next_20261013_bridge_verify.py --verify
PYTHONPATH=. python scripts/fresh_h9_remote_next_20261013_bridge_frontier.py --verify
PYTHONPATH=. python scripts/fresh_h9_remote_next_20261013_projection_obstruction.py --verify
```

The [86-wall review](fresh_h9_remote_next_20261013_source_review.json) independently checks the symbolic frame and \(H_4\) identities and repeats the full replay at 160 digits. The [88-wall bridge review](fresh_h9_remote_next_20261013_bridge_review.json) independently checks the symbolic projection/norm/orthogonality identities and reruns all 448 enclosures at 160 digits, nested inside the archived 80-digit bounds. These reviews reuse the interval implementation; neither is claimed to be a second field-arithmetic implementation. The construction scripts are [the endpoint producer](../scripts/fresh_h9_remote_next_20261013.py) and [the bridge producer](../scripts/fresh_h9_remote_next_20261013_bridge.py).

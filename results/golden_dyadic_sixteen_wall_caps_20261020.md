# Two integral caps and a finite sixteen-wall compactness gate

20 October 2026. New paper calculation and exact constant list only.
The finite sixteen-wall cone is NONCOMPACT by the exact witness in
Section5. No worker or cone enumeration has run in this note; the
full integral chamber remains open.

Use the frozen four-edge/M6 packet
`results/golden_dyadic_d4_research_20261020_freeze.json`.
In particular import the exact sealed-face splitting from
`golden_dyadic_d4_fourth_vertex_and_sealed_h3_face_20261020.md`,
SHA256 `d6b148954b09e259d0e0598cc4f754904ba1977b373a3f491196ef87f93f9501`,
and the corrected actual M6 edge/reversal from
`golden_dyadic_m6_actual_short_edge_preparation_20261020.md`,
SHA256 `2afd545c3c3166cd07381949b972fb012280b18ed673a090025b6aac30537dd2`,
with its generation and separate replay certificates. The native compact
M6 lineage remains an explicit inherited input. No frozen file is changed.

## 1. Coordinates and the fourteen known walls

Put O=Z[phi], phi^2=phi+1. Work in the orthogonal decomposition

    (V0; e1,...,e6; a3,s,r),
    V0^2=-phi, e_i.e_j=delta_ij,
    G3=Gram(a3,s,r)=[[4,-2,0],[-2,4,-2phi],[0,-2phi,4]].

Here e1=-c and e2,...,e6 are the original I5 coordinates. All original
coordinates and every normal below are recorded exactly in
`results/golden_dyadic_sixteen_wall_caps_20261020.json`, with coefficients
encoded as rational pairs a+b*phi. The original lattice is
`R5=<-phi>+I5+D4(O)`, where the four D4 coordinates have sum in2O.

Let beta_i=-e_i+e_(i+1) for i<6 and beta6=-e6. The nine M6 walls are

    beta1,...,beta6,
    u=(phi-1;phi,0,0,0,0,0),
    w=(phi^2;phi,phi,phi,phi,phi,0),
    h=(2phi^3;2phi^2,phi^3,phi^3,phi^3,phi^3,phi^3).

These are now ACTUAL walls of the full ambient chamber. The actual
order-two reversal A fixes the H3 summand and interchanges u,w and
beta6,h. Set

    wP=a2-V0=-2omega1
      =-(phi+2)a3/2-phi^2*s-phi^3*r/2.

The remaining five known actual walls are a3,s,r,a2,A(a2), where

    a2=(1;0^6;wP),
    A(a2)=(12phi+7;(6phi+4)^6;wP).

The H3 fundamental weights are defined by omega_i.alpha_j=delta_ij;
they are not claimed integral separately. Their inverse Gram is

    G3^-1=(phi^2/4)*[[3-phi,2,phi],[2,4,2phi],[phi,2phi,3]].

## 2. New caps

The vector M0=(2;1^6) in M6 has square -2phi^-3 and is fixed by A.
In the H3 summand define

    theta=phi*a3+2phi*s+phi^2*r,
    theta^2=4, (theta.a3,theta.s,theta.r)=(0,2/phi,0).

The two new integral reflecting normals are

    W=phi^3*M0-phi*theta,      W^2=2,
    Z=phi^3*M0+phi*wP,        Z^2=1,       W.Z=0.       (1)

Their original coordinates, in order(time;I5;D4), are

    W=(3phi^4;(phi^3)^5;-5phi-2,-4phi-2,0,-3phi^4),
    Z=(2phi^5;(phi^3)^5;-phi^5,-4phi-3,-phi,-2phi^5).

The D4 sums lie in2O, and unit I-coordinates prove primitivity.
For a primitive integral vector of norm1 or2 its reflection preserves
the integral lattice, so these are full-arrangement reflectors.
They have not been proved simple facets of the full chamber, and such
firstness is unnecessary for the finite compactness gate below.

Here is the complete contact ledger beyond the old corner. The old
corner is B6+H3, in the order(beta1,...,beta6,a3,s,r). Omitted entries
are zero:

|new wall|nonzero old-corner products|
|---|---|
|u|beta1:-phi|
|w|beta5:-phi|
|h|beta1:-1; beta6:-phi^3|
|a2|a3:-2|
|A(a2)|beta6:-2phi^4; a3:-2|
|W|beta6:-phi^3; s:-2|
|Z|beta6:-phi^3; a3:-2phi|

For the seven-wall order(u,w,h,a2,A(a2),W,Z), the upper triangular
off-diagonal products are

    u: 0,0,-1,-phi^3,-phi,-phi;
    w: 0,-phi^3,-1,-phi,-phi;
    h: -2phi^4,0,-phi^3,-phi^3;
    a2: -18phi-10,-2phi^3,-3phi^2;
    A(a2): -2phi^3,-3phi^2;
    W: 0.

Together with the norms(2,2,1,2,2,2,1), this checks every one of the
120 unordered wall pairs. Every nonzero normalized product is a native
Coxeter angle3,4,5 or strictly ultraparallel. Exact recomputation from
the coordinates, including all inequalities defining these alternatives,
is a mandatory checker obligation; the ledger is not a numerical test.

## 3. Strict feasibility and bounded affine section

Put epsilon=1/100 and take

    q=(1;epsilon*(6,5,4,3,2,1);
       -epsilon*(5phi^2/4,(5phi+4)/2,(9phi+6)/4)).

All nine old-corner products equal -epsilon. Moreover

    q^2=-phi+epsilon^2*(383+24phi)/4<0.

The seven further wall values, in the order above, are

    -1+6phi*epsilon,
    -phi^3+20phi*epsilon,
    -2phi^4+(12phi^2+15phi^3)*epsilon,
    -phi+(5/2)phi^2*epsilon,
    -(19phi+12)+(257phi+173)*epsilon/2,
    -2phi^4+(47phi+25)*epsilon,
    -2phi^4+(47/2)phi^3*epsilon.

They are all strictly negative, for instance directly using
3/2<phi<5/3. Thus the sixteen halfspaces have a full-dimensional
hyperbolic interior and a consistent future sheet.

There is no hidden affine boundedness assumption. In original time
coordinate a, write D_i=-y_i. The original F4 walls (u=a1,a2,a3,s)
and r imply

    0<=y3<=y2<=y1,
    y1+y2+y3<=y4<=a.

The old I5 walls imply I1>=...>=I5>=0. The original W inequality is

    phi^3*sum(I_i)+(5phi+2)y1+(4phi+2)y2+3phi^4*y4
        <=3phi^5*a,

hence sum(I_i)<=3phi^2*a. Therefore a=0 forces the zero vector,
and the a=1 section of the homogeneous sixteen-halfspace cone is a
bounded Euclidean polytope. This proves pointedness and supplies a
finite exact compactness test, not the sign of every extreme ray.

## 4. What the caps remove, and what remains undecided

The fourteen-wall prefix alone has the exact spacelike escape

    X_M=M0/2, X_P=-omega3/2,
    X^2=(27-13phi)/16>0.

Its original integral multiple is
`4X=(3phi+1;2^5;-3,1-2phi,0,-2phi^2)`.
The new cap cuts it strictly: `W.X=-1+phi^3/2=sqrt(5)/2>0`.

There is a stronger finite section check. At fixed M-coordinate M0/2,
write the H3-negative cone as `-sum(lambda_i*omega_i)`, lambda_i>=0.
The W and Z inequalities become

    phi^2*lambda1+2phi^2*lambda2+phi^3*lambda3<=1,
    (3phi+1)*lambda1/2+phi^3*lambda2+phi^4*lambda3/2<=1.

The nonzero vertices of this three-variable polytope are

    (2/(3phi+1),0,0), (0,1/(2phi^2),0), (0,0,phi^-3),
    (phi^-3,phi^-4/2,0), (phi^-3,0,phi^-5).

Convexity of the positive H3 norm gives its maximum `2phi^-6` at
the last listed vertex. Consequently every point of this fixed section
has total square at most `55/2-17phi<0`. The two caps remove the
entire symmetric midpoint escape section, rather than just X.
This section check is not a substitute for a full cone calculation.

The complete cone gate must control transverse M6 directions too.
The following exact counterexample settles that gate without an
extreme-ray enumeration or a further prefix search.

## 5. A spacelike escape for all sixteen walls

Set

    E_M=(1;phi^-1,phi^-1,phi^-1,phi^-2,phi^-2,phi^-2),
    E_P=-phi^-3*omega1-phi^-5*omega3
       =(phi^-1*(-1/2),-1/2,(4-3phi)/2)

in the H3 simple-root coordinates. The six spatial coordinates sum
to3, just as in the symmetric midpoint section, but their square sum
is larger. The exact M6 and positive H3 squares are

    E_M^2=21-13phi,
    E_P^2=2phi^-6,
    E^2=47-29phi>0.                                  (2)

For example phi<81/50<47/29 proves the last strict sign rationally.
Every signed product is recorded here, in the constant-list order:

|wall|E.wall|
|---|---|
|beta1,beta2,beta4,beta5|0|
|beta3,a3|-phi^-3|
|beta6,h|-phi^-2|
|s,u,w,W,Z|0|
|r|-phi^-5|
|a2,A(a2)|-1|

In particular both new cap inequalities are equalities. The original
integral vector is

    2E=(8-2phi; (2phi-2)^2,(4-2phi)^3;
         1-2phi,1-2phi,3-2phi,5-6phi),
    (2E)^2=188-116phi>0.

Its D4 sum is10-12phi in2O. Its positive original time and its
M6 time coefficient1 fix the sheet convention. Together with the
strict timelike interior point q, its feasibility proves noncompactness:
the feasible segment q+tE reaches the future light cone at a finite
parameter and approaches an ideal direction. Thus the hyperbolic
polyhedron is not compact. No claim about its finite volume is needed.

Indeed E is an extreme ray. Its active walls are
`beta1,beta2,beta4,beta5,s,u,w,W,Z`. The four beta equalities reduce
the M6 spatial coordinates to(A,A,A,B,B,B); u,w then give
A=phi^-1*time and B=phi^-2*time. The s equality gives lambda2=0.
The remaining W,Z coefficient determinant in(lambda1,lambda3) is
`-phi^4/2`, so the nine active rows are independent. This explicit
rank argument is optional for noncompactness but gives a small exact
extreme-ray certificate without any DD enumeration.

Only a tiny fixed-identity replay, rather than complete DD enumeration,
is needed to reproduce this negative answer. It should check all
coordinates, lattice membership/reflections,120 products, q, and the
sixteen inequalities and strict square of E. This closes this finite
prefix. It does not exclude the full root chamber of R5 and does not
authorize an unbounded sequence of added caps.

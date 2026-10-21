# A nonpolar mirror cutting the twenty-five-wall H9 escape

19 October 2026. New exact paper construction, awaiting independent review.
No scientific computation, enumeration, optimizer, Lean job or build ran.
The pinned twenty-five-wall construction and escape reports are unchanged.

**Result.** An explicit primitive norm-two lattice root cuts the escape X
although X's polar line contains no lattice reflector. It passes every
retained wall, preserves both complete compact H8 facets, and supplies one
further actual finite vertex. The resulting twenty-six-wall seed has at
least 84 finite vertices. Finite volume is not established here.

## 1. Fixed inputs and the rank-two pencil

Use the native form and notation of
`compact_h9_discriminant19_polar_caps_research_20261019.md`, SHA256
`792634e935866ed78d95d0a93c04fdd7a8484f3605b5c50a66c1d376c9599f1b`,
and the complete escape table in
`compact_h9_discriminant19_twentyfive_escape_20261019.md`, SHA256
`1cc74516c9e22b7150c64a95b4c4251dec9d3c02f238dbb4c946ef937e37c096`.
In particular

    phi²=phi+1,  b=4+phi,  D=1+5phi=phi*b,
    G=[[-2phi,phi],[phi,2]] + C_E8,
    omega1=C_E8^{-1}e1=(7,14,20,16,12,8,4,10),
    omega1²=14,  M=8phi-5,
    p=(2,-phi,0),  lambda=D/M,
    X=(2,-phi,lambda*omega1).

The E8 index 1 is zero based. The new root is sought in the exact pencil

    s=(2u,-phi*u,t*omega1).

Its norm is 2(7t²-Du²), and its covector is

    Gs=(-D*u,0;t*e1).

The following integral solution requires no finite search:

    u=10+17phi,     t=12+19phi,
    u²=389+629phi,  t²=505+817phi,
    D*u²=3534+5719phi,
    7t²-Du²=1.                                      (1)

Both u and t are positive at the physical embedding. Thus s²=2 and
the norm-two reflection v -> v-(v,s)s preserves the entire integral
native lattice. It is primitive: its coordinates generate both u and t
(phi is a unit and gcd(7,20)=1 in omega1), and (1) implies (u,t)=O.

Explicitly its binary coordinates and sparse covector are

    (x_s,y_s)=(20+34phi,-17-27phi),
    Gs=(-95-152phi,0;t*e1).

This is a different direction from X: s² has the canonical norm-two
class, whereas the reviewed prime-59 argument excludes every integral
reflector on X's polar line.

## 2. Every retained wall is compatible

The essential additional exact identity is

    delta=D*u-M*t=3-phi>0.                           (2)

For any retained positive-control wall h=(x_h,y_h,z_h), the complete
escape table gives

    R_h=z_h,1-M*x_h <= 0.

Using (2), rather than assuming a new numerical payment or root census,

    (s,h)=-D*u*x_h+t*z_h,1
         =t*R_h-delta*x_h < 0.                       (3)

This covers all sixteen positive-control walls of the twenty-five-wall
seed, including the six later walls r,k,w,a,d,z. The remaining nine
height-zero walls are old0 and the eight negative E8 simple roots.
Their products with s are zero except

    (s,old2)=-t<0.

Thus no wall, including an E8 spectator, has been omitted.

All nonzero pairs are strictly ultraparallel, rather than merely acute:

* For a retained norm-two positive-control wall other than old9, the
  literal prefix and the displayed later roots have x_h >= phi. Hence
  (3) is at most -(3-phi)phi=-(2phi-1)=-sqrt(5)<-2.
* The remaining norm-two wall old9 has product -19-28phi<-2. The
  height-zero old2 has product -t<-2.
* For the long wall r, (s,r)=-D*u=-95-152phi. For the long walls w,z,
  their x_h=2phi^4>13 gives (s,h)<-13 by (3), since delta>1.
  These magnitudes exceed the mixed-norm threshold 2sqrt(b)<5.

All other new products are zero. In particular no forbidden finite
angle, parallel pair, or unexamined large label occurs.

The inherited strict timelike center is

    c=p+(0,1/10,C_E8^{-1}1/10).

Since (C_E8^{-1}1)_1=91, its product is exactly

    (c,s)=-2D*u+91t/10=-(808+1311phi)/10<0.          (4)

The old center therefore remains strict on all twenty-six inequalities.
The standard acute projection

    c_s=c-((c,s)/2)*s

is timelike, lies on the new wall, and is strict against every old wall:
all old products with s are nonpositive. Thus the new inequality is an
actual facet. Applying the same projection to any other normal proves
that every retained inequality remains an actual facet as well.

## 3. It removes X and is genuinely additional

Another exact identity is

    7t-M*u=2(phi-1)>0.

Consequently

    (X,s)=-2D*u+14lambda*t
         =2lambda*(7t-M*u)
         =4lambda*(phi-1)>0.                        (5)

The old escape is cut strictly. Equations (3)--(4) also show that the
new mirror passes through the interior of the previous seed chamber.
The acute Coxeter seed is a fundamental chamber for its own generated
reflection subgroup, whose mirrors cannot cross that interior. Thus s
is genuinely outside that previously generated subgroup; it is not an
old-wall reflection word presented as a new cap.

## 4. The two full H8 facets and a further finite vertex survive

The new s is strictly ultraparallel to both cap normals r and w. The
acute projection of the common strict center onto either cap is a
timelike point of its hyperplane in s's retained halfspace. Since s
does not meet that hyperplane, the whole hyperplane lies in that
halfspace. Therefore each complete compact Bugaenko H8 facet, with its
41 actual finite vertices, is unchanged. Their vertex sets are disjoint.
The original finite vertex p satisfies (p,s)=-2D*u<0 and also survives.

There is one additional actual finite vertex, not just an elliptic
subdiagram. Put

    p_s=p+D*u*s.

Then (p_s,s)=0 and

    p_s²=-2D-2D²u²<0.

The point p_s is orthogonal exactly to s, old0, and the seven negative
E8 simple walls other than old2. Every retained positive-control wall
h is strict there, since

    (p_s,h)=-D*x_h+D*u*(s,h)<0.

The omitted simple old2 is also strict: (p_s,old2)=-D*u*t<0.
The eight old normals form 2A1+A6: deleting node1 from the native E8
tree leaves the isolated node0 and the chain7--2--3--4--5--6.
Adjoining the orthogonal s gives the positive-definite rank-nine
diagram 3A1+A6. These nine independent actual walls meet at the
timelike point p_s and all remaining inequalities are strict, so this
is an actual finite vertex. It lies on neither cap and differs from p.

The justified lower bound is therefore 41+41+1+1=84 finite vertices.
It is not a complete ambient vertex census.

## 5. Precise scope

This solves the next-mirror problem for the previously obstructed polar
line by a different integral root, using exact pencil equations and all
twenty-five inherited inequalities. It does not establish that all
spacelike escapes are removed or that the twenty-six-wall chamber is
finite-volume. No full mixed-norm Vinberg prefix, global cone census,
K3 nef equality, nonreflectivity, compact H9 example or dimension record
is claimed. The inherited pair ledger and compact H8 completeness
certificate remain explicit dependencies; no scientific job ran here.

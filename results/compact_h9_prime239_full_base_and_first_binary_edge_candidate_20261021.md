# Golden norm239: the unique first binary mirror and its actual finite germ

21 October2026. **Complete manual supporting candidate, pending full root
and outside-author mathematical review.** No worker, scientific import,
enumeration program, test, Lean or TeX process ran. The finite rows below
come from explicit coefficient intervals and hand arithmetic. No frozen
paper, source, guide or receipt was changed. This proves a first actual
edge of one literal module; it proves neither compactness nor
nonreflectivity and supplies no new polytope example.

## 1. Actual module and inherited boundary

Let O=Z[phi], phi²=phi+1, and take

    D=1+16phi, b=phi D=16+17phi,
    L=E8(O) perpendicular-sum [[-2phi,1],[1,8]].       (1)

Use the ORIGINAL binary basis X,Y. The native E8 basis has diagonal2
and products-1 on02,13,23,34,45,56,67. Physically D>0, D'<0 and
b,b'>0; N(D)=-239, N(b)=239. The latter is prime and coprime to30.
This is the exact c=4 member of the accepted good-prime odd-glue family.

| Direct input | SHA-256 | Interface used |
| --- | --- | --- |
| `results/golden_odd_glue_full_base_and_discriminant_reduction_20261021.md` | `b107f0de4f76da942100e0c82867cd6d4eeee40d833304fee67cebd84d3e41ce` | Full actual Bugaenko base, complete primitive norm alphabet and original binary profiles for this c=4 module |
| `results/bugaenko_facet_binary_extension_and_collar_audit_20261020.md` | `774e87f7abb33680a06879a17a5be9f06e2218bf12a95023bcff5c7046ba1227` | Section5: no mixed reflector meets the whole E8-orthogonal hyperbolic subspace; not merely a local collar width |
| `data/known_examples/bugaenko_compact_h8.json` | `c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c` | Actual Bugaenko calibration and native E8 chamber |

The literal module is admissible of signature(9,1), with positive
definite conjugate. Discreteness/local finiteness of the full integral
reflection arrangement and its ordinary spherical germ recognition are
inherited geometric/arithmetic inputs. Finite covolume of its reflection
group is **not** an input. No transport from a similar form or another
invariant lattice is made.

The full base wall is

    r=X+2phi Y, r²=2b, Gr=(0,D;0).

At its actual E8+r corner use the future ray p0=X. Release r and retain
the native E8 chamber. The outgoing binary axis is

    p(t)=X-t r,
    p(t)²=-2phi(1-Dt²), 0<t<1/sqrt(D).              (2)

The full-base theorem supplies the starting chamber C. The global
mixed-E8 statement means every mirror meeting the open axis is binary.
Primitivity, evenness, the golden PID/unit input and the prime norm239
leave precisely primitive squares2 and2b after unit rescaling.

## 2. All binary profiles, with original reflection covectors

Orient a cutting mirror m so that (p0,m)<0 and its derivative along
the axis is positive. Writing m=(x,y;0), these conditions give y<0.
For the binary block the exact identity is

    2phi m²=D y²-(2phi x-y)².                      (3)

For norm2 put v=-y>0 and u=2phi x-y>0. Then

    m=((u-v)/(2phi),-v;0),
    u²=Dv²-4phi, u-v in2O,
    t=u/(Dv).                                     (4)

Conversely (4) gives an integral norm2 vector, and therefore an integral
reflection on the original even lattice. Its first covector is -u.

For norm2b put the first covector equal to -b w and y=-phi z,
where w,z>0. Reflection integrality and (3) give exactly

    m=((Dw-z)/2,-phi z;0),
    z²-Dw²=4, z-w in2O,
    t=w/z.                                        (5)

Indeed the second covector divided by b is `(w-z)/(2phi)`;
D=1 modulo2 also makes the displayed first coordinate integral.
Thus (5) checks both the original-module coordinates and both reflection
covectors. These are full binary root classes, not just native
finite-label candidates. Mirrors with nonpositive derivative do not cut
the positive axis; the starting wall r is already at t=0.

## 3. An exact integral upper mirror

Set

    v*=3+5phi, u*=17+25phi,
    m*=(3+7phi,-3-5phi;0).                         (6)

Direct expansion gives

    v*²=34+55phi, Dv*²-4phi=914+1475phi=u*²,
    u*-v*=14+20phi,
    Gm*=(-17-25phi,-21-33phi;0), m*²=2.

Its original coordinates and reflection covectors are integral; norm2
also makes it primitive. The collision is

    t*=u*/(Dv*)=(17+25phi)/(83+133phi),
    1-Dt*²=4phi/(Dv*²)>0.                          (7)

Thus it is a genuine finite mirror upper on (2). No firstness has yet
been assumed to obtain this upper.

## 4. Complete norm2 interval: exactly six possible v

For any norm2 cut with t<=t*, equations(4) and(7) give

    t²=1/D-4phi/(D²v²), hence 0<v<=v*<111/10.

At the other embedding,

    |v'|²<=4(phi-1)/(16phi-17)<12/43<(53/100)².     (8)

The middle inequality uses phi>8/5: that rational function is
decreasing. Put v=a+kphi, a,k integers. Since
`k=(v-v')/sqrt5`, the bounds above and sqrt5>11/5 give
`0<=k<=5`. Moreover `a=k(phi-1)+v'` and
`3/5<phi-1<31/50`. The complete integer intervals are therefore

| k | Possible a |
| ---: | --- |
|0|0, which is excluded by v>0|
|1|1|
|2|1|
|3|2|
|4|2 or3|
|5|3|

This is a necessary coefficient box derived from the actual finite
upper, with no all-height cutoff. For its six nonzero entries:

| v | Dv²-4phi | Rational norm | Exact square exclusion or root |
| --- | --- | ---: | --- |
|1+phi|50+79phi|209|14²<209<15²|
|1+2phi|133+212phi|941|30²<941<31²|
|2+3phi|349+561phi|2869|53²<2869<54²|
|2+4phi|532+860phi|944|30²<944<31²|
|3+4phi|665+1076phi|-11|A field square has nonnegative rational norm|
|3+5phi|914+1475phi|7921|Exactly (17+25phi)²; 7921=89²|

If an element of O is a square in O, its rational norm is the square
of an integer. The first five rows are impossible. In the last row
the only field roots are +/-u*, and positivity chooses u*. Its original
parity was already checked in(6). Consequently the sole norm2 cutting
hyperplane at or before t* is m*.

## 5. Complete norm2b interval: exactly five possible w

For a norm2b cut with t<=t*, (5) gives

    t²=w²/(Dw²+4),
    0<w²<=4t*²/(1-Dt*²)=u*²/b<81.                 (9)

For the last strict inequality,
`81b-u*²=382-98phi>0`. Also

    |w'|²<=4/(16phi-17)<20/43<(71/100)².           (10)

Write w=a+kphi. Equations(9)–(10) imply `0<=k<=4`.
Using the same bounds `3/5<phi-1<31/50`, the possibilities are
k1:a0,1; k2:a1; k3:a2; k4:a2,3. At k0 only w0 occurs and is
excluded. At k4,a3, w=3+4phi>9, so it is excluded by(9).
Thus the five remaining rows, without even needing their parity test,
are

| w | Dw²+4 | Rational norm | Exact square exclusion |
| --- | --- | ---: | --- |
|phi|21+33phi|45|6²<45<7²|
|1+phi|54+83phi|509|22²<509<23²|
|1+2phi|137+216phi|1705|41²<1705<42²|
|2+3phi|353+565phi|4829|69²<4829<70²|
|2+4phi|536+864phi|3904|62²<3904<63²|

Every row fails the necessary integer-norm-square test. There is no
norm2b mirror at or before t*, and in particular none tied with m*.

## 6. Firstness, actual endpoint and entire germ

The finite interval proof exhausts every binary mirror that can cut
before or at the known finite upper. Mixed mirrors never meet the
binary hyperbolic subspace by the pinned all-root theorem; pure E8
mirrors contain the axis. Hence no other root class can give an earlier
collision. This proves firstness for the FULL integral arrangement,
not an arbitrary selected reflection subgroup.

More explicitly, let s=21+33phi and

    q=s X-u* Y.

Then q²=-2D, (q,m*)=0 and
`p(t*)=[2phi/(Dv*)]q`. It is future, since `(X,q)=-2phi s-u*<0`.
The old wall r is strict: `(q,r)=-D u*<0`. Every other old base
facet extends orthogonally to r, so its product with p(t) is unchanged
along the axis and stays nonpositive; its strict products stay strict.
The E8 walls remain zero. Every remaining oriented chamber wall has
the same side as at the starting germ until it cuts the axis. The
global mixed exclusion and the complete two-profile interval show
that none cuts earlier, so q is feasible against EVERY original wall.

No additional mixed mirror passes through q. The binary orthogonal
line to a future point is one-dimensional, and the only primitive
binary hyperplane at t* is m*. Therefore the entire local arrangement
at q is the E8 root arrangement perpendicular to this A1, and its
actual simple inward facet germ is exactly **E8+A1(m*)**. Usual
local finiteness and spherical germ recognition identify q as a genuine
finite vertex of the same full chamber. There is no finite-volume
assumption in this local conclusion.

For future continuation, the orthogonal pair q,m* has original
coordinate determinant -2, not an O-unit. Precisely,

    X=(v* q-u* m*)/2,
    Y=((3+7phi)q-s m*)/2.

Thus a later germ map or transverse root must respect this original
index2 glue. The norm79 or unit-coupled formulas cannot be transported
by treating q,m* as a unimodular basis.

**Scope.** The literal norm239 host now has a manually proved unique
first binary edge and full finite endpoint germ, subject to the named
actual-base/arrangement inputs. No executable census was required or
performed. Transverse edges and whole-chamber finite volume remain
unresolved; no period, nonreflectivity, new H9 example, form-class
exclusion or dimension-bound change is claimed.

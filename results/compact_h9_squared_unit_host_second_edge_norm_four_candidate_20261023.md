# A second actual edge in the squared-unit glued nine-dimensional host

23 October2026. Complete source-only mathematical candidate by root,
pending full different-author review and root acceptance. No worker,
enumeration, optimizer, symbolic computation, Lean or TeX job was used.
The first-edge input below is itself pending review. This note proves
no reflectivity, finite volume, dimension bound or new example.

## 1. Original module and precise input

The literal input is
`results/compact_h9_glued_e7_a2_squared_unit_coupling_and_first_edge_candidate_20261023.md`,
SHA256 `f11c9dbbb4eabaca146b47d3ef8188131b0f78649e3fa1c32cb2d6c5e98ff664`.
Its Sections3--7 supply the original module, primitive alphabet and
first full actual germ. They are not transported from the closed
coupling-one or coupling-phi chambers.

Write O=Z[phi], phi^2=phi+1, phi>0, K=Q(sqrt5), and b=2+phi.
On the ORIGINAL module

    L=E7(O) perpendicular <e1,e2,g>,
    Gram(e1,e2,g)=[[2,-1,phi^2],[-1,2,0],[phi^2,0,-2phi]],

put n=phi*g+2e1+e2. The input proves that the chamber has a full
finite vertex on the ray of g, with outward normals

    -f1,...,-f7, -e2, n,
    n^2=4b, n.g=n.e2=0, n.e1=2b.              (1)

The original walls -e1 and -e2 remain chamber walls. Every primitive
reflector, after positive-unit scaling, has square2d with

    d in {1,2,b,2b},    (r.L)/d contained in O. (2)

Conversely, the divisibility condition together with a primitive
positive root makes its reflection integral. The K-form has signatures
(9,1) and (10,0), and the full original integral reflection arrangement
is discrete and locally finite. These are the input boundary.

## 2. A short explicit mirror upper

Release -e2, retaining E7 and n. Parameterize the physical edge by

    Q(tau)=g+tau*e2,  0<tau<sqrt(phi).

Its square is -2phi+2tau^2. It stays in the future component since
its product with g is -2phi. The primitive ORIGINAL vector

    r0=(phi-1)g+phi*e2                         (3)

has square4. Its original covector is

    (r0.e1,r0.e2,r0.g)=(0,2phi,-2).

Dividing by d=2 is integral. Its g coefficient phi-1 is a unit,
so it is primitive. It is orthogonal to E7 and n, and acute to
every already known original simple wall, including -e1 and -e2.
It meets Q at

    tau0=phi-1,                               (4)

strictly before the light point. Thus r0 is a genuine integral
mirror upper, not yet assumed a simple wall. The initial full germ
and local finiteness guarantee an actual first wall no later than
tau0. Every such first collision is finite.

## 3. Complete finite profiles against the retained walls

A first simple wall cannot touch E7. Indeed a finite rank9 diagram
containing E7 and an extra touching node must have E8 as that
component, with n isolated. Its new root has the old norm2 scale
and its E7 projection is the minuscule fundamental weight of square
3/2. That weight is not in E7(O). The direct original E7 summand
admits no ternary cancellation. This is the same explicit
original-module obstruction as at the first edge; it does not
assume a new full stabilizer or an inherited chamber symmetry.

Consequently write every candidate as

    r=x*e1+y*e2+z*g, x,y,z in O.

Retained-wall acuteness gives x<=0 because r.n=2b*x. Put

    A=2z-phi*x, B=2y-x.

Since the collision is at positive tau,

    A>0, B>0, tau=phi*A/B,
    r.g=-phi*A, r.e2=B,
    B^2-phi*A^2+2b*x^2=4d.                   (5)

Crucially, ORIGINAL reflection integrality in (2) gives

    A/d in O, B/d in O.                      (6)

We retain these even though e2 is the released wall.

At a finite endpoint the retained n and arriving r have a finite
native angle. If M=4cos^2(pi/m), then

    x^2=d*M/(2b).                            (7)

The complete M-list in K is {0,1,2,phi^2,3,b}, with labels
2,3,4,5,6,10. To see completeness, the degree of cos(2pi/m)
is at most2, so Euler's phi-function is at most4; the remaining
labels8 and12 would require sqrt2 and sqrt3, outside K.
Equation(7), physical x<=0 and the K-square test give precisely:

| d | all possible x |
|---|---|
|1|0|
|2|0,-1|
|b|0,-1|
|2b|0,-1,-phi|

Here the nonzero entries have labels10,4,3,5 respectively. This
table does not discard a native10 or dyadic class. For completeness,
the excluded square tests reduce either to a nonsquare rational norm,
or to the elementary statement that a positive rational number is
a K-square only if it is a rational square or five times one.
For the latter statement expand (a+c*sqrt5)^2, a,c rational: its
irrational coefficient vanishes only when a*c=0. This excludes,
in particular, 1/2, 2, 3 and 3/2 when they occur in (7).

If d=2 or2b and x is nonzero in the table, (6) makes A divisible
by2. But A=2z-phi*x and each listed nonzero x is a unit modulo2.
This is impossible. If x=0 and d=b or2b, (6) makes BOTH y and z
divisible by b: 2 and phi are units modulo b. Then all original
coordinates of r are divisible by b, contrary to primitivity.

The ENTIRE first-wall domain is therefore reduced to

    x=0, d=1 or2; OR x=-1, d=b.              (8)

This reduction uses original covectors, not only a normalized
Gram matrix or a guessed coefficient box.

## 4. Excluding the remaining nonorthogonal profile

Suppose x=-1 and d=b. From tau<=tau0=1/phi and (5), B>=phi^2*A.
Thus

    B^2-phi*A^2=2b,
    A^2<=b/phi^2=3-phi=b'.                  (9)

Conjugating the equality, using phi'=-1/phi, gives

    A'^2<=2phi*b'.                           (10)

By (6), write A=b*C with C in O. These two bounds imply

    0<C<=sqrt(b')/b<1/2,
    |C'|<=sqrt(2phi/b')<5/3.                 (11)

The first strict comparison follows from 4b'<b^2, equivalent
to 7<9phi. The second follows from phi<5/3:
2phi/(3-phi)<5/2<25/9. These are strict rational comparisons,
not numerical observations.

Write C=a+j*phi with integers a,j. Then

    -5/3 < C-C'=j*sqrt5 < 13/6.

Since 5/3<sqrt5 and 13/6<sqrt5, j=0. But then C=a is an
integer strictly between0 and1/2, a contradiction. This excludes
the only possible nonorthogonal first profile, including its
closed upper endpoint and both conjugate signs.

## 5. Complete orthogonal domain and the unique survivor

It remains to consider x=0, d=1 or2. Here

    r=y*e2+z*g, y>0,z>0,
    y^2-phi*z^2=d,
    tau=phi*z/y<=1/phi.

Hence y>=phi^2*z, so

    2phi^2*z^2<=d,  z'^2<=phi*d.             (12)

For d=1 this gives 0<z<1/2 and |z'|<3/2. Writing z=a+j*phi,
we obtain -3/2<j*sqrt5<2, so j=0. There is no positive integer
z<1/2. Thus no norm2 reflector lies in the closed first domain.

For d=2, (12) gives 0<z<=phi-1<1 and |z'|<2. Consequently
-2<j*sqrt5<3, so j is0 or1. The integer case j=0 is impossible.
For j=1 the physical bound gives -phi<a<=-1, hence a=-1 and

    z=phi-1.

Now y^2=2+phi*(phi-1)^2=phi^2. Positivity gives y=phi, exactly
the original vector r0 in (3). This is a complete terminating
two-embedding coefficient argument. It does not use a census,
square-root implementation or any floating-point rejection.

## 6. Actual firstness, full endpoint and scope

Every possible first simple wall up to the explicit upper was
included, and r0 is the sole survivor. Thus r0 is the unique
actual opposite wall of this edge. Unknown walls do not cut the
segment earlier, by the complete first-wall classification.

Its endpoint has the primitive original representative

    v=phi*g+e2=phi*Q(phi-1), v^2=-4phi.

The nine retained/arriving outward normals are

    -f1,...,-f7, n, r0,
    Gram=E7 perpendicular[4b]perpendicular[4], det=32b. (13)

Their Gram is positive definite. Finite acute simplicity makes
this the FULL actual germ. The released wall is strict since
Q(phi-1).(-e2)=-2(phi-1)<0; the older -e1 is strict since
Q(phi-1).(-e1)=-2. The known wall signs and both actual norm4
scales are retained. No normalization silently turns either
root into the norm2 roots of a closed earlier host.

The new endpoint is an actual ambient nine-dimensional vertex,
conditional only on the precise first-edge input and the stated
external arrangement/finite-diagram foundations. This proves no
compact chamber, finite volume, period or new example. The next
edge must start from (13) and keep all original norm classes,
known walls and actual firstness obligations.

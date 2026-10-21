# Two cubic-seven simple-end quadrilaterals at B19/B20 vertices

19 October 2026. **New exact local countermodels, not globally completed
compact polytopes.** Paper/algebra only; no scientific job, enumeration,
Lean, TeX or build was run. Existing sources are unchanged.

The eight-field B simple-end distance-two deletion expressly excludes
K7 before its final six-row argument. Here the missing finite label7
produces TWO genuine compact local quads. The first has a slack charge
recipient. The second has two whole-B source corners and requires a
separate simultaneous-consumer rule.

## 1. Exact field and old vertex

Use K=Q(cos(2pi/7)), k=cos(pi/7), and

    delta=4k^2-3,
    delta^3+4delta^2+3delta-1=0,
    0<delta<1/4,
    Norm(delta)=1,
    delta^-1=delta^2+4delta+3,
    (delta+2)^2=1+delta^-1,
    Norm(delta+2)=-1.                            (1)

The last norm follows by evaluating the cubic at-2. At both
nonphysical embeddings k_i^2<1/2 and delta_i<-1. These exact
field identities and angle comparisons are also derived in
`arithmetic_k7_an_far_bridge_quad_survivor_20261019.md`.

Let n=19 or20 and start with the entire old B_n path

    a-b-c-t_1-...-t_L, L=n-3,

whose only4 bond is at the remote terminal t_L. In unit roots
every early path product is-1/2 and the terminal4 is-sqrt2/2.
Take quad sides a,c,E,D and common whites

    W={b}+{t_1,...,t_L}=A1+B_(n-3).

In both models E has ONLY the old contact E.a=-k and E.D=0.
The native opposite product c.D=-y will be specified below.

## 2. Model I,I: a single source with a slack recipient

Let D touch ONLY c among the old roots, and put

    y=(delta+2)/2>1.

The corank relation is exactly

    (4k^2-3)*(4y^2-1)=delta*delta^-1=1.          (2)

Its complete corner diagrams are

|corner|entire finite diagram|baseline at the mark|
|---|---|---:|
|(a,c)|B_n|1, simple-end first/third pair|
|(c,E)|B_(n-1)+A1|0, second node versus isolate|
|(E,D)|B_(n-3)+3A1|0, two isolates|
|(D,a)|A2+B_(n-3)+A1|0, A2 root versus isolate|

Thus an unbalanced deletion of even half the old unit loses
the required face payment in this local model. At n=20 the
last recipient forest has baseline1+31=32 and has only TWO
possible A2-root/isolate marks. A full unit coupon at the
actual receiving role is therefore compatible with allowance35
for that forest even without a sharper geometric capacity.
This observation addresses this model only, not every consumer.

## 3. Model I,J: adjacent simultaneous sources

Now let D.b=D.t_1=-1/2, D.c=-y, with no other old contact,
and put

    2y=delta^-1-1>3.

The exact corank relation is

    delta*(2y+1)=1.                              (3)

The two corners(a,c),(D,a) are whole B_n simple-end
distance-two sources. The other two corners(c,E),(E,D)
are BOTH whole B_(n-1)+A1 with the marked isolate E and
the SECOND node from the simple end of the B component.
They have baseline zero. The old balanced-quarter recipient
on the simple ENDPOINT versus the isolate does not pay
these different marked roles.

At n=20 these two recipient forests have baseline35. A new
half-coupon there must therefore be accompanied by a funded
refund or a new global obstruction. The slack recipient of
Model I,I does not occur in this model.

The simple four-cycle b-c-t_1-D-b has the dotted c.D
diagonal. It is NOT an induced affine diagram. Dropping that
actual chord would manufacture a false compactness obstruction.

## 4. Integral norm-root realizations and full-support radicals

Use norm-two roots on the long B path and norm-one at its
short terminal. All old Gram entries are integers: the remote4
has product-1 with norms2 and1. Take E,D of norm2. The only
nonintegral-looking entries are

    a.E=-2k,
    c.D=-(delta+2) [I,I],
    c.D=-(delta^-1-1) [I,J],

which are all algebraic integers by(1). There is no assertion
that sqrt2 itself lies in K.

In these norm-root coordinates the literal full-support kernel
relations are respectively

    2a+2kE+(1-delta)b-2delta*(c+sum_i t_i)
      -delta*(delta+2)D=0,                       [I,I]

    2a+2kE+(1-delta)b-delta*(c+D)
      -2delta*sum_i t_i=0.                       [I,J]          (4)

Every row cancels directly using(1)--(3). For example the
Model I,I c row uses delta*(delta+2)^2=delta+1; the Model
I,J c row uses delta*(2y+1)=1. The terminal short-root row
has diagonal1 rather than2, which is included in(4).

The coefficients of D are units. Hence D belongs to the
integral module on the independent old B_n roots plus E.
All displayed norm-root reflections preserve that module.
The extra norm2 roots have the required integral pairings;
these are stronger data than a list of integral cyclic products.

The old B_n plus E has determinant ratio1-2k^2<0 physically
and1-2k_i^2>0 otherwise. Its exact inertia is(n,1) physically
and positive definite at both nonphysical embeddings. Relations
(4), with nonzero D coefficient, make the full Gram congruent
to that block plus zero. Thus both models have physical inertia
(n positive,1 negative,1 zero), and conjugate inertia
(n+1 positive,0 negative,1 zero).

All kernel coefficients remain nonzero at every embedding.
Every proper nonphysical principal block is consequently
positive definite; no proper principal affine singularity
can occur physically. These facts do NOT prove reflectivity
or finite-volume compact completion of the local wall system.

## 5. Actual compact projected quadrilaterals

After eliminating W the old a,c Gram is

    a_bar^2=3/4,c_bar^2=1/4,a_bar.c_bar=-1/4.

Model I,I has E_bar^2=D_bar^2=1, products a.E=-k,
c.D=-y, and other new projected products zero. Eliminating
E,D leaves the two-by-two matrix

    [-delta/4   -1/4      ]
    [-1/4       -1/(4delta)],

negative semidefinite of rank1. Thus the projected Gram
has inertia(2,1,1). Each cyclic adjacent pair is positive
definite: its determinants are1/8,1/4,1,3/4. Both opposite
pairs are indefinite, with determinants-delta/4 and
-1/(4delta).

For Model I,J, c,D have equal projected square1/4,
product-y-3/4, and products-1/4 with a. They are both
orthogonal to E, whose norm is1. Eliminate E and split
c,D into their symmetric and antisymmetric directions.
The antisymmetric direction is positive; on a,c+D the
remaining matrix has diagonals-delta/4 and-delta^-1,
off-diagonal-1/2, and determinant0. It has one negative
and one zero direction. Again the projected inertia is
(2,1,1). Its cyclic adjacent determinants are1/8,1/4,1/4,1/8,
and its opposite pairs are indefinite.

In both cases all projected products are nonpositive and the
contact graph is connected. Perron's negative eigenvector
gives a strict timelike point of all four halfplanes. Acute
positive-definite projection gives each of the four stated
finite vertices inside the other two halfplanes. Every side
has two finite endpoints, so the intersection is a compact
quadrilateral. The full H^n wall intersection may have infinite
volume; no compact n-polytope is constructed here.

## 6. Exact audit of the new finite label7 branches

In the whole-B simple-end proof, E=I and finite x give
q_E=2x^2, K=Q(x^2), with all other x_i^2<1/2. In K7
native label7 survives; label14 fails because a conjugate
cos^2(3pi/14)>1/2. Rational labels force K=Q. Thus x=k
is the only missing finite entry. The E=L finite test has
K=Q(x) and all other x_i<0; k has a positive conjugate,
so that branch has no new finite K7 case.

At x=k the old six-row identities give exactly:

    I,I: 4y^2=1+delta^-1=(delta+2)^2;
    I,T: 2(2y+1)^2=(delta+2)^2;
    I,J: 2y+1=delta^-1.

In I,T, the simple actual cycle puts2y in K7, so the middle
identity would put sqrt2 in K7, impossible. The other two
are exactly the models above. This is the relevant finite
entry audit, not a claim that all dotted or mixed consumers
of a future funded rule have already been checked.

An exact local arithmetic survivor does not refute a theorem
about globally compact chambers; global walls might exclude
it. It does refute an attempted proof that dismisses it by
the current corank, finite-corner, conjugate, integrality or
proper-affine tests. A correct K7 charge needs to replace
the face payment or introduce that additional global argument.

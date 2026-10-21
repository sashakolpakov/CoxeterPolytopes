# A prime-three golden extension: complete G2 edge and full diagonal facet

20 October 2026. Own paper proof; no mathematical program or search was
run. This note proves an actual edge and a uniform full-facet identity.
It does **not** claim reflectivity or nonreflectivity of the resulting
eight-dimensional diagonal facet, nor a new compact example.

Put `phi=(1+sqrt(5))/2`, `O=Z[phi]`, and

    L_n = <-phi> + I_n + A2,                    n>=7,
    A2 Gram C = [[2,-1],[-1,2]].

Coordinates are `(x; y_1,...,y_n; z_1,z_2)`. The form is Lorentzian at
phi>0 and positive definite at its conjugate. Its dimension is n+2.
The determinant ratio to `<-phi>+I_(n+2)` is 3, a nonsquare in K because
3 is inert and has odd valuation. The rational norm `N(3)=9` is square
and supplies NO such distinction. In even rank, scalar similarity
preserves this K-valued determinant squareclass. No classification of
other integral modules or commensurability classes is asserted.

## 1. Complete primitive norm pool and zero germ

The ring O is a PID (its ideal Minkowski bound is `sqrt(5)/2<2`), with
units `+-phi^j`; totally positive units are squares. Both 2 and 3 are
inert. For a primitive reflective vector v, the pairing ideal contains
`det(L_n)` by the adjugate identity and primitivity. Thus `v^2` divides
`2 det(L_n)`, and unit rescaling gives the complete necessary pool

    q in {1,2,3,6}.

For q=3 or 6, reflection integrality forces x and every y_i into 3O,
and forces `Cz` into `(3O)^2`. At time zero a nonzero such y_i contributes
at least 18 to the trace of the positive norm, exceeding `Tr(q)<=12`.
So these zero roots lie entirely in A2.

Write any spatial vector as `a+phi b`, with a,b in the corresponding
ordinary integer lattice. For a rational norm q its two equations are

    a^2+b^2=q,                  2(a,b)+b^2=0.       (1)

The ordinary A2 lattice is even, has minimum 2, and has no vector of
norm4. Equations(1) therefore give no A2 norm3 vector; norm6 vectors
are precisely the ordinary A2 norm6 vectors. Norm1 roots are coordinate
roots of I_n. At norm2, a nonzero A2 component already consumes the
entire trace4; equations(1) give an ordinary A2 root. In I_n the norm1
and norm2 equations likewise give the ordinary coordinate roots and
`+-e_i+-e_j`; the case `a^2=b^2=1` is excluded by the second equation.
There is no mixed zero root.

Consequently the FULL zero arrangement at `p=(1;0)` is `B_n+G2`.
For the last block choose actual outward simples

    rho=-a1,          sigma=a1-a2,
    rho^2=2,          sigma^2=6,          rho.sigma=-3.

The G2 sector is `z_1<=z_2<=2z_1`. The other simples are the standard
`-e_i+e_(i+1), -e_n`. Local finiteness of integral mirrors follows by
bounding their reflection matrices in both embeddings; it gives the
usual actual-first-wall interpretation of a finite edge interval.

For orientation only: the face retaining the entire G2 block has the
full `<-phi>+I_n` arrangement. A finite rank-three spherical extension
cannot connect to an I2(6) block, so every neighboring wall is
orthogonal to it. Conversely all intrinsic norm1/2 reflections extend;
orthogonal primitive norms3/6 are impossible by the divisibility above.

## 2. A finite upper and all finite continuation profiles

Release sigma, retaining rho and B_n. The ray has points

    P(tau)=(1;0^n;tau,2tau).

The integral primitive reflector

    c=(3phi;0^n;phi^2,2phi^2)

has norm3 and covector `(-3phi^2;0^n;0,3phi^2)`. It is orthogonal to
rho and B_n, acute to sigma, and meets this ray at tau=1/2. This is
finite because `P(1/2)^2=3/2-phi<0`.

Let w be the first wall, normalized to its canonical q. Write
`b=-w.rho>=0`. The finite products with rho, for q=1,2,3,6 respectively,
are exactly

    b in {0,1}, {0,1,phi}, {0}, {0,3}.           (2)

Indeed the possible native labels over K are 2,3,4,5,6,10. Substitution
in `b^2=2q cos^2(pi/m)` leaves precisely (2); nonintegral values and
the nonsquare classes 2,3,6 and phi+2 exclude the others.

The finite spherical classification permits an added root either to
leave B_n separate, or to extend it to B_(n+1) at its long end. The
latter root has q=2 and I_n coordinates `(1,0,...)`. It is either
orthogonal to rho or joins it by a 3-bond, producing B_(n+2).
No other connected finite extension of B_n exists for n>=7.
Thus all ten profiles are:

| q | I_n contribution epsilon | b | R=2q-2epsilon-b^2 |
|---|---|---|---|
|1|0|0|2|
|2|0|0|4|
|3|0|0|6|
|6|0|0|12|
|1|0|1|1|
|2|0|1|3|
|2|0|phi|3-phi|
|6|0|3|3|
|2|1|0|2|
|2|1|1|1|

Set `Z=2z_1-b=z_2`, so `z_1=(Z+b)/2`. The norm and intersection
conditions imply

    3Z^2-2phi*x^2=R,
    0<x, Z>0,                 phi*x/(3Z)<=1/2,
    x^2 <= (3phi^2/2)R,       x'^2 <= (phi/2)R'. (3)

The last bound is conjugate positive definiteness. All tests below
overcount by allowing the boundary and by omitting some integrality
and old-wall inequalities; hence they cannot lose a first wall.

## 3. Complete tiny interval census

Modulo3 the nonzero squares in F9 are
`{1,2,1+phi,2-phi}`. Equation(3) requires
`x^2=R(phi-1)` in F9. Therefore R=1,2,4 are impossible. For R=3,
`x=3u`; its bounds give `0<u<7/6`, `|u'|<8/15`.
Writing `u=a+d phi` forces d=0, then no positive integer u survives.
For `R=3-phi`, the safe bounds `x<7/3`, `|x'|<7/4` give exactly
`x in {1,phi-1,phi}`; none has the required residue `x^2=2 mod3`.

The remaining isolated q=3 and q=6 rows have x=3u by reflection
integrality. Direct coefficient strips in (3) give:

| row | exact bounds | complete positive u list |
|---|---|---|
|q3,R6|`u<=phi, u'^2<=phi/3`|`phi`|
|q6,R12|`u<=sqrt(2)phi, u'^2<=2phi/3`|`1,phi`|

For clarity, writing `u=a+d phi` in either row forces d=0 or1.
In the first row d=0 gives no integer; d=1 gives only u=phi because
u=phi-1 violates the conjugate bound. In the second d=0 gives u=1,
and d=1 again leaves only u=phi. Thus there is no unbounded Pell or
unit-orbit premise.

At q6, u=1 gives `Z^2=2phi^4`, impossible since 2 is nonsquare in K.
At u=phi it gives `Z^2=10+12phi`, whose rational norm is76, not a
square. At q3 the sole solution is `x=3phi, Z=2phi^2`, reconstructing
exactly c. It passes all integrality and actual-wall inequalities.
It is consequently the unique first wall of the full arrangement.
The endpoint has the full finite germ `B_n+A1(rho)+A1(c)`.

## 4. The actual c facet is a full unimodular diagonal chamber

The pairing ideal of c is exactly(3), and c^2=3. Hence orthogonal
projection along c is integral and

    L_n = M + <3>,                  M=c-perp intersect L_n.

The equation of M is `z_2=x`. Its binary Gram in x,z_1 is

    [[2-phi,-1],[-1,2]].

The vectors `(phi,phi)` and `(phi^3,phi^2)` are orthogonal, have norms
1 and -phi, and have determinant `-phi^2`, a unit. Therefore

    M is integrally <-phi>+I_(n+1).             (4)

This identifies the FULL induced arrangement, not just its lattice.
Write an ambient primitive root as `v+k c`, v in the unimodular M.
For its mirror to cross c, `v^2=q-3k^2` must be positive in both
embeddings. For q=1,2,3 a nonzero k contradicts its nonzero integral
field norm. For q=6 the two bounds `|k|,|k'|<sqrt(2)` leave only
k=+-1 (apart from0). Then v^2=3, whereas reflection integrality and
unimodularity force v in3M, impossible since its norm would lie in9O.
Thus every crossing is orthogonal. Orthogonal primitive ambient
norms3/6 are likewise impossible; every intrinsic norm1/2 mirror
extends to L_n. The actual c facet is exactly a full M chamber.

In particular, compactness of this literal L_n chamber would imply
compactness of the literal golden diagonal chamber of dimension n+1.
No commensurability or finite-index reflection transfer is substituted
for this full-facet statement.

## Scope and remaining published-input gate

For n=7 this gives an actual dimension8 facet with lattice
`<-phi>+I8`. A purported published cutoff for that literal unimodular
family would settle this route uniformly for n>=7. It is NOT imported
as an established exclusion here without an exact verified primary
statement. Allcock's accessible §2 explicitly gives the SILVER
`<-1-sqrt(2)>+I_n` cutoff n<=6; it must not be cited as the golden
cutoff. The golden n<=7 statements found in later tables and Vinberg's
1985 discussion point to Bugaenko1984, whose exact original theorem
text remains the external-input check. No job, new example, global
dimension bound or whole commensurability-class exclusion is claimed.

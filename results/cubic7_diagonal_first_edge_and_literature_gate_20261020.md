# Cubic-seven diagonal lattice: complete first edge, primary-source gate

20 October 2026; paper only, no scientific execution. Put

    f(T)=T³+T²−2T−1,    t=2cos(2π/7),    u=t²+t−1=2cos(π/7),
    K=Q(t),             O=Z[t],          L_n=<-t>+I_n.

The distinguished embedding has t>1; the other two are negative. Thus
L_n is admissible. This note proves one complete actual edge for n≥5
and records a decisive unresolved literature check. It does **not**
claim reflectivity, nonreflectivity, novelty, or a new H9 example.

## Arithmetic and the complete initial germ

The polynomial discriminant is49, and O is the full ring of integers:
an index7 would force field discriminant1, impossible by Minkowski.
The ideal-class bound is14/9<2, so O is a PID. The units −1,t,t+1
have independent real sign patterns. Since a real cubic field has
eight unit squareclasses, every totally positive unit is a square.
Modulo2, f is irreducible, so2 is inert. Primitivity and unimodularity
then make the complete primitive reflector norm alphabet `{1,2}`.
Every integral vector with one of these norms is reflective.

For nonzero z∈O, `Tr(z²)≥3` by its positive integral norm and
arithmetic–geometric mean; equality forces z=±1. At time zero a norm1
root is a coordinate vector. A norm2 root has at most two nonzero
coordinates. Two require both to be ±1; one would require √2∈K,
impossible in degree3. Thus the *entire* initial mirror arrangement
at p=(1;0^n) is Bn, with outward simple walls

    a_i=(0;−e_i+e_(i+1)), i<n,       a_n=(0;−e_n).

## One complete actual edge for every n≥5

Release a1, retaining B_(n−1). The explicit upper root

    r=(1;u,0,…,0),       r²=2,

meets this ray at `p_r=(u;t,0,…,0)`, whose norm is −2t. Also
`r.a1=−u`, hence the native angle isπ/7. Its other old products vanish.

The complete finite continuation inventory for retained rank≥4 consists
of an isolated root or a Bn extension at the long end. An isolated root
has spatial coordinates `(y,0,…)` and residual
`y²−t x²=q`, q=1 or2. The connected root has norm2 and spatial
coordinates `(y,1,0,…)`, with the same residual equation for q=1.
It must be included; assuming orthogonality would omit this profile.

Any wall no later than the displayed upper satisfies

    x>0,       x²≤q/2,       σ_i(x)²≤q/(−σ_i(t)), i=1,2.

Indeed its collision parameter is `t x/y≤t/u`, so y≥ux and
`q≥(u²−t)x²=2x²`. The other bounds are conjugate positivity.
Since `Σσ(1/t)=−2`, these imply

    Tr(x²)≤q(5/2+1/t)<7q/2.

For q=1 the integer trace is≤3, so x=1, contradicting x²≤1/2.
For q=2 it is≤6. Write x=a+bt+ct². The exact trace form is

    Tr(x²)=[(3a−b+5c)²+14(b²−bc+c²)]/3.

Consequently `b²−bc+c²` is0 or1. The zero case gives x=0,±1.
The six unit hexagon pairs `(b,c)` give precisely

    ±t, ±(t²−2), ±u, ±(t+1), ±(t²−1), ±(u−1).

This is a hand coefficient census, not a height search. Put v=t²−2.
The three real conjugates of t, in cyclic order, are t,v,−u. Physical
positivity and x≤1 immediately remove ±t,±u,±(t+1). Apart from x=1,
the only physical possibilities are −v,v+1,u−1. They fail one of
the conjugate inequalities respectively by the exact identities

    u t²=u+1>2,
    u(t+1)²>2,
    (−v)(t+1)²=t+1>2.

Thus x=1 is the only survivor, and y²=t+2=u² forces y=u.
Local finiteness of the integral arrangement supplies a first wall on
the compact timelike segment to the explicit upper. Every first wall
has appeared in the finite continuation inventory. Hence r is the
unique first simple wall of the full arrangement, and p_r is an actual
`B_(n−1)+A1` vertex. No finite compatible prefix has been mistaken for
a compact chamber.

## Lower-face information and the nonunit-glue boundary

The same trace census proves that an algebraic integer b all of whose
conjugates have absolute value<√2 is0 or±1. Therefore a norm2 mirror
crossing a norm1 coordinate facet has product0 or±1 with its normal;
its nonorthogonal projection is an integral norm1 root. Norm1 mirrors
cross such a facet only orthogonally. Conversely every intrinsic norm1
or2 root extends. Thus the actual coordinate facet carries the full
L_(n−1) arrangement, uniformly in n.

The actual pair (a1,r) has positive Gram determinant β=2−t,
`N(β)=7`. No finite irreducible rank3 diagram extends I2(7), so a
mirror crossing its orthogonal hyperbolic face is orthogonal to the
pair. The integral complement has basis

    (u;t,t,0,…,0),e3,…,en

and form `<-tβ>+I_(n−2)`. The sum of this complement and the I2(7)
lattice has determinant ideal index(β), not1. Its induced arrangement
contains exactly the complementary norm1/2 mirrors. One must not
replace it by the *full* reflection arrangement of that nonunimodular
complement: extra norm classes β and2β can require a separate glue
argument. No such unproved full-complement descent is used here.

## Repository and primary-literature screen

The existing [cubic E8 H9 scope audit](fresh_cubic7_h9_independent_scope_audit_20261013.md)
and [dihedral audit](fresh_cubic7_dihedral_independent_audit_20261013.md)
explicitly certify bounded prefixes only. Their escapes do not exclude
the full reflection group or decide this literal diagonal lattice.

The bibliographic red flag is Bugaenko, *Arithmetic crystallographic
groups generated by reflections, and reflective hyperbolic lattices*,
Advances in Soviet Mathematics8 (1992), pp33–56,
[AMS DOI10.1090/advsov/008/02](https://doi.org/10.1090/advsov/008/02).
Bogachev's author-maintained bibliography attributes a classification
of unimodular cubic-seven reflective lattices to that chapter and gives
the phrase “reflective for all n≤4”. The inspected
[primary Bogachev paper, Section5.1](https://nvbogachev.netlify.app/Bogachev-GeomArithmCoxeter-2020.pdf)
confirms that Bugaenko classified the unimodular case, but does not
state the exact dimension cutoff and quantifiers needed here.

The original AMS chapter was not retrievable in this bounded check:
its landing/PDF endpoints returned403 or inaccessible. The related
Bielefeld1992 booklet exceeded the web extractor's size limit, and
the Guglielmetti thesis endpoint timed out. A Springer DOI sometimes
attached to the citation, `10.1007/BF00883263`, is an unrelated
mechanics article and is not evidence for this theorem.

**Stop condition:** park this route before any further edge effort or
worker until the exact primary theorem is read. A bibliography phrase
does not establish that all n≥5, all unimodular lattices, or this
specific literal form are excluded. If that scope is confirmed, the
route is already published dead and should not be reopened. The
first-edge calculation above remains a small independent calibration,
not a new exclusion claim.

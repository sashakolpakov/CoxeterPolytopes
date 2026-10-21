# Exact local obstructions for the two remaining initial-edge pencils

19 October2026. New root written argument, awaiting independent review.
No enumeration or proof program was run to derive it. This is a constraint
on specified root pencils and one finite-vertex continuation, not a proof
that the entire H19 lattice is nonreflective.

Use the reviewed H19 lattice G=[[-2phi,phi],[phi,2]]+E8 over O=Z[phi],
phi²=phi+1, D=1+5phi=phi(4+phi), b=4+phi. Its canonical primitive
reflection norms are2 and2b. The exact26-wall check leaves initial
spacelike edges in directions omega2 and omega3, where

    omega_i=C_E8^{-1}e_i, omega2²=30, omega3²=20.

The pure pencil orthogonal to all eight fixed edge walls is

    r=(2u,-phi*u,t*omega_i),    u,t in O.

Indeed the old0 product forces the displayed binary relation; phi is a
unit, so u is integral. The seven fixed E8-simple products force the
spatial vector to be t*omega_i, and omega_i is primitive integral,
so t is integral. Write k=omega_i²/2. Then

    r²=2*(k*t²-D*u²).                               (1)

## 1. The omega2 pencil contains no lattice-reflection direction

First take canonical norm2. Equation(1) becomes15t²-Du²=1. Modulo3,
O/3O is F9, because X²-X-1 has no root in F3. The class of D is
1+2phi and its inverse is2phi; hence

    u²=-D^{-1}=phi in F9.

But the complete set of squares in F9 is

    {0,1,2,1+phi,2+2phi}.

Thus there is no solution, independent of height or coefficient size.
Equivalently the F9/F3 norm of phi is2, which is not a square in F3.

For canonical norm2b, reflection integrality gives t=b*tau: the spatial
covector is t*e2 and must lie in b times the dual lattice. Divide(1)
by b to obtain15b*tau²-phi*u²=1. Reducing modulo3 gives

    u²=-phi^{-1}=1+2phi,

which is also absent from the same square set. Thus neither canonical
norm class occurs. Rescaling a primitive reflector by an allowed unit
does not evade this argument; its canonical representative would still
lie in the same pencil.

## 2. The omega3 long pencil is also obstructed

For norm2b and k=10, the same divisibility gives

    10b*tau²-phi*u²=1.

At the ramified prime over5, phi reduces to3 in F5. The equation gives
-3u²=1, hence u²=3. The square classes in F5 are0,1,4, a contradiction.
The short equation10t²-Du²=1 is not ruled out by this calculation and
remains a legitimate bounded-search target.

## 3. A long root cannot close the omega2 edge through a finite corner

There is no finite nonright angle between a canonical long root of norm2b
and a short root of norm2. If m were its Coxeter label and g their
product, then g²=b*q_m with q_m=4cos²(pi/m) in K=Q(sqrt5).
The elementary cyclotomic-degree list gives m=3,4,5,6,10; labels8,12
have quadratic fields sqrt2,sqrt3 instead. The rational norms are

|m|q_m|N(b*q_m)|
|---|---|---:|
|3|1|19|
|4|2|76|
|5|phi²|19|
|6|3|171|
|10|phi+2|95|

None is a rational square, whereas N(g²)=N(g)² must be one. The
nonright angle is impossible. At a finite corner with eight fixed short
walls every product of a new long wall with those walls would therefore
have to be0. For the omega2 edge this is precisely the pencil excluded
in Section1.

Consequently every possible new mirror that closes that edge at a
finite vertex must be SHORT and contact at least one of its eight fixed
walls nonorthogonally. For short/short canonical roots the only finite
nonright products in K are-1 and-phi (labels3 and5); labels4,6,10 would
require a product outside K. The finite continuation problem has thus
been narrowed to those explicit short contact profiles.

This does not exclude such profiles. Nor does it establish a complete
mixed-norm Vinberg prefix, an infinite-order chamber symmetry, or a
whole-lattice nonreflectivity theorem. A finite-sided finite-volume
completion would be compact because the conjugate lattice form is
positive definite and forbids a K-rational null cusp. The unexcluded
short-root finite continuations remain essential work.

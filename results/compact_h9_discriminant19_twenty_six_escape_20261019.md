# An exact initial-edge escape in the26-wall discriminant19 seed

19 October 2026. **New written algebraic obstruction, pending independent
review.** No scientific computation, checker body, root search, optimizer,
Lean process or build ran. The26-wall checker was only prepared and syntax
parsed separately. Every earlier constructive packet is unchanged.

The26-wall seed defined in
`compact_h9_discriminant19_nonpolar_pell_mirror_20261019.md` has an exact
spacelike feasible vector. Therefore its volume is still infinite. Its
26 actual facets, two full compact H8 facets and at least84 finite
vertices are not invalidated. No global nonreflectivity or whole
completion obstruction is proved.

## Vector and norm

Use the same native E8 indexing, phi^2=phi+1 and D=1+5phi. Put

    omega2=(10,20,30,24,18,12,6,15),
    C omega2=e2,          omega2^2=30,
    A=13+8phi,
    Y=(66,-33phi,A*omega2).

The displayed vector verifies C omega2=e2 directly. The binary part of
Y is33 times the old timelike control p=(2,-phi,0). Since
A^2=233+272phi, its square is

    Y^2=-2178D+30A^2=4812-2730phi>0.

The strict inequality follows already from phi<5/3, giving lower
bound262. Also (p,Y)=-66D<0, so the relevant light-cone crossing is
in the future component selected by p.

## Feasibility against all26 actual walls

For a positive-control root h=(x_h,y_h,z_h), the exact product is

    (Y,h)=A*z_h,2-33D*x_h.

The subscript2 denotes the THIRD native E8 coordinate, zero based,
not either binary coordinate. Reading the literal prefix coordinates
and the displayed later root coordinates gives the following complete
sixteen-plus-one positive-control ledger:

|Wall|Product with Y|
|---|---:|
|old9|15-39phi|
|old11|-85+12phi|
|old12|0|
|old13|28-31phi|
|old14|-72+20phi|
|old15|15-39phi|
|old16|-41+23phi|
|old17|5-13phi|
|old18,old19|-8-21phi|
|r|-33D|
|k|-87-150phi|
|w|-126-174phi|
|a|30-78phi|
|d|-152-190phi|
|z|8-188phi|
|s|-110-132phi|

All nonzero entries are strictly negative, using3/2<phi<5/3 where
needed. Among the nine zero-control old walls, old3=-e2 has product
-A<0; old0 and the other seven negative E8 simples have product0.
This accounts for every one of the26 inequalities.

For an independent view of the equality, the projectively rescaled
point is on the fixed direction from the initial vertex:

    Y/33=p+tau*(0,0,omega2),
    tau=(13+8phi)/33
       =D*(1+phi)/(6+15phi).

The latter denominator and numerator are the third E8 coordinate
and D times the control coordinate of old12. The complete ledger
shows that old12 is its first other-wall intersection. All other
positive-control walls are strict there. No extremality or complete
cone enumeration is required for the volume argument.

## Infinite volume and scope

The preserved center c is strict timelike for all26 walls by the
reviewed construction. A sufficiently small positive perturbation of
Y towards c stays spacelike and is strict against every wall. The
segment from c to this point crosses the future light cone with all
26 inequalities strict. It exposes an open ideal patch and therefore
proves infinite hyperbolic volume.

The nonpolar norm-two root s cuts the previous escape X; it does not
cut this new Y. Thus solving that exact next-mirror problem did not
settle the whole cone. A different integral wall or another construction
could still change this chamber. No polar-line classification, complete
mixed-norm Vinberg search, K3 marking, nonreflectivity or dimension
record is claimed here.

The native prefix and all constructive source pins are those in
`h19_twenty_six_wall_seed_checker_preparation_20261019.md`. The prepared
checker's nine initial-p direction diagnostics recompute this kind of
first-hit norm test, but no execution or saved replay of that source
is asserted by the present written proof.

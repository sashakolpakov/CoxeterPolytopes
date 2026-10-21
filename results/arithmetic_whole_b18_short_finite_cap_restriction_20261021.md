# Whole B18: the short-edge finite-own cap has one possible type

21 October2026. New source/manual supporting lemma, pending full root and
outside-author review. No scientific, symbolic, numerical, enumeration,
Lean or TeX job ran. Frozen sources, receipts and bounds are unchanged.
The author wrote the inherited D19 finite-angle argument below; this is
new B18 arithmetic and not an independent review of that earlier input.

**Proposed theorem.** At a vertex with WHOLE diagram B18 of a compact
arithmetic first-type hyperbolic18-polytope, release the short node18.
If its actual opposite has a positive FINITE own contact, its only
possible profile is own3 and one white3 contact at node17. The ACTUAL
ground field is then Q(sqrt2). In the rational old-B normalization its
new-root norm is2 and its squared time component is24+17sqrt2.

This is a necessary literal-cap restriction, valid without a weak-height
hypothesis or a farthest-vertex assumption. It does not exclude the
remaining profile, dotted short-edge caps, whole B18 or dimension18.

## Inputs and actual root normalization

| Input | SHA256 | Used interface |
| --- | --- | --- |
| results/arithmetic_uniform_classical_weak_helper_obstruction_20261021.md | 535ac0d0fa92d28ccab2976fa147a4ccff0bd3bfcd1e04ca4e134c6ae6857905 | Field-valued B long-square2/short-square1 normalization, cyclic integrality, strict conjugate positivity |
| results/arithmetic_d19_nonisolated_finite_fork_exclusion_20261021.md | 6b824daa86d75c97a04e3f429bdce56f3d01c93a7d79d40b8cd0b8510b2006da | Section4's elementary finite-angle degree<=2 calculation only |
| manuscripts/compact_arithmetic_dimension_bound.tex | 1f09df2976085fd1ad1aa56ed35f1ab182ff95a8332dcfa0b0b80d4849b82203 | Lines1050–1052: rational first-type compactness obstruction via Meyer |

External inputs are compact Coxeter simplicity, finite spherical
classification, and admissible first-type root lines/cyclic integrality.
Meyer's theorem is used in its usual general form: a rational indefinite
form in at least five variables is isotropic, so its simplest-type
arithmetic quotient cannot be compact. Here the form has19 variables.
No larger-rank whole-diagram theorem, field list, prefix estimate,
rank19 quadraticity, host nonreflectivity or lattice classification is used.

Write the old field-valued roots as

    r_i=e_i-e_(i+1) (1<=i<=17), r18=e18.

Their Gram is rational, with long square2, short square1 and consecutive
products-1. Its raw duals are Omega_j=P_j, including j18, and

    J_jj=j, J_18,j=j, J_18,18=18.

The endpoint opposite c and the17 retained old roots form a positive
definite rank18 white Gram: the endpoint is finite by compactness.
The complete one-node finite extensions of A17 are therefore exactly:
no contact; one3 at1,2,16,17; or one4 at1,17. A connected finite
completion of rank18 can only be A18, B18 or D18; two contacts would
create a cycle. There are no other components or foreign contacts.

Put gamma=cos(pi/m), m>=3, at the own short node. For a white3,
rescale c OVER K to square2 and helper product-1. Its short product
is-t/2, with physically t=2sqrt2*gamma. Traverse the cycle through
the old B chain: its directed Cartan product is plus or minus t.
Thus t belongs to O_K without adjoining sqrt2. For white4 instead
normalize c OVER K to square1 and helper product-1; then its short
product is-t/2 with t=2gamma in K, integral by the native finite angle.
This norm-one scaling is valid because the original helper product a
satisfies a²=c² when the helper long root has square2.

In each case the full19-root Gram is Lorentzian and has entries in Q(t).
The new root is independent of the old18 because its time component is
nonzero. Any nonphysical embedding of a proper K/Q(t) fixing Q(t)
would preserve this Lorentzian Gram, contradicting the positive conjugate
ambient form. Hence ACTUAL K=Q(t). Rational t is excluded by Meyer.
The same argument applies with B in place of t in the no-white case.

## No-white finite own contact

Put B=4gamma², a nonzero algebraic integer in K. Scaling the new root
over K gives square B and short product-B/2, so the full Gram is over
Q(B). At every nonphysical embedding its Schur complement gives

    0<B_sigma<4/18=2/9.

Physically0<B<4. A minimal degree d>=2 would give
0<|Norm(B)|<4*(2/9)^(d-1)<=8/9, impossible. Rational B forces
K=Q and is excluded by Meyer. Thus every isolated finite-own short
opposite is impossible, including all nonnative finite labels.

## The six nonisolated roles and exact conjugate strips

If N is the normalized new-root square, its squared time component is

    H²=(9/2)t²+j*t+j-N,

where j is the helper position. Physically H²>0; every other conjugate
has H_sigma²<0. Set

    a=(sqrt19-1)/9 <3/8,
    b=(sqrt19+1)/9 <3/5.

The exact intervals are:

| White role | N | Required interval for every other t_sigma |
| --- | ---: | --- |
| 3 at1 | 2 | (-b,a) |
| 3 at2 | 2 | (-4/9,0) |
| 3 at16 | 2 | (-2,-14/9) |
| 3 at17 | 2 | (-2-a,-2+b) |
| 4 at1 | 1 | (-2/9,0) |
| 4 at17 | 1 | (-2,-16/9) |

For white4, physically1<=t<2. At1 the integral norm of t is
<2*(2/9)^(d-1)<1 for d>=2. At17 use z=t+2: physically3<=z<4
and0<z_sigma<2/9, giving the same contradiction. Thus both roles
are rational and impossible by Meyer.

For white3 at2, physically sqrt2<=t<2sqrt2<3; d>=3 would give
|Norm(t)|<3*(4/9)²=16/27. At16 instead use z=t+2<5 and
0<z_sigma<4/9, giving Norm(z)<5*(4/9)²=80/81 for d>=3.
Thus these roles require degree<=2.

For white3 at1, d>=4 gives |Norm(t)|<3*(3/5)³<1. If d3,
its norm is a nonzero integer of absolute value<27/25<2. If either
other conjugate is positive, its absolute value is <3/8, making
|Norm(t)|<3*(3/8)*(3/5)=27/40<1. Hence both are negative and
Norm(t)=1. Then t>25/9, and their sum is>-6/5, so Trace(t)>71/45>1.
AM–GM gives the sum of their magnitudes >=2/sqrt(t)>1 because t<3,
so Trace(t)<t-1<2. This contradicts its integral trace. Degree<=2
therefore holds here too.

For white3 at17 use z=t+2. Physically z<5 and every other z_sigma
lies in(-a,b). A degree d>=5 gives |Norm(z)|<5*(3/5)^4<1.
In degree4 the norm has absolute value<27/25<2. If any of the
three other conjugates is negative, its absolute value is<3/8,
giving |Norm(z)|<5*(3/8)*(3/5)²=27/40<1. Otherwise Norm(z)=1
and all three are positive. Then z>125/27 and their sum>3/2 by
AM–GM (z<5<8). Thus6<Trace(z)<5+9/5<7, again impossible.
This leaves degree<=3 at the near end.

## Complete finite-angle reduction

For white3 put B=t²/2=4cos²(pi/m)=2+2cos(2pi/m). If t had
degree3, Q(B) would have degree1 or3. Degree1 gives m3,4,6 and
t=sqrt2,2,sqrt6, none cubic. Degree3 gives exactly m7,9,14,18,
by phi(m)/2=3. Their respective integral norms Norm(B) are1,1,7,3.
But in Q(t)=Q(B), Norm(B)=Norm(t)²/8, whose2-adic valuation
cannot be zero. Therefore finite-own cubic t is impossible.

At degree<=2 the complete finite list is

    m3: t=sqrt2; m4: t=2; m6: t=sqrt6;
    m12: t=1+sqrt3.

For completeness, degree(B)<=2 first gives m3,4,5,6,8,10,12.
The extra t² values for m5,8,10 are3+sqrt5,4+2sqrt2,5+sqrt5.
The latter two have nonsquare rational norms8,20; the former would
require a rational square1/2 or5/2. Thus each extra t has degree4.

Compare the three irrational conjugates -sqrt2,-sqrt6,1-sqrt3
with the table. None belongs to the strips at1,2 or16. At17,
only-sqrt2 lies in(-2-a,-2+b); -sqrt6<-12/5<-2-a, whereas
1-sqrt3>-1>-2+b. Rational t2 is the already excluded rational
ground-field case. Consequently only own3/white3-at17 remains.

## Exact surviving type and the remaining obstruction

For this type t=sqrt2 and c has square2. Its spatial projection is

    V=Omega17+(sqrt2/2)Omega18,
    H²=V²-2=24+17sqrt2.

The nonphysical squared time is24-17sqrt2<0, so this type is not
eliminated by conjugate positivity. Its actual field is Q(sqrt2),
and its time norm is-2. Moreover

    24+17sqrt2 = sqrt2*(3+2sqrt2)².

Thus, after a K-rescaling of the time coordinate, the actual ambient
quadratic form is proportional to <-sqrt2>+<1>^18. This identifies
a rational-form class only. No compatible integral lattice has been
identified, so the preserved nonreflectivity proof for a literal
ramified-silver lattice cannot be transferred without an additional
invariant-lattice or actual-mirror argument. It is not a proof that
this surviving cap occurs in a compact chamber.

All dotted own contacts remain outside this theorem. The restriction
above is a fresh B18 finite-own inventory, not a D19 quadratic-field
conclusion, a global B18 exclusion or a new arithmetic bound17.

# Dimension nineteen: the actual field is quadratic

21 October 2026. New manual supporting theorem by root, pending full
outside-author review. No scientific, symbolic, enumeration, numerical,
TeX or Lean job was used. This belongs to the existing arithmetic paper
family. It does not exclude dimension nineteen.

## Statement and inputs

Every hypothetical compact arithmetic first-type Coxeter polytope of
dimension19 has an ACTUAL quadratic ground field. No prescribed list
of fields is assumed. More locally, this conclusion holds whenever
one of its farthest vertices has whole diagram D19.

The first statement follows from the accepted whole-D19 frontier.
The new step excludes the cubic fields without any hypothesis about
an outgoing fork wall. In particular, it does not apply the old
isolated-fork conclusion to a nonisolated wall.

All paths in the input table are under results/.

| Source | SHA256 | Interface used |
|---|---|---|
| arithmetic_dimension19_only_whole_d19_frontier_20261021.md | bbced9c52d498e371db2220704353e6ed89f77ca92b9d4cfbb0138390ffdb8ad | Every farthest vertex has whole D19 |
| arithmetic_d19_isolated_dotted_fork_continuation_20261021.md | 59e6c68bd60bce749ac7da60e0ba80ad82f91e0aa431dbfb297369125871fa5e | Section4 COMPLETE fresh .48 prefix inventory; first part of Section5 cubic prefix analysis BEFORE its original-fork interpolation; scalar E16(.48)>1.009327 |
| arithmetic_d19_centre_complete_actual_profile_floor_20261021.md | 3d88932839dbd25efa3d2f4ce0d433d0b779ca15a79721021579f6669abcdb14 | Complete actual centre17 floor, with arbitrary white contacts |
| compact_farthest_vertex_mixed_replacement_inequalities_20261020.md | 3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84 | Actual future feasible-point and farthest comparison |

These input interfaces have full internal reviews in the predecessor
acceptance ef3f26e2. The cited59e source also proves a B-bound and a
quadratic conclusion for an ISOLATED outgoing fork. Neither conclusion
is used here. Its prefix analysis has no isolated-fork premise: this
distinction is made explicit below rather than inferred from its title.

External boundaries are compact acute simplicity, finite spherical
classification, first-type K-valued normalizations/cyclic integrality
and positivity at nonphysical embeddings, together with Meyer's
rational isotropy theorem. No new external classification is imported.

Use unit old D19 roots, with omega_i=sqrt2*P_i for i<=17,
P_i=e1+...+e_i. Let x be the ONE original farthest vector,
||x||^2<1, x.n_j>0. Write S_i=x.P_i and h_i=sqrt2*S_i.

## 1. Degree at most three, without a fork anchor

Call i in5..16 weak when h_i<=sqrt(2i)-12/25. At least one index
is weak: otherwise the ordinary positive prefix reference through
5..16 has energy E16(.48)>1.009327>1 on this same x.

The complete changed-threshold inventory in59e Section4 treats
deleting i as A_(i-1)+D_(19-i), including D4 at15 and D3 at16.
It retains all finite-white endpoint/D/E/fork/end4/end5 roles,
both full bridges and isolated caps. Its arguments use only this
actual prefix opposite and the old D19 basis; none uses an original
outgoing fork or its parameter B. At every weak direction that proof
gives the following necessary alternatives:

* Nonisolated finite native/helper cases lock the actual field to
  degree at most2, or are excluded/strong by their affine, conjugate
  or actual H4 comparisons.
* A weak near-left dotted helper has2<t<2.3 and all other
  conjugates in(-1.3,0). The integer t(t+1) has physical absolute
  value<7.59 and other absolute values<.39. Hence degree>=4
  would give a nonzero integral norm<7.59*.39^3<1.
* A weak far simple3 helper has2<t<3, other conjugates in(-.7,.3).
  The same polynomial has physical value<12 and other absolute
  values<.39, again excluding degree>=4. D-position and right
  fork/end3 dotted cases already have degree<=2 by their stronger
  strips. The remaining terminal and bridge roles are strong.
* An isolated finite weak cap has degree<=2. For an isolated
  dotted weak cap, alpha=4gamma^2 is integral, K=Q(alpha),
  4<alpha<16 and0<alpha_sigma<2/i<=2/5. Degree>=5 gives
  |Norm(alpha)|<16*(2/5)^4<1. In degree4 its norm must be1,
  so every nonphysical conjugate exceeds25/64. Then
  |Norm(1-2alpha)|<31*(7/32)^3<1, also impossible.

Here the field locks refer to the ACTUAL field: the full old19-plus-cap
Gram has rank20 and is physically Lorentzian. An extra embedding fixing
its coefficient field would preserve that Lorentzian Gram at a
nonphysical place. Thus no unspecified field extension evades the
minimal-field norm arguments. These are precisely the independently
usable conclusions of the fresh inventory, not a restriction to
previously named quadratic fields.

Since some weak direction exists, [K:Q]<=3.

## 2. What a cubic weak direction can be

Assume [K:Q]=3. The two possible cubic nonisolated helper roles
are also impossible, as follows.

For a near-left helper let the other conjugates be -a,-b,
0<a,b<1.3. If either is at most1, the norm of t(t+1) has absolute
value<7.59*.25*.39<1. Thus both exceed1. The integer trace
t-a-b is in(-.6,.3), hence0, and the integer norm tab is3.
Consequently a+b=t, ab=3/t and t^3>=12. The integral second
coefficient 3/t-t^2 is strictly between-4 and-3, impossible.
Indeed its lower bound at t<2.3 is greater than-9167/2300>-4,
and its upper bound is at most-3t^2/4<-3.

For a far helper, a positive other conjugate would give
|Norm(t)|<3*.3*.7<1. Both others are therefore negative, and
|Norm(t(t+1))|<12/16<1. All remaining nonisolated weak roles
already force degree<=2. Every cubic weak cap is consequently
ISOLATED DOTTED at its prefix index, regardless of the fork opposites.

Write its conjugates as alpha=x in(4,16), a,b in(0,2/5), and
its minimal polynomial as X^3-TX^2+SX-N. The norm N is1 or2;
the trace ranges are7<=T<=16 if N1,13<=T<=16 if N2. The bounds

    2sqrt(N(T-4/5)) < S < (2/5)T+(5/2)N-4/25

follow from ab=N/x and evaluation of the polynomial at2/5.
Checking the finite integer intervals by squaring gives exactly:

| N | T | S | Discriminant |
|---:|---:|---:|---:|
|1|7|5|-44|
|1|12|7|257|
|1|13|7|-268|
|1|15|8|985|
|1|16|8|229|
|2|13|10|-104|

This is the full coefficient computation in59e before its fork
interpolation. The three negative discriminants are impossible for
three real roots. Thus only(T,S,N)=(12,7,1),(15,8,1),(16,8,1)
remain. In particular EVERY cubic weak alpha has x>T-4/5>=56/5.

For the first two polynomials their value at1/3 is1/27>0.
Their largest root exceeds4, so1/3 lies between their small roots.
The larger small root is>1/3; its Schur bound<2/i implies i<6.
For the last polynomial its value at2/7 is1/343>0, so the larger
small root is>2/7 and i<7. It follows, WITHOUT a fork anchor, that

    every cubic weak index is5 or6;
    S_i>sqrt(i)-.48/sqrt2 for every i7..16.          (1)

No assertion that every candidate polynomial is realized is needed.

## 3. The weaker floor at five

If5 is not weak, its .48 floor implies the .54 floor. If it is
weak, the actual isolated dotted wall there has alpha>56/5.
Put u=4/alpha<5/14. Its unit normal is c=A p-gamma*omega5,
gamma>1. The positive combination xi=c+gamma^-1*n5 is an ACTUAL
future feasible timelike vector with xi^2=-(1-u). It is orthogonal
to c and nonpositive against every actual wall. Farthest comparison
therefore gives

    h5-u*x.n5 >= sqrt(10-u)-sqrt(u(1-u)).

Since x.n5>0, its height loss is strictly less than

    u/(sqrt10+sqrt(10-u))+sqrt(u(1-u)).

The reciprocal is<.163 for0<u<1, because sqrt10>3.16 and
.163*(3.16+3)>1. Also u<5/14<1/2, whence
sqrt(u(1-u))<3sqrt5/14<.48, the last comparison following by
squaring positive quantities. Thus the loss is less than

    .163*(5/14)+.48 < .54.

In BOTH cases we have the strict floor

    S5>sqrt5-.54/sqrt2.                             (2)

The complete centre theorem independently supplies

    S17>sqrt17-.525/sqrt2.                          (3)

All floors(1)--(3) refer to this same original x.

## 4. One concave reference excludes the cubic case

Let f_i=sqrt(i)-.48/sqrt2, delta=.06/sqrt2, and
d=sqrt17-4-.045/sqrt2. Define w by

    w1=...=w5=(f5-delta)/5,
    w6=w7=(sqrt7-sqrt5+delta)/2,
    wj=sqrt(j)-sqrt(j-1), j8..16,
    w17=d, w18=w19=0.

This reference uses knots5,7,8,...,16,17. It omits6 completely.
Its coordinates are nonnegative and nonincreasing. The first slope
is>1/4 while the averaged slope is<1/4. The averaged slope is at
least the ordinary7th slope, hence exceeds the8th slope. The usual
square-root increments decrease. Finally .09<d<.1, whereas
4-sqrt15>1/8. Each comparison follows by squaring positive rational
bounds, such as4.123<sqrt17<4.124 and.031<.045/sqrt2<.032.

The dual coefficient w6-w7 is zero. Every positive dual coefficient
occurs at a floor in(1)--(3), where w.P_i equals that floor. Therefore
summation by parts and Cauchy give

    x.w>||w||^2,            ||x||^2>||w||^2.

Write A=2f5/5-(sqrt7-sqrt5) and
g=2sqrt6-sqrt5-sqrt7. Directly averaging the6th and7th slopes
of the full .48 reference, with its changed initial value, yields
the EXACT identity

    ||w||^2=E16(.48)-A*delta+(7/10)delta^2-g^2/2+d^2. (4)

Here A<.35, delta<.043,0<g<.02 and d>.09.
For the first bound use sqrt5<2.2361,sqrt7>2.6457 and
1/sqrt2>.7071; for g use sqrt6<2.45,sqrt5>2.236,
sqrt7>2.645. Each rational bound is verified by its positive
square. Strict concavity gives g>0. Omitting only the positive
delta^2 term in(4) gives

    ||x||^2>||w||^2
       >1.009327-.35*.043-.0002+.0081
        =1.002177>1.

This contradicts ||x||^2<1 and excludes every cubic actual field.
No separately optimized vectors or overlapping norm gains are added.

## 5. Conclusion and limits

The actual degree is at most2. Degree1 would give an indefinite
rational first-type form of rank20; Meyer's theorem makes it
isotropic, incompatible with a compact arithmetic quotient.
Consequently the actual field is quadratic.

This is a reduction over ALL actual fields, not a theorem restricted
to a proposed field list. It does not identify which quadratic fields
remain, exclude all nonisolated dotted forks, or prove ceiling18.
It uses no new computed proposition or script output. The frozen
six manuscripts and their selected Lean kernels are unchanged.

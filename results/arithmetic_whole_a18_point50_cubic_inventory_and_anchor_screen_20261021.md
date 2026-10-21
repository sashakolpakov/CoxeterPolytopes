# Whole A18: the fresh half-loss cubic inventory and its anchor comparison

21 October2026. Root supporting proof and finite-algorithm specification,
pending outside-author review and bounded execution of the new domain.
No whole-A18 or global dimension exclusion is asserted here.

## Inputs and purpose

Use the reviewed actual-field theorem
`results/arithmetic_whole_a18_point52_actual_field_degree_reduction_20261021.md`
(`b57aa75adbde5c1a1b6f78354c58bba47d03b2047a923438b133be9aee5f2707`)
and the reviewed isolated-cap algorithm specification
`results/arithmetic_whole_a18_isolated_cap_domain_and_norm_screen_20261021.md`
(`06172cc3d1f1e218a29d7e20916f105366d10344f4e086a4aeb4d4eed45e75ac`).
The latter's bounded execution is
`results/a18_isolated_cap_domain_norm_screen_execution_20261021.json`;
its complete saved output is
`results/a18_isolated_cap_domain_norm_screen_20261021.json`
(`671550c23549a9f55af602cd5f26578c5901e16646deecd540b328d8010386fc`).

That execution visits420 cubic and63 quadratic rows. It retains six
cubic polynomials and30 cubic index cells. All fifteen cubic pair groups
have exactly one polynomial and one k=i(19-i). The proved reflected
Cartan-integrality obstruction therefore shows that an actual cubic
whole-A18 case has at most one .52-weak central index. At least one
exists byb57. These are necessary-domain conclusions, not realizability.

Here we freshly classify the larger set of .50-weak central indices
under the hypothesis that ACTUAL K is cubic. The old .52 list is not
silently reused at a changed threshold. Unit old normals, duals omega_i,
heights h_i, D_i and C_ij retain the N19 definitions ofb57, i=5..14.
An index is .50-weak when h_i<=sqrt(D_i)-1/2.

## 1. The new edge threshold and complete white inventory

Writing q=sqrt(D_i), the sufficient edge-ratio threshold at loss1/2 is

    T(q)=2-(q-1/2)/(q(2q-1/2)).

It increases on the present range sqrt(140/19)<=q<31/10. Thus

    T(q)<3274/1767<1853/1000,

the last cleared difference being251. A weak index therefore has
R=sigma/sqrt(delta)<1853/1000. Only this necessary implication is used.

The complete finite white extensions are the same geometric TYPES as
inb57: isolated; one3 at a path end, D-position or E-position; one4
at an end; or two3 contacts at the ends of the two retained paths.
This follows from finite spherical classification with both white
path ranks>=4 and connected join rank18>8, independently of the
loss. There is no foreign component at a whole-A18 vertex.
Own0/parallel and the own3/6 affine configurations retain their actual
geometric exclusions fromb57; none depends on the chosen loss.

For a native3 helper at side position j, with m=i or19-i, use

    delta=1-j(m-j)/(2m), sigma=gamma+j/(2m), t=2gamma.

Normalizing overK to square2 gives integral t and ACTUAL K=Q(t).
All nonphysical conjugates satisfy the full old-plus-new Schur test.

## 2. Far endpoint3: no cubic parameter at the new loss

Here j=1. The identity R²>=4gamma-2 gives

    phi<t<68/25,

since ((1853/1000)²+2)/2<68/25. The strict lower bound phi may be
replaced by t>=phi; only t>202/125 is needed below.
The complete N19 inverse calculation inb57 gives each nonphysical
conjugate in(-a,b), a=193/250, b=12/25.

There is an additional useful bound for the SAME actual helper.
Let -a_i,b_i be the two roots of its Schur quadratic
D_i y²+2C_ie*y+D_e-4, where e is1 or18. Since D_e=36/19,

    a_i*b_i=(4-D_e)/D_i<=2/7.                       (1)

This couples the interval endpoints; independently maximized a and b
must not replace it in a mixed-sign norm calculation.

Suppose t has degree3. If both other conjugates are positive, its
absolute norm is<(68/25)b²<1. If one is positive and one negative,
(1) gives absolute norm<(68/25)*(2/7)<1. Both contradict nonzero
integral norm. Hence the conjugates are-a1,-a2, 0<a_j<a. The norm
is positive and<(68/25)a²<2, so it is1. The integral trace is1 or2:
it is>202/125-2a>0 and<68/25<3.

Trace2 would give a1+a2=t-2<18/25 while a1*a2=1/t>25/68;
the arithmetic-geometric mean inequality makes this impossible.
For trace1, t<1+2a=318/125<51/20, and
t(t-1)²>=4 forces t>9/4. The second integral coefficient is

    S=t-t²+1/t,              -4<S<-2.

The bounds follow by monotonicity on9/4<t<51/20. Therefore S=-3,
but X³-X²-3X-1=(X+1)(X²-2X-1) is reducible. This excludes the
far endpoint3 role over an actual cubic field at the NEW loss.

## 3. Near endpoint3 and all other helper roles

For a DOTTED near endpoint3, gamma>1 and j=m-1. Since m>=5,
j/(2m)>=2/5 and delta<=3/5. At gamma>=21/20 this gives

    R²>=(29/20)²/(3/5)=841/240>(1853/1000)².

Thus a weak dotted case has2<t<21/10. The independent Schur
calculation inb57 gives every other conjugate in(-3/2,0).
In degree3 write them-a,-b. The integral trace T=t-a-b is0,1 or2,
because-1<T<21/10. If T=1 or2,

    0<Norm(t)=t*a*b<=t(t-T)²/4<1.

If T=0, the positive integral norm M<=t³/4<3, so M=1 or2.
The second coefficient S=M/t-t² obeys-4<S<-3 in BOTH cases:
for M=1 its lower bound is10/21-441/100>-4 and its upper
bound is1/2-4<-3; for M=2 the corresponding bounds are
20/21-441/100>-4 and1-4=-3, with strict inequality. No integral
S exists. Hence the near dotted case is not cubic.

The remaining complete role checks at loss1/2 are:

* D-position3, near or far: dotted R>=(m+1)/sqrt(2m)>=6/sqrt10
  exceeds1853/1000. The finite own>=7 conjugate exclusions ofb57
  Section5 do not depend on loss; own4/5 locks the field to a
  quadratic field. The same finite-label statement treats near
  endpoint3 after the preceding new dotted argument.
* E-position3: finite white classification gives exactly
  (m,j)=(6,3),(7,3),(7,4),(8,3),(8,5). For own>=4,
  gamma>=1/sqrt2. At m6,j3, R=2gamma+1/2>=sqrt2+1/2>1853/1000.
  At m7 the bound is sqrt7*(1/sqrt2+3/14), and at m8 it is
  4*(1/sqrt2+3/16), each larger. Own3 is already affine.
* Endpoint4: near own3 has R>2 and far own3 is affine. Every
  own>=4 has R>=6/sqrt10>1853/1000.
* Two white endpoint3 contacts: own3 is affine; every own>=4
  has R>=q/sqrt2+1/q>2>1853/1000.

Finally the complete isolated FINITE cyclotomic argument ofb57
Section6 leaves only5,8,12, all quadratic. It is independent of
weakness. Thus in an ACTUAL CUBIC field every .50-weak central
index has an actual isolated DOTTED opposite. This is the complete
fresh inventory needed for the new computation.

## 4. The expanded isolated domain

For an isolated dotted cap put beta=4gamma² and u=4/beta. The actual
foot gives height loss<.191u+sqrt(u(1-u)). For0<u<=1/4 this is at
most.04775+sqrt3/4<.5. Weakness at the new threshold therefore forces

    4<beta<16,       0<beta_sigma<19/35.

The raw time, actual field, integral beta and per-index Schur condition
are exactly those proved in the earlier specification; their formulas
are independent of the new loss. In cubic degree, trace and symmetric
sums give the COMPLETE finite box

    p=X³-TX²+SX-M,
    T=5..17, S=1..17, M=1..4: 884 rows.

Indeed T<16+38/35<18, S<32*(19/35)+(19/35)²<18 and
M<16*(19/35)²<5. Keep precisely two roots in(0,19/35), one in(4,16),
and at index i require the two small roots in(0,38/[i(19-i)]).

For completeness the worker also records the isolated QUADRATIC
outer box T=5..16,M=1..8:96 rows, with one small root and one
physical root in the same intervals. This does NOT enumerate
nonisolated quadratic weak caps or prove that those are absent.

## 5. Algorithm and the exact proposition to be checked

The new worker imports the unchanged rational polynomial primitives
of `scripts/check_a18_isolated_cap_domain_and_norm_screen_20261021.py`
(`29da7434634c80fb001784d9d5f0d69c3d6247a7b6de1e55e2a0eb6b7ff22932`).
It uses the exact same irreducibility tests, Sturm root counts,
multiplication-matrix norm verification, squarefree witnesses and
necessary discriminant/norm groups. Only the coefficient boxes and
the physical root upper endpoint change as explicitly stated above.

The proof of termination and complete domain coverage is the same
finite-loop and rational-Sturm argument, now for980 rows. The worker
records all kept polynomials, index cells and pair groups. It also
matches every saved .52 cubic anchor against every expanded .50
cell by the SAME necessary signature. For each anchor it records
every surviving target cell, including any with a different polynomial
or k. No expected compatibility pattern is hardcoded.

If an anchor's entire match list has the same polynomial and k as
the anchor, then the reflected integrality lemma excludes every
distinct weak target: the anchor itself is the only possible .50-weak
index in that case. This implication will be used ONLY for anchors
whose actual output establishes that pattern. No norm pass is called
a K-square, and no retained polynomial is called realizable.

The new worker, its imported primitive source, both specifications,
b57 and the old output are pinned before and after execution. Use
one new output and one new guarded59-second/1GiB receipt. No field
list, floating-point optimizer or algebra package is needed.

The energy consequence requires the following separate legal same-x
reference and is not inferred merely from these finite domains.

## 6. Conditional energy closure, with every coefficient paid

Assume the computation shows that EVERY saved .52 cubic anchor matches
only the same polynomial and k in the expanded domain. The reflected
integrality argument then proves that its index a is the ONLY possible
.50-weak central index. This implication uses the complete fresh cubic
inventory above, as well as the computed necessary comparison.

Set the exact rational floors

    H5=2.214,H6=2.365,H7=2.473,H8=2.543,H9=H10=2.577,
    H_(19-i)=H_i.

They are strictly below sqrt(D_i)-.50: adding.50 gives exactly
the oldb57 floor plus.52, whose five positive residual numerators
are displayed there. Interpolate linearly to heights0 at indices0,19.
The successive slopes are.4428,.151,.108,.070,.034,0 and their
reflections. Thus the representing dual coefficients are nonnegative
and supported on5..14. Its exact energy is

    Efull=2.214²/5+.151²+.108²+.070²+.034²
         =1.0208802.

If6<=a<=13, omit just that knot and join its two neighbors linearly.
Concavity is preserved and the omitted index gets coefficient zero.
The energy debit is one quarter of the squared difference of the two
old neighboring slopes. Their largest difference here is.043.
The resulting reference energy is therefore at least

    1.0208802-.043²/4=1.02041795>1.

Every remaining positive coefficient has an actual strict strong
floor. Completing the coordinate square contradicts||x||<1.

Suppose instead a=5; a=14 follows by reversing the old chain. Use the
reference w0 with the first SIX heights linear from0 toH6, then the
knots6..14 just displayed, then linear to19. It is concave, has no
coefficient at5, and has exact energy and endpoint coefficient

    E0=59124131/60000000>.9854,
    b6=(H6/6-(H7-H6))/2=1717/12000>2/15.

The six possible .52 cubic anchor polynomials have coefficients(T,S,M)
(9,6,1),(11,7,1),(12,7,1),(12,8,1),(13,8,1),(14,9,1).
Their values at8 are respectively-17,-137,-201,-193,-257,-313.
All their other roots are<19/35, so their unique physical beta exceeds8.
The already proved .52 upper bound beta<40/3 then gives

    3/10<u=4/beta<1/2.

For this ACTUAL isolated cap set

    U=omega5-u*n5,
    Y=U-(5/6)*omega6,
    g=sqrt(D5-u)-sqrt(u(1-u)).

The actual feasible foot gives x.U>=g. Since w0's first six coordinates
are constant, w0.Y=0. The raw dual identities give

    Y²=5/3-2u+u²<347/300,
    g>sqrt(D5-1/2)-1/2>53/25.

For the last strict comparison use D5-1/2=261/38>(131/50)²,
so its square root exceeds2.62. Choose mu=4/25. The vector

    W=w0+mu*Y

uses a new actual-foot coefficient mu and subtracts5mu/6=2/15
from b6, leaving every old coefficient nonnegative. The SAME original
x therefore satisfies

    x.W>E0+mu*(g-5H6/6),
    W²=E0+mu²Y².

No unknown separate energy is added. Completing ONE square yields

    ||x||²>E0+2mu*(179/1200)-mu²*(347/300)
           =60211491/60000000=1.00352485>1.

This closes the endpoint as well. Consequently, IF the new complete
anchor comparison has the stated property, no actual cubic whole-A18
farthest vertex exists. Combiningb57 with the explicit rational Meyer
obstruction would then force ACTUAL field degree exactly2 at such a
vertex. This is not a whole-A18 exclusion or an arithmetic bound17;
the quadratic branch is separate. Pending execution and full review,
the conditional inference is not promoted to an accepted theorem.

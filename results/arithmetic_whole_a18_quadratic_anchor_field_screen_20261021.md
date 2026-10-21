# Whole A18: quadratic anchor matching and a universal endpoint foot

21 October2026. Root supporting proof and algorithm specification,
pending full outside-author review and bounded execution. The target
is a necessary actual-field restriction, not a chamber construction or
a new dimension bound. The proposed finite output is not assumed true.

## 1. Already proved input domains

The fresh whole-A18 theoremb57 and its two complete isolated domains
are the inputs:

* `results/arithmetic_whole_a18_point52_actual_field_degree_reduction_20261021.md`
  SHA `b57aa75adbde5c1a1b6f78354c58bba47d03b2047a923438b133be9aee5f2707`;
* `results/arithmetic_whole_a18_isolated_cap_domain_and_norm_screen_20261021.md`
  SHA `06172cc3d1f1e218a29d7e20916f105366d10344f4e086a4aeb4d4eed45e75ac`;
* `results/arithmetic_whole_a18_point50_cubic_inventory_and_anchor_screen_20261021.md`
  SHA `bcefd51d39d814b11a37eb04dd08e17962a1e71caee932ebc0e35af10dd8d125`.

Their saved exact outputs, each with an actual guarded execution, are
`results/a18_isolated_cap_domain_norm_screen_20261021.json`
(SHA `671550c23549a9f55af602cd5f26578c5901e16646deecd540b328d8010386fc`)
and `results/a18_point50_isolated_domain_anchor_screen_20261021.json`
(SHA `8c3b9ffb07f7a6aef39463b464c3a5c18afe5fbd07929eb2cdf0b0e3b820beae`).

The first contains37 quadratic polynomials and314 .52 index cells;
the second contains59 quadratic polynomials and504 .50 index cells.
These are complete necessary domains for ISOLATED weak caps. They
do not enumerate the nonisolated quadratic types. The new computation
only joins these saved finite domains; it reruns no root enumeration.

The fresh .50 cubic inventory and its actual30-anchor comparison,
together with its written single-reference argument, exclude degree3.
Rational degree is excluded by full-Gram locking and Meyer. Thus a
whole-A18 farthest vertex now has ACTUAL quadratic field. The next
two sections are a fresh complete reduction when that field is neither
Q(sqrt2) norQ(sqrt5).

## 2. Every nonisolated .50-weak quadratic cap is silver or golden

Use the N19 unit-normal notation, i=5..14. At loss1/2 a weak cap
has edge ratio R<1853/1000. The complete finite white list and actual
own0/parallel/affine exclusions are those proved independently of
loss inb57 and rechecked in the fresh .50 cubic source.

For a far endpoint3 helper normalize overK to square2, and put
t=2gamma. Its actual field isQ(t), t is integral, t>=phi>202/125,
and the identity R²>=2t-2 gives t<68/25. The threshold-independent
conjugate strip is-193/250<t'<12/25. Thus its integral trace is1,2
or3. These three cases can be checked completely without a field list:

* Trace1: t<443/250. The nonzero integral norm t(1-t) is strictly
  between-2 and0, hence is-1. Thus t=phi and K=Q(sqrt5).
* Trace2: t(2-t)<1, and t<68/25 makes it>-2. Its only possible
  nonzero integral value is-1. Thus t=1+sqrt2 and K=Q(sqrt2).
* Trace3: t>3-12/25=63/25 and t<68/25. Hence
  0<t(3-t)<(63/25)(12/25)=756/625<2. Its integral norm is1.
  Thus t=phi² and K=Q(sqrt5).

In trace2 the zero norm t=2 is impossible for a primitive quadratic
parameter, as well as being an excluded parallel own contact. No
positive integral norm is missed because its upper bound is strictly1.

For a near dotted endpoint3 helper the fresh weak bound gives
2<t<21/10, and its nonphysical conjugate is in(-3/2,0). Its trace
is1 or2. Norm t(1-t) lies in(-2.31,-2), or norm t(2-t) lies
in(-.21,0), respectively. Both intervals contain no integer.

For near/D-position FINITE helpers, the entire finite own>=7 exclusion
inb57 Section5 depends only on the conjugate Schur interval, not on
the loss. The remaining own4/5 values generate exactly the silver
or golden field. Dotted D-positions are strong sinceR>=6/sqrt10.
E-positions, endpoint4 and two-endpoint3 roles are all .50-strong
by the explicit complete checks in the fresh .50 source. Every role
has therefore been covered, including the new far-helper phi² type.

For an isolated finite cap, the complete cyclotomic argument leaves
labels5,8,12 after affine exclusions. They generate Q(sqrt5),Q(sqrt2),
Q(sqrt3), respectively. Label8 is incompatible with the raw diagonal
>=70/19. The remaining native12 case is strong by an ACTUAL finite
two-wall projection, which we give explicitly. Put gamma=cos(pi/12),
c=A*p-gamma*omega_i, A²=gamma²D_i-1 and delta=1-gamma²>0.
The spacelike vector xi=c+gamma*n_i has square delta, p.xi=-A,
and nonpositive products against every wall outside this finite core.
Projection of p to the intersection of the two core walls is

    q=(p+(A/delta)*xi)/sqrt(1+A²/delta).

It is future unit timelike, annihilates both core normals, and is
feasible against every other actual wall because its coefficients on
p and the two core normals are positive. Farthest comparison gives

    x.(omega_i-n_i)>=sqrt(D_i-1)-tan(pi/12).

Since x.n_i>0 and
sqrt(D_i)-sqrt(D_i-1)<.191, its height loss is
<.191+tan(pi/12)=.191+2-sqrt3<.459<.50.
Thus no isolated finite cap outside the silver/golden fields is weak.
No dotted-foot formula with gamma^(-2)>1 is used.

Consequently, if the actual quadratic field is neither silver nor
golden, EVERY .50-weak central cap is isolated dotted. The .52-weak
index forced byb57 is also isolated dotted, and belongs to the saved
314-cell domain. For this anchor beta<40/3, so u=4/beta>3/10.

## 3. A universal legal endpoint reference

The following elementary bound works for EVERY isolated dotted anchor
at index5, without a field or polynomial restriction:

    g(u)=sqrt(140/19-u)-sqrt(u(1-u))>211/100,
    0<u<1.                                             (1)

For u<=4/5, the sum of denominators
sqrt(140/19)+sqrt(140/19-u) exceeds2.714+2.562=5.276>100/19.
Thus sqrt(140/19)-sqrt(140/19-u)<.19u. The elementary circle identity
for s=sqrt(u(1-u)) gives

    .19u+s<(.19+1.018)/2=.604,

because1+.19²<1.018². Since sqrt(140/19)>2.714, (1) follows.
For u>=4/5, sqrt(140/19-u)>sqrt(121/19)>2.523 and
sqrt(u(1-u))<=.4; this gives the stronger g>2.123.
All displayed decimal comparisons are strict rational square tests.
The circle identity follows by expanding
(2u-1)²+4s²=1 and Cauchy--Schwarz; it is also proved in the
already checked continuous-loss kernel.

Suppose the anchor is the ONLY .50-weak central index. If it is
internal, the one-omitted .50 reference of the fresh cubic source
has energy>=1.02041795>1. For endpoint5, use exactly that source's
concave w0 with H6=2.365,H7=2.473,H8=2.543,H9=H10=2.577,
reflected H11..13, H14=2.214 and linear end tails. It has

    E0=59124131/60000000,
    b6=1717/12000>2/15.

For the actual cap put U=omega5-u*n5 and Y=U-(5/6)omega6.
The same-vector actual-foot pressure gives x.U>=g(u), while
w0.Y=0 and Y²=5/3-2u+u²<347/300, since3/10<u<1.
Choose mu=4/25. The subtraction5mu/6=2/15 is paid by b6, and
every other old coefficient has its valid strict strong floor.
Since g(u)-5H6/6>167/1200, one completed square gives

    ||x||²>E0+2*(4/25)*(167/1200)
                   -(4/25)²*(347/300)
           =60019491/60000000=1.00032485>1.

Reflection treats endpoint14. Thus ONE isolated .52 anchor with no
other .50-weak central index is impossible, regardless of its actual
quadratic field or the precise continuous parameter u.

## 4. The complete saved-domain comparison and conditional conclusion

For every one of the314 saved .52 QUADRATIC cells, form its necessary
signature(degree, squarefree polynomial discriminant, squarefree
absolute time-norm numerator). List EVERY one of the504 expanded
quadratic cells with the same signature. Keep all matches, including
different polynomial/k types. The fixed input arrays, integer signatures
and loops make this a terminating complete join of the two saved
necessary domains; their generation proofs are the two earlier sources.

Every actual same-field pair must pass this test. In quadratic degree
the discriminant squarefree part is exactly the squarefree D inQ(sqrtD).
The new checker reports, without assuming it, whether every old anchor
outside D2,5 matches ONLY its own polynomial and k. If this predicate
holds, its only possible distinct target is the reflected index. The
already proved Cartan integrality obstruction excludes that target
(the OLD anchor M<=7 is not divisible by19). Therefore the anchor
is the only .50-weak central index. Section3 contradicts its existence.

IF the exact saved-domain comparison establishes that predicate, a
whole-A18 farthest vertex must therefore have ACTUAL field Q(sqrt2)
orQ(sqrt5). This is a necessary restriction over ALL initial actual
admissible fields. It does not exclude those two fields, wholeA18,
any other retained forest, or dimension18. An unsuccessful predicate
would instead leave the explicitly recorded additional matching groups.

The two saved outputs, their two worker sources, all three proof inputs,
this new proof and the new checker are hashed before and after the
bounded run. Use a new compact output and receipt under the unchanged
59-second/1GiB/one-thread/exclusive-lock/resource gates. No saved
mathematical output or frozen manuscript is overwritten.

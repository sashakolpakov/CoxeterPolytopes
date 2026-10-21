# Cubic-seven extension of the two-component D fork half-refund

19 October 2026. **New written theorem candidate, awaiting independent
review.** Paper/source work only: no enumeration, scientific job, Lean,
TeX or build was run. No earlier packet has been changed or silently
broadened from its eight-field statement.

## 1. New field, new baseline, and exact rule

Let P be a compact arithmetic Coxeter20-polytope with ACTUAL ground
field K7=Q(cos(2pi/7)). Start from the activated A4+A1 transfer in

    arithmetic_k7_active_a4_exception_transfer_20261019.md,
    SHA256 cb912246464e30badcb82d6927dbf1fc20f375e40335d196e7f5d5fd6b42575d.

Its validity and capacity are a separate new proof dependency, not
assumed to have completed independent review. The charge has no D5
bonuses. It is the distance-at-most-two baseline plus an entire
I2(7/14) unit and the activated half-coupon on an A4 endpoint and
an entire isolate. It pays all triangles and quadrilaterals.

At a vertex whose ENTIRE forest is one of

    D_r+D_(20-r), r=10,...,16,
    D12+E8, D13+E7, D14+E6,

select the indicated rank-r D component, choosing either component
in D10+D10. Replace its two fork/short-leaf weights1 by1/2.
The following is a complete standalone consumer proof over K7.
It does NOT use any later whole-A/B/D20 refund or apply a published
eight-field result to a ninth field by declaration.

Use the notation a for the selected fork, b,e for its leaves,
a-l_1-...-l_L for its long arm, N=r-2>=8,L=N-1, and T_s for
the complete other factor. The white graph of an affected face is

    W=A1(e)+A_L+T_s.                              (1)

It has no connected A3 component. An active A4+A1 coupon has
complete white graph A3 plus its orthogonal spectators, by the
sealed-core theorem. Hence there is NO active coupon at any corner
of a face with white graph(1). This is whole-component separation,
not an assertion about whether a large component contains an A3
subdiagram.

## 2. Triangles and simultaneous half-sources

The complete K7 Lanner/triangle theorem supplies a7 or14 black
bond in every triangle. At its finite corner it is an entire
I2 component and retains weight2. A fork/leaf pair in a large
D component cannot be that high-I2 role. Thus every triangle
keeps its payment despite all simultaneous half reductions.

In a quadrilateral, two half-sources already pay1. If exactly
one is present, any other positive original baseline corner is
unchanged and pays at least1. It remains a quad at(a,b) whose
three other ORIGINAL weights vanish. With cyclic sides a,b,c,d,

    a.d=b.c=c.d=0, d.e=d.l_1=0.                   (2)

The latter conditions follow from the marked distance-two test
at(a,d). All other restrictions from the original zeros are kept.
Finite W+a+d has the same COMPLETE field-independent first-d
inventory as before:

    I no white contact; F far endpoint3; B far endpoint4;
    T far-next3; J spectator-only; bridge far3+spectator-long3.

This is a finite Coxeter classification at rank20, including the
entire spectator. A connected path-plus-spectator continuation
must be D20. Its bridge is simple and uses a D long endpoint;
no E bridge or unlisted interior exceptional extension exists.

We now prove that all six rows are excluded over K7. The old
field-independent Schur identities and actual finite inventories
are reused with their original hypotheses. Every field-dependent
step is explicitly audited below.

## 3. Isolated d, no long or spectator contact by c

This is the only no-long row needing a new native field entry.
Use the exact identities in
`arithmetic_rank20_isolated_side_no_long_contact_20261019.md`,
SHA256 `af2ff3c9c4667b037ea93ed1afa8f1dcaf364c0ac41e894e3541b2bcde1df879`.
Let x=-a.c, gamma=-e.c, s=r/(4N),

    R=(x+gamma/2)^2/(1-gamma^2),
    delta=4(-b.d)^2-4=1/(R-s)>0.

Here d touches only b, so beta=-b.d>1 follows from the ACTUAL
compact quad's common perpendicular side, not an ambient-wall
assumption. All nonphysical delta conjugates lie in(-4,0).
The source's full old-star strict positivity and cyclic identities
are unchanged in K7.

If gamma=0, its finite-x norm argument and dotted-x denominator
ideal argument use no field exclusion except K!=Q. Thus they
apply directly to K7.

For gamma>0 and x=0, K=Q(gamma^2), and every nonphysical
gamma-square is less than2/r<=1/5. The K7 native e-labels are
3,4,6,7,14. The first three squares are rational, forcing K=Q.
For label7, its conjugate cos^2(2pi/7)>1/4; for label14,
cos^2(3pi/14)>1/2. Both violate the bound1/5.

For finite x=cos(pi/j), j>=4, the source proves the SAME
field equality K=Q(x^2), conjugate interval(x^2)'<r/(4N)<=5/16,
and complete trace/totient elimination. It leaves only
j=4,5,6,8,12, whose square fields are Q,Q(sqrt5),Q(sqrt2)
or Q(sqrt3), never K7. In particular it explicitly excludes
j=7 and14 BEFORE the old eight-field check. If j=3, the
source gives K=Q(gamma), with every other gamma conjugate in
(-1,-1+4/r) contained in(-1,-3/5). Only gamma=cos(pi/7)
could generate K7, and its positive conjugate cos(3pi/7)
contradicts that interval. This covers all finite/zero x.

For dotted x>=1, retain the source's GENERAL ideal norm lemma:
with u=4gamma^2-4 and U=|Norm(u)|,

    delta>4(N/8)^3/U,  delta<=48/85.              (3)

The full actual field has degree3. Native label14 is not a
prime power, so U=1; label6 has U=1 also. These contradict(3).
Native label7 has U=Phi_7(1)=7, so

    delta>4/7>48/85,

since340>336. This is the needed NEW label7 step. The rational
labels3,4 are handled by the source's denominator/Bezout arguments,
which are totally real field-independent and cover every N=8,...,14.
Thus the whole isolated/no-long row is excluded over K7.

## 4. Isolated d, long contact by c

The complete finite inventory and proper-block identities in
`arithmetic_rank20_isolated_long_finite_product_20261019.md`,
SHA256 `9644c76993b26d97803686787c5a87c284051184aa08cd8d2d7d6250a06e1ece`,
remain valid. This includes the N=8 finite E8 contacts at l_3,l_5.

For a simple long contact and e contact absent or simple, the
source's exact centered moment has degree at most3 and proves
K=Q(x) for finite x. Only j=7 can now generate K7. Write
eta=cos(3pi/7)>1/8. To prove the last inequality without a
numerical angle estimate, eta is in(0,1/2) and satisfies
P(X)=8X^3-4X^2-4X+1=0. This polynomial is strictly decreasing
on[0,1/2] and P(1/8)=29/64>0, so eta>1/8.

The required nonphysical interval is |x'+kappa|<rho. In every
simple-long row except the far endpoint, kappa>=rho (strict
for the additional e3 rows and both finite E8 rows). For the
far endpoint without e contact,

    rho-kappa=(sqrt(N+1)-1)/(2N)
              =1/[2(sqrt(N+1)+1)]<=1/8.

Thus the positive conjugate eta violates EVERY such interval.
This excludes the new j=7 case; the source's zero-x affine and
positive-index witnesses are field-independent.

For endpoint4 and no e contact, let w=2sqrt2*x in O_K7,
an actual cycle scalar. The far-end degree-three norm exclusion
in the source is field-independent. At the near end it proves
K=Q(w), so Q(w^2)=Q(8x^2) is Q or K7. The COMPLETE native
finite j list is3,4,6,7,14. Their exact norms in K7 are

|j|3|4|6|7|14|
|---|---:|---:|---:|---:|---:|
|Norm(8cos^2(pi/j))|8|64|216|8|56|

The last two follow from Norm(4cos^2(pi/j))=Phi_j(-1), equal
to1 and7 respectively. Every entry except64 is not the square
of an integer, whereas it must equal Norm(w)^2. Label4 gives
w=2 rational, contradicting K=Q(w). This supplies the explicit
K7 replacement for the old K9-specific line. It does NOT assume
sqrt2 lies in K7.

A simple endpoint plus an e4 bond contains an actual old-path
cycle with scalar sqrt2. Since sqrt2 is not in K7, this entire
attachment is impossible independently of x. Finally all dotted
long contacts, including E8, are excluded by the separately
proved universal totally-real norm/unit packet
`arithmetic_rank20_isolated_d_dotted_long_unit_closure_20261019.md`,
SHA256 `5c583ce01db154908d5d0e139c61f4986cf1491cd4f3f991f0c7e26c8a5a92c6`.
That proof has NO native field-list or integral-basis restriction.

## 5. Isolated d, spectator contacts by c

The complete spectator inventory is unchanged and simple. Its
one-root residual t_T is at most1/2. For spectator-only c, the
full-star integral squares w=4x^2,z=4beta^2 generate the same
actual K. The finite-x norm bound is below4*(1/8)^2<1.
For dotted x, the existing coupled norm inequality at degree3 is
at most256/507<1, again a contradiction. No K7-specific finite
label has escaped these continuous inequalities.

When c also contacts e simply or a path endpoint simply, the
complete residual/odd-unit arguments in
`arithmetic_rank20_isolated_d_spectator_contact_obstruction_20261019.md`,
SHA256 `4dd1934e502251339b019b116ecc854bac495cf4a686f371ac6015afa542b13f`,
exclude every x over EVERY totally real field. They include all
three e3 residuals1/4,1/8,1/12 and both path/spectator bridges.
The path/spectator finite diagram is D18+A1, not D19+A1.
The underlying spectator-only calculation is Section3 of
`arithmetic_rank20_fork_spectator_continuation_research_20261019.md`,
SHA256 `95b0102731bb12a43416197260d65ff2998cc0a92db189722bee8ad43bf2d754`.
This finishes ALL isolated-d contacts over K7.

## 6. The five nonisolated first-d rows

Here is the full field audit of the remaining sources. Their
finite continuation and Lorentz determinant calculations are
field-independent; only the listed last field checks change.

* F and T, finite beta: their actual simple cycle gives
  K=Q(beta), and every other beta conjugate lies in(-1/2,1/4).
  In K7 only beta=cos(pi/7) is native, and cos(5pi/7)<-1/2
  violates that interval. Native label14 has degree6 and cannot
  be this cyclic scalar in the cubic field. Beta=0 already fails
  the physical full-star cost.
* T, beta>=1: its COMPLETE c inventory, including spectator
  contacts, and exact positive-root bounds are unchanged. Every
  e-label is physically cut to at most5, or to rational/golden
  full-star determinants; these force Q or Q(sqrt5), not K7.
  All remaining cases have the stated principal affine witnesses.
* F, beta>=1: the near and bridge bounds are field-independent.
  No-path c leaves native labels at most6 and forces Q or
  Q(sqrt5), or an affine witness. In the shared-far-end case,
  A=2x has2<=A<3 and all other |A_i|<1/2; its degree-three
  integer norm is below3/4. Finite x=cos(pi/7) fails the same
  tiny conjugate window. Thus the complete F proof applies.
* B, finite beta: the ACTUAL scalar A=2sqrt2*beta satisfies
  K=Q(A)=Q(A^2) and -8/r<A_i<0. Its complete cosine-square
  table explicitly retains labels7 AND14; each has a nonphysical
  A_i^2>2, contrary to A_i^2<16/25. The only old exceptional
  label12 forces Q(sqrt3), not K7, so that corner profile never
  arises. No unsupported beta-in-K claim is used.
* B, beta>=1: its complete finite c inventory and Q<1/2 bound
  leave only affine witnesses or a golden full-star field. These
  exclude K7 with no new native label.
* J: finite beta has an integral square with physical value<4
  and the two other conjugates<2/5; its norm is below16/25.
  Thus the old quadratic label12 profile is absent in K7.
  At beta>=1, the complete spectator-orthogonal and joint
  E6/E7/E8 extension proofs leave only native labels at most5,
  rational/golden exact fields, or the old affine witnesses.
  Its REQUIRED E8 coverage addendum remains included.
* Bridge: finite beta is killed by the same small norm in degree3.
  Its beta>=1 proof is pure Lorentz/affine geometry with the
  complete finite D18 continuation. It has no additional field
  premise except exclusion of Q in the degree-one finite case.

These statements reuse the exact source identities, not their
old eight-field theorem conclusions. The relevant source hashes
are the local ledger in Section7 of
`arithmetic_rank20_two_component_d_fork_half_refund_20261019.md`,
SHA256 `38654093cee3ee7419fb92c91d860545a614fda8fd55ef222a8369b300488f01`.
In particular the J theorem is bound together with its separate
E8 addendum. Equations and inventories were reread in full for
this audit; no new source execution is asserted.

## 7. Consumer conclusion and remaining K7 forests

All six possible first-d rows in Section2 are impossible when
the other three original corner weights vanish. Thus a quad
with a single new half-source still receives at least1. Quads
with several new half-sources already paid1, and Section2 paid
every triangle. All larger-face requirements follow from
nonnegativity. The two fork half-refunds are therefore consumer
safe relative to the new K7 activated-A4 charge.

At each of these ten whole D/E forests there is no A4+A1 coupon;
its original charge is36. Two half reductions give charge35.
D10+D10 needs only an arbitrary choice of one of its components;
omitting reductions at another corner can only help payment.

After this standalone theorem the still heavy K7 rank20 forests
for this charge are EXACTLY

    A20,B20:37; D20:38;
    D19+A1:36;
    D18+I2(7),D18+I2(14):36.

Every other forest has upper35, by the complete frontier in the
activated-A4 transfer plus the ten reductions here. This is a
precise six-labelled-type frontier, not a dimension19 theorem.
Whole-A/B/D20 and the two small-spectator refunds must still
receive complete K7 consumer proofs. Existing eight-field source
statements are not automatically those proofs. New-rule composition
must also retain the actual active-coupon white graph and any
shared A/D zero-source interfaces. No global dimension ceiling
or priority claim is made here.

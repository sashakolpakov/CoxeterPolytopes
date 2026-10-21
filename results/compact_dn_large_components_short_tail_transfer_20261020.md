# Large D components: arbitrary spectators and the two short tails

20 October 2026. Complete paper candidate, awaiting full root and
outside-author review. No scientific computation, optimizer, Lean, TeX,
or accepted-source edit. The two averaged tail constructions were
developed with root; resource_cleanup independently checked their
algebra before this source was written.

**Theorem.** A vertex farthest from a strict interior point of a compact
hyperbolic Coxeter polytope cannot have a whole finite component D_r
with24<=r<=29, even with arbitrary orthogonal finite spectator
components. This is a component theorem, not a new ambient ceiling.

## 1. Normalization, inputs, and what is transferred

Number the D_r chain1,...,r-2, with fork leaves r-1,r at r-2. Use unit
old normals and Euclidean coordinates with

    n_i=(e_i-e_(i+1))/sqrt(2) (i<r),
    n_r=(e_(r-1)+e_r)/sqrt(2).

The spatial farthest vector satisfies

    ||x_D||<1, x1>...>x_(r-1)>|xr|,
    omega_i=sqrt(2)*prefix_i, D_i=2i (i<=r-2),
    omega_(r-1)=(prefix_(r-1)-e_r)/sqrt(2),
    omega_r=(prefix_(r-1)+e_r)/sqrt(2).

Thus C_ij=2min(i,j) for long indices, C_i,f=i for a fork f, and
C_f,f=r/2. Put h_i=x_D.omega_i and S_i=h_i/sqrt(2). Along the long
arm the S_i are strictly concave. Put phi=(1+sqrt5)/2. Call i weak if
h_i<=sqrt(2i)-61/100, and strong otherwise.

We use the following unchanged sources, all under results/:

* `compact_whole_b_bulk_weak_packet_and_endpoint_pressure_20261020.md`,
  SHA256 `e5a3c51858bb7b1e392b2828c874387be8d2044d5ae1be583a546b67eaddb54c`;
* `compact_whole_b_short_tail_extension_20261020.md`, SHA256
  `b6adc03a5950a60a1e263aeb79869785268f1bb7905e74ec47a784446c989c5c`;
* `compact_bn_large_components_arbitrary_spectator_transfer_20261020.md`,
  SHA256 `9852635a5271c1a733999fb9bd49cd4170d019ae3f37f2ae8bb4af246ebc79ed`;
* the foreign-far comparison
  `compact_classical_foreign_far_packet_strong_20261020.md`, SHA256
  `3766485d02300a78951e26e7a162086aa4ff74db9dc3aaa0842584ffb5fb1520`;
* the general isolated-clearing construction, also written in Section2
  of `compact_an_large_components_arbitrary_spectator_transfer_20261020.md`,
  SHA256 `8d3b999761be7dd46e7990973b3ec4d4a746d4d682b044d53695f4aae405e3e5`.

The generic actual-point input is
`compact_farthest_vertex_mixed_replacement_inequalities_20261020.md`,
SHA256 `3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84`;
reciprocal shadow is SHA256
`c913ec3bd8e0f26a8c679a1ede2b98b251ad3f36efeae60486213ddfa3e8badc`;
the orthogonal-component opposite obstruction is SHA256
`4216d06c2e03fd336898b8dac289b9969fa096c8a1d08bd67ddeb8ab5ec2f23d`.

For an actual outside unit wall c=A*p-v, A>0, write
gamma=-c.n_i and Q for the complete spatial inverse-Gram cost of its
contacts in all foreign components. A nonzero Q is at least1/4.
Deleting direction i gives an ACTUAL positive definite full Gram,
including every spectator. Each finite connected component is a tree,
so c has at most one contact with each retained old component. We use
the ordinary finite Coxeter classification with its actual native
labels; no projected wall is presumed native or actual.

For a finite old set J, if w is the positive old-normal correction
and delta=1-||w||^2, actual finite projection (delta>0) gives

    x_D.(v-w)>=H(A,delta),
    H(A,delta)=A/[1+sqrt(delta/(A^2+delta))]>A-sqrt(delta).       (1)

Here and below the displayed spatial vector is pure D whenever this
formula is used without foreign coordinates. If delta<0, the positive
combination xi=c+w is future timelike, actual feasible, and gives

    x_D.(v-w)>=A-sqrt(-delta).                               (2)

Feasibility holds against ALL unknown walls: the selected products
are nonpositive by the displayed calculation, and every other wall
has nonpositive product with each term of this positive combination
of actual normals. Formula(1) uses the analogous positive acute
projection. These are actual points, not formal finite completions.

## 2. Isolated contacts and the foreign cost

Suppose the only D contact is gamma at i>=5. The exact affine
principals exclude own3 (two simple forks using the old right fork),
own4 (a terminal4 joined to the old fork), and own6 (the6,3 path).
Parallel actual normals give a positive feasible null sum and are
also impossible in a compact polytope. Hence gamma>=cos(pi/5),
finite or dotted. These principal obstructions ignore foreign nodes
and therefore remain valid with spectators.

For Q>0 write v=gamma*omega_i+vS and u=(1-Q)/gamma^2. Then
A^2=gamma^2*(D_i-u). For u<1 the actual point

    xi=c+vS+[(1-Q)/gamma]*n_i

has xi.c=xi.foreign=0, xi.n_i=gamma*(u-1)<0 and
xi^2=-gamma^2*u*(1-u). Its pure-D pressure is

    h_i>sqrt(D_i-u)-sqrt(u*(1-u)).                            (3)

The case u=1 gives a nonzero feasible future null vector and is
forbidden. For u>1 the full finite set c,n_i,all foreign old walls
has residual gamma^2*(u-1)>0, and(1) gives

    h_i>sqrt(D_i-u)-sqrt(u-1).                               (4)

For D_i>=10 the loss in(3) is less than
u/6+sqrt(u*(1-u))<=(1+sqrt37)/12<3/5.
For(4), Q>=1/4 gives u<=3/phi^2=6-3phi<7/6,
sqrt(u-1)<=2-phi<2/5, and
sqrt(D_i)-sqrt(D_i-u)<1/5. Thus all foreign isolated caps are strong.

At Q=0 the accepted isolated actual-point estimates are unchanged.
In particular finite5 needs only the LEFT H4 old segment. Its vector
is gamma*(3omega_(i+1)/4+omega_(i-3)/4), whose x-product is at most
gamma*h_i by long-arm concavity. This remains valid through i=r-3:
i+1<=r-2 is still a long node. Native7 uses the old left A2 negative
core, native>=8 the finite old singleton, and dotted caps the actual
opposite Lanner foot. The accepted bounds require only D_i>=85/11,
whereas here D_i>=10. Every Q0 isolated cap is therefore strong too.

## 3. Complete bulk inventory, 5<=i<=r-5

The retained D part is A_(i-1)+D_s, s=r-i>=5. On the left the
complete one-node finite extensions are endpoint3 (A_i), endpoint4
(B_i), next-to-end3 (D_i), or the admissible E6/E7/E8 positions3.
Their Schur residuals are, respectively,

    (i+1)/(2i), 1/i, 2/i, and 1/4,1/7,1/16.

For a left simple3 at j counted from the far end1,
sigma=gamma+j/(2i). For a left4 at an endpoint,
sigma=gamma+j/(sqrt(2)*i). In each case cosh L_i is
sqrt(2i)*sigma/sqrt(delta), and foreign cost subtracts Q from delta.

On the right, a simple3 at its long leaf i+1 makes D_(s+1).
A simple3 at a fork leaf makes E_(s+1), only when s=5,6,7.
These are ALL right extensions. At s>=8 no fork attachment is finite.
If c contacts both old pieces, the connected full white rank is r>=24.
Finite classification forces a simple3 endpoint bridge into D_r;
a fork attachment or a left internal/4 attachment cannot make a large
finite diagram. In this bridge delta=1/(2i) and
sigma>=gamma+1/2+1/(2i), so it is strong. At gamma0 the full all3
cycle is affine instead. No foreign extension of this white D_r at c
is finite: c is an internal node of the long arm, and an extra arm
would create two forks. Thus no foreign contact was dropped here.

The zero-own branches also remain covered. A left internal3 contains
the affine two-fork D principal; left4 contains the terminal4 affine
fork; a right long3 gives the affine two-fork diagram. A right fork3
gives the exact affine E8 principal with arms1,2,5. Two old-piece
contacts give an all3 affine cycle. The only zero-own branch lacking
such a principal is far3@1 with foreign contacts. Its height obeys

    A^2=Q-1/2>0, Q<(i+1)/(2i), hence A^2<1/(2i)<=1/10.

For any contacted foreign node j, its squared derivative is at least
1/4>A^2. Reciprocal shadow makes c an actual opposite in that foreign
component too, contradicting the orthogonal-component opposite
obstruction. With no D contact the derivative in direction i is zero.

When there is no right contact, own4 is ALWAYS forbidden, even with
a near left3: select c and the entire right D_(r-i+1) starting at i,
and OMIT the left helper. Its long-leaf inverse unit diagonal is2,
so the4 attachment has cost1 and is an exact connected affine
principal with positive null vector. Own6 gives the local6,3 affine
path; own3 with a left3 gives an all3 cycle. Far left4 with own3
gives the exact affine F4 path old3,old2,old1,c,oldi. Near left4 with
own3 is allowed, but its actual edge has
cosh L_i=(i-1)+i/sqrt(2), strong for i>=5. All remaining near3,
D-position and E-position estimates are exactly the accepted B
ones with the preceding left residuals; subtracting Q only increases
the edge length. They are strong.

For right long3, its cost is1/2 and sigma=gamma+1/2.
Own0 and3 are the affine fork/triangle; gamma>=1/sqrt2 gives
sigma/sqrt(delta)>=1+1/sqrt2>100/61, so it is strong.
For right fork3, its cost is s/8 and sigma=gamma+1/4.
Own0 gives affine E8; own3 gives the all3 cycle; own4 gives the
exact affine F4 path old(i-2),old(i-1),oldi,c,fork; own6 gives6,3.
Thus gamma>=cos(pi/5), and

    sigma/sqrt(delta)>=(cos(pi/5)+1/4)/sqrt(3/8)>1.72>100/61.

Extra foreign contacts only decrease delta. A finite exceptional
extension could include foreign nodes at c; this inequality includes
them and does not assume they are sealed.

It follows that the ONLY weak bulk profile is

    v=gamma*omega_i+omega1/2,
    A^2=2i*gamma^2+2gamma-1/2,
    Q=0, own5, finite>=7, or dotted,
    101/125<gamma<6/5.                                      (5)

Indeed a foreign far3@1 packet is strong by the pinned foreign-far
lemma, valid here for5<=i<=29. The gamma upper bound and all small
left-position edge estimates are the accepted B estimates; they use
the identical prefix inverse matrix, not a B right-tail assumption.

## 4. The right D4 and D3 finite inventories, including bridges

First let i=r-4, so the right old part is D4. Its one-node finite
extensions are simple3 at ANY leaf, making D5. Long-leaf3 is already
handled. With a left contact, finite rank r forces a simple3 endpoint
bridge into D_r; its delta=1/(2i) and
sigma>=gamma+1/4+1/(2i) make it strong. Gamma0 gives an all3 cycle.
The white D_r bridge has no finite foreign extension at its internal
c, as above. With no left contact, a fork3 has

    delta=1/2-Q, sigma=gamma+1/4.                            (6)

Its own0/3/4/6 are the exact affine principals already listed. If
Q>0, delta<=1/4 and gamma>=cos(pi/5), immediately strong. At Q0,
gamma>=cos(pi/7) gives sigma/sqrt(delta)>1.61, whereas the exact
edge target T(sqrt(2i)) is<1.577 for i<=29 (the pinned foreign-far
lemma gives this target bound). Thus only own5/fork3,Q0 needs a
new actual-point argument, supplied in Section5.

Next let i=r-3, so the right old diagram D3 is the A3 path with its
middle node the old center and its ends the two forks. Its COMPLETE
one-node extensions are: simple3 at the middle (D4), simple3 at an
end (A4), end4 (B4), or end5 (H4). More than one contact to this
old A3 would make a cycle and is not finite.

Middle3 has delta=1/2-Q and sigma=gamma+1/2. Own0 is the affine
four-leaf star at the old center; own3 is the affine triangle;
own>=4 is strong. If it also joins the left path, the white diagram
is D_r with a simple3 at a left endpoint, delta=1/(2i), and the
larger bridge numerator, strong; the gamma0 bridge has an all3 cycle.

End4 has delta=1/4-Q and sigma=gamma+1/(2sqrt2). With gamma0,
the old center has two simple leaves (oldi and the other fork) and
the third arm is the contacted fork--4--c. Its cost at the center
is1/4+1/4+1/2=1: an exact affine principal. Thus gamma>=1/2,
giving ratio>=1+1/sqrt2, strong. Adding a left path through c would
put the4 internally in a long finite diagram and is impossible.

End5 is H4, sealed against any further white contact. Importantly,
gamma0 here need NOT be affine and is NOT discarded. Even then

    delta=1-3phi^2/8=(7-3sqrt5)/16,
    sigma>=phi/4,
    sigma/sqrt(delta)>100/61,

so the direction is strong for every gamma>=0.

For end3 the exact data, with NO left contact, are

    delta=5/8-Q, sigma=gamma+1/4.                            (7)

If it also contacts the left path, finite classification at rank r
forces simple3 at a left endpoint, making A_r. The root c is the
fourth node from the right end and has a long left arm, so adding a
foreign arm cannot make D or any finite exceptional diagram. The
bridge is sealed. Its residual and numerator are

    delta=1/8+1/(2i), sigma>=gamma+1/4+1/(2i).

Gamma0 gives the full all3 cycle. Otherwise gamma>=1/2; at i>=21,
delta<.16 and sigma>.75, so the ratio is>1.875, strong.

For the sole end3 case(7), own0/3/4/6 are the same affine E8,
cycle, F4 and6,3 principals as before. If Q>0, delta<=3/8 and
gamma>=cos(pi/5), again strong. At Q0, gamma>=1 gives ratio at
least sqrt(5/2)>1.58>1.577, strong. It remains exactly own5, or
finite own>=7 with gamma<1, fork3, Q0, no left contact. These are
the only new cases in Sections5--6.

## 5. Own5 plus fork3: two actual H4 projections

This applies to BOTH short layers, i=r-4>=20 or i=r-3>=21;
the estimates below even work for i>=19. Put gamma=phi/2.
The full D projection is v=gamma*omega_i+omega_f/2, where f is
the contacted fork. Its height satisfies

    A^2=i*(phi+5/8)-1/2  (right D4),
    A^2=i*(phi+5/8)-5/8  (right D3).                         (8)

On the left select oldi,old(i-1),old(i-2). On the right select
oldi, its next chain neighbor, and the following chain neighbor;
at the D3 layer the last one is the OTHER, uncontacted fork.
Both selected old paths have only the own5 contact with c, so the
actual H4 correction and residual are

    w_+-=gamma*(3n_i/2+n_side+n_second/2),
    delta5=1-3gamma^2/2=(7-3sqrt5)/16>0.

Every correction coefficient is positive. Formula(1) applies to
each ACTUAL finite projection, retaining all foreign old walls.
Since Q0 these points have pure D displacement. Average to obtain

    y=v-(w_++w_-)/2,
    x_D.y>=H(A,delta5)>A-sqrt(delta5).

The exact cross product w_+.w_-=3gamma^2/4 gives

    ||y||^2=A^2+1-15gamma^2/8<A^2.                          (9)

The old corrections have zero first5 coordinates. Consequently
y.u5=sqrt(10)*(gamma+1/4), where u5=omega5/sqrt(10).
Using 1.6<phi, gamma<.81, sqrt10<3.17, and i>=19,

    A^2>41.65, A>6,
    ||y||^2>A^2-1/4>41.4>6.4^2,
    sqrt(delta5)<1/7,
    (A-sqrt(delta5))/||y||>41/42>.97,
    rho=(y/||y||).u5<3.3602/6.4<.55.                       (10)

Section3 already gives the complete bulk weak inventory. Its actual
mixed points imply nonadjacency in that bulk. The accepted weak5
argument uses strong6 and one of16/17; all are bulk even at r24,
since r-5>=19. All its selected old normals and pure prefix duals
exist in D. Thus5 is strong, independently of either short-tail
argument, and x_D.u5>1-.61/sqrt10>.8.

Equations(10) would give x_D.(u5+y/||y||)>1.77, but
||u5+y/||y||||^2<3.1 and 1.77^2=3.1329>3.1. This contradicts
||x_D||<1. Both own5/fork3 short-tail cases are excluded.

## 6. D3 tail, finite own>=7: two actual negative A2 cores

Now i=r-3>=21, cos(pi/7)<=gamma<1, Q0, and the sole helper is
fork3. Put

    A^2=i*(2gamma^2+gamma+1/8)-5/8,
    d=4gamma^2/3-1>0,
    w_+-=gamma*(4n_i/3+2n_side/3),

where side is old(i-1) or the old center i+1. The points
xi_+-=c+w_+- have selected old products0, product with c equal
-d, norm -d, and positive future time coefficient A. The fork
helper is outside each old A2 selection; its negative product only
helps. All coefficients are positive, hence BOTH xi are feasible
against every actual wall, including all spectators. Formula(2)
therefore gives, on averaging,

    y=v-(w_++w_-)/2,
    x_D.y>=A-sqrt(d).

Here w_+.w_-=8gamma^2/9, whence

    ||y||^2=A^2+1-14gamma^2/9<A^2.                         (11)

Again y.u5=sqrt10*(gamma+1/4). At i21 exactly,

    ||y||^2-40*(gamma+1/4)^2=4gamma^2/9+gamma+1/2>0;

larger i adds a positive multiple of 2gamma^2+gamma+1/8.
Thus rho=(y/||y||).u5<1/2. Also d/A^2 increases with gamma:
the derivative numerator, after positive denominator factors, is

    [i*(4gamma^2+13gamma+3)-5gamma]/3>0.

It decreases with i. Its maximum on the enlarged interval gamma<=1,
i>=21 is therefore 1/195 (at gamma1, A^2=65,d=1/3). Consequently

    x_D.(y/||y||)>1-1/sqrt195>.928.                        (12)

The normalized vector is well defined: already gamma>.9 makes
A^2>54, while (11) subtracts less than1. Strong5 gives
x_D.u5>1-.61/sqrt10>.807, since 10*(.193)^2>.61^2.
Together with(12) this implies

    x_D.(u5+y/||y||)>1.735,
    ||u5+y/||y||||^2<3,
    1.735^2=3.010225>3,

again contradicting ||x_D||<1. This excludes ALL remaining finite
fork3 cases together, with no per-label table or field restriction.

## 7. Prefix closure through21, and the limited frontier corollary

Sections2--6 prove the complete weak inventory(5) for EVERY
5<=i<=r-3. Its accepted native-gap plus actual mixed-replacement
proof makes weak indices nonadjacent. All prefix inverse products
are C_ij=2min(i,j); the common-helper Schur matrix, positive
displacements, and far-helper actual finite/negative-core points
are therefore exactly those of the accepted B proof. Retaining all
spectator walls adds an orthogonal positive definite block. None
of these actualness assertions is inferred merely from a formal
completion or from a metric facet.

For r>=24, index21 is now in the complete range. If21 is strong,
the accepted concave-prefix reference has energy E21>1.009. At most
eight nonadjacent internal weak holes reduce it by less than8/3200;
its nonnegative dual coefficients give ||x_D||^2>1.

If21 is weak,20 is strong and the same actual far-helper pressure
gives x1>.424. Its finite5 H4 uses old1,old2,oldi; its7/8 core
uses old1,oldi; its other finite/dotted feet use oldi. The dual
identity needed to expand them uses simple neighbors i-1,i+1,
both still available at i21 even when r24. Concavity is along long
prefixes, so the accepted estimates apply without a new constant:

    E20>.99752, f5/5<.361,
    ||x_D||^2>.99752+(5/4)*(.424-.361)^2-7/3200
               =1.00029375>1.                               (13)

There are at most seven nonadjacent holes in6..19. The fixed x1
constraint preserves the decreasing slopes and positive duals of
the reference curve; each isolated omitted floor loses<1/3200.
Both alternatives contradict ||x_D||<1, proving the theorem.

For clarity, the index5 argument was established in Section5 using
only the bulk through17; it does not use21 or either averaged
tail contradiction. This makes the assembly noncircular.

One accepted-input corollary records the resulting frontier, WITHOUT
an unconditional dimension reduction. The source
`compact_dn_actual_projection_nine_fourteenths_20261020.md`, SHA256
`99bde3780f9d720d9809f5d8bea6795e4a77fcaacf8ede5871f0f831fb7823e9`,
proves that in ambient dimension n>=15, if a farthest finite forest
has no A/B component of rank>=5, then its only remaining large
component is a unique D_r with r>=max(10,n-5). At n29 this forces
r>=24, contradicted here. Thus EVERY farthest vertex of any compact
29-dimensional candidate must contain A_r with5<=r<=20 or B_r
with5<=r<=22, using the pinned accepted A/B spectator transfers for
the upper limits. The restricted class admitting a farthest forest
without such A/B components consequently has ceiling28 after the
existing ambient ceiling29 is imported. The unrestricted compact
and arithmetic ceilings are not changed by this corollary.

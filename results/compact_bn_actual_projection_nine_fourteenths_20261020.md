# Uniform actual B projection above 9/14, with arbitrary spectators

20 October 2026. Complete paper candidate, awaiting independent review.
No scientific execution, optimizer, host classification, Lean, TeX or
accepted-source edit. Derived in discussion with root.

**Theorem.** At a farthest vertex of a compact Coxeter polytope, every
WHOLE B_r component, r>=9, satisfies sum_(j=1)^7 x_j^2>9/14 in the
coordinates below, hence also ||x_B||^2>9/14, with arbitrary orthogonal
finite spectator components. This is not an
edge-only covering rate or an ambient dimension bound. In particular
the proof does not drop spectators from the whole-B23 inventory.

## 1. Conventions and actual-geometry inputs

Use unit roots n_i=(e_i-e_(i+1))/sqrt(2), i<r, and n_r=e_r. For long
indices i,j<r the dual entries and directions are

    omega_i=sqrt(2)*prefix_i, C_ij=2min(i,j),
    u_i=prefix_i/sqrt(i), h_i=x_B.u_i.

Put x7 for the first-seven-coordinate projection and suppose for
contradiction ||x7||^2<=e=9/14. For every i<=7, 0<h_i<1, and h_i
majorizes the actual edge half-tangent. Distinct actual opposites obey
F(h_i)F(h_j)>=u_i.u_j, F(h)=2h/(1+h^2).

Write Q for a wall's FULL spectator inverse-Gram cost. Each nonzero
component contribution is at least1/4; one below1/3 can only be an
A1_3 cost1/4. Indeed native contact squares are>=1/4 and inverse
unit-Gram diagonals are>=4/3 in a nontrivial component. Positive acute
projection and the actual reciprocal-shadow criterion are used below.
Affine principal diagrams are excluded by compactness. The external
classical structural inputs are finite spherical classification and
the Lanner order bound5. Actual opposite stars are indefinite; every
minimal nonpositive principal diagram is Lanner.

## 2. Complete common internal-pair inventory

Let4<=i<=r-4 and suppose an actual c closes i and i+1. Deleting i
retains right B_(r-i), rank>=4. Its finite one-node extensions permit
no contact or just one simple3 at LONG endpoint i+1. Internal/multiple
contacts or bonds>=4 cannot extend it finitely. F4 does not contain
B4 as a proper principal diagram. Deleting i+1 retains left A_i,
allowing at most one contact:3, or4 at an endpoint. Here i>=4 excludes
the smaller H/F path exceptions.

If c contacts i+1 by3 but not i, the principal fork with leaves
(old i,c) and right arm terminating in4 is exact affine B. Hence c
contacts i. Own3 gives the affine triangle, leaving own4. No second
left contact is possible. Any spectator contact puts a nonpositive
path with bonds at least3,4,3,3 in the ACTUAL finite replacement
deleting i+1. Thus no spectator contact is possible.

With no right contact, an internal left3 at j>=2 gives the same affine
fork ending in4. A left4 gives an affine4,3,3 star or the double4
terminal path. The sole remaining left contact is3 at1, possibly
with spectators. Every Lanner in the full star contains i+1, because
the actual replacement deleting i+1 is positive definite. Connectedness
then requires the full path c,1,...,i+1 of i+2>=6 vertices. Spectators
attach to B only through c and cannot shorten it. This contradicts
the Lanner bound. With no left contact either, c has zero B derivative
and cannot close i+1.

Therefore the ONLY common profile is4@i,3@(i+1), all other old products
zero. The retained A_i cost is i/(i+1), residual1/(i+1), giving

    cosh L_(i+1)=i+(i+1)/sqrt(2)>7.                 (1)

So4/5 are distinct under e, and6/7 are distinct for r>=10.

### The complete r9 exception for6/7

Deleting6 retains B3. Besides long-end3@7 it permits short-end3@9,
making F4. The long-end case gives the preceding4@6,3@7 profile.
The F4 branch seals all left and spectator contacts while6 is deleted:
no connected finite rank>=5 diagram extends F4. Thus c can contact
only9 and the removed6. Own0 gives the affine F4 path6,7,8,9,c.
Deleting7 retains A6, so own6 can be only3 or4;4 gives the affine
4,3,3 star at6. The sole extra common profile is consequently3@6,3@9.
Its actual direction6 has

    delta=1/4, sigma=1/2+1/(2sqrt(2)),
    cosh L6=(1+1/sqrt(2))*sqrt(12)>11/2,
    tau6^2>9/13>9/14.                              (2)

This also contradicts the projection upper. Thus actual6/7 are
distinct for EVERY r>=9, including all exceptional right-tail roles.

## 3. Fixed low-prefix pressures and their checked B hypotheses

The computations imported below use only the actual4/5 and6/7
distinctness just proved, C_ij=2min(i,j) for i,j<=7, and the local
unit roots n1,n2,n3. Their actual feasible points use positive
combinations of c and these roots, so every unknown wall and every
spectator inequality is respected. No D-tail inventory is imported.
Every dual and Cauchy vector in these fixed arguments is supported on
the first seven coordinates, so the Cauchy upper needs only ||x7||^2<=e,
not a bound on the uncontrolled remaining B coordinates.

First, ANY actual edge2 with cosh L2>=1+sqrt(2) is impossible under e:
it gives h2^2>=sqrt(2)-1. Distinct6/7 requires h6+h7>3/2. For
w=u6+u7, orthogonal Cauchy about u2 has

    rho=1/sqrt(3)+sqrt(2/7)<1.114,
    S=29/21+4sqrt(2/21)<2.62,

and its decreasing envelope gives x_B.w<1.494416<3/2. The exact
seals are Sections2--3 of
`compact_dn_long_end_native_eight_plus_energy_20261020.md`, SHA256
`f103e5378a6997fc97e10f802502cf882ce7e8b9c59a724567d868a82af8606b`.
After its initial floor that argument has precisely the B hypotheses
established here, with no dependence on a D tail.

Three other fixed profiles are also excluded under e:

* Isolated5@1,Q0: actual shadow makes c close3 and4, with
  h4^2>=phi-1, phi=(1+sqrt(5))/2. Orthogonal Cauchy about u4 gives
  h6+h7<362087/250000<3/2.
* Isolated5@2,Q0, actual common2/3: shadow1 and actual edge3 give
  x1>31/50 and S3>131/100. Distinct4/5 gives h4+h5>57/40. The
  fixed positive dual using x0=(31/50,69/200,69/200,0,...), coefficient
  1/6 and the residual of u4+u5 has exact lower Cauchy margin3/700.
* Contacts4@2,3@3: actual deletion3 has white cost2/3, so Q=0 or
  one A1_3 cost1/4. At Q=1/4, cosh L3=4+3sqrt(2)>7. At Q0 put
  a^2=(5+4sqrt(2))/2. For J={2,3}, w=sum(B_J^-1 beta)_j*n_j
  has residual delta=-sqrt(2)/3. Its ACTUAL future point c+w yields

      x_B.U>=a-sqrt(sqrt(2)/3),
      U=((2+sqrt(2))/2)e1+((1+sqrt(2))/3)(e2+e3+e4).

  The selected products are (c+w).n2=(c+w).n3=0,
  (c+w).c=delta<0; positive coefficients verify all other actual
  inequalities. The edge2/3 floors and this point, with positive-dual
  coefficient1/16, have Cauchy margin8825113/11200000000>0.

These exact computations and their positive-coefficient/rational seals
are pinned in
`compact_dn_low_index_three_profile_farthest_obstruction_20261020.md`,
SHA256 `b03a84986519cf9a9c16b1c84461dd0ed52806d9ada505a1066f8f750f345c01`.
Every derivative, shadow length, dual and residual used there lies in
the identical first-seven-coordinate subspace just verified. Actual
closure is supplied by the generic shadow theorem, not a D-specific
finite extension. This explicitly checks the transfer's hypotheses.

Finally isolated7@1,Q0 is excluded by an actual point. Put
t=cos(2pi/7),g=cos(pi/7),c=sqrt(t)*p-g*sqrt(2)*e1. Then

    xi=c+n1/g+(1-t)*n2/g
       =sqrt(t)*p-V/sqrt(1+t),
    V=t(e1+e2)+(1-t)e3,
    xi.c=xi.n1=0, xi.n2=(1/2-t)/g<0,
    xi^2=-(1-t)*(2t-1)/(1+t)<0.

It is future, and positive coefficients prove all other actual wall
inequalities. Thus

    x_B.V>=sqrt(t(1+t))-sqrt((1-t)(2t-1)),
    x_B.(V/||V||)>.73.

Distinct6/7 and decreasing Cauchy about V/||V|| give x_B.w<591/400,
contradicting >3/2. The complete cubic/rational seals are pinned in
`compact_dn_native_seven_actual_dual_foot_energy_20261020.md`, SHA256
`a39306eb2637acf75959530fb5cf5d37dd77d28efdcf6cf2060708a63620d769`.
All local normals, heights and actual feasibility above are explicitly
the same in B; no D conclusion is assumed merely from similarity.

## 4. Complete common1/2 inventory with all spectators

Deleting1 retains B_(r-1), rank>=8. The right contact is none or3
at LONG leaf2. Deleting2 retains old1 as A1, so gamma=-c.n1 is zero
or finite native<1. If c2=3, own0 is the affine fork ending in4 and
own3 the affine triangle. For gamma>=1/sqrt(2), the full Q gives

    cosh L2=(1+gamma)/sqrt(1-gamma^2-Q)>=1+sqrt(2),

contradicting Section3. This exhausts the right-contact case.

With no right contact, c has just its node1 B contact. Its full height
is a^2=2gamma^2+Q-1>0 and actual deletion2 requires gamma^2+Q<1.
Own0 cannot close a B direction. For own3, a^2=Q-1/2<1/4; a contacted
spectator has squared derivative>=1/4 and hence is also closed by c.
This violates the farthest orthogonal-component theorem. Own4 gives
the double4 affine path, own6 the6,3 path. For own5,Q>0 implies
Q>=1/4 and cosh L2>=phi^2>1+sqrt(2). At Q0 the isolated5 pressure
above applies. For m>=7, Q<1-gamma^2<1/4 forces Q0. Native m>=8
gives h2^2>=sqrt(2)-1; native7 is excluded by its actual xi point.
Thus actual1/2 opposites are distinct.

## 5. Complete common2/3 inventory with all spectators

Deleting2 retains B_(r-2), rank>=7: right contact none or3 at LONG
leaf3. Deleting3 retains A2, permitting one contact3/4/5 at1 or2;
two contacts give a nonpositive triangle, label>=6 a3,6 path.

If c3=3 but c2=0, the fork with leaves(old2,c) and right arm ending
in4 is affine. Thus c2 is3/4/5 and c1=0. Own3 is an affine triangle;
own4 is the4+3 profile above; own5 has A2 white cost phi^2/3>3/4,
forcing Q0, and cosh L3=(2+5phi)/sqrt(2)>7. There is no fork or
exceptional E extension of this retained B component of rank>=7.

With no right contact, no left contact gives zero B derivative. The
remaining cases are exactly3@1 with spectators,4@1,5@1,3@2,4@2,5@2.
The first has a^2=Q-1/2<1/6, since actual deletion3 retains A2 cost
1/3; a contacted spectator therefore also closes, forbidden. The4@1
case is the double4 affine path;3@2 the fork ending in4;4@2 the
affine4,3,3 star. Both5 cases are excluded by Section3, with Q0
forced by the actual deletion3 white cost phi^2/3. Thus actual2/3
opposites are distinct. No spectator profile was dropped.

## 6. Fixed three-pair dual and scope

The actual pairs1/2,2/3,6/7 are all distinct. Their correlations are
sqrt(1/2),sqrt(2/3),sqrt(6/7). The exact positive-dual calculation in
`compact_dn_actual_projection_nine_fourteenths_20261020.md`, SHA256
`99bde3780f9d720d9809f5d8bea6795e4a77fcaacf8ede5871f0f831fb7823e9`,
is directly applicable as ALGEBRA on these same prefix directions.
Its independently sealed pair sums are >1.091,>1.265,>1.512. For

    W=.1*(u1+u2)+.055*(u2+u3)+.305*(u6+u7),

they give x_B.W>.639835 and W^2<.636435825, with exact gap

    .639835^2-(9/14)*.636435825
       =70423123/280000000000>0.

Cauchy gives (x_B.W)^2=(x7.W)^2<=||x7||^2 W^2, contradicting the
assumption ||x7||^2<=9/14. This completes the stronger localized
theorem and its component-projection consequence.
No B5..B8 result, rank/14 rate for every B_r, or ambient dimension
ceiling is claimed. The whole-B23 proof and its absent-spectator
inventory were not used.

The generic actual-shadow/point inputs are SHA256
`c913ec3bd8e0f26a8c679a1ede2b98b251ad3f36efeae60486213ddfa3e8badc`,
`3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84`,
and farthest orthogonal-component obstruction SHA256
`4216d06c2e03fd336898b8dac289b9969fa096c8a1d08bd67ddeb8ab5ec2f23d`.
All geometry is actual; no arithmetic or literal-host cutoff is used.

# D8: three consecutive directions force covering energy greater than4/7

20 October2026. Written proof candidate, awaiting outside-author review.
No scientific execution, optimizer, Lean or TeX. Other finite components
at the actual compact Coxeter vertex are arbitrary. No arithmeticity.

Use the component covering energy from
`compact_universal_e6_e7_e8_component_energy_20261020.md`, SHA256
`acf261d97f6eb1a11a64502e17e079ef7f418b9746fbf85dbdf6915750fe9a18`.
For normalized outgoing directions u_i and actual edge lengths L_i,
write tau_i=tanh(L_i/2). The energy is the minimum squared norm of a
component vector x satisfying x.u_i>=tau_i. Distinct actual opposites
obey the full-tangent law tanh(L_i)tanh(L_j)>=u_i.u_j. These are pinned
geometric inputs; finite spherical classification and Lanner order<=5
are explicit classical inputs as in that source.

**Claim.** Every D8 component has covering energy greater than4/7.
This would add D8 to the rank/14 forest estimate. It does not assert
a new ambient dimension bound or estimates for D9 and higher.

## 1. A common opposite for directions4 and5

Label the D8 chain1-2-3-4-5-6, with fork leaves7,8 at6; all old bonds
are3. Suppose one actual wall c is opposite to both directions4,5.
Deleting4 retains a right D4 on5,6,7,8 and a left A3 on1,2,3.
Finite classification allows c no right contact, or just one simple3
contact at one of the three right leaves5,7,8. A middle contact makes
the affine degree-four star; a leaf4 or higher is already nonpositive.
Two right contacts make a nonpositive cycle. Deleting5 retains the
left A4 on1,2,3,4, so c has at most one left contact there.

If the right contact is3 at5, c.4 must be negative: otherwise
c,4,5,6,7,8 form the exact affine two-fork D diagram. In the actual
replacement deleting5, its contact at endpoint4 is3 or4. A3 would
make the exact affine triangle c,4,5, so it is4. Every spectator
contact is excluded by the retained path x,c,4,3,2 with bonds at
least3,4,3,3, which is nonpositive. Thus the complete profile is

    c.4=-1/sqrt2, c.5=-1/2, all other old products0.       (1)

If the right contact is3 at7 (or8), any left contact must also be3.
Indeed an internal4 contact gives the retained4,3,3 star, and an
end4 gives an internal4 in the retained path of rank8, containing
the affine F4 path3,3,4,3. Larger contacts only worsen its positive
nonpositive witness. A left3 creates an exact all3 cycle through
c,7,6,5 and the left chain, impossible in a compact polytope. With
no left contact, the full D8 plus c is the exact affine E8 tree with
arms1,2,5. Spectator contacts cannot remove that principal affine
subdiagram. Thus no such right-fork contact is possible.

With no right contact, a left3 at an internal node j>=2 creates the
exact two-fork affine D tree, using leaves c,j-1 and7,8. An internal4
is forbidden in the retained A4 replacement; an end4 at4 gives the
full4,3,3 star, and an end4 at1 gives the terminal-4 affine D extension.
The sole possible left contact would be3 at1, with arbitrary spectator
contacts. A minimal nonpositive principal diagram must contain c and5,
since the actual replacement deleting5 is positive definite. It is
connected and has the path c,1,2,3,4,5 of six vertices; spectators
cannot shorten this path. This contradicts the classical Lanner bound.
No left or right contact gives zero derivative in direction5 and cannot
close that edge. This exhausts all cases and proves(1).

The unit inverse-D entries are C44=8,C55=10,C45=8. At the actual5
endpoint the white residual is1/5, so

    cosh L5=(4sqrt2+5)/sqrt2=4+5/sqrt2>7.

Thus a common4/5 opposite forces tau5^2>3/4. A common5/6 opposite
also forces energy>3/4 by the already reviewed D_r center/neighbor
classification:
`compact_dn_center_neighbor_common_opposite_energy_20261020.md`, SHA256
`e371611c4903ad362654ee0d4656457140cdc5e4374b4e19a2697718c2d679b1`.

## 2. Consequences of energy at most4/7

Suppose the covering energy is at most e0=4/7 and choose a minimizing
x. Put h_i=x.u_i. Then0<h_i<1. By Section1 both pairs4/5 and5/6 have
distinct actual opposites. The full-tangent law and monotonicity of
F(h)=2h/(1+h^2) on[0,1] give

    F(h4)F(h5)>=q=sqrt(4/5),
    F(h5)F(h6)>=k=sqrt(5/6).                         (2)

For any h,h' in(0,1), F(h)F(h')>=t implies

    hh'>=p(t)=t/(1+sqrt(1-t))^2.                    (3)

This follows by putting p=hh', using h^2+h'^2>=2p, and solving
4p/t>=(1+p)^2. In particular h4*h5>=p(q).

Let lambda=(1+k)/2, c0=sqrt(lambda), s0=sqrt(1-lambda). The orthonormal
vectors w=(u5+u6)/(2c0) and v=(u5-u6)/(2s0) satisfy

    u5=c0*w+s0*v, u6=c0*w-s0*v,
    u4=q*u5+(1/sqrt5)*t0,

where t0 is unit and orthogonal to w,v. The last identity follows from
u4.u5=q and u4.u6=q*k, using the standard inverse-D entries. Write
x=a*w+b*v+z*t0+y with y perpendicular to these three vectors. Hence
a^2+b^2+z^2+||y||^2<=e0 and

    p=h5*h6=lambda*a^2-(1-lambda)b^2
      <=lambda*e0-b^2-lambda*z^2.                  (4)

The second inequality in(2) is equivalently

    psi(p)=4p/k-(1+p)^2>=4s0^2*b^2.

The concave quadratic psi vanishes at p(k), with derivative there
4sqrt(1-k)/k. Its tangent line is an upper bound, so

    p>=p(k)+alpha*b^2, alpha=k*sqrt(1-k)/2.

Combining with(4), set

    T=(1+alpha)b^2+lambda*z^2
      <=delta=lambda*e0-p(k).                      (5)

The forthcoming rational seals give delta<1/400. This confines the
vector near the maximizing direction for5/6, with an additional
penalty for their difference; this penalty is needed.

## 3. An incompatible upper bound for h4*h5

Since h5<=sqrt(e0), the exact identity for u4 gives

    h4*h5=q*h5^2+(z/sqrt5)*h5
      <=q*lambda*e0 + A*|b|+B*|z|
                          -q*k*b^2-q*lambda*z^2,  (6)

where A^2=8/105 and B^2=4/35. For h5^2 use
lambda*a^2+2c0*s0*a*b+s0^2*b^2, bound |a|<=sqrt(e0), and retain
the negative squared terms from a^2<=e0-b^2-z^2.

Weighted Cauchy gives

    A|b|+B|z|<=sqrt(C*T),
    C=(8/105)/(1+alpha)+(4/35)/lambda<47/250.

Also q*k/(1+alpha)>71/100 and q>71/100. The quadratic loss in(6)
is consequently at least(71/100)T. Since q*lambda*e0<489/1000,

    h4*h5 <489/1000+sqrt((47/250)T)-(71/100)T.

This last function is increasing for0<=T<=1/400. Thus(5) gives

    h4*h5 <489/1000+sqrt(47/100000)-71/40000
      <489/1000+271/12500-71/40000
      =101781/200000<509/1000<p(q),                (7)

contradicting(3). This proves the claim, subject to independent review.

## 4. Complete rational seals

Squaring gives

    .8944<q<.8945, .9128<k<.9129,
    .29<sqrt(1-k)<.2953<.296.

Consequently .13<alpha<.14 and .95<lambda<.95645. These imply

    C<(8/105)/(113/100)+(4/35)/(19/20)<47/250,
    q*k/(1+alpha)>(.89*.91)/1.14>.71,
    q*lambda*e0<.8945*.95645*(4/7)<.489.

For(5), p(k)>.9128/(1.2953)^2 and lambda*e0<.95645*(4/7).
Their difference is less than1/400. One exact cross multiplication
for the only tight comparison is

    91280000*70000-167780209*38083=26300653>0.

Here38083/70000=.95645*(4/7)-1/400. Next sqrt(1-q)<.325, so
p(q)>.8944/(1.325)^2>.509. Finally(271/12500)^2>47/100000,
and the increasing-function test in(7) follows from
10sqrt(47/250)>71/100. Every decimal here denotes an exact
terminating rational. No sampled optimization or field approximation
is used.

If accepted, the pinned forest aggregation extends to D8 because
4/7=8/14. All other exceptional-defect and short-A1 branches remain
unchanged. In n>=15 a farthest vertex would then need an A/B component
of rank>=5 or a D component of rank>=9. This still leaves large
component families and does not prove a uniform dimension ceiling.

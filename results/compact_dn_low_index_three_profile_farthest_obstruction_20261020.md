# Three low-index D profiles force actual farthest projection above 9/14

20 October 2026. Paper proof candidate, awaiting independent review.
No arithmetic, scientific execution, optimization, Lean or TeX. Accepted
sources are unchanged. Let p be a farthest vertex of a compact Coxeter
polytope, x its common spatial farthest vector, and x_D its projection
onto a WHOLE D_r component, r>=9. Other finite components are arbitrary.

Use unit roots `(e_i-e_(i+1))/sqrt2` and
`(e_(r-1)+e_r)/sqrt2`. For i<=r-2,
`omega_i=sqrt2*(e1+...+ei)`, `u_i=(e1+...+ei)/sqrt(i)`.
Write phi=(1+sqrt5)/2, s=sqrt2, and e=9/14.

**Theorem.** Suppose one actual wall c closes directions 2 and 3,
and its D contacts are one of

1. a single 5 at node 1;
2. a single 5 at node 2;
3. 4 at node 2 and 3 at node 3.

Then `||x_D||^2>9/14`. The first two arguments actually give the same
conditional lower bound for the edge-covering energy. The third uses
an additional ACTUAL negative-core point; it does not assert an
edge-only covering-energy bound. This proves no unconditional D_r rate
or ambient dimension ceiling.

## Inputs and spectator costs

Use the actual reciprocal-shadow theorem SHA256
`c913ec3bd8e0f26a8c679a1ede2b98b251ad3f36efeae60486213ddfa3e8badc`,
the mixed actual-point theorem SHA256
`3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84`,
and the accepted internal common-opposite theorem SHA256
`ff1aaa687171eec9ebf2669746bdc37cd753dcbcef62890ebf49be20a791bd11`.
The r=9 center/neighbor pair additionally uses SHA256
`e371611c4903ad362654ee0d4656457140cdc5e4374b4e19a2697718c2d679b1`.
Their paths are the respective named actual-shadow, mixed-replacement,
internal-D and center/neighbor sources under `results/`.

Let Q be the full spectator inverse-Gram cost of c. It is nonnegative.
In an actual finite replacement c has at most one contact per old
spectator component, since two give a nonpositive cycle. A nonzero
component cost is at least 1/4; one below 1/3 must be a single A1_3.
Indeed native contact squares are at least 1/4, and an inverse unit
Gram diagonal in a connected component of rank at least two is at
least 4/3, by Schur complementation against any neighbor.

For cases 1 and 2 the replacement deleting 3 retains left A2 and has
white cost `phi^2/3>3/4`. Thus Q=0. In case 3 the left white cost is
2/3, so Q is either zero or a single A1_3 cost 1/4. If the latter
occurs, its actual direction-3 residual is 1/12 and
`cosh L3=4+3sqrt2>7`. Then `tau3^2>3/4>9/14` and the conclusion
already follows. It remains to consider Q=0 in case 3.

Assume for contradiction `||x_D||^2<=e` in each remaining case.
Put h_i=x_D.u_i. These satisfy the actual half-tangent floors.
Every distinct pair obeys `F(h_i)F(h_j)>=u_i.u_j`, where
`F(h)=2h/(1+h^2)`; positivity and h_i<1 justify monotonicity.
For correlation t this implies

    h_i*h_j >= p(t)=t/(1+sqrt(1-t))^2.

The derivation is `4p/t>=(1+p)^2`, using
`h_i^2+h_j^2>=2h_i*h_j`.

## 1. Isolated 5 at node 1

Put gamma=phi/2. Its unit height has square `(phi-1)/2`, and its
direction-4 shadow satisfies `T4^2=phi^3/4>1`, so it is the actual
opposite there too. The reciprocal actual length gives exactly
`tau4^2=phi-1`: square
`(sqrt(phi^3)-sqrt(phi^3-4))/2` and use
`phi^3=2phi+1`, `phi^-3=2phi-3`.

The 6/7 opposites are distinct under the assumed e bound. For r=9
this is the center/long-neighbor theorem; for r>=10 it is the
internal theorem at i=6. A common wall in either case forces energy
above 3/4. The 6/7 correlation is k=sqrt(6/7), and
`p(k)>9/16`, since `k>576/625`; the last squared comparison has
positive cleared difference `6*390625-7*331776=21318`.
Consequently `h6+h7>3/2`.

Let w=u6+u7, rho=w.u4 and S be its squared orthogonal residual:

    rho=2*(1/sqrt6+1/sqrt7)<787/500,
    S=16/21+4/sqrt42<69/50.

For the first bound use `1/sqrt6<409/1000`,
`1/sqrt7<378/1000`. For the second use `sqrt42>162/25` and
`16/21+50/81=782/567<69/50` (cleared difference 23).
Since `phi-1>309/500`, `sqrt(e)<401/500`, and
`e-(phi-1)<1/40`, orthogonal Cauchy gives

    x_D.w <= rho*sqrt(e)+sqrt(S*(e-(phi-1)))
          < (787/500)*(401/500)+93/500
          =362087/250000 <3/2.

Here `(93/500)^2>69/2000`; every radical bound follows by squaring.
This contradicts h6+h7>3/2. Only actual edge constraints were needed.

## 2. Isolated 5 at node 2

The height square is phi. Its direction-1 shadow has
`T1^2=phi/2<1`, hence the other actual opposite gives
`tanh^2 L1>=phi/2`. This yields `x1>31/50`:
`F(31/50)=3100/3461` and
`(3100/3461)^2<161/200<phi/2`.
The first cleared difference is 6541881.

The actual direction-3 edge has `cosh L3=sqrt2*phi^2>37/10`,
so its partial sum S3=x1+x2+x3 is greater than 131/100.
Indeed `3*tau3^2>81/47>(131/100)^2`, with cleared difference
3433. The cosh comparison uses
`sqrt2>707/500`, `phi>809/500`, giving product
`(707/500)*(1309/500)>37/10`.

Directions 4/5 have distinct actual opposites: a common one has
`cosh L5=4+5/sqrt2>7`, impossible under e. Their correlation is
q=sqrt(4/5). The accepted exact comparison `p(q)>.509` gives

    x_D.w=h4+h5>57/40,       w=u4+u5.

This also follows directly by squaring the formula for p(q), so
the terminating constant is not a numerical approximation.
Let v=e1+e2+e3 and k0=1/2+1/sqrt5<237/250. The projection of w
onto the first three coordinate subspace is k0*v, and its residual
square is `S=13/20+1/sqrt5<11/10`.

Choose the fixed positive-dual vector

    x0=(31/50,69/200,69/200,0,...),
    E0=||x0||^2=12449/20000>31/50,
    x0=(11/40)e1+(69/200)v.

The lower constraints above give x_D.x0>=E0. Put
`Delta=57/40-k0*(131/100)>9/50` and mu=1/6.
For `W=x0+mu*(w-k0*v)`, its coefficients of e1,v,w are
`11/40`, `69/200-mu*k0>0`, and mu. Therefore

    x_D.W > E0+mu*Delta,       W^2=E0+mu^2*S.

Cauchy under the assumed e would give the reverse of
`(E0+mu*Delta)^2>e*(E0+mu^2*S)`. The difference increases with
E0 in this range. Its lower rational seal is

    (31/50+(1/6)*(9/50))^2
      -(9/14)*(31/50+(1/36)*(11/10))=3/700>0.

This contradiction again uses only actual edge/shadow constraints.

## 3. Contacts 4 at node 2 and 3 at node 3, Q=0

The two actual edge floors are

    tau2^2=(1+2s)/7, tau3^2=(2s-1)/3,
    tau2*tau3=1/sqrt3.

Both inverse-Gram weights are positive. Their equality vector x0
in span(u2,u3) has
`E0=(6s-4)/7>.6407`; this alone is below 9/14.

The old A2 core J={2,3} has inverse
`[[4/3,2/3],[2/3,4/3]]` and contact vector `(1/s,1/2)`.
Its Schur residual is `delta=-s/3`. With the positive old correction
the actual vector xi=c+correction is future timelike, has square
delta, and satisfies every actual wall inequality. Thus its normalized
point lies in P. Farthest comparison gives

    x_D.U >= T=A-sqrt(s/3),       A^2=(5+4s)/2,
    U=((2+s)/2)e1+((1+s)/3)(e2+e3+e4).

This is an actual-point inequality, not a condition on all edge-cover
vectors. Project U onto span(u2,u3):

    rho2=(2s+1)/6, rho3=(1+s)/sqrt3,
    S=||U-proj U||^2=(7+4s)/12<211/200,
    Delta=T-rho2*tau2-rho3*tau3>3/50.

Here are complete simple rational seals. From 1.4142<s<1.4143,
we have A>2.3083, sqrt(s/3)<.6867,
tau2 in(.7395,.7396), tau3 in(.7806,.7808),
rho2<.6381 and rho3<1.394. Hence
`Delta>1.6216-.6381*.7396-1.394*.7808=.06122604>.06`.
All terminating decimals denote exact rationals; each root bound is
checked by squaring. The two old dual coefficients exceed .3 and
.52 (use u2.u3<.8166), whereas subtracting (1/16)rho removes less
than .04 and .088, respectively.

Set mu=1/16 and `W=x0+mu*(U-proj U)`. Its old coefficients remain
positive, so `x_D.W>=E0+mu*Delta` and `W^2=E0+mu^2*S`.
The difference in the ensuing Cauchy test is increasing in E0.
The exact lower seal is

    (.6407+(1/16)*.06)^2
      -(9/14)*(.6407+(1/256)*1.055)
      =8825113/11200000000>0.

This proves the actual farthest-projection conclusion in the final
case, with arbitrary spectators already handled at the start.

## Boundary

These three branches are completely excluded if the actual farthest
D projection has squared norm at most 9/14. They do not by themselves
exhaust all common-2/3 profiles or establish an unconditional rank/14
estimate. The complete small-left-chain classification is a separate
obligation. In particular no arithmetic whole-D long-end theorem or
noncofinite-subspace implication was imported.

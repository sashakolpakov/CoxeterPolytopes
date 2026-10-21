# An actual Lanner dual foot excludes the D native7 long-end cap

20 October 2026. Paper only. The decisive dual-foot calculation was
supplied by finite57 and independently reconstructed by root and guard.
No scientific execution, host classification, Lean, TeX or accepted
source edit. The actual-component statement below is universal.

**Theorem.** Let p be a farthest vertex of a compact Coxeter polytope
from a strict interior point. Suppose its finite diagram has a WHOLE
D_r component, r>=9. If an actual wall c has exactly one nonzero old
contact, a native7 at D node1, and is perpendicular to every other
D node and every spectator, then ||x_D||^2>9/14.

Consequently, under ||x_D||^2<=9/14 the actual opposites for directions
1 and2 are DISTINCT, by the complete arbitrary-spectator reduction
`compact_dn_common_one_two_native_seven_reduction_20261020.md`, SHA256
`ccea4a223a5c0dafc3f83166d19a4d2f5f1f18ae622860b3d2bc5dafd6eac95c`.
The proof does not assume another reflection normal is an actual wall.

## 1. The exact future timelike point is feasible in the full polytope

Use unit D roots n1=(e1-e2)/sqrt(2), n2=(e2-e3)/sqrt(2), with
the usual chain and fork numbering. Put

    g=cos(pi/7), t=cos(2pi/7), g^2=(1+t)/2,
    c=sqrt(t)*p-g*sqrt(2)*e1,
    xi=c+g^-1*n1+(1-t)*g^-1*n2.                       (1)

All coefficients in (1) are positive. Direct calculation gives

    xi.c=0, xi.n1=0, xi.n2=(1/2-t)/g<0,
    xi^2=-(1-t)*(2t-1)/(1+t)<0,
    xi=sqrt(t)*p-V/sqrt(1+t),
    V=t*(e1+e2)+(1-t)*e3.                            (2)

Its positive p coefficient makes xi future timelike. Every other
actual wall has nonpositive product with each of c,n1,n2, so also
with their positive combination xi. The three exceptional products
were checked explicitly in (2). Therefore the normalization of xi is
an ACTUAL point of the full polytope, including all unknown walls.
This is where the known Lanner dual foot supplies more than one-cap
edge lengths or a hypothetical second mirror.

Write the strict interior comparison point as z=E(p-x), with ||x||<1.
Comparison with p's maximum distance gives

    x_D.V>=L,
    L=sqrt(t*(1+t))-sqrt((1-t)*(2t-1)).              (3)

Only the D projection occurs, because V has no spectator component.

## 2. A rationally sealed normalized pressure

The native value t satisfies 8t^3+4t^2-4t-1=0 and lies in(1/2,1).
This polynomial is strictly increasing there. Its values at .6234
and .6235 are respectively -.000926392768 and .000105223. Hence

    .6234<t<.6235.

Set Q=||V||^2=3t^2-2t+1. The preceding interval gives

    sqrt(t*(1+t))>1.0059,
    sqrt((1-t)*(2t-1))<.305,
    .95<sqrt(Q)<.9588<.959.

For example t*(1+t)>1.01202756>1.0059^2;
(1-t)*(2t-1)<.0929955<.305^2; and
Q<.91925675<.9588^2. All endpoints are positive rationals.
Thus for v=V/||V||, (3) implies

    h=x_D.v>=L/||V||>.7009/.959>.73.                 (4)

## 3. A distinct actual pair contradicts that pressure

Suppose ||x_D||^2<=e=9/14. The actual opposites for directions6 and7
are distinct: a common opposite forces an actual squared half tangent
greater than3/4, by the center/neighbor theorem when r=9 and the
internal consecutive theorem at i=6 when r>=10. Their correlation is
k=sqrt(6/7). The full-tangent distinct-pair law consequently gives

    h6*h7>=k/(1+sqrt(1-k))^2>9/16,
    x_D.w=h6+h7>3/2,              w=u6+u7.           (5)

The last strict comparison follows from k>576/625, whose squared
cleared difference is 21318>0. Here h_i=x_D.u_i majorizes the actual
half-tangent floor; the pair law uses F(h_i)=2h_i/(1+h_i^2).
No full tangent is identified with a half tangent.

Write w=rho*v+w_perp and S=||w_perp||^2. With
a=1/sqrt(6)+1/sqrt(7), the exact projection is

    rho=a*(1+t)/sqrt(Q).

The elementary bounds .78<a<.787 and Section2 give
1.3<rho<1.35. Also ||w||^2=2+2sqrt(6/7)<31/8, hence

    S<31/8-(1.3)^2=437/200<11/5.

Orthogonal Cauchy now bounds

    x_D.w<=f(h)=rho*h+sqrt(S)*sqrt(e-h^2).

The function is decreasing for h>=.73: its derivative is negative
exactly when rho^2*e<||w||^2*h^2. Here the left side is less than
6561/5600, while the right side exceeds3*(.73)^2, a larger rational.
Thus (4), evaluated at .73, gives

    x_D.w<1.35*.73+sqrt((11/5)*(9/14-.73^2))
          <1.35*.73+.492=591/400<3/2,

contradicting (5). The squared final radical margin is exactly
.492^2-(11/5)*(9/14-.73^2)=277/1750000>0.

## Inputs and boundary

The full actual farthest/projection conventions and the distinct6/7
argument are pinned in
`compact_dn_long_end_native_eight_plus_energy_20261020.md`, SHA256
`f103e5378a6997fc97e10f802502cf882ce7e8b9c59a724567d868a82af8606b`.
Its complete common-opposite inputs are the internal theorem SHA256
`ff1aaa687171eec9ebf2669746bdc37cd753dcbcef62890ebf49be20a791bd11`
and center/neighbor theorem SHA256
`e371611c4903ad362654ee0d4656457140cdc5e4374b4e19a2697718c2d679b1`.
All these allow arbitrary finite spectator components.

The result closes the remaining common1/2 native profile under the
stated actual projection upper. It does not by itself establish a
rank/14 rate for all D components, a uniform ambient dimension ceiling,
or a classification of cubic reflection lattices. No literal-family
nonreflectivity input, candidate second wall, or new finite search is
used.

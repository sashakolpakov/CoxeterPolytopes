# A universal A3 energy above thirteen fiftieths

22 October 2026. **Complete manual candidate, pending full different-author
review and root acceptance.** Root authored this refinement. No worker,
enumeration, symbolic computation, mathematical import, test, compiler,
Lean or TeX job was used. All fractions and decimals are exact.

**Theorem.** At an actual farthest vertex of a compact hyperbolic Coxeter
polytope, every whole A3 component has projection energy greater than
13/50, with arbitrary finite spectators. Arithmeticity is not needed.

**Consumer.** Together with the accepted arithmetic A11 energy above
13/20, this excludes the exact arithmetic16 forest A11+O7+A3, where
O7 means OLD I2(7) with its ACTUAL NEW isolated-three opposite. A bare
old seven is not O7. This is a component improvement and one necessary
forest exclusion, not a new global dimension bound or a polytope example.

## 1. Exact interfaces and one original vector

All mathematical paths below are under `results/`.

| Input | SHA256 | Exact interface used |
|---|---|---|
|compact_universal_small_a_component_energy_20261020.md|801581c4b855289aca7982b2d0f4bb6f5849f151a78f40b86c9d403e84f415e6|Actual A3 inverse Gram, distinct-opposite tangent-product law, pair/log-concavity method, complete contact inventory, native gap, positive foreign cost strictly below height and actual affine exclusions; not its older quarter conclusion at a changed target.|
|compact_farthest_vertex_orthogonal_component_opposite_obstruction_20261020.md|4216d06c2e03fd336898b8dac289b9969fa096c8a1d08bd67ddeb8ab5ec2f23d|No actual opposite shared across whole orthogonal components; used with the preceding source to derive the positive-foreign and small-height isolation statements.|
|arithmetic16_a11_paid_helper_and_all_isolated_feet_uniform_point65_candidate_20261022.md|2e69ab9ccdc2fee7700ca12f89574ca59f93caef6f033fff990e4951663a290e|ONLY the accepted spectator-safe arithmetic A11 energy above13/20, for the final consumer.|
|compact_farthest_rank_four_forest_exclusion_20261020.md|68ad62a619603eb7e6c2acc5689c901dde3f9d0897a7c1580b3deb948f5832b6|ONLY the actual old-seven/new-three defect energy above9/100, for the final consumer.|

Compact acute simplicity, spherical/affine classification and the actual
native-wall gap are the inherited geometric boundaries. The selected
vertex p maximizes distance from a fixed strict interior point z, written
z=E(p-x), ||x||²<1. Every component energy and actual comparison below
concerns this SAME original vector; all whole components are orthogonal.

For the whole A3 component let C be the inverse UNIT Gram, in path order,

    C=[[3/2,1,1/2],[1,2,1],[1/2,1,3/2]],
    d_i=C_ii, u_i=omega_i/sqrt(d_i), e=||x_A3||².

For actual edge length L_i, farthest comparison gives

    x_A3.u_i>=tau_i=tanh(L_i/2)>0.

We assume e<=13/50 and derive a contradiction. Neither a rank-
decremented theorem nor an assumption about spectators is used.

## 2. The changed target still forces ONE actual opposite

Put H(t)=2t/(1+t²), increasing on (0,1). Its logarithm is concave,
with second derivative -1/t²-2(1-t²)/(1+t²)²<0. If two old edge
directions have correlation k, then Cauchy and farthest comparison give

    tau_i+tau_j<=sqrt(2e(1+k)),
    tanh(L_i)tanh(L_j)<=H(sqrt(e(1+k)/2))².

For adjacent A3 directions k=1/sqrt3. At e=13/50 the latter bound is

    5200(1+k)/(113+13k)² < k.                         (1)

Indeed the positive denominator clears to

    k(113+13k)²-5200(1+k)=(22876k-12662)/3>0.

Here k>4/7, since 3<49/16, and 22876*(4/7)=13072>12662.
The distinct-actual-opposite inequality would instead require the
product of tangents to be at least k. Thus adjacent directions share
their actual opposite. Connectivity gives ONE wall c closing all three.

Individually tau_i<=sqrt(e), so EVERY actual tangent square satisfies

    tanh²(L_i)<=4e/(1+e)²<=2600/3969=:r.              (2)

This is a new bound at the changed target, not the older16/25 value.

Write the UNIT wall as c=A p-v_A-v_S, h=A²>0, and
v_A=sum beta_i omega_i, beta_i=-(c,n_i)>=0. Let

    q=beta^T C beta, Q=||v_S||²>=0,
    h=q+Q-1, v_i=(C beta)_i/sqrt(d_i).

Each v_i is positive because c closes its actual edge. Equations (2)
and tanh(L_i)=sqrt(h)/v_i give

    q-1<=h<=r*min_i v_i².                            (3)

Every positive contacted whole spectator costs at least1/4 and is
STRICTLY below h, by the actual orthogonal-opposite theorem. Thus
h<1/4 forces Q0. All spectators, rather than chosen contacts, enter Q.

## 3. Complete contact coverage and the two-contact contradiction

If adjacent beta_i,beta_j are positive, the ACTUAL finite endpoint
of the third edge contains c,n_i,n_j. Its white triangle already has
three contacts at least1/2 and is nonpositive, affine at the minimum.
This contradicts the positive finite endpoint Gram. Hence two positive
contacts can occur only at endpoints1 and3. There cannot be three.
There is at least one contact because all v_i are positive.

Suppose both endpoint contacts u=beta1 and w=beta3 are positive.
If both are1/2, the actual principal c,n1,n2,n3 is affine A3:
its positive root sum is nonzero future and feasible against every
other original wall by acuteness. This contradicts compactness.
In every other case native discreteness gives u,w>=1/2, with one
at least1/sqrt2. The exact own cost and middle derivative are

    q=(3/2)(u²+w²)+u w, v2²=(u+w)²/2.

Therefore

    q-r*v2²=((3-r)/2)(u²+w²)+(1-r)u w
      >=3(3-r)/8+(1-r)sqrt2/4
      =(27921+2738sqrt2)/31752 >1.                    (4)

For the last strict comparison sqrt2>7/5 gives the stronger lower
bound158771/158760=1+11/158760. On the other hand (3) implies
q-r*v2²<=1. This excludes EVERY two-contact alternative, including
larger finite, parallel and dotted contacts. No arithmetic restriction
or upper native label cutoff was used.

## 4. Every single-contact case is impossible

For a single endpoint contact beta, direct multiplication of C gives

    q=(3/2)beta², min_i v_i²=q/9.

From (3),

    q<=35721/33121,
    beta²<=23814/33121<3/4,
    h<=2600/33121<1/4.

The complete actual own-contact possibilities below3/4 are native3/4/5.
Each has q<1; for native five this follows from
cos²(pi/5)=(3+sqrt5)/8<2/3, using sqrt5<7/3. Small-height
isolation gives Q0, so h=q-1<0, contradiction. Native6 and every
larger finite/parallel/dotted value were excluded by beta²<3/4.

For a single middle contact beta,

    q=2beta², min_i v_i²=q/3.

Equation (3) gives

    q<=11907/9307,
    beta²<=11907/18614<16/25<cos²(pi/5).

Here 25*11907=297675<297824=16*18614; the last strict comparison
follows from sqrt5>11/5. Only native3 or4 can remain, and hence
q<=1. Using THIS finite-label conclusion in (3),

    h<=r/3=2600/11907<1/4.

Small-height isolation again gives Q0 and h=q-1<=0, contradiction.
The generic earlier bound2600/9307 is NOT claimed below1/4; the
native restriction is established before taking the stronger height
bound. Thus no foreign case is lost through a reversed inequality.

The complete zero/one/two-contact alternatives have been exhausted.
This contradicts e<=13/50 and proves the unconditional compact A3
component theorem with arbitrary spectators.

## 5. Exact arithmetic16 consumer and limits

For A11+O7+A3 at the same original farthest vertex, the accepted
arithmetic A11 theorem, actual typed-defect floor and the new universal
A3 theorem are disjoint whole-component estimates. They give

    ||x||²=e_A11+e_O7+e_A3
      >13/20+9/100+13/50=1,

contradicting ||x||²<1. The ACTUAL old-seven/new-three condition is
retained; no bare I2(7) component is substituted and no defect energy
is counted twice. This closes original group5 of the fresh arithmetic16
list only after full outside-author/root acceptance.

The component theorem itself is nonarithmetic and rank-free. No other
whole forest, general compact ceiling, arithmetic ceiling fifteen,
new example, worker execution, Lean theorem or manuscript revision is
claimed by this candidate.

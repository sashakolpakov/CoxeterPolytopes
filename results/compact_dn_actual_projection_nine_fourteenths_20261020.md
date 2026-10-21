# Uniform actual D projection above 9/14 and the five-rank spectator frontier

20 October 2026. Complete paper candidate, awaiting full independent
review. No scientific execution, optimization, arithmeticity, Lean or
TeX. Accepted sources are unchanged. The native7 input below is a
fixed paper proof with completed root and outside-author reviews;
this assembly still awaits its own full review.

**Theorem.** At a farthest vertex p of any compact Coxeter polytope,
the actual common spatial vector x satisfies

    ||x_D||^2>9/14

for every WHOLE D_r component, r>=9, with arbitrary finite spectator
components. This is an ACTUAL farthest-projection bound. It does not
claim that every vector satisfying only edge half-tangent floors has
this norm, or that an edge-only covering energy H_K exceeds9/14.

**Forest consequence.** In dimension n>=15, if the farthest diagram
has no A or B component of rank at least five, then it has exactly one
large D component D_r, with

    r>=max(10,n-5),    spectator total rank<=5.

This strengthens a structural frontier; it does not give a new
unconditional dimension ceiling or exclude these remaining forests.

## 1. Pinned actual-wall inputs and noncircularity

All paths in this paragraph are under `results/`.

* Complete common1/2 reduction:
  `compact_dn_common_one_two_native_seven_reduction_20261020.md`, SHA256
  `ccea4a223a5c0dafc3f83166d19a4d2f5f1f18ae622860b3d2bc5dafd6eac95c`.
* Its remaining native7 case:
  `compact_dn_native_seven_actual_dual_foot_energy_20261020.md`, SHA256
  `a39306eb2637acf75959530fb5cf5d37dd77d28efdcf6cf2060708a63620d769`.
  This constructs an ACTUAL future timelike point as a positive
  combination of the known wall and old nodes1,2. Its farthest
  pressure, rather than a hypothetical new opposite, closes the case.
* Complete common2/3 inventory:
  `compact_dn_low_index_common_pair_inventory_and_relaxation_20261020.md`,
  SHA256 `2bd35ff5bd18599c5f4624128e055c3669373eea88d4b9edf8cb7f830c9a99f3`.
* Internal common pairs:
  `compact_dn_internal_consecutive_common_opposite_20261020.md`, SHA256
  `ff1aaa687171eec9ebf2669746bdc37cd753dcbcef62890ebf49be20a791bd11`;
  the r9 center/neighbor alternative is
  `compact_dn_center_neighbor_common_opposite_energy_20261020.md`, SHA256
  `e371611c4903ad362654ee0d4656457140cdc5e4374b4e19a2697718c2d679b1`.
* Actual farthest/shadow conventions:
  `compact_farthest_vertex_energy_obstruction_20261020.md`, SHA256
  `50c64e5a521d3e4e0f277a8ece02d3a29aa4ae8b71485cf3088211e587f45bff`,
  and `compact_actual_wall_reciprocal_shadow_20261020.md`, SHA256
  `c913ec3bd8e0f26a8c679a1ede2b98b251ad3f36efeae60486213ddfa3e8badc`.

The local low-index arguments use only independently established
common4/5 or6/7 alternatives and actual feasible-point comparisons;
none assumes this new uniform9/14 conclusion. Thus the assembly
below introduces no circular component estimate.

## 2. Three distinct actual pairs

Assume `||x_D||^2<=e=9/14`. Number the old D chain1,...,r-2,
fork leaves r-1,r. For i<=r-2 write

    u_i=(e1+...+ei)/sqrt(i), h_i=x_D.u_i,
    F(h)=2h/(1+h^2).

These have 0<h_i<1. The pinned complete common1/2 reduction and
native7 theorem make those two opposites distinct. The common2/3
theorem makes that pair distinct. For6/7 a common opposite forces
an actual squared half tangent above3/4: use the center/neighbor
theorem at r9, and the internal theorem at i6 for r>=10. Thus those
opposites are distinct as well.

The actual distinct-opposite shadow law gives, respectively,

    F(h1)F(h2)>=sqrt(1/2),
    F(h2)F(h3)>=sqrt(2/3),
    F(h6)F(h7)>=sqrt(6/7).                          (1)

All ranks use precisely the same first seven coordinate directions.
No cap-product label or spectator profile is omitted by imposing a
desired finite continuation; (1) applies to the actual opposite walls.

## 3. Exact lower bounds for the three pair sums

On (0,1), F is increasing and log(F) is strictly concave, since

    (log F)''=-1/h^2-2(1-h^2)/(1+h^2)^2<0.

Therefore if h+h'<=s then
`F(h)F(h')<=F(s/2)^2=16s^2/(4+s^2)^2`.
The following rational seals show that this maximum is below the
corresponding correlation in (1), giving the strict floors listed.

| Pair | Strict sum floor s=m/1000 | Correlation lower bound p/d |
| --- | ---: | ---: |
| 1/2 | 1091/1000 | 7071/10000 |
| 2/3 | 1265/1000 | 81649/100000 |
| 6/7 | 1512/1000 | 9258/10000 |

The squared comparisons placing p/d below the exact correlations
have positive cleared differences959,322397,26052, respectively:
`10000^2/2-7071^2`, `2*100000^2-3*81649^2`, and
`6*10000^2-7*9258^2`.
For the other side, the positive integers

    p*(4000000+m^2)^2-16*d*m^2*1000000

are respectively

    40828209713231,
    358399613480625,
    52443946917888.

This proves all three pair-sum floors by rational cross multiplication;
no approximate minimization or numerical root computation is used.

## 4. A fixed positive dual vector contradicts Cauchy

Put w12=u1+u2, w23=u2+u3, w67=u6+u7 and choose

    W=(1/10)w12+(11/200)w23+(61/200)w67.

Every coefficient is positive, so the preceding floors give

    x_D.W > .1*1.091+.055*1.265+.305*1.512=.639835. (2)

Here and below every terminating decimal is an exact rational. The
Gram entries have the following upper bounds:

| Entry | Exact expression | Strict upper |
| --- | --- | ---: |
| w12^2 | 2+sqrt2 | 3.415 |
| w23^2 | 2+2sqrt(2/3) | 3.633 |
| w67^2 | 2+2sqrt(6/7) | 3.852 |
| w12.w23 | 1+1/sqrt2+1/sqrt3+sqrt(2/3) | 3.102 |
| w12.w67 | 1/sqrt6+1/sqrt7+1/sqrt3+sqrt(2/7) | 1.899 |
| w23.w67 | 1/sqrt3+sqrt(2/7)+1/sqrt2+sqrt(3/7) | 2.474 |

They follow by positive squares from
`sqrt2<1.415`, `sqrt(2/3)<.8165`, `sqrt(6/7)<.926`,
`1/sqrt2<.7072`, `1/sqrt3<.5774`, `1/sqrt6<.4083`,
`1/sqrt7<.378`, `sqrt(2/7)<.5346`, `sqrt(3/7)<.6547`.
Consequently

    W^2 < .01*3.415+.003025*3.633+.093025*3.852
          +.011*3.102+.061*1.899+.03355*2.474
        = .636435825.                              (3)

The exact rational Cauchy margin is

    .639835^2-(9/14)*.636435825
        =70423123/280000000000>0.                  (4)

Equations (2)--(4) contradict
`(x_D.W)^2<=||x_D||^2 W^2<=e W^2`. This proves the universal actual
projection theorem. In particular the negative-core actual-point
inputs have not been replaced by edge-only energy assumptions.

## 5. Conditional forest frontier

Use the accepted arbitrary-spectator accounting in
`compact_dn_uniform_energy_and_chain_relaxation_20261020.md`, SHA256
`1bd7c4eb372083973ca059b34121cdf0e3f575945a1552d771d7b09c84617804`.
For n>=15 and no A/B rank>=5, it already supplies a unique D_r with
r>=9, and every spectator component is ordinary with cost at least
rank/14, except for the following mutually controlled defects:
at most one short A1, or at most one isolated rank-two7 defect of
cost greater than9/100. The rank-two defects are OLD A2/NEW7 and
OLD I2(7)/NEW3; they are different from the long-end7 just excluded.
They cannot coexist with a short A1, by the accepted orthogonal
support/native product-gap argument. Other covered components retain
their full rates; no spectator rank is discarded except that single
short A1 when present.

Since the actual orthogonal projections sum to ||x||^2<1, the new
D bound gives the complete table, writing S for spectator total rank:

* No defect: S>=5 gives total greater than9/14+5/14=1; hence S<=4.
* One short A1: S>=6 leaves ordinary rank at least5 and gives the
  same contradiction; hence S<=5.
* One rank-two7 defect: S>=6 leaves ordinary rank at least4, giving
  total greater than9/14+4/14+9/100=713/700>1; hence S<=5.

Thus always S<=5 and r=n-S>=n-5>=10. The older4/7 bound still gives
uniqueness of the large D block, so the stated frontier follows.
This is an actual-vector aggregation, not a claim that the new bound
is a purely local edge-covering rate.

## Scope and contributions

The three-pair dual weights were supplied by root. Finite57 independently
derived the exact Gram/sum seals and rational margin and assembled this
note. The decisive native7 actual dual foot was supplied by finite57
and independently reconstructed and written by root/guard; its source
and all earlier profile inputs are pinned above. This is a collaborative
proof candidate, not an outside-author review. A new uniform ambient
dimension bound remains open.

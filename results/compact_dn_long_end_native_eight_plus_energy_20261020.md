# D long-end native labels at least eight: actual projection above 9/14

20 October 2026. Paper only, independently derived with root. No
scientific execution, lattice classification, Lean or TeX job. No
accepted-source edits or ambient dimension claim.

**Theorem.** Let p be a farthest vertex from a strict interior point
of a compact Coxeter polytope. Suppose its finite diagram has a WHOLE
component D_r, r>=9, with arbitrary orthogonal finite spectators.
Number the long chain1,...,r-2 and the fork leaves r-1,r. If an actual
wall c closes directions1 and2, has just one nonzero D contact at
node1, and that native label is a finite m>=8, then the actual
farthest projection satisfies ||x_D||^2>9/14.

This is a complete conditional branch exclusion. It does not classify
all common1/2 opposites, exclude native7, or establish a rate for every
D component.

## 1. Spectators are forced to vanish

Use unit old roots. For long-arm nodes i,j<=r-2,

    omega_i=sqrt(2)*(e1+...+ei), C_ij=2min(i,j),
    u_i=(e1+...+ei)/sqrt(i).

Put beta=cos(pi/m), theta=pi/m, and let Q be c's total spectator
inverse-Gram cost. At the actual endpoint deleting2, old node1 is an
isolated A1 in the retained D part. Positive definiteness therefore
requires beta^2+Q<1. Any contacted spectator contributes at least1/4:
its native contact square is at least1/4 and the relevant inverse
unit-Gram diagonal is at least1. Since beta^2>=cos^2(pi/8)>3/4,
Q<1/4, and consequently Q=0.

Thus c=a*p-v with v=beta*sqrt(2)*e1 and a^2=2beta^2-1. Its actual
direction2 full tangent is a/beta. The half-tangent floor is exactly

    tau2^2=(beta-sin(theta))/(beta+sin(theta))
          =tan(pi/4-theta)>=sqrt(2)-1.                (1)

In particular h2=x_D.u2>=sqrt(sqrt(2)-1). This calculation does not
identify the full tangent with the half tangent.

## 2. Directions6 and7 require a large sum

Assume for contradiction ||x_D||^2<=e=9/14. Their actual opposites
must be distinct. For r=9 the common6/7 theorem is the center/long-
neighbor theorem; for r>=10 it is the internal consecutive theorem
at i=6. In either case a common wall forces an actual squared half
tangent greater than3/4, incompatible with the assumed projection.

The correlation of u6,u7 is k=sqrt(6/7). The distinct-opposite shadow
law and F(h)=2h/(1+h^2) give

    F(h6)*F(h7)>=k,
    h6*h7>=p(k)=k/(1+sqrt(1-k))^2>9/16.

Here the product bound follows from h6^2+h7^2>=2h6*h7. Its strict
last comparison is equivalent to k>576/625, whose squared comparison
has numerator6*390625-7*331776=21318>0. Hence

    x_D.(u6+u7)=h6+h7>3/2.                          (2)

## 3. Orthogonal Cauchy contradicts that sum

Let w=u6+u7 and write w=rho*u2+w_perp. Directly from the prefix
directions,

    rho=1/sqrt(3)+sqrt(2/7)<1.114,
    S=||w_perp||^2=29/21+4sqrt(2/21)<2.62.

For t=h2, orthogonal Cauchy gives

    x_D.w<=f(t)=rho*t+sqrt(S)*sqrt(e-t^2).

This function is decreasing for t>=q=sqrt(sqrt(2)-1). Indeed its
derivative is negative precisely when rho^2*e<||w||^2*t^2. The
left side is less than(5/4)*(13/20)=13/16, whereas the right side
is greater than3*(2/5)=6/5. Thus (1) permits evaluation at q.

For explicit rational seals, q<.644 and e-q^2<.23. Also
sqrt(2.62*.23)<.777. Therefore

    x_D.w<=f(q)<1.114*.644+.777
                  =1.494416<3/2,

contrary to (2). The radical bounds are elementary positive squares:
1/sqrt(3)<.578, sqrt(2/7)<.535, sqrt(2/21)<.309;
1.414<sqrt(2)<1.4143; and .777^2>.6026.

## Pinned actual-geometry inputs and boundary

The farthest projection inequality and its strict norm convention are
`compact_farthest_vertex_energy_obstruction_20261020.md`, SHA256
`50c64e5a521d3e4e0f277a8ece02d3a29aa4ae8b71485cf3088211e587f45bff`.
The shadow input is `compact_actual_wall_reciprocal_shadow_20261020.md`,
SHA256 `c913ec3bd8e0f26a8c679a1ede2b98b251ad3f36efeae60486213ddfa3e8badc`.
The complete arbitrary-spectator common-opposite inputs are
`compact_dn_internal_consecutive_common_opposite_20261020.md`, SHA256
`ff1aaa687171eec9ebf2669746bdc37cd753dcbcef62890ebf49be20a791bd11`,
and `compact_dn_center_neighbor_common_opposite_energy_20261020.md`,
SHA256 `e371611c4903ad362654ee0d4656457140cdc5e4374b4e19a2697718c2d679b1`.

Every closure assertion above concerns an actual opposite. No formal
finite completion is promoted to a wall, no spectator is assumed
harmless, and no nonreflective subspace or published literal-family
cutoff is imported. Native7 remains the explicit open branch here.

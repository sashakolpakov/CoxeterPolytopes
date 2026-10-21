# Independent audit of the cubic-seven dihedral hosts

The independent replay passes for all eleven saved searches: the H2 calibration is the known compact triangle with angles π/2, π/3, π/14; all ten H3, H4, and H9 prefixes have exact infinite-volume certificates. The extra Pell cap is a valid compatible reflection wall beyond the height cutoff, but the saved H9 escape survives it. This produces no new finite-volume polytope and proves no nonreflectivity statement about a full lattice.

Replay with standard-library Python only:

```sh
python scripts/fresh_cubic7_dihedral_independent_audit_20261013.py --verify
```

The [verifier](../scripts/fresh_cubic7_dihedral_independent_audit_20261013.py) imports no project arithmetic or search code. Its [JSON record](fresh_cubic7_dihedral_independent_audit_20261013.json) includes every case, exact strictly spacelike witnesses, the binary root census, and the cap check. Signs use rational isolating intervals for the cubic's three real roots. Forged level heights and timelike replacements of escape witnesses are rejected.

| Host | Height | Walls | Independent result |
|---|---:|---:|---|
| H2, α=2τ, one block, no scalar | 10 | 3 | Compact known (2,3,14) triangle |
| H3, α=τ / 2τ / τβ | 30 | 6 / 4 / 6 | Exact escapes |
| H4, α=τ / 2τ / τβ | 30 | 4 / 8 / 4 | Exact escapes |
| H9, α=2τ | 10 | 13 | Exact escape |
| H9, α=τ / 2τ / τβ | 30 | 18 / 25 / 18 | Exact escapes |

## Norm alphabet

Let τ³+τ²−2τ−1=0, with τ=2cos(2π/7) at the distinguished real place. Set C=τ²+τ−1, β=2−τ and D=2−τ². Exact identities are C²=2+τ, D²=2−C, N(β)=N(2+C)=7, N(2−C)=1 and N(C)=−1. The block B=[[2,−C],[−C,2]] has determinant β and Smith factors (1,β), since C is a unit. The host is ⟨−α⟩⊕B^k⊕⟨1⟩, omitting the last summand in even hyperbolic dimensions.

The field has integer ring Z[τ], discriminant 49 and class number one. Indeed a possible index 7 would give discriminant 1, excluded by Minkowski; the ideal-class Minkowski bound is 14/9<2. The signatures of −1, τ and C generate all eight sign patterns. Since a real cubic field has eight unit square classes, every totally positive unit is a square. Also 2 is inert: X³+X²+1 is irreducible over F₂.

For a primitive reflection vector r of norm q, integrality says 2Gr/q is integral. If eG⁻¹ is integral, then (2e/q)r is integral; primitivity forces q to divide 2e. Here e=β for α=τ or τβ, giving normalized positive norm representatives {1,2,β,2β}. For α=2τ, e=2β gives the safe producer pool {1,2,4,β,2β,4β}.

The latter pool can be sharpened for primitive vectors: norms 4 and 4β are impossible. Reflection integrality first forces every spatial coordinate to be divisible by 2, because the spatial Gram matrix is invertible modulo 2. The norm equation then forces the time coordinate to be divisible by 2 as well. The producer's larger pool is conservative and removes parallel duplicates; it does not miss primitive mirrors.

## Complete zero-height chamber

Write F(x,y)=x²−Cxy+y². For a nonzero integral pair, F is a totally positive algebraic integer, so Tr(F)≥3. This follows from its positive integral norm and the arithmetic–geometric mean inequality. Apply this separately to every block of a zero-height vector.

For q=1 the equation 2ΣF+z²=1 forces all block coordinates to vanish, leaving z=±1. For q=2, reduction modulo the inert prime 2 makes z even; after division by 2, the trace inequality forces exactly one nonzero block, with F=1, and z=0.

For q=β or 2β, reflection integrality forces Br to be divisible by β in each block and z=βs. Write r=adj(B)t. Dividing the norm equation by β gives 2ΣF₊(t)+βs²=1 or 2, where F₊(t)=t₁²+Ct₁t₂+t₂². In the first case all blocks vanish and βs²=1 is impossible by field norms. In the second, s is even; dividing by 2 and using traces leaves precisely one F₊=1 block and s=0. The conservative q=4 roots are only ±2 times the scalar root; q=4β gives none.

The verifier independently enumerates every binary F=1 solution. Completing the square gives x²,y²≤4/β at all embeddings. Trace-dual coefficient bounds enclose each coordinate in |a|≤7, |b|,|c|≤4. Exact enumeration finds fourteen solutions, exactly the orbit of (1,0) under (x,y)↦(Cx−y,x). The β-family gives another fourteen oriented roots, the orbit of D(−1,1). Thus H9 has 114 oriented primitive zero-height roots, or 57 mirrors, and the full chamber is I₂(14)^4×A₁. The dominant census has nine roots, with an additional scalar duplicate only in the conservative norm-4 pool.

The two block mirrors e₁ and D(e₂−e₁) have 4cos²(π/14)=2+C, and reflection in the second literally exchanges the coordinates. This angle cannot occur between roots of normalized norm classes 1 or 2 in the older unimodular cubic hosts: taking field norms in p²=q₁q₂(2+C)/4 gives an odd 7-adic valuation on the right and an even one on the left. This distinguishes mirror families, not quadratic spaces or commensurability classes.

## Geometry and search scope

The replay checks every saved root's norm, integrality and reflection integrality; every wall pair's sign and finite Coxeter angle; the exact control heights and their chronology; the complete zero-height census; and every strict timelike center. Positive-height candidate enumeration has been reviewed but is not independently rerun by this verifier. The producer bounds dual block coordinates by wᵢ²≤2(q+αx₀²) at every embedding, reconstructs integral coordinates using B⁻¹, and enumerates all sorted block multisets and their permutations. Its equal-height packets are checked for mutual acuteness.

For the triangle, all three affine vertices are timelike, and the three spatial normals have a strictly positive dependence with coefficients 2+4τ+2τ², 2+5τ+2τ², 1. Their span has full rank, so the affine section is bounded; the vertex census therefore proves compactness.

Every higher-dimensional escape has positive time coordinate, positive Lorentz square and all wall pairings nonpositive. The verifier adds a small exact positive multiple of the strict center until the witness has strictly negative pairings while remaining spacelike. The segment from the center to that point crosses the future light cone strictly inside all inequalities, yielding an open ideal-boundary patch and infinite volume for the saved prefix.

All four saved H9 escapes vanish in the scalar coordinate and two spatial blocks. This suggests targeted searches on the two active blocks; the saved H4 calibrations test that mechanism at the same cutoff. It is not a proof that future roots cannot cut those directions. In particular, no positive-cone projection obstruction can use a slice containing the spacelike witness itself.

The extra norm-one cap has time X=6+18τ+8τ² and scalar Z=11+28τ+12τ², satisfying Z²−2τX²=1. Its height is 648+1816τ+808τ². All 25 original pairings are nonpositive, all finite added angles are right angles, and it contains the old strict center. Its pairing with the saved H9 escape is −16−44τ−20τ²<0, so the escape remains. This explicitly demonstrates that the height-30 prefix is not the full collection of available reflection walls.

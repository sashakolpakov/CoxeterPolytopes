# Norm79 A9 bridges: an exact integral reversal gate and a conditional next edge

21 October 2026. **New manual supporting candidate; not yet independently reviewed.** No worker, import, enumeration, symbolic calculation, Lean, TeX or test ran. No frozen source, publication, guide or checkpoint was changed. Neither bridge has been observed, and the interval census remains unexecuted.

**Result.** At either possible full actual A9 bridge germ, its nontrivial diagram reversal preserves the ORIGINAL lattice if and only if both bridge controls are divisible by `lambda=sqrt5=2phi-1`. In that branch the reversal gives a complete actual transverse edge with a proved finite upper and full finite endpoint germ. In the other branch the reversal is not an integral isometry. This is not a bridge exclusion, a chamber period, finite-volume completion, new example or form-class statement.

## 1. Exact inherited boundary

All paths below are under `results/`; the displayed hashes match current bytes.

| Input | SHA256 |
|---|---|
| compact_h9_prime79_seven_model_continuation_candidate_20261021.md | `343e8422d2705a9dec3a82fea3658367213ffd2d978fc6346604de7d3130348c` |
| compact_h9_prime79_a7_c_finite_upper_and_complete_interval_contract_20261021.md | `b4de91853f6fd625831d591554cb47f7a61400fde49c233017ef905e8b24c3fa` |
| compact_h9_prime79_next_edge_geometric_continuation_candidate_20261021.md | `cabea5bfbaa7d45ecd8fe461000c18f67ee6c601353dfa89376010ce26f6bfca` |
| review_a16_lower_a_arithmetic17_and_prime79_root_20261021.md | `9f1dd75a06b6826821a6aa186c0b8183ec782d2bba00e9c1d0a31d80957ae272` |
| review_prime79_seven_model_continuation_full_compact_reviewer_20261021.md | `9dd4c9f09239d5d0aea6dc935797cacdaa455d81adfb258203ae56fa729f0e70` |

Use `O=Z[phi]`, `phi²=phi+1`, `a=3phi²`, `D=18+31phi`, and the literal module

    L=E8(O) orthogonal-sum [[-2phi,a],[a,2]].

Its original binary basis is `Bx,By`; put `c=-By`, `p=2Bx+a*c` and `Q=omega1`. Thus `p²=-2D`, `Q²=8`, and `p,c,Q` are pairwise orthogonal. The native E8 edges are `02,13,23,34,45,56,67`; `Q.e1=1`, all other `Q.ei=0`. The full integral-arrangement chamber has the inherited actual corner `p` of germ E8+A1(c). Release `t=-e1` and retain the A7 path `(-e0,-e2,-e3,-e4,-e5,-e6,-e7)` and `c`.

A bridge normal is

    s=(U/2)p-c/2+omega_j+YQ,   j=0 or7,
    z=8Y+I_j,   I_0=5, I_7=3,
    z²-4DU²=5,   U=phi mod2,   U>0, Y>2.

The original-module integrality and all first-simple coverage hypotheses are inherited, not reproved by this equation. The corresponding finite ray is

    P=z p+DUQ,       P²=-10D,       p.P=-2Dz<0.             (1)

For the geometric conclusions below, **assume this is the FULL ACTUAL first A9 germ**, with `s` an actual simple wall. A norm-equation solution or an unverified minimum does not supply this premise.

## 2. Exact integrality criterion for the corner point symmetry

Let `r_P(w)=w-2(w.P)P/P²`; this is reflection in a TIMELIKE vector, not a hyperbolic wall reflection. Its negative `J=-r_P` is the geodesic point symmetry fixing the ray of `P`. We claim

    J(L)=L  iff  U,z belong to lambda O,  lambda=2phi-1.     (2)

Here `lambda²=5` and `O/(lambda)=F5`, with `phi=3` and `D=1`. Reducing the bridge equation gives `z²=4U²`, so U and z are either both zero or both nonzero modulo lambda. If they are nonzero, then

    r_P(e1)=e1+(U/5)P

has Bx coefficient `2Uz/5`, which is not integral at the prime `(lambda)`. Neither r_P nor J preserves L.

Conversely write `U=lambda u`, `z=lambda z0` with `u,z0 in O`. Then

    q=P/lambda=z0p+DuQ,   q²=-2D,
    z0²-4Du²=1.                                         (3)

On the ORIGINAL basis, using `Bx.p=-D` and `By.p=0`,

    r_q(e1)=e1+u q,      r_q(ei)=ei  (i!=1),
    r_q(Bx)=Bx-z0 q,     r_q(By)=By.                     (4)

Every image is integral. As an involution, r_q has integral inverse, proving equality of modules. Consequently J is an integral isometry, is future because it fixes the future timelike q, and acts as minus the identity on `q`-perp. This proves (2) at both bridge endpoints without a height restriction. The timelike reflection itself is not subject to the primitive POSITIVE mirror alphabet and is not an additional arrangement wall.

## 3. The unique A9 germ reversal

Order the nine inward norm2 simple roots as follows:

    j0: (beta0,...,beta8)=(c,s,-e0,-e2,-e3,-e4,-e5,-e6,-e7),
    j7: (beta0,...,beta8)=(-e0,-e2,-e3,-e4,-e5,-e6,-e7,s,c).

Let w0 be the longest A9 Weyl element. In the standard A9 realization it reverses ten coordinates, hence `w0 beta_i=-beta_(8-i)` and `w0²=1`. It is a finite word in the nine norm2 integral reflections, so preserves the original L, fixes P, and is future. Define

    T=J w0.

Then T fixes P, sends `beta_i` to `beta_(8-i)`, and is an involution. Since w0 is an integral automorphism, T is integral exactly under condition (2). Any future isometry fixing this corner ray and reversing its inward A9 chain is T: the nine roots and P span the whole ten-dimensional space. Thus the nondivisible branch cannot be repaired by choosing another extension of the same diagram reversal.

Under the FULL ACTUAL germ premise, integral arrangement invariance and the permutation of every inward incident facet show `T(C)=C`: the two chambers have the same germ at the finite vertex and hence are the same arrangement chamber. No arbitrary reflection subgroup is substituted for the full arrangement. Within the K-linear group `O^+(L)`, the chamber stabilizer fixing this vertex is therefore of order2 in the divisible branch and trivial otherwise: the only A9 diagram automorphisms are identity and reversal, and fixing all nine normals and the future ray determines the identity.

## 4. A complete actual transverse edge, not another guessed upper

Suppose (2) and the full actual germ hold. The inherited complete old edge joins p to q. Its eight retained roots are all beta_i except s. Applying T gives a complete edge of the SAME chamber:

| Bridge | Released wall at q on the new edge | Other endpoint | Full endpoint germ |
|---|---|---|---|
| j0 | `-e6=beta7` | `T(p)` | `T(E8+A1(c))` |
| j7 | `-e2=beta1` | `T(p)` | `T(E8+A1(c))` |

The actual wall `T(t)` is a finite upper for this transverse edge. Indeed `(q,Tt)=(q,t)=-Du<0`, `(Tp,Tt)=(p,t)=0`, and `(Tp)²=-2D`. Firstness, absence of intervening walls and the complete endpoint germ follow by the isometry of the old COMPLETE actual edge, not from finite Gram positivity. The original backward wall t, ultraparallel to s, is not being reused as an upper.

For an exact distance check, let Lambda_i be the A9 duals and k the position of s (1 for j0,7 for j7). Since `p.beta_i=-DU` only at i=k,

    p=z0q-DU Lambda_k,      Tp=z0q-DU Lambda_(8-k).

The inverse A9 entries are `Lambda_k²=8/5` and `Lambda_k.Lambda_(8-k)=2/5`. Thus

    cosh d(p,Tp)=-(p,Tp)/(2D)=1+3Du²>1.                 (5)

These are distinct finite endpoints. No trace or infinite-order conclusion follows just from this distance.

## 5. Remaining bottleneck and execution boundary

The explicit finite upper in b4de has `z*=lambda Z0`, `U*=lambda W0/2`, with even W0. It lies in the divisible branch of (2). If it or any other divisible bridge is the full first A9 germ, Section4 now supplies its next complete actual edge without a second interval census.

The bridge norm equation does NOT force condition (2): modulo lambda it also permits the two nonzero branches `z=+/-2U`. No global O-solution, actual minimum or germ is asserted in those branches. Conversely the failure of this particular integral reversal does not exclude a genuine asymmetric chamber or every possible chamber period. Even in the divisible branch, combining T with the inherited base involution still needs a proved infinite-order argument before any noncofiniteness claim.

The unresolved tasks are actual firstness/all ties of the two bridges, their divisible or nondivisible controls, and an unconditional continuation or period. No observed census result, finite-volume polytope, full host exclusion, field-wide theorem, commensurability-class conclusion or change in a published/accepted dimension bound is supplied here.

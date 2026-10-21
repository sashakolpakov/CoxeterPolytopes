# Exact low-Z8 gate and complete fixed-prefix low-RZ union

Status: FINAL/FROZEN. Generation14783 and fresh saved replay99542 terminal PASS. Every one of the392 low-Z8 norm-frame recipes is excluded on both exact branches. Together with the frozen [continuous high-Z8 gate](fresh_joint_rz_high_z8_gate_20261017.md), this closes all finite Z8 and ZX labels in the specified **RZ≤5** polar-Y incidence. It does not exclude high RZ≥6 or other prefix geometry.

The fixed prefix is C rightS6+7+30 with C10=9, A rightS6+7+8 with AC=3, and X with labels8/A/C=2/3/2; all three have zero A2 profile. R is zero-profile, rightS6+8+63, with native R44 label2..15. The latter is the complete necessary circle domain from the [native circle theorem](fresh_r_circle_admissibility_20261017.md). The exact construction uses four three-contact/norm recipes perR, verifies positive norm remainders, and checks R·10<0 explicitly to bind the correct permitted sheet. The numerical seed selects only an algebraic branch sign.

Z must be an actual side of the target S6 face and have the three finite corners8/X/Z, Y/8/Z and Y/R/Z. H4 finite-extension isolation and the inherited A2 extension classifier give exactly seven Z profiles. Y is stipulated to have norm2, zero A2 profile, and be right8/R/Z. In the low range Z8,ZR,ZX∈2..5, exhaustive finite-corner filtering leaves126 recipes. Since8 andR are orthogonal norm-two vectors, the exact Schur remainder of the projected8/R/Z Gram is

`qZ − (theta(Z8)^2 + theta(ZR)^2)/2`.

The signature table is independent ofR:96 remainders are strictly positive,2 are exactly zero, and28 are strictly negative. The positive cases have a positive-definite3-space with a timelike orthogonal line, so cannot contain a spacelike Y in its complement.

Both zero cases have Z profile0, Z8=ZR4 and ZX2 or3. Their literal Gram is

```
[ 2       0       −sqrt2 ]
[ 0       2       −sqrt2 ]
[ −sqrt2  −sqrt2   2      ]
```

Its positive kernel is(1,1,sqrt2). Let ν=8+R+sqrt2 Z. In any coherent candidate chamber, a strict common intrinsic center pairs negatively withν, soν is nonzero. A linear dependence of the three roots would have to lie in their one-dimensional Gram kernel, contradictingν≠0. Thus the roots span a three-dimensional degenerate subspace with null radical Rν. Its orthogonal complement in Lorentz4 is exactly that null line, containing no norm2Y. This handles the28 zero-limit instances across14 R values; it does not assume ideal vertices are simple.

For the28 Lorentzian recipes at each R, the checker solves Z in the frame8/R/X, certifying the positive complementary line before applying the norm formula. Both square-root branches are retained. Every native product restores the A2 offset against the retained old wall, including old walls whose own projections require the fullS6 correction. Every branch has a strict old-native obtuse or consecutive-angle-gap witness.

| First rejecting wall | Strictly obtuse | Strict consecutive-angle gap |
|---|---:|---:|
|4|392|0|
|7|38|65|
|10|0|116|
|30|66|107|
|Total|496|288|

This is392 exact recipes and784 rejected branches, with no unresolved case. It replaces the numerical evidence for these recipes in the [frozen joint scout](fresh_joint_rz_circle_scout_20261017.md) by a separate exact certificate; that frozen scout is unchanged.

To complete the Z-label domain under RZ≤5, a high Z8≥6 in the finite8/X/Z corner forces zero profile and ZX right. The Lorentzian polar-frame condition leaves RZ3,4,5, giving the42 entire high-Z8 families already excluded by the frozen continuous gate. A high ZX≥6 instead forces Z8 right and zero profile; the8/R/Z span is then positive definite for every finite RZ≤5, so it cannot supply this Y. These reductions cover every finite Z8/ZX label without enlarging a numerical cutoff. The lower endpoint of the high tail is relaxed conservatively, including the affine/null limit.

The result is an exact necessary-condition obstruction for this fixed C/A/X prefix, zero-profile R circle, RZ≤5 and zero-profile polar Y. High RZ, changing A/C/X, nonzero-profile R, nonpolar Y and other disk incidences remain outside the theorem. The parent numerical prefixes are not promoted to exact chambers. The intended deletion set12,47,52,64,65 still debits11 old finite-vertex and six compact2-face certificates; nothing here reconstructs those faces. The exact63-wall/199-finite/10-open-end infinite-volume control remains unchanged, with no new finite-volume example claimed.

The finite-diagram and outward-rational physical helpers are inherited. Fresh replay is not an independent physical re-enclosure. Frontier independently checked the positive/affine-frame orthogonal-complement argument, without checking the literal frame census or784 arithmetic witnesses.

A separate Fraction-only certificate walk checked all126 saved signature entries, the14 permitted R-sheet signs, positive exact frame-complement/remainder bounds and all784 saved strict inequalities. This recount does not independently reconstruct the physical intervals. All18 direct source/data manifests, Python syntax and local report links pass.

[Source](../scripts/fresh_joint_rz_low_z8_gate_20261017.py), [data](fresh_joint_rz_low_z8_gate_20261017.json), [high-Z8 certificate](fresh_joint_rz_high_z8_gate_20261017.md), [R-circle certificate](fresh_r_circle_admissibility_20261017.md), [bounded scout](fresh_joint_rz_circle_scout_20261017.md).

Source SHA256 `f5b4d224a99252729b496d452f953a189786fd74d46edc8ab77b087fdf2febda`. Data SHA256 `44def83e6d3b098388ea24f5a9359169bca054037c455ab0c4584eb32cb7a233`.

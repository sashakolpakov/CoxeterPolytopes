# Independent review of the S61 one-cap gate

The [separate checker](../scripts/fresh_h20_s61_single_cap_20261014_review.py) and [certificate](fresh_h20_s61_single_cap_20261014_review.json) verify the frame, complete label bound and all 32 real branches of the [primary gate](fresh_h20_s61_single_cap_20261014.md). Every branch has a positive pairing greater than `99/100` with an old retained wall: sixteen with wall 0 and sixteen with wall 95. No issue was found in the formula, normalization or stated scope.

The check imports no arithmetic from the primary gate or its field helper. It reuses the earlier independent sparse-radical implementation to reconstruct the physical `X,W,r95` frame and its support orthogonality. It checks `X²=−A`, `W²=w>2`, `r95²=2`, the removed weight square `2/3`, and all frame cross products. The two serialized old side normals are compared coefficientwise with `A(−ω−r95/2±W)+√Δ X`, where `Δ=A(w−5/6)`. Thus the normalized projected frame is tied to the actual old walls, not just to an abstract diagram.

After dividing the old side normals by the positive number `A`, their ambient squares are two and their S61 projected squares are `4/3`. A prospective square-two cap perpendicular to S61 has, at a finite endpoint with either side, final Schur pivot `2−3cos²(π/m)`. This is positive exactly for integer labels `m=2,3,4,5`; every `m≥6` gives a pivot at most `−1/4`. No arbitrary label cutoff is used.

The independent calculation solves the norm equation in the **other** coordinate. Write a candidate in the unit Lorentz frame as `(a,b,c)`, put `s=cos(π/m)+cos(π/n)`, `d=cos(π/n)−cos(π/m)`, `v=w−5/6` and `L=v−1/2>0`. The two linear contact equations give

`a=(s−b/√2)/√v`, `c=d/√w`.

Direct symbolic elimination gives

`L b² + √2 s b + v d²/w − s² − 2v = 0`,

and hence

`b = [−s/√2 ± √(v(s²+L(2−d²/w)))]/L`.

The discriminants are strictly positive. The two signs cover every real candidate for each of the sixteen ordered label pairs. The branch correspondence with the primary quadratic is checked, but the root formula and scalar evaluation are independently reconstructed. Symbolic identities establish the roots; interval overlap is only a consistency check.

New outward-rounded dyadic interval arithmetic at 256 bits evaluates the two offending wall pairings, using independently reconstructed physical frame products with wall 0 and the exact identity `r95·N0=−1`. All 32 positive lower bounds exceed `99/100`. The certificate records those strict inequalities without repeating the enormous unreduced rational interval arrays. Generation and saved replay pass; source manifests and links are checked.

A single perpendicular cap making the retained S61 three-sided region compact must meet both old outgoing side walls at finite endpoints. The existing two finite corners survive any coherent acute extension, so the test does not omit a compact quadrilateral obtained by silently removing an old side. Both conic branches are tested even if one has an inadmissible orientation. Excluding this larger candidate set proves the stated one-cap obstruction.

This does not exclude a polygon using several new walls, nonright contacts to S61, or deletion/replacement of retained walls. It does not assert global H20 nonexistence. No frozen predecessor was modified and no old 107-wall pair census was repeated.

# Antipodal cycle16 payments from two finite half-coupons

**Complete local theorem; verification metadata below.** Let Q be an actual quadrilateral of a finite-sided finite-volume hyperbolic Coxeter polytope. Suppose one ideal corner has a varying affine A factor on16 nodes, with the two marked nodes antipodal. Then either some corner is paid by finite marked distance at most9 or ideal marked distance at most7, or Q has two finite corners whose marked roots lie in separate, whole components A1 and A15, with the A15 mark its unique middle node. Assigning one half to each of these two finite corner types pays Q.

Consequently the16-cycle factor can use just its112 pairs of distance at most7, dropping its8 antipodal pairs, provided this finite half-coupon is included. This is a universal local statement, including nonsimple cusps and arbitrarily many supporting whites. It supplies no global finite-vertex budget or ambient dimension bound by itself.

The [source](../scripts/fresh_cycle16_half_coupon_transfer_20261017.py) and [exact certificate](fresh_cycle16_half_coupon_transfer_20261017.json) bind the frozen [mixed finite-nine/ideal-seven lemma](fresh_ideal_distance7_mixed_replacement_20261017.md), its [weak-profile and spectator argument](fresh_distance8_ideal_corner_replacement_20261017.md), and the [actual factor-fiber theorem](fresh_cusp_factor_payment_framework_20261017.md). The existing infinite-volume order16 example is used only as the literal equality core, not as a finite-volume construction.

## Complete antipodal profile reduction

Write the face star as W+a+b+c+d, with opposite pairs {a,c} and {b,d}. Assume no finite-distance9 or ideal-distance7 corner is already paid. Then every adjacent black product is zero, since a nonzero such product supplies a length-one marked path. At the selected ideal corner,

    W+a+b = F orthogonal-sum Z,

where F is the16-node affine cycle, a,b are antipodal, and Z is positive definite. The two old arc lengths are P=Q=8. The factor-fiber theorem makes this decomposition valid even when the ambient cusp has several affine factors: the other factors contribute spherical one-node deletions to Z.

The complete one-root positive-semidefinite extensions of the old A15 path give the same eight profiles as in the frozen weak-profile proof: isolated; a simple endpoint attachment on either arc; a doubled endpoint attachment on either arc; a simple second-node attachment on either arc; and simple attachments at both endpoints. Completeness here applies directly at rank15: connected finite possibilities are A16, B16 and D16, and the connected affine possibility is the16-cycle. Exceptional finite or affine diagrams are too small. A root touching several old nodes cannot produce another finite or affine type containing this induced A15 path.

For each of the six single-attachment profiles, the adjacent corner containing the replacement root and the other old mark has a finite A16/B16/D16 principal component. The opposite endpoint of the old path remains a simple leaf whose neighbor has degree2. Supporting spectators are orthogonal to every old F root, so they cannot change those two degrees. No large connected affine A/B/C/D diagram has that configuration; exceptional affine diagrams are too small. Hence the full component stays finite. The actual factor-fiber theorem then excludes an ideal corner with these two marks in that finite component.

The marked distance is8 for an endpoint attachment and7 for a second-node attachment. Thus each single profile already gives a finite-distance9 payment, regardless of how many spectator roots enlarge the finite component. The checker verifies all six principal positive-definite diagrams, distances and fixed-leaf configurations explicitly. Under the unpaidness assumption, only the isolated and both-end3 profiles remain for c and d.

## The unique determinant equality

Let h=-c.nu and k=-d.nu, where nu is the positive affine null vector with unit coefficients on F. These heights are strictly positive because the opposite sides do not contain the selected ideal point. Put u=1/h and v=1/k. An isolated profile has h>=1 and constant c0=2; a both-end3 profile has h>=2, c0=0 and equal moments1 on the two arcs. The bounds use only the least nonzero acute Coxeter contact magnitude1; they impose no upper cutoff on opposite labels or dotted products.

At P=Q=8 the inherited normalized determinant becomes

    F = 1 - (l_c*u+l_d*v)/2 - (c0*u^2+d0*v^2)/4,

with (l,c0)=(0,2) for isolated and (1,0) for both-end3. Its four profile-pair minima on the allowed height rectangles are respectively

| c profile | d profile | Minimum |
|---|---|---:|
| isolated | isolated | 0 |
| isolated | both-end3 | 1/4 |
| both-end3 | isolated | 1/4 |
| both-end3 | both-end3 | 1/2 |

All16 rectangle values are exact rational checks. The only possible zero is the isolated pair, since

    F = ((1-u^2)+(1-v^2))/2,

and 0<u,v<=1. Therefore h=k=1. The native products are a.c=b.d=-1, and c,d have no other contact with the old cycle. This forces the literal simple16-cycle with a simple leaf at each antipodal mark.

Eliminating arbitrary spectator whites subtracts a positive-semidefinite matrix [[p,q],[q,r]] from the c,d block. The source independently verifies the formal determinant change, after dividing by h^2*k^2:

    Delta = (p/h^2 - 2q/(h*k) + r/k^2)/4 >= 0.

The full face star has corank1, so F+Delta=0. It follows that the isolated equality holds and Delta=0. If the spectator contact vectors of c,d are x,y and its Gram is G_Z, then at h=k=1 this equality is

    (x-y)^T G_Z^(-1) (x-y) = 0.

Positive definiteness forces x=y. Any common nonzero coordinate gives a white-interior path c-e-d of length2, immediately paying that actual corner, whether finite or ideal. Under the unpaidness assumption x=y=0. Thus every supporting spectator is orthogonal to the complete18-root equality core. This is a statement about the entire spectator block, not a bounded-rank attachment test.

## Sealed finite corners and the coupon

Number the cycle0,...,15 and its leaves16,17 attached at0,8. Supporting whites in the core are1,...,7 and9,...,15. The four adjacent marked pairs are (0,8), (8,16), (16,17), (17,0). The rational core has signature(16,1,1), with radical

    (-4,-3,-2,-1,0,1,2,3,4,3,2,1,0,-1,-2,-3,-2,2).

The selected corner is affine A15. At (8,16) and (17,0), the marked components are A15 and A1, and the A15 mark is7 edges from either endpoint, hence its unique middle. The last corner has its two marks in different A1 components, with two unmarked A7 components. All three corners are finite. Orthogonality of every extra spectator seals the displayed components as whole actual vertex components; it cannot silently turn A15 or A1 into a larger diagram.

The two half-coupons therefore contribute exactly1. They are defined using whole finite components, with the A15 marked node its unique center. No coupon is assigned merely because an A15 subdiagram happens to occur inside a larger component. The source checks every literal corner, both middle-node identifications, and the total payment.

This proves the stated local transfer. It also explains why the earlier infinite-volume local countermodel is compatible with the new rule: its two previously unpaid finite corners now pay one half each. Neither that example nor this transfer asserts the existence of a finite-volume chamber containing the equality configuration.

## Verification and review scope

Generation and fresh saved replay both terminated successfully. The certificate checks six finite single-profile diagrams,16 rational rectangle values, the formal spectator determinant correction, the literal18-root signature and radical, all four corner diagrams, two half-coupons, and the112 retained/8 removed ideal pairs. Syntax, all8 direct source hashes and all5 local report links passed. A separate internal mathematical/source review read the complete source and proof, finding no gap in the rank15 extension completeness, fixed-leaf finite-corner argument, unique rectangle equality, arbitrary-spectator equality and sealing, or whole-component half-coupon. That review did not independently rerun the rectangle or core arithmetic. The inherited complete extension classification and actual factor-fiber geometry are stated proof inputs; this package does not repeat the old full cone census or published classification tables. Source, data and report are final/frozen. No global dimension conclusion, external refereeing or proof-assistant formalization is claimed.

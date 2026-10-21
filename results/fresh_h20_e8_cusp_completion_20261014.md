# Two bounded cusp links in the 79-wall H20 extension

**Both deficient cusp links of the completed H12 face can be repaired coherently by two norm-two walls.** The resulting chamber has 79 actual facets and retains the 35 certified finite vertices of that face. Its previous strict spacelike escape survives, so the ambient H20 chamber still has infinite volume. This checkpoint stops after the requested two-cap test.

The [exact checker](../scripts/fresh_h20_e8_cusp_completion_20261014.py) and [certificate](fresh_h20_e8_cusp_completion_20261014.json) use the frozen [77-wall predecessor](fresh_h20_contact_closure_20261014_h12_replacement.json). They preserve every predecessor wall and do not alter earlier certificates. Acceptance uses exact arithmetic in `Q(sqrt(2))(sqrt(79))`; square comparisons determine signs without numerical angle fits.

## The two caps

The E8 support in the 77-wall numbering is

```text
[12,13,14,15,18,26,33,67].
```

Its Cartan inverse has exactly one diagonal entry equal to two, at wall 12. Let `omega` be the corresponding fundamental weight. It satisfies `omega^2=2`, pairs one with wall 12, and pairs zero with the other seven E8 simple roots.

Let `nu1,nu2` be the two future null rays obtained from the complete [14-wall H12 source](fresh_h12_a12_source_20261014.md). They are orthogonal to E8 and to their incident H12 cusp walls. For each ray set

```text
r_i = -omega + 8 nu_i.
```

Then `r_i^2=2`; its only nonright contact in its old cusp link is pi/3 with wall 12. All remaining nonincident predecessor walls have negative pairing with `nu_i`, so a sufficiently large positive shift makes those pairs dotted. The checker tries the common integer shifts 1, 2, 4, and 8, recording failures until the first passing value. At 8 every required strict gate holds, including the mutual dotted relation of the two new caps.

All **155 new pairs** are checked exactly: 39 right, two pi/3, and 114 strictly dotted. The common predecessor center remains strict and timelike. Projection to each of the 79 normal hyperplanes produces an explicitly verified future timelike facet point, strict on every other inequality. Thus the facet count refers to actual facets.

## Full bounded links, not rank alone

The checker recomputes each complete incident set against all 79 walls. The first cusp has 22 incident facets, in affine blocks of sizes 4, 9, and 9; the second has 21 incident facets, in affine blocks of sizes 12 and 9. In both cases the newly added nine-node block is affine E8.

For every connected block it verifies a positive Gram kernel, a positive-definite principal minor of corank one, and a weighted vector identity

```text
sum_j m_j r_j = mu nu_i,   all m_j>0, mu>0.
```

The block quotient ranks sum to **19**, the full Euclidean horosphere dimension. Consequently any link recession vector, whose pairings with the facet directions are nonpositive, has zero pairing with every direction by the positive relations. Full quotient rank forces the recession vector to vanish. The links are therefore bounded. This verifies more than the rank condition that the incomplete predecessor links already satisfied.

Each new cap projects to `8 nu_i` on the E8-perpendicular H12 space. That vector pairs strictly negatively with every future timelike point of the face. Thus the intrinsic finite-volume H12 face remains unchanged. The checker additionally verifies strict new-cap pairings at all 35 saved ambient finite vertices. It does not claim a complete new ambient corner or cusp census.

## Remaining global failure

The same saved future spacelike vector remains strict against all 79 facets, with exact norm

```text
269344645863312077 / 40000000000000000 > 0.
```

Together with the strict timelike center, this gives an open feasible patch at infinity and proves infinite volume. Completing these two links is a concrete local improvement, not a global finite-volume completion. The inherited conditional non-quasi-arithmetic statement concerns any hypothetical finite-volume completion retaining the stated old walls; no lattice is claimed for this infinite-volume chamber.

```sh
python scripts/fresh_h20_e8_cusp_completion_20261014.py --verify
```

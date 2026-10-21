# Independent check of the H20 two-angle extension

The three new caps and the resulting **85-facet H20 chamber pass** the independent exact check. All three have squared norm \(2/27\), preserve the common strict timelike center, and retain every one of the preceding 82 walls. The chamber still has infinite volume: the saved strict future spacelike point of square \(3849/32768\) satisfies all 85 inequalities strictly.

The [checker](../scripts/fresh_dimension_bridge_two_angle_20261013_independent.py) imports no construction producer or numerical package. It uses the separate Fraction implementation of \(\mathbf Q(\sqrt6)\) from the prior independent verifier and a new recursive quadratic-extension implementation. For the three cross pairs, it reverses the producer's order of adjoining radicals; sign-based equality also handles dependent radicals.

The [saved replay result](fresh_dimension_bridge_two_angle_20261013_independent.json) independently checks:

- All three explicit vectors, positive radicands and exact norm identities.
- All 246 new/old pairs: 54 right, nine \(\pi/6\), and 183 dotted; the three new/new pairs are dotted.
- Each cap's 20 prescribed pairing equations, with coefficient rank 20: 18 right angles and two prescribed \(\pi/6\) anchors. Each also has a third, unprescribed \(\pi/6\) relation.
- Strict center and escape inequalities, including the exact saved pairing values.
- Positive cutting of the three saved target rays. Each target is future spacelike, satisfies all 82 old inequalities, and has 20 independent active walls.
- Rejection of a forged normal, forged norm and forged target pairing.

The original 82-wall angle and cusp certificate is an explicitly hash-anchored dependency; this checker does not repeat its full cusp census. The new pair checks preserve global acuteness. Orthogonally projecting the common strict center onto any wall therefore proves that every retained and added inequality is a true facet. No completeness of the cap search, finite-volume completion, or new dimension is claimed.

Reproduce from the repository root:

```sh
python scripts/fresh_dimension_bridge_two_angle_20261013_independent.py --verify
```

Both the initial run and saved-result replay passed. The audited input is [the exact three-cap archive](../data/enumerations/fresh_dimension_bridge_two_angle_20261013_exact.json).

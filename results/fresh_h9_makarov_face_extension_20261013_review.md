# Review of the compact-H5 seed and changed-center cap

**PASS:** source and mathematical review, followed by a fresh saved-result
replay of the [standalone checker](../scripts/fresh_h9_makarov_face_extension_20261013.py)
and its [certificate](fresh_h9_makarov_face_extension_20261013.json).
This review uses the checker's exact arithmetic for replay; it is not a second
implementation of every root calculation.

The compactness argument for the seven-facet H5 seed is complete. The checker
enumerates every possible rank-five extreme ray, finds ten timelike rays,
and verifies that they have a common future orientation. Its timelike center
is a strictly positive combination of all seven wall normals. Thus its
normalized affine slice has no nonzero recession direction: such a direction
would pair to zero with every wall, and the first six walls span. The ray
census therefore proves compactness, rather than merely checking a bounded
list of finite vertices. Direct sum with the positive-definite H4 Gram gives
signature (9,1); all ten lifted points are actual rank-nine finite vertices.

The forty-family obstruction is also sound. For an old timelike vertex v,
H4 dual vector omega_j and transverse pairing magnitude t, the normalized
cap is

```text
r(t)/t = sqrt((B_jj - 2/t^2)/(-v^2)) v - omega_j,  t >= 1.
```

The saved escape E has (v,E)<0 for all ten vertices. The square-root
coefficient is increasing in t, so it suffices to check t=1. The forty
exact squared inequalities prove strict negativity there. Together with
the exact norm, future sign and eleven strict seed inequalities, this
excludes every finite addition from those forty families, including all
finite labels m>=3 and the larger real pairing magnitudes. It does not
exclude arbitrary cap centers or profiles.

The changed-center cap leaves that family. The inverse-Gram construction
prescribes its five old pairings, and the displayed quadratic gives exact
norm two. The final field identity is checked exactly:
`sqrt(RAD) = (16*sqrt(11)+4*sqrt(55))/11`, with positive sign. The form and
seed lie in `Q(sqrt(5))`; the new cap lies in `Q(sqrt(5),sqrt(11))`.
All eleven relations are checked: two label-three contacts,
seven right contacts and two dotted pairs. Its pairing with the former
common escape is positive. It is strict at the center and at all ten old
vertices, so the compact H5 face remains intact. Projecting the center onto
the new wall supplies an actual facet point with every other inequality
strict. The new saved future spacelike point is strict for all twelve
walls; consequently this chamber still has infinite volume.

The companion [H4 transverse gate](fresh_h9_h4_transverse_gate_20261013.md)
also passes review and saved replay. Positivity of the H4 inverse gives
`b^T B b >= 3+sqrt(5)>2` whenever a norm-two acute Coxeter normal has a
nonzero H4 pairing. Its orthogonal component is timelike, so that wall
cannot meet the H4 face, even at infinity. The conclusion is a restriction
on walls meeting this face, not a completion obstruction for the seed.

Replayed commands:

```sh
python3 scripts/fresh_h9_makarov_face_extension_20261013.py --verify
python3 scripts/fresh_h9_h4_transverse_gate_20261013.py --verify
```

The first replay reports 152 accepted and 48 rejected proposals in the
separate bounded 200-cap screen. Those counts do not assert compatibility
among the accepted proposals. The all-label family obstruction does not
depend on that screen.

Reviewed SHA-256 bindings:

| File | SHA-256 |
| --- | --- |
| Seed/cap checker | `2135abd59a7acea1cf4901fe0e91a4bf878a9413d8d7eff4b4022f6806186d3a` |
| Seed/cap certificate | `d7f43b1e5402dcbbf85094c63f56778b22908e4f08c1169c53ab72c04d50e0d9` |
| H4 gate checker | `276efa894ff60e94cb30a80a4c79f4d51475eccb58b4981f8ed0a9cb845c8f43` |
| H4 gate certificate | `c696daadde753726a4d6599c6255ff27ae599cd158fd61bedbd4088e82dbb9ee` |

No mathematical or scope issue was found. This is a new geometric search
seed and a certified cap within this experiment, not a newly established
finite-volume H9 example or a novelty claim for the known H5 prism.

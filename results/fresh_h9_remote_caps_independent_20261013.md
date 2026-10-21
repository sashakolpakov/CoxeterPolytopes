# Review and higher-precision replay of the remote H9 insertions

**PASS.** The four-cap extension to 83 walls and the subsequent two-cap extension to 85 walls have sound exact norm/angle identities, genuine finite vertices, and certified strict spacelike escapes. Both stages still have infinite volume. The second insertion gives an actual two-edge path from an existing finite vertex through a new H3+6A1 vertex to a new 9A1 vertex.

The separate [replay script](../scripts/fresh_h9_remote_caps_independent_20261013.py) reruns the constructors' outward-rounded integer interval checker at **100 decimal digits**, then checks that all 949 resulting enclosures lie inside the corresponding archived 80-digit enclosures. All non-interval metadata is also compared exactly. Its [saved result](fresh_h9_remote_caps_independent_20261013.json) records source hashes, counts, and exact escape norms. This is a source-level proof review and a different-precision replay of the same arithmetic implementation, **not a second independent field implementation**. The separately verified, hash-anchored 79-wall input remains a dependency.

| Stage | New walls | New/old right pairs | New/old dotted pairs | Additional finite relation | Enclosures replayed |
|---|---:|---:|---:|---|---:|
| 83 walls | 4 | 20 | 296 | none | 457 |
| 85 walls | 2 | 14 | 151 | one π/3 | 492 |

For the first insertion, the exact Q(φ) frame has one negative and four positive orthogonal directions. Rescaling gives squares (−1/2,1/2,1/2,1/2,1/2). A rational boost satisfies `ct²−st²=1`; four Hadamard sign patterns then give four norm-two, mutually orthogonal normals. They are exactly orthogonal to five retained orthogonal walls. The displayed point of square −1/2 is orthogonal to those nine independent normals and strict on all other walls, proving an actual 9A1 finite vertex. Every new/old relation outside those five walls is certified strictly dotted, and all four new normals cut the targeted old escape. The saved replacement escape is strict on all 83 walls.

For the second insertion, the rescaled frame has squares (−1,1,1). With `k>1`, the two normals are

```
r± = ((k−1/k)/2)e0 + ((k+1/k)/2)e1 ± e2.
```

They have norm two and are mutually orthogonal. If α,β,γ are the anchor pairings of e0,e1,e2, the condition `(r+,anchor)=−1` is exactly

```
((α+β)/2) k² + (γ+1) k + ((β−α)/2) = 0.
```

The checker defines k by the quadratic radical formula and proves that the selected branch is the unique root greater than one; the other branch is less than one. Thus the π/3 identity is symbolic, not an approximate interval-equality test. All remaining nonright relations are strictly dotted.

The joint 9A1 vertex is explicit. The H3+6A1 vertex and the retained H2+7A1 vertex are obtained by exact orthogonal projection of a strict timelike center to positive-definite nine-wall spans. Their physical Grams are anchored to the original exact pair certificate; inverse-Gram and final strict-inequality checks verify the actual vertices. Each claimed edge has eight common spherical walls and two distinct feasible finite endpoints. Convexity and strictness of the other inequalities give a true finite edge, rather than just a formal diagram adjacency. The final strict future spacelike escape proves infinite volume of the 85-wall stage.

The sign orientation was also checked: the saved positive covector h corresponds to coefficient center `−G⁻¹h`, so positive root/h products mean strict negative halfspace evaluations at that center. The same convention is handled consistently for the saved target and escape covectors.

The preliminary rank-seven/rank-six/rank-five projection census is only a census of projected-target signs for the specified seed subsets. It is not a nonexistence theorem for other cap constructions. This audit does not claim a complete 85-wall vertex census, an exhaustive cap search, or finite-volume impossibility for further extensions.

Primary artifacts are the [83-wall input](../data/enumerations/fresh_h9_remote_caps_20261013_caps.json), [83-wall certificate](../data/enumerations/fresh_h9_remote_caps_20261013_verify.json), [85-wall input](../data/enumerations/fresh_h9_remote_caps_20261013_coupled_caps.json), and [85-wall certificate](../data/enumerations/fresh_h9_remote_caps_20261013_coupled_verify.json). The reviewed checker sources are [the four-cap checker](../scripts/fresh_h9_remote_caps_20261013_verify.py) and [the two-cap checker](../scripts/fresh_h9_remote_caps_20261013_coupled_verify.py).

Reproduce from the repository root:

```sh
python scripts/fresh_h9_remote_caps_independent_20261013.py --verify
```

The replay passed after the final source/hash, strict-center, and preliminary-census checks were added.

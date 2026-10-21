# Independent review of the H22 face-orbit family

**PASS.** I reviewed the [exact checker](../scripts/fresh_h22_gap_face_family_20261013.py) and [family report](fresh_h22_gap_face_family_20261013.md), then ran

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h22_gap_face_family_20261013.py --verify
```

The replay passes for 241,920 new orbit directions and 2,736,428 actual facets. This is an independent source/proof review plus a fresh replay of the existing arithmetic, not a second arithmetic implementation. The previous chamber and orbit certificates are explicit hash-anchored dependencies; their nested source manifests are also checked.

The 22 independent linear equations determine `p` exactly, with `p²=−31/3` and `(h,p)=−36`. Integral generator checks preserve the form, center and old walls. Array closure and an earlier-predecessor argument prove both connectivity from the representative and completeness under the saved generators. The already certified invariance of each prior cap family justifies transfer of every representative cross-pair check to the whole new orbit.

For `q=1/(2cos²(pi/k))`, `k>=6`, one has `1/2<q<=2/3`. The nearest distinct old projections pair by `−35/3`, so the corresponding same-sign new walls pair by `q−4/3<=−q`; all other distinct pairs are more negative. Opposite-height pairs are dotted. The strict dotted bounds against old timelike caps, old null caps and center caps correctly include their different norms and height factors. Thus the proof covers every integer `k>=6` and every retained center label `m>=3`, without a bounded label enumeration.

The complete incidence at the displayed `k=6` null ray is justified: the 18 common original normals and one affine `A1` supply quotient rank 19, while every other cap is strict there. This leaves two unbounded directions in the 21-dimensional horosphere. The stronger volume conclusion is the saved strict spacelike point, whose norm is `647099999951/35000000000000000>0`. I checked the factors in the squared comparisons: `704` for the center family, `440` for the earlier timelike cap orbit, and `4v²>11t²d²` for null caps stored with four-times height. The new-family bound `R·cov+132d<0` is exactly the worst-height inequality. These bounds prove the same escape works uniformly in both labels.

The old escape is cut even at the limiting smallest new height. The test `390*z²>(6p·x)²`, with negative old projection pairing, uses precisely `36*(31/3+1/2)=390`; choosing the appropriate height sign gives a strict cut for every finite `k`.

## Small labels and the actual-facet premise

The report's small-label algebra is correct. The nearest pair at `k=3` has positive pairing `2/3`; at `k=4` its normalized negative pairing is `1/3`; at `k=5` it is `sqrt(5)/3`, strictly between `cos(pi/4)` and `cos(pi/5)`. The latter two are forbidden Coxeter angles.

The following direct projection supplies the actual-facet premise even at `k=3`, where new/new acuteness fails. For every `k>=3`, `q<=2` and `H=sqrt(31/3+q)<4`. All new/old pairings remain nonpositive: the earlier bounds give `−8+4sqrt(5/2)<0`, `−5+4<0`, and `−36+16<0`, while original-wall pairings were already nonpositive. Every new wall has norm `q` and center pairing `−36`. For a selected new normal `r`, project the common center onto its wall:

`x = h + (36/q) r`.

This point is timelike. Every old inequality is strict at `x`. For any other new normal `r'`, the orbit calculation gives `(r,r')<q`, including opposite-height pairs, and therefore

`(x,r') = −36 + (36/q)(r,r') < 0`.

Thus every new inequality is a genuine facet for the small labels as well. The positive pair at `k=3` violates the necessary acute Gram condition for a Coxeter chamber; the forbidden pairs at `k=4,5` likewise cannot be dismissed as redundant inequalities. This closes those small labels only for the stated full signed orbit ansatz.

No issue was found in the main conclusions. They do not exclude proper orbit subsets, different projections, further caps or other H22 constructions.

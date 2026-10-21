# Narrow review of the 94-wall compact H3 face

The four-cap formula and the compact tetrahedral face argument pass source review. The frozen producer-free verifier also passes, and a second run at 160 decimal places, twice its archived precision, passes all reconstructed inequalities and the strict escape check.

For the projected anchor \(a^2=q=(4-2\phi)/3>0\), with \(P^2=-1\) and \((a,P)=0\), the cap
\[
r=-a/q+\sqrt{1/q-2}\,P
\]
has norm two and pairs −1 with the anchor. Its other three companions lie in the Euclidean space perpendicular to both \(a\) and \(P\); their explicitly given coordinates have Gram matrix \(2I_3\). Consequently all four new caps have norm two and are mutually orthogonal. The remaining claimed strict inequalities are certified by outward rational intervals.

The exact four-wall 3–5–3 Gram matrix is Lannér, perpendicular to the rank-six common face diagram. Its four future timelike dual rays satisfy every other wall strictly. Their projectivized convex hull is therefore an actual compact H3 tetrahedron in the 94-wall chamber. The saved strict spacelike witness still proves that the entire H9 chamber has infinite volume.

This review is independent algebraic/source reasoning plus a higher-precision replay of the existing interval implementation. It is **not** a second implementation of all 94 roots or their angle table.

- [Review certificate](fresh_h9_compact_next_20261013_source_review.json)
- [Review/replay script](../scripts/fresh_h9_compact_next_20261013_source_review.py)
- [Frozen source archive](../data/enumerations/fresh_h9_compact_next_20261013_caps.json)

```sh
PYTHONPATH=src:. python scripts/fresh_h9_compact_next_20261013_source_review.py --verify
```

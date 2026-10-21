# Polar truncations and two changed-core probes

The [numerical scout](../scripts/fresh_h20_polar_frontier_20261014.py) samples projective extreme rays of the 98-wall chamber and three variants deleting wall 8, wall 15, or both. Two simple rational polar normals pass its angle screen on the full chamber. Their subsequent [exact construction and independent review](fresh_h20_polar_pair_20261014.md) certify a 100-wall chamber with sixteen additional finite vertices. The strict escape survives; no finite-volume example is established.

The deletion choices come from the [exact escape dual](fresh_h20_s75_completion_20261014.md): any escape cutter retaining the old walls must have a nonright contact with wall 8 or wall 15. Removing those walls tests a different core rather than repeating the same face-contact family.

| Deleted walls | LP objectives | Distinct sampled rays | Spacelike rays | Polar angle passes |
| --- | ---: | ---: | ---: | ---: |
| None | 5,322 | 4,242 | 4,240 | 2 |
| 8 | 3,322 | 897 | 897 | 1 |
| 15 | 3,322 | 1,329 | 1,329 | 2 |
| 8 and 15 | 3,322 | 520 | 520 | 2 |

The objectives include fixed random directions and directions concentrated near the current escape. The changed-core passes numerically recover removed walls or one of the same full-core polars. None of the passing proposals cuts the saved escape. Many sampled rays do cut it, but their polar hyperplanes fail finite-angle conditions. These are bounded numerical observations, not exclusions of polar completion or changed cores.

The original physical coordinates are badly conditioned: the newest walls have very large boosts, and the rational escape has coordinates of order `10^10` although its square is about 20. The script decodes the exact radical towers and constructs an orthonormal Lorentz chart at 90 decimal digits before running floating-point LP and angle tests. It verifies the escape's exact saved square numerically at that precision. Ordinary direct double-precision norm subtraction was rejected as unreliable. Only the separate exact construction accepts walls.

The four saved numerical transcripts are [full core](fresh_h20_polar_frontier_20261014.json), [delete 8](fresh_h20_polar_frontier_20261014_remove_8.json), [delete 15](fresh_h20_polar_frontier_20261014_remove_15.json), and [delete both](fresh_h20_polar_frontier_20261014_remove_8_15.json). All four generation and saved-result replays pass with the recorded seed and sample counts.

```
python3 scripts/fresh_h20_polar_frontier_20261014.py --samples 5000 --verify
python3 scripts/fresh_h20_polar_frontier_20261014.py --samples 3000 --remove 8 --verify
python3 scripts/fresh_h20_polar_frontier_20261014.py --samples 3000 --remove 15 --verify
python3 scripts/fresh_h20_polar_frontier_20261014.py --samples 3000 --remove 8 15 --verify
```

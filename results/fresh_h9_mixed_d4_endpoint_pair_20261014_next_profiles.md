# Finite-endpoint profiles for the fourteen-wall mixed-D4 chamber

The [exact profile table](fresh_h9_mixed_d4_endpoint_pair_20261014_next_profiles.json) covers all **71 open directions** in the frozen fourteen-wall chamber's135-direction finite-corner queue. It contains **2,221 finite spherical profiles and92 isolated I2(m) tails**, from13,239 finite trials. These are complete necessary local endpoint possibilities, not globally compatible caps or a finite-volume completion.

The JSON preserves the original queue indices: `profiles[str(direction_index)]` contains `support`, `finite_profiles` and `dihedral_tail_profiles`. The21 directions starting at the three new finite vertices contribute459 finite profiles and21 tails. The other50 open supports reproduce their frozen twelve-wall profile lists exactly.

| New finite ray | Incident type | Open direction indices | Finite profiles | Tails |
|---|---|---|---:|---:|
| 18 | A1+A8 | 9,11,12,13,14,15,16 | 177 | 8 |
| 67 | I2(5)+D7 | 90,92,93,94,95,96,97 | 106 | 4 |
| 74 | A1+E8 | 108,109,110,111,112,113,114 | 176 | 9 |

Direction108 retains E8 and permits only an all-right finite extension. Direction112 retains A1+A4+A2+A1 and has72 finite profiles plus two unbounded dihedral tails, giving a broader constructive target. Every tail allows a single contact of every integer label m≥6 to an isolated old A1, with all other support contacts right.

The [checker](../scripts/fresh_h9_mixed_d4_endpoint_pair_20261014_next_profiles.py) rechecks all91 physical old-wall pairs with the separate quartic Fraction arithmetic. It then normalizes the eight support roots to square2 and performs exact inverse/Schur tests in Q(√5). The proposed new node is14, avoiding collision with the actual new walls12 and13. Every retained edge support is checked against the bound finite start and spacelike cone endpoint. The complete cone census and the physical Lorentz frames are inherited rather than recomputed.

The all-label completeness argument is the previously [reviewed cycle/star/path reduction](fresh_h9_mixed_d4_endpoint_profiles_20261014_scope_review.md): at most one attachment per old component, at most three attachments overall, three attachments necessarily all label3, and every label at least6 isolated in a rank-two component. All remaining3/4/5 profiles are accepted exactly when their positive Schur complement is nonzero. This excludes no affine endpoint or vertex-free construction; those lie outside the stated finite-endpoint task.

Generation and saved replay pass:

```sh
python scripts/fresh_h9_mixed_d4_endpoint_pair_20261014_next_profiles.py --verify
```

# New frontiers of the 29-wall chamber

The [complete-line census](fresh_h9_octic_base_lines_20261015.json), produced by its [exact checker](../scripts/fresh_h9_octic_base_lines_20261015.py), finds 72 complete geodesic one-faces without finite endpoints. Of 3,477 relevant exterior-ray pairs, 1,810 fail the clique test and 1,595 fail positive definiteness of the support. Two supports are new relative to the preceding chamber:

* `[0,1,2,5,7,10,27,28]`;
* `[0,1,2,6,7,10,27,28]`.

The full extreme-ray census is inherited from the independently reviewed 29-wall chamber. This new derived frontier has exact generation and saved replay; no separate complete-line implementation has been run on this checkpoint.

The [new finite-corner frontier](fresh_h9_octic_base_germs_20261015.json), with its [checker](../scripts/fresh_h9_octic_base_germs_20261015.py), contains all nine open germs at the sole new vertex introduced by wall 29. It records 122 finite spherical contact profiles and five symbolic dihedral tails. It does not repeat the older 285 open germs. The limiting null points have only their eight support walls active and unbounded links.

Two bounded numerical proposal searches target the exact escape surviving both octic alternatives:

| Frontier | Distinct numerical proposals | Escape cutters |
|---|---:|---:|
| 72 complete lines | 68,526 | 0 |
| Nine new finite-corner germs | 2,173 | 60 |

The [line scout](fresh_h9_octic_base_scout_20261015.json) and [germ scout](fresh_h9_octic_germ_scout_20261015.json) bind their full compressed banks. They use bounded anchor labels, sampled symbolic tails and floating-point pair gates. Their zero counts are not all-parameter exclusions. The germ scout's inherited `scope` string incorrectly says “24” germs from “3” vertices; its actual inputs, nine edge records and enumeration cover the nine germs described here. This note corrects that text without changing the frozen certificate.

All 60 germ cutters come from support `[0,1,2,5,6,8,10,28]`. One nine-contact recipe has now been [exactified in the original quartic field](fresh_h9_nine_contact_cap_20261015.md); it does not require an octic extension. Generation and saved replay pass for both frontiers and both numerical scouts.

```sh
python scripts/fresh_h9_octic_base_lines_20261015.py --verify
python scripts/fresh_h9_octic_base_germs_20261015.py --verify
python scripts/fresh_h9_octic_base_scout_20261015.py --verify
python scripts/fresh_h9_octic_germ_scout_20261015.py --verify
```

# Independent review of the π/7 incidence branch

**PASS.** Source review found no issue in the primary checker's norm formula, root isolation, finite-angle identities, outward-rounded interval gates, strict centers, or target cuts. A separate implementation independently verified two representative caps against all 812 parent normals.

Files: [separate Fraction checker](../scripts/fresh_h22_incidence_prism_20261013_independent.py), [audit output](fresh_h22_incidence_prism_20261013_independent.json), and [primary checker](../scripts/fresh_h22_incidence_prism_20261013_verify.py). The independent checker imports no producer or project arithmetic module. It reconstructs the old form, incidence bipartition, and all 812 parent normals from raw input and binary profiles. It isolates T and the positive square roots by rational bisection, independently of the primary checker's fixed-point intervals.

| Recipe | Candidate | New–parent checks |
|---|---:|---|
| 0 | 3 | 21 right, two π/7, one π/4, 788 dotted |
| 22 | 52 | 21 right, one π/7, 790 dotted |

Both norms are exactly two. Both centers are strictly negative, both prescribed positive-square parent rays are feasible and strictly cut, and the two new caps have a strictly dotted mutual relation. All **1,624** new–parent checks passed.

For the primary twice-squaring test, write the actual pairing as `p=LT+B√d`, with T>0 and d>0, and let E be the required squared pairing. Set `Y=E−L²T²−B²d`. The checked polynomial identity `Y²=4L²B²T²d`, together with `sign(Y)=sign(LB)`, proves `Y=2LBT√d`, and hence `p²=E`. If LB=0, checking Y=0 proves the same identity directly. When E=0, this proves p=0; no separate orientation choice is needed. When E>0, the strictly negative pairing interval selects the Coxeter angle. Thus the squaring introduces no extraneous accepted branch. In these two representative recipes all finite contacts have LB=0, so the nonzero-cross-term case was reviewed algebraically rather than exercised by those caps.

The primary fixed-point addition, multiplication and square-root functions round the lower endpoint down and the upper endpoint up. The cubic is strictly increasing on [1.8,1.9], so sign bisection selects the intended `T=2cos(π/7)`. The norm formula follows from u⊥v and `d=(q−q²u²T²/4)/v²`. The separate checker verifies this identity in the quotient polynomial ring.

Scope: this is a complete source/formula review plus a separate arithmetic replay of **two** caps. It does not independently replay the other 21 caps, the full 835-wall volume calculation, or any later extensions. The primary all-cap certificate remains the source for those pair counts. No finite-volume conclusion follows from this audit.

```sh
python scripts/fresh_h22_incidence_prism_20261013_independent.py --verify
```

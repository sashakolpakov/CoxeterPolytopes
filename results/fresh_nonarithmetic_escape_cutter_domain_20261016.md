# Every escape cutter misses all three selected caps

Let E be the renewed exact strict spacelike escape of the [filtered52-wall chamber](fresh_nonarithmetic_affine_filtered_20261016.md). **Every spacelike normal r acute-compatible with that chamber and satisfying `(r,E)>=0` is strictly dotted to each of its new caps49,50,51.** This is an all-real geometric statement: it imposes no Coxeter-angle alphabet, field, norm normalization, or finite-vertex premise.

The [exact certificate](fresh_nonarithmetic_escape_cutter_domain_20261016.json) and [checker](../scripts/fresh_nonarithmetic_escape_cutter_domain_20261016.py) consider normals as vectors in the polyhedral cone

```
C_E = {r : (r,n_i)<=0 for all52 old walls, (r,E)>=0}.
```

The additional inequality is a diagnostic on possible normals. The checker does not claim that adding the hyperplane with normal `-E` produces a Coxeter chamber.

It reproduces the complete frozen52 cone and clips it by this additional halfspace. Of the2,457 old extreme rays,146 have `(r,E)>0` and2,311 have `(r,E)<0`; none has zero pairing. The196 crossing edges give a complete resulting cone of342 extreme rays:341 spacelike and one timelike. All directions have a strict future chart relative to the inherited H. Exact positive crossing formulas and rational outward intervals certify the signs and complete common-zero-face adjacency.

For each resulting ray v and each cap n in49,50,51, the checker verifies

```
p = v - (v,n)*n/2,
p² = v² - (v,n)²/2 < 0,
(p, H - (H,n)*n/2) < 0.
```

Thus all342 projections lie strictly inside the same future timelike cone in `n-perp`. Every nonzero vector in `C_E` is a nonnegative combination of the complete extreme-ray list. Its projection is therefore future timelike as well: positive sums of future timelike vectors remain future timelike, and cannot vanish. If r itself is spacelike, this gives `(r,n)²>2*r²`. Since r is acute-compatible with n, its pairing is negative, proving strict dottedness. The closed boundary `(r,E)=0` is included in the certificate, so there is no omitted equality case or zero projection.

This proves that **one added wall cannot both cut E and meet any selected cap**. It does not exclude a simultaneous collection with separate roles: noncutting walls can supply the finite vertices required on the new caps' pi/10 ridges, while other walls cut E. Any finite-volume completion retaining the52 chamber must contain both kinds, hence at least two additional walls. A single noncutting wall might meet several caps; no larger lower bound is asserted.

The certificate is a complete producer cone replay with exact interval signs, using frozen algebraic cap recipes. It is distinct from any independent second implementation of the52 cone. Generation and a fresh saved `--verify` replay both passed, as did source hashes, syntax and links. The four [bounded bridge pilots](fresh_nonarithmetic_ridge_bridge_20261016.md) are numerical diagnostics and are not used in this proof.

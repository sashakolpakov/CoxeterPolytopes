# Independent review of the perpendicular S61 escape obstruction

**PASS.** The [gate](fresh_h20_s61_multicap_20261014_perpendicular_gate.md) correctly excludes finite-volume completion by any finite coherent family of additional acute walls perpendicular to S61, retaining the frozen 133-wall chamber. There is no restriction on real coefficients, norm, field, radius or Coxeter label. Nonright contacts to S61 and replacement of retained walls are outside the conclusion.

The [independent checker](../scripts/fresh_h20_s61_multicap_20261014_perpendicular_gate_review.py) and its [certificate](fresh_h20_s61_multicap_20261014_perpendicular_gate_review.json) bind the primary source/data and the preceding [independent 133-wall review](fresh_h20_s61_multicap_20261014_exact_review.md). They import no producer arithmetic. Using the prior independent sparse radical field and 512-bit outward dyadic intervals, this small replay reconstructs the physical Lorentz frame, the rational escape and its orthogonal projection, all 29 polygon-side signs, and all 29 rational active-frame inverses. It inherits the complete compact polygon and ambient geometry from that preceding review; it does not repeat the full 3,107-pair or polygon-intersection census.

The escape has square `38017800175279/500000000000`. Its projected frame coordinates have positive time coordinate, projected square less than `−76`, and pairings less than `−15` with every polygon side. The frame is checked directly against the ambient rational form and is perpendicular to all 18 S61 roots. The sign change from the time pairing to the time coordinate is necessary and is implemented correctly.

All 29 active 20-root Grams are positive definite with entrywise nonnegative inverses: two old corners involving wall 95, two chain endpoints, and 25 internal chain corners. The physical normalization and actual incidence of these frames were independently certified in the preceding review.

The proof uses the **proposed extension's own common future timelike strict interior point**, not an assumption that the old center survives. For any actual finite vertex `v` and its active normals `s_i`, decompose that point as `h = αv + Σλ_i s_i`. Future orientation gives `α > 0`, while `λ = G⁻¹((h,s_i)) ≤ 0`. Acuteness gives `(n,s_i) ≤ 0`, hence `Σλ_i(n,s_i) ≥ 0`. Since `(n,h) < 0`, necessarily `(n,v) < 0`. Thus every proposed wall strictly retains every polygon vertex.

The future projection `E_F` lies strictly inside the complete compact Klein polygon. After positive-time normalization it is a convex combination of its 29 vertices, so `(n,E_F) < 0`. Perpendicularity to S61 then gives `(n,E) = (n,E_F) < 0`. This is the only step requiring perpendicularity of the added wall.

For a finite family, the original spacelike `E` remains strict for every wall. The segment from a common future strict interior point toward `E` meets the future light cone while all inequalities are strict. A neighborhood of that null direction remains feasible, giving an open ideal-boundary region and infinite volume. This does not require the old center to survive or assume the pairing of `E` with the new center has a prescribed sign. Infinite limits of wall families are not covered.

Generation and saved replay pass:

```sh
python scripts/fresh_h20_s61_multicap_20261014_perpendicular_gate_review.py --verify
```

No mathematical or scope issue was found. Only source/data files are hash-bound, so the primary report may add this review link without invalidating the certificate.

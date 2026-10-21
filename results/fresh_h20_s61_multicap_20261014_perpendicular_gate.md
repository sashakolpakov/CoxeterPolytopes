# Further perpendicular S61 walls cannot remove the global escape

For the [133-wall chamber](fresh_h20_s61_multicap_20261014.md), **every coherent acute addition perpendicular to S61 preserves the saved strict spacelike escape**. This holds for arbitrary real spacelike normals, without a norm, field, angle-label, radius, or parameter bound. Consequently no finite collection of such additions can give finite volume while retaining this chamber. Nonright S61 contacts and replacement of retained walls remain outside this theorem.

The [small exact certificate](fresh_h20_s61_multicap_20261014_perpendicular_gate.json) is produced by [this checker](../scripts/fresh_h20_s61_multicap_20261014_perpendicular_gate.py). It binds the frozen 133-wall construction and reuses its field and interval helpers; it is not a separate arithmetic implementation.

Let `E` be the saved rational escape, whose square is `38017800175279/500000000000`. Its orthogonal projection `E_F` to the S61 Lorentz three-plane is future timelike. The checker reconstructs this projection from the physical old-field vectors and proves strict negative pairings with all 29 polygon sides. Thus the positive-time normalization of `E_F` lies in the compact polygon, which is the convex hull of its 29 finite vertices. This step uses the complete bounded polygon census, rather than a finite sample of interior points.

The finite-vertex preservation argument applies to any new acute wall. Let `v` be one of these actual finite vertices, let its 20 active normals be `s_i`, and let their positive-definite Gram matrix be `G`. The checker verifies entrywise nonnegativity of `G⁻¹` for all 29 active frames. If `h` is a common future timelike strict interior point of the retained chamber and the proposed new wall, write

`h = αv + Σ λ_i s_i`, where `α=(h,v)/(v,v)>0` and `λ=G⁻¹((h,s_i))≤0`.

For an acute new normal `n`, every `(n,s_i)≤0`. Hence `Σλ_i(n,s_i)≥0`. Since `(n,h)<0`, it follows that `(n,v)<0`. This proves strict preservation of every polygon vertex, even when the common interior point depends on the proposed added walls.

Linearity then gives `(n,E_F)<0`: normalize all future points in the same Klein chart and use the convex-hull description of the compact polygon. When `n` is perpendicular to S61, orthogonal projection gives `(n,E)=(n,E_F)<0`. The old escape therefore remains strict at every proposed wall. Together with a common future interior point, it still yields an open ideal-boundary region and infinite volume.

The hypothesis “coherent acute addition” includes a common future timelike strict interior point. It is not enough merely to list pairwise nonpositive normals without checking that they define the intended chamber. The conclusion concerns finite additions, as appropriate for a finite-sided Coxeter polytope; it does not assert anything about limits of infinitely many walls.

Replay:

```sh
python scripts/fresh_h20_s61_multicap_20261014_perpendicular_gate.py --verify
```

Generation and saved replay both pass, including all 29 projection inequalities and all 29 active-frame inverse checks. Source manifests, syntax, and report links also pass.

# Independent review of the two-distance polar cap

**PASS: the new cap is valid and cuts the targeted escape, but the resulting 45-wall H9 chamber has infinite volume.** The [independent checker](../scripts/fresh_h9_quartic_completion_20261013_polar_independent.py) reconstructs the cap directly from its displayed equations, then checks the [saved construction](../data/enumerations/fresh_h9_quartic_completion_20261013_polar.json). It imports no search or producer module. See the [exact PASS record](fresh_h9_quartic_completion_20261013_polar_independent.json).

With `phi^2=phi+1`, set

```
A=-20+25phi/2, B=-25+15phi, C=-42/5+24phi/5,
D=B^2-4AC=-62+39phi,
s=(-B+sqrt(D))/(2A), delta=sqrt(2+phi).
```

Exact signs give `A>0`, `C<0`, `D>0`, and `3<s<4`. Therefore the displayed branch is the unique positive root of the norm equation. The only nonzero seed pairings are `-delta*s` and `-delta*(s+1)` at zero-based indices 4 and 8. Multiplication by the authoritative seed Gram inverse yields norm exactly 2 and the correct strict center orientation.

All 44 pairings are checked independently. Eight seed relations are right angles, two are dotted; old cap 348 has angle pi/3, the other 21 retained old caps are dotted, and all 12 preceding new caps are dotted. The checker embeds the earlier quadratic towers in explicit composita over the new field, so this comparison does not assume that the new radicand equals the earlier anchor radicand. It uses exact zero, pi/3, and dotted inequalities, with no numerical angle tolerance or finite-label truncation.

The cap cuts exactly saved target 3. Its seed zero pattern also differs from every one of the 3,030 literal preceding entries in the cited dictionaries, proving distinctness of its mirror regardless of normalization. This is a new mirror within this search, with no claim of literature priority.

The saved volume witness has all ten seed inequalities strict, all 35 nonseed inequalities strict, positive Lorentz square, and negative pairing with the common future timelike center. These exact facts prove that an open set of ideal directions remains in the 45-wall cone, hence its volume is infinite. Numerical halfspace vertices only proposed the witness; their enumeration is not used to certify the conclusion. Both a forged cap coefficient and a forged escape are rejected.

```sh
python -m scripts.fresh_h9_quartic_completion_20261013_polar_independent
```

No issue was found in the construction. This escape applies to the literal 45-wall chamber and its wall subsets. A subsequent [four-wall cusp obstruction](fresh_h9_fixed_core_cusp_obstruction_20261013.md) additionally excludes finite-volume completion by adding acute-compatible walls while retaining the inherited core. Replacing inherited walls remains outside that theorem.

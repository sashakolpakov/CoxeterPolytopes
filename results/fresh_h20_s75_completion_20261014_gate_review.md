# Independent review of the S75 escape dual

**PASS.** The [independent checker](../scripts/fresh_h20_s75_completion_20261014_gate_review.py) and its [saved certificate](fresh_h20_s75_completion_20261014_gate_review.json) verify the exact identity behind the [S75 completion gate](fresh_h20_s75_completion_20261014.md). Generation and saved replay pass. No mathematical issue was found.

The checker imports no producer or project arithmetic. It multiplies squarefree radicals using their integer gcd, encloses coefficient signs with outward rational square-root bounds, and reconstructs all 21 physical coordinates of the saved rational escape. An exact finite-field specialization proves the 21 physical root columns independent. This is a separate rank check from the producer's 336-dimensional rational regular representation. All denominators are units in the chosen finite field, and the specialization is checked against the entire radical multiplication table.

The support inverse is checked by rational matrix multiplication. Its nonnegative entries and the diagonal threshold `J_ii >= 2` force the ten stated right contacts for any Coxeter normal whose projection gives a geodesic boundary on S75. This applies to every finite label and also rules out parallel or dotted contacts at those nodes. It does not depend on the bounded numerical profile scout.

The physical identity has exactly two negative coefficients, at walls **8 and 15**, and 19 strictly positive coefficients. Thus a nonzero normal right to those two walls and acute to the other 19 has strictly negative pairing with the escape. Equality would force it perpendicular to all 21 independent roots. The original ten-right-contact formulation follows immediately.

Consequently every acute escape cutter must be nonright to wall 8 or 15. Their S75 inverse diagonals are `24/11` and `4`. For a Coxeter normal of square `q>0`, the first such contact already gives projected square at most `-q/11`. In a coherent extension the projection is future timelike: projecting the extension's own strict center off the acute positive-definite support preserves strictness. Such a wall is strictly negative everywhere on S75 and adds no boundary there. Cutting this escape and bounding S75 therefore require different walls while the premises are retained.

The independently checked escape square is

`20181174579467302741 / 10^18 > 0`.

A finite family consisting only of S75 boundary-producing caps keeps this same strict escape and has infinite volume. This is not a proof against arbitrary mixed additions or changed retained walls. It also does not prohibit a compact polygon on S75 whose walls fail to cut this escape.

The review binds the primary source/data, the inherited chamber, and the complete support-contact gate. It does not duplicate the old 98-wall pair or vertex census, and it does not certify the numerical scout counts.

```sh
python3 scripts/fresh_h20_s75_completion_20261014_gate_review.py --verify
```

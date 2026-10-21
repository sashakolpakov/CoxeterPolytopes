# Independent review of the E8+A1 chamber obstruction

Source and mathematical review **PASS**. I also reran both frozen global checks successfully:

```sh
PYTHONPATH=src:. python scripts/fresh_h9_e8_next_20261013_global.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_e8_next_20261013_global_independent.py --verify
```

The [global certificate](fresh_h9_e8_next_20261013_global.json) and [separate Fraction reconstruction](fresh_h9_e8_next_20261013_global_independent.json) support the full-host conclusion in the [main report](fresh_h9_e8_next_20261013.md). This review did not repeat the complete norm-orbit enumeration or independently implement PARI.

Both matrices are reconstructed from oriented simple frames and timelike corner points. The actual E8 lattice basis is used to check each matrix and its inverse, so rational changes of coordinates are not confused with lattice automorphisms. Positive scales, form preservation and future orientation are explicit. The complete first-edge certificates and their independent replays are hash-bound; the three oriented paths connect the initial E8+A1 sector to the D6+3A1 and second E8+A1 sectors of the same full reflection chamber. The finite reflection closures have 121, 33 and 121 mirrors. Closure counts alone would not establish completeness; here that premise comes from the complete starting stabilizer and first-edge/endpoint certificates, including the simultaneous-root cone test on the D7 edge.

A lattice automorphism permutes the full reflective arrangement. Its positive map of a complete corner sector therefore sends the chamber to itself. This applies to both involutions H and J. Their product has exact trace greater than ten and hence infinite order. Discreteness of the admissible integral orthogonal group makes an infinite chamber stabilizer incompatible with finite chamber volume: a small generic interior ball has infinitely many disjoint stabilizer translates. Thus the specified host `&lt;-3(T-1)&gt;+E8+&lt;1&gt;` is nonreflective.

The finite fourteen-wall escape is a separate, weaker observation. The global proof does not claim a lattice identification with D9, nor nonexistence of arbitrary H9 Coxeter polytopes. Completeness of the relative norm/unit interface remains the explicitly stated shared arithmetic dependency. The report's five-edge, 71-profile and 46-interval-point totals agree with the saved global output; only the three displayed path segments are needed for the chamber-symmetry argument.

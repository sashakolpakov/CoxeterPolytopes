# Bounded numerical frontier after the fourth H9 cap

The [replay script](../scripts/fresh_h9_after_fourth_rescore_20261015.py) and [saved certificate](fresh_h9_after_fourth_rescore_20261015.json) rescore the existing bounded bank of 68,111 proposals against the frozen [29-wall chamber](fresh_h9_fourth_escape_extension_20261015.json). They reproduce **67,056 numerically compatible proposals and 152 numerical cutters of its saved escape**. No new cap is certified against the full chamber, and no extension or finite-volume claim follows from these counts.

The script binds the chamber, original compressed bank, bank metadata and producer sources by SHA-256. It has no temporary-file or project-arithmetic dependency. It selects an independent physical covector frame exactly, reconstructs the three new linear functionals at 100 decimal digits, and then evaluates the saved numerical proposal coordinates using binary64 arithmetic. The three outputs are the pairing with normalized wall 28, the current escape and the current strict center. The JSON keeps all 152 cutter recipes with their original bank indices, plus hashes of the full acceptance masks; it does not duplicate the 68,111-row bank.

The added-wall gate accepts a fit to a finite label 2 through 256 within `1e-7`, or a pairing below `-2-1e-7`. It also requires center pairing below `-1e-8`; cutters have escape pairing above `1e-7`. The old 28-wall gates are inherited numerical tests. The originating bank used finitely many anchor labels and sampled its dihedral tails, as recorded in the JSON. Neither this bank nor its rescore is a complete all-label or all-parameter search, and these floating-point gates are not exact angle certificates.

Two promising rows also admit independently reconstructed exact nine-contact recipes. Let

`F = Q(sqrt(5))`, `D = 6 + (16/5)sqrt(5)`, and `K = F(sqrt(D))`.

Both use the eight-wall support `S = [0,1,3,5,6,7,10,27]`, are right to the first seven entries, and meet wall 2 at label 3. Their contact with wall 27 has label 3 or 5 respectively. The script solves these nine physical linear equations over `K`, chooses their kernel vector `v` with coordinate 9 equal to 1, and computes the unique particular solution `p` perpendicular to `v`. The norm-two solutions are

`r = p ± sqrt(delta) v`, where `delta = (2 - p²)/v²`.

The complete exact coordinates of `p` and `v` are saved. The branch with the minus sign matches each numerical proposal to less than `6e-13` in its old normalized pairings.

| Original bank index | Label at wall 27 | Exact delta | Norm from F of delta | Norm from F of delta/D |
|---:|---:|---|---:|---:|
| 1596 | 3 | `50 + 22 sqrt(5)` | `80` | `-100/19` |
| 1437 | 5 | `245/4 + (111/4)sqrt(5)` | `-395/4` | `1975/304` |

Neither delta is a square in `K`. Indeed, `Norm_F(D)=-76/5`, so `D` is not a square in `F`. If `(x+y sqrt(D))²` belongs to `F`, then `xy=0`; therefore delta can be a square in `K` only if delta or delta/D is a square in `F`. The norm of an `F`-square is a rational square, whereas each of the four displayed norms is not a rational square. This is a field-membership obstruction only: both positive deltas define legitimate further quadratic extensions. Their full new-wall angle gates, actual endpoint geometry and joint compatibility remain unproved.

Generation and saved replay pass:

```sh
python scripts/fresh_h9_after_fourth_rescore_20261015.py
python scripts/fresh_h9_after_fourth_rescore_20261015.py --verify
```

The useful next exact experiment would work in one of these quadratic extensions and certify all 29 old-wall relations and the proposed endpoint. This archive performs no such extension search.

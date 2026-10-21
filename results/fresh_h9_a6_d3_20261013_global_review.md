# Independent review of the A6+D3 full-chamber obstruction

**PASS.** I reviewed the [lane report](fresh_h9_a6_d3_20261013.md), the initial/last-edge construction, the [global transport producer](../scripts/fresh_h9_a6_d3_20261013_global.py), and its [separate Fraction checker](../scripts/fresh_h9_a6_d3_20261013_global_independent.py). Both global `--verify` invocations passed against the frozen sources. This review does not repeat the complete norm enumeration or independently implement PARI's norm/unit algorithms.

The initial chamber argument is sound. The trace identity for a vector with four coefficients in the rational even lattice A6+D3, together with its minimum norm two, makes every norm-two spatial vector rational. Such roots generate the rational lattice, so every integral spatial isometry is rational. The two irreducible root components have different ranks and cannot be exchanged. In A6, an outer element `-w`, with `w` in S7, cannot be a reflection: it would require five negative eigenvalues for the involution w, whereas at most three are possible. The D3 factor has the B3 reflection extension. Thus A6+B3, with30 mirrors, is the full zero-height chamber, rather than merely a displayed subgroup.

The safe primitive norm argument has the correct lattice and primes. The actual dual exponent is28; primitive reflection integrality consequently gives a norm dividing56. The base field's certified principal-ideal and full unit-signature results let totally positive unit factors be removed by square scaling. There are13 dyadic valuations and two distinct degree-two primes over seven, each of exponent at most one. Their product gives the52 saved classes. The report correctly separates these base-field facts from the last degree-eight relative field, whose certified class number is eight. The absolute-norm interface uses `bnfisintnorm` and does not assume that latter class number is one.

For the decisive A5+B3 edge, the fixed-projection formula reduces the residual norm to `Y²-42*beta*X²=42*(q-r0²)`. The source checks the necessary integrality of X and Y on the actual lattice basis. Enlarging to all integral solutions is safe because every solution is subsequently tested in the original lattice. The complete finite-angle profile and norm-orbit replay, including previous/next orbit endpoints, supplies firstness. The40-power proposal loop only supplies an upper endpoint and is not asserted to exhaust solutions. The global theorem depends on this certified first edge—71 profiles, two interval points and one primitive first mirror—not merely on the13-wall prefix's volume test.

The global matrix is newly reconstructed from the two oriented full corner frames. It sends the initial A6 chain `(0,1,2,3,4,5)` to `(1,2,3,4,5,12)`, fixes the D3 coordinates, and sends the initial future timelike point to a positive multiple of the endpoint. Both engines check the Lorentz form, inverse, integrality on the actual A6+D3 module in both directions, positive simple-normal scales, primitive reflective normals, and both30-mirror closures. Retained-eight-wall orthogonalities and the strict arrival/release signs reconnect the corners along the complete certified first edge.

Consequently these are two sectors of the same full reflection chamber. An actual lattice isometry permutes the complete root arrangement; carrying one full oriented sector to the other therefore preserves that chamber. Its trace exceeds ten, so it has infinite order. The arithmetic isometry group is discrete because the form is definite at every other real embedding. An infinite chamber stabilizer supplies infinitely many disjoint translates of a sufficiently small interior ball, forcing infinite volume. This conclusion is stronger than, and logically independent of, the saved finite-prefix escape.

The determinant comparison is also appropriately scoped: in even rank ten, scalar similarity changes the determinant by a square. The saved nonsquare ratios distinguish the named rational quadratic spaces. They do not classify all lattices or establish worldwide novelty.

No issue was found. The theorem excludes finite volume for the full root reflection group of this exact host `⟨-(T-1)(T+2)⟩+A6+D3`, with `T⁴-4T²+2=0`. It does not exclude arbitrary compact H9 polytopes or other hosts. The exact edge census and the certified PARI norm/unit interface remain explicit dependencies.

Fresh replays:

```sh
PYTHONPATH=src:. python scripts/fresh_h9_a6_d3_20261013_global.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_a6_d3_20261013_global_independent.py --verify
```

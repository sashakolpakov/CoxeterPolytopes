# Independent review of the all-real quartic-pencil cap obstruction

**PASS.** The [independent checker](../scripts/fresh_h20_quartic_rational_cap_gate_20261014_review.py) re-enumerates the rational endpoint constraints using the ANP algebraic-field arithmetic from our preceding independent review. It imports no producer arithmetic. Its [certificate](fresh_h20_quartic_rational_cap_gate_20261014_review.json) binds the frozen [primary certificate](fresh_h20_quartic_rational_cap_gate_20261014.json), source files and [separate seed-geometry audit](fresh_h20_quartic_orbit_cap_gate_20261014_review.md). Generation and saved replay pass.

The reviewed conclusion is: **for local spatial direction 0, every real shift `u>0` is incompatible with a finite-volume Coxeter completion that retains the 22 cusp facets and uses projections of rational ambient vectors for all added facets.** Ambient norms are unrestricted. Other spatial directions, genuinely nonrational ambient cap directions, and replacement of seed facets remain outside this result. Rational ambient directions alone are not an assertion of arithmeticity.

The independent replay reproduces all 21 rational-norm polynomial calculations, all twelve recorded necessary cases, the three anchor survivor records, and all three identities forcing a rational shift. It also checks the field premise in a different form: the quartic minimal polynomial remains irreducible over `Q(sqrt5)`. Its compositum with that quadratic field therefore has degree eight, so `sqrt5` is absent from the quartic field. The producer instead factors `x²-5` over the quartic field. Both checks use SymPy as a CAS, but the producer's field operations use its own Fraction quotient-ring implementation.

The following points are essential to completeness.

1. **The endpoint constraint is independent of the shift.** The seed triangle retains label3 at all three pairs. Its positive-semidefinite three-root endpoint block forces either one right contact and a second label at most6, or two label3 contacts. A right contact fixes `T=q*g_i/a_i` in the quartic field. For two label3 contacts, the positive kernel fixes their length ratio and gives `T=q*(k_i*g_i-k_j*g_j)/(k_i*a_i-k_j*a_j)`; all denominators are nonzero. Thus `T` belongs to that field even when the shift is transcendental. Only after this observation is it valid to exclude label5 using field membership.

2. **Every rational normalization is covered.** For a selected nonzero old pairing, rational rescaling makes it1; zero cases are handled separately. Both orientations are retained by requiring the two projected pairings to have a common sign, rather than assuming normalization preserves the outward orientation. The later test only requires nonzero cusp pairing, another safe relaxation. The family with both contacts right and free positive norm is included.

3. **There is no hidden continuous rational-norm branch.** Each nonright angle determines the ambient norm as a quadratic polynomial in the remaining rational pairing over the quartic field. The checker reconstructs that polynomial, then requires all three nonconstant field coordinates to vanish. It independently computes their rational gcd and obtains all rational roots from its linear irreducible factors. None of the 21 systems has all three coordinate polynomials zero. Constants, zero roots and the case where the first pairing vanishes are explicitly retained.

4. **Singular identities are retained.** Equality of the four field coordinates gives a rational linear system for the unused old pairing and the cusp pairing. The review uses augmented-matrix ranks and reduced row-echelon form. It checks consistency, every archived parameterized solution, and whether the cusp pairing is identically zero. In particular it retains the one-dimensional solution families rather than treating a singular system as failure.

At the anchor `u=1/2`, the pairs `{23,200}` and `{50,200}` have no necessary endpoint candidate with nonzero cusp pairing. This does not rule out every individual cap: three literal records, representing two projective contact types, survive on `{23,50}`.

For **rational** `u>0`, substitute `h*=2*u*h` into the normal identity. This is a bijective rescaling of nonzero rational cusp pairings, so the anchor obstruction still blocks those two edge classes, totaling **168** germs.

For **irrational** `u>0`, inspect the angle-stage records on `{23,50}` before imposing the anchor identity. The only possibilities are both contacts right, or both labels3 with normalized pairings `g23=0,g50=1`. The exact identities are respectively

```text
u*h = g200,
u*h = g200 + 2.
```

Without normalization the second right side is `g200+2*g50`. All these pairings are rational and a genuine added facet has `h!=0`. Either identity forces u rational, a contradiction. Hence the **84** germs in this third edge class are unclosable for irrational u. The proof does not assert that the other 168 germs are individually blocked at irrational shifts.

The geometric premises extend to every positive shift: all seed Gram entries are unchanged, the target positive affine kernel sums to `u*L*D`, and the other two kernels remain positive multiples of D. Their orthogonal quotient ranks are still 11, 2 and 6. Thus the link is a bounded product of simplices with 252 actual edge germs. Its full-support positive relations forbid a genuinely new acute spacelike facet through D, so every additional facet is strict there. The previously audited edge-endpoint argument applies unchanged: any finite-volume completion must close every germ at a finite or ideal endpoint, while at least one full edge class cannot acquire an endpoint cap. A finite set of walls cannot repair that failure through a limiting accumulation.

The [primary report](fresh_h20_quartic_rational_cap_gate_20261014.md) states this rational/irrational distinction and the retained-facet scope correctly. No proof-sensitive issue was found.

```sh
python3 scripts/fresh_h20_quartic_rational_cap_gate_20261014_review.py
python3 scripts/fresh_h20_quartic_rational_cap_gate_20261014_review.py --verify
```

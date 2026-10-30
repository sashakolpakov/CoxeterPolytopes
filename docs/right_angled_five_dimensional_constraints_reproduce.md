# Five-dimensional attachment and density constraints: reproduction and lemma map

For current reading, see the [current sharp right-angled manuscript](../arxiv/sharp_dimension_bounds_8_and_4/main.tex) and the
[submission guide](arxiv_submission_editions_20261028.md). The supporting
arguments and recorded verification scopes below remain distinct from those
current submission builds. Standalone draft typesetting routes have been retired.

The supporting argument and four-page PDF translate two separate packets:

1. The frozen [one-layer attachment theorem](../results/ideal_ra5_positive_carrier_one_layer_20261018.md) constructs an infinite-volume right-angled H5 chamber having the positive-kappa sixteen-cell Q as an actual facet. All 457 facets are genuine and all 258 original cusps have bounded rectangular four-box links. This does not solve finite-volume completion.
2. The [density theorem](../results/right_angled_mixed_five_face_density_20261018.md) gives ideal RA5 at least 56 facets and consequent necessary mixed-dimensional face/vertex inequalities. It does not give a dimension exclusion or assert optimality against all published bounds.

## Reproduction

Run from the repository root:

```sh
python scripts/fresh_ideal_ra5_positive_carrier_one_layer_20261018.py --verify
python scripts/check_ideal_ra5_one_layer_saved_20261018.py
python scripts/fresh_right_angled_mixed_five_face_density_20261018.py --verify
```

Python uses only the standard library. Default generation refuses to overwrite the saved JSON; `--verify` reconstructs and compares it. The attachment helper imports no producer and uses reverse-column elimination. No final 457-wall cone enumeration is needed: explicit strict facet witnesses prove genuineness, and an explicit strictly feasible ideal direction proves infinite volume.

The attachment generation `57676/160d1b` and separate fresh/helper process `17764/d89074` passed. Density generation and its separate-process replay passed in `32cd3c`. The four-page PDF's final build passed in `aab5e4`, with no warnings, overfull boxes or undefined references in the final log. Review scopes and final integrity hashes are recorded in the reports and their audit JSONs. Typesetting only renders the written proof.

## Precise dependencies

|Statement|Obligation and status|
|---|---|
|Original sixteen-cell ideal RA4 Q|Frozen exact doubling certificate and its written halfspace-union proof. The 198-facet original geometry is inherited.|
|Opposite cusp wall has form `lambda*v-n`|Internal orthogonal-complement and orientation argument; not an externally assumed new lemma.|
|Simultaneous native wall products|Internal Lorentz identities; checked on all 104,196 pairs of the example.|
|All 457 walls genuine|Internal explicit timelike relative-interior witnesses; each checked against every wall.|
|Original Q is the exact actual facet|Internal restriction of all inequalities to `n`-orthogonal space.|
|258 bounded original cusp links|Internal rank-five and four opposite-pair argument with positive common null sums; all saved links checked.|
|Open ideal patch implies infinite volume|Internal continuity and Klein volume-density argument. No finite-volume completion exclusion.|
|Ideal RA5 sum kappa48|Frozen internal carrier theorem derived from cusp incidences and Euler.|
|Positive carrier has integer delta at least5|Frozen internal second-value and equality classification; not an external axiom or a consequence of this arithmetic script.|
|Ideal RA4 minimum24 facets, equality only24-cell|Published Kolpakov2012 theorem.|
|Distinct ridge neighbors and cusp opposites|Internal actual-hyperplane incidence injection. This is needed for the new56 bound and is not replaced by a diagram count.|
|Faces finite-volume/right-angled; simple at edges|Standard geometric input stated in the cited primary literature. Used to apply the published averages and incidence counts.|
|No compact RA5; general RA5 at least16 facets|Published Potyagailo–Vinberg and Dufour inputs. Dufour’s exact statements were read through Alexandrov/Egorov; the subscription original proof was not independently reviewed.|
|Nikulin–Khovanskii bound|Published external theorem for `0 <= l < k <= ceil(n/2)` and polytopes simple at edges. Only `(k,l)=(5,4),(5,1)` in dimensions9–12 are used.|
|Density fractions and finite/ideal vertex ratios|Internal double counting and rational elimination, checked exactly by the small density producer.|

The one-layer lemma is valid without a field restriction. The rational example does not establish an arithmetic lattice. The new density theorem and the universal attachment remain written geometric proofs with internal mathematical reviews; this package makes no full Lean or external human-review claim. Formalization should prove the Lorentz identities, actual facet/cusp arguments, incidence injections and algebra internally while exposing the published NK/Dufour hypotheses explicitly.

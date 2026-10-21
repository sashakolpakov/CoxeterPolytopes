# Restricted RA5 cusp compatibility: reproduction and lemma map

The [manuscript](../manuscripts/ideal_ra5_cusp_compatibility.tex) and [three-page PDF](../build/ideal_ra5_cusp_compatibility/ideal_ra5_cusp_compatibility.pdf) translate the frozen [coupled cusp proof](../results/ideal_ra5_known_carrier_cusp_compatibility_20261018.md). **Every theorem in that paper assumes the explicit three-type four-facet spectrum**: the regular24-cell, its reflected facet double, and the particular sixteen-cell T from the frozen [positive-kappa construction](../results/ideal_ra4_positive_kappa_doubling_20261018.md). It does not classify all ideal RA4 facets or exclude unrestricted RA5.

The distinct [ideal6 paper](../manuscripts/ideal_right_angled_six_face_constraints.tex) contains universal necessary bounds without this three-type restriction. The two scopes must not be conflated.

## Commands

```sh
python scripts/fresh_ideal_ra5_known_carrier_cusp_compatibility_20261018.py --verify
python scripts/check_ideal_ra5_carrier_cusps_saved_20261018.py
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build/ideal_ra5_cusp_compatibility manuscripts/ideal_ra5_cusp_compatibility.tex
```

Generation and separate-process replay passed in `46215/f8ef0a`. The final standalone checker passed in `5302/ffa62e`. Integrity-only audit `1e815a` verifies four input bindings, three report links and both Python ASTs. The audit JSON retains exact source/data/helper/report hashes. The three-page manuscript built in `82493/f956f5`, with no warnings, overfull boxes or undefined references. Python uses only standard-library exact fractions, and default generation refuses to overwrite data.

The producer reconstructs all324 actual cusp stars and all four finite necessary-box domains. Negative cases use complete backtracking; the independent helper imports no producer, reconstructs the reverse coordinate mappings, and uses arc consistency plus complete branching. Positive output is just a labeling checked against every necessary constraint. Neither path solves for a five-dimensional chamber.

## Lemma map

|Statement|Status and precise obligation|
|---|---|
|Actual regular24-cell, double and T|Frozen exact construction certificates and written halfspace-union completeness proofs.|
|Cusp width proportional to opposite-normal sum coefficient|Internal upper-half-space unit-normal calculation, not an externally assumed new theorem.|
|324 complete stars and ridge-size histograms|Fresh exact reconstruction from all inherited normals/vertices, independently checked.|
|Complete ambient four-box shape list|Internal triple-projection argument; all possible length pairs and triples occur.|
|Carrier sum48|Frozen universal internal theorem, valid for any ideal RA5 boundary.|
|Divisibility6 and double-count relation|Internal actual cusp occurrence arithmetic, with no orientation multiplier.|
|22 distinct T neighbors|Internal actual facet intersection argument plus the explicit restricted facet spectrum.|
|No six-clique|Internal Lorentz positive-index-five argument for mutually orthogonal facet normals.|
|Weighted clique edge bound|Internal finite-simplex maximization and nonadjacent-weight transfer proof.|
|At least30 T and66 double facets|Internal consequence of the preceding count and orthogonality lemmas.|
|Complete four-case local domain|Internal width restriction and actual cusp-record permutations/flips. Permitted transformations form a necessary overdomain; no global facet symmetry is assumed.|
|Three excluded local cases|Complete exact necessary-constraint searches, replayed by a different search mechanism.|
|One surviving label and457/2397 aggregate assignment|Explicit necessary-data witnesses only; no Gram, manifold, boundary-complex or chamber realization.|

The universal proof and producer/helper were fully read by the root agent, without execution or reproof of the frozen polytopes/carrier theorem. The package makes no full Lean or external human-review claim. A formal translation must keep the actual facet/cusp geometry as internal proof obligations and must retain the three-type hypothesis in every conditional counting theorem.

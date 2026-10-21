# Right-angled five-face carriers: reproduction and lemma map

The [four-page manuscript](../manuscripts/right_angled_five_face_carriers.tex) contains two separate results:

1. Every entirely ideal right-angled hyperbolic five-polytope has facet sum $\sum_Q(3f_0(Q)-4f_3(Q)+24)=48$. Combined with the internally proved RA4 second-value and equality theorems, some facet has volume at least $3\pi^2$. There are also necessary ideal6 and mixed-dimensional incidence inequalities.
2. An exact sixteen-cell reflected union is an ideal right-angled four-polytope with f-vector $(258,1032,972,198)$ and positive $\kappa=6$. Positive kappa is unbounded under further doubling. This disproves the suggested universal nonpositive-kappa route; it is not an ideal5 example or a dimension bound.

The frozen [carrier proof](../results/right_angled_five_face_carrier_20261018.md) precedes the [positive example](../results/ideal_ra4_positive_kappa_doubling_20261018.md). Its final research suggestion about nonpositive kappa is historical: the successor construction refutes it. The identity and its proved corollaries remain valid.

## Commands

Run from the repository root:

```sh
python scripts/fresh_right_angled_five_face_carrier_20261018.py --verify
python scripts/fresh_ideal_ra4_positive_kappa_doubling_20261018.py --verify
python scripts/check_ideal_ra4_positive_kappa_saved_20261018.py
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build/right_angled_five_face_carriers manuscripts/right_angled_five_face_carriers.tex
```

Python uses only the standard library. Default generation refuses to overwrite saved data. The standalone doubling checker imports no producer; it uses reverse-column elimination for rank checks. Typesetting does not verify the mathematical proof.

The carrier generation and separate-process replay passed in `b132c2`. The positive example generation `6ba641`, fresh replay `1d5a60`, and independent saved checker `5738ef` all passed. The four-page manuscript final build passed in `43d1cc`; its final log has no warnings, overfull boxes or undefined references. Frozen reports retain precise final review and integrity status.

## Published, inherited and internal obligations

|Statement|Status and dependency|
|---|---|
|Ideal RA4 volume formula and regular 24-cell|Published Kolpakov2012 input.|
|Second distinct value and equality classification|Internal universal proofs in the separate frozen packets; not external axioms or conclusions of the identity script.|
|Cusp link is a rectangular box; finite vertices simple|Standard right-angled geometric input. The actual local incidence counts are derived in the carrier report.|
|Three-face identity $f_2=f_0+2$|Euler and four-valent ideal RA3 vertices.|
|Relative-interior coface multiplicity|Internal use of the orthogonal normal diagram: a four-face has exactly $n-4$ five-cofaces.|
|Ideal5 sum kappa48|Internal exact elimination of the cusp and Euler equations. Four integer linear-combination certificates are saved.|
|Positive carrier beyond the two classified types|Internal application of the second-value/equality theorem and integer delta.|
|Ideal6 and mixed inequalities|Internal double counting, with finite-vertex and positive-carrier alternatives kept separate.|
|Complete original 24-cell and double normal/vertex sets|Inherited exact construction certificate; the final 198-facet basis census is not rerun.|
|Facet-double halfspace-union identity|Internal Lorentz reflection and acuteness argument. It transports actual completeness and finite volume.|
|No unrecorded cross-copy wall coincidences|Internal sign argument: $m=s_n\ell$ and both products with $n$ nonpositive force orthogonality and $m=\ell$.|
|Exact sixteen-cell positive example|New rational reflection, complete transformed sets, genuine facets, strict centers, cusp links and merging checks.|
|First positive whole-doubling level16|Bounded family statement using the two facet types of the first double; not global volume or polytope minimality.|
|Unbounded positive kappa|Internal induction using $b\ge6$ and $\kappa'=2\kappa+b-8$.|

The scripts do not formalize the universal geometric statements, and this package makes no end-to-end Lean claim. Any Lean translation should prove the elimination, recurrence and finite incidence arithmetic internally while exposing the actual geometric coface and doubling lemmas as still-to-be-formalized internal obligations.

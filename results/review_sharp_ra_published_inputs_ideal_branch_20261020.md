# Sharp right-angled paper: ideal-branch published-input audit

20 October 2026. **PASS for the cited scopes; no gap found.** I first
read the current TeX, SHA-256
`3c3ddd8540593004bcd5fb57832afc77ffc639268e66cb0ab29302b286654310`,
then checked the original journal PDFs below. No scientific, Lean or TeX
job ran; no canonical source changed. Only this small report was saved;
no downloaded PDF copies or other build products were created.

|Current-paper input|Original primary source and precise location|Verified scope/application|
|---|---|---|
|Musin spherical bound|O. R. Musin, *The kissing number in four dimensions*, Annals of Mathematics 168 (2008), 1–32, [journal PDF](https://annals.math.princeton.edu/wp-content/uploads/annals-v168-n1-p01.pdf): spherical-code definition p. 2; Main theorem and proof p. 4, Section 2. DOI [10.4007/annals.2008.168.1](https://doi.org/10.4007/annals.2008.168.1).|The Main theorem gives k(4)=24. The source identifies the kissing problem with finite subsets of S³ having angular separation at least π/3, equivalently distinct-pair inner products ≤1/2. Thus the manuscript's proposition has the right sphere/ambient dimensions and weak inequality. Its stricter code qualifies immediately. No uniqueness assertion is required.|
|Kolpakov facet count and metric equality|A. Kolpakov, *On the optimality of the ideal right-angled 24-cell*, Algebraic & Geometric Topology 12 (2012), 1941–1960, [journal PDF, p. 1945](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-p.pdf#page=5): Theorem 4, Section 4. DOI [10.2140/agt.2012.12.1941](https://doi.org/10.2140/agt.2012.12.1941).|Theorem 4 directly supplies f₃(P)≥24 and says every four-dimensional ideal right-angled P with f₃(P)=24 is isometric to C. Section 2, p. 1942, defines ideal as all vertices ideal and f₃ as the number of facets. This is metric rigidity, beyond combinatorial recognition.|
|Actual finite-volume right-angled H⁸ example|L. Potyagailo and E. B. Vinberg, *On right-angled reflection groups in hyperbolic spaces*, Commentarii Mathematici Helvetici 80 (2005), 63–73, [publisher PDF](https://ems.press/content/serial-article-files/43032): Section 3, pp. 66–67; Proposition 3.1, p. 67. DOI [10.4171/CMH/4](https://doi.org/10.4171/CMH/4).|The construction uses finite orbits of finite-volume simplices in dimensions 3–8 to give actual right-angled finite-volume Pₙ. Proposition 3.1 explicitly includes P₈, with 240 facets and 2160 cusps. This supports the general bound's attainment; entire ideality is not claimed for this example.|
|Actual entirely ideal H⁴ example|[Kolpakov journal PDF](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-p.pdf), introduction p. 1941 and Section 2 p. 1942.|C is explicitly the regular ideal hyperbolic 24-cell with right dihedral angles and 24 octahedral facets. It supplies the entirely ideal bound's attained dimension 4.|

The facet/volume distinction is correctly handled. Kolpakov's **Theorem 3,
p. 1944**, is the separate minimum-volume result. The manuscript cites
**Theorem 4, p. 1945** for facet number and its equality case, both of which
appear explicitly in that theorem. Its proof, Section 4.3, pp. 1949–1958,
ends with the isometry conclusion for 24 facets. No volume-to-facet inference
is needed. [Original theorem and proof](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-p.pdf#page=5).

In the manuscript's ideal proof, the projected four-facet has a code to
which Musin applies. Its upper facet count and Kolpakov's lower count
meet at 24; Kolpakov then gives the actual metric 24-cell needed by the
manuscript's opposite-normal height calculation. Uniqueness of a kissing
configuration is neither supplied nor used. The cited page/theorem labels
and bibliographic data match the published sources.

This audit checks those published statements, definitions and the relevant
proof endpoints against their use. It does not re-referee the full cited
proofs or audit the separate Dufour/Vinberg inputs. No source correction,
new numerical result or end-to-end formal-verification claim follows.

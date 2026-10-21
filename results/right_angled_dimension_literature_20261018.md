# Right-angled dimension bounds: primary-source audit

Checked online on 18 October 2026. This is a literature/status note, not a new proof or research certificate. “Ideal” means **all vertices ideal**; “mixed” means **at least one finite and at least one ideal vertex**. Noncompact finite volume includes both ideal and mixed cases.

| Right-angled class | Published dimension upper bound | Equivalent nonexistence statement | Highest known example dimension | Repository improvement |
| --- | ---: | --- | ---: | --- |
| Compact | 4 | None for `n >= 5` | 4: compact right-angled 120-cell | None; published bound remains stronger |
| All-ideal, finite volume | 6 | None for `n >= 7` | 4: ideal right-angled 24-cell | None; published bound remains stronger |
| Mixed finite and ideal vertices, finite volume | 12 | None for `n >= 13` | 8: `P^8` in the Potyagailo–Vinberg family | None; published bound remains stronger |
| Finite volume, without a vertex-type restriction | 12 | None for `n >= 13` | 8 | None; published bound remains stronger |

The compact bound is sharp. The ideal bounds leave dimensions 5 and 6 open, and the mixed/general finite-volume bounds leave dimensions 9 through 12 open in the checked literature. An upper bound of 12 does **not** assert that examples exist in dimension 12.

1. **Compact bound and sharpness.** Potyagailo–Vinberg, *On right-angled reflection groups in hyperbolic spaces*, Comment. Math. Helv. 80 (2005), 63–73, §2, gives a proof of the known Vinberg proposition excluding compact right-angled polyhedra above dimension four. Its immediately following remarks identify the compact regular right-angled 120-cell in dimension four. The publisher lists publication on 31 March 2005. [Publisher record](https://ems.press/journals/cmh/articles/457), [full paper, printed p.65](https://ems.press/content/serial-article-files/43032).

2. **Finite-volume bound.** Guillaume Dufour, *Notes on right-angled Coxeter polyhedra in hyperbolic spaces*, Geometriae Dedicata 147 (2010), 277–282, proves the upper bound 12. The primary publisher abstract states the dimension bound directly. Published online 20 January 2010; August 2010 issue. The full article was subscription-only through the browser, so this audit verified the statement from its primary abstract and the later original papers below, not by reading Dufour's entire proof. [Publisher abstract and bibliographic record](https://link.springer.com/article/10.1007/s10711-009-9454-2).

3. **All-ideal bound and dimension-four example.** Alexander Kolpakov, *On the optimality of the ideal right-angled 24-cell*, Algebraic & Geometric Topology 12 (2012), 1941–1960. Corollary 7 on printed p.1959 excludes ideal right-angled hyperbolic polytopes in every dimension `n >= 7`. The introduction describes the ideal right-angled 24-cell in dimension four; the paper proves its optimal volume/facet properties. [Full primary paper](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-s.pdf).

4. **Dimension-eight examples really are mixed.** Potyagailo–Vinberg §3 constructs finite-volume right-angled `P^n` for `3 <= n <= 8`; Proposition 3.1 lists `P^8` with 240 facets and 2,160 cusps. The original construction paper by Giovanni Italiano, Bruno Martelli and Matteo Migliorini, *Hyperbolic manifolds that fibre algebraically up to dimension 8*, explicitly states in its “Outline of the construction” that these `P^n` have both ideal and real vertices. It appeared online 10 November 2022 and in J. Inst. Math. Jussieu 23(2) (2024), 609–646. This prevents treating the dimension-eight examples as all-ideal. [Primary paper](https://doi.org/10.1017/S1474748022000536).

5. **Recent status checks.** Stepan Alexandrov's original 2023 preprint, §1.1, gives new proofs of the known bounds and explicitly records examples only through dimension four in the ideal case and dimension eight for finite volume, with sharpness unknown outside the compact case. Andrey Egorov's 31 March 2026 original preprint, §1, repeats the three nonexistence thresholds and the dimension-eight example limit while improving vertex-count estimates. Its abstract's loose description of 12 as a maximal dimension should not be read as an existence theorem: the introduction explicitly distinguishes upper bound 12 from examples through 8. [Alexandrov 2023](https://arxiv.org/pdf/2303.09533), [Egorov 2026](https://arxiv.org/pdf/2603.29925).

The two preprints in item 5 are used as current-status checks, not mislabeled as peer-reviewed publications. The dimension bounds in the table already have the published sources in items 1–3. No stronger right-angled dimension theorem or higher-dimensional right-angled example was found in the targeted current searches. No claim of exhaustive bibliographic coverage is made.

## Inputs for a distinct right-angled research lane

Egorov's 2026 preprint claims the following necessary vertex counts (Theorems 3.1, 3.2, 4.1, 4.2):

| Dimension | Ideal vertices at least | Finite vertices at least |
| --- | ---: | ---: |
| 5 | 3 | — |
| 6 | 9 | — |
| 7 | 35 | 4 |
| 8 | 205 | 41 |
| 9 | 3,690 | 1,312 |

Its dimension-nine calculation also gives at least 288 facets. These are lower bounds, not exact minima or existence claims. With `F` facets and facet-level minima `c,f`, its incidence recurrences are `v_ideal >= ceil(F*c/[2(n-1)])` and `v_finite >= ceil(F*f/n)`. They count the cubical ideal link and simple finite link respectively. This audit checked the statements, not the new gluing proofs. Before using the preprint as a theorem input, perform a separate proof review. [Egorov, §§3–4](https://arxiv.org/pdf/2603.29925).

Alexandrov's 2023 paper supplies a complementary starting point: weighted finite/ideal vertex incidence strengthens the usual average-face inequality and yields short proofs of the known bounds. It is useful for locating the exact inequality deficit in a lower target dimension; reproducing its existing exclusions would not improve the table. [Alexandrov, §§3–6](https://arxiv.org/pdf/2303.09533).

Concrete proposed tasks, not results:

- **All-ideal dimension six.** Fix the whole cubical five-dimensional link at a cusp, then derive complete compatibility rules for adjacent ideal stars and their actual three- and four-faces. Every four-face is an ideal right-angled hyperbolic four-polytope, so the published 24-facet minimum is available; the equality case is the 24-cell. Audit the equality and next-incidence cases with exact Gram rank/signature and actual-face recovery, seeking an inequality stronger than the current all-ideal face average. A complete dimension-six contradiction would improve the ideal upper bound from six to five. A bounded facet search alone would not do so.
- **All-ideal dimension five.** Use the cubical four-dimensional cusp link and shared ideal four-faces as a smaller geometric pilot. The immediate target is a new forced-incidence or bounded-fiber inequality, or a completely verified example. If a contradiction covered all such polytopes, faces would also exclude dimension six and make the ideal upper bound four, sharp. Restricting to 24-cell facets is only a subclass unless an independent argument forces that restriction.
- **Mixed dimension nine.** All finite-volume right-angled examples in dimension nine would necessarily be mixed, by the compact and ideal exclusions. Start with whole finite orthogonal nine-stars and cubical eight-dimensional cusp links; test paired actual-face/coface incidences, then combine their exact counts with weighted finite/ideal inequalities. The preprint's large necessary counts favor a local incidence certificate over an unstructured full-facet enumeration. A complete exclusion in dimension nine would propagate through faces and improve the finite-volume upper bound to eight; a finite-volume example would improve the known-example dimension instead. Neither outcome follows from vertex lower bounds alone.

For all three tasks, right-angled Gram entries at meeting facets are fixed, but distances between disjoint facets remain variables. Restricting those variables to `Q(sqrt(5))` or another arithmetic field would address only that restricted model, not the full right-angled dimension problem. Any future exact certificate must retain this distinction and prove actual face incidences and occurrence multiplicities before exporting counts. No new right-angled bound, completed search or Lean formalization is claimed here.

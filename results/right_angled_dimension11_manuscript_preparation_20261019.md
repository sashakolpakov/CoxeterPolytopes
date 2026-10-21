# Focused mixed dimension-eleven manuscript: source handoff

19 October 2026. **Finished source, awaiting full source review and a
coordinator-scheduled build.** No TeX build, scientific checker/import,
Lean execution or numerical job was run by the author of this source.
No central documentation, catalogue, predecessor manuscript/PDF or
frozen proof packet was changed.

Primary source: `manuscripts/right_angled_dimension11.tex`, SHA256
`6c455e9e36420581a1987d4fab8e24e2bfa745c8c4648b4c81ff73f89576ddcf`.
It has32189 bytes,759 source lines and three inline TikZ vector figures.
Estimated typeset length is10--16 pages; that is not a build result.
Static label/reference/citation, environment-stack and brace checks passed on the initial mathematically identical source. Root subsequently polished numeral spacing only and independently reviewed the complete source; the final TeX is pinned above.
There is no new PDF yet. The existing mixed-corner and opposition-transport
papers remain unmodified supporting predecessors in the same paper family.

## Theorem and source-proof map

The main theorem is mixed dimension<=11. A separate corollary applies
the published compact<=4 and entirely-ideal<=6 facts to get the same
upper bound for arbitrary finite-volume right-angled polytopes. The
new proof itself excludes dimension12. The earlier Dufour<=12 theorem
is the explicit higher-dimensional tail. No dimension11 example or
dimension10 exclusion is asserted.

|Source section|Complete argument included|
|---|---|
|1, scope|Actual Lorentz/face conventions; published12/4/6 tails and82-facet input; precise f,c,B counts.|
|2, actual three-face oppositions|Two-normal acute projection/product gap; bad-triangle wall; P3 opposite and every-other-wall feasibility; exact two-finite/three-ideal hull; distinct triples and undirected B-edge graph; h invariance; one fixed future null Z on each component.|
|3, formal triples|No feasibility claim for absent triples; complete positive ratio domain; subset inclusion probabilities and all four moment identities with cubic sign.|
|4, quartic bound|Both P>=0 and reciprocal majorization; convex lower parabola; exact remainder identity; full degree-six expansion; all17 Bernstein coefficient bounds and central derivative chain on the entire interval.|
|5, global upper|Positive-part actual-subset deletion, ordinary unweighted vertex sum, reciprocal actual-edge pairing, exact1420561/12500 coefficient.|
|6, strict weighted polygon row|Complete simple-truncation convention, DS generating polynomial and two coefficient rows, positive-dimensional coface counts, published82 and elementary11-facet incidence bounds; strict a2>15f+4246c.|
|7, small faces and exception budget|Euler/area identities; own F6/nu9 P3 recognition; full ideal-edge triangle attachment proof; Q8 zero-excess/quadrilateral structure and exactly one ideal edge at each finite corner; k=0..12 exceptional-triple bound; full nu-case proof of p budget.|
|8, elimination|Actual E/R and three-face Euler incidences;101B coefficient row; two strict lower inequalities; positive cusp cancellation; exact5861741 separation from transport upper.|
|Final verification record|Exact scalar generation/replay status, conditional standard-library Lean final count logic and its actual premise boundary; no claim of end-to-end geometric formalization.|

The figures show(1) the actual P3 incidence object,(2) the rational
interval proof domains, and(3) formal-row deletion versus actual
reciprocal edge pairing. None is presented as a metric projection,
polynomial approximation or proposed ambient polytope.

## Exact proof and verification pins

The principal written proof and its independent complete review are

* `results/right_angled_mixed_positive_part_opposition_cap_research_20261019.md`,
  SHA256 `19ec3d814759cd2d9b9b1df8ef128f1df7d7afee68590aa3310479c954a65c9f`;
* `results/review_mixed_positive_part_opposition_cap_resource_20261019.md`,
  SHA256 `358969924c4ce3681084ee38755468e0917cf3fc2cdda5b9747dcd33758abf56`.

The independent review rederived all geometry, moment/majorant steps,
the entire inherited lower bridge and strict comparison, without
executing the scientific checker or reading the root's review first.
The source-bound supporting18October proofs are pinned in that review.
Their necessary new arguments have been integrated into the manuscript,
rather than left as unnamed local-lemma assumptions.

The root separately scheduled exact generation and read-only replay:

* `scripts/check_right_angled_mixed_positive_part_opposition_cap_20261019.py`,
  SHA256 `50627366713d0a8f2417fc0ed7251a62a8572a702b79a5389d3a2130eb663892`;
* `results/right_angled_mixed_positive_part_opposition_cap_research_20261019.json`,
  SHA256 `2b1e08dd3b2ce25efec25e5d3947ba581ed8ef3f06582454d8776296fb13c7e1`.

The settled conditional final logic is
`formal/CoxeterProofs/RightAngledMixed11.lean`, SHA256
`f0b1cffede6210dbddf7cc7d95f29d14baa5984bcd2e19641e029afb1e7f864d`.
Root reports compilation and separate printed-axiom audit PASS, with
only propext/Quot.sound and an axiom-free exact comparison. It imports
only Std. Its real transport/majorant upper interface, actual incidence
data and published higher tail remain premises, not a formalization
of the geometric proof now included in this source.

## Preserved parallel work

The separately authorized arithmetic19 integration has now resumed:
`manuscripts/appendices/compact_arithmetic_eight_field19.tex` is included
by the same primary `manuscripts/compact_arithmetic20.tex`. Its source
handoff is recorded separately in
`results/compact_arithmetic_eight_field19_manuscript_preparation_20261019.md`;
its existing compiled PDF remains unchanged. Before that integration,
the old primary source was byte-preserved at
`manuscripts/source_archives/compact_arithmetic20_before_eight_field19_20261019.tex.gz`
with verified decompression identity. The main coordinator owns all
central catalogue/documentation updates and future serial builds.

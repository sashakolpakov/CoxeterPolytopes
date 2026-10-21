# Ideal5/6 primary-manuscript integration, source-only

19 October 2026. The SAME primary source,
[ideal_right_angled_six_face_constraints.tex](../manuscripts/ideal_right_angled_six_face_constraints.tex),
now presents the current independently reviewed necessary bounds:

| Class | Facets | Vertices |
|---|---:|---:|
| Entirely ideal right-angled dimension5 | ≥63 | ≥211 |
| Entirely ideal right-angled dimension6 | ≥275 | ≥6106 |

The dimension6 relation is `272V≥6039F`. For mixed ambient dimensions11/12,
the ideal five-face fraction is `<4/47` and the ideal six-face fraction
is `<15/248`, conditional on the precisely stated published NK/Dufour
premises. No dimension ceiling or existence claim is made.

The source has866 lines,39,288 bytes and four inline TikZ vector figures.
Its preparation SHA256 is
`bf05cb7b634bc4e5c89822774d32f2c18bbefe24338737d41cef0aa67003184f`.
No TeX, Lean, checker, census, optimizer or scientific body ran during
this integration. Page count and actual figure layout await a serial
build and visual review; the intended scale is12–16 pages, without
truncating the human arguments to reach a target.

## Source preservation and stale PDF

The catalogue's `ideal6` entry reported no explicit TeX freeze. The
reviewed transport/reciprocal freeze is
`results/ideal_ra6_cusp_transport_reviewed_20261018_freeze.json`; it binds
the older scientific PY/JSON/MD packets and reviews, not this primary TeX.
None of those packets, their PDFs, or the companion ideal4 source was
edited. The central bounds and result-correspondence documents were not
edited by this integration agent.

Before replacing the primary source, its old bytes were archived at
[the source archive](../manuscripts/source_archives/ideal_right_angled_six_face_constraints_preintegration_20261019.tex.gz).
The old source was12,795 bytes, SHA256
`26744b042f8919e0687f5070e2cd0ca341f8c87f54c8755be940ffc5b714563d`.
The gzip is5,576 bytes; streaming decompression and SHA256 comparison
proved it byte-identical to the original, whose hash was checked again
immediately before atomic replacement. The temporary replacement path
was consumed; it was not retained as clutter.

The existing [PDF](../manuscripts/build/ideal_right_angled_six_face_constraints.pdf)
was preserved. It is298,518 bytes, SHA256
`efbd6e1fb230aec6e1a0c65b515f4d594569fcdec5c4af62f966abc1c8ae56e1`,
and is **stale**: it records the older five-page bounds. The archived
source reproduces that historical text. No revised PDF or new compilation
is claimed here. The existing reproduction page and catalogue still
describe the predecessor package until the coordinator updates them;
they must not be read as execution evidence for the revised source.

## Theorem and dependency map

| Primary label / section | Full content and dependency |
|---|---|
| `lem:projection`, §1 | Actual all-wall acute projection, nonright wall-product gap, exclusion of mutually orthogonal finite-vertex walls. The foundational acute-normal and rectangular-cusp facts are explicit. |
| `lem:injection`, §1 | Distinct ridge neighbors and actual opposite-cusp facet injection. |
| `lem:capacity`, §2 | Complete hemisphere/box proof that an outside cusp has at most two facet-neighbors; no untested hemisphere apex is assumed feasible. |
| `lem:polygon`, `lem:quad`, §3 | Fixed null scales, spinor recurrence including positive wrap signs, two-small-corner rank proof, triangle coefficient1, and exact `8/B²` quadrilateral identity. |
| `lem:partition`, §4 | Exact `N1=2q` via the actual N2/triangle bijection, shared for facets of dimensions3/4. |
| `lem:goodwitness`, §4 | One common proof of nonorthogonality, hemisphere tangency, convexity against all walls, and the actual edge; `G≤2N0` and `9G≤2N0`. The separate extension without goodness requires independently proved nonorthogonality. |
| `thm:plus3`, §5 | Complete polygon-corner count and integer defect `V(Q)≥3V(R)+3`; no withdrawn +6 claim. |
| `lem:octasquare`, `thm:eleven`, §6 | Complete equality exclusion: saturation, six-vertex opposite facets, own rectangle/hemisphere square-cusp proof, exact quad coefficients, two4-of5 quad-facet sets, and only then the all11-facet witness contradiction. |
| `thm:fiveprop`, §7 | Complete `3V(P)≥13V(Q)+6F(Q)`, with six cusp facets/three coefficient pairs and four corner occurrences each. |
| `prop:bridge`, §8 | Explicitly cited INTERNAL companion theorem: through10 facet spectrum implies an opposite-restricted24-cell chain, hence nonpositive κ. Its external3D classification and realization uniqueness inputs are named. The arbitrary four-polytope is not assumed to be a garland. |
| `lem:carrier`, §8 | Human cusp/Euler proof `Σκ=48`, positive carrier, bridge to ≥11, local equality exclusion to37 vertices, and SEPARATE injection to25 facets. |
| `thm:fivebounds`, §8 | `F5≥63,V5≥211,F6≥275` by injection and the complete propagation inequality. |
| `lem:goodedge`, §9 | Actual triangle wall, feasible endpoint, inverse widths, noncubic matching bound and missing-link-edge count. |
| §9 and `lem:reciprocal`, §10 | Canonical cubic invariant, positive formal energy ratios, compact/all-moment balancing proof of the sharp103/132 average including boundary/equality cases, and unweighted reciprocal-pair edge summation. |
| `thm:sixbounds`, §10 | Complete face identities, retained triangle factor3, universal136/99 transport,63-facet input, ratio6039/272 and6106 rounding. |
| §11 | Published NK range checked for(5,4),(6,5) in11/12; general minima16/27 and noncompact/compact scope explicit; fractions4/47,15/248 proved. |
| §12 | Mathematical review, small script and conditional Lean scope are distinguished; no new execution or full geometric Lean claim. |

The valid weaker208 extra-credit argument remains preserved in
`results/right_angled_ideal5_carrier_equality_strictness_20261019.md`.
It is not needed by the stronger211 proof and was omitted from the
primary narrative. Historical intermediate constants are likewise not
presented as successive main results.

## Vector figures

1. `fig:box`: box-edge contact capacity and the actual parallel-wall
   tangency section, used by both witness proofs.
2. `fig:quad`: fixed-scale quadrilateral coefficients and the exact
   opposite-null product, rather than an unlabelled combinatorial picture.
3. `fig:equality`: six-good/three-N0 incidence saturation, followed by
   the4-of5 obstruction; its caption separates the good-only and all-facet
   stages of the proof.
4. `fig:transport`: actual inverse-width good-edge transport, invariant
   cubic component and the reciprocal-pair summation.

All figures are repo-native TikZ; no raster assets or external rendering
steps were introduced.

## Written evidence and verification boundary

The primary human proofs translate the already reviewed reports:

* `right_angled_outside_cusp_capacity_20261019.md`;
* `right_angled_ideal4_cusp_defect_research_20261019.md`;
* `right_angled_ideal5_facet_cusp_propagation_20261019.md`;
* `right_angled_ideal4_eleven_vertex_equality_exclusion_20261019.md`;
* `right_angled_five_face_carrier_20261018.md`;
* `ideal_ra6_cusp_transport_research_20261018.md`;
* `ideal_ra6_sharp_reciprocal_20261018.md`.

The new equality source is pinned at
`b9af9b497c69da69a25faced7067b854f0e3b8332f89da57daa12a5a968c2121`;
its independent root and finite57 mathematical reviews PASS. The
small-facet bridge additionally uses
`ideal_ra4_ten_vertex_facet_exclusion_research_20261019.md` and the frozen
unified ideal4 structural proof. These are internal dependencies,
not newly declared external axioms.

Existing exact frame, reciprocal and carrier scripts retain their
historical inputs and stated coverage. No finite frame list is treated
as a proof for all widths. Older conditional Lean results do not imply
formal coverage of the new equality geometry or the current constants.
The coordinator must record any subsequent checking, source review,
typesetting, visual inspection or formal extension in separate ledgers.

The only local checks in this preparation were source-metadata/static
checks: no control characters, balanced braces/environments and dollar
delimiters,52 unique labels, no undefined local references, and four
figure environments. These are not a TeX parse, build or mathematical
checker. The coordinator owns the next resource-gated serial build.

Root fully read the preceding865-line draft and reported mathematical
PASS subject to a face-inheritance scope clarification. That requested
clarification is now applied: faces of dimension at least three inherit
right angles, whereas two-faces in the entirely ideal class are ideal
polygons. The second requested edit keeps execution/typesetting records
separate rather than embedding a stale-PDF notice in the mathematical
text. The final866-line bytes await root's handoff/build; the PDF-stale
status above remains the actual preparation status.

## New-ledger builder source, not executed

Prepared
`scripts/build_illustrated_ra_manuscript_with_ledger_20261019.py`, SHA256
`ae521225307be682f120860c1e0a90279b206b24191a084be19ccedec0c12c8e`.
It is an exact source copy of the stable20261018 illustrated builder
except for the required `--ledger results/NEW.json` argument, validation
that the ledger is a new direct JSON child of `results`, and exclusive
`open('x')` creation before PDF replacement. Existing ledger records
cannot be replaced even with `--replace-unfrozen`.

The59s/1GiB caps, shared lock, resource checks, parent-owner check,
output bounds, stop-group cleanup, temporary directory handling and
watchdog loop are unchanged. The prior builder was preserved at SHA256
`27f14501b52a38c92121266aceb58eb2e1d36bfc21c2615ab0b4e9be4f9dee15`.
The new source passed AST parsing only; neither its main function nor a
TeX engine was run. Root owns source review and any scheduled execution.

Root also archived the historical PDF at
`manuscripts/source_archives/ideal_right_angled_six_face_constraints_preintegration_20261019.pdf.gz`
(292,587 bytes). A fresh bounded streaming decompression/hash check here
confirmed the old PDF SHA256 `efbd6e1fb230aec6e1a0c65b515f4d594569fcdec5c4af62f966abc1c8ae56e1`.
Thus both historical source and PDF bytes remain recoverable without
overwriting their scientific proof packets.

## Coordinator completion, superseding preparation-only status

The final primary source is782cea54d88970b70e1d6692934f48f5727a46b262bf6b49e342bd9cd68493ad.
The bounded final build in`results/ideal56_figure_spacing_build_20261019.json`
passes:12 pages, four figures, two passes and zero warnings. All four
figure pages were inspected; Figure4 received one spacing correction
and its final rendering was rechecked. The old source/PDF archives remain
unchanged. The current conditional Lean successor also compiles and
passes its separate axiom audit; see
[right_angled_eleven_vertex_counts_verification_20261019.md](right_angled_eleven_vertex_counts_verification_20261019.md).
Actual hyperbolic geometry remains explicit written premises.

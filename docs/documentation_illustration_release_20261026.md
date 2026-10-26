# Documentation and illustrations, 26 October 2026

This revision keeps six paper families and changes no dimension ceiling. All 14
tracked READMEs are shortened; detailed explanations are in the
[Sphinx site](https://sashakolpakov.github.io/CoxeterPolytopes/).

The complete illustrated reading copies have 50 pages (finite volume 57),
286 pages (general compact 17) and 217 pages (compact arithmetic 15).
Each adds two TikZ figures. Removing all marked illustration blocks restores
the reviewed consolidated source byte for byte. The six frozen editions,
three consolidated bases and three right-angled reading copies are preserved.
The sharp right-angled, ideal RA4 and compact RA4 papers already have one,
five and four figures respectively.

## What was audited

Independent source review checked each new figure against the exact proposition
it illustrates. Root review checked the source changes and all six newly rendered
figure pages. Build receipts identify the final source and PDF hashes. None has
an overfull box or undefined reference; inherited underfull spacing warnings
remain. This is sampled layout inspection, not an all-page visual audit of the
long appendices. The [visual record](../results/illustrated_manuscripts_six_new_figures_visual_acceptance_20261026.json)
and [independent review](../results/review_illustrated_manuscripts_diagram_semantics_20261026.md)
retain those boundaries.

The Sphinx book explains common arguments, case splits, finite domains,
termination and completeness, computed propositions, and conditional Lean
interfaces. Five original SVG diagrams show these dependencies. The short
READMEs point into the book rather than duplicating the full evidence inventory.
Current theorem scopes were checked against the reviewed sources and their
existing maps. This revision does not rerun successful scientific checks or
Lean modules merely to refresh documentation.

## Historical records and current documentation

The old right-angled release record has 116 recursively present path/hash
bindings: 111 still match current files, while five occurrences bind the former
`formal/README.md`. That long guide has a concise successor and a full dated
text in [the formal records](formal_verification_records_20261025.md). The exact
old README bytes survive in the byte-verified preservation archive.

The third-party input manifest retains its original 37-file payload hashes and
historical provenance-README hashes. The nine current provenance READMEs are
concise successors, not claims that the old hashes match their new text.
The [successor manifest](../results/readme_documentation_provenance_successor_20261026.json)
binds all 14 new READMEs to their preserved originals. Old manifests, reviews
and execution receipts remain immutable historical evidence.

## Publication checks and availability

The Pages workflow audits explicit local Sphinx inputs, transitive source
includes and all static assets before building with warnings as errors. The
input guard rejects manifest-listed unuploaded scientific files and excluded
third-party originals by path and exact byte hash. It also checks source/index
membership in CI. It is a publication check, not mathematical verification or
an analysis of arbitrary extensions and modified/compressed unknown payloads.
The actual Git index is audited separately before pushing.

No bulk asset is uploaded by this revision. The manifest still describes 936
scientific files in 33 additional archives. Original third-party papers/data,
local compiler caches and generated sites are not included. A clean source
checkout is not a promise that every historical replay can run without separately
obtained inputs. The detailed site distinguishes these availability limits from
logical assumptions and actual execution credit.

# Explicit prism replacement diagnostic

The two algebraic branches from the [one-cap obstruction](fresh_h20_opposite_prism_gate_20261017.md) bound a compact target H3 prism only after **142 of the 168 old walls are deleted**. The resulting 28-wall branch is a separate diagnostic. It does not repair the retained168-wall construction and has no established whole-H20 volume conclusion.

The [source](../scripts/fresh_h20_prism_replacement_20261017.py) and [saved data](fresh_h20_prism_replacement_20261017.json) retain the exact core recipes, then check every old native pairing numerically. In the notation of the exact obstruction, the two norm-two normals are

\[
r_\pm=p\pm\sqrt{(2-p^2)/u^2}\,u,
\qquad p^2=(-5+\sqrt2)/23,\quad u^2=1/30.
\]

Both are right to T,159,161 and have product -1 with160. Their mutual product is exactly2p^2-2<-2. The exact recipes do not turn the remaining numerical old-wall tests into an exact ambient certificate.

| Retained branch | Deleted old walls | Total walls | Target H3 finite vertices | Target H3 exterior vertices | Maximum fixed-chart squared radius |
|---|---:|---:|---:|---:|---:|
| Both new branches |142|28|6|0|0.9999799758428216|
| Branch0 only |142|27|3|1|1.0006340642367981|
| Branch1 only |1: wall158|168|10|14|2.2284093476008735|

The complete142-wall deletion list and both contacts responsible for each deletion are saved. Exactly26 old walls survive the two-cap branch: the17 roots of T, the three partners159,160,161, and the six cusp-completion walls162 through167. The compact prism is therefore not evidence that a small modification closes the original chamber.

Branch1 is a different, smaller modification: it replaces only158, retaining the other167 old walls. It removes the old E8 attachment of the triangular base because the replacement is right to T. Its target face still has14 exterior vertices and the same largest squared radius as the old target. This supplies a possible alternative base for future coordinated finite-contact work, not a successful closure.

Each target cone is intersected with one fixed Klein affine chart by a strict-center linear program and a full numerical halfspace vertex census. Coincident numerical vertices are deduplicated at distance1e-7; the raw count is also saved. The finite/exterior test uses a1e-6 distance from squared radius1. There are no near-sphere vertices in these three outputs. Old native gate classification uses1e-7 tolerances; no all-label or exact sign claim is inferred from those tolerances. Exactification would still need to certify every retained gate, coherent ambient interior, actual facets, and any claimed complete ambient census.

Generation and a fresh saved replay pass. The literal source, source hashes, and local report links are checked. This diagnostic is final/frozen; its scope remains numerical.

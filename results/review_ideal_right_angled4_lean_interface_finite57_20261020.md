# Independent source review of the ideal-four conditional Lean interface

20 October 2026. Full source read only; no Lean invocation or axiom audit
was run in this lane.

Reviewed `formal/CoxeterProofs/IdealRightAngled4.lean`, SHA256
`3aad406da960989e2fe7a1e7a5e60c700243ba011d9427e333571ae32ebc4fad`.

**PASS for the stated conditional logic and assumption boundary.** This
does not extend the formalized geometry beyond the source's explicit
claims, and compilation remains the coordinating agent's separate job.

`Height` is an arbitrary type. There is no integer-height proxy or hidden
discretization. The existence and positivity of a real maximum are not
proved: geometry supplies a maximum value, an attaining actual facet, a
zero-height cusp facet, and their distinctness. This is sufficient for the
equality-propagation argument the kernel is intended to check.

The published inputs explicitly contain Musin's cardinality bound,
Kolpakov's facet lower bound, and Kolpakov's **metric** equality case. They
do not assume uniqueness of a kissing configuration. The actual-foot,
normalized-bound, spherical-code and center predicates are opaque to Lean;
their local constructions and implications are fields of `WrittenGeometry`.
In particular all-wall acuteness, Lorentz normalization, injection of
neighbors into the code, opposite-normal metric rigidity and the final
neighbor-height formula remain actual written geometric inputs.

No field assumes that all facets are maximal or that a five-polytope is
impossible. Those conclusions are obtained by the cardinality squeeze,
metric equality, local neighbor closure and induction along the explicitly
supplied adjacency paths. The final contradiction compares the actual
chosen cusp facet to the nonzero maximum. The higher-dimensional theorem
exposes actual five-face inheritance as a hypothesis producing the full
local certificate; it does not purport to formalize the face lattice.

The code imports only `Std`, includes no `sorry` or new `axiom`
declaration, and asks for audits of its eight public conclusions. This is
source inspection, not a substitute for the separately scheduled actual
kernel compilation and printed-axiom audit.

# Known-example provenance

This directory contains exact, machine-readable Gram matrices, root lists,
diagrams, and source locators. An entry is usable only to the extent stated by
its `artifact_status`; source-audited Gram data are not thereby an exhaustive
Vinberg reconstruction.

The manifest distinguishes exact polyhedral geometry from root-search
completeness.  The Bugaenko artifact is marked `verified`: a clean exact
Vinberg run reconstructs the complete eleven-wall system from the lattice and
control vector.  The separate chamber verifier checks the Coxeter Gram data
and Vinberg's endpoint criterion: all 164 elliptic rank-seven edge diagrams
have exactly two elliptic rank-eight endpoint extensions, while no parabolic
rank-seven cusp diagram occurs.  This certifies compactness; the enumeration
of all 41 finite vertices is corroborating geometry, not by itself the volume
argument.  Direct access to the 1992 chapter remains a provenance desideratum,
not a mathematical dependency of the reproduced certificate.

The `P_{6,10}` entry is deliberately split at the theorem boundary. Its
existence, compactness, exact Gram, and index-two relation to Bugaenko's
nine-wall `H6` chamber are frozen; its mixed norm-two/norm-four golden root
scaling is also exact. The uniqueness claim of arXiv:2608.02894v1 is not
imported because its 387-to-11 reduction uses the false general assertion
that every facet of a compact Coxeter polytope is Coxeter. The complete
source and lineage audit, including the `[5,3,5]` facet counterexample and the
hidden parent-wall warning for surface carriers, is
`results/p6_10_bugaenko_double_lineage_and_classification_audit.md`.

The five cases in `k3_high_rank_nef_chambers.json` have status
`source_audited_gram_verified`. They transcribe the complete published
negative-curve matrices
for full K3 nef chambers of Picard ranks 15 through 19, giving known
noncompact finite-volume Coxeter polytopes in hyperbolic dimensions 14 through
18. The verifier checks every entry, exact rank and inertia, basis determinant,
lattice discriminant, and the discriminant factorization and generic fibre
configuration of the derived rank-19 elliptic equation. Its generic Picard-rank
conclusion additionally uses the Noether--Lefschetz argument recorded in the
report. Provenance and one corrected redraw ambiguity are documented in
`results/k3_high_rank_nef_chambers.md`.

`enriques_finite_automorphism_nef_chambers.json` is the exact seven-image
transfer certificate for complex Enriques surfaces with finite automorphism
group.  Kondō--Martin supply the iff classification by the complete graph of
all smooth rational curves; Shimada supplies the exact `L10` wall vectors and
intersection matrices.  The project independently replays all seven matrices,
enumerates all 2,508 projective cone rays with GMP rational arithmetic, and
applies Vinberg's full edge-endpoint and affine-zero-neighborhood test.  The
result is two 12-facet and five 20-facet full-nef Coxeter polytopes in `H9`,
all finite-volume noncompact.  This is a strong realization/classification of
known configurations, not a claim of worldwide polytope novelty.  A separate
no-import hostile certifier reconstructs the same 2,508 rays and every
Vinberg gate directly from Shimada's raw arrays; its exact certificate is
`data/enumerations/enriques_finite_aut_nef_chambers_hostile_audit.json`.

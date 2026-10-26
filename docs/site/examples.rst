Known examples and construction records
========================================

The construction lane has not produced a new hyperbolic Coxeter polytope.
It has reconstructed known examples, transferred published surface data to
exact wall systems, and recorded obstructions to proposed constructions.
These achievements have a different role from the dimension upper bounds.
A necessary diagram list is not an example, and an exact finite wall Gram
matrix alone does not prove that a root-search algorithm found every wall.

The machine-readable starting point is
:download:`the example index <../../data/known_examples/index.json>`.
Each artifact's status records its actual scope. The concise
:download:`directory README <../../data/known_examples/README.md>` links the
main certificates; :doc:`sources` documents their third-party inputs.

Three questions for an alleged example
--------------------------------------

A useful certificate separates three propositions.

#. **Exact local data:** roots, scalar products, labels, rank and Lorentz
   signature agree exactly with the proposed diagram.
#. **The chamber:** the chosen halfspaces have the claimed finite or ideal
   vertices and satisfy an appropriate finite-volume/compactness criterion.
#. **Search completeness or geometric identification:** a stated root search
   found every necessary wall, or a published geometric theorem identifies
   the chamber as the intended full nef or reflection chamber.

A signature check addresses the first question. Enumerating some vertices
can support the second, but does not by itself prove the endpoint criterion.
A source-transcribed wall matrix does not become an exhaustive Vinberg run.
All three questions require their own witnesses and interpretation.

Bugaenko's compact eight-dimensional chamber
---------------------------------------------

The :download:`H8 artifact <../../data/known_examples/bugaenko_compact_h8.json>`
and :download:`Vinberg transcript <../../data/known_examples/bugaenko_vinberg_transcript.json>`
record an exact eleven-wall reconstruction from the lattice and control vector.
The chamber verifier separately checks the Coxeter Gram data and the complete
edge-endpoint criterion: all 164 elliptic rank-seven edge diagrams have exactly
two elliptic rank-eight endpoint extensions, and no parabolic rank-seven cusp
diagram occurs. The 41 finite vertices corroborate that geometry; their count
alone is not the compactness proof.

The certificate concerns a known example. Access to the historical 1992 chapter
remains a provenance desideratum, but the recorded reconstructed certificate
has its own explicit mathematical input boundary. No dimension-nine compact
example is inferred from the eight-dimensional reconstruction.

The P6,10 lineage record
------------------------

The :download:`exact P6,10 certificate <../../data/known_examples/p6_10_exact_certificate.json>`
keeps existence, compactness, exact Gram data and its index-two relation to
Bugaenko's nine-wall H6 chamber separate from a broader uniqueness claim.
Its mixed norm-two/norm-four golden scaling is exact data, not permission to
identify all root scales or all chamber facets.

The :download:`source and lineage audit <../../results/p6_10_bugaenko_double_lineage_and_classification_audit.md>`
explains why the project does not import the claimed classification reduction.
It records the facet-heredity issue and an explicit facet counterexample.
This is a statement about the pinned source argument examined by that audit;
it is not a fresh survey of every subsequent version of the literature.

K3 nef chambers in dimensions fourteen through eighteen
-------------------------------------------------------

Five :download:`source-audited chamber matrices <../../data/known_examples/k3_high_rank_nef_chambers.json>`
transcribe published complete negative-curve configurations for Picard ranks 15
through 19. Their hyperbolic chamber dimensions are 14 through 18. The recorded
status is ``source_audited_gram_verified``: exact entry, rank, inertia,
basis-determinant and lattice-discriminant checks are distinguished from a
root-by-root reconstruction of an entire reflection system.

The :download:`source report <../../results/k3_high_rank_nef_chambers.md>`
records a corrected redraw ambiguity and the separate Noether–Lefschetz input
used for the generic Picard-rank conclusion of the rank 19 elliptic model.
The :download:`volume record <../../data/known_examples/k3_high_rank_nef_vinberg_volume.json>`
retains its own endpoint and geometric interpretation. None of these records
makes a known configuration a new example.

Seven Enriques chamber transfers
--------------------------------

The :download:`finite-automorphism Enriques certificate <../../data/known_examples/enriques_finite_automorphism_nef_chambers.json>`
uses Kondō–Martin's classification by the complete graph of smooth rational
curves and Shimada's exact L10 wall/configuration arrays. The project records
exact reconstruction of all seven matrices, 2508 projective cone rays, and
Vinberg edge-endpoint plus affine-zero-neighborhood tests. This gives two
12-facet and five 20-facet known finite-volume noncompact chambers in H9.

The independent no-import certifier reconstructs the same rays and gate data
from the raw arrays. Its saved certificate,
``data/enumerations/enriques_finite_aut_nef_chambers_hostile_audit.json``,
is a scientific bulk input rather than a promise that a source checkout
contains that file. See :doc:`sources` for its actual availability. The
important distinction is between a strong realization/classification of known
configurations and worldwide novelty of a polytope.

Other chamber stocks and parked searches
----------------------------------------

The directory also retains Apéry–Fermi, K3 quotient and Jacobian–Kummer stocks.
Their individual statuses determine whether an entry is a source-audited Gram,
a chamber candidate, an exact geometric certificate or a search intermediate.
The project does not promote the existence of such a file to an assertion of
finite covolume or complete wall enumeration.

Historical searches can require large third-party supplements and are not part
of the six paper families' selected verification plans. Reproducing one begins
with its own frozen source, input manifest and algorithm contract, rather than
an unbounded invocation copied from an old note. The current work adds
explanation and illustrations; it does not restart those searches.

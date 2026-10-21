# Compact-H9 incidence repair: exact-family review

The [exact checker](../scripts/fresh_h9_incidence_repair_20261014_verify.py)
passed a fresh root-agent saved-result replay against its
[certificate](fresh_h9_incidence_repair_20261014_verify.json). Its source
and the following two mathematical arguments were reviewed. This is a
source/mathematics review and producer replay, not an independent
implementation of the whole enumeration.

For the72 starting right/dotted patterns, the finite-only15-node matrix
is positive definite. Every clique avoiding dotted pairs therefore has a
positive-definite Gram block, regardless of dotted distances. Changing
additional old right pairs to dotted pairs only deletes cliques, leaving
a flag nerve. A compact H9 chamber would have a simplicial8-sphere nerve,
and any flag8-sphere has at least18 vertices. Indeed, each of the nine
ridges of one maximal simplex needs a second maximal simplex; their nine
outside vertices are distinct, since a repeat would create a10-clique.
Thus the15-facet family is impossible for arbitrary extra right-to-dotted
changes and arbitrary allowed dotted distances.

The other exact family changes seven originally right pairs to labels3
or4, preserves old finite labels, and leaves original parallel pairs
dotted with arbitrary distances. The checker exhausts the eleven minimum
seven-pair covers of the relevant old affine circuits and all128 labelings
per cover. Its spherical-block predicate is exact: a nonright cycle
cannot be positive definite, and each remaining tree admits a rational
diagonal rescaling with root norms2 or4, checked by Fraction LDL pivots.
The rooted weighted-tree cache identifies only Gram blocks with the same
normalized diagram.

For each of the1,408 profiles, every potential nine-wall vertex is removed
by iterative unique-ridge peeling. An eight-wall ridge in a compact simple
H9 polytope has two endpoints. A candidate vertex with no second potential
endpoint cannot survive; deleting it may expose more such defects. Extra
right-to-dotted changes only remove potential vertices, so cannot restore
a closed nerve after complete peeling. This argument does not choose
dotted distances or assume numerical rank tests are exact.

The review does not certify the separate SAT or numerical metric scouts.
Changing more finite contacts, using other labels, replacing old parallel
incidences by finite contacts, adding facets, or changing the source
remains outside these two exact family exclusions.

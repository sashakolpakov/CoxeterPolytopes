# Full-norm Vinberg search and compact-star contact extensions

**The discriminant-11 arithmetic host is now proved nonreflective.** An exact
infinite-order chamber symmetry rules out every finite-volume Coxeter polytope
whose wall reflections preserve this lattice. The separate combinatorial
search advances through new contact-family exclusions. No new finite-volume
polytope or dimension has been established.

## Arithmetic construction

The explicit discriminant-11 host from the
[previous checkpoint](disjoint_cap_edge_extensions_and_star_symmetry_2026_10_14.md)
is isometric to the repository's older H11 host (the name refers to its
discriminant, not hyperbolic dimension 11; this is an H9 search). The older investigation
enumerated norm-two reflections; its accepted walls are not a prefix of the
full reflection group. The complete primitive norm alphabet is **2 and
2(3+phi)**. The new search exhausts both classes in increasing geometric height,
starting from the [complete E8+A1 zero sector](fresh_h9_beta11_zero_sector_20261014.md).

The [height-100 enumeration](fresh_h9_beta11_vinberg_20261014.md) has 21 walls.
The [height-300 continuation](fresh_h9_beta11_vinberg_h300_20261014.md) has
35 walls, after all 80 positive levels and 1,082,417 enumeration nodes.
[Independent height-100](fresh_h9_beta11_vinberg_20261014_review.md) and
[height-300 reviews](fresh_h9_beta11_vinberg_h300_20261014_review.md) reproduce
the complete root censuses using separate arithmetic and a reversed search
order. Every wall has a checked relative interior point.

The [21-wall geometry](fresh_h9_beta11_geometry_20261014.md) has 260 finite
vertices. A new long root removes its escape. Nevertheless the
[35-wall geometry](fresh_h9_beta11_geometry_h300_20261014.md), with 880 finite
vertices, still has infinite volume. Its escape and the earlier escape have a
[separate arithmetic review](fresh_h9_beta11_escape_20261014_review.md).

The decisive step is an
[integral corner-sector isometry](fresh_h9_beta11_corner_symmetries_20261014.md)
with determinant 1 and trace **4+8phi > 10**. It sends the initial corner to
finite vertex 735, mapping all nine outward incident roots exactly. Every
finite corner of an acute Vinberg prefix persists in the full chamber:
later walls cannot cut or touch it. Thus this is a symmetry of the full
chamber, at all heights. Its infinite order proves that chamber has infinite
volume. The [independent proof review](fresh_h9_beta11_corner_symmetries_20261014_review.md)
checks the matrix, inverse, signatures, orientation, corner signs and physical
preservation decomposition with different arithmetic.

This settles the full norm-two/norm-2(3+phi) group, not just the old norm-two
subgroup. Do not continue integral reflection searches in this host. The
[height-1000 experiment was stopped](fresh_h9_beta11_vinberg_h1000_20261014_stopped.md)
after this decisive proof; it supplies no completed prefix or additional
result used here.

A [universal mixed-norm restriction](fresh_h9_beta11_mixed_norm_gap_20261014.md)
shows that a norm-two wall and a norm-2(3+phi) wall must be orthogonal or
ultraparallel. Its [independent review](fresh_h9_beta11_mixed_norm_gap_20261014_review.md)
checks the number-theoretic argument, including the parallel boundary case.

The [literature and lattice comparison](fresh_h9_beta11_literature_20261014.md)
finds no verified classification excluding this host. McLeod's related
diagonal-lattice nonreflectivity result does not transfer automatically across
the dyadic lattice change. Indeed, his transported infinite-order isometry is
a [verified product of 94 reflections already present here](fresh_h9_beta11_mcleod_word_20261014.md).
It supplies no chamber-stabilizer obstruction for this lattice.

## Combinatorial search

The [mixed-single contact extension](fresh_h9_star_mixed_singles_20261014.md)
allows every noncore own contact to vary in its specified star family.
Its [independent review](fresh_h9_star_mixed_singles_20261014_review.md) checks
448 new core witnesses, 27 core guards, three budgets and the high-label
argument; the 7,840-case affine arithmetic remains an explicitly inherited
dependency. This first extension advanced the search to 1,271 visits.

The [opposite-single core theorem](fresh_h9_star_opposite_core_20261014.md)
then removes the noncore hypotheses from the next family. Its
[independent review](fresh_h9_star_opposite_core_20261014_review.md) checks
3,612 new core witnesses, while the 47,040-case polynomial bank remains an
explicitly inherited, previously reviewed calculation. Both core theorems
transfer to [twelve literal incidence placements](fresh_h9_star_core_transfers_20261014.md).
The [separate transfer review](fresh_h9_star_core_transfers_20261014_review.md)
checks 26,624 induced-subset decisions, 48 predicates, 108 maximal-face
witnesses and 4,123 occurrences of source-proof witnesses. No global
automorphism of the full incidence is assumed.

The combined binder now has **253 reviewed predicates**. Its continuation
visits 19,566 additional states and is saved at **20,837 cumulative visits,
with nine pending frames**, in
[batch 011](fresh_h9_star_core_transfers_20261014_search_batch_011.json).
The new local transfers alone prune 3,150 partial states.

The next assignment has internal labels **(4,at least 7,at least 7,at least 7)**.
The first positive leaf has own contacts (3,3), and its cut leaf 14 has (3,4).
Each of the other three positive/cut leaf pairs has shared single contacts
(2,3). The third leaf's changed (3,4) profile lies outside the proved (3,3)
core theorems. This is a necessary combinatorial state, not a geometric
realization; omitted contacts and the Lorentzian metric are unresolved.

## Verification and continuation

All frozen predecessor artifacts remain unchanged. The
[archive audit](full_norm_vinberg_and_star_contacts_2026_10_14_audit.json)
checks recursive source hashes, report links and preservation of the previous
checkpoint. Archive integrity is separate from mathematical verification.

Primary checks, independent reviews and saved replays pass. The final
combinatorial traversal replay takes 89.812 seconds and preserves all nine
pending frames. The archive audit binds 22 reports, 129 local links and 299
unique source files through 163 recursive manifests, preserving the frozen
predecessor. The stopped height-1000 experiment remains explicitly incomplete.

```sh
python scripts/fresh_h9_beta11_corner_symmetries_20261014_review.py --verify
python scripts/fresh_h9_star_core_transfers_20261014_search.py --batch 11 --verify
python scripts/audit_full_norm_star_contacts_checkpoint_20261014.py --verify
```

Preserve the arithmetic certificates and move to different hosts or walls
that do not preserve this lattice. In the compact-star search, preserve the
saved cursor and address its next unresolved core rather than restarting a
bounded label search. General nonintegral Coxeter constructions remain outside
the host-wide exclusion.

# Intersecting block families add four exact packing constraints

Four vertex types violate valid clique capacities even after every subset
integer-packing inequality has been imposed. The repaired exact witness is

\[
 \frac{2790968068837079040}{217745752925063023}
 =12.817554562350466\ldots>12.
\]

Thus [all-subset floor inequalities](fresh_sqrt5_all_subset_port_packing_20261017.md)
do not exhaust the useful local packing constraints in this model. The new
cuts still do not remove its obstruction. No optimizer or dimension
exclusion is claimed.

At a finite vertex, actual maximal-simplex port blocks are disjoint.
Consequently a family of possible blocks in which every two intersect
contains at most one actual block. Its corner indicator has vertex cost 1.
To make this invariant under the labelled-forest automorphisms, take a union
of complete block orbits. Every orbit must be internally pairwise
intersecting, and every literal block in one selected orbit must intersect
every literal block in each other selected orbit. These universal conditions
are essential; a pair of representative blocks is insufficient.

The [source](../scripts/fresh_sqrt5_intersecting_block_cliques_20261017.py)
tests all 6,011 positive block orbits at all 740 positive-simplex forests.
Exactly 3,862 orbits are internally intersecting. At each forest it builds
the graph of universal cross-orbit intersection and enumerates every maximal
clique, 2,002 in total. With nonnegative masses every clique lies in a maximal
one of at least the same load. Zero-load orbits cannot improve the maximum.
Disjoint witness pairs are recorded for every rejected orbit and every
missing cross-orbit edge.

There are exactly four violating maximal cliques, at vertex indices 879,
904, 1451 and 2948. Their selected block-orbit indices are respectively
`[342,2325]`, `[1287,1484,4436]`, `[427,874,2181,4685]` and
`[168,1312,1319,3206,3207]`. The
[saved result](fresh_sqrt5_intersecting_block_cliques_20261017.json)
exports their integer coefficients on the positive simplex rows, together
with every tested clique graph and its maximal cliques. All four capacities
equal 1. Quadrilateral rows receive zero payment from this new family.

The face masses stay fixed and only the four necessary vertex masses
increase. Hence every previously checked positive-cost inequality remains
valid and all signed transfer loads are unchanged. The numerator is the same
as in the full-model and all-subset certificates; the larger normalizer gives
the ratio above.

Generation and a fresh saved replay both terminated with exit 0 and
`PASS_EXACT_ALL_POSITIVE_INTERSECTING_ORBIT_CLIQUES`. The construction agent
separately reviewed the proof, literal-orbit interpretation, universal
intersection tests and monotone repair, with no issue found. This is a narrow
source and mathematical review, not a second implementation of maximal-clique
enumeration. The enumeration uses NetworkX's standard `find_cliques` routine.

This exhausts the stated invariant clique family at this fixed witness.
General joint packing inequalities, new face classes, other fields and an
ambient dimension theorem remain separate. **FINAL/FROZEN.**

```sh
python3 scripts/fresh_sqrt5_intersecting_block_cliques_20261017.py
python3 scripts/fresh_sqrt5_intersecting_block_cliques_20261017.py --verify
```

# TOOT: a 112-regular finite-subtree factory in `H12`

Status: **exact theorem; hostile audit PASS / FREEZE** (2026-08-08); see the
[independent raw-Gram audit](jacobian_kummer_h12_maximal_disjoint_doubling_family_hostile_audit.md).

## Positive theorem

Let `Q=JK-H12-172` be the certified 172-facet finite-volume noncompact
Coxeter polytope in `H12`.  Among its exact 124 doubling walls there is a
unique largest labelled set of pairwise disjoint walls:

```text
W61,W62,...,W172                         (112 walls).
```

Here **disjoint** has Allcock's precise meaning: the hyperplanes do not meet
in hyperbolic space, so both parallel and ultraparallel pairs qualify.  The
relations inside this set are

```text
  664 parallel pairs,
5,552 ultraparallel pairs,
    0 intersecting pairs.
```

There is also a unique largest labelled set for the strictly stronger
relation “ultraparallel”:

```text
W77,W78,...,W172                          (96 walls),
```

and all `binom(96,2)=4,560` pairs in it are ultraparallel.

Let `R_i` be the reflections in either selected stock.  They generate

```text
*_{i=1}^112 C2        or        *_{i=1}^96 C2,
```

whose Cayley graph in the displayed generators is respectively the
112-regular or 96-regular tree.  For every finite subtree `T`, the union
`Q_T` of the corresponding copies of `Q` is a finite-sided, finite-volume,
noncompact hyperbolic Coxeter polytope in `H12`, with

```text
[W(Q):W(Q_T)] = |V(T)|,
Vol(Q_T)       = |V(T)| Vol(Q).
```

This is a much larger branching construction than the previously isolated
path, cube, and direct-product subfamilies.  It does not provide a closed
facet formula for arbitrary tree shapes.

## Why the maxima and uniqueness are exact

The independent verifier reparses Dolgachev--Shimada's 316 published wall
vectors, reconstructs the saturated trace-nine fixed lattice, and rebuilds
all 172 primitive facet roots.  It imports no Jacobian--Kummer scout or
production certifier and uses no graph library.

All 14,706 facet-pair relations are reclassified from integer Gram data.  The
complete 124-wall doubling-stock relation census is

```text
  264 orthogonal pairs,
1,138 parallel pairs,
6,224 ultraparallel pairs.
```

A deterministic Bron--Kerbosch enumeration gives every maximal clique for
each of the two relations:

| relation | exact maximal-clique size census | unique maximum |
|---|---|---:|
| parallel or ultraparallel | `36^9, 55^36, 74^42, 93^12, 112^1` | 112 |
| strictly ultraparallel | `32^32, 57^12, 61^16, 96^1` | 96 |

There are respectively 100 and 61 maximal cliques.  Every clique and its
maximality are replayed directly.  Uniqueness is not inferred from the
single displayed construction: every doubling wall excluded from the
112-set lies in cliques of size at most 93, and every wall excluded from the
96-set lies in cliques of size at most 61.

Separate upper-bound certificates avoid relying only on clique search.  The
complement of the first graph contains an explicit 12-edge vertex-disjoint
matching, so a clique has size at most `124-12=112`.  The complement of the
second contains an explicit 28-edge matching, giving `124-28=96`.  The full
matchings and all maximal cliques are frozen in the JSON archive.

The words “maximum” and “unique” are scoped to these exact 124 seed doubling
walls.  They do not classify roots outside the seed diagram, reflected
descendant walls, other K3 chambers, or arbitrary redoubling constructions.

## Exact applicability of Allcock's theorems

Allcock calls a wall doubling when every wall it meets does so at an even
submultiple of `pi`.  In this seed, all 124 displayed walls pass exactly:
each is a primitive square-minus-four root, no one has an incident `m=3`
edge, and every reflection is an integral involutive isometry of the full
rank-13 lattice.

The hypotheses of Allcock's Theorem 4.1 are therefore literal:

1. `Q` is a Coxeter polytope;
2. the selected walls are doubling walls; and
3. they are pairwise disjoint.

The theorem makes every subtree union a Coxeter polyhedron.  If `T` is
finite, only finitely many seed chambers and facets occur.  Their interiors
are disjoint and their overlaps have measure zero, proving the displayed
volume formula and finite volume.  Since a translated copy of the
noncompact seed is contained in `Q_T`, every such finite union is
noncompact.  This conclusion does not infer volume from the lattice: the
seed's literal Vinberg certificate has zero bad edges and 6,214 cusps, and
the independent exact cone has zero hyperideal rays.

Allcock's Theorem 4.3 requires a finite-volume seed with at least three
pairwise disjoint doubling walls.  It applies overwhelmingly here; one may
even choose three directly from the 96 ultraparallel walls, so the
preliminary doubling step used in the general proof is unnecessary.  Its
conclusion is that the number of index-`I` reflection subgroups, up to
conjugacy by isometries of `H12`, has an exponential lower bound in `I`.
Equivalently, Allcock's argument selects exponentially many nonisometric
finite-tree polytopes of volume `I Vol(Q)`.  This does **not** say that every
pair of nonisomorphic finite subtrees gives nonisometric polytopes.

## Full-Picard K3 inheritance

The same seven-dimensional very-general K3 family carries every finite-tree
descendant, with

```text
NS(X) = U(2) + D5(-1)^2 + <-4>.
```

This assertion has a fixed boundary-core proof.  After translating `T`, take
it to contain the identity seed chamber.  Only selected generator walls are
ever glued.  Hence all other identity-chamber walls remain exterior walls.
For the 112-generator construction these are exactly `W1,...,W60`; their
primitive coordinates have rank 13 and Smith diagonal

```text
1,1,1,1,1,1,1,1,1,1,1,1,1.
```

Thus they span the full Picard lattice with index one.  For the
96-generator construction, the larger core `W1,...,W76` has the same Smith
diagonal.  Every other descendant boundary normal is the image of a
primitive seed root under integral Picard isometries.  Consequently every
finite-tree boundary consists of primitive integral reflection roots and
contains an index-one Picard-spanning core.

This remains a **weak divisor-reflection carrier**.  The reflections are not
claimed to be surface automorphisms or Hodge isometries; the roots are not
claimed effective; and no descendant chamber is claimed to be the nef cone
of a K3 surface.

## Lineage and claim boundary

The finite-subtree and exponential-growth mechanism is classical, from
Daniel Allcock, *Infinitely many hyperbolic Coxeter groups through dimension
19*, Geometry & Topology 10 (2006), 737--758,
[DOI 10.2140/gt.2006.10.737](https://doi.org/10.2140/gt.2006.10.737),
[arXiv:0903.0138](https://arxiv.org/abs/0903.0138).  In particular, Allcock
already proved infinitely and exponentially many finite-covolume Coxeter
examples in `H12`, using a Conway-face seed.  This theorem is not a first in
dimension 12 and not a new arithmetic commensurability class.

Dolgachev--Shimada supplied the ambient Jacobian--Kummer surface, divisor
lattice, chamber and wall vectors; their source does not state this
trace-nine 172-wall fold or its maximal doubling stocks.  The exact project
gain is the unique `112/96` graph theorem for this surface-derived seed and
the full-Picard inheritance certificate for all finite subtrees.  Whether
these exact maxima or descendant diagrams appeared elsewhere has not been
exhaustively audited, so worldwide priority remains open.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/certify_jacobian_kummer_h12_maximal_disjoint_doubling_family.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_jacobian_kummer_h12_maximal_disjoint_doubling_family.py
```

Artifacts:

- [`jacobian_kummer_h12_maximal_disjoint_doubling_family.json`](../data/enumerations/jacobian_kummer_h12_maximal_disjoint_doubling_family.json)
- [`certify_jacobian_kummer_h12_maximal_disjoint_doubling_family.py`](../scripts/certify_jacobian_kummer_h12_maximal_disjoint_doubling_family.py)
- [`test_jacobian_kummer_h12_maximal_disjoint_doubling_family.py`](../tests/test_jacobian_kummer_h12_maximal_disjoint_doubling_family.py)

Frozen JSON SHA-256:

`47d930ba3ed70d0ba4f9a77a0c10db04605c76fd8e63b5450ba000233e6cd4c4`

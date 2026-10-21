# Hostile audit: the JK-H12 112-wall finite-subtree factory

Status: **PASS / FREEZE** (2026-08-08). No mathematical correction was
required.

## Independent replay

The hostile verifier imports neither the production certifier nor its JSON
output. It reads the frozen 172-by-172 integer facet Gram directly and
reclassifies all 14,706 pairs. It recovers

```text
whole seed:       2820 m2, 144 m3, 800 m4,
                  1750 parallel, 9192 ultraparallel;
doubling stock:   264 m2, 1138 parallel, 6224 ultraparallel;
doubling walls:   124.
```

Maximal cliques are not replayed with the production Bron--Kerbosch code.
The audit instead enumerates maximal independent sets of each conflict graph
by the exact recurrence obtained by separating sets that contain a chosen
vertex from those that omit and dominate it. The result is

| relation | edges | maximal-clique census | unique maximum |
|---|---:|---|---|
| parallel or ultraparallel | 7,362 | `36^9, 55^36, 74^42, 93^12, 112^1` | W61--W172 |
| strictly ultraparallel | 6,224 | `32^32, 57^12, 61^16, 96^1` | W77--W172 |

The complete independently generated clique families have SHA-256 hashes

```text
disjoint:       54fdfc5849c88f8e26368c82bd6f526da82a3dde461094c0a50be3c4433d45bf
ultraparallel:  4e7761d1b5a486e19e8f9d79b079b47254d39170241b45ace36f6df17ead5793
```

The separate conflict matchings replay exactly. Their 12 and 28
vertex-disjoint nonedges prove the independent upper bounds
`124-12=112` and `124-28=96`; hence the enumeration is not the only
maximality witness.

## Primary-theorem transfer

Section 4 of Daniel Allcock,
[*Infinitely many hyperbolic Coxeter groups through dimension
19*](https://doi.org/10.2140/gt.2006.10.737), states that reflections in
`k` pairwise disjoint doubling walls generate the free product of `k` copies
of `C2`, whose Cayley graph is the `k`-regular tree. Theorem 4.1 makes the
union indexed by every subtree a Coxeter polyhedron; Corollary 4.2 gives
index equal to the number of tree vertices; and Theorem 4.3 gives an
exponential lower bound for index-`I` reflection subgroups up to conjugacy by
hyperbolic isometries when the finite-volume seed has at least three such
walls.

The hypotheses hold literally here. The 96-stock already supplies three
pairwise ultraparallel walls, so Allcock's preliminary doublings used in the
general proof of Theorem 4.3 are unnecessary. For every finite subtree `T`,

```text
[W(Q):W(Q_T)] = |V(T)|,
Vol(Q_T)       = |V(T)| Vol(Q).
```

Finite volume follows by additivity over the finitely many seed chambers.
Noncompactness is not guessed from the lattice: the frozen literal Vinberg
certificate has 6,214 cusps and zero failures, the independent exact cone has
zero hyperideal rays, and `Q_T` contains an unbounded isometric copy of the
noncompact seed. Allcock's exponential conclusion does not say that every
pair of nonisomorphic finite subtrees produces nonisometric polytopes.

## Full-Picard inheritance

The audit separately reads the 172 primitive wall coordinates in

```text
NS(X) = U(2) + D5(-1)^2 + <-4>
```

and reconstructs their complete relation matrix from the rank-13 carrier
Gram. It agrees pair-for-pair with the raw facet Gram. Every one of the 124
doubling reflections is independently checked to be an integral involutive
isometry of the full Picard lattice.

Smith reduction gives

```text
W1--W60:  rank 13, Smith 1^13, index 1;
W1--W76:  rank 13, Smith 1^13, index 1.
```

After translating a finite subtree to contain the identity chamber, these
unglued walls remain exterior. Thus every descendant boundary contains a
full-Picard index-one core; every other boundary normal is a primitive root
transported by an integral Picard isometry.

This is only a weak full-Picard divisor-reflection carrier. It proves no nef
equality, effectivity, Hodge extension, or surface-automorphism action.
Maximality is confined to the exact 124 seed doubling walls. No arbitrary-tree
facet formula, new dimension, new commensurability class, or worldwide
priority is claimed; Allcock's exponential `H12` abundance is classical.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/audit_jacobian_kummer_h12_maximal_disjoint_doubling_family_hostile.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_jacobian_kummer_h12_maximal_disjoint_doubling_family_hostile_audit.py
```

Artifacts:

- [`hostile certificate`](../data/enumerations/jacobian_kummer_h12_maximal_disjoint_doubling_family_hostile_audit.json)
- [`independent verifier`](../scripts/audit_jacobian_kummer_h12_maximal_disjoint_doubling_family_hostile.py)
- [`focused tests`](../tests/test_jacobian_kummer_h12_maximal_disjoint_doubling_family_hostile_audit.py)
- [`production theorem`](jacobian_kummer_h12_maximal_disjoint_doubling_family.md)

Frozen hostile JSON SHA-256:

`6df5a785d4682c9afa8ac33e5728a2378205e22aa73a1b646fdf2552d2085d89`

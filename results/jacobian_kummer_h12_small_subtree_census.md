# TOOT: 35 explicit three-chamber Jacobian--Kummer polytopes in `H12`

Status: **exact positive census; independent hostile audit PASS / FREEZE**
(2026-08-08).  Worldwide priority remains open pending a
publication-complete lineage audit.

## The positive result

Let `Q=JK-H12/172` be the certified finite-volume noncompact Coxeter
polytope carried by the rank-13 Picard lattice

\[
 U(2)\oplus D_5(-1)^2\oplus\langle-4\rangle .
\]

Its walls `W61,...,W172` are 112 pairwise disjoint doubling walls in
Allcock's sense: every pair is parallel or ultraparallel.  The complete
small-finite-subtree census in this stock is:

| finite subtree | raw choices | exact diagram types | facet counts | volume multiple |
|---|---:|---:|---|---:|
| one edge / two seed chambers | 112 | 2 | `317` (one type), `326` (one type) | 2 |
| colored path / three seed chambers | 6,216 | **35** | `462` (4 types), `471` (7), `480` (24) | 3 |

All 35 three-chamber diagrams are pairwise nonisomorphic as abstract
labelled Coxeter diagrams.  Thus the census adds three explicit
finite-volume noncompact facet counts to the project ledger:

```text
H12/462, H12/471, H12/480.
```

More precisely, 120 raw generator pairs yield the four `H12/462` types,
1,536 yield the seven `H12/471` types, and 4,560 yield the 24 `H12/480`
types.  This is a classification of every one-edge double and every
three-vertex colored path inside the exact 112-wall stock, not a sample.

The one-edge result sharpens the earlier representative theorem: the
16 walls `W61,...,W76` form one diagram type with 317 facets, while the
96 walls `W77,...,W172` form one type with 326 facets.  Those two examples
were already in the project ledger; the three-chamber facet counts and 35
diagram types were not.

## Exact orbit and nonisomorphism proof

The verifier reparses Dolgachev--Shimada's published 316-wall data.  Their
order-23,040 chamber stabilizer has a trace-nine centralizer of order 768;
the induced faithful action on the 172 folded facets has order 384.  It has
two orbits of sizes 16 and 96 on `W61,...,W172`, and exactly 35 orbits on the
6,216 unordered pairs.  The pair-orbit size census is

```text
8^1, 32^2, 48^4, 96^6, 192^16, 384^6.
```

Six pair orbits are parallel and 29 are ultraparallel.  Every raw pair in
one orbit gives an isometric colored-path union by the seed isometry.
Treating the pair as unordered loses no type: in the free product,
left translation by `ba` carries the chamber set `{1,a,ab}` to
`{1,b,ba}`, reversing the three-vertex path.

Different pair orbits are separated without using the Picard decoration.
For each descendant vertex of the Coxeter diagram, form the histogram of
all exact normalized Coxeter labels incident to it, including the exact
rational normalized square on every ultraparallel edge.  Sort the resulting
multiset of vertex histograms canonically.  The 35 exact multisets are all
different.  This is a facet-permutation invariant of the abstract labelled
Coxeter diagram, so the 35 representatives are pairwise nonisomorphic.  The
JSON freezes a SHA-256 of each full canonical multiset; the script compares
the underlying exact tuples before hashing.

In fact the separation is stronger than the metric labels: all 35 multisets
remain different after first forgetting every ultraparallel distance, and
again after collapsing both parallel and ultraparallel pairs to the single
ordinary Coxeter-matrix entry `infinity`.  Hence the result cannot be an
artifact of Picard root types or fine dotted-edge weights.

The complete representative table is:

| # | walls | pair orbit | generator relation | facets | profile hash |
|---:|---|---:|---|---:|---|
| 1 | `61,62` | 8 | parallel | 462 | `a064112135a657b4` |
| 2 | `61,63` | 48 | parallel | 462 | `290d7472a96b9f19` |
| 3 | `61,69` | 32 | parallel | 462 | `232c8ef1bce8ec35` |
| 4 | `61,73` | 32 | ultraparallel `4` | 462 | `ea2bb93e9614c03d` |
| 5 | `61,77` | 192 | parallel | 471 | `4932a8dee83e1117` |
| 6 | `61,78` | 192 | parallel | 471 | `1e3583af23fcd447` |
| 7 | `61,92` | 192 | parallel | 471 | `103dcfd2b883d2ce` |
| 8 | `61,113` | 384 | ultraparallel `9` | 471 | `3c8f986ffb76646d` |
| 9 | `61,125` | 192 | ultraparallel `9` | 471 | `f23948987121c1a7` |
| 10 | `61,127` | 192 | ultraparallel `9` | 471 | `afa546b70ceea513` |
| 11 | `61,146` | 192 | ultraparallel `9` | 471 | `3e94a2b219da9adb` |
| 12 | `77,78` | 192 | ultraparallel `4` | 480 | `9c2647c448eb8595` |
| 13 | `77,79` | 48 | ultraparallel `9` | 480 | `14e57b6cc0cc19ea` |
| 14 | `77,80` | 384 | ultraparallel `16` | 480 | `8b4594e4356b953d` |
| 15 | `77,83` | 192 | ultraparallel `4` | 480 | `c687bee99ed9b0a7` |
| 16 | `77,84` | 384 | ultraparallel `25` | 480 | `cad89804016c66d7` |
| 17 | `77,85` | 384 | ultraparallel `16` | 480 | `59fe5591c7749f22` |
| 18 | `77,86` | 192 | ultraparallel `36` | 480 | `a622ee292c3fce70` |
| 19 | `77,87` | 96 | ultraparallel `9` | 480 | `383945ab8de68fd2` |
| 20 | `77,90` | 192 | ultraparallel `4` | 480 | `d742dd495523613a` |
| 21 | `77,92` | 384 | ultraparallel `36` | 480 | `fc8680cda65ecce4` |
| 22 | `77,93` | 192 | ultraparallel `16` | 480 | `3fc6eb78d60d1320` |
| 23 | `77,94` | 192 | ultraparallel `9` | 480 | `da82f67758414088` |
| 24 | `77,109` | 96 | ultraparallel `9` | 480 | `b5ebc22c28263675` |
| 25 | `77,110` | 96 | ultraparallel `49` | 480 | `15677970cb7bda9c` |
| 26 | `77,111` | 192 | ultraparallel `25` | 480 | `d00946f8e1ada9c8` |
| 27 | `77,113` | 384 | ultraparallel `36` | 480 | `dce0ff5e923dacc5` |
| 28 | `77,114` | 192 | ultraparallel `16` | 480 | `866454252990c9ca` |
| 29 | `77,115` | 96 | ultraparallel `25` | 480 | `841ec3b8ccd777cb` |
| 30 | `77,117` | 96 | ultraparallel `25` | 480 | `da97b5cb7b5d85e8` |
| 31 | `77,125` | 48 | ultraparallel `9` | 480 | `c0bb3198e3fdb85f` |
| 32 | `77,126` | 48 | ultraparallel `49` | 480 | `69aa7e61fa7a9c09` |
| 33 | `77,132` | 192 | ultraparallel `64` | 480 | `5cc377142cbf691e` |
| 34 | `77,138` | 96 | ultraparallel `49` | 480 | `16fcb65947591206` |
| 35 | `77,143` | 192 | ultraparallel `16` | 480 | `76c698a84f802ffb` |

Here the ultraparallel number is the normalized squared Gram ratio
`(r,s)^2/(r^2 s^2)`.

## Coxeter and volume certification

For every representative, the program forms the integral reflection
matrices and the chamber set

```text
{1, R_a}                    (one edge),
{1, R_a, R_a R_b}           (three-vertex path).
```

It deletes the one or two internal mirrors, reconstructs all transported
seed walls, and deduplicates exact primitive hyperplanes.  A summed strict
interior orients the surviving roots.  Direct integer replay then checks:

1. every boundary root is primitive and its reflection is integral;
2. every boundary pair is acute and has an exact Coxeter, parallel, or
   ultraparallel relation;
3. the boundary has full rank 13; and
4. the unchanged exterior core `W1,...,W60` has Smith diagonal `1^13`.

Allcock's finite-subtree theorem applies because the seed is a finite-volume
Coxeter polytope and all selected walls are pairwise disjoint doubling walls.
Therefore each reconstructed union is a finite-sided Coxeter polytope with
reflection index and volume multiplier two or three.  The seed group is a
nonuniform lattice, and a finite-index subgroup remains nonuniform, so every
descendant is finite-volume **noncompact**, not merely finite-sided.

The fixed `W1,...,W60` core also proves that every descendant still spans the
full Picard lattice with index one.  Hence the same seven-dimensional
very-general K3 family weakly carries all 37 diagram types as integral
divisor-reflection chambers.  This is not a nef-cone, effectivity, Hodge, or
surface-automorphism claim.

## Lineage boundary

The construction mechanism is classical: Daniel Allcock's finite-subtree
theorem already gives infinitely and exponentially many finite-volume
Coxeter polytopes in `H12`.  The new project-level content is the exhaustive
small-tree classification for this particular surface-derived 112-wall
stock, the 35 explicit exact diagrams, the new ledger facet counts
`462/471/480`, and their common full-Picard K3 carrier.

The pinned CoxIter lanes contain 1,214 ordinary files and 2,959 index-two
files; their `H12` facet headers are respectively `13,14,16,18` and
`18,20,25,35`.  Thus `462,471,480` are absent there.  This is only bounded
same-header evidence, not a publication-complete novelty proof.  Worldwide
priority for the exact diagrams remains open.

## Reproduction

```bash
PYTHONPATH=.:src python3 scripts/certify_jacobian_kummer_h12_small_subtree_census.py --check-archive
PYTHONPATH=.:src pytest -q tests/test_jacobian_kummer_h12_small_subtree_census.py
```

Artifacts:

- [`jacobian_kummer_h12_small_subtree_census.json`](../data/enumerations/jacobian_kummer_h12_small_subtree_census.json)
- [`certify_jacobian_kummer_h12_small_subtree_census.py`](../scripts/certify_jacobian_kummer_h12_small_subtree_census.py)
- [`test_jacobian_kummer_h12_small_subtree_census.py`](../tests/test_jacobian_kummer_h12_small_subtree_census.py)

Frozen JSON SHA-256:

`796575ebbaaf8ec0ba3054911397d83834d8939e1eb4c20cd24aae1b56aee56d`

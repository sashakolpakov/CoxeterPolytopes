# Independent full-prefix and nonreflectivity review for the silver E8 host

The literal lattice `〈−1−√2〉 ⊕ E8 ⊕ 〈1〉` is nonreflective: its full
integral reflection chamber has infinite volume. The separate
[Fraction-only verifier](../scripts/fresh_silver_e8_full_chamber_20261016_review.py)
and [certificate](fresh_silver_e8_full_chamber_20261016_review.json) pass
generation and a fresh saved replay. They import no producer arithmetic.
This excludes this specific arithmetic H9 host, without excluding other
lattices or proving any dimension bound.

## Independent enumeration and witness arithmetic

The verifier reconstructs the E8 Cartan matrix and host Gram matrix from
literal entries. It checks unimodularity, physical signature and positive
definiteness at the other embedding. It independently enumerates all three
primitive norm classes `1, 2+√2, 2`, the complete zero-height dominant census,
and all 29 positive-height controls through height 100.

The dominant-covector recursion differs from the producer. It uses only the
physical quadratic energy of the current coordinates in the entrywise
positive inverse E8 matrix. It does not use the producer's leading-Cartan
Schur bounds. Both-embedding coordinate bounds give a finite alphabet;
the conjugate energy is checked at every complete tuple. This independently
visits 132,442 positive-height prefixes and reproduces every candidate count,
selected root count and parallel deduplication. The wall set agrees after
every height packet, ending with 48 walls.

All 48 reflections are integral. The verifier checks strict interior and
actual-facet points, three complete outward E8+A1 corner sectors, and the
basepoint decompositions needed for persistence. It independently checks
both displayed integral isometries and their integral inverses, their
positive future corner scales, and their images of all nine outward corner
normals. Their product squared has determinant one and trace

`10 + 4√2 > 10`.

A finite-order real matrix of size ten has trace at most ten. Thus this
integral full-chamber stabilizer has infinite order. The replacement strict
escape for the finite prefix is also checked, although the full-chamber
symmetry is the stronger host-wide obstruction. The separate complete
100-vertex/182-edge census is not independently rerun here; only the three
corners needed by the symmetry proof are checked.

## Mathematical inputs and scope

Because the lattice is unimodular over the PID Z[√2], a primitive integral
reflection norm divides 2 up to a unit. Totally positive units are squares;
the totally positive representative of the ramified prime is 2+√2.
This gives exactly the three norm classes used above. The finite E8 Weyl
group and scalar reflection place each zero-height root in the dominant
sector, whose complete census is checked directly.

For each retained corner, the basepoint is a positive multiple of its
future timelike vertex plus a nonpositive combination of its outward
simple normals. A later simple root is nonpositive against those normals
and the basepoint. This forces a strictly negative pairing with the vertex,
so the complete corner survives in the full chamber. An integral future
isometry mapping one complete outward corner sector to another therefore
maps the full root-arrangement chamber to itself.

The integral isometry group is discrete because the other real embedding
is positive definite. An infinite discrete stabilizer would give infinitely
many disjoint translates of a small interior ball inside the chamber, so
the chamber cannot have finite volume. Consequently no finite-volume
Coxeter polytope can have all its wall reflections preserving this literal
lattice. The proof does not exclude other commensurable lattices or walls
whose reflections fail that integrality condition.

The [primary continuation and symmetry](fresh_silver_e8_h9_continuation_20261016.md)
remain separately source-bound. This is internal mathematical and
computational verification, not external peer review or a novelty claim.

```sh
python3 scripts/fresh_silver_e8_full_chamber_20261016_review.py
python3 scripts/fresh_silver_e8_full_chamber_20261016_review.py --verify
```

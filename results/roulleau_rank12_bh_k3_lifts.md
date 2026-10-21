# Three rank-12 Roulleau folds now have genuine K3-automorphism realizations

## Positive theorem

There are complex projective K3 surfaces with the following generic
Néron--Severi lattices and effective automorphisms.  After marking their
Kähler chambers as Roulleau's complete nef chambers, the indicated
automorphisms induce the exact archived wall actions.

| polytope | generic `NS(X)` | Brandhorst--Hofmann row | effective group / Picard image | Vinberg census |
|---|---|---|---|---|
| `R12-E8A1x2-C2-H10-12` | `U+E8(-1)+A1(-1)^2` | `0.4.4.1` | `C4 / C2` | 12 facets, 17 finite vertices, 2 cusps, 99 elliptic edges |
| `R12-D8A1x2-C2a-H10-15` | `U+D8(-1)+A1(-1)^2` | `0.4.4.2` | `C4 / C2` | 15 facets, 42 finite vertices, 10 cusps, 294 elliptic edges |
| `R12-D4x2A1x2-p67-H10-38` | `U+D4(-1)^2+A1(-1)^2` | `0.4.4.6` | `C4 / C2` | 38 facets, 526 finite vertices, 141 cusps, 4,518 elliptic edges |

All three polytopes are finite-volume and noncompact.  Every literal
Vinberg edge and parabolic-neighborhood audit is green.

There is also a genuinely different realization of the third polytope.
Full-database row `1.4.4.2` has effective group `C4 x C2`, symplectic kernel
`C2`, generic `NS=U+D4(-1)^2+A1(-1)^2`, and Picard image `C2 x C2`.  Its
distinguished order-four automorphism restricts to the same archived `p67`
Picard involution.  The common algebraic coinvariant of the full Picard image
has rank 8, determinant 256, and a Gram equal to twice an even
negative-definite integral Gram; consequently it contains no `(-2)`-vectors.

Thus the theorem gives **three newly geometric chamber realizations and four
effective K3-group realizations**.  It does not claim three new abstract
Coxeter polytopes: the chambers were already certified, and the `H10/12`
abstract chamber is classical.  Nor is a new arithmetic commensurability
class claimed.

## Why these are exact wall actions

For each of the four rows, the certificate reconstructs the primitive
cyclotomic transcendental lattice inside the full even unimodular K3 lattice
and takes its integral orthogonal complement.  This gives the generic Picard
lattice, not a rational approximation or a same-determinant guess.

The selected Picard action has characteristic polynomial

`(lambda - 1)^11 (lambda + 1)`.

Its fixed/coinvariant decomposition has rank-one coinvariant `<-4>` of
divisibility 2 and primitive-glue index 2.  On the fixed lattice the
certificate supplies an explicit discriminant-form isometry:

| row | fixed Smith factors | discriminant elements checked |
|---|---:|---:|
| `0.4.4.1` | `4` | 4 |
| `0.4.4.2` | `2,2,4` | 16 |
| `0.4.4.6` | `2,2,2,2,4` | 64 |
| `1.4.4.2` | `2,2,2,2,4` | 64 |

Each displayed map is bijective, preserves the quadratic form modulo `2Z`,
and carries every primitive glue pair.  The fixed lattices have signature
`(1,10)` and rank at least discriminant length plus two, so the indefinite
lattice lifting theorem promotes these finite-form maps to integral
isometries.  Together with the identity intertwiner on `<-4>`, this proves
integral conjugacy of the **full Néron--Severi action pairs**.

The database rows are effective characters, hence come with an invariant
Kähler class.  Transport that Kähler chamber to the archived complete root
chamber by the `(-2)`-Weyl group.  The resulting automorphism lies in the
exhaustively computed weighted chamber group.  For `U+E8+A1^2` its
nontrivial element is unique; for the other two lattices the fixed
discriminant form separates the two possible high-rank involution classes.
Thus a final chamber symmetry conjugates to the literal archived generator.
This closes the Kähler/Torelli and chamber-conjugacy gates, rather than merely
matching fixed-lattice invariants.

## Finite positive filter

The same exact scanner was run over both Brandhorst--Hofmann lanes:

- all 587 purely non-symplectic rows: 53 had generic Picard rank 12 or 13,
  11 had one of the five target Picard lattices, and exactly the three rows
  above matched an archived target action;
- all 3,561 full effective-group rows: 92 had generic Picard rank 12 or 13,
  four had a target Picard lattice, and exactly `1.4.4.2` matched a target
  action;
- no rank-13 target action survived.

The empty rank-13 output is retained only as a reusable finite search
boundary.  The positive rows are the theorem.

## Reproducibility

The offline certificate is
[`scripts/certify_roulleau_rank12_bh_k3_lifts.py`](../scripts/certify_roulleau_rank12_bh_k3_lifts.py),
with machine artifact
[`data/enumerations/roulleau_rank12_bh_k3_lifts.json`](../data/enumerations/roulleau_rank12_bh_k3_lifts.json).

Run the pinned local proof with:

```bash
PYTHONPATH=.:src python3 scripts/certify_roulleau_rank12_bh_k3_lifts.py --check-archive
```

Replay all four complete 22-dimensional database rows with:

```bash
PYTHONPATH=.:src python3 scripts/certify_roulleau_rank12_bh_k3_lifts.py \
  --brandhorst-hofmann-zip /path/to/K3Groups.zip --stdout
```

The external archive is Brandhorst--Hofmann database version 2.1,
DOI `10.5281/zenodo.10410278`, pinned SHA-256
`c1272db2c0fe2ae0a788aed7745201fccd2a1337b6c2d85555880c434cab4d7c`.

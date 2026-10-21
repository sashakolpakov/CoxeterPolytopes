# Two positive `H14` chambers from the Type-VI Coble divisor lattice

**Status (2026-08-09): exact positive theorem; worldwide priority open.**

The complete `H14/60` nef chamber of Kondō's characteristic-three,
five-boundary Type-VI Coble surface has exactly `35` doubling walls.  Its
weighted Gram automorphism group has order `120` and splits them into two
orbits, of sizes `15` and `20`.  Doubling once gives exactly two
nonisomorphic finite-volume noncompact Coxeter polytopes:

| chamber | marked choices | facets | finite vertices | ideal vertices | hyperideal rays |
|---|---:|---:|---:|---:|---:|
| `COBLE-VI3-D2A-H14-100` | 15 | 100 | 54,318 | 1,993 | **0** |
| `COBLE-VI3-D2B-H14-101` | 20 | 101 | 55,130 | 2,043 | **0** |

These are genuine additional abstract Coxeter polytopes and exact integral
divisor-lattice realizations, not merely proposed wall systems.  They are
new to this project's positive ledger.  They are not claimed worldwide-new:
Allcock already proved infinitely many finite-volume examples in `H14`, and
the two pinned CoxIter lanes are not a publication census.

Both catches are also **redoublable**.  Consequently each seeds an infinite
family of pairwise nonisometric finite-volume `H14` Coxeter polytopes on the
same divisor lattice.  This family statement uses Allcock's doubling lemma;
no closed facet-count formula is claimed below.

## The algebraic starting object

The input is not a random Coxeter matrix.  It is the full Néron--Severi nef
cone of the Type-VI Coble surface, with complete physical negative-curve
stock

```text
5*(-4) + 10*(-2) + 45*(-1).
```

After adjoining the five boundary half-classes and doubling the form, the
arithmetic divisor lattice is

```text
2 Pic~(S) = U(2) + E8(-1) + A1(-1)^5,
```

of signature `(1,14)`.  The sixty primitive mirror normals span it with
index one and have split

```text
15*(-2,div1) + 35*(-2,div2) + 10*(-4,div2).
```

The `35` middle roots are precisely the one-boundary `(-1)`-curves recovered
from primitive Coble--Mukai cusps.  Each has no odd finite-angle neighbor,
so each is a doubling wall.  The other 25 walls have an `m=3` neighbor and
fail the gate.

For contrast, the other sharp Coble endpoint, Type MII `H17/48`, has no
doubling wall: its 48 facets have respectively `2`, `4`, or `8` incident
`m=3` walls, with histogram

```text
2^16, 4^8, 8^24.
```

Thus the Type-MII branch closes exactly for this factory, while Type VI is
productive.

## Exact doubling and irredundancy

Let `P` be the Type-VI chamber and let `e` be one of the 35 selected
square-minus-one physical divisor walls.  Every finite dihedral angle at
`e` is an even submultiple of `pi`.  Therefore

```text
D_e(P) = P union s_e(P)
```

is convex and is again a Coxeter polytope.  Its boundary consists of the
walls of the two copies other than `e`, modulo the walls fixed by `s_e`.

The exact seed weighted-Gram automorphism group has order `120`.  It has
two orbits on the 35 eligible walls:

- 15 walls are orthogonal to 18 other seed facets and give
  `2*(60-1)-18 = 100` facets;
- 20 walls are orthogonal to 17 other seed facets and give
  `2*(60-1)-17 = 101` facets.

An automorphism carrying one gluing wall to another carries the entire
doubled boundary with it, so these are at most two types.  Their different
facet counts prove that they are not isomorphic, hence there are exactly
two.

The certificate archives every primitive doubled wall vector, the full raw
wall Gram, and the full exact normalized-squared Gram.  It separately
repeats the construction in the physical lattice `NS(S)=I_(1,14)`: the
glued class has square `-1`, its reflection is an integral `NS` isometry,
the resulting 100/101 primitive divisor vectors span `NS` with index one,
and their normalized Gram agrees entry by entry with the arithmetic Gram.
The results are:

| chamber | primitive mirror types `(-2,1);(-2,2);(-4,2)` | `m2;m3;m4;parallel;ultraparallel` | wall-Gram inertia |
|---|---:|---:|---:|
| `H14/100` | `18;68;14` | `1566;53;128;611;2592` | `(1,14,85)` |
| `H14/101` | `20;68;13` | `1578;55;137;609;2671` | `(1,14,86)` |

All off-diagonal wall pairings are nonnegative.  Every normalized ratio
below one is exactly `0`, `1/4`, or `1/2`, so every intersecting pair has
angle `pi/2`, `pi/3`, or `pi/4`; ratio one is parallel and larger ratios are
ultraparallel.  Every wall reflection is a primitive integral involutive
isometry, and each doubled wall set still spans the arithmetic lattice with
index one.

Irredundancy is not inferred from the wall count.  A strict seed interior
`h` has `h^2=70` and minimum wall pairing one.  For every seed facet `w`,

```text
h_w = h - (h.w/w^2) w
```

is timelike, lies on `w`, and is strict at every other seed wall.  Transport
`h_w` into the relevant copy.  The replay checks it against every doubled
inequality and obtains zero at its target and strict positivity everywhere
else.  Thus all 100 and all 101 inequalities support genuine facets.

## Literal finite-volume gate and the flat-vertex correction

Finite volume is true geometrically because each double is a union of two
copies of the known finite-volume chamber and its reflection group has index
two.  The certificate nevertheless performs the requested zero-hyperideal
test from the complete exact seed ray archive.

The seed has

```text
27,845 finite + 1,027 ideal + 0 hyperideal rays.
```

A seed ray off the glued facet occurs once in each copy.  A ray on the glued
facet needs care: after deleting the glued normal, it remains extreme exactly
when another active facet is nonorthogonal to that normal.  The reflected
copy then restores the missing normal direction.  If every other active wall
is orthogonal, reflection fixes the whole active set and the putative vertex
flattens.  This is the correction that a naive `2R-R_on_wall` count misses.

The exact accounting is:

| chamber | seed rays on glue `(finite,ideal)` | retained on glue | flattened on glue | doubled result `(finite,ideal,hyperideal)` |
|---|---:|---:|---:|---:|
| `H14/100` | `(1052,61)` | `(732,61)` | `(320,0)` | `(54318,1993,0)` |
| `H14/101` | `(280,11)` | `(0,11)` | `(280,0)` | `(55130,2043,0)` |

Thus both chambers pass the direct cone criterion with exactly zero
hyperideal extreme rays.  Their ideal rays prove noncompactness.

## What is algebraic here, and what is not

Every new wall is an integral reflection image of a divisor class in the
same Type-VI Picard/arithmetic lattice.  This is therefore a genuine
algebro-geometric divisor-lattice factory, and not a search through arbitrary
Coxeter matrices.

But doubling deliberately crosses one wall of the nef cone.  The reflected
divisor classes are not asserted effective on the original Coble surface,
and the doubled chambers are not asserted to be nef cones or fundamental
domains of surface automorphism groups.  The exact gain is:

```text
complete surface nef chamber
    -> exact integral divisor reflections
    -> two complete finite-index Coxeter chambers.
```

It is a positive `SURFACE -> POLYTOPE` theorem, not the still-missing broad
`POLYTOPE -> SURFACE -> classification -> theorem on P` leverage theorem.

## Exact redoubling corollary: two infinite `H14` families

The catch is not isolated.  Each new chamber has exactly `68` doubling
walls.  Among them the `H14/100` chamber has `316` parallel pairs and the
`H14/101` chamber has `324`.

More concretely, both chambers contain the same archived pair

```text
a = (0,0,0,0,0,2,0,2,-1,3,1,0,0,0,-2),
b = (0,0,1,0,2,1,0,0, 1,1,1,0,0,0,-2).
```

In the arithmetic divisor lattice,

```text
a^2 = b^2 = -2,   a.b = 2,   (a+b)^2 = 0.
```

Their normalized squared ratio is exactly one, so the two walls are
parallel/asymptotic and do not meet in `H14`.  The full relation replay also
shows that both are doubling walls.  Hence each new chamber is redoublable
in Allcock's precise sense.

Allcock's Lemma 1.3 says that doubling a redoublable Coxeter polytope again
produces a redoublable Coxeter polytope.  His finite-subtree theorem and
Corollary 4.2 identify every finite path in the two-generator Cayley line
with a reflection-subgroup fundamental chamber.  Therefore, for every
`N >= 1`, a path of `N` copies of either new chamber is a finite-volume
noncompact Coxeter polytope with

```text
index in the H14/100 or H14/101 group = N,
index in the parent H14/60 group       = 2N,
volume / volume(parent H14/60)         = 2N.
```

Different `N` have different volume and are therefore nonisometric.  Every
boundary wall remains an integral reflection image of a divisor class in
the same lattice.  The exact existence of two infinite algebro-geometric
divisor-lattice families is thus proved.  What remains deliberately open is
a closed facet formula and a collision classification for all path lengths,
and whether the two families overlap across different seeds.

The relevant general doubling mechanism is Allcock's
[*Infinitely many hyperbolic Coxeter groups through dimension 19*](https://arxiv.org/abs/0903.0138),
Theorem 4.1 and Corollary 4.2.

## Lineage boundary

- `H14` was already populated, indeed infinitely, before this project.
- Neither `H14/100` nor `H14/101` occurs as a header in the pinned 1,214
  ordinary CoxIter records or the pinned 2,959 index-two-doubling records.
- That absence is bounded-corpus evidence only.  The same corpus omits known
  controls, including the parent `H14/60` chamber.
- The allowed wording is therefore **new-to-project abstract `H14/100` and
  `H14/101` divisor-lattice chambers; worldwide priority open**.

## Reproducibility

- [certifier](../scripts/certify_coble_type_vi_h14_single_doublings.py)
- [certificate](../data/enumerations/coble_type_vi_h14_single_doublings.json)
- [test](../tests/test_coble_type_vi_h14_single_doublings.py)
- [parent Type-VI theorem](coble_type_vi_char3_five_boundary_full_ns_h14.md)

Replay with

```bash
PYTHONPATH=src:. python3 \
  scripts/certify_coble_type_vi_h14_single_doublings.py --check-archive

# Optional bounded-lineage replay, when the pinned local CoxIter corpus exists:
PYTHONPATH=src:. python3 \
  scripts/certify_coble_type_vi_h14_single_doublings.py \
  --check-archive --check-corpus

PYTHONPATH=src:. pytest -q \
  tests/test_coble_type_vi_h14_single_doublings.py
```

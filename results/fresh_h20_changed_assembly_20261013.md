# Doubling the H19 source before a changed H20 lift

This round constructs an exact **179-facet Coxeter chamber in H20**, replacing the previous long-wall core and removing its deficient last facet. It cuts the previous 176-wall escape. A new exact strict spacelike escape nevertheless proves **infinite volume**, and seventeen inherited cusp links are deficient. This is constructive progress outside the earlier fixed-core obstructions, not a new finite-volume example or a claim of global novelty.

## Changed finite-volume source and lift

Start with the archived [Vinberg X3 H19 chamber](../data/enumerations/vinberg_x3_h19.json), using the coordinates and form reconstructed in the [original dimension-bridge certificate](../data/enumerations/fresh_dimension_bridge_20261013.json). The 36 source walls have squared norms 2 and 6. Write `a=r24`, with `a²=6`.

The exact check finds 21 right-angle relations, three π/6 relations, four parallel pairs and seven dotted pairs at `a`. Thus `a` is an even doubling wall. The union with its reflection across `a` is a finite-volume Coxeter chamber by the doubling lemma: [Allcock, *Infinitely many hyperbolic Coxeter groups through dimension 19*, Lemma 1.3](https://web.ma.utexas.edu/users/allcock/research/doubling.pdf). Reflection is explicitly

```
s_a(r) = r − (r,a) a/3.
```

Remove the shared wall and merge the identical boundary planes. This gives **49 H19 walls: 27 short and 22 long**. The exact full pair census is 771 right, 36 of angle π/3, 66 of angle π/6, 259 dotted and 44 parallel.

Adjoin an orthogonal coordinate `z²=3`. The graph on the 22 long normals with edges `(r_i,r_j)=−6` is connected and bipartite. Fix its first node's sign to +1 and write its signs `ε_i`. Keep each short wall at height zero; replace each long wall by

```
(r_i, 2 ε_i),   square 18,
(r_i, −ε_i),   square 9.
```

The resulting **71-wall core** has 1,393 right pairs, 102 π/3 pairs, 106 π/4 pairs, 704 dotted pairs and 180 parallel pairs. These numbers sum to all 2,485 unordered pairs. The form still has signature (20,1).

This phase is opposite to the original phase on the surviving original long directions. All 24 original long-wall lifts are removed. Retain exactly 102 compatible later walls from the frozen [176-wall certificate](../data/enumerations/fresh_h20_coupled_completion_20261013_second_exact.json), obtaining **173 walls**. In particular, the previous last facet 175 is absent. This changes the core rather than adding dotted walls to that excluded fixed facet.

The strict center is the old center projected to `a⊥`:

```
h' = h − (h,a) a/6,     (h')² = −7599/2.
```

All 7,242 core–retained pairs pass exactly. The previous 176-wall strict escape is cut by new core wall 59 with pairing **12171/4096 > 0**.

## Six coherent polar cuts and the remaining escape

A bounded numerical search proposed rational polars. Six were accepted in exact arithmetic; each passes every pair gate against all preceding walls and cuts an explicit strict spacelike point of the preceding chamber.

| New wall | Square | Right | π/3 | π/4 | Parallel | Dotted |
|---|---:|---:|---:|---:|---:|---:|
| 173 | 1 | 22 | 0 | 3 | 3 | 145 |
| 174 | 18 | 21 | 2 | 0 | 7 | 144 |
| 175 | 4 | 20 | 0 | 2 | 0 | 153 |
| 176 | 8 | 22 | 2 | 0 | 4 | 148 |
| 177 | 72 | 21 | 2 | 0 | 6 | 148 |
| 178 | 144 | 20 | 0 | 1 | 1 | 156 |

These labels are the new list's indices; its wall 175 is unrelated to the removed old wall 175. The first new polar, in the existing 21-coordinate basis, is

```
[52,-22,-14,7,8,13,29,36,39,-13,-20,-11,-7,-8,-13,18,22,26,-3,-6,-1].
```

The [exact archive](../data/enumerations/fresh_h20_changed_assembly_20261013_exact.json) supplies all six vectors, every new pairing, strict target cuts, center and final escape. It verifies 10,780 changed/new unordered pairs; the remaining 5,151 pairs are inherited from the retained old subset. Forty recursively hash-bound inputs tie that old-old dependency to the files actually reconstructed.

Every retained inequality is an actual facet. Project the strict timelike center orthogonally to its spacelike normal. The projected point remains timelike, and the complete nonpositive pair gate makes every other wall inequality strict. The same projection to a positive-definite two-wall span certifies the finite-angle ridges. Consequently this is a Coxeter chamber, not merely a collection of formally compatible Gram entries.

The final exact witness satisfies all 179 inequalities strictly and has

```
x² = 10640902782393308787 / 20000000000000000000 > 0,
(x,h') = −308390396781 / 5000000000 < 0.
```

Together with the strict timelike center, this gives an open set of ideal directions and proves infinite volume. The finite polar search does not exclude further compatible additions or another replacement.

## Seventeen deficient inherited cusp links

The separate [cusp diagnostic](../scripts/fresh_h20_changed_assembly_20261013_cusps.py) rebuilds the union of the 107 original source cusp rays and their reflected rays. There are **204 distinct primitive rays**, and all are feasible for the final 179 walls. All 36,516 signs are checked exactly in the [cusp archive](../data/enumerations/fresh_h20_changed_assembly_20261013_cusps.json).

Seventeen rays have only short, height-zero incident normals: eleven reflected copies of the original short-only cusps and six cusps lying on the removed doubling wall. At each of these seventeen rays, the actual normal span has rank 19, and the Euclidean quotient has affine rank **18 instead of 19**. The checker decomposes its negative-pair graph and verifies a strictly positive kernel for every component; each weighted normal sum equals a positive multiple of the null ray. Thus the rank statement describes actual positive affine components, not just a numerical rank deficiency.

The explicit vector `(0,...,0,1)` has square 3, is perpendicular to the null ray and to every incident normal, and gives a free horospherical direction. The links are therefore unbounded. More directly, a sufficiently small nonzero displacement of the null ray in that vertical direction is spacelike and leaves all nonincident inequalities strict; a still smaller displacement toward `h'` makes the incident inequalities strict as well. This is an independent local explanation of infinite volume.

This census covers the 204 inherited rays only. It does not enumerate every cusp of the new chamber, nor does this diagnostic claim that every other inherited link is bounded. Completing these seventeen links alone would not invalidate the separately saved remote escape without checking its new inequalities.

## Scope and replay

The changed source and all 179 accepted normals are exact. Numerical optimization only proposes indices and rational polar vectors; it does not certify angle compatibility or volume. Earlier bounded reflection, copied-neighbor and freed-anchor scouts are exploratory inputs, not universal exclusions. Their failure did not enter the final construction proof.

The [main checker](../scripts/fresh_h20_changed_assembly_20261013_exact.py) and cusp checker use the existing exact rational/radical implementations. Fresh generation and deterministic replay passed:

```sh
python scripts/fresh_h20_changed_assembly_20261013_exact.py --verify
python scripts/fresh_h20_changed_assembly_20261013_cusps.py --verify
```

The [independent root review](fresh_h20_changed_assembly_20261013_root_review.md) separately reconstructs the rational doubled source, 71-wall lift, six rational polars and their 77-wall subset using Fraction arithmetic, and independently checks the old target cut. It also reruns the full primary 179-wall checker; its scope distinguishes that shared radical replay from the separate rational arithmetic.

A separate [cusp source review and fresh replay](fresh_h20_changed_assembly_20261013_cusp_review.md) passed as well. It checks the distinction between the normal span and horospherical quotient, the positive affine relations, and the explicit free vertical direction; it shares the frozen chamber's exact arithmetic and signature/center dependencies.

The exact construction depends on the frozen prior 176-wall pair certificate for old-old relations. It does not infer finite H20 volume from the doubled source's finite H19 volume. Further work must address the actual missing cusp direction and remaining escape while preserving the full angle gate.

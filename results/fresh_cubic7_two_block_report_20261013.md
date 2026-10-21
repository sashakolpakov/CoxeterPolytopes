# Two-block cubic-seven completion at the unchanged height bound

The complete two-block search through height 1000 produced two exact Coxeter completions, with 158 and 122 walls after a separately verified scalar Pell cap. Both retain a spacelike escape. Reconsidering the previous wall selection found a useful cut that would have been lost by retaining all 73 earlier walls, but it did not close the volume.

Write τ=2cos(2π/7), C=τ²+τ−1, β=2−τ and D=2−τ². The host is <-2τ>+B⁴+<1>, where B has diagonal 2 and off-diagonal -C. The starting 25-wall height-30 cone and the earlier 73-wall asymmetric completion are linked directly from the artifacts below.

## Complete bounded support search

The new mechanism allows two fully variable dihedral blocks, with final scalar coordinate zero. In each block the full zero-height chamber forces covector coordinates 0≤u≤w. Equality and the zero block are included. The integral block vector and its energy are

```
y = ((2u+Cw)/β, (Cu+2w)/β),
E = u*y1+w*y2.
```

All scalar-zero lattice vectors have even norm, so the complete primitive norm classes in this slice are 2 and 2β. These classes are squarefree as ideals and hence cannot give imprimitive integral normals.

For height H=1000, the time controls obey x²≤Hq/(2τ) in the chosen embedding and x'²≤-q'/(2τ') at both other embeddings. Each block covector satisfies u²,w²≤2q(H+1) in the chosen embedding and ≤2q at the other embeddings. These follow from Cauchy bounds for the positive blocks. The exact thin-box enumerator covers all these alphabets.

There are 80 admissible block vectors for each norm class. Exact energy joins solve E1+E2=q+2τx², and all 12 ordered embeddings into the four spatial blocks are tested. Root integrality, reflection integrality, all retained-wall pairings, and every finite Coxeter angle are checked. The finite-angle gate rejects the arithmetic acute angle 2π/7, whose value 4cos² is τ². Ordered embedding counters can include duplicates; the saved dictionaries identify rays exactly.

This is a complete search for the stated support and height. It is not a uniform ambient Vinberg search and omits walls with a nonzero final scalar, apart from the separately appended Pell cap.

## Retaining the previous 73 walls

The original 73-wall cone admits 204 new compatible two-block normals: 114 of norm 2 and 90 of norm 2β. None cuts its saved escape. A height-ordered cross-compatible selection adds 84 walls. The norm-one Pell cap from the earlier certificate is then appended after checking all its pairings and angles, giving **158 walls**.

The complete dictionary, selected indices, final walls, strict timelike interior and exact escape are in [the retained-base artifact](fresh_cubic7_two_block_h1000_20261013.json). The cap lies beyond the bounded two-block census and is explicitly marked `additional_pell_cap`.

## Replacing premature caps

Retaining all 73 walls would exclude six further norm-two candidates compatible with the original 25-wall cone. One of those six cuts the saved 73-wall escape. It has

```
x = 1+2τ+τ²,
y1=y2=y5=y6=1+4τ+2τ²,
all other spatial coordinates zero,
q = 2,
height = 5+14τ+6τ².
```

Its pairing with that escape is (10-20τ+12τ²)/29>0. The full original-25-compatible dictionary has 258 roots, including the earlier 48 asymmetric caps. Rebuilding in height order selects this useful cut, retains 36 of those 48 caps, and selects 96 dictionary roots overall. With the original 25 walls and the Pell cap this gives **122 walls**.

The [replacement artifact](fresh_cubic7_two_block_replacement_h1000_20261013.json) records the entire original-25 dictionary, membership in the previous 73-wall cone, the rejected-new-wall flags, selected indices and final geometry. The selected set is a height-ordered compatible subset; no maximum-cardinality claim is made.

Both final cones admit the same certified spacelike vector:

```
x=1,
y5=y6=1+τ,
z=(32+36τ-20τ²)/29,
all other spatial coordinates zero,
norm=(866-2658τ+2656τ²)/841 > 0.
```

The norm-one cap is active at this vector, so the certificate already accounts for that cap. This proves infinite volume of these two finite cones, not nonreflectivity of the full host or impossibility of other two-block selections outside the stated height bound.

## Files and verification

The producer is [fresh_cubic7_two_block_search_20261013.py](../scripts/fresh_cubic7_two_block_search_20261013.py). It writes both artifacts and checks every final cone angle, strict interior and exact geometry certificate. The independent follow-up verifier is maintained separately by the coordinating agent; this report does not claim a second full candidate enumeration.

```sh
PYTHONPATH=src:. python3 scripts/fresh_cubic7_two_block_search_20261013.py --height 1000 --output results/fresh_cubic7_two_block_h1000_20261013.json
```

The surviving vector is supported on one equal-coordinate block and the scalar line. Writing that block as (a/D,a/D) gives norm equation 2a²+z²−2τx²=q. In this saturated slice the complete primitive classes reduce to 1 and 2: reflection integrality at β would force x,a,z all divisible by β, contradicting the norm valuation. That is a different mechanism, pursued separately from this two-block run.

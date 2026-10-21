# Published precursor of the all-ideal seven-simplex example

The ten-root precursor of the exact twelve-facet H9 chamber is already
published. It is literally the bottom-right diagram in Figure 14, printed
page 24, of Hao Chen and Jean-Philippe Labbé, *Lorentzian Coxeter systems
and Boyd–Maxwell ball packings*, Geometriae Dedicata 174 (2015), 43–73
([primary PDF, version 6](https://arxiv.org/pdf/1310.8608v6),
[journal DOI](https://doi.org/10.1007/s10711-014-0004-1)).
The repository chamber is its canonical truncation by the two polar
hyperplanes. This identification supplies provenance, not a novelty claim.

## Literal comparison

The figure was checked visually in a rendered copy of the primary PDF;
the [saved source copy](../data/references/fresh_all_ideal_seven_simplex_literature_20261017_chen_labbe.pdf)
is hash-bound by the [exact checker](../scripts/fresh_all_ideal_seven_simplex_literature_20261017.py).
Number the figure's octagon clockwise from its top vertex, and then label
the interior left leaf 8 and the exterior right leaf 9. Its edges are
the eight-cycle and (6,8), (2,9), all with Coxeter label 3. The map from
these figure labels to the labels in the
[exact realization](fresh_all_ideal_seven_simplex_realization_20261017.md) is

`[1,0,7,6,5,4,3,2,9,8]`.

The checker compares every entry of the corresponding norm-two Gram
matrices. It also verifies that all 45 two-node deletions are positive
semidefinite and that a one-leaf deletion remains indefinite: this is
indeed a level-two diagram. The transcription of the printed figure is
a human-checked input; the literal graph comparison is exact.

Let G be the repository's original ten-root Gram matrix and J its inverse.
The two added roots are exactly `-2 J e8` and `-2 J e9`. These are the polar
normals of the two hyperideal simplex vertices, with the orientation used
in the frozen realization. Thus the twelve-wall completion follows from
the published ten-root diagram without an additional free parameter.
The complete finite-volume proof, all 23 vertices and proper cusp links,
and the all-ideal seven-simplex face remain certified by the realization
package; they are not inferred merely from the picture.

## Limits of the identification

This targeted check did not locate an explicit source naming this
twelve-wall completion or highlighting its all-ideal seven-simplex face.
It therefore makes neither a novelty claim for those observations nor
an exhaustive classification claim. Adding the polar reflections changes
the generating reflection group; the completed group is not being
identified with the original ten-generator group studied in the paper.

There are useful scope checks against nearby classification results.
The chamber is noncompact, so a compact n+3-facet classification does not
apply. It has three nonsimple ideal vertices and is not a pyramid, as
checked directly from the frozen vertex and facet lists. Consequently it
is outside the exactly-one-nonsimple-vertex hypothesis of
[Roberts' classification](https://arxiv.org/pdf/1511.08451).
[Tumarkin's short primary note](https://arxiv.org/pdf/math/0311272)
establishes the n+3 dimension bound and special classifications; no
literal match to its full n+3 list is asserted here.

One must also distinguish different meanings of “truncated simplex.”
The original edge between the two hyperideal vertices has square
`(t-u)^2/2` in nonnegative homogeneous coordinates. Only its midpoint is
ideal, and no part of that edge is inside hyperbolic space. Thus an
edge-finite truncated-simplex classification requiring every original
edge to meet hyperbolic space cannot be applied to this precursor.

## Verification and status

The [saved certificate](fresh_all_ideal_seven_simplex_literature_20261017.json)
binds the frozen realization, the primary PDF and the checker. Generation
and a fresh saved replay pass, as do syntax, source-hash and local-link
checks. This is a literal primary-source identification plus exact
matrix/combinatorial verification, not an independent repetition of the
entire twelve-wall cone computation. **FINAL/FROZEN.**

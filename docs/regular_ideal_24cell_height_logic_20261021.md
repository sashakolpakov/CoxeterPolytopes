# Regular ideal24-cell: opposite heights and the coordinate center

**Bounded compilation and separate actual axiom audit PASS.**
[RegularIdeal24CellHeightKernel.lean](../formal/CoxeterProofs/RegularIdeal24CellHeightKernel.lean),
SHA256 `633ae6fef8c5d0f0d0cafb2363d99f38e6922405d61c25c0e3a569e7bceaecec`
(108 lines,4720 bytes), has five actual audited declarations, matching its
five `#print axioms` commands. The
[completion receipt](../results/lean_publication_completion_20261021/manifest.json)
binds the unchanged source, object/import/compiler ledgers and actual stdout.
The compile took1.13s and the separate audit0.76s. All five outputs use only
`propext`, `Classical.choice` and `Quot.sound`, with no `sorryAx`. Two
unused-section-variable warnings are retained. No successful dependency was
rerun. The six papers and their proof maps are unchanged; the source's
preparation header is historical and retained to preserve its compiled hash.

The only import is the existing projection/height kernel
`2fac8028e432a2104a58e423c65fc38c9e214fed5f52210e408aabe3f15f0e90`.
Its finite-coordinate Lorentz definition and the cubic kernel's proved
nonnegativity of a coordinate norm are reused. No mathlib, geometry
library, saved certificate, scientific worker or extra framework is added.
The scalar type is an arbitrary linearly ordered commutative ring R;
the intended geometric instance is the real field. Four natural-number
indices address coordinates; no real height or norm is replaced by an
integer variable.

The target is Lemma `lem:center` and the equality endpoint of the ideal
proof in the current [five-page sharpRA paper](../manuscripts/right_angled_sharp_dimension_bounds_publication.tex),
source `0c29a7869fa46e4e1a542d8ae6f4dfc7e616136c7a1aad3984d941a95cf89722`.
The [current projection map](right_angled_projection_height_logic_20261020.md)
already lists this opposite-normal calculation as outside its coverage.

Write q=(a,w), with only coordinates0,...,3 used by `lorentz 4`.
The new source defines the literal two normals

    m+=(1,(1,1,0,0)), m-=(1,(-1,-1,0,0)).

These are an opposite spatial pair among the paper's normals
`(1,epsilon*e_i+delta*e_j)`. The primitive point hypotheses are ONLY
`a>0` and `lorentz 4 q q=-1`. The cap hypotheses are the actual
coordinate heights `height 4 q m+<=1` and `height 4 q m-<=1`, or their
strict variants. Unit-normal status and the height substitution are
derived from the explicit definitions, rather than supplied conclusions.

| Declaration | Derived content |
|---|---|
| `opposite_pair_identities` | Both literal normals have Lorentz square1; their heights are a-w0-w1 and a+w0+w1, whose sum is2a |
| `future_unit_time_lower_bound` | Coordinate norm nonnegativity and the primitive future-unit equation imply a>=1 |
| `opposite_height_caps_force_center` | The two weak caps imply a=1 and w_i=0 for EACH i<4 |
| `center_has_unit_heights` | At that coordinate center every normal with time coordinate1 has height1 |
| `no_strict_opposite_height_caps` | Strict caps on this pair contradict the derived time lower bound |

For the center theorem, unit timelikeness at a=1 gives a sum of four
nonnegative squares equal to0. Each square vanishes; strict positivity
of a nonzero square then forces its coordinate to vanish. The theorem
does not claim anything about the unused tail of a Nat-to-R function.
It assumes no lower bound a>=1, no center predicate, no desired norm
conclusion and no height-sum or zero-coordinate conclusion. The strict
contradiction needs only this one pair, not a spanning/uniqueness premise.

The following remain external geometric or published application inputs:

* Kolpakov's metric equality theorem, and the hyperbolic isometry placing
  an actual equality facet in the displayed24-cell normal model;
* identification of the two literal vectors with actual outward UNIT
  facet normals and of a future unit coordinate vector with the selected
  intrinsic finite point;
* transport of actual heights and the existing strict height bound through
  that isometry, and the time-coordinate1 form of the remaining normals;
* actual facet projection/inheritance, Musin's bound and the minimum facet
  count that force metric equality in the first place.

The coordinate results supply the internal argument behind the old
`WrittenRegularCellGeometry` fields. They do NOT automatically instantiate
that abstract record, construct the metric isometry or all24 actual
facets, classify kissing codes, or prove the ideal dimension bound.
No dimension wrapper, new theorem or external input is added to the paper.
The recorded compile and separate audit used the usual shared lock, one
thread,2048MiB setting and unchanged40/6GiB admission and8/3GiB abort gates.
This is actual verification of the five conditional coordinate declarations;
it does not discharge the geometric or published application inputs above.

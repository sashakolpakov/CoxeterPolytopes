# Regular ideal24-cell height kernel: outside-author source review

21 October2026. **FULL mathematical SOURCE/INTERFACE PASS before
compilation.** No mathematical worker, Lean module or compiler was
executed/imported, and no reviewed source or existing paper was changed.
This is not compilation, an actual axiom audit or an end-to-end geometric
formalization.

| Artifact | SHA256 | Bytes |
|---|---|---:|
| formal/CoxeterProofs/RegularIdeal24CellHeightKernel.lean |633ae6fef8c5d0f0d0cafb2363d99f38e6922405d61c25c0e3a569e7bceaecec|4720|
| docs/regular_ideal_24cell_height_logic_20261021.md |f57a82535c5a2cc79af40bd1a2df0e7de519597629737877bbaaaf5b890e2450|4386|
| manuscripts/right_angled_sharp_dimension_bounds_publication.tex |0c29a7869fa46e4e1a542d8ae6f4dfc7e616136c7a1aad3984d941a95cf89722|15738|

I authored neither new source nor map. I share earlier source audits and
mathematical foundations of the sharpRA packet, but independently derived
this coordinate argument before reading the new kernel. I read all108
Lean lines, all73 map lines, the current publication's normal model and
center lemma (lines295--318), and its ideal-branch application
(lines320--335). Relevant imported coordinate definitions and positivity
proofs were read, rather than assumed from earlier PASS labels.

## Coordinate mathematics and every public theorem

The source uses an arbitrary linearly ordered commutative ring R. Heights
and coordinates live in R; Nat indexes only the finite coordinate sum.
There is no integer replacement for real height, division, reciprocal or
square-root premise. The established Lorentz form is
`-q.time*m.time + sum_(i<4) q.spatial_i*m.spatial_i`.

1. **`opposite_pair_identities`, lines25--37.** The literal normals are
   `(1,(1,1,0,0))` and `(1,(-1,-1,0,0))`, an actual opposite pair in the
   paper's displayed24-normal model after shifting coordinate labels.
   Their squares are `-1+1+1=1`; no missing normalization factor occurs.
   Expanding the finite sum gives heights `a-w0-w1` and `a+w0+w1`,
   whose sum is `2a`. These identities and unit norms are derived,
   not input fields. The explicit normal tails vanish by definition.
2. **`future_unit_time_lower_bound`, lines39--48.** Primitive inputs are
   only `a>0` and Lorentz square-1. The imported coordinate positivity
   theorem gives `sum_(i<4) w_i*w_i>=0`, hence `a*a>=1`. If `a<1`,
   multiplication by positive a gives `a*a<a<1`, a contradiction.
   Thus `a>=1` is proved; no equivalent time-bound premise is supplied.
3. **`opposite_height_caps_force_center`, lines50--78.** The two weak
   caps and the proved sum give `2a<=2`; combined with the derived
   lower bound this forces `a=1`. Substitution in the unit equation
   gives the sum of FOUR nonnegative squares equal0. Each square
   vanishes, and a nonzero coordinate would have a strictly positive
   square. The finite Nat case split covers exactly0,1,2,3. Its
   conclusion explicitly restricts spatial zero to `i<4`; it does not
   wrongly identify the unused tail of a Nat-to-R function with zero.
4. **`center_has_unit_heights`, lines80--91.** Given the coordinate
   center obtained above, direct substitution makes the height of
   EVERY normal with time coordinate1 equal1. Unit/spatial conditions
   on that normal are unnecessary for this stronger algebraic statement.
   The theorem is a center implication, not a supplied geometric-center
   conclusion replacing the preceding derivation.
5. **`no_strict_opposite_height_caps`, lines93--100.** Strict caps on
   just this pair imply `2a<2`, contradicting the derived `a>=1`.
   It needs neither the center theorem nor a spanning/uniqueness premise.

All five statements match the written lemma, which needs only one
opposite pair. The future-unit and cap hypotheses are elementary point
and wall-coordinate data, rather than the desired center or contradiction.
The statements remain valid without assuming the point lies inside the
cell, so omission of a polytope-membership premise is a legitimate
strengthening. The proof uses only positive multiplication, finite
sum expansion and ordered-ring square positivity. No sign or equality
case is lost. No mathematical or source-level API blocker was found;
actual tactic elaboration must still be established by compilation.

## Imported logic, requested prints and geometric boundary

The sole direct import is `RightAngledProjectionHeightKernel`, source
2fac8028e432a2104a58e423c65fc38c9e214fed5f52210e408aabe3f15f0e90.
Its Vector/lorentz/height definitions are used literally. The reused
`coordinateDot_self_nonneg` is proved from the transparent finite sum
and square nonnegativity in `RightAngledCubicScalarKernel`, source
f55568035b8852a6284f1bd404c897758f30b301b1d85f38f09e4d9de365e6ee.
Both current source hashes match. This review does not recompile their
already successful modules or assert a fresh transitive Std object audit.

Exactly the five public declarations listed above have requested
`#print axioms` commands at lines104--108. There is no `sorry`, added
axiom, native decision oracle or conclusion-shaped structure. Before
root execution these remain REQUESTED prints; this review records no
actual axiom lists, source/object/compiler binding or runtime PASS.

The map accurately separates the remaining application obligations:
Kolpakov's metric equality and an isometry identifying the actual
equality facet with the displayed24-cell; actual outward unit facet
incidence; future-unit coordinates for the selected intrinsic point;
transport of heights/strict caps through the isometry; time-one form
of the other normals; and projection, face inheritance and published
facet/kissing inputs. Their role is visible in the paper: Musin and
Kolpakov force an actual regular equality facet, whose transported
strict heights contradict this coordinate lemma.

The kernel constructs neither the isometry nor the24-wall geometric
object. It does not instantiate `WrittenRegularCellGeometry` or prove
ideal dimension<=4 by itself. The map explicitly states those limits
and correctly marks source preparation uncompiled. No correction is
required before a separately admitted single-module bounded compile
and five-print audit. The six papers and current dimensions are unchanged.

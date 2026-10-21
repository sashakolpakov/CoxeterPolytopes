# Right-angled projection and height descent: assumptions map

**Bounded compilation and all eight actual axiom prints pass.** [RightAngledProjectionHeightKernel.lean](../formal/CoxeterProofs/RightAngledProjectionHeightKernel.lean) has SHA-256 `2fac8028e432a2104a58e423c65fc38c9e214fed5f52210e408aabe3f15f0e90` (282 lines,12333 bytes). Each printed declaration uses only `propext`, `Classical.choice` and `Quot.sound`, with no `sorryAx`. The [verification receipt](../results/publication_logic_verification_20261020.json) now records seven selected kernels and 62 actual axiom prints: the preserved six/53 plus nine for the separate [dimension58 shape kernel](finite_volume58_shape_kernel_assumptions_20261020.md), whose geometry/normalization inputs remain explicit. All seven have raw output and fourteen immutable compile/audit snapshots. SharpRA's core/full selections remain six jobs each with PDF; its actual geometric inputs remain explicit.

The source imports the verified [RightAngledCubicScalarKernel](../formal/CoxeterProofs/RightAngledCubicScalarKernel.lean), SHA-256 `f55568035b8852a6284f1bd404c897758f30b301b1d85f38f09e4d9de365e6ee`, for transparent finite sums and coordinate products. It uses the same generic linearly ordered commutative-ring interface, with no division, inverse, square-root or real-analysis library. No frozen module or manuscript was changed. Its written target is Lemma `lem:projection` and the strict height thresholds following equation `eq:iterate` in the [sharp right-angled manuscript](../manuscripts/right_angled_sharp_dimension_bounds.tex), SHA-256 `3c3ddd8540593004bcd5fb57832afc77ffc639268e66cb0ab29302b286654310`.

## Coordinate statements

A vector is a time coordinate paired with a spatial function `Nat → R`; only the first `d` spatial coordinates enter

```text
lorentz d x y = -(x.time*y.time) + coordinateDot d x.spatial y.spatial.
```

Symmetry and bilinearity are proved from the finite coordinate definition, rather than supplied as certificate fields. Set `v=q+t*n` and `w=lambda*v` using the explicit coordinate addition/scaling definitions.

| Declaration | Elementary inputs | Derived conclusions |
|---|---|---|
| `raw_projection_identities` | `q.q=-1`, `n.n=1`, `q.n=-t` | `v.n=0`, `v.v=q.v=-1-t*t` |
| `raw_wall_strict` | `t>0`, `q.m<0`, `n.m<=0` | `v.m<0` |
| `normalized_projection` | The three unit/incidence equations, `lambda>0`, `lambda*lambda*(1+t*t)=1` | `w.w=-1`, `w.n=0`, `q.w<0` |
| `normalized_wall_strict` | `t>0`, `lambda>0`, `q.m<0`, `n.m<=0` | `w.m<0` |
| `orthogonal_height_identity` | `n.m=0` and the normalizing equation | With `h=-q.m`, `h'=-w.m`, derive `h'=lambda*h` and `(1+t*t)*h'*h'=h*h` |
| `maximal_height_square_step` | Orthogonality, positive normalizer, its equation, and `0<h<=t` | `h'>0` and `(1+t*t)*h'*h'<=t*t` |

The norm equations and wall signs above are primitive coordinate data. No projected unit norm, all-wall feasibility, future-point predicate, improved height bound or dimension conclusion is an input field. The maximum condition is needed only for the last row. The wall statements quantify over an arbitrary wall, so a finite actual wall family supplies them pointwise without losing non-neighbor walls.

## Finite descent

`descent_step_strict` takes `a>=0`, the local recurrence `(1+a)*b<=a` and the previous invariant `k*a<1`, with natural `k`. It derives `(k+1)*b<1` from the polynomial identity

```text
(1+a)*(1-(k+1)*b)
  = (1-k*a) + (k+1)*(a-(1+a)*b).
```

The right side is positive. If the desired bound failed, both `1+a>=0` and `1-(k+1)*b<=0` would make the left side nonpositive. No division or reciprocal monotonicity is required. Nonnegativity of `b` is unnecessary for this stronger one-step implication.

`descent_strict` handles an arbitrary finite length `steps`. It assumes nonnegative `a j` only for `j<=steps`, and the local recurrence only for `j<steps`. Ordinary natural-number induction derives `j*a j<1` for every `j<=steps`; its initial case is `0<1`. Thus no infinite projection chain or desired terminal invariant is assumed. Taking `a j` to be the squared maximum height at each actual face, the preceding coordinate theorem supplies every local recurrence. One and three steps give the strict thresholds used by the two manuscript exclusions. The stronger displayed formula involving the initial height is not separately formalized here.

## Scalar radial gap

`radial_gap` takes an arbitrary positive scalar `p`, nonnegative heights `hi,hj`, positive spatial lengths `ri,rj`, the unit-wall equations `ri*ri=1+hi*hi` and `rj*rj=1+hj*hj`, the bounds `(p-1)*hi*hi<1` and `(p-1)*hj*hj<1`, and the all-wall inequality `ri*rj*z<=hi*hj`. It derives `p*z<1`. Multiplying the two strict height bounds proves `(p*hi*hj)^2<(ri*rj)^2`. Nonnegativity removes the squares, and the positive product `ri*rj` gives the strict conclusion. No desired radial gap, cardinality or dimension bound is an input. Taking `p=4` after three descents gives `4z<1`; taking `p=2` after one gives `2z<1`.

The [parsed actual axiom output](../results/RightAngledProjectionHeightKernel_axiom_summary_20261020.json) covers exactly eight qualified declarations: `raw_projection_identities`, `normalized_projection`, `normalized_wall_strict`, `orthogonal_height_identity`, `maximal_height_square_step`, `descent_step_strict`, `descent_strict` and `radial_gap`. The raw-wall strictness helper is a dependency of the printed normalized-wall theorem. The saved1057-byte [Lean output](../results/RightAngledProjectionHeightKernel_axiom_stdout_20261020.txt) matches both execution-ledger output hashes. First compilation took3.47s and the initial separate audit3.15s, each under the one-thread/2048MiB runner. The current six-job core run repeats the audit with the same1057-byte child output; immutable copies of the execution records are bound in the publication logic receipt. The object hash is `b2e265c8445767735f7ea149cfcf8344df074229bf62bafe158b0b482d7ef2a3`.

The ledgers bind the complete project-local object closure: the directly imported cubic kernel and its `RightAngledSharpDimensions` dependency. Their current objects and source identities were checked. Compiler/configuration/shared-library fingerprints bind the pinned Lean4.32.2 toolchain; a separately hashed transitive standard-library object closure is not claimed. The [independent source/API review](../results/review_right_angled_projection_height_kernel_resource_20261020.md) predates compilation and remains distinct from these actual execution records.

## Remaining geometric boundary

The real application must still supply positive normalizers (over the reals, `lambda=1/sqrt(1+t*t)`), actual finite unit facet normals, a relative-interior point and a maximum-height facet. It must justify the all-wall nonpositive Gram condition from right angles, and identify all intrinsic facet walls with the orthogonal ridge neighbors used in the height calculation.

The interpretation of unit timelike norm and negative pairing with the original future point as the future sheet remains a Lorentz-model foundation. The interpretation of the derived strict halfspace signs as relative interior of the actual facet, plus finite-volume/right-angle/entire-ideality face inheritance, remains geometric. These inputs contain no improved-height assertion. The source neither builds actual hyperbolic polytopes nor constructs the old `ProjectionGeometry` record, whose geometric incidence fields still require this bridge.

Radial unit-vector realization and its actual-facet indexing/injectivity, the tensor contractions, the regular24-cell opposite-normal argument and the named published facet/kissing/rigidity results are outside this module. The scalar height-to-product comparison is derived here, but no real radial vectors or finite code objects are constructed. No new dimension implication or end-to-end formalization is asserted. This metadata update ran no scientific worker, Lean or TeX job; the earlier root-owned compilation and audit are evidenced separately above.

The [fresh independent foundations review](../results/review_right_angled_sharp_foundations_adversarial_20261020.md) audits the written projection, face inheritance and code composition at these explicit geometric and published-input boundaries. The [facet/acuteness source audit](../results/review_sharp_ra_published_inputs_dufour_20261020.md) checks the mixed 27 consequence in Alexandrov as declared (Dufour's original remains paywalled) and Vinberg's original all-wall implication; the [ideal-branch audit](../results/review_sharp_ra_published_inputs_ideal_branch_20261020.md) checks Musin's bound and Kolpakov's facet minimum/metric equality in their originals. These are source reviews, not new Lean instantiations or executions.

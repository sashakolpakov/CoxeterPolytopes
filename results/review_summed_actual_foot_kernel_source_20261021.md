# Mathematical source review: SummedActualFootKernel

21 October 2026. **PASS for the conditional finite-coordinate kernel.**
This is a source/manual mathematical review, not a compilation or actual
axiom-output audit. No Lean, scientific or TeX job was run. No Lean,
manuscript or canonical source was edited.

Reviewed source: `formal/CoxeterProofs/SummedActualFootKernel.lean`,
318 lines, SHA256
`f4dcb86bc2fe13d310fa22d030ae6b5f4f99cbb30376e0127cb7873601e3d3ec`.
The inspected direct import is `RightAngledProjectionHeightKernel.lean`,
SHA256 `2fac8028e432a2104a58e423c65fc38c9e214fed5f52210e408aabe3f15f0e90`;
its finite-sum/coordinate dependency `RightAngledCubicScalarKernel.lean`
has SHA256 `f55568035b8852a6284f1bd404c897758f30b301b1d85f38f09e4d9de365e6ee`.
These are source bindings, not a transitive object/compiler receipt.

## Derived algebra and norm bound

`weightedSum` (lines31--32) and the imported `coordinateDot` are actual
finite coordinate sums: indices i<n and coordinates c<d. Weighted
pairings, bilinearity, subtraction squares and the full double Gram
expansion (lines34--121) are derived from those definitions. They are
not fields of an abstract inner-product record.

The essential `obtuse_weighted_norm_bound` (lines130--168) derives

    ||sum a_i n_i||^2 <= sum a_i^2

from a_i>=0, unit old normals and their nonpositive distinct products.
Its induction bounds the accumulated cross pairing by a finite sum of
nonpositive terms, scales it by the next nonnegative coefficient, and
adds the exact unit diagonal. Thus it covers ANY finite family,
including adjacent normals; no orthogonality or nonadjacency is assumed.
The diagonal upper bound is not supplied as an application premise.

The imported coordinate-square nonnegativity is itself proved by summing
nonnegative coordinate squares in a linearly ordered commutative ring.
Consequently `coordinate_energy_lower` (lines238--242) genuinely derives
2<x,W>-||W||^2<=||x||^2 by expanding ||x-W||^2>=0.

## Exact pressure, shift and final energy

Lines170--180 derive E=||w||^2=sum beta_i<w,omega_i> for the actual
weighted reference w=sum beta_i omega_i. No free scalar E or desired
energy value is assumed. Lines181--207 derive the single vector

    W=sum beta_i(omega_i-u_i n_i)=w-Z,
    sum beta_i T_i <= <x,W>, Z=sum beta_i u_i n_i,

by summing the individual foot comparisons with nonnegative beta_i.
There is no independence requirement between the individual feet.

The `dual` premise in lines209--236 is the linear identity
<w,n_i>=beta_i on the summed range. It gives exactly
<w,Z>=sum beta_i^2 u_i and then the shifted-square identity.
This is a reference/coordinate compatibility assumption, not the
desired shift norm, effective energy or whole-A20 conclusion. A concrete
application may derive it from old-root/dual-root biorthogonality; this
module does not instantiate that application.

For `summed_actual_foot_energy` (lines244--305), all primitive conditions
are explicit: beta>=0,u>=0, unit/obtuse old normals, the reference pairing,
and each individual foot inequality. Applying the independently derived
norm bound to a_i=beta_i u_i gives ||Z||^2<=sum beta_i^2 u_i^2.
The pressure, exact shifted square, reference energy identity and
coordinate-square inequality then give

    E+sum [2 beta_i(T_i-<w,omega_i>)
                       +(2u_i-u_i^2)beta_i^2] <= ||x||^2.

I re-expanded the final sum identity in lines277--293 and the concluding
algebra; all coefficients and signs match this formula. In particular,
dropping nonpositive off-diagonal contributions to ||Z||^2 improves the
energy lower bound. Adjacent feet cannot create an omitted negative
overlap. Neither u_i<=1 nor a nonadjacency hypothesis is needed for this
generic theorem; the expression may simply have negative quadratic
terms outside the application's 0<=u_i<=1 range.

All index hypotheses are restricted to i<n and, for obtuseness, j<n.
Values outside the summed range have no effect. No positivity of n or d
is needed: empty sums and coordinate sums retain their ordinary meaning.

## Match to the accepted whole-A20 proof and remaining boundary

The theorem matches equation(12) in accepted written source
`results/arithmetic_whole_a20_summed_actual_feet_exclusion_candidate_20261021.md`,
SHA256 `457debf95ed29488a90159ac2b9f9f8fe8780ebdd79e2ba8d50e6a72a07208a0`.
For example, an application can take the twelve central indices5..16,
their old unit normals and spatial duals, and coordinate dimension21.
For dotted indices set u_i=gamma_i^-2 and T_i to the actual dotted-foot
lower bound; for strong indices set u_i=0 and T_i=H_i. Then
H_i=<w,omega_i>, so the displayed kernel expression is precisely the
paper's summed effective energy, including its strong branches.

The application must still establish the coordinate realization,
concrete dual reference, strong-or-dotted inventory, actual feasible
future points and farthest comparisons. The radical bounds, explicit
beta coefficients, E=1.007658, internal error .005511, endpoint
positivity, ||x||^2<1 and the resulting numerical contradiction are
also outside this module. No A20 exclusion, field argument, curvature
estimate or dimension conclusion is silently formalized by it.

The source contains ten distinct `#print axioms` requests for its named
derived declarations, including the main theorem. I found no `sorry`,
`admit`, new axiom declaration or unsafe bypass in this source. Actual
elaboration and the printed dependency/axiom sets require the separate
bounded compile and axiom-audit receipt; this review supplies no runtime
or no-sorryAx execution claim. No mathematical source blocker was found.

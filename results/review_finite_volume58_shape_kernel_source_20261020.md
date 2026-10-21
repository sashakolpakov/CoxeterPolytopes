# Independent review of the dimension-58 shape assembler

20 October 2026. **Mathematical source and proposed primary-selection
replacement review: PASS.** No Lean/scientific/TeX job was executed by
this reviewer. No canonical source was edited. Compilation and actual
axiom-output verification remain separate.

Reviewed source:
`formal/CoxeterProofs/FiniteVolume58ShapeKernel.lean`, 12401 bytes/339 lines,
SHA-256
`40d12b5761a69470eee963a57ae9b6baf8c53ee120b83f05d11a4d98e1eebcec`.
Assumptions map:
`docs/finite_volume58_shape_kernel_assumptions_20261020.md`,
`48079a367758003e4a47fb2c61bee81b6c3858290717edbe9a7ce05cd97e81c3`.
The manuscript remains `c300db69…ca82a0`; current catalogue inspected was
`932269ff…5eec7b`. The mathematical target is manuscript lines 801–827,
as specified in `results/finite57_next_logic_obligation_audit_20261020.md`.

## Primitive data and derived local inequalities

`Link` contains only natural index multiplicities and a below/above split
`L(j)+R(j)=h(j)` for j<58. The functions G, E, V and U are actual finite
sums of those values, not independent fields (lines 100–124). Their
coefficients use natural subtraction only for `57-j` inside the bounded
range. Curvatures and subsequent differences use Int throughout.

`minimum_uncut` derives G<=U termwise from both minimum-coefficient
inequalities and the split (143–157). It does not take this inequality as
a hypothesis. `coefficient_gap` is valid for arbitrary signed coefficients
a<=b and multiplicity h>=1: it multiplies `1<=h` by the nonnegative gap
b-a and rearranges. The two computed gap sums are exactly

```text
sum(j-min(j,57-j)) = 841
sum(57-j-min(j,57-j)) = 841,
```

over j=0,...,57. Both are the odd sum 1+3+...+57. Thus positive
multiplicities give the full 841 bound, not merely a constant identity.

At an all-below extremum, U=E follows directly from L=h,R=0, and
`edge_gap_from_positive` supplies the bound (167–192). At an all-above
extremum, `all_above_gap` directly uses the **other** gap coefficients
and L=0,R=h (194–211). It neither assumes nor derives U=E. Consequently
there is no hidden symmetry or degree-balance requirement in either
extremal lower bound.

## Two extrema, finite/ideal split, and global incidence

`two_extrema_curvature` requires both positions inside the complete
0,...,n-1 domain and explicitly requires lo!=hi (214–249). Its pointwise
delta comparison distinguishes lo, hi and all other positions. Each
selected position contributes 841 and all other slacks are nonnegative.
The two delta sums are individually one; hence their total is 1682.
Multiplication by five gives the exact doubled correction 8410. This
does not count either extremum a second time or require it to be finite.

For a finite link, h=1 gives G=812, E=1653, V=58 and
`5G-2E=754=2*377` (251–274). The Boolean tag is used once at each bounded
vertex position: `finiteCount` counts true positions, while
`idealCurvature2` sums precisely the complementary positions. The split
theorem derives the finite term from h=1 rather than assuming the final
curvature value (276–298). Empty ideal domains are allowed. Curvature
at other links can be negative; no truncation or positivity premise is
silently imposed.

The main theorem takes exactly the two primitive global equalities

```text
I = sum E_v,
2*f2 = sum U_v.
```

It derives `sum(5U-2E)=10f2-2I` and combines the preceding theorems to
obtain

```text
754*finiteCount + idealCurvature2 + 8410 <= 10*f2-2*I.
```

Neither a shape lower bound, the positive correction, a charge budget,
terminal surplus nor any dimension exclusion is a supplied field
(300–321). The proof therefore adds internal local-to-global substance
absent from the earlier conditional `WrittenGeometry.shapeLower`.

`balanced_curvature` separately requires `2E=57V` and only then identifies
`5G-2E` with `5G-57V` (323–327). This matches the paper's doubled ideal
curvature. It does not contaminate the all-above proof with a false
arbitrary-data equality.

## Assumptions and integration limits

The map accurately retains actual oriented-link/index realization, uncut
edge assignment, complete vertex enumeration, two actual generic extrema,
finite/ideal tags, product-link coefficient positivity, and the published
incidence/index identifications as explicit application obligations.
These counts must be instantiated from the actual chamber. Compilation
alone will not instantiate `ManuscriptFinite57.WrittenGeometry.shapeLower`.
The finite-count arithmetic is proved, but actual geometric count
interpretation and any Nat/common-denominator adapter remain separate.
The nine literal axiom printers match the intended source declarations;
this review makes no claim about their future output.

The proposed integration replaces **only the primary finite57** standalone
`CoxeterProofs.ShapeCurvatureArithmetic` root with this new module after
successful compilation/audit. Independent source search found no Lean
module importing the old file or referencing its declarations. The old
file's 870 slack, finite 783/2 and correction 4350 are ambient-59/link-58
arithmetic; they do not supply the current rank58/377/4205 proof. Its
global contradiction already assumes its shape inequality. No claimed
rank58 verification is lost by replacing this isolated primary slot.

Keep the old source/evidence and the separate historical `finite57_atlas`
catalogue selection unchanged. The replacement then preserves the primary
paper's nine-root/11-core/123-full job counts while adding the actual
dimension-58 assembler. This source review authorizes no claim of combined
runtime PASS, end-to-end geometry, or whole transitive regeneration.

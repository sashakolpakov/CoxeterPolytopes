# Next finite57 internal logic obligation: derive the shape lower bound

20 October 2026. Fresh source reading only; no scientific, Lean or TeX job,
and no canonical source edit. This note proposes one bounded formalization
target, not a new manuscript or a claim of implementation.

Current sources read:

* `manuscripts/finite_volume57.tex`:
  `c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0`.
* `docs/finite_volume57_lemma_map.md`:
  `24dbd895d7607f833505b608919c601e57501bff9e80feca88f36c8acae9c948`.
* `formal/CoxeterProofs/ManuscriptFinite57.lean`: `c680d5d0…70002dc`;
  `ShapeCurvatureArithmetic.lean`: `f54ff77d…52d1f2`;
  `FiniteVolumeTailKernel.lean`: `b8f03fd4…39daa`.

## Recommendation and missing substance

The highest-value small next theorem is an **actual shape-curvature
assembler for dimension 58**, deriving the strict 4205 correction from
primitive finite oriented-link counts. Target manuscript lines 801–827 and
lemma-map item 18, currently lines 407–410.

`ManuscriptFinite57.WrittenGeometry.shapeLower` still takes the entire
377/4205 lower bound as a field (lines 45–47). `FiniteVolumeTailKernel`
also consumes a shape inequality. Existing `ShapeCurvatureArithmetic`
has the older ambient-59 constants (link dimension 58, slack 870,
correction 4350), and its final contradiction likewise assumes a global
shape inequality. None of the six newer kernels derives the actual
dimension-58 local-index-to-global-shape implication.

Another final subtraction theorem would repeat existing work. The proposed
assembler instead discharges a specific substantial internal field while
leaving honest geometric incidence/realization inputs. A single small
module is enough; do not expand it into a general polytope library.

## Concrete minimal data and theorem

Use the 58 indices j=0,...,57. For each actual vertex position v supply
natural counts `h_v(j), L_v(j), R_v(j)` with

```text
L_v(j) + R_v(j) = h_v(j).
```

Here h counts oriented link vertices of index j, L counts those below the
induced level section, and R those above it. Define, rather than assume,

```text
E_v  = Σ_j j*h_v(j)
G_v  = Σ_j min(j,57-j)*h_v(j)
E0_v = Σ_j [j*L_v(j) + (57-j)*R_v(j)]
K2_v = 5*(G_v : Int) - 2*(E_v : Int).
```

The finite global vertex domain must distinguish positions, with a
finite/ideal tag. Require `h_v(j)=1` at every finite vertex. Choose two
**distinct** positions p,q satisfying

```text
L_p(j)=0        for every j       -- global minimum: all above
R_q(j)=0        for every j       -- global maximum: all below
h_p(j)>=1, h_q(j)>=1 for every j.
```

Only the selected two links need positivity for this scalar theorem;
actual simplex/product links have it at every vertex. Do not require
either selected position to be finite.

The one global geometric/index input is the exact, doubled identity

```text
D2 = Σ_v [5*(E0_v : Int) - 2*(E_v : Int)],   D2 = 2*D5(P).
```

Then derive, with v_f the number of finite-tagged positions,

```text
754*v_f + Σ_(v ideal) K2_v + 8410 <= D2.
```

This is exactly twice the manuscript's shape lower bound after the
curvature interpretation below. There must be no field named or
equivalent to `shapeLower`, `G+841<=E0`, `4205 correction`, or the desired
dimension exclusion.

## Proof that supplies the new internal logic

1. From the split identity and `min(j,57-j)<=j,57-j`, prove termwise and
   then by finite summation that `G_v<=E0_v` at every vertex. This is the
   orientation-assignment minimum estimate, not a supplied inequality.
2. At q (all below), derive
   `E0_q-G_q = Σ_(j=29..57)(2j-57)h_q(j) >= 841`.
   At p (all above), derive
   `E0_p-G_p = Σ_(j=0..28)(57-2j)h_p(j) >= 841`.
   Both numerical coefficient sums are `1+3+...+57=841`; kernel `decide`
   can check these small constants. Positivity and the finite-sum bound
   prove the inequalities for arbitrary multiplicities.
3. Sum nonnegative local slack and reserve the two distinct positions
   exactly once. The total slack is at least 1682, even when an extremum
   is ideal. Consequently the difference between D2 and `Σ_v K2_v` is
   at least `5*1682=8410`.
4. For a finite link h=1, compute `E=1653`, `G=812`, and
   `K2=5*812-2*1653=754`. Split the finite/ideal tagged sum to obtain the
   displayed theorem. The theorem also covers an empty ideal domain.

An important normalization subtlety: for arbitrary h, an all-above link
has `E0=Σ(57-j)h_j`, which need not equal `E=Σj h_j`. Therefore one must
not infer `E0=E` merely from its sign. The direct lower-half estimate in
step 2 avoids any hidden symmetry or balance premise. Alternatively actual
regular-link balance would justify equality, but it is unnecessary for
the strict correction proof.

## Precise geometric, script and internal boundaries

**Published geometric inputs.** Khovanskii's exact incidence/index
identities give the global expression and identify the oriented index
distribution; the finite/cusp structural classification identifies the
simplex/product-of-simplices links. These are the paper's stated inputs at
lines 804–806, 813–814 and 857–859. This source-only review has not checked
the cited publications themselves or their full hypothesis applicability.

**Actual geometric incidence bridge.** The application must identify E
with the actual link edge count (each oriented edge has one upper endpoint),
and E0 with the uncut edge count via the below-upper/above-lower assignment.
It must supply a complete actual vertex domain, its tags, and two distinct
generic global extrema with uniform signs. These are primitive incidence
and sign data, not numeric shape conclusions. The assignment argument at
paper lines 809–814 and the generic-extrema identification remain internal
geometric interpretation; do not relabel them as externally proved merely
because the scalar theorem is conditional.

**Cusp curvature interpretation.** For an actual 57-dimensional simple
link, `2E=57V`. Together with `V=Σh=∏(r_i+1)`, this identifies
`K2=5G-57V=2κ` with the manuscript's cusp curvature (lines 649–654).
Regular-link degree counting and the product h-vector identification must
be explicit in this adapter. Symmetry/positivity of product coefficients
are elementary internal facts of that identification; the scalar theorem
does not establish actual product-link realization. No numeric κ lower
bound is required to prove the Int theorem.

**Allowed script facts.** The finite charge cap 377 and terminal surplus
36 may remain the existing exact-script inputs when applying the resulting
shape theorem to payment subtraction. They are not needed in the proposed
shape proof itself. Do not supply the global shape bound as a script fact:
the endpoint's checks of 812/841/4205 do not prove the universal assembler.

**New internal proof content.** The finite split-to-minimum argument,
both extremal weighted bounds, distinct-position summation, finite-link
constant and tagged global assembly are all derived. This replaces a
conclusion-shaped field in the current conditional interface.

Keep the main theorem in Int to avoid truncated curvature subtraction.
If an adapter to the existing Nat interface is added, choose a common
positive charge scale `2s`: multiply the doubled theorem by s. Encoded
deficit is `s*D2=2s*D5`, ideal curvature is `s*K2=2s*κ`, finite coefficient
is `754s=377*(2s)`, and correction is `8410s=4205*(2s)`. Justify
nonnegative Nat encodings separately. This adapter is optional, not a
reason to broaden the initial module.

## Stopping rule

One such assembler advances the requested internal proof logic beyond the
six kernels. It does not prove actual Lorentz/polytope realization, the
whole simultaneous local atlas, or the full higher tail. After this one
step, prioritize applying its primitive-data interface and auditing the
published Khovanskii/Nikulin hypothesis boundaries rather than adding
another scalar wrapper or unrelated manuscript/kernel family.

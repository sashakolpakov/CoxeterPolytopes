# Dimension-58 shape counting: scope of the formal bridge

The target is the shape-curvature lemma in
[finite_volume57.tex](../manuscripts/finite_volume57.tex), lines 795–831,
source `c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0`.
The [new source](../formal/CoxeterProofs/FiniteVolume58ShapeKernel.lean),
SHA-256 `40d12b5761a69470eee963a57ae9b6baf8c53ee120b83f05d11a4d98e1eebcec`,
passed its first bounded compilation in **0.59 seconds**, followed by a
**0.55-second audit of all nine printed declarations**. Only standard axioms
occur, with no `sorryAx`. The [actual output summary](../results/FiniteVolume58ShapeKernel_axiom_summary_20261020.json)
and [execution receipt](../results/finite_volume58_shape_kernel_initial_execution_20261020.json)
bind the source, object, immutable ledgers and inspectable raw child output.
The unchanged source emits eight harmless deprecation/unused-simp warnings.
This is supporting verification within the existing paper, not another
manuscript or a new dimension bound.

## Counts supplied by geometry

An actual vertex has a simple 57-dimensional link. Fix a generic linear
function on that link, with the induced section at level zero. For each
index `j=0,...,57`, let `h(j)` count link vertices of index `j`, and split
that count into `L(j)` below the section and `R(j)` above it. The primitive
count relation is `L(j)+R(j)=h(j)`.

The scalar source defines the following sums rather than assuming bounds
on them:

```text
E  = sum j*h(j)
G  = sum min(j,57-j)*h(j)
E0 = sum [j*L(j)+(57-j)*R(j)].
```

The application identifies `E` with the actual link edge count. Every
oriented edge has exactly one upper endpoint. To identify `E0` with uncut
edges, assign an uncut edge below the section to its upper endpoint, and
one above the section to its lower endpoint. A below vertex of index `j`
receives exactly its `j` incoming edges; an above vertex receives exactly
its `57-j` outgoing edges. A cut edge is assigned to neither endpoint.
These are incidence interpretations of the supplied counts; the scalar
kernel does not construct a geometric link or prove transversality.

The complete global vertex domain has distinct positions for the generic
minimum and maximum. At the minimum every link vertex is above the
section (`L=0`); at the maximum every link vertex is below it (`R=0`).
Only those two links need the premise `h(j)>=1` for the strict correction.
For actual Coxeter links this follows from the simplex or
product-of-simplices h-polynomial. Identification of that polynomial and
its positive coefficients remains an explicit application obligation.

## Internal deductions

The split and the two inequalities for `min` imply `G<=E0` term by term.
At the all-below link, its slack is the weighted sum with weights
`j-min(j,57-j)`; at the all-above link, use
`57-j-min(j,57-j)`. Each sum of weights is `841`. Positivity of the
corresponding `h(j)` therefore proves `G+841<=E0` at either extremum.
This argument needs no h-vector symmetry and does not assume the desired
841 slack or infer `E0=E` from an all-above sign pattern.

Distinct positions allow the two contributions to be reserved exactly
once in the global sum. In doubled integral normalization the result is

```text
sum_v (5*G_v-2*E_v) + 8410
  <= sum_v (5*E0_v-2*E_v).
```

All differences are in `Int`, so a negative local curvature cannot be
silently truncated. At a finite vertex `h(j)=1`, the source checks
`G=812`, `E=1653` and `5G-2E=754=2*377`. Splitting the complete domain into
finite and ideal positions gives

```text
754*v_f + sum_(v ideal) (5*G_v-2*E_v) + 8410 <= D2.
```

This includes the case with no ideal positions. No finite charge cap,
terminal surplus or final dimension exclusion is required for this shape
counting argument.

## Exact application to the manuscript

The geometric incidence identities are
`I=sum E_v` and `2*f_2=sum E0_v`. They identify
`D2=2*(5*f_2-I)` with the right-hand side above. For a simple
57-dimensional link with `V` vertices, degree counting gives `2E=57V`.
Its h-index distribution identifies `G` with the manuscript's `G`, so
`5G-2E=5G-57V=2*kappa` at an ideal vertex. These are the precise
interpretation obligations; an arbitrary natural index row need not
satisfy degree balance.

The supplied doubled identity then gives exactly
`D5>=377*v_f+sum kappa+4205`. A bridge to a natural-number charge interface
must also supply a common denominator and justify the nonnegative
encodings. The unscaled `Int` theorem supplies no such encoding by fiat.
The older `ManuscriptFinite57.WrittenGeometry.shapeLower` field is not
automatically instantiated merely by compiling the new source.

The [mathematical specification](../results/finite57_next_logic_obligation_audit_20261020.md)
records why this derives an internal assumption missing from the preceding
six kernels. The [independent source review](../results/review_finite_volume58_shape_kernel_source_20261020.md)
binds the unchanged compiled source and the archived pre-execution version
of this map. The [primary-source audit](../results/review_finite57_published_shape_inputs_20261020.md)
checked Khovanskii's original Theorems 2, 6 and the proof of Theorem 9,
including the finite-volume ideal-cusp scope in Section 6. The source
review and bibliographic audit make no end-to-end geometry claim.

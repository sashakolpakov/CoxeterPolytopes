# Full-Picard carrier ingestion queue after the `X2 H10/37` transport

## Scope

This is a bounded audit of the maintained higher-dimensional catch ledger,
not a literature-wide census.  It asks which certified finite-volume rows in
`H^n`, `n>=10`, already have exact local wall/Gram data but still lack a
machine-frozen **very-general full-Picard** algebraic carrier.  A special
surface with the chamber lattice only as an invariant sublattice does not by
itself close that stronger gate.

Rows already carrying a complete-nef surface whose full Picard lattice is
the chamber lattice, an equation-level full Picard lattice, or one of the
recently frozen generic K3 embeddings were not put back into the queue.  Odd integral carrier lattices were not silently called
K3 Picard lattices, and rank above 20 was not treated as an ordinary complex
K3 problem.

## Closed in this pass: `X2-RAT-C2-H10-37`

The highest-yield row was the 37-wall rational-involution section on
Shimada's `X_2`.  Its exact source package already proved

```text
dimension / facets       H10 / 37
fixed lattice            U+E8(-2)+<-2>
signature / determinant  (1,10) / 512
wall span                index 1
mirror types             6*(-2,div1)+21*(-2,div2)+10*(-4,div2)
Vinberg                   3236 finite + 261 ideal, 24350 edges, 0 failures
```

The missing step was not another polytope computation.  The existing
quadratic-base-change theorem gives an explicit nine-parameter K3 family
whose very-general full Picard lattice is the same
`F=U+E8(-2)+<-2>`.  The new composition certificate replays the pinned
determinant-one isometry from the geometric standard `F` basis to
`NS(X_2)^iota`, transports all 37 walls, and checks every primitive
reflection.  Their 37-by-11 coordinate matrix has Smith diagonal `1^11`.

The result is frozen in the
[`H10/37` full-Picard carrier theorem](x2_h10_37_ue8m2a1_full_picard_carrier.md),
with [certificate](../data/enumerations/x2_h10_37_ue8m2a1_full_picard_carrier.json),
[certifier](../scripts/certify_x2_h10_37_ue8m2a1_full_picard_carrier.py), and
[tests](../tests/test_x2_h10_37_ue8m2a1_full_picard_carrier.py).

## Second closure in this pass

### `K3-D4-H12/17`: coordinate genericization closed

This is the cleanest next row.  The local archive
`k3_rank15_order3_h12_fold.json` is already green and contains the complete
17-wall orbit Gram, the integral orbit reflections, and two independent
finite-volume proofs:

```text
fixed lattice       U+E8(-1)+A2(-1)+A1(-1)
rank / determinant  13 / 6
volume              93 finite + 3 ideal + 0 hyperideal
Vinberg edges        623, zero endpoint or cusp failures
maximality           full reflection chamber, index 1
```

The lattice is the primitive kernel of the order-three triality action on an
actual K3 Picard lattice.  Therefore it is already primitively K3-embeddable.
What is still missing is one composition certificate that reconstructs a
saturated fixed basis, expresses all 17 orbit normals in it, verifies their
span, and then takes a very-general period with full `NS` equal to that
rank-13 lattice.  This does **not** require realizing triality as an
automorphism of the returned surface: the output is a divisor-reflection
carrier.  An equation-level upgrade would additionally require a seven-moduli
family, naturally suggested by an elliptic frame with root stock
`E8+A2+A1`.

Clear completion test:

1. freeze the 15-to-13 saturated kernel matrix;
2. transport all 17 oriented walls and prove reflection integrality on the
   full rank-13 lattice;
3. prove index-one wall span or record the exact index;
4. give an explicit primitive embedding in `Lambda_K3`;
5. replay `93+3+0` and the literal 623-edge Vinberg gate;
6. state no triality-automorphism or full-nef claim unless separately proved.

Pinned source SHA-256:
`54fea0698d48241cd50d6654ccc7ad687e227e9fa9d8908e011e7455782e097b`.

All six completion tests now pass.  The determinant-one standard
decomposition, 17 transported walls, index-one span, visible-unit-minor K3
embedding, and both volume lanes are frozen in the
[`K3-D4 H12/17` full-Picard carrier](k3_d4_h12_full_picard_k3_carrier.md).
The output is a seven-dimensional period family, not yet an explicit
equation family.

## Third closure in this pass

### `K3-S3-C2-H10/13`: third chamber on the same equation family

The complete local survivor in
`k3_rank18_19_all_subgroup_completion.json` has 13 facets, 34 finite
vertices, four cusps, zero hyperideal rays, fixed lattice
`U+E8(-2)+<-2>`, and facet mirrors generating the saturated lattice.  It is
the known full `W_r(11,9,1)` chamber.

The missing machine gate is now closed.  An explicit determinant-one
`GL(11,Z)` matrix identifies the archived fixed Gram with the geometric
standard `F` Gram.  Its construction splits the divisibility-one isotropic
pair

```text
v=[0,0,0,-1,-1,-2,-2,-3,-2,-1,1]
w=[0,-1,0,-1,-1,-2,-2,-3,-2,-1,1]
```

and solves the remaining exact rank-nine definite isometry.  All 13
transported roots are primitive and integral, span `F` with index one, and
have type split `2*(-2,div1)+2*(-2,div2)+9*(-4,div2)`.  A square-358 class
strictly orients every wall.

The volume statement was freshly recomputed from the transported facet Gram,
not merely copied: exact cone enumeration gives `34+4+0`; the literal audit
has 208 rank-nine edges, all with two endpoints; and all seven optimized
affine-neighborhood audits pass.  Thus the same explicit family now carries
three distinct chambers with 13, 21, and 37 facets.  The special `S3` action
and quotient-nef interpretation are still kept separate from the
very-general family.

Pinned source SHA-256:
`148ce106b968709c94b31ac8bed13b17269861ea6051eeba68b5bc8b89e7e5c4`.

The frozen package is the
[`H10/13` full-Picard carrier](k3_s3_c2_h10_13_ue8m2a1_full_picard_carrier.md),
with [certificate](../data/enumerations/k3_s3_c2_h10_13_ue8m2a1_full_picard_carrier.json),
[certifier](../scripts/certify_k3_s3_c2_h10_13_ue8m2a1_full_picard_carrier.py),
and [tests](../tests/test_k3_s3_c2_h10_13_ue8m2a1_full_picard_carrier.py).

## Fourth closure in this pass

### `GHESS-C2-H10/27`: generic period carrier closed

The Hessian-quartic transposition archive already gives a saturated even
rank-11 fixed divisor lattice of determinant 256, 27 primitive integral
facet reflections, and literal Vinberg data
`895 finite + 67 ideal`, 6,221 elliptic edges, with zero failures.  Its
normalized chamber is the known `D16-D4-C7-H10/27` polytope.

That carrier step is now closed by an explicit discriminant-glue theorem.
The fixed lattice has

```text
A_S=(Z/4)+(Z/2)^6
q_S=<3/4>+u(2)^2+v(2).
```

The rank-16 Hessian Picard Gram and
`T_H=U+U(2)+A2(-2)` have an explicit anti-isometry on all 48 discriminant
classes.  Graph gluing gives a displayed even unimodular signature-`(3,19)`
lattice, and the induced 22-by-11 fixed-lattice embedding is integral,
isometric, and primitive by Smith diagonal `1^11`.  Thus a nine-dimensional
period family has very-general full `NS=S_H`.

All 27 oriented mirrors reconstruct in this full Picard lattice, every
reflection is integral, and their coordinate matrix also has Smith diagonal
`1^11`.  Fresh literal Vinberg gives `895+67`, 6,221 edges, and zero edge or
cusp failures.  Independent GMP double description gives exactly
`895+67+0`, so no hyperideal ray is being hidden by the diagram test.

Pinned source SHA-256:
`c637c54385b33f1f116fff8b9d3cb833c08f48e22b78b8cd809b6d25a04532fb`.

The frozen package is the
[`GHESS H10/27` full-Picard carrier](ghess_h10_27_full_picard_k3_carrier.md),
with [certificate](../data/enumerations/ghess_h10_27_full_picard_k3_carrier.json),
[certifier](../scripts/certify_ghess_h10_27_full_picard_k3_carrier.py), and
[tests](../tests/test_ghess_h10_27_full_picard_k3_carrier.py).  It is
period-level, not equation-level, full-nef, or a generic surface-automorphism
theorem; the abstract chamber remains the known `D16-D4-C7-H10/27` polytope.

## Fifth closure in this pass

### The `A20` pair: `H12/19` and `H11/20`

The last ordinary ingestion job is now closed jointly.  Starting from the
same singular Picard-rank-20 K3 witness, the exact glue reconstruction gives
primitive embeddings of both saturated fixed kernels in one marked K3
lattice.  Very-general periods therefore give full-Picard families of
dimensions seven and eight:

```text
                         H12/19                    H11/20
rank / determinant       13 / 256                  12 / -1024
wall count / span        19 / index 1              20 / index 1
finite + ideal + hyper   527 + 25 + 0              311 + 29 + 0
Vinberg edges             3846                      2229
```

Every primitive wall reflection is integral on the full Picard lattice.  An
exhaustive discriminant-module audit gives stable integral subgroup orders
`[1,2]` in each case, but stable even orders only `[1]`: the unique proper
index-two overclosure is odd.  Hence each displayed kernel is the unique even
reflection-stable carrier, so the denominator is forced before the period
theorem is invoked.

The frozen package is the
[`A20` paired full-Picard theorem](a20_pair_full_picard_k3_carriers.md), with
[certificate](../data/enumerations/a20_pair_full_picard_k3_carriers.json),
[certifier](../scripts/certify_a20_pair_full_picard_k3_carriers.py), and
[tests](../tests/test_a20_pair_full_picard_k3_carriers.py).  The source
order-four automorphism and involution are primitive-embedding witnesses; no
automorphism is claimed on a very-general returned surface.  These are mixed
divisor-reflection subsystems, not full-nef or new abstract-polytopal claims.

Frozen payload SHA-256:
`2d6c68b8aed6fbe7c39b527d94a23a1b7e466504d3a6097b13536f9edc96cbfe`.

## Structural blockers, not ordinary ingestion jobs

- `BOR-RM-H17/960`: the complete chamber is machine-known, but the complex
  K3 Picard realization is excluded, and the maintained `K3^[2]` and
  `K3^[3]` full-wall interpretations are also excluded.  Progress requires a
  genuinely different algebraic carrier class or a weaker subsystem notion,
  not another Gram ingestion.
- `BOR-D4-H21/210`: rank 22 is beyond the complex-K3 Picard bound.  It already
  occurs as the full Picard lattice of the characteristic-two
  Artin-invariant-one supersingular K3, so the live question is a different
  characteristic-zero or higher-dimensional carrier, not existence of an
  algebraic carrier simpliciter.
- Odd reflection-stable maxima in several rows cannot be Neron--Severi
  lattices of complex K3 surfaces.  Their even source lattices, even
  overclosures, or scaled duals must be tested individually; passing to an
  odd maximal group is not a legitimate shortcut.

## Operational order

With `K3-D4-H12/17`, the 13-wall common-family transport, the Hessian
27-wall genericization, and both high `A20` sections now closed, this bounded
ordinary ingestion queue is empty.  The next positive-first work is no
longer Gram transport: it must either find a new complete chamber, produce an
equation/effectivity upgrade of a period carrier, or change carrier category
for the structural blockers above.

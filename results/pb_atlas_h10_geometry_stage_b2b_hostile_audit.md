# Hostile audit: PB-Atlas H10 generated geometry, Stage B2b

**Verdict (2026-08-10): PASS after two substantive certificate repairs and
one scope repair.**  On the one structurally selected plane packet, the
incidence component, collision exclusion, exact strict transforms, all thirty
physical facets, and strong-Coble conclusion survive an independent replay.
The resulting surface is over `Q/C`, so its canonical K3 cover gives the
carrier-class bound `n <= 19`.  This is neither a second positive packet nor a
global raw-polytope existence theorem.

Audited result:
[`pb_atlas_h10_geometry_stage_b2b.md`](pb_atlas_h10_geometry_stage_b2b.md).

## Pinned audited artifacts

| artifact | SHA-256 |
|---|---|
| `data/enumerations/pb_atlas_numerical_contraction_stage_a.json` | `6eeb5726cf991e3e7f014429e7a5eb88a7ad4101376b1ced1ef32d30f4b6dcc4` |
| `data/enumerations/pb_atlas_numerical_cluster_stage_b1.json` | `c645ed67182b765717fd2fad4c1d89c7c97a332ad2d4b10ef31c04d084fa28c1` |
| `scripts/pb_atlas_cluster_chart_stage_b2a.py` | `2b0d534c6e014f3e15e85a200ceaad029160aa33a10ba8b5d6f80cde76904cbf` |
| `scripts/pb_atlas_geometry_stage_b2b.py` | `c60b92990e0f10d3e8bbf6646614288246265b089d3e93d4f3846db273fb4302` |
| `scripts/certify_pb_atlas_h10_geometry_stage_b2b.py` | `a47113454ea79aa7fbf47f15481aadbd87b9a6c4b0c02a6fa6abe0cf1e966177` |
| `data/enumerations/pb_atlas_h10_geometry_stage_b2b.json` | `cc0e88a3a4eba07873f69b9b5b3372b60ca69cbadece812a870d4d7d96bc19b0` |
| `results/pb_atlas_h10_geometry_stage_b2b.md` | `884f6755ca2202a241dfd38c90f9ce780ae0e65b7e578cda47d9a48448419053` |
| `tests/test_pb_atlas_h10_geometry_stage_b2b.py` | `0d9e572b69d8ad109afeb4ab94cb1b70046d8d9ef73a5174a5c0eb49365e810d` |
| `tests/test_pb_atlas_h10_geometry_stage_b2b_hostile_audit.py` | `74650e884a165fb2c90fd70d964e936dfa69064648cb6449fa48502d8c46ef53` |
| `results/canonical_orthogonal_k3_cover_dimension_theorem.md` | `d57ed2005df38b76dcc26c7876eb9e6532ba60b158632b21605bb5af1be3ea5c` |
| `results/reflection_maximal_coble_exact_nef_return.md` | `cfd6e175165244b5a4df6e94d888fbe7803d10b0ae7223dfd0c111a663eb67be` |

## 1. Incidence was rebuilt from the cluster, not replayed from B2b

The hostile test starts from the Stage-A endpoint programs, recompiles their
Stage-B1 packets, and independently finds exactly one `P2` packet with
proximity-component sizes

```text
[2,2,2,1,1,1,1].
```

It does not import either B2b module.  It fixes the deterministic four-point
frame and writes a line as `aX+bY+cZ`.  A proper multiplicity contributes
evaluation at the proper point.  A free-center multiplicity contributes the
directional row

```text
A chart: du + t dv,       B chart: t du + dv.
```

These direct rows regenerate all nine line-root determinants, up to a nonzero
rational unit.  Exhausting the `2^3=8` products of the two direction charts
leaves exactly

```text
center 3: A,  center 6: B,  center 9: B.
```

This is completeness only over the selected proper affine chart and its
declared general-position `PGL3` open.  It does not claim to exhaust the other
proper affine patches or gauge quadruples.  On that stated open, the two
blowup charts cover each exceptional `P1`, so the eight-way direction census
is complete.

The independent Gröbner basis is

```text
q5_0-q8_1, q5_1-1, q7_0-1, q7_1-q8_1, q8_0-q8_1,
(q8_1-1)(q8_1+1), t3, t6, t9.
```

It has precisely the two reduced closed points in the certificate.  Direct
projective normalization finds no proper-center collision at `q8_1=-1`.  At
`q8_1=+1`, positions `4,5,7,8` coincide and give the six recorded collision
pairs.  Both points have nonzero selected affine coordinates and all four
gauge minors; the positive component fails only the ordered-configuration
collision open.  There are no satellite centers, and the three free centers
have distinct parents, so no untested same-exceptional collision remains.

## 2. Every plane section survived an independent strict-transform replay

The hostile jet oracle expands homogeneous monomials directly in the seven
proper local charts.  It takes Taylor coefficients below the declared order,
then substitutes each of the three first-near blowup maps and divides by the
declared exceptional power.  It imports neither the production low/high
splitter nor the B2b exact-order routine.

This gives, independently,

```text
|-K|:   10 x 10, rank 10, nullity 0;
|-2K|:  30 x 28, rank 27, nullity 1.
```

For the boundary and all nineteen nonexceptional, nonboundary facets, the
archived polynomial spans the independently computed one-dimensional kernel.
Transporting each whole polynomial through the chronological blowups gives
the declared order exactly at every one of the ten centers.  In particular,
a zero declared multiplicity is also checked to be an actual zero order; an
unrecorded exceptional component cannot hide there.

Production uses Singular's `absFactorize`.  The hostile suite reruns every
nonlinear polynomial with the alternate `absFactorizeBCG` algorithm.  Each
has one nontrivial absolute factor of multiplicity one.  As a negative
control, `(X+Y)^2` has multiplicity two and is rejected.

## 3. Facet partition, smoothness, and the strong-Coble conclusion

The ten contracted wall indices equal the ten reverse-blowup source indices.
More strongly, for each center the hostile replay verifies

```text
wall exceptional coefficients = strict exceptional coordinates,
wall multiplicities            = their coefficientwise negatives.
```

Those ten indices, the nineteen other generated curve indices, and the one
boundary index are pairwise disjoint and partition `0,...,29`.  The boundary
root has class

```text
3H-E1-...-E10,
```

whereas the constructed physical boundary has exactly twice that class.  Its
intersection data are `B^2=-4`, `K.B=2`, and hence `p_a(B)=0`.  Direct use of
`d^2-sum(m_i^2)` and `-3d+sum(m_i)` similarly gives arithmetic genus zero for
all nineteen other plane strict transforms.

Absolute integrality plus the exact strict-transform class makes each such
curve an integral Cartier curve on the smooth blowup.  In characteristic
zero,

```text
p_a = genus(normalization) + sum of delta invariants.
```

Thus `p_a=0` forces smoothness and rationality; a separate Jacobian guess is
not needed.  The reverse-blowup exceptional curves are already smooth
rational curves.  Combined with `h0(-K)=0`, `h0(-2K)=1`, and the reduced
smooth boundary, this proves the stated strong-Coble surface.

## 4. K3 pushback and exact scope

For the double cover branched in `B=-2K`, the canonical bundle formula gives
`K_X=0`.  The branch is smooth and nonempty, so the cover is a smooth connected
surface.  Moreover `pi_* O_X = O_S + O_S(K_S)`; rationality and
Riemann--Roch give `chi(O_X)=2`, and Serre duality then gives `q(X)=0`.
Thus it is K3.  Pullback of `Pic(S)` is injective because push-pull is
multiplication by two, and it doubles the intersection form.  Therefore a
rank-`n+1` marked carrier injects into `NS(X)`.  Over `C`, `rho(X)<=20`, giving

```text
n <= 19.
```

This argument needs carrier existence but not equality with the nef cone.
The current chart and factorization compiler is over `Q`, with its geometric
conclusion over `C`.  The valid `n<=21` characteristic-not-two analogue is
conditional on an independently constructed positive-characteristic strong
Coble carrier; it is no longer presented as an output of this certificate.

Likewise, no literal same-polytope conclusion is asserted for this
calibration.  Such a conclusion is available only after adding the separate
reflection-maximal exact-return hypothesis.  The audit therefore confirms a
genuine one-packet AG existence gate and a class-relative dimension theorem,
not universal essential surjectivity or global L4.

## 5. Repairs forced before PASS

The hostile pass found three issues in the first candidate freeze.

1. The report claimed coefficientwise exceptional matching while the code
   compared only wall-index sets.  The certifier now performs and archives
   the coefficientwise comparison for every center.
2. The absolute-factor gate counted components but ignored factor
   multiplicity, so a nonreduced power could masquerade as an integral curve.
   It now requires the sole nontrivial factor to have multiplicity one and
   includes a square negative control.
3. The odd-characteristic rank line could be read as an output of a compiler
   implemented over `QQ`.  It is now separated explicitly from the proved
   `Q/C` calibration.

After those repairs, the production replay and the independent hostile suite
both pass.

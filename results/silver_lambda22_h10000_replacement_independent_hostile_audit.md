# Silver Lambda22 H10000 replacement: independent hostile audit

**Verdict: PASS.** A clean-room exact replay confirms the bounded
necessary-wall obstruction for

`z=(1,(28/87)d)`, with `d^2=4`,

in `<-(sqrt(2)-1)> + Lambda22`. The replay imports no project module. It
reimplements `Q(sqrt(2))` arithmetic, control enumeration, the binary Gram
census, lattice gates, integral positive-semidefinite extensions, and exact
ray signs.

## Complete bounded frontier

The interval is open at `2786+1970sqrt(2)` and closed at the rational field
element `10000`, with comparisons made at the identity embedding. From
height at most 10000 one gets `|x_0|<705`; conjugate target positivity gives
`|x_1|<3`. Thus, for `x=m+n sqrt(2)`, `|m|<=353` and `|n|<=250`. The replay
exhausts the larger symmetric coefficient box 354.

It independently obtains:

| quantity | exact count |
|---|---:|
| new reflective controls | 129 |
| retained gate-compatible necessary Gram types | 1,555 |
| types with nonnegative relaxed ray maximum | 1 |

The unique type is the production record at height
`3196+2260sqrt(2)`, norm 2, control `88+62sqrt(2)`, and coefficient Gram

`(u^2,v^2,u.v)=(3194,1600,2260)`.

## Independent 665-extension theorem

Put `p=d.u`, `q=d.v`, and `d^2=4`. The binary determinant is

`3194*1600-2260^2=2800`.

The two `2x2` Cauchy bounds give `-113<=p<=113` and `-80<=q<=80`.
Completing the square in the exact `3x3` determinant gives the equivalent
integer ellipse

`(3194q-2260p)^2 <= 2800(12776-p^2)`.

The hostile replay uses this ellipse to derive the complete integer `q`
interval for each `p`, independently of the owner's rectangular scan. It
finds exactly **665 integral PSD extensions**. Every one also passes the
declared local gate, which tests all nonzero `(a,b,c)` in `[-3,3]^3` for a
norm-1, -2, or -3 vector. The canonical extension stream has SHA-256

`3bc3dc04c50765b15ba1be158bdefc6d6bbbbae0e8e90b6d721a659b85556594`.

The finite gate is logically safe for this necessary-superset certificate:
if it rejects an extension, it has exhibited an actual nonzero integral
combination `au+bv+cd` of norm 1--3, impossible in `Lambda22` because its
minimum is four. The gate is not advertised as a complete ternary minimum
algorithm. Indeed, the final contradiction uses the larger coefficient 28.

## Unique nonnegative extension and contradiction

Exact multiplication gives

`(sqrt(2)-1)(88+62sqrt(2))=36+26sqrt(2)`.

Consequently

`87 <z,r> = (28p-3132)+(28q-2262)sqrt(2)`.

Exact sign checks on all 665 extensions leave exactly one nonnegative row:

`(p,q)=(113,80)`.

Its extension determinant is zero and its margin is

`(32-22sqrt(2))/87 > 0`,

because `32^2-2*22^2=56>0`. There are four determinant-zero extensions in
the full 665-set, so determinant zero alone did not select this record; the
ray-sign census did.

For the unique row,

`(v-20d)^2=1600+400*4-40*80=0`.

Positive definiteness of `Lambda22` forces `v=20d` (and
`u.(v-20d)=2260-20*113=0` is consistent). Then

`(u-28d)^2=3194-56*113+28^2*4=2`.

This vector is nonzero and belongs to `Lambda22`, contradicting the
independently certified minimum four. Hence the sole nonnegative necessary
extension is arithmetically empty.

## Continuity and scope

The audit pins the H2786 artifact, verifies the same scale and spacelike ray
norm `10705/7569-sqrt(2)>0`, and checks that the ray is strict on all four
previously accepted cutters. It separately pins the earlier independent
`Lambda22` minimum certificate before using the norm-two contradiction.

Therefore the ray is strict against every actual wall contained in the
complete necessary superset through identity height 10000. This is a
**bounded-prefix** conclusion. It is neither Vinberg termination nor a
proof of global nonreflectivity, and it says nothing about walls beyond the
cutoff.

Reproduce with:

```sh
python3 scripts/independent_hostile_audit_silver_lambda22_h10000_replacement.py --verify
PYTHONPATH=src python3 scripts/audit_silver_lambda22_h10000_replacement.py
PYTHONPATH=src pytest -q tests/test_silver_lambda22_h10000_replacement.py tests/test_silver_lambda22_h10000_replacement_independent_hostile_audit.py
```

# Actual K7 compact arithmetic dimension at most nineteen

20 October 2026. **Complete internally reviewed written proof; conditional
Lean final logic verified.** This is a geometric proof based on the named
internal local theorems, not a numerical conjecture. It has not been
externally refereed or formalized end to end.

A compact arithmetic hyperbolic Coxeter polytope whose **actual** ground
field is K7=Q(2cos(2pi/7)) has dimension at most19. The earlier universal
compact arithmetic bound20 supplies the higher-dimensional tail; the new
proof excludes dimension20. No arithmetic inheritance by arbitrary faces
or containment in a convenient larger field is assumed.

Together with the previous eight-field theorem, this proves ceiling19
for the nine actual fields Q(sqrt3), Q(sqrt6), K7, K9, K11, K15,
K16, K20 and K24, where Km=Q(2cos(2pi/m)) and
K24=Q(sqrt2,sqrt3). The unrestricted compact arithmetic ceiling stays20.
No ceiling19 is proved for Q(sqrt2), Q(sqrt5) or Q(sqrt2,sqrt5).
There is no new polytope example or sharpness claim.

## Written proof and verification

The [complete charge composition](compact_arithmetic_k7_dimension19_charge_composition_20261020.md)
has SHA256 `e6a28d7e12013be9431b4e9add5bd9b243ea8fed4da59fa82f9955d9038c9f65`.
It combines seven complete local consumer theorems into one simultaneous
nonnegative angle charge, including every shared decreasing white support
and receiving role. Every triangle pays2, every quadrilateral pays1,
and every finite rank20 vertex costs at most35. Complete forest coverage
is a deficit argument from finite-type classification, not an optimizer
support or a new enumerated census.

For total charge C, the nonnegative symmetric dual h-vector gives

    35 f0 >= C >= 5 f2 - 190 f0
       = 35 f0 + (5/2) sum_(i=0)^20 (i-10)^2 h_i
      >= 35 f0 + 500,

since h0=h20=1. This contradiction excludes dimension20. The
[root review](review_k7_dimension19_charge_composition_root_20261020.md)
and [independent review](review_k7_dimension19_charge_composition_resource_20261020.md)
read the full composition and its complete consumer interfaces.
The [nine-field proof map](../docs/compact_arithmetic_nine_field19_proof_map_20261020.md)
identifies all proof dependencies and reproduction boundaries.

The [20-page primary paper](../manuscripts/build/compact_arithmetic20.pdf)
includes both the eight-field and K7 composition appendices. Its
[source-bound build](compact_arithmetic_nine_field19_included_build_20261020.json)
passed with zero warnings and all76 input pins checked before and after
execution. PDF SHA256:
`aeb248f077bd49860705838cb8d74a6a1077c9d33380171ad52932fb32f2698b`.
The complete local geometric derivations remain an explicitly named
internal supplement; the short PDF does not replace those derivations.

`ArithmeticK7Field19.lean`, SHA256
`08f1398fd3d80b5176bb4e8f59d844ddd9a9e6f8beae089716f9d92230d4fad6`,
passed one capped compilation and a separate six-theorem axiom audit.
The kernel checks the thirteen quarter-unit family budgets, total charge,
endpoint surplus and dimension logic. Actual finite-forest coverage,
geometric payments and h-vector/face realization are explicit premises.
Only propext and Quot.sound occur in the printed audit; no sorryAx.
The older imported CompactHVector object retains its documented legacy
cache provenance.

Reproduce from the repository root:

```sh
python3 -B scripts/verify_manuscript.py compact20 --mode core --pdf
```

The current declared roots include K7, which imports and checks the cached
eight-field dependency. Its earlier source-bound audit record is retained;
a redundant standalone audit would overwrite that historical input record.
The [final combined run](compact20_nine_field_core_verification_final_20261020.json)
passes all eleven jobs: three exact scripts, seven Lean roots and the
76-input source-bound temporary PDF rebuild. The earlier run's memory-gate
deferral is preserved separately; its gate was respected. The temporary
PDF was discarded and the published20-page PDF remains unchanged.

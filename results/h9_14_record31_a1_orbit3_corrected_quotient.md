# Record 31 `a=1`, source orbit 3: corrected terminal quotient

## Outcome

The corrected production `--quotient` traversal for positive-`A4` source
orbit 3 (mask `8220`) terminated normally.  It finds exactly **11**
nonisomorphic incidence targets satisfying the necessary 120-facet and
ridge-degree-two gate.  All eleven then pass the connected mod-two
homology-sphere audit, but all eleven have zero states in the uncapped
unrestricted compact Coxeter-label join.  An independent hostile-order join
replay agrees.

Thus source orbit 3 supplies **no compact Coxeter labeling, no Gram state, and
no hyperbolic example**.  This closes this source orbit only.  It is not yet a
complete Record-31 `a=1` schedule theorem, and it says nothing about
finite-volume noncompact chambers.

## Terminal production counts

The corrected schedule is

`(14,87,319,770,1281,1477,1142,540,120)`.

| count | exact value |
|---|---:|
| search nodes | 469,990,577 |
| depth-node partition | `[1; 196; 1,554; 25,284; 691,780; 29,702,118; 439,208,906; 360,738]` |
| canonical states by depth | `[1; 9; 154; 4,388; 235,670; 14,715,798; 295,348,921; 162,031]` |
| quotient-isomorphic prunes | 159,523,605 |
| dead uncoverable | 39,299,444 |
| dead five-facet bound | 0 |
| dead slot bound | 265,737,338 |
| states with no transversal below order five | 13,273 |
| order-five count exactly 120 | 6,193 |
| pure 120-facet targets | 4,445 |
| literal targets with ridge histogram `{2:540}` | 11 |
| target orbits | 11 |

The order-five-transversal histogram is

`120:6193, 121:3552, 122:1896, 123:612, 124:547, 125:290, 126:60, 127:60, 128:35, 129:28`.

For the 6,193 states with 120 order-five transversals, the complete blocker-
size histogram is

`120:4445, 121:554, 122:79, 123:480, 124:400, 125:87, 126:79, 127:4, 128:28, 129:33, 130:3, 132:1`.

The full 4,445-entry ridge-signature histogram is retained verbatim in the
frozen JSON transcript.  Its corrected target entry is exactly
`"2:540": 11`; the histogram totals 4,445.

## Exact target families and downstream decision

Each row contains the unique `A4` mask followed by seven added `A5` masks.
Every literal multiplicity is one.  “Topology” means the exact connected
mod-two homology-sphere audit, not convex polytopality.  “Labels” is the
terminal state count of the uncapped compact Coxeter-label join.

| orbit | minimal-nonface masks | multiplicity | topology | labels |
|---:|---|---:|:---:|---:|
| 0 | `8220,121,122,1144,1924,7040,7264,9092` | 1 | pass | 0 |
| 1 | `8220,481,482,1504,1564,6680,7744,8612` | 1 | pass | 0 |
| 2 | `8220,481,482,1504,1564,6680,6848,9508` | 1 | pass | 0 |
| 3 | `8220,233,234,248,1812,6368,7936,9988` | 1 | pass | 0 |
| 4 | `8220,233,248,1812,3844,4328,6368,7936` | 1 | pass | 0 |
| 5 | `8220,481,992,1564,3100,7192,7744,8612` | 1 | pass | 0 |
| 6 | `8220,481,992,1564,3100,7192,7360,8996` | 1 | pass | 0 |
| 7 | `8220,233,248,1812,3844,6368,7936,9988` | 1 | pass | 0 |
| 8 | `8220,121,248,1924,3844,6368,7936,9988` | 1 | pass | 0 |
| 9 | `8220,481,488,1564,3604,4576,7744,8616` | 1 | pass | 0 |
| 10 | `8220,121,124,1924,3968,4216,6496,9860` | 1 | pass | 0 |

The production label joins terminate after three or four decisive support
stages.  The independent replay changes the support order, reconstructs the
classical order-four and order-five Lanner tables with 108 and 420 entries,
and again obtains terminal state count zero on every row.  Since no finite
compact Coxeter labeling exists, there is no downstream Lorentzian Gram
candidate to test.

## Integrity and scope

The raw terminal transcript is 25,556 bytes with SHA-256
`b4ea149c1259d713dd8c47910c5052d2f00c51a8b4b55c988bf22c22ec035345`.
The verifier recomputes all count identities, every literal blocker, and every
ridge histogram.  The equivalence-certified generated source has SHA-256
`186418cf82872f3145c8ebcf4e1fe953a61cb9836eb04c4734ee0e1c35fd8443`;
the production binary had SHA-256
`b740f0d8ba7fc357df6db1e50b70629533b9375a585ea2fc62f09846907a86cc`.

The frozen downstream artifacts have SHA-256
`1562231aedd9f3fc138aeb31c6af48b01cd70eee30fa32e12bab1987c9c5184b`
(production topology/label audit) and
`efd2f975f79bc3d28215e952269165b934d3744fd6e9e7a084810b93f07542af`
(independent hostile replay).

Artifacts:

- `data/enumerations/h9_14_record31_a1_orbit3_corrected_quotient.json`
- `data/enumerations/h9_14_record31_a1_orbit3_downstream.json`
- `data/enumerations/h9_14_record31_a1_orbit3_independent_downstream.json`
- `scripts/verify_h9_14_record31_a1_orbit3_transcript.py`
- `scripts/verify_h9_14_record31_a1_orbit3_downstream.py`
- `tests/test_h9_14_record31_a1_orbit3_transcript.py`
- `tests/test_h9_14_record31_a1_orbit3_downstream.py`

Replay with:

```bash
PYTHONPATH=src:. python3 scripts/verify_h9_14_record31_a1_orbit3_transcript.py
PYTHONPATH=src:. python3 scripts/verify_h9_14_record31_a1_orbit3_downstream.py
PYTHONPATH=src:. python3 scripts/explore_h9_14_record31_a1_downstream.py \
  data/enumerations/h9_14_record31_a1_orbit3_corrected_quotient.json
PYTHONPATH=src:. python3 scripts/independent_h9_14_record31_a1_downstream_audit.py \
  data/enumerations/h9_14_record31_a1_orbit3_corrected_quotient.json
```

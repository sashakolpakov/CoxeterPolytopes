# Fourteen high schedules in seven open H9/14 shadows are impossible

## Result

The top two Dehn--Sommerville schedules in each of graph-shadow records
`0,1,3,4,5,7,9` are impossible for a **compact** hyperbolic Coxeter
9-polytope.  This closes fourteen schedule branches independently of any
completion enumeration:

| record | excluded missing-four counts (a) | summed vertex-link (g_4) |
|---:|---|---|
| 0 | `17,18` | `12,6` |
| 1 | `14,15` | `12,6` |
| 3 | `14,15` | `12,6` |
| 4 | `11,12` | `12,6` |
| 5 | `9,10` | `12,6` |
| 7 | `8,9` | `12,6` |
| 9 | `8,9` | `12,6` |

No whole graph shadow is closed by this result; all lower schedules remain in
the relevant enumeration queues.  No noncompact finite-volume conclusion is
made.

## Proof

For every vertex (v), its link is the boundary of a nonsimplex simplicial
8-polytope and inherits the Coxeter--Lanner missing-face order bound five.  If
(g_4(\operatorname{lk}v)=0), the generalized lower-bound equality theorem
would give a nontrivial 3-stacked triangulation and hence a missing link face
of order at least six.  Therefore

\[
g_4(\operatorname{lk}v)\geq 1
\]

at all fourteen vertices, so their sum is at least 14.

For the boundary (K) of the polar simplicial 9-polytope, incidence
double-counting gives

\[
\sum_v g_4(\operatorname{lk}v)
=5f_4(K)-24f_3(K)+63f_2(K)-112f_1(K)+14\cdot126.
\]

Substitution of every frozen Dehn--Sommerville row gives an arithmetic
progression with step (-6).  In each of the seven shadows the last two rows
have sums 12 and 6, contradicting the required lower bound 14.

## Replay

Production uses the established link-(g_4) double-count implementation.  A
no-project-import replay reads the support archive directly, derives (h_3)
and (h_4) from the 8-dimensional link h-vector formula, and separately
checks the expanded polynomial above.  It reproduces all 95 schedule rows and
all fourteen exclusions.

- canonical certificate SHA-256:
  `a7f9d6c4ca31c91d6d456555703a273a67a70bf16ad74ab54b26afeb735156f5`
- independent summary SHA-256:
  `787e3680c788620538c108fdec894662b674315b8df2dd30d3725b64364c5b22`
- certificate file SHA-256:
  `8dc8fe7399a24a130323bbb9e0c7b2ec8865acc12bbaee948c24b88b7d38180e`

Artifacts:

- `data/enumerations/h9_14_open_high_schedule_g4_no_go.json`
- `scripts/certify_h9_14_open_high_schedule_g4.py`
- `scripts/independent_h9_14_open_high_schedule_g4.py`
- `tests/test_h9_14_open_high_schedule_g4.py`

```bash
PYTHONPATH=src:. python3 scripts/certify_h9_14_open_high_schedule_g4.py --verify
PYTHONPATH=src:. python3 -m pytest -q tests/test_h9_14_open_high_schedule_g4.py
```

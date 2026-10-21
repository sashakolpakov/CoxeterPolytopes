# Complete edge walking closes the second conductor-sixteen host

The changed host is **not reflective**. No new finite-volume polytope was found. Two independently certified first-edge steps now connect its initial corner to a corner with an integral symmetry; together with the initial F4 reversal, this gives an infinite-order stabilizer of the **full arithmetic reflection chamber**. This supersedes the open-host status in the [previous report](fresh_h9_new_host_20261013.md).

Let `K=Q(T)`, `T^4−4T²+2=0`, `T=2cos(pi/8)`, `pi=2+T`, and `alpha=(T−1)(2+T)`. The lattices are

```
M = <-alpha> + F4 + <1>             (H5),
L = <-alpha> + F4 + F4 + <1>        (H9).
```

Their determinant squareclass differs from the previously excluded coefficient `T−1` by the nonsquare `pi`. Thus the result closes a genuinely different quadratic space among the compared hosts. It is not a classification of other arithmetic spaces or Coxeter polytopes.

## Complete local edge censuses

The new search starts at the actual H5 corner `[0,1,3,4,8]` of the saved 10-wall face. This is the F4 reversal of the previously certified first endpoint. For each released wall, a positive rational parameter gives an exact timelike probe on the outgoing edge. The probe is chosen by a bound on the squared hyperbolic cosine of its distance from the starting corner; this is a finite edge interval, not a uniform Vinberg height increase.

Every primitive root norm class `1, pi, T², 2/pi, 2` is included. Folding by the finite edge reflection group fixes the entire edge. Cauchy inequalities bound each integral pairing to the fixed edge walls at all four embeddings. For each such profile, a shifted Cauchy bound controls the pairing to the released wall. An integral translation and a unit-balanced trace ball produce an exhaustive finite alphabet. The remaining quadratic is solved for **all** field square roots; no residual squareclass or nonsquare norm orbit is discarded. Integral coordinates, primitivity, reflection integrality, all embeddings, and the full finite reflection orbit are checked.

| Released wall | Probe cosh² bound | Profiles | Quadratics | Square cases | Primitive signed solutions | Result |
|---|---:|---:|---:|---:|---:|---|
| 0 | 4 | 19 | 32 | 6 | 6 | first mirror; 19-mirror H5 endpoint |
| 1 | 4 | 19 | 32 | 6 | 6 | no new mirror in this interval |
| 3 | 4 | 24 | 32 | 16 | 16 | no new mirror in this interval |
| 4 | 4 | 13 | 25 | 6 | 2 | no new mirror in this interval |
| 1 | 64 | 19 | 138 | 38 | 38 | first mirror; 33-mirror H5 endpoint |
| 3 | 64 | 24 | 38 | 16 | 16 | no new mirror in this interval |
| 4 | 64 | 13 | 53 | 6 | 2 | no new mirror in this interval |

Both successful edges reach the same new mirror, with norm `2/pi` and vector

```
(-40−21T+68T²+36T³,
 -82−44T+140T²+76T³,
 -165−88T+280T²+152T³,
 -172−96T+296T²+162T³,
 -86−48T+148T²+81T³,
 0).
```

All seven censuses pass the [separate Fraction replay](fresh_h9_new_host_next_20261013_all_census_independent.json). That implementation chooses trace-balancing units independently and obtains the quadratic anchor by inversion of the entire positive rank-five Gram, rather than the producer's residual-axis formula. It also replays the previously complete second-F4 census. For every allowed norm, that block has no intermediate energy `0<S<q` at all embeddings, so the full H9 first-edge statements follow. The 33-mirror endpoint has 57 incident mirrors in H9, including the other F4.

## The exact global obstruction

Use the released-wall-1 edge with support `[0,3,4,8]`. Its endpoint has spherical type `I2(16)²+A1` in H5. In the incident order `[w0,w3,w4,w8,new]`, the permutation `[1,0,2,4,3]` exchanges the two I2(16) components with exactly matching root norms and preserves the isolated A1. The corresponding form isometry `H`, fixing the timelike corner line, is constructed by exact incident-Gram inversion. Both `H` and its inverse are integral, and `H²=1`.

Let `J` reverse the initial F4 chain and fix time and the scalar coordinate. The [global certificate](fresh_h9_new_host_next_20261013_global.json) checks:

- `J²=H²=1`, exact form preservation, integral inverses, and preservation of the future component;
- the exact transport by `J` of the previously certified first endpoint and every oriented incident simple root to the start of the new edge;
- the complete initial, intermediate, and final H9 stabilizers, with 49, 44, and 57 mirrors respectively;
- the positive oriented permutation of the final simple roots by `H`, including an explicit fixed inward tangent;
- extension of both maps by the identity on the second F4;
- the exact product trace
  `tr(HJ) = −56022−30320T+95640T²+51760T³ > 6` on H5, and `tr(HJ)+4 > 10` on H9.

The complete nearest-mirror path identifies the same full reflection chamber at all three corners. A lattice isometry permuting the oriented simple roots of the complete finite stabilizer preserves the local chamber sector and hence this chamber. Thus both `J` and `H` stabilize it. Their product has infinite order by the trace bound. Discreteness of the arithmetic lattice isometry group then gives infinitely many disjoint translates of a small interior ball, proving infinite volume. No finite-facet hypothesis is used.

This is a full-host nonreflectivity theorem, rather than an escape from a finite prefix. It also excludes a finite-volume reflection subgroup generated by mirrors of this lattice: such a subgroup chamber contains a chamber of the full reflection group. It says nothing against changing the lattice or using unrelated mirrors.

The frontier agent independently source-reviewed and reran the global certificate, including the oriented-sector joins and H9 extension: **PASS**. This review is separate from the independent first-edge enumerations; it does not claim a third enumeration.

## Reproduction

```sh
PYTHONPATH=src:. python scripts/fresh_h9_new_host_next_20261013_global.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_new_host_next_20261013_independent.py \
  results/fresh_h9_new_host_next_20261013_census_e0_c4.json \
  results/fresh_h9_new_host_next_20261013_census_e1_c4.json \
  results/fresh_h9_new_host_next_20261013_census_e2_c4.json \
  results/fresh_h9_new_host_next_20261013_census_e3_c4.json \
  results/fresh_h9_new_host_next_20261013_census_e1_c64.json \
  results/fresh_h9_new_host_next_20261013_census_e2_c64.json \
  results/fresh_h9_new_host_next_20261013_census_e3_c64.json \
  --output results/fresh_h9_new_host_next_20261013_all_census_independent.json --verify
```

Producer: [shifted interval census](../scripts/fresh_h9_new_host_next_20261013_census.py). Independent checker: [Fraction census](../scripts/fresh_h9_new_host_next_20261013_independent.py). Global proof: [exact integral involutions](../scripts/fresh_h9_new_host_next_20261013_global.py). All seven jobs are complete; no search remains active on this excluded host.

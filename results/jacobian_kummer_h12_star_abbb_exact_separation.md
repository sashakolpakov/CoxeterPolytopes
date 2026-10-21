# TOOT: exact star--`ABBB` separation in the Jacobian--Kummer `H12` factory

Status: **PASS/FREEZE; complete literal-byte production replay and independent
archive replay green** (2026-08-08).

## Exact project-scoped theorem

Inside the frozen `W61,...,W172` 112-wall Jacobian--Kummer stock, the complete
five-vertex star layer has 16,872 pairwise nonisomorphic ordinary Coxeter
matrices and the complete `ABBB` path shard has 36,247.  Their only common
facet header is 779:

```text
K1,4 with 779 facets:   6,117 types
ABBB P5 with 779 facets: 36,247 types
```

The new certificate rebuilds all 42,364 boundaries in that common header and
compares their exact ordinary degree-profile byte strings.  The intersection
is empty.  The other 10,755 star types have 752, 761, 770, or 788 facets and
are therefore automatically distinct from every 779-facet `ABBB` type.
Consequently the two frozen layers together contain exactly

```text
16,872 + 36,247 = 53,119
```

pairwise nonisomorphic finite-volume noncompact Coxeter polytopes in `H12`,
all of reflection index and volume multiplier five.

## The literal-byte invariant

For each facet `v` of an ordinary Coxeter matrix, form the five-vector

```text
(deg_m2(v), deg_m3(v), deg_m4(v), deg_m6(v), deg_infinity(v)).
```

Sort the 779 vectors lexicographically and serialize their unsigned entries as
little-endian 16-bit integers.  The result is exactly 7,790 bytes.  An
ordinary Coxeter-matrix isomorphism only permutes facets, so it leaves this
sorted byte string unchanged.  Equality is not a complete isomorphism test,
but unequal strings are an exact obstruction to isomorphism.

The star side has 5,949 distinct byte keys among its 6,117 rows.  The
certifier retains those actual `bytes` objects and tests every rebuilt
`ABBB` key by literal set membership.  It finds

```text
cross_shape_exact_byte_match_count = 0.
```

No hash equality or collision-resistance assumption enters this comparison.
SHA-256 has only three archival roles:

1. pinning the immutable source builders and archives;
2. checking that each newly rebuilt byte string is exactly the row represented
   by its frozen source archive; and
3. recording compact digests of the complete archive-order byte streams.

Those stream digests are

```text
star-779: 98496abf4ea3c653eb98730be30ebf1d3af7afed44484534d65902ddac364dc1
ABBB:     9c4811d7fafb06b343982c72e99836346b2ee9e8a7bd9218cfca1477feb2381e
```

The source star and `ABBB` theorems still supply their respective internal
pairwise-nonisomorphism proofs; the present result closes only the missing
cross-shape comparison.

## Hyperbolicity and carrier boundary

Finite volume is not inferred from these byte strings or from the divisor
lattice.  It remains transported from the independently certified nonuniform
172-facet seed by Allcock's finite-subtree theorem.  Each row is a five-copy
fundamental polytope, hence has reflection index and volume multiplier five;
a finite-index subgroup of a nonuniform lattice remains nonuniform.

All 53,119 rows retain the same weak full-Picard K3 divisor-reflection carrier

```text
U(2) + D5(-1)^2 + <-4>.
```

This is a `SURFACE -> POLYTOPE` carrier-preserving redoubling result.  It is
not a nef-cone or effectivity theorem, not a surface-automorphism theorem, and
not a `POLYTOPE -> SURFACE` recognizer.

## Novelty and completion boundary

This is an exact promotion of **project cross-shape accounting**.  It does not
make 779 a new facet header, does not prove that these literal polytopes were
absent from all prior literature, and makes no worldwide novelty claim.
Worldwide priority remains open.

Nor does it complete the full `P5` or all-shape five-vertex census.  The
`BABB` and `BBBB` path shards still have a combined 144,807 color-orbit
denominator awaiting matrix census, and the degree-three fork remains open.

## Reproduction and frozen artifacts

```bash
python3 scripts/certify_jacobian_kummer_h12_star_abbb_exact_separation.py --check-archive
PYTHONPATH=src:. pytest -q \
  tests/test_jacobian_kummer_h12_star_abbb_exact_separation.py
```

Frozen SHA-256 values:

```text
certifier: f27f18c36531c48f75e93f5f72b36da694ef0de6f64175d7e84b38186991b926
archive:   7a328edb6cbedafa13c32140c2f2f1262a8b1159617a1fc59d5ad8d2165c3fe6
tests:     33ab916e7338483f8e1ecb55e1e9f8942d38453e6fc5532e5d1662e6d2ea3e5a
```

Artifacts:

- [`certify_jacobian_kummer_h12_star_abbb_exact_separation.py`](../scripts/certify_jacobian_kummer_h12_star_abbb_exact_separation.py)
- [`jacobian_kummer_h12_star_abbb_exact_separation.json`](../data/enumerations/jacobian_kummer_h12_star_abbb_exact_separation.json)
- [`test_jacobian_kummer_h12_star_abbb_exact_separation.py`](../tests/test_jacobian_kummer_h12_star_abbb_exact_separation.py)

# Hostile audit: exact four-black `WW` no-go

**Status (2026-08-10): PASS.**

This audit covers the theorem that no degree-`29` four-black system can put
both normalized endpoints in the white fibre.  It includes all endpoint
statuses `SS,SC,CS,CC`, but it does not close the remaining `BW+WB`
one-black/white-simple block or the global degree-`29` problem.

## Pinned artifacts

```text
6528cf01fcd16f90ef65389226dfd1e094794d4f48490e97b6eaa86b8beeb8b1  scripts/certify_gram_belyi_m5_degree29_four_black_ww_simple.py
b6d3d547a9862144867ce1c871ea101c83a399d2f2c5907b7a830650f9ff6ecf  data/enumerations/gram_belyi_m5_degree29_four_black_ww_simple.json
7d83fc06d564b52ec21739c310efcdace2af6c9222021dfd7333734ca9a81d10  results/gram_belyi_m5_degree29_four_black_ww_simple.md
142e30eb15898b0c231c5085f64cf2268df355b8d714f7e8f86556a3772119e3  tests/test_gram_belyi_m5_degree29_four_black_ww_simple.py
```

The canonical replayed row-list digest is

```text
679979b1d3544e86b8c2ab4eaa34029435d2f16aaa2f2dd318f0abf16fe1ad6
```

## Independent checks

The hostile audit independently verified:

1. The affine normalization sending the marked root to zero and another
   distinct black root to one, including the endpoint ratio
   `q=-(1-lambda)/lambda=29 lambda-13`.
2. The scheme-theoretic equivalence between a constant remainder modulo the
   residual cubic and `R_1=R_2=0`, including repeated residual roots.
3. Total degrees `28,27`, bidegrees `(27,27)`, and the sharp degree-`756`
   resultant bound and preservation.
4. Coverage of every irreducible `u` factor and every factor of the
   specialized `v` gcd.  Weighted `u`-factor degrees sum to `756` for every
   row.
5. The local-integrality chain for `u,v,y`; no characteristic-zero branch can
   disappear at infinity after reduction.
6. Tower-field arithmetic by static review and a separate pure-Python replay
   over `F_101[T]/(T^2+61)`, which independently returned endpoint gcd one.
7. The strengthened scope.  No simplicity inequation is used, so the same
   `640` broad systems exclude all `2,560` status systems.  The marked-record
   multiplier is `3+4+4+3=14`, giving `8,960` records.

## Replay and resources

The independent full replay matched all `640` rows and all `10,696` endpoint
orbits exactly.  It completed in `441` seconds with peak RSS `88.6 MB`, below
the `0.5 GiB` cap.  Targeted tests gave

```text
6 passed in 0.95s
```

The archive stores every monic degree-`756` resultant and every component and
orbit degree.  It does not serialize all factor polynomials or Bezout
witnesses; independent validation of those objects relies on the deterministic
exact replay, which passed.

## Claim boundary

```text
four-black WW black-degree master systems eliminated = 640
four-black WW endpoint-status systems eliminated = 2560
four-black WW marked passports eliminated = 8960
four-black systems remaining = 3458
four-black marked passports remaining = 10374
global marked-passport remainder = 42572709
four-black status = PARTIAL
global degree-29 status = OPEN
```

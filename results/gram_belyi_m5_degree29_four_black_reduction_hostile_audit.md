# Hostile audit: degree-29 four-black reduction

**Status (2026-08-10): PASS.**

This audit covers the exact elimination of `6,734` four-black master systems,
representing `23,660` marked passports.  It does not claim that the remaining
four-black systems, or degree `29` globally, are impossible.

## Pinned artifacts

```text
c816a4dbc66df9b826eb29f28ff1d49a0f1265b65811b8053159bf9f23c2685f  scripts/certify_gram_belyi_m5_degree29_four_black_reduction.py
abe0e36528737331eafbcb3eef5f3b9b088ba3fff9ca0d0038dbb0a50189f641  results/gram_belyi_m5_degree29_four_black_reduction.md
279622e36dd47756e76e342decb7da2583f63a9a33e6e9df37ca2e9265066cb1  data/enumerations/gram_belyi_m5_degree29_four_black_reduction.json
7a1e5e93a43573b8c512b45a56f1b163c0e38a86f9894602c9c694488344088a  tests/test_gram_belyi_m5_degree29_four_black_reduction.py
```

## Independent algebra audit

The hostile audit independently checked the following points.

1. If `D` is the residual derivative cubic and
   `F mod D=R_2 x^2+R_1 x+R_0`, equality of all non-root critical values is
   exactly `R_1=R_2=0`, including every repeated-root factorization of `D`.
2. In the two-black-endpoint block, homogeneity gives preliminary parameter
   degrees `28,27`.  The pure degree-`28` term of `R_1` vanishes when the
   other three black roots coalesce at zero: both `F` and `D` are divisible
   by `x^2`, so the remainder has only an `x^2` term.  The certified degree
   bounds `27,27` are therefore correct.
3. In the one-black/critical-white block, `D(1)=C(r)s+K(r)` with affine
   `C,K`.  The stored nonzero linear resultant excludes the complete
   exceptional branch `C=K=0`; division by `C` loses no solution.
4. Substitution `s=-K/C` gives the cleared equations `C^28 R_1` and
   `C^27 R_2`.  Their bounds are `55,54`: the only nominal degree-`56` term
   is again the vanishing pure `r^28` coefficient of `R_1`.
5. Newton interpolation is exact.  Fifty-six distinct values determine a
   polynomial of degree at most `55`, and the producer verifies the result at
   every remaining residue value with `C(r)!=0`.  Every used prime exceeds
   `56`.
6. The reflection `x -> 1-x` correctly sends the `WB` case to a black
   endpoint at zero, a critical endpoint at one, and marked root `1-lambda`.
7. All primes split in `Q(sqrt(5))` and avoid `29`.  Preserved exact degrees,
   nonzero exceptional resultants, and full finite-field gcd one therefore
   prove nonzero characteristic-zero resultants.

The separate quotient audit also recovered all census numbers and emphasized
the correct scope: `12,752` counts necessary simple/critical master systems,
not multiplicity-specific constructible passport schemes.

## Replay

The independent replay checked all

```text
3276 + 3458 = 6734
```

rows, compared their keys against the expected system set, and passed.  It
returned archive SHA-256

```text
279622e36dd47756e76e342decb7da2583f63a9a33e6e9df37ca2e9265066cb1
```

Targeted tests gave `6 passed`; the combined Gram--Belyi suite gave
`25 passed`.

## Claim boundary

```text
eliminated four-black systems = 6734
eliminated four-black marked passports = 23660
remaining four-black systems = 6018
remaining four-black marked passports = 19334
global marked-passport remainder = 42581669
four-black status = PARTIAL
global degree-29 status = OPEN
```

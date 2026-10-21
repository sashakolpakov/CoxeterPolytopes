# ADEG24 remaining-package K3 stock promotion ledger

This is the frozen row ledger accompanying the
[promotion theorem](bredon_adeg24_remaining_k3_stock_promotions.md) and
[exact archive](../data/enumerations/bredon_adeg24_remaining_k3_stock_promotions.json).
It contains all 30 stable-range packages left after the already closed `N11`
row and separately audited `P_star` are removed.

Status codes:

- `P`: exact complete-nef divisor-lattice fixed-section promotion;
- `R`: every locally frozen stock with the same rank/root-count header is
  rejected by Smith form or exact weighted Gram;
- `O`: open relative to this search because the 55-stock local denominator
  has no matching header.

`O` is not a nonexistence result.  The denominator consists of 17
target-header-selected complete Roulleau stocks in ranks 5--8 and 38 locally
frozen complete K3 nef stocks in ranks 9--19; it is not an all-rank or
worldwide stock classification.

| Bredon row | orientation | unfolded rank/roots | matching local headers | status |
|---|---:|---:|---|---:|
| `ADEG-N2-01-path-66` | 0 | `7/7` | `R7-UA5-NEF-H6-7`, `R7-UD5-NEF-H6-7` | `R` |
| `ADEG-N2-01-path-66` | 1 | `5/5` | `R5-UA3-NEF-H4-5`, **`R5-P4D4-NEF-H4-5`** | **`P`** |
| `ADEG-N2-03-triangle-663` | 0 | `7/7` | `R7-UA5-NEF-H6-7`, `R7-UD5-NEF-H6-7` | `R` |
| `ADEG-N2-03-triangle-663` | 1 | `5/5` | `R5-UA3-NEF-H4-5`, `R5-P4D4-NEF-H4-5` | `R` |
| `ADEG-N3-01-path-633` | 0 | `6/6` | `R6-UA4-NEF-H5-6`, `R6-U2D4-NEF-H5-6`, **`R6-UD4-NEF-H5-6`** | **`P`** |
| `ADEG-N3-01-path-633` | 1 | `10/10` | `R10-UE8-NEF-H9-10`, `R10-UD8-NEF-H9-10` | `R` |
| `ADEG-N3-02-path-636` | 0 | `8/8` | `R8-UD6-NEF-H7-8`, `R8-UA6-NEF-H7-8`, `R8-UE6-NEF-H7-8` | `R` |
| `ADEG-N3-02-path-636` | 1 | `8/8` | `R8-UD6-NEF-H7-8`, `R8-UA6-NEF-H7-8`, `R8-UE6-NEF-H7-8` | `R` |
| `ADEG-N3-03-path-363` | 0 | `8/8` | `R8-UD6-NEF-H7-8`, `R8-UA6-NEF-H7-8`, `R8-UE6-NEF-H7-8` | `R` |
| `ADEG-N3-03-path-363` | 1 | `8/8` | `R8-UD6-NEF-H7-8`, `R8-UA6-NEF-H7-8`, `R8-UE6-NEF-H7-8` | `R` |
| `ADEG-N3-04-star-6` | 0 | `6/6` | `R6-UA4-NEF-H5-6`, **`R6-U2D4-NEF-H5-6`**, `R6-UD4-NEF-H5-6` | **`P`** |
| `ADEG-N3-04-star-6` | 1 | `10/10` | `R10-UE8-NEF-H9-10`, `R10-UD8-NEF-H9-10` | `R` |
| `ADEG-N3-05-triangle-leaf-6` | 0 | `6/6` | `R6-UA4-NEF-H5-6`, `R6-U2D4-NEF-H5-6`, `R6-UD4-NEF-H5-6` | `R` |
| `ADEG-N3-05-triangle-leaf-6` | 1 | `8/10` | `R8-UD4A1x2-NEF-H7-10`, **`R8-UA2x3-NEF-H7-10`** | **`P`** |
| `ADEG-N3-07-square-opposite-6` | 0 | `8/8` | `R8-UD6-NEF-H7-8`, `R8-UA6-NEF-H7-8`, `R8-UE6-NEF-H7-8` | `R` |
| `ADEG-N3-07-square-opposite-6` | 1 | `8/8` | `R8-UD6-NEF-H7-8`, `R8-UA6-NEF-H7-8`, `R8-UE6-NEF-H7-8` | `R` |
| `ADEG-N5-01-path-633336` | 0 | `10/11` | `R10-UE7A1-NEF-H9-11`, **`R10-UD4D4-NEF-H9-11`**, `R10-UA2E6-NEF-H9-11` | **`P`** |
| `ADEG-N5-01-path-633336` | 1 | `16/17` | none | `O` |
| `ADEG-N5-02-path-triangle` | 0 | `8/9` | `R8-UA3x2-NEF-H7-9`, `R8-UA2A4-NEF-H7-9`, `R8-UA1A5-NEF-H7-9`, **`R8-UA2D4-NEF-H7-9`**, `R8-UA1D5-NEF-H7-9` | **`P`** |
| `ADEG-N5-03-octagon-alternating` | 0 | `14/16` | none | `O` |
| `ADEG-N5-03-octagon-alternating` | 1 | `14/16` | none | `O` |
| `ADEG-N5-04-octagon-clustered` | 0 | `14/16` | none | `O` |
| `ADEG-N5-04-octagon-clustered` | 1 | `14/16` | none | `O` |
| `ADEG-N6-01-diamond-tail` | 1 | `9/10` | `R9-UD6A1-NEF-H8-10`, **`R9-UD4A3-NEF-H8-10`**, `R9-UD5A2-NEF-H8-10`, `R9-UA1E6-NEF-H8-10` | **`P`** |
| `ADEG-N7-01-pentagon-tail` | 1 | `10/11` | `R10-UE7A1-NEF-H9-11`, `R10-UD4D4-NEF-H9-11`, `R10-UA2E6-NEF-H9-11` | `R` |
| `ADEG-N7-02-path-star` | 0 | `10/11` | `R10-UE7A1-NEF-H9-11`, **`R10-UD4D4-NEF-H9-11`**, `R10-UA2E6-NEF-H9-11` | **`P`** |
| `ADEG-N7-03-three-armed-path` | 0 | `14/16` | none | `O` |
| `ADEG-N7-04-triangle-branch-path` | 0 | `12/14` | `R12_UE8A1x2` | `R` |
| `ADEG-N7-05-two-triangles-branch` | 1 | `10/12` | none | `O` |
| `ADEG-N9-01-branched-path` | 0 | `12/13` | `R12_UA2E8` | `R` |

The census is therefore

```text
P = 8,    R = 15,    O = 7.
```

For every `P` row, the exact archive freezes the complete wall permutation,
unimodular basis change, triality matrix, invariant ample class, full-nef
provenance, and orbit-sum same-polytope certificate.  The order-three K3
automorphism/Hodge lift remains open in all eight rows.

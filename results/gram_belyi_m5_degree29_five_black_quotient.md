# Exact degree-29 five-black passport/master quotient

**Status (2026-08-10): the five-black marked-passport census of `190,047`
records maps to exactly `33,079` necessary binary-status coefficient
systems.  This is a combinatorial quotient, not a five-black no-go.**

The convention is `B(lambda)=0`.  The letters `B/W` retain the ordered
colours of the normalized source points `x=0,x=1`.  Cycles with equal
valency are indistinguishable, while endpoint order is retained.  As in the
four-black reduction, `S/C` records only whether a white endpoint is simple
or critical; a black endpoint's exact degree is already part of its black
degree record.

## 1. White partitions and endpoint statuses

For a degree-`29` polynomial, the black and white vertex counts sum to `30`.
Five black vertices therefore force `25` white vertices, so the white
ramification excess is `29-25=4`.  There are exactly five white partitions:

| white partition | `S` | `C` | `SS` | `SC` | `CS` | `CC` |
|---|---:|---:|---:|---:|---:|---:|
| `(5,1^24)` | 1 | 1 | 1 | 1 | 1 | 0 |
| `(4,2,1^23)` | 1 | 2 | 1 | 2 | 2 | 2 |
| `(3^2,1^23)` | 1 | 1 | 1 | 1 | 1 | 1 |
| `(3,2^2,1^22)` | 1 | 2 | 1 | 2 | 2 | 3 |
| `(2^4,1^21)` | 1 | 1 | 1 | 1 | 1 | 1 |
| **total** | **5** | **7** | **5** | **7** | **7** | **7** |

Thus one marked white endpoint has `5+7=12` exact valency choices across
the five partitions.  Two ordered white endpoints have
`5+7+7+7=26` choices.  For example, the `CC=3` entry for
`(3,2^2,1^22)` consists of the ordered pairs `(3,2)`, `(2,3)`, and `(2,2)`.

## 2. Black degree-coordinate records

Three black-side counts are needed:

| symbol | retained data | canonical degree form | count |
|---|---|---|---:|
| `M0` | marked `lambda` root | `(a,b,c,d,e)`, `b<=c<=d<=e` | 1,292 |
| `M1` | marked root and one black endpoint | `(a,b,c,d,e)`, `c<=d<=e` | 4,225 |
| `M2` | marked root and two ordered black endpoints | `(a,b,c,d,e)`, `d<=e` | 11,011 |

All entries are positive and sum to `29`.  The producer obtains these counts
in two independent ways: first by marking parts in every length-five
partition of `29`, and second by directly enumerating the displayed
ordered/unordered degree records.  Both give `(1292,4225,11011)`.

The canonical record-list SHA-256 digests are:

```text
M0 = fc18d3c8fc308091820ef68d9ff6642cff96d4660eed39ba684ba6afffeccea0
M1 = 37dc0743a6ba1a3132f54f6a08168ad1983f69c4f92a2134cffaaaa8f54e206a
M2 = 0b4978f7cbdf4c4c6b883e689e22d73950f9b45818ef3f3e1c3e835ac4596d2a
```

## 3. Raw marked-passport census

Restoring exact white partitions and endpoint valencies gives:

| endpoint colours | factorization | raw marked records |
|---|---:|---:|
| `BB` | `M2*5 = 11011*5` | 55,055 |
| `BW` | `M1*12 = 4225*12` | 50,700 |
| `WB` | `M1*12 = 4225*12` | 50,700 |
| `WW` | `M0*26 = 1292*26` | 33,592 |
| **total** |  | **190,047** |

This independently reproduces the five-black row of the global marked
passport census.  These are Riemann--Hurwitz-compatible marked records, not
a claim that `190,047` dessins exist.

## 4. Necessary binary-status master systems

For no-go calculations, exact white ramification indices can be forgotten
while the `S/C` status is retained.  The exact quotient, including its count
multiplier back to raw marked passports, is:

| refined block | master systems | marked multiplier | raw records |
|---|---:|---:|---:|
| `BB` | 11,011 | 5 | 55,055 |
| `BW-S` | 4,225 | 5 | 21,125 |
| `BW-C` | 4,225 | 7 | 29,575 |
| `WB-S` | 4,225 | 5 | 21,125 |
| `WB-C` | 4,225 | 7 | 29,575 |
| `WW-SS` | 1,292 | 5 | 6,460 |
| `WW-SC` | 1,292 | 7 | 9,044 |
| `WW-CS` | 1,292 | 7 | 9,044 |
| `WW-CC` | 1,292 | 7 | 9,044 |
| **total** | **33,079** |  | **190,047** |

Aggregating only by endpoint colours gives:

| endpoint colours | binary-status master systems |
|---|---:|
| `BB` | 11,011 |
| `BW` | 8,450 |
| `WB` | 8,450 |
| `WW` | 5,168 |
| **total** | **33,079** |

This is a necessary-system quotient only.  It does not assert that every
binary-status system realizes each represented exact passport: exact white
indices impose additional multiplicity and inequation conditions on the
residual derivative quartic.

## 5. Reproducibility and claim boundary

The standalone producer and canonical archive are:

- [`certify_gram_belyi_m5_degree29_five_black_quotient.py`](../scripts/certify_gram_belyi_m5_degree29_five_black_quotient.py), SHA-256
  `deaaed9f90e550ecd8a2256dc65a2bdc4ca57aea7e18354b7dabfff0d344936f`;
- [`gram_belyi_m5_degree29_five_black_quotient.json`](../data/enumerations/gram_belyi_m5_degree29_five_black_quotient.json), SHA-256
  `89e2d1fafc10ac5ad91edd6e20ed5e96f6c238b404b8bb303a054f330cbc0e87`.

Replay with:

```text
python3 scripts/certify_gram_belyi_m5_degree29_five_black_quotient.py --check
```

The calculation uses only exact integer partitions and counts.  It launches
no coefficient resultant, eliminates no five-black system, and leaves the
degree-`29` minimum open.

```text
five_black_raw_marked_passports = 190047
five_black_binary_status_master_systems = 33079
coefficient_resultants_run = false
five_black_no_go_claimed = false
global_degree_29_status = OPEN
```

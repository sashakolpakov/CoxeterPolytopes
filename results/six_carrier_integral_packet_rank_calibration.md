# Six-carrier faithful-packet calibration

**Status (2026-08-10): exact lightweight experiment on six frozen carriers.**
This report reads only the direct SHA-pinned numerical inputs listed in the
machine ledger.
It performs no chamber census and proves no bound for general maximal
arithmetic chambers.

## Result

The trace module is not the missing object in any of these six controls.
Every row has a genuine global injection into divisor cohomology.  In the five
rational-field full-nef rows it is simply the whole rational Neron--Severi
space; for Bugaenko it is a rank-18 integral sublattice of a rank-20 rational
surface Picard lattice and a rank-18 divisor-Tate Chow summand.

What is missing is a **uniform** AG ceiling.  The rowwise values below are
certificates for six fixed inputs, not constants forced by reflection
maximality.

Here `B,C,D,t` have exactly the meanings in the bounded faithful-packet
theorem: packet rank, integral-coordinate height, total ample facet degree,
and number of interval generators.  A dash in `D` means “not certified,” not
“does not exist.”

| carrier | `n/f` | `k`, `d` | carrier `rho=b2` | trace rank `d(n+1)` / `B` | actual cohomological return | `C` | certified `D` | `t` | missing instance gate |
|---|---:|---|---:|---:|---|---:|---:|---:|---|
| Bugaenko rational eigencarrier | `8/11` | `Q(sqrt(5))`, 2 | 20 | `18/18` | rank-18 lattice embeds in `NS=I_(1,19)`; Chow Tate summand | 48 | -- | 1 | `D` |
| toric full-nef endpoint | `6/9` | `Q`, 1 | 7 | `7/7` | identity with full `NS`; integral wall span | 1 | 12 | 0 | none |
| Coble Type II | `10/30` | `Q`, 1 | 11 | `11/11` | identity with full `NS_Q`; integral wall span | 17 | -- | 0 | `D` |
| Coble Type I | `11/26` | `Q`, 1 | 12 | `12/12` | identity with full `NS_Q`; physical rank 12 | 21 | -- | 0 | `D` |
| Coble Type VI, char. 3 | `14/60` | `Q`, 1 | 15 | `15/15` | identity with `NS_Q=I_(1,14)` | 8 | -- | 0 | `D` |
| Coble Type MII, char. 3 | `17/48` | `Q`, 1 | 18 | `18/18` | identity with `Pic_Q=I_(1,17)`; index-one wall span | 1 | 168 | 0 | none |

`rho=b2` is literal for these rational surfaces.  In the Bugaenko row it is
an invariant of the chosen weak rational eigencarrier, not an invariant of
the abstract chamber or a claim of rank-minimality of the ambient surface.

## 1. The two exact degree packets

For the toric endpoint, the archived integral interior class

\[
 A=7H-2(E_1+E_2+E_3)-(E_4+E_5+E_6)
\]

has degree pattern `(1,1,2)^3` on the nine invariant facet curves.  The
complete-nef theorem makes this strict integral nef interior class ample, so

\[
                         \sum_{i=1}^9 A\cdot D_i=12.
\]

For Type MII, the archived class

\[
 H=5F_1+5F_2-\sum_{i,j}E_{ij}
\]

lies strictly inside the certified full nef cone.  It has degree one on the
eight boundary and sixteen exceptional facets, and degree six on the twenty-
four graph facets.  Hence

\[
                    \sum_{i=1}^{48}H\cdot D_i
                    =8+16+24\cdot6=168.
\]

No analogous number is promoted for the other rows.  Bugaenko's eleven
transported divisor classes are not asserted effective or irreducible.  For
the two plane Coble models, raw plane degree is not ample after blowing up.
The Type-II, Type-I, and Type-VI archives do not pin an integral ample class
together with the total degree of every facet curve.

## 2. The exact field and coordinate packets

The five rational-field rows need no spectral generator: `t=0`.  For the
golden row use

\[
 \alpha=2\phi=1+\sqrt5,\qquad
 \alpha^2-2\alpha-4=0.
\]

Since `2<sqrt(5)<3`, the distinguished value lies in `[2,14]` and the other
conjugate `1-sqrt(5)` lies in `[-2,1]`.  Multiplication by `alpha` is integral
on the frozen rank-18 trace lattice.  Thus this row really passes the
project-native interval gate with `t=1`; one need not use the ambient
integral extension whose eigenvalue is `4phi`.

The displayed `C` values are exact ceilings in the stated pinned bases, not
minimal heights:

```text
Bugaenko trace basis: 48   toric NS basis: 1
Type-II NS basis:     17   Type-I NS basis: 21
Type-VI NS basis:      8   Type-MII NS basis: 1.
```

They bound the packet pairing/action/facet coordinates for each fixed row.
Taking their sample maximum `48` is not a class theorem.

## 3. No monodromy number is hidden in the table

No row has a certified faithful monodromy representation whose rank can be
entered as a replacement for `B`.  The Bugaenko carrier has a rank-18 Chow
correspondence action, with an integral ambient extension of 40 divisor-
product terms and height 9; this is not monodromy.  In the toric row,
`W(A2)^3` of order 216 is the deck group of the product simultaneous-
resolution marking cover.  The frozen theorem explicitly declines to call
it the monodromy image of a single smoothing pencil.  The four Coble rows
archive no monodromy packet.

## 4. What the experiment decides

Instancewise, `B`, `C`, and `t` pass on every row; `D` passes only on the
toric `H6` and Type-MII `H17` rows.  This is already enough to reject the
idea that the known surface examples fail merely because their trace modules
remain external.

It does **not** establish any of the four uniform gates for general maximal
arithmetic finite-volume chambers.  The descriptive sample maxima

```text
B=18, C=48, t=1, carrier rho=20
```

are not global bounds, and `D` is missing on four controls.  The unresolved
theorem remains exactly the maximality-to-bounded-packet implication: one
must force `B,D,C,t` on an entire declared class, not read them after the
chamber and carrier have already been fixed.

## 5. Machine ledger

The compact ledger is
[`six_carrier_integral_packet_rank_calibration.json`](../data/enumerations/six_carrier_integral_packet_rank_calibration.json).
Its focused test hashes every declared direct numerical source, rechecks all rank identities,
recomputes the two total degrees and all six `C` ceilings, and verifies that
no monodromy rank was manufactured:
[`test_six_carrier_integral_packet_rank_calibration.py`](../tests/test_six_carrier_integral_packet_rank_calibration.py).

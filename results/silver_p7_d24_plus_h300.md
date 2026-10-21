# Silver prime-seven `D24+`: exact H=300 wall census

## Outcome

For

\[
L=\langle-(2\sqrt2-1)\rangle\oplus D_{24}^{+},
\]

the glue-complete Vinberg enumeration has now been extended from algebraic
height 100 to height 300.  It gives an exact census of 324 accepted
positive-height walls, hence 348 walls after adjoining the 24 height-zero
`D24` walls.

This complete prefix is **not of finite volume**.  The old H=100
hyperideal witness survives the first new wall but is cut by the next
packet, and all 109 hyperideal rays of the first complete packet are cut at
that stage.  Nevertheless a separate exact calculation produces a new
rank-24 hyperideal extreme ray satisfying all 348 wall inequalities.  This
settles the H=300 prefix only: later mirrors may cut that ray, so the eventual
chamber remains unresolved.

## H=100 lineage

The independently reverified H=100 freeze was pinned before extending the
backend and tests:

| file | SHA-256 |
|---|---|
| backend | `e9d9359209e980eb7e1d76371c3de344e2ac8d1a9276ee5bb84bcb73f367cd81` |
| tests | `719e37b48b2edd2a55a99ff21e3f28fd5efa0af2f581c0cc12caf0a8c8650533` |
| audit | `67997a7ba00214b573b5574ebc74b3cd6ff3a40d7c262e8707d84c1fe36da2c8` |
| artifact | `ec01a747418b4fb0ecbfa25040cbae6d3f6369746f533c5b5447ff9ac11c5fbd` |
| report | `5d88f4454ff5a6d23c51d89ddb04034b6c8be760b0bea621df06d6131740afea` |

That result remains valid: the 55-wall H=100 prefix has an exact rank-24
hyperideal extreme ray and therefore infinite volume.  The present report
asks what happens when later Vinberg walls are added.

## Exact H=300 census

The canonical primitive mirror norms remain the proved-complete list

```text
1, 2, 3+sqrt(2), 6+2sqrt(2).
```

The search exhausts 48 control levels grouped into 33 distinct heights.  It
uses the D24 root coset and all three nonzero spinor residues in
`O_K/2O_K`, imposes the required `(3+sqrt(2))` spatial divisibility before
enumerating prime-seven norms, and removes vectors divisible by `sqrt(2)`
in the integral `D24+` glue basis.

Only six height packets are nonempty:

| height | raw reflective | imprimitive | primitive | accepted | accepted spinor |
|---|---:|---:|---:|---:|---:|
| `5+4sqrt(2)` | 4 | 1 | 3 | 3 | 0 |
| `31+22sqrt(2)` | 30 | 1 | 29 | 28 | 5 |
| `56+40sqrt(2)` | 1 | 0 | 1 | 1 | 0 |
| `62+44sqrt(2)` | 100 | 0 | 100 | 87 | 27 |
| `82+58sqrt(2)` | 5 | 0 | 5 | 5 | 1 |
| `124+88sqrt(2)` | 232 | 0 | 232 | 200 | 69 |

The first new nonempty packet after H=100 is therefore

```text
h = 56+40sqrt(2),
q = 6+2sqrt(2),
x = 12+8sqrt(2),
y = (13+9sqrt(2), 8+5sqrt(2), 5+4sqrt(2), 0^21).
```

It contains one D24-coset wall and is followed by the much denser 87-wall
packet at `62+44sqrt(2)`.

## What happens to the H=100 witness

The frozen H=100 ray is

```text
z = ((20+12sqrt(2))/7;
     2+2sqrt(2), 2+2sqrt(2), 2, 2, 0^20),

z^2 = (48-16sqrt(2))/7 > 0.
```

Its exact packetwise pairings are summarized below.  A positive pairing
means the corresponding wall cuts the ray.

| height | walls | negative | zero | positive | minimum | maximum |
|---|---:|---:|---:|---:|---|---|
| `31+22sqrt(2)` | 28 | 23 | 5 | 0 | `-12-8sqrt(2)` | `0` |
| `56+40sqrt(2)` | 1 | 1 | 0 | 0 | `-4-2sqrt(2)` | `-4-2sqrt(2)` |
| `62+44sqrt(2)` | 87 | 79 | 4 | 4 | `-28-20sqrt(2)` | `4` |
| `82+58sqrt(2)` | 5 | 4 | 0 | 1 | `-28-18sqrt(2)` | `4+2sqrt(2)` |
| `124+88sqrt(2)` | 200 | 183 | 14 | 3 | `-40-28sqrt(2)` | `2sqrt(2)` |

Thus the old witness survives through `56+40sqrt(2)` and is first cut at
`62+44sqrt(2)` by four walls.

The stronger census checks every one of the first packet's 109 hyperideal
extreme rays against cumulative later walls:

| through height | cumulative later walls | surviving old hyperideal rays |
|---|---:|---:|
| `31+22sqrt(2)` | 28 | 2 |
| `56+40sqrt(2)` | 29 | 2 |
| `62+44sqrt(2)` | 116 | 0 |
| `82+58sqrt(2)` | 121 | 0 |
| `124+88sqrt(2)` | 321 | 0 |

The table by itself is only a witness-cut result: new extreme rays can be
created during facet insertion.  In this case one such replacement ray can
be certified directly, without enumerating the full ray fan.

## A replacement hyperideal ray through H=300

Put

\[
c=\frac{172-120\sqrt2}{7}
\]

and, in the orthonormal `D24+` coordinates, set

```text
v = (1; c^19, 0^5).
```

Direct exact pairing with every wall in the 348-wall archive gives

```text
(v,r) <= 0  for all archived walls r.
```

Exactly 24 pairings vanish.  Their wall functionals have rank 24 over
`Q(sqrt(2))`, so `v` is an extreme ray, not merely a feasible point.  Its
Lorentzian norm is

\[
(v,v)=\frac{1109345-784418\sqrt2}{49}>0.
\]

Thus it is hyperideal and proves that the complete H=300 prefix has infinite
volume.  The active facets are 23 height-zero `D24` walls together with the
fourth accepted wall of the `62+44sqrt(2)` packet (global wall index 59).

## Exact rank acceleration

The dedicated double-description backend now accelerates full-row-rank
checks by specializing `sqrt(2)` to both roots modulo the split primes 17,
23, and 31:

```text
sqrt(2) = 6 or 11 mod 17,
sqrt(2) = 5 or 18 mod 23,
sqrt(2) = 8 or 23 mod 31.
```

A nonzero modular minor proves that the same minor is nonzero over
`Q(sqrt(2))`; this is an exact certificate, not a floating-point heuristic.
If every specialization is inconclusive, the code falls back to exact
Gaussian elimination over the field.  The frozen first-packet ray stream is
unchanged:

```text
b5066aaa308ccfc7a48fb8c5897f4fe4613e4683364c5d0fd2def8d81158d24d
```

Even with that acceleration, unrestricted refinement through the 28-wall
H=31 packet expanded beyond the bounded run used for this lane.  The direct
rank-24 certificate above bypasses that expansion for the single question of
whether the H=300 prefix has finite volume; it is not a full ray census.

## Reproducibility and scope

- Backend: `src/coxeter_polytopes/search/silver_p7_d24_plus.py`
- H=300 audit: `scripts/audit_silver_p7_d24_plus_h300.py`
- Tests: `tests/test_silver_p7_d24_plus.py`
- H=300 archive: `data/enumerations/silver_p7_d24_plus_h300.json`
- Pinned H=100 archive: `data/enumerations/silver_p7_d24_plus_h100.json`

The H=300 archive contains every control, candidate count, accepted wall,
the exact old-witness pairings, survivor counts, the replacement ray and its
rank certificate, and the H=100 SHA lineage.  It proves that the H=300 cone
is neither finite-volume nor compact.  It does not prove nonreflectivity:
the eventual Coxeter chamber after all later heights remains unresolved.

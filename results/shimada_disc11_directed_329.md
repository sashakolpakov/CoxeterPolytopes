# A 329-wall K3-generated Coxeter subsystem in dimension nine

**Status (2026-08-08).** Exact theorem, sealed resumable checkpoint, and full
external-atlas replay.  This advances the discriminant-11 golden K3
eigen-nef wall lower bound from 309 to 329.  One exact hyperideal ray still
survives, so this is not yet a compactness or finite-volume theorem for the
complete chamber.

## Result

Start with the frozen 309-wall subsystem on the singular discriminant-11 K3
surface, its self-adjoint golden action, and the enhanced automorphism group

\[
V_4=\{1,A,B_+,B_-\}\subset\operatorname{Aut}(X_{11},\tau).
\]

Five further escape-directed searches find five physical Néron--Severi
classes of square `-2`.  Each restriction is an integral norm-two root over
`Z[phi]`, meets the previous outer cone in a strict timelike face, and has a
four-element `V4` orbit.  Every orbit cuts the selected hyperideal-ray orbit
and remains pairwise acute with all earlier walls.  Thus

```text
309 -> 313 -> 317 -> 321 -> 325 -> 329
```

is a chain of certified effective divisor-wall lower bounds.  The five exact
Shimada gallery lengths are

```text
138, 105, 64, 108, 55.
```

The resulting 329 roots have canonical ordered hash

```text
eb11488e15a822878b2d1f9b41a0039fdc289b6d887d8fb4c829206a3d6b6f16
```

and all `329*328/2 = 53,956` pairings have recognized Coxeter type:

| relation | count |
|---|---:|
| `m2` | 2,816 |
| `m3` | 816 |
| `m5` | 472 |
| ultraparallel | 49,852 |

There are no unrecognized finite angles and no positive off-diagonal root
pairings.

## The surviving exact escape

The terminal 329-wall outer cone still contains the exact projective ray

```text
(-1,
 (113-101phi)/183,
 (-367+226phi)/183,
 (-707+450phi)/183,
 (-651+421phi)/183,
 (-620+369phi)/183,
 (-363+176phi)/183,
 (-5+2phi)/3,
 (-48-71phi)/183,
 (-388+153phi)/183).
```

Its norm is

\[
\frac{297890-184100\phi}{33489}>0,
\]

its active facets are

```text
4, 9, 13, 16, 20, 38, 58, 62, 193,
```

and every one of its 329 wall pairings is nonpositive.  Their exact hash is

```text
8c031d3c5f1fedd21e50f4fb080a11cbfb1bf8a57a79340ca6646590d0074578
```

Consequently the finite wall system is still an infinite-volume outer
approximation.  Further genuine curve walls may cut this direction; the
result does not imply that the complete eigen-nef chamber has infinite
volume.

## Independent replay and checkpoint

Before resumption, the old 309-wall archive was replayed through all 30
orbits and 2,217 Shimada atlas transitions.  The extension then generated the
five new galleries from scratch.  A second full replay of the resulting
checkpoint verifies all 35 orbits and all **2,687** atlas transitions,
reconstructs every physical `(-2)` class and field root, and recomputes the
complete Coxeter relation histogram above.  Its explicit
`frozen_309_prefix_match` flag is true; `terminal_equals_frozen_309` is false,
as it must be for the 329-wall extension.

The sealed checkpoint payload SHA-256 is

```text
c527c7eed01f2defbd3f46adbc5c7df8384e1616e9b0bbdb05500f6adec336e4
```

The serialized extension artifact has SHA-256

```text
ab44dfd516239667b729a995ca0ea9a2f6899a8ee9e2d619b5270d387abe472c
```

The separately regenerated full-replay audit has SHA-256

```text
de9cce2683ad3247c028e98c67416217dcfc5a1ebcddefb0b0bef8777e84d31e
```

Artifacts:

- [`shimada_disc11_directed_309_plus5.json`](../data/candidates/shimada_disc11_directed_309_plus5.json)
- [`shimada_disc11_directed_329_replay_audit.json`](../data/candidates/shimada_disc11_directed_329_replay_audit.json)
- [`probe_shimada_disc11_hyperideal_escape.py`](../scripts/probe_shimada_disc11_hyperideal_escape.py)
- [`test_shimada_disc11_directed_329.py`](../tests/test_shimada_disc11_directed_329.py)

With the SHA-pinned Shimada `disc11.txt`, reproduce the full replay by

```bash
PYTHONPATH=src:. python3 scripts/probe_shimada_disc11_hyperideal_escape.py \
  /path/to/disc11.txt \
  --verify-checkpoint data/candidates/shimada_disc11_directed_309_plus5.json
```

Resume the constructive boundary by replacing `N` with the desired number of
additional `V4` wall orbits:

```bash
PYTHONPATH=src:. python3 scripts/probe_shimada_disc11_hyperideal_escape.py \
  /path/to/disc11.txt \
  --resume-checkpoint data/candidates/shimada_disc11_directed_309_plus5.json \
  --additional-orbits N --skip-gallery-replay --output /path/to/resume.json
```

Skipping the old gallery replay is legitimate only after the full checkpoint
replay above has passed; every newly generated gallery is still checked
exactly.

## Scope

This is positive algebro-geometric construction progress: twenty additional
Coxeter walls arise from actual effective divisor classes on one K3 surface,
not from a free Lorentzian-vector search.  It is not a new finite-volume
polytope and does not settle the enhanced automorphism centralizer.  The two
terminal routes remain unchanged: either enough further divisor walls make a
compact outer cone, or an infinite-order enhanced automorphism proves the
complete chamber has infinite volume.

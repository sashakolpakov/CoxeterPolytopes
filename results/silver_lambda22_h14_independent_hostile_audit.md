# Silver Lambda22 H14 independent hostile audit

**Verdict: PASS.**  A clean-room replay confirms the complete simultaneous
Vinberg packet at `H=14+10sqrt(2)` for
`<-(sqrt(2)-1)> + Lambda22` in hyperbolic dimension 22.

The replay imports neither the production silver-Lambda22 search nor its
number-field implementation.  It reparses the vendored Gram matrix, obtains
determinant 12 and Smith form `1^20,2,6`, reconstructs the pairing radicals
modulo 2 and 3, and runs exact PARI/GP count-only shell censuses after LLL
reduction.  The PARI call is `qfminim(Q,b,m=1)` with the default algorithm
flag 0: `m=1` limits stored representatives but the returned total `N` is
the exact complete count.

The ordinary Lambda22 cumulative signed counts through norms 4, 6, 8, 10
are respectively

`49,896; 2,871,528; 53,204,118; 520,800,918`.

Thus the exact norm-10 shell contains **467,596,800** vectors.  For
`S2={u:G*u=0 mod 2}`, the exact shell histogram is 8,448 at norm 12,
49,896 at norm 16, and **1,368,576** at norm 20.  Consequently the complete
H14 packet has **468,965,376** roots.

The disputed `q=4` component is complete.  Its five necessary Gram types
have `v^2=10,12,14,16,20`.  Since both coefficients lie in `S2`, so does
`u-v`, and

`(u-v)^2 = 20-v^2 = 10,8,6,4,0`.

The independently certified minimum `min(S2)=12` excludes the first four;
only `u=v` remains, giving exactly the full `S2` norm-20 shell.  The ternary
gate minimum is **24**, correcting the obsolete value 36; all earlier
ternary candidate differences have norm at most 14 and remain excluded.

The ordinary, dyadic, and cross-packet pairings at their sharp proved dot
bounds are respectively

`-4-4sqrt(2), -14-12sqrt(2), -1-2sqrt(2)`,

so the simultaneous packet is mutually Vinberg-acute.  Finally, for every
minimal `d in Lambda22`, the ray `z=(1,d/3)` has positive square
`13/9-sqrt(2)` and worst packet pairings
`-(1+sqrt(2))/3` and `-4/3+2sqrt(2)/3`, both strictly negative.  Hence the
complete H14 prefix has infinite volume.  This remains a bounded-prefix
statement; the eventual chamber is unresolved.  The corrected next frontier
is the H36 mixed packet, not the obsolete H41 target.

Reproduce with:

```sh
PYTHONPATH=src python3 scripts/independent_hostile_audit_silver_lambda22_h14.py
PYTHONPATH=src pytest -q tests/test_silver_lambda22.py tests/test_silver_lambda22_h14_independent_hostile_audit.py
```

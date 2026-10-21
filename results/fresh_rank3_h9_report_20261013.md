# Rank-three cuts in the golden H9 cone

The bounded search found seven reflective walls outside the previously tested
rank-two plane. Pairing each with its exact block-swap image gives seven
compatible two-wall orbits. Their compatibility graph has exactly four maximal
cliques, each with five orbits. All four resulting **382-wall Coxeter cones have
an exact spacelike escape**. No finite-volume or compact polytope was found.

A stronger certificate closes this particular slice, independent of bounds:
one surviving escape cannot be cut by **any** additional wall normal from the
slice or its block-swap image that is compatible with the original 372 walls.
This is a restriction on a pair of slices, not on the full ambient reflection
group.

## Slice and bounded search

The host is `B=<-2phi>+E8+<1>` in the standard E8 basis, with
`phi^2=phi+1`. Put `z=(7,14,20,16,12,8,4,10)`, so `(z,z)=14`.
The slice contains roots `(x,t*z,u)` with

`14*t^2+u^2-2phi*x^2=q`, for `q in {1,2}`.

It is saturated: `3*z[0]-z[2]=1`, so integral E8 coordinates force integral t.
Freeing u/t lets this search leave the previously saturated binary plane.
The exact box was `0<t<=4000`, `0<=u<=6000`, with conjugate condition
`14*t'^2+u'^2<=q`. Reduction modulo two requires `u=a+b*phi` with
`a=q mod2`, `b=0 mod2`.

The search tested 3,847,150 parameter pairs in 17.184 seconds. An exact
quadratic inequality discarded pairs that could not cut the selected escape
before taking square roots. It found seven roots satisfying integrality,
reflection integrality, all 372 acute pairings, the strict interior, the escape
cut condition, and every finite Coxeter-angle gate. This is a complete search
of the stated parameter box, not an ambient root census. The independent
certificate audit replays the seven accepted roots; it does not independently
repeat all 3.8 million parameter pairs.

The gate is `4*(r,s)^2/(r^2*s^2) in {0,1,2,3,1+phi,2+phi}` below 4.
Values at least 4 represent parallel or disjoint hyperplanes. In particular,
an acute integral pair with angle `2pi/5` is rejected.

Numbering the seven cuts 0 through 6 in increasing old-basepoint height, only
orbit pairs (0,4) and (1,3) are incompatible. The four maximal selections and
exact positive escape norms are:

| Orbits | Walls | Escape norm |
| --- | ---: | --- |
| 0,1,2,5,6 | 382 | `(128-79phi)/12` |
| 0,2,3,5,6 | 382 | `(4507400-2767458phi)/1794005` |
| 1,2,4,5,6 | 382 | `61152-37794phi` |
| 2,3,4,5,6 | 382 | `(10071-6224phi)/121` |

Every wall pair in every cone passes the exact Coxeter gate; the original
372-wall checks are reused across the four cones. Each cone has a strict
timelike interior and a feasible spacelike extreme ray with nine independent
active facets. These are exact infinite-volume certificates for these finite
cones.

## Bound-independent slice obstruction

Let S be the ternary slice, U the stored integral block-swap involution, and
v the escape in the first row above. Represent orthogonal projection onto S
by `(time,(z,w)/14,last_scalar)`. Let `w_i` denote the zero-based original
372-wall list. The exact identities are

```
P_S(v) = (-8+5phi) P_S(w_107)
       + (-16/3+10phi/3) P_S(w_1)
       + (-11/2+7phi/2) P_S(w_8),

P_S(Uv) = (9-11phi/2) P_S(w_101)
        + (-21/2+13phi/2) P_S(w_8).
```

All displayed coefficients are positive. Also U permutes the original wall
set. If r belongs to S and `(r,w_i)<=0` for every original wall, the first
identity gives `(r,v)<=0`. If r belongs to U(S), apply the second identity to
Ur, using the permutation property. Thus v survives every such compatible
cut in either slice, without any arithmetic or height restriction.

## Reproduction and validation

- [Bounded search](../scripts/fresh_rank3_h9_search.py) and [seven cuts](fresh_rank3_h9_search_20261013.json).
- [Compatible-cone construction](../scripts/fresh_rank3_h9_cones.py) and [four cones plus obstruction](fresh_rank3_h9_cones_20261013.json).
- [Independent Fraction verifier](../scripts/fresh_rank3_h9_verify.py) and [passing audit](fresh_rank3_h9_audit_20261013.json).
- [Regression tests](../tests/test_fresh_rank3_h9.py) compare the bounded integer alphabet against an independent coefficient box and reject a forged obstruction coefficient.

```sh
PYTHONPATH=src:. python3 scripts/fresh_rank3_h9_search.py --t-bound 4000 --u-bound 6000
PYTHONPATH=src:. python3 scripts/fresh_rank3_h9_cones.py
PYTHONPATH=. python3 scripts/fresh_rank3_h9_verify.py \
  results/fresh_rank3_h9_search_20261013.json results/fresh_rank3_h9_cones_20261013.json \
  --output results/fresh_rank3_h9_audit_20261013.json
PYTHONPATH=src:. python3 -m pytest -q tests/test_fresh_rank3_h9.py
```

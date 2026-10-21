# An explicit abstract completion of H9/14 record 47

## Subsequent whole-shadow theorem

This selected completion is now known to be the unique residual incidence
orbit under the compact H9/14 missing-face necessities.  The exhaustive
all-schedule census and the already frozen unrestricted label obstruction
therefore exclude the entire record-47 graph shadow in the compact case; see
[`h9_14_record47_whole_shadow_compact_no_go.md`](h9_14_record47_whole_shadow_compact_no_go.md).
No finite-volume noncompact conclusion follows.

H9/14 second-layer survivor record 47 now has an explicit abstract
completion. The missing-edge graph is

```text
M?????????O?o?g?_
```

with missing edges

\[
(0,11),(0,12),(1,12),(0,13),(2,13),(12,13),
\]

and four additional minimal A5s

\[
\begin{aligned}
&(3,4,5,6,7),\\
&(1,2,8,9,10),\\
&(1,8,9,10,11),\\
&(2,8,9,10,11).
\end{aligned}
\]

The complete blocker, reconstructed directly over all `2^14` subsets, is
pure of order five with exactly 90 members. Their complements generate a
complex with face vector

\[
(14,85,300,690,1088,1190,880,405,90).
\]

Every one of its 405 ridges has degree two, the dual graph is connected, and
the generated face table has exactly the declared six missing edges and four
A5s. All strict local `g`, Macaulay, Euler, and mod-two link-homology audits
pass.

The facet-complement and packed-face-table digests are, respectively,

```text
8306d9c6c39a90aae519bd5ae587743b3737371f9efbe01ebd3757105b0b56ed
4d07f78cba8551c9d0fe9fc4d6041f3df8ae35c5017e9942e10005686c628aba
```

A uniform rank-four Gale search also succeeds. The archived 1,001-sign word
passes all 45,045 Grassmann--Pluecker sign relations and has exactly the 90
facet complements as its positive five-circuits. This is a necessary
oriented-matroid condition only and does not prove that this particular word
is realizable.

A subsequent exact multiwedge factorization independently proves that this
face lattice is the nerve of an ordinary convex simple 9-polytope; see
`results/h9_14_multiwedge_polytopality.md`. It supplies neither Coxeter labels
nor a Lorentzian Gram realization. The current theorem-facing partition is
`121 + 6 + 15`: records 43, 44, 47, and 35 join the two product fixtures in
the ordinary-polytopal class, while the later exact record-34 exclusion leaves
15 bounded-timeout shadows undecided.

The frozen certificate digest is

```text
960ef6564c6b8013caf0989c95a46cc8401c846d32d47c66c057091340301546
```

and the Gale sign digest is

```text
b8a9b3c5a8b4712e908892a53e75934e49db2e6a2fa6856c1428a38c8d0749b2
```

Reproduce the solver-independent archive verification and hostile replay with

```sh
PYTHONPATH=.:src python3 scripts/certify_h9_14_record47.py
PYTHONPATH=.:src python3 scripts/independent_h9_14_record47_audit.py
PYTHONPATH=.:src python3 -m pytest -q tests/test_h9_14_record47.py
```

Only `--write-certificate` reruns the Gale search.

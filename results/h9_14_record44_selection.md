# An explicit abstract completion of H9/14 record 44

Subsequent work proves that this is the unique residual incidence orbit of
the record-44 graph shadow under the compact missing-face necessities and
then transfers the unrestricted selected-nerve label obstruction to the
whole shadow; see
[`h9_14_record44_whole_shadow_compact_no_go.md`](h9_14_record44_whole_shadow_compact_no_go.md).
The upgrade is compact-only and makes no finite-volume noncompact claim.

## Result

H9/14 second-layer survivor record 44 has an explicit abstract facet
completion. Its missing-edge graph is

```text
M?????????O?o?o@?
```

with missing edges

\[
(0,11),(0,12),(1,12),(0,13),(1,13),(11,13),
\]

and four additional minimal nonfaces

\[
\begin{aligned}
&(2,3,4,5,6),\\
&(1,7,8,9,10),\\
&(7,8,9,10,11),\\
&(7,8,9,10,12).
\end{aligned}
\tag{1}
\]

The complete blocker of these ten nonfaces is pure of order five and has
exactly 90 members. Their complements define an eight-dimensional simplicial
complex on fourteen vertices with face vector

\[
(14,85,300,690,1088,1190,880,405,90).
\tag{2}
\]

This moves record 44 from an unresolved bounded timeout to an **explicit
globally valid abstract completion**. A subsequent exact multiwedge
factorization proves that this face lattice is the nerve of an ordinary
convex simple 9-polytope; see
`results/h9_14_multiwedge_polytopality.md`. No Coxeter labels or Lorentzian
Gram matrix are supplied.

## Global audit

The certificate independently enumerates all `2^14` potential transversals;
there are no blocker members of any order other than five. It then generates
the full face table and verifies:

- precisely the six missing edges and four A5s in (1) are the minimal
  nonfaces;
- all 90 facet complements lie in the original support core;
- all 405 ridges have degree two;
- the dual graph is connected;
- the face vector is exactly (2), matching the `a=0` Dehn--Sommerville
  schedule;
- every strict local `g` and Macaulay constraint holds;
- every link has spherical Euler characteristic;
- every link has the expected mod-two homology of a sphere.

The facet-complement digest is

```text
a706233135e61171db1da19172f4cb5234191f0a0a703a3388b4a4d9a9877917
```

and the packed full-face-table digest is

```text
449e514bd9f090a3c52dba38db833014e696c7d9cb20753a6b8f4f7390cb6ad6
```

These tests establish a connected mod-two homology-manifold candidate. They
do not by themselves prove that the complex is a PL sphere; the later
multiwedge realization independently proves that the face lattice is a
polytopal sphere.

## Uniform rank-four Gale gate

Because a putative nine-polytope has fourteen facets, the dual realization
would have fourteen vertices and Gale rank four. An exact Z3 search found a
uniform rank-four chirotope with 1,001 basis signs. The archived word is then
checked without a solver:

- all 45,045 rank-four Grassmann--Pluecker sign relations hold;
- the positive five-circuits are exactly the 90 facet complements;
- every ground-set element lies in, and is avoided by, a positive circuit.

Thus record 44 passes the necessary uniform oriented-matroid gate. This does
**not** prove that this oriented matroid is realizable over the reals, so the
Gale word itself is not a polytopality certificate. The independent
multiwedge construction proves polytopality of the face lattice without
claiming to realize that word.

The current theorem-facing partition is

\[
121\text{ independently impossible}
\; +\;6\text{ explicit ordinary-polytopal selections}
\; +\;15\text{ theorem-undecided}.
\]

Records 43, 44, and 47 are the three Delzant multiwedge selections added to
the two product fixtures; record 35 is the later rational multiwedge. They are
outside the undecided class; after the later record-34 exclusion, the remaining
15 are unresolved bounded-timeout shadows.

## Reproduction

The frozen archive is
`data/enumerations/h9_14_record44_selection.json`, with certificate digest

```text
e66e493d0c68f7ae28ffe03e63f96209276a8d83a5ceb27c498f35b42c4a1867
```

and Gale sign digest

```text
2e7c05277d8f34686093fc27b5c3a170be321edefad31c8eff586943aaa39444
```

Replay the frozen certificate and the independently hard-coded hostile audit:

```sh
PYTHONPATH=.:src python3 scripts/certify_h9_14_record44.py
PYTHONPATH=.:src python3 scripts/independent_h9_14_record44_audit.py
PYTHONPATH=.:src python3 -m pytest -q tests/test_h9_14_record44.py
```

Only `--write-certificate` invokes Z3 to search for another Gale sign word.
Ordinary verification imports no search backend.

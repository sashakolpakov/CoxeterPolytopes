# Independent common-escape audit of the H9 cutter dictionary

An exact spacelike vector satisfies **all 382 old walls and all 670 direct or block-swap cutter roots**. Consequently no selection from this finite dictionary can complete a finite-volume polytope, regardless of how its pair conflicts are resolved. This does not exclude additional roots outside the bank.

Write `φ=(1+√5)/2`. The common escape is

\[
v=(1,16/5,32/5,48/5,39/5,6,4,2,5,1/3-\phi/12),
\qquad
v^2=2389/720-295\phi/144>0.
\]

The archived Gram form and coordinates are those of `results/fresh_rank4_h9_dictionary_20261013.json`. Independent `Fraction` arithmetic checks every one of the 1,052 nonpositive wall pairings, the nine independent active normals, the positive time coordinate, and the displayed positive squared norm. The dictionary file is anchored by its SHA-256 digest.

The common strict timelike center `c` satisfies

\[
c^2=621/10000-2\phi<0,
\qquad \max_r(c,r)=-1/100.
\]

For an explicit strictly feasible escape, set `x=v+c/2048`. Exact replay gives

\[
\begin{aligned}
x^2&=1252687467989/377487360000
      -967578209\phi/471859200>0,\\
(x,c)&=27240263/61440000-153739\phi/76800<0,\\
\max_r(x,r)&=-1/204800<0.
\end{aligned}
\]

The segment from `c` to `x` crosses the future null cone while every wall inequality remains strict. Thus the intersection contains an open set of ideal directions and has infinite volume. Deleting any dictionary inequalities preserves this same witness.

A numerical LP proposed the original escape. The saved certificate and strict perturbation are checked without numerical optimization, project `Golden` arithmetic, external solvers, or network access:

```sh
PYTHONPATH=. python3 scripts/fresh_rank4_h9_selection_audit_20261013.py --verify
```

The replay passed. The exact artifact is `results/fresh_rank4_h9_selection_audit_20261013.json`. A corrupted copy replacing the spacelike vector with the timelike center was rejected.

A separate integer audit reconstructs all 224,115 potential pairs of the 670-node graph and finds exactly the archived 340 conflicts. It checks that all 12 saved greedy selections are pairwise compatible and maximal under adding one root, and independently replays all their escape certificates. Its script and artifact are `scripts/fresh_rank4_h9_selection_audit_graph_20261013.py` and `results/fresh_rank4_h9_selection_audit_graph_20261013.json`.

This finite-dictionary result is separate from the stronger two-slice projection obstruction in `results/fresh_rank3_h9_rank4_obstruction_20261013.json`. That replay also passed here; an additional direct check using the inverse E8 Gram confirms that its four restricted coordinates really are pairings against time, fundamental weights 4 and 7, and the scalar basis vector. Both strictly positive cone identities were then checked against those explicit basis vectors. The stronger proof excludes completion using the two specified slices at any height; neither result establishes nonreflectivity of the ambient arithmetic host.

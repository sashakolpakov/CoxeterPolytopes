# Exact bound-10 norm-two Vinberg prefix on the golden OG10 carrier

> **Subsequent extension.**  The search now continues through control
> coordinate `4+5phi`, with 1,144 walls and another exact spacelike
> replacement ray.  See
> [`og10_golden_vinberg_through_4_plus_5phi.md`](og10_golden_vinberg_through_4_plus_5phi.md).
> A later symmetry theorem, independent of bounded enumeration, proves the
> complete carrier chamber infinite-volume; see
> [`og10_golden_h9_infinite_stabilizer.md`](og10_golden_h9_infinite_stabilizer.md).

## Result

For

\[
H_{OG}=\begin{pmatrix}-2\phi&1\\1&2\end{pmatrix}
       \oplus E_8\oplus\langle2\rangle,
\]

the norm-two Vinberg search is now exhaustive through identity control bound

\[
0<x\leq 10.
\]

The bound-eight search had 171 walls.  The only two further admissible
control coordinates are

\[
x=2+4\phi,\qquad x=3+4\phi.
\]

Exact enumeration gives respectively 150 and 136 new simple roots, hence a
457-wall prefix.  All \(\binom{457}{2}=104{,}196\) pairs are recognized:

| relation | pair count |
|---|---:|
| `m=2` | 4,413 |
| `m=3` | 1,810 |
| `m=5` | 3,266 |
| ultraparallel | 94,707 |

The 457-wall outer cone still has an exact spacelike extreme ray.  Therefore
this recorded cone is neither compact nor finite-volume.  This remains a
bounded-prefix theorem: a norm-two root at a later height may cut the ray.

## Exhaustive extension transcript

The positive-definite reduction and height-zero dominance conditions are the
same as in the bound-eight certificate.  Starting from its 161 positive-height
roots, the exact Schur-complement recursion records:

| control coordinate | new roots | visited nodes |
|---|---:|---:|
| \(2+4\phi\) | 150 | 237,278 |
| \(3+4\phi\) | 136 | 319,736 |

Thus the extension visits 557,014 nodes.  Including the archived bound-eight
run, the complete bound-ten computation visits 717,248 nodes.  No
floating-point filter is used to accept or reject a root.  Every recorded
root is integral, has field norm two, and is acute with all earlier simple
roots.  The ordered 457-root digest is

```text
e86b65bb5348f5c0e645e482b26eb6006c20245a578249fea8c8707c7bbdcd9b
```

## What happens to the bound-eight escape ray

Of the 150 walls at \(x=2+4\phi\), exactly 147 pair negatively with the old
bound-eight hyperideal ray, two contain it, and one cuts it.  Thus merely
reusing the old certificate would be invalid.

The exact local cone update gives a replacement extreme ray with active
facets

\[
(0,1,2,3,4,5,6,7,11,225).
\]

In the standard field basis it is

\[
\begin{aligned}
v={}&(1,-1/2,
135/158+50\phi/79,
135/79+100\phi/79,\\
&405/158+150\phi/79,
162/79+120\phi/79,
243/158+90\phi/79,\\
&81/79+60\phi/79,
81/158+30\phi/79,
108/79+80\phi/79,
3/4+\phi/4).
\end{aligned}
\]

Its exact norm is

\[
(v,v)=\frac{27755}{24964}-\frac{32489}{49928}\phi>0.
\]

Every one of the subsequent 136 walls at \(x=3+4\phi\) pairs strictly
negatively with this ray.  Hence the same replacement ray certifies the full
bound-ten prefix.  Its ordered 457-pairing digest is

```text
8700c4bf56ee8e7619a68916992cf8fe6d1fd7231b4af2e230bcae0f86e4f629
```

The independent verifier reconstructs the ray from the ten active facets,
orients it using all 457 inequalities, confirms that no additional facet is
active, and checks its positive norm.

## Geometric scope

This calculation concerns the norm-two reflection arrangement, which is the
movable Coxeter arrangement on the golden OG10 eigenspace.  The two additional
nonreflective square-`-4` mirror orbits are Kähler-chamber subdivision data;
they are not missing roots of this norm-two Vinberg run.  Consequently the
precise conclusion here is:

- the norm-two simple-root list is complete through \(x\leq10\);
- its 457-wall outer cone is not finite-volume; and
- neither the completed norm-two chamber nor the Kähler subdivision has yet
  been proved finite-volume or infinite-volume.

## Reproducibility

The generator is
[`scripts/probe_og10_golden_vinberg_extension.py`](../scripts/probe_og10_golden_vinberg_extension.py).
It first freezes a one-level checkpoint and then resumes the second level
without repeating the expensive first enumeration:

```sh
PYTHONPATH=.:src python3 scripts/probe_og10_golden_vinberg_extension.py \
  --max-new-levels 1 \
  --output data/candidates/og10_golden_vinberg_through_2_plus_4phi.json \
  --summary-only

PYTHONPATH=.:src python3 scripts/probe_og10_golden_vinberg_extension.py \
  --max-new-levels 2 \
  --output data/candidates/og10_golden_vinberg_through_3_plus_4phi.json \
  --summary-only
```

The final exact artifact is
[`data/candidates/og10_golden_vinberg_through_3_plus_4phi.json`](../data/candidates/og10_golden_vinberg_through_3_plus_4phi.json).

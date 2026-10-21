# Exact OG10 norm-two Vinberg prefix through `4+5phi`

> **Subsequent decision.**  A direct fixed-H9 recursion has exhausted
> `3+6phi` and the hyperideal ray survives all 71 fixed walls.  More strongly,
> an infinite complete-chamber stabilizer proves the H9 and H10 chambers
> infinite-volume without any wall-completion claim.  See
> [`og10_golden_h9_vinberg_through_3_plus_6phi.md`](og10_golden_h9_vinberg_through_3_plus_6phi.md)
> and
> [`og10_golden_h9_infinite_stabilizer.md`](og10_golden_h9_infinite_stabilizer.md).

## Scope first

This is an exhaustive theorem through the control coordinate

\[
                         x=4+5\phi,
\]

not through every identity value below `13`.  The candidate generator was
called with identity bound `13`, but its next admissible value is

\[
3+6\phi\approx12.7082<13,
\]

and that fifth level was deliberately not computed.  Roots at
\(3+6\phi\) and later may cut the final ray below.

Within that scope, the enumeration and the non-finite-volume conclusion for
the recorded outer cone are exact.

## Four completed levels

Starting from the 457-wall checkpoint through \(3+4\phi\), the next four
Vinberg levels are:

| control coordinate | roots | exact recursion nodes | level-root SHA-256 |
|---|---:|---:|---|
| \(2+5\phi\) | 71 | 131,946 | `b7cb33523501cc78cba2a6ef23933a78d87386c63ecff58644795a44afe1bf8d` |
| \(4+4\phi\) | 13 | 6,732 | `89c0432bc56791664d39e1a426a2f10da56c4a8f7f1e7eed2cd96d6d25c250bc` |
| \(3+5\phi\) | 371 | 1,143,573 | `ed2a48f27fe024be5cf6dbeba7d4ff8da9dde6563fd19622ca731ae94659cb55` |
| \(4+5\phi\) | 232 | 650,813 | `4b2189fa0c04e321ac4294d11d08da3af236ef31426fdcfd24ca92ac81c49cc4` |

Thus the extension adds 687 roots in 1,933,064 exact recursion nodes, for a
total of

\[
                              1144
\]

oriented norm-two simple roots.  The full ordered-root digest is

```text
d4e898a0c4310108c95ac790bc9f7860a239b266b6a9bd051bdca91e944b9987
```

The artifact records, at every depth of every level, the visited count,
Schur-complement pruning count, dominant-wall pruning count, parity pruning,
and leaf norm mismatches.  The regression test checks the ordered nesting,
each level digest, every root's integrality and norm, and the ray history.

All \(653,796\) pairs have an allowed norm-two Coxeter relation or are
ultraparallel:

| relation | pairs |
|---|---:|
| \(m=2\) | 11,629 |
| \(m=3\) | 5,459 |
| \(m=5\) | 10,897 |
| ultraparallel | 625,811 |

## The ray is cut twice and survives twice

At every level, all new walls were tested exactly against the current
spacelike extreme ray.  A cutting event triggers a local reconstruction from
the cutting wall and nine inherited active facets; the result is then replayed
independently against every wall in the enlarged prefix.

| completed through | new-wall signs at previous ray | resulting active facets | ray norm status |
|---|---|---|---|
| \(2+5\phi\) | 70 negative, 1 positive | `(0,1,2,3,4,5,6,7,11,463)` | spacelike |
| \(4+4\phi\) | 13 negative | unchanged | spacelike |
| \(3+5\phi\) | 371 negative | unchanged | spacelike |
| \(4+5\phi\) | 231 negative, 1 positive | `(0,1,2,3,4,5,6,7,11,912)` | spacelike |

After the first cut the replacement norm is

\[
\frac{27091}{6724}-\frac{33049}{13448}\phi
\approx0.05261709591>0.
\]

The final replacement ray is

\[
\begin{aligned}
(&1,-\tfrac12,
-\tfrac{115}{41}+\tfrac{235}{82}\phi,
-\tfrac{230}{41}+\tfrac{235}{41}\phi,
-\tfrac{345}{41}+\tfrac{705}{82}\phi,\\
&-\tfrac{276}{41}+\tfrac{282}{41}\phi,
-\tfrac{207}{41}+\tfrac{423}{82}\phi,
-\tfrac{138}{41}+\tfrac{141}{41}\phi,
-\tfrac{69}{41}+\tfrac{141}{82}\phi,\\
&-\tfrac{184}{41}+\tfrac{188}{41}\phi,
\tfrac34+\tfrac14\phi).
\end{aligned}
\]

Its exact norm is

\[
\boxed{
\frac{39643}{6724}-\frac{48969}{13448}\phi
\approx0.003903450692>0.}
\]

Every one of the 1,144 wall pairings is nonpositive, exactly ten are zero,
and the active normals have full rank ten.  Hence this is an exact spacelike
extreme ray.  The 1,144-wall outer cone is therefore not finite-volume.

The much smaller positive norm is useful boundary evidence: successive cuts
are pushing this particular escape direction toward the light cone.  It is
not evidence that later walls must close it; the next level may cut it again
or produce a different hyperideal direction.

## Amplification by the complete-chamber involution

The complete-chamber involution from
[`og10_golden_complete_chamber_involution.md`](og10_golden_complete_chamber_involution.md)
maps every one of the original 457 walls into this 1,144-wall archive.  On the
larger archive it retains 771 images and transports 373 walls farther out:

| transported control coordinate | walls |
|---|---:|
| \(3+6\phi\) | 106 |
| \(4+6\phi\) | 171 |
| \(5+6\phi\) | 5 |
| \(4+7\phi\) | 91 |

Adjoining those images gives an involution-invariant set of 1,517 distinct
norm-two complete-chamber walls.  All 1,149,886 pairs in that closed set are
exactly acute.  The new and closed-set digests are

```text
373 transported roots:
10cf5a825b4e2700283c63c6f9b65b9d894cf912fabc73f540c98e42fddbd585

1517-root involution closure:
ecc3ef936d3c8766899420d5917ba239e642608fa6b13e07a97f3a62b0663039
```

These 373 roots are genuine complete-chamber walls because they are symmetry
images.  They do **not** exhaust any of the four newly reached control levels.
In particular, the 106 roots at the next level \(3+6\phi\) are a certified
subset, not the result of a completed Vinberg enumeration at that height.

### Independent check on the `h`-fixed `H9` restriction

As a separate bounded cross-check, restrict the 1,144 inequalities to the
fixed hyperplane of \(h\), namely `y=0`.  Its field lattice is

\[
\langle-2\phi\rangle\oplus E_8\oplus\langle2\rangle.
\]

Exact positive-proportional reduction gives 789 projected classes.  Of
these, 728 have timelike normals and are inactive on the selected hyperbolic
component; the remaining 61 have spacelike norm-two normals.  Independently
of the full-slice analysis, the nine restricted walls represented by original
wall indices

\[
(2,3,4,5,6,7,14,29,118)
\]

give the exact ray

\[
\begin{aligned}
(&1,2\phi,
\tfrac{10}{19}+\tfrac{69}{19}\phi,
\tfrac{20}{19}+\tfrac{100}{19}\phi,
\tfrac{16}{19}+\tfrac{80}{19}\phi,
\tfrac{12}{19}+\tfrac{60}{19}\phi,\\
&\tfrac{8}{19}+\tfrac{40}{19}\phi,
\tfrac{4}{19}+\tfrac{20}{19}\phi,
\tfrac{14}{19}+\tfrac{51}{19}\phi,
-\tfrac{13}{38}+\tfrac{11}{38}\phi),
\end{aligned}
\]

of norm

\[
\frac{625}{361}-\frac{761}{722}\phi
\approx0.0258672>0.
\]

Its other 52 restricted-wall pairings are strictly negative; the complete
pairing digest is
`f1b44caea280f363e133e62f7addd85bafc2f6af7c6d0c5f75d9ed14fe99c109`.
Thus the fixed slice of this **bounded** prefix is also not finite-volume.
This cross-check neither replaces the global analysis in
[`og10_golden_h9_fixed_slice.md`](og10_golden_h9_fixed_slice.md) nor asserts
by itself that later walls cannot cut the ray.

## Exact conclusion and nonclaims

The proved bounded statement is:

> The norm-two Vinberg prefix through \(4+5\phi\) has 1,144 walls and an
> independently verified spacelike extreme ray.  Its outer cone is not
> finite-volume.

It does not prove that the complete movable chamber is infinite-volume.  It
does not decide whether the norm-two reflection subgroup has finite index,
whether the enhanced birational group is finite, or whether later walls make
the complete chamber compact.  It contains no claim about the two additional
nonreflective square-`-4` Kähler-wall families.

## Reproducibility

The generator is
[`scripts/probe_og10_golden_vinberg_extension_4_plus_5phi.py`](../scripts/probe_og10_golden_vinberg_extension_4_plus_5phi.py),
the 1,144-root frozen artifact is
[`data/candidates/og10_golden_vinberg_through_4_plus_5phi.json`](../data/candidates/og10_golden_vinberg_through_4_plus_5phi.json),
and the independent replay is
[`tests/test_og10_golden_vinberg_extension_4_plus_5phi.py`](../tests/test_og10_golden_vinberg_extension_4_plus_5phi.py).

```sh
PYTHONPATH=.:src python3 -u \
  scripts/probe_og10_golden_vinberg_extension_4_plus_5phi.py \
  --progress \
  --output data/candidates/og10_golden_vinberg_through_4_plus_5phi.json \
  --summary-only

PYTHONPATH=.:src pytest -q \
  tests/test_og10_golden_vinberg_extension_4_plus_5phi.py
```

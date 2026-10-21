# Three new silver walls at `x=3+sqrt(2)`

The intrinsic norm-two Vinberg search on the designed silver K3 carrier has
now been exhausted through the next admissible coordinate

\[
x=3+\sqrt2.
\]

This shell is nonempty: it adds exactly three walls relative to the complete
lower-height prefix, taking the acute prefix from 34 to 37 walls.  None of
the three is the missing second finite endpoint of the targeted
`A2+F4+B2` edge.

The frozen manifest is
`data/candidates/silver_k3_norm2_vinberg_extension_3_plus_sqrt2.json`; its
canonical digest is

```text
172301cae3723d15b80425c7303ab12638574a84e69a5cc1ba28dd19179cf4a1
```

An independent reversed-coordinate Schur replay has digest

```text
5a2da5690108b4cbf3c90ada01577c5d37b74b9ca1405b159a32e025689634f9
```

## Exact shell enumeration

Write `s=sqrt(2)`.  The norm equation in finite-simple-root weight
coordinates is

\[
z^T K^{-1}z=2+(1+s)(3+s)^2=25+17s.
\]

The exact two-embedding square bounds give the ten-element alphabet

```text
0, 1, 1+s, 2+s, 2+2s, 3+2s, 4+2s, 3+3s, 4+3s, 5+3s.
```

The production left-to-right recursion gives

```text
visited by depth:         1,10,100,200,250,250,1170,2540,1550,570
Schur-pruned by depth:    0,0,80,175,225,133,916,2385,1493,345
exact weight solutions:  12
nonintegral roots:       0
rejected by old walls:   9
accepted at this shell:  3
```

The hostile replay fixes the nine weight coordinates in the opposite order,
rebuilds all prefix Schur complements, and independently obtains the same
12 reconstructed roots, the same nine blocker witnesses, and the same three
accepted roots.  Its traversal is

```text
visited by depth:         1,10,100,690,1930,1210,400,1970,4350,6960
Schur-pruned by depth:    0,0,31,497,1809,1170,203,1535,3654,6735
```

Thus the result is an exact finite enumeration, not a numerical search or a
timeout.

## The three roots and simultaneous admission

In the frozen free `Z[sqrt(2)]` basis, the new roots at global indices
34, 35, and 36 are

```text
r34 = (3+s, 4+2s, 4+4s, 8+6s, 8+6s,
       2+s, 4+2s, 3+3s, 2+2s, 1+s)

r35 = (3+s, 4+2s, 4+4s, 8+6s, 8+6s,
       3+2s, 6+4s, 6+4s, 3+2s, 0)

r36 = (3+s, 4+2s, 4+4s, 8+6s, 8+6s,
       2+2s, 3+3s, 4+2s, 2+s, 1+s).
```

Every pairing with the 34 lower-height walls is nonpositive.  The
same-height Gram matrix is

\[
\begin{pmatrix}
2&-3-3s&-s\\
-3-3s&2&-3-3s\\
-s&-3-3s&2
\end{pmatrix}.
\]

Hence sequential versus simultaneous admission makes no difference: no one
of the three cuts off either of the others.  They are certified simple walls
relative to the exhausted lower-height prefix and same-height cohort.  This
terminology does **not** assert that the 37-wall list is complete.

The full wall Gram matrix remains acute and has exact inertias

```text
identity embedding:   (9,1,27)
conjugate embedding:  (10,0,27).
```

## The targeted edge still has one finite endpoint

The edge remains

```text
E = (1,2,4,5,6,7,9,12),   type A2+F4+B2.
```

The pairings of the three new roots with `E`, in that order, are

```text
r34: (0,0,0,0,0,-1-s,0,-2-s)
r35: (0,0,0,-1,-1,0,0,-8-5s)
r36: (0,0,-1-s,0,0,0,0,-2-s).
```

For `E union {r34}`, `E union {r35}`, and `E union {r36}`, respectively,
the identity-embedding determinants are

```text
-60-48s,  -744-528s,  -60-48s,
```

and all three identity inertias are `(8,1,0)`.  Therefore none is elliptic
of rank nine.  Wall 33 remains the only certified finite endpoint, and the
second endpoint status is still

```text
not_found_through_this_shell.
```

No ideal endpoint can be certified from a bounded prefix fragment: Vinberg's
criterion requires its full componentwise-affine zero-neighbourhood in the
complete wall system.

## Literal status

The next admissible control coordinate is

\[
x=2+2\sqrt2.
\]

The formal status remains

```text
undecidable_incomplete_simple_roots.
```

There is still no proof of simple-root completeness, no global audit of all
elliptic edges and affine neighbourhoods, and therefore no conclusion of
compactness, finite volume, or infinite volume.  Reflective norms other than
two also remain outside this enumeration.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/certify_silver_k3_vinberg_3_plus_sqrt2.py
PYTHONPATH=src:. python3 scripts/independent_silver_k3_vinberg_3_plus_sqrt2.py
PYTHONPATH=src:. pytest -q tests/test_silver_k3_vinberg_3_plus_sqrt2.py
```

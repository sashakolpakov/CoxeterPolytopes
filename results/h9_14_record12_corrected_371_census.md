# Record 12: corrected exact 371-orbit abstract census

> **Subsequent compact-label theorem.**  Complete Lanner-placement and
> actual-vertex joins leave zero compact Coxeter label states on all 371
> residual orbits, with two no-import replays.  Thus record 12 is now excluded
> as a whole compact Coxeter graph shadow; see
> [`h9_14_record12_all_coxeter_label_gate.md`](h9_14_record12_all_coxeter_label_gate.md).
> The abstract PL/convex status and finite-volume noncompact problem remain
> separate.

## Certified result

For the compact-`H^9`, 14-facet support row whose missing-edge graph is
`K_{2,2}+10I` (survivor record 12), the complete graph-shadow search leaves
exactly **371 abstract residual incidence-complex orbits**.  If `a` denotes
the number of missing four-faces, their distribution is

| `a` | abstract residual orbits | facets |
|---:|---:|---:|
| 0 | 48 | 130 |
| 1 | 95 | 128 |
| 2 | 144 | 126 |
| 3 | 63 | 124 |
| 4 | 21 | 122 |

There are no survivors in the scheduled branches `a=5,6,7`; the branches
`a=8,9` were already excluded by the strict summed vertex-link `g_4` bound.

This is an exact **abstract combinatorial residual census** under the
published compact-`H^9` missing-face reductions.  Every one of the 371 rows
passes the complete blocker, face-vector, ridge-degree, dual-connectivity,
local `g`/Macaulay, all-face Euler-link, and all-face mod-two homology-link
checks.  A mod-two homology-sphere link audit is not a PL-sphere proof.

In particular, this result does **not** prove ordinary convex-polytopal
realizability, provide Coxeter labels or a Lorentzian Gram matrix, or produce
a compact or finite-volume hyperbolic Coxeter polytope.  A separate rational
certificate realizes one `a=2` member as an ordinary simplicial polytope;
ordinary polytopality of the other 370 rows is open here.

## Exhaustive pipeline

The corrected search has the following exact stages:

| stage | count |
|---|---:|
| triangle-free graphs `X` on ten isolate vertices | 12,172 |
| normalized switching cuts | 56,310 |
| locally ridge-feasible `(X,p,y)` states | 1,391,382 |
| states surviving the proved one-sided filter | 6,184 |
| colour-preserving branch orbits | 1,308 |
| profile/facet-compatible ordered type pairs | 22,858 |
| distinct relative embeddings, using both cut orientations | 2,725,719 |
| embeddings passing the complete coupled blocker schedule | 42,000 |
| full `D8 x S10` orbits | 371 |

All 42,000 global passes use the same ten balanced branch types, with a
`5+5` cut and `y=0`.

## Why the one-sided filter is sound

Write the two facet-complement branches as

```text
P = { A union t : t in T_P },
Q = { B union u : u in T_Q }.
```

If a minimal transversal `r` of `T_P` meets `B`, it hits every member of `Q`
through `B`.  It hits `P` by definition.  For every vertex of `r`, minimality
inside `T_P` supplies a member of `P` missed after that vertex is removed;
because `r` avoids `A`, the same witness proves global minimality.  Thus `r`
is already a minimal transversal of `P union Q`.  In a valid scheduled
completion it must consequently have order four or five.

Only branch-blocker rows meeting a special vertex are filtered this way.
Rows supported solely on the ten isolates are retained until both branches
have been coupled and the complete blocker has been formed.

## Correction of the provisional 344 count

The older 344-row archive was explicitly quarantined because it imposed a
one-branch purity condition before coupling, which was not proved necessary
for isolate-only blocker rows.  The corrected enumeration removes that
condition and also enumerates both normalizations of the second cut.

The latter normalization acts as relative transposition of the second
`5 x 5` branch matrix.  Including it adds 27 full orbits:

```text
a=2: +15
a=3:  +9
a=4:  +3
```

As a compact cross-check, adjoining the transposed second layer to the old
2,418 fixed-orientation representatives and applying the unchanged full
canonicalizer gives the same 371 keys as the new 42,000-row census.

The withdrawal record for 344 is preserved in
[`h9_14_record12_residual_hostile_audit.md`](h9_14_record12_residual_hostile_audit.md).
The independent replay of the corrected census is documented in
[`h9_14_record12_corrected_371_hostile_audit.md`](h9_14_record12_corrected_371_hostile_audit.md).

## Frozen artifacts

- Archive: `data/enumerations/h9_14_record12_residual.json`
  - file SHA-256:
    `faca367738bba58afc0ec78a43e9d4aa4b02721802c5d9e77b7d155c8e7faa05`
  - internal certificate SHA-256:
    `50ca012334e7a4973a4c29bde4226d4460bc8ade2d644a3b9dbf56b1bf92e8f4`
- Certificate: `data/enumerations/h9_14_record12_residual_certificate.json`
  - internal certificate SHA-256:
    `c759980c4a18cbc9c5c5cdf892b7fe68a142bec472f7b8be5f7e9abf9b6c9a98`
- Production verifier:
  `src/coxeter_polytopes/combinatorics/h9_14_record12_residual.py`
- Independent semantic topology digest:
  `6a58b02d39324916d24de5d43a8b1091393f83cadc7963de8581023266c66c84`

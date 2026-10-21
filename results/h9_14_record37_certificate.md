# Record 37: exact blocker exclusion

## Result

H9/14 second-layer survivor 37, with graph6 word

```text
M?????????O?O?o@?
```

cannot be the nerve of a compact hyperbolic Coxeter 9-polytope.  The proof is
a finite exact enumeration using local link (g)-vectors, blocker duality,
and the ridge-degree condition.  It does not trust the archived bounded Z3
timeout.

Consequently the current theorem-facing partition of the 142 graph shadows is

\[
121\text{ independently impossible}
+6\text{ explicit ordinary-polytopal}
+15\text{ undecided}.
\]

The later sixth ordinary-polytopal realization is record 35, and the later
record-34 certificate supplies the additional exclusion. This certificate
excludes one graph shadow; it is not an H9/14 nonexistence theorem.

## The five Dehn--Sommerville schedules

The missing-edge graph is

\[
\bigl\{\{0,11\},\{1,12\},\{0,13\},\{1,13\},\{11,13\}\bigr\}.
\]

If
(a=|\mathcal A_4|), the possible schedules have (a=0,1,2,3,4), facet
counts (102,100,98,96,94), and summed vertex-link (g_4) values
(30,24,18,12,6), respectively.  Strict (g_4\geq1) in all fourteen
vertex links therefore excludes (a=3,4).

There are exactly 462 four-sets that individually satisfy all local incidence
intervals.

- For (a=1), all 462 choices satisfy the family intervals, but every forced
  five-set shadow makes at least one vertex-link (g_4) nonpositive.  The
  shadow sizes are (8,9,10), requiring respectively (3,2,1) further
  minimal five-faces.

- For (a=2), the search checks all
  \(\binom{462}{2}=106{,}491\) pairs; 96,033 satisfy the family intervals.
  Their forced-shadow size distribution is
  
  \[
  15:168,\quad16:21{,}168,\quad17:23{,}268,\quad
  18:33{,}474,\quad19:11{,}340,\quad20:6{,}615.
  \]
  
  The schedule permits total five-face offset 17.  Thus 51,429 families
  already overrun the count, while each of the other 44,604 forces a
  nonpositive vertex-link (g_4).

Hence only (a=0) remains.

## Complete (a=0) census

The (a=0) schedule requires exactly five missing five-faces and 102 facets.
Local intervals allow activity at vertices (0,1,11,12), forbid vertex 13,
and treat vertices (2,\ldots,10) as nine interchangeable isolates.  Every
independent active pair has upper incidence zero, so a missing five-face has
at most one active vertex.  This yields an exact five-row, isolate-column
quotient rather than a labelled brute-force search.

The complete counts are:

| census layer | count |
|---|---:|
| row-normalized isolate-quotient states | 2,865,450 |
| distinct-row states | 2,774,810 |
| states with no transversal of order at most four | 522,427 |
| states with the required 102 minimal five-transversals | 12 |
| distinct five-face families after deduplication | 5 |

Across the 522,427 blocker-valid states, the number of minimal
five-transversals ranges from 40 to 105.  The full histogram, not merely its
target bin, is frozen in the certificate.

For each of the five target families the complete minimal-transversal clutter
is pure of order five and has exactly 102 members, as required.  None is a
closed pseudomanifold:

| number of families | ridge-degree histogram | dual components using degree-two ridges |
|---:|---:|---:|
| 2 | (1:28,\ 2:403,\ 3:28) | 3 |
| 3 | (1:24,\ 2:411,\ 3:24) | 2 |

In particular, every candidate has ridges of degrees one and three, whereas
every ridge in a simplicial polytope boundary has degree two.  This completes
the exclusion of (a=0) and therefore all five schedules.

## Trust and reproducibility

The frozen manifest is
`data/enumerations/h9_14_record37_certificate.json`.

- file SHA-256:
  `cdb7df45b1a886cc1d689c34a2a2ac0cc786afd8edd19c2770eee1ddaacccdad`;
- canonical internal certificate SHA-256:
  `4da8617bd4dfc77889c3c3f987e7988efd87d6d366d7efd29c7403370536cc5c`;
- full (a=0) census SHA-256:
  `8b794a8f4419907acb3cb2e29b6a034403fdc2f8b69d745f1116142f188664ae`;
- target-family audit SHA-256:
  `661aef1cd1cb62f2938f347cfca441a6ff323adc290596ae23f8277e6f3fdf2a`.

The full census rebuilt and verified in about 157 seconds on the current
machine.  A separate replay reconstructs every blocker by scanning all
(2^{14}) masks, rebuilds the packed face tables, and recomputes ridge degrees
and dual components without calling the census routines.

```text
PYTHONPATH=src:. python3 scripts/certify_h9_14_record37.py --verify-only
PYTHONPATH=src:. python3 scripts/certify_h9_14_record37.py --verify-only --recompute
PYTHONPATH=src:. python3 scripts/independent_h9_14_record37_audit.py
PYTHONPATH=src:. pytest -q tests/test_h9_14_record37_certificate.py
```

The certificate makes no assertion about Coxeter labels, Lorentzian Gram
matrices, or the remaining seventeen graph shadows.

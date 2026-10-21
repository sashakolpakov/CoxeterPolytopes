# Record 0, compact schedule `a=2`: stronger audit of the four open sources

## Outcome

The residual source orbits `24, 29, 30, 32` remain genuinely **OPEN**.  No
timeout has been promoted to an exclusion, and no compact-incidence witness
has been promoted to a hyperbolic example.  Consequently the frozen frontier
does not change:

| status | source orbits | count |
|---|---|---:|
| independent two-solver incidence `UNSAT` | `0,8,9,10,12,23,31,34` | 8 |
| explicit incidence-overcover `SAT` | the previously frozen 24 sources | 24 |
| open after the stronger audit | `24,29,30,32` | 4 |

Thus the partition remains `8 + 24 + 4 = 36`.

## Stronger exact formulation

For each fixed pair of missing four-sets, the new primal model selects the
additional independent missing five-sets.  It imposes all of the following
simultaneously:

- the complete Dehn--Sommerville vector
  `(14,89,341,878,1579,1967,1618,792,176)`;
- exact extension of every allowed face of orders `1,...,8` to a selected
  rank-nine facet;
- closure of every allowed set of more than nine vertices by a missing face
  of order at most five;
- simplicial ridge degree `0` or `2`;
- all strict local `g` and Macaulay inequalities for the `14` vertex links,
  `89` edge links, and `341` triple links.

There are `444` local figures and `902` exact two-entry Macaulay tables.  This
strictly contains the earlier compact-incidence scout.  A second composition
fixes the same sources inside the production Z3 facet-selection system, so
the strengthened formulation is not tied to the primal CP-SAT encoding.

## Bounded dispositions

The lossless source-stabilizer lex slice was removed to check whether it was
responsible for the residual hardness.  It was not:

| source | CP-SAT status | observed wall time | branches | conflicts |
|---:|---|---:|---:|---:|
| 24 | `UNKNOWN` | 60.910 s | 0 | 0 |
| 29 | `UNKNOWN` | 52.112 s | 0 | 0 |
| 30 | `UNKNOWN` | 51.889 s | 0 | 0 |
| 32 | `UNKNOWN` | 48.513 s | 0 | 0 |

These runs spent their windows in model presolve.  For source 24, retaining
the lex slice reached 13,978 branches in 60.091 seconds but was still
`UNKNOWN`.  The independently composed production Z3 model for source 24 has
22,687 assertions and returned `unknown` with reason `timeout` after its
30-second check window (following 98.530 seconds of generic construction).

All timings and branch counts are noncanonical run metadata.  The canonical
result is only that no terminal SAT or UNSAT certificate was obtained.

## Positive shadow result for source 29

There is one exact positive result below the facet-completion gate.  An
explicit family of 21 additional A5 sets for source 29 realizes all three
scheduled missing-shadow counts:

| order | allowed ground | missing under the witness | scheduled faces |
|---:|---:|---:|---:|
| 6 | 2,114 | 147 | 1,967 |
| 7 | 2,038 | 420 | 1,618 |
| 8 | 1,421 | 629 | 792 |

It also satisfies every vertex-link incidence interval.  The archive records
the 21 masks and their vertex sets, and the certifier recomputes these counts
directly.

This is **not** a compact-incidence solution.  Its complete blocker has order
histogram

```text
4: 10, 5: 90, 6: 27,
```

so purity and the scheduled 176 facets already fail.  Independently, it has
48 strict-link failures: two edge figures and 46 triple figures.  The frozen
failure-transcript digest is

```text
c87d6941843fca0a43d64039e87e2149431b8bbbf0caa4bf9c8c09e99d618102.
```

The useful conclusion is precise: the global six/seven/eight
Dehn--Sommerville shadow counts alone do **not** kill source 29.  The next
obstruction must use strict links, blocker purity/facet extension, ridge
incidence, or Coxeter/Gram data jointly.

## Subsequent exact localization of source 29

The later
[`triple-012 cut`](h9_14_record0_a2_source29_triple012_cut.md) proves an
analytic strict-link reduction while leaving the status **OPEN**.  The
source stabilizer is `S2 x S2 x S8` of order 161,280.  The link at `012`
deletes one 28-element A5 orbit and confines two further orbits, of total
size 168, to at most one free-triple fibre.  Hence the remaining exact search
splits into one zero-support branch and five normalized active-support
branches, each with 1,404 undecided candidates.

That certificate also replaces the 48-failure staging witness by an exact
scheduled-shadow seed with one strict triple defect of size one.  On the
proved post-cut ground, exhaustive canonical neighbourhood replay finds 596
shadow-exact one-replacement families and 515,354 shadow-exact
two-replacement families, but no strict-link seed.  Thus any repair requires
at least three A5 replacements.  This is a localization, not a source-29
exclusion or a compact/hyperbolic witness.

## Reproduction

- Archive: [`h9_14_record0_a2_open_strong_audit.json`](../data/enumerations/h9_14_record0_a2_open_strong_audit.json)
- Certifier: [`certify_h9_14_record0_a2_open_strong_audit.py`](../scripts/certify_h9_14_record0_a2_open_strong_audit.py)
- Specialized primal model: [`explore_h9_14_record0_a2_open_primal_full.py`](../scripts/explore_h9_14_record0_a2_open_primal_full.py)
- Production Z3 composition: [`explore_h9_14_record0_a2_open_full_links_z3.py`](../scripts/explore_h9_14_record0_a2_open_full_links_z3.py)
- Staged shadow model: [`explore_h9_14_record0_a2_open_shadow_design.py`](../scripts/explore_h9_14_record0_a2_open_shadow_design.py)

## Next exact completion route

1. Solve the 144 pair configurations from the subsequent
   [`triple-123 fibre-budget theorem`](h9_14_record0_a2_source29_triple123_cut.md),
   or search repairs at least three A5 replacements from the one-defect near
   seed.
2. For every surviving A5 family, impose blocker purity, all 176 facets, and
   ridge degree `0/2`.
3. Replay any terminal UNSAT result in an independent SAT/Z3/MILP encoding;
   replay any SAT target by direct face, blocker, ridge, link, Coxeter-label,
   and Gram audits.

No finite-volume noncompact chamber, cusp condition, Coxeter labeling, or
Lorentzian Gram realization is asserted here.

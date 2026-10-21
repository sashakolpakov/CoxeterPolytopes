# Record 31 `a=1`: ridge-count correction and the first incidence survivor

Status: **correction recorded; full `a=1` census still running**.  This note
withdraws the preliminary zero-target count produced by an invalid port of the
`a=2` ridge gate.  It does not claim a hyperbolic Coxeter polytope.

## The caught error

The first generated `a=1` engine correctly changed the facet count from 118
to 120, but accidentally retained the `a=2` expected ridge count 531.  The
preliminary orbit-0 transcript therefore rejected every pure blocker, even
though its complete ridge histogram contained

```text
"2:540": 1.
```

That raw invalid transcript is frozen at
`data/enumerations/h9_14_record31_a1_invalid_ridge531_transcript.json`.  Its
zero target field is explicitly invalid and must not be used as evidence.
All longer searches launched with that binary were stopped, and all their
target counts are withdrawn.

## Why the correct number is 540

The corrected engine does not replace one hard-coded ridge count by another.
It reconstructs the graph-shadow face prefix

```text
f_{-1},f_0,f_1,f_2,f_3 = 1,14,87,319,770,
```

where the unique missing A4 removes one of the 771 independent four-sets.  It
then computes

```text
h_0,...,h_4 = 1,5,11,18,25,
```

applies Dehn--Sommerville symmetry, and inverts the h-vector transform.  The
result is

```text
f = (14,87,319,770,1281,1477,1142,540,120).
```

Thus a ridge pseudomanifold target must have histogram `{2:540}`.  A focused
regression test derives the number independently and forbids the stale
`{2:531}` condition in the generated source.

## The corrected orbit-0 survivor

Both corrected quotient implementations find exactly one target orbit over
the A4 `{2,3,4,5}`.  One representative has the five added A5s

```text
{6,7,8,9,10}
{3,4,7,11,12}
{7,8,9,11,12}
{2,5,6,10,13}
{6,8,9,10,13}.
```

Its literal blocker has histogram `{5:120}` and its ridge histogram is
`{2:540}`.  It passes the exact scheduled face vector, connected dual graph,
every-link Euler characteristic, and every-link mod-two homology audit.

This is a genuine topology-valid abstract incidence survivor.  It is not a
hyperbolic example.  The uncapped compact-label join has stage counts

```text
420 -> 96 -> 8 -> 0.
```

An independent implementation regenerates all 108 order-four and 420
order-five classical Lanner placements, uses the opposite support order, and
again obtains `420 -> 96 -> 8 -> 0`.  Perron--Frobenius monotonicity makes the
zero unrestricted over all omitted finite non-right labels.  Consequently no
Gram state survives and no Lorentzian, discreteness, compactness, or
finite-volume claim arises from this target.

The corrected orbit-0 source and downstream records are frozen at

- `data/enumerations/h9_14_record31_a1_orbit0_corrected.json`;
- `data/enumerations/h9_14_record31_a1_orbit0_downstream.json`; and
- `data/enumerations/h9_14_record31_a1_orbit0_independent_downstream.json`.

The other corrected `a=1` source-orbit searches remain in progress.  Until
they terminate, neither a complete `a=1` no-go nor a complete list of its
incidence survivors is claimed.  Nothing here addresses noncompact
finite-volume chambers or Vinberg endpoint completeness.

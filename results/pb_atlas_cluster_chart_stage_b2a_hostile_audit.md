# Hostile audit: PB-Atlas Stage B2a exact plane-cluster charts

**Verdict (2026-08-10): PASS after strict public-boundary repairs.**  On the
declared coordinate-compiler scope, the proper, free, and satellite chart
substitutions and every tested jet coefficient are exact.  The audit found
and caused repair of several coercive or silently ignored API paths before
the archive was frozen.  It proves no incidence nonemptiness, rank,
irreducibility, hidden-wall, full-nef, or same-polytope statement.

Audited result:
[`pb_atlas_cluster_chart_stage_b2a.md`](pb_atlas_cluster_chart_stage_b2a.md).

## Pinned audited artifacts

| artifact | SHA-256 |
|---|---|
| `data/enumerations/pb_atlas_numerical_cluster_stage_b1.json` | `c645ed67182b765717fd2fad4c1d89c7c97a332ad2d4b10ef31c04d084fa28c1` |
| `scripts/pb_atlas_cluster_chart_stage_b2a.py` | `2b0d534c6e014f3e15e85a200ceaad029160aa33a10ba8b5d6f80cde76904cbf` |
| `scripts/certify_pb_atlas_cluster_chart_stage_b2a.py` | `e288c7cbee8ea7454809da4108c67ae15f433493def7bd8afeb57fac30955213` |
| `schemas/pb-atlas-cluster-chart-stage-b2a.schema.json` | `f6fb32928aeaf82f8b57aed1fd2424c2feb2498b8c7b68ad98d19cd89e3e5be9` |
| `results/pb_atlas_cluster_chart_stage_b2a.md` | `167ef792dc5bf9e2f8f08fe9c76ee46a274d4d10371b6943013900b15f490e51` |
| `data/enumerations/pb_atlas_cluster_chart_stage_b2a.json` | `edf4b93e5023a46a950c006d440c47586453c00530a576235ae3772c5672c2f4` |
| `tests/test_pb_atlas_cluster_chart_stage_b2a.py` | `cbf1b8758e69d1325f6f31fb308601bfe012b9fa4b8a50cfab8d960b3497fe0b` |
| `tests/test_pb_atlas_cluster_chart_stage_b2a_hostile_audit.py` | `68565dc09ff1c97c47caf1a87023392df476a2064289d31031662259b52d9eb3` |

## 1. Direct chart and strict-transform oracle

The hostile suite does not import the production low/high splitter or
strict-transform routine.  For a degree-four monomial it starts with a
nontrivial rational local frame at a proper point, substitutes each displayed
blowup chart into the monomial itself, and reads the required coefficient
directly.

For the free (A)-chart and its forced satellite (B)-chart, the composite
map is

\[
 (u_0,v_0)=
 (u_2v_2, u_2v_2(t+v_2)).                         \tag{1.1}
\]

For the free (B)-chart and its forced satellite (A)-chart, it is

\[
 (u_0,v_0)=
 (u_2v_2(t+u_2), u_2v_2).                         \tag{1.2}
\]

If the parent and free-center multiplicities are (m_0,m_1), division in
(1.1) shifts the satellite coefficient by

\[
 (a,b)\longmapsto(m_0+a, m_0+m_1+b),
\]

whereas (1.2) shifts it by

\[
 (a,b)\longmapsto(m_0+m_1+a, m_0+b).
\]

The independently extracted rows agree entry by entry with Stage B2a for
both orientations, not merely up to rank or row span.  The same test checks
the proper-point Taylor rows and the exact parent-axis labels at the
satellite.

There is also a coefficientwise reason no local Laurent denominator can
occur.  In chart (A), a parent monomial (u^av^b) becomes

\[
 u^{a+b}(t+v)^b;
\]

in chart (B), it becomes

\[
 v^{a+b}(t+u)^a.
\]

After terms of total degree below (m) are removed, every retained term is
divisible by the selected exceptional coordinate to at least the (m)-th
power.  The hostile rows additionally contain none of the local symbols
(u_i,v_i) after coefficient extraction.

## 2. Commuting centers and chart options

The two frozen free-forest fixtures encode the same proximity DAG with two
independent proper roots in opposite chronological orders.  The audit uses
the asymmetric multiplicity row

```text
[2, 1, 2, 1]
```

so an accidental symmetric comparison cannot pass.  After transporting the
center bijection, multiplicities, jet descriptors, and the (q)-variables,
the two degree-four matrices agree coefficientwise.

Every declared option was also checked for actual dataflow.  Ungauged proper
patches, rational local frames, both free-chart modes, explicit gauge
quadruples, automatic gauges, and disabled gauges all change or annotate the
claimed object as specified.  A proper-patch choice on a gauged center was
initially accepted and silently ignored; it is now rejected.

The four determinant conditions attached to a gauge are exactly the four
three-by-three minors of its selected projective frame.  Alternate proper
quadruples, three affine patches, two free blowup charts, and the explicit
ungauged mode give the stated finite coordinate-cover boundary.  This is not
promoted to a claim that an incidence locus meets the default chart.
Sibling collisions beyond the factors explicitly emitted remain part of
Stage B2b, as the primary report says.

## 3. Independent H10 and H11 calibration

The hostile calibration imports neither
`_cluster_condition_matrix` from the H10 certificate nor
`_quintic_condition_matrix` from the H11 certificate.  The old archives are
used only to match the contracted raw-wall set to its historical point labels.
The expected rows are then rebuilt from literal point, frame, and slope data.

For H10, homogeneous monomials are evaluated at the seven proper points and
pulled through the three displayed near-point maps.  The keyed center/jet
rows reproduce both the (10\times10) cubic and (30\times28) sextic
matrices exactly.

For H11, the four root charts are independently taken as

\[
 (u,v),\quad(1+u,v),\quad(v,1+u),\quad(-1+u,2-2u+v).
\]

Successive (A)-chart pullbacks are composed directly, with second slopes
(0,a,b) on the three length-three chains.  Shifted monomial coefficients
reproduce the complete (23\times21) quintic matrix over
(\mathbf Q(a,b)), keyed by center and multiindex.  Thus the hostile check
cannot pass because the new and old certifiers happen to share a matrix
helper or a row-multiset convention.

These remain entry calibrations only.  No legacy rank, kernel, curve, or
surface conclusion is re-imported.

## 4. Public-boundary defects repaired before freeze

The first implementation silently converted some inputs with `int(...)`.
It therefore accepted a degree of `3.75` as degree three, `True` as degree
one, degree `-1` as an empty-column system, and fractional exponent rows as
truncated integer monomials.  A duplicate satellite parent reached a bare
`StopIteration`.

The second hostile pass found normalized option-key aliases such as `"0"`
and `"00"`, scalar frame rows that leaked `TypeError`, non-string unknown
keys that broke error formatting, approximate SymPy floats admitted as huge
rationals, and the ignored gauged-center patch choice.  The repaired boundary
now:

1. requires exact nonnegative Python integers for degrees, multiplicities,
   and monomial exponents;
2. requires a nonempty list of distinct homogeneous monomial columns;
3. requires distinct, chronological proximity parents and rejects malformed
   satellite data before chart construction;
4. accepts only canonical ASCII decimal option keys and rejects collisions
   after normalization;
5. accepts exact rational frames but rejects Boolean or approximate entries;
6. returns clean `ValueError` failures for malformed packet, option, frame,
   multiplicity, and column containers; and
7. rejects any proper-patch override that would otherwise be ignored by a
   selected PGL gauge.

The frozen certifier replays these hostile probes in addition to the owner
and independent hostile tests.

## 5. Exact boundary of the PASS

The archive covers all 42 frozen plane representatives, representing 630
labelled endpoints and 433 center occurrences.  The current real packets
have no satellite centers; the separate two-parent fixture is therefore
essential evidence rather than a census statistic.

Stage B2a supplies exact coordinate and evaluation data.  It does not prove:

- that any selected chart-open or incidence stratum is nonempty;
- a rank or independence statement for any evaluation matrix;
- collision, Jacobian, irreducibility, singular-scheme, or Hilbert-factor
  completeness;
- a reduced strong-Coble boundary;
- hidden-wall exclusion or full-Mori/full-nef generation;
- literal same-(P) return or L4.

Those are the Stage-B2b geometric gates.  In particular, this PASS is not a
dimension-bound existence theorem by itself.

## Replay

```bash
/Users/sasha/.venvs/global/bin/python \
  scripts/certify_pb_atlas_cluster_chart_stage_b2a.py

/Users/sasha/.venvs/global/bin/python -m pytest -q \
  tests/test_pb_atlas_cluster_chart_stage_b2a.py \
  tests/test_pb_atlas_cluster_chart_stage_b2a_hostile_audit.py
```

The final combined count is 15 passing tests: eight owner tests and seven
independent hostile tests.

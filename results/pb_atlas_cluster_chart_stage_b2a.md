# PB-Atlas Stage B2a: exact cluster charts for the \(\mathbf P^2\) packets

**Status (2026-08-10): PASS on the declared Stage-B2a scope.**  The compiler
turns every frozen Stage-B1 \(\mathbf P^2\) endpoint packet into one explicit
ordered-cluster chart and replaces all formal \(-K\) and \(-2K\) jet tokens by
exact rational-polynomial matrix entries.  It is a coordinate compiler, not an
incidence-existence or same-polytope theorem.

## Frozen interface

The production entry points are

```python
compile_p2_cluster_chart(stage_B1_packet, chart_options=None)
compile_p2_divisor_evaluation_matrix(
    stage_B1_packet,
    degree=d,
    multiplicities=(m_1, ..., m_r),
    chart_options=None,
    section_columns=None,
)
```

They consume the frozen
`coxeter.pb_atlas.numerical_cluster_stage_b1.v1` contract and accept only its
\(\mathbf P^2\) endpoint packets.  The compiler reads the chronological center
records, proximity parents, and indexed section/jet arrays.  It does not read a
row name, surface type, calibration alias, Coxeter type, or historical point
label.  Ruled \(F_e\) packets are rejected at this stage rather than silently
fed through plane formulas.

The general divisor entry point is genuinely row agnostic: \(d\), the complete
chronological multiplicity vector, and (optionally) the degree-\(d\) monomial
columns determine the matrix.  The special entry point calls that same engine
for the two Stage-B1 rows.

The public boundary is strict.  Boolean, fractional, approximate, or negative
degrees/exponents are not coerced; monomial columns must be nonempty and
distinct; proximity parents must be distinct; option-map indices have one
canonical nonnegative decimal spelling; and local frames accept exact rational
entries but reject floating approximations.  These controls prevent two
different hostile inputs from collapsing silently to the same alleged exact
chart.  A proper-patch override aimed at a \(\mathrm{PGL}_3\)-gauged center is
also rejected, since the canonical gauge target already fixes its affine patch.

## The cluster chart

For a proper center \(p_i\), select a homogeneous coordinate that is nonzero,
normalize it to one, and write the other two coordinates as
\((q_{i0},q_{i1})\).  An optional invertible rational \(2\times2\) frame fixes
local coordinates \((u_i,v_i)\).  Thus the default \(Z\ne0\) chart is

\[
 p_i=(q_{i0}:q_{i1}:1),\qquad
 (X,Y,Z)=(q_{i0}+u_i,q_{i1}+v_i,1).
\]

For a free center over \(p_j\), the two standard blowup charts are

\[
 \begin{array}{ll}
 A:&u_j=u_i,\quad v_j=u_i(t_i+v_i),\\[2mm]
 B:&u_j=v_i(t_i+u_i),\quad v_j=v_i.
 \end{array}
\]

The parameter \(t_i\) is the selected point of the exceptional line in that
chart.  For a satellite center, its two proximity parents must appear as the
two exceptional axes in a parent chart.  The required blowup chart is then
forced and the intersection direction is \(t_i=0\).  The output records which
parent divisor is the \(u_i\)-axis and which is the \(v_i\)-axis.  The hostile
three-center fixture `proper -> free -> satellite` certifies this transport
without relying on the current atlas census, whose 42 representative plane
packets happen to contain no satellite center.

When four proper centers exist, the default chart applies the deterministic
\(\mathrm{PGL}_3\) gauge to the first four chronological proper centers:

\[
 (1:0:0),\quad(0:1:0),\quad(0:0:1),\quad(1:1:1).
\]

This is valid only where the four corresponding triple determinants are
nonzero; those factors are emitted explicitly.  The API also accepts an
alternate proper quadruple or the empty gauge.  The empty-gauge mode covers
configurations with no general-position quadruple.  Affine patch choices and
the two blowup charts give a finite coordinate cover, not a proof that the
desired incidence locus meets the selected default member of that cover.

## Exact jet transport

Let \(f_i(u_i,v_i)\) be the current strict transform of a section at center
\(p_i\), and let its required multiplicity be \(m_i\).  Stage B2a performs the
following exact operation:

1. expand \(f_i\) in \(u_i,v_i\);
2. emit the coefficient of every \(u_i^a v_i^b\) with \(a+b<m_i\), in the
   frozen Stage-B1 total-order/multiindex order;
3. delete those low terms;
4. substitute the chosen blowup chart and divide the remaining expression by
   the exceptional coordinate to the power \(m_i\);
5. pass that strict transform to every child for which \(p_i\) is the selected
   primary parent.

Because the low terms are removed before division, no illicit local Laurent
term is introduced.  Every serialized entry is an exact numerator/denominator
pair over \(\mathbf Q\) in the selected \(q\)- and \(t\)-variables.  The
normalization is the Taylor coefficient, rather than a factorial-scaled
derivative.  The \(-K\) multiplicity-one and \(-2K\) multiplicity-two rows are
therefore unchanged under the legacy derivative convention.

## Census and independent plane calibration

The frozen replay covers all **42** Stage-A \(\mathbf P^2\) orbit
representatives, representing **630** labelled endpoints and **433** center
occurrences.  Every Stage-B1 \(-K\) and \(-2K\) formal matrix has the same row
descriptors, monomial columns, and dimensions after exact replacement.

The certifier selects the two pre-existing plane inverse examples without a
type/name route:

* The unique proximity forest with component sizes
  \([2,2,2,1,1,1,1]\) is matched to the old plane branch by equality of its
  contracted-wall set.  After specializing the universal coordinates to the
  old seven proper and three near points, the new \(10\times10\) cubic and
  \(30\times28\) sextic matrices equal the old matrices row-for-row as
  multisets.
* The unique proximity forest with component sizes \([3,3,3,2]\) is likewise
  matched by its contracted-wall set.  The generic divisor API—not a
  \(-K\)/\(-2K\) special case—uses the frozen degree-five boundary row and its
  chronological multiplicities.  Its specialized \(23\times21\) matrix equals
  the old quintic matrix over \(\mathbf Q(a,b)\), again as an exact row
  multiset.

These are **jet-entry equalities**.  They do not re-import, or newly assert,
the legacy rank, kernel, nonemptiness, irreducibility, or surface conclusions.

## Exact boundary

Stage B2a does **not** prove any of the following:

* that a declared chart-open set or incidence locus is nonempty;
* independence or rank of the jet conditions;
* completeness of collision, Jacobian, reducibility, singular-scheme, or
  Hilbert bad factors;
* existence of a strong-Coble boundary;
* hidden-wall or full-Mori/full-nef completeness;
* literal same-\(P\) return, or L4.

Those are Stage-B2b geometric tasks.  In particular, an exact matrix on a
universal coordinate chart is data for the incidence problem, not a solution
of it.

There is a useful downstream separation: after Stage B2b has independently
proved the strong-Coble and effective-facet gates, the project-native
reflection-maximal exact-return lemma in
`results/reflection_maximal_coble_exact_nef_return.md` can remove a separate
hidden-wall/full-Mori enumeration.  Stage B2a proves neither gate and therefore
does not invoke that lemma as part of this PASS.

## Artifacts and replay

* Core: `scripts/pb_atlas_cluster_chart_stage_b2a.py`
* Certifier: `scripts/certify_pb_atlas_cluster_chart_stage_b2a.py`
* Packet schema: `schemas/pb-atlas-cluster-chart-stage-b2a.schema.json`
* Frozen compact census: `data/enumerations/pb_atlas_cluster_chart_stage_b2a.json`
* Tests: `tests/test_pb_atlas_cluster_chart_stage_b2a.py`

Replay with the stable project interpreter:

```bash
/Users/sasha/.venvs/global/bin/python scripts/certify_pb_atlas_cluster_chart_stage_b2a.py
/Users/sasha/.venvs/global/bin/python -m pytest -q tests/test_pb_atlas_cluster_chart_stage_b2a.py
```

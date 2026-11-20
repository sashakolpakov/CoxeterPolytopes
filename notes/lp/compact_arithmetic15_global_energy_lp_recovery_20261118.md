# Arithmetic 15 from the proved local inputs: a small energy LP

This is a source-only recovery model, prepared on 18 November 2026. It
retains the existing proof's local theorems and complete necessary covers.
It is not the unfinished programme deriving those theorems from primitive
Gram conditions. No worker, solver, constructor or scientific body was run
to prepare it. The proposed LP and its hand dual still require independent
review and root acceptance; this note credits no new execution or bound.

The fixed mathematical source is
[`compact_arithmetic_dimension_bound_consolidated.tex`](../../manuscripts/revisions/compact_arithmetic_dimension_bound_consolidated.tex),
SHA-256 `d0f040b253146f318bafd2391c8809c9d995f53fcd920477f78f42038d424295`.
The complete covers are `prop:f17` / Table `tab:f17` (lines 467--500) and
`prop:f16` / Table `tab:f16` (lines 1000--1030). Their terminal consumers are
`app:arithrankseventeenclosure` (lines 20197--20316) and
`app:arithranksixteenclosure` (lines 20320--20476).

## What is retained as a proved input

Fix one actual compact arithmetic polytope, an original farthest vertex
and the source's one spatial vector x, with squared norm less than one.
For each original orthogonal component j let E_j be its projection energy.
Thus E_j is nonnegative and the sum of all E_j is less than one.

The initial recovery retains these named input families:

* **GEO:** the original-vector geometry, actual opposite-wall meanings,
  orthogonal component partition and strict energy budget.
* **COST:** the proved component, joint-component and fully paid local
  reference inequalities below, with their actual rank, spectator and
  arithmetic hypotheses. A joint inequality owns every component in its
  support once; its terms are not independent extra energy slots.
* **COVER16/COVER17:** the two complete necessary covers, including all
  original defect and field hypotheses and all baseline theorems used to
  prune them. These are proved classification outputs, not unexplained
  physical axioms. Ordinary unbounded dihedral labels remain symbolic.
* **EXIT:** the exact local exclusion theorems listed below. They retain
  actual domain inclusions, null/point/acuteness arguments, field descent,
  integral cycles, all nonphysical embeddings and full field-square tests
  used in those theorems. No local exclusion is inferred from a source name.
* **DATA1--DATA4:** the four computed-output premises of
  `prop:pair-output`, at their complete proved input domains: quadratic A17
  emptiness; the two cubic A17 finite-anchor empty domains; the shared
  quadratic B/D17 two-survivor output; and the whole-B16 one-survivor output.
  B17 and D17 retain distinct actual-domain inclusion theorems.
* **LATTICE:** wherever a retained local proof uses the golden-lattice
  theorem, retain its literal invariant lattice and full reflection-group
  containment bridge. A shared golden field alone is insufficient.
* **COMPACT18--29:** the same fixed compact branch/cost catalogue in each
  actual rank 18--29, valid also for arithmetic polytopes because they are
  compact. The compact extraction is prepared separately. Append those
  branches unchanged to this rank-16/17 catalogue; do not replace them by
  an unexplained dimension-at-most-17 row.
* **CEILING29:** the explicit published general compact ceiling 29,
  separately excluding all dimensions at least 30. The arithmetic master
  does not use a site-6 coefficient construction to establish this input.

These inputs are already mathematical parts of the existing proof. Their
new low-level typed transcription is not a prerequisite for this initial
proof recovery. In particular no list of 454 semantic keys or sixty prefix
consumers is an input to this master.

## Local costs and notation

All printed decimals below are exact rational numbers. Every strict source
floor is weakened to a closed lower bound in the LP; strictness is retained
in the geometric upper budget instead. No positive epsilon is introduced.

T2 means an ordinary irreducible I2(m), m at least 3; T3 means A3/B3/H3;
T4 means A4/B4/D4/F4/H4. I, J, N7 and O7 have the actual-wall meanings
in source lines 151--176, not just old diagram types. In particular I owns
one rank-one component, J owns two, and N7/O7 own a rank-two component.
An unmarked seven bond is not an O7 defect. Unmarked small components in
the cover tables are ordinary. Every alternative in a table cell remains.

The following reusable cost rows are proved before the terminal dispatch:

| Cost key | Support and closed LP floor | Source anchor |
|---|---|---|
| O1 | E(A1) >= 1/10 | `subsec:shared-seals`, lines 174--197 |
| O2 | E(T2) >= 1/7; ordinary T2 >= 3/19 when needed, including A2 by its separate 319/2000 floor | shared seals; `prop:ordinarydihedralthreenineteenths`, lines 2867--3070 |
| O3/O4/O5 | E(T3) >= 1/4; E(T4) >= 2/7; each retained ordinary rank-five forest >= 5/14 | shared seals and the explicit cover partitions |
| DEF | E(N7), E(O7) >= 9/100 | source lines 164--176 |
| PROP | Beside an actual typed defect, E(A1) >= 2/15 and E(T2) >= 1/5; use only at its proved propagated domain | cover proofs and `app:arithcoverinputs` |
| I2/I3/I4 | E(I+T2) >= 21/100; E(I+T3) >= 2/7; E(I+T4) >= 23/70 | `lem:Ijoint`, lines 212--267 |
| I1 | E(I+A1) >= 5051/25000 (= .20204) | `lem:Ijoint` |
| A10d | E(A10) >= 404/625 beside actual I or either actual typed defect | `app:arithA10consumers`; final dispatch lines 20335--20342 |
| A11 | E(A11) >= 13/20 | `app:arithA11uniform` |
| A11I | E(A11+I) >= 18/25, jointly | `app:arithA11Ijoint` |
| A12/A13/A14 | E(A12) >= 5/7; E(A13) >= 53/70; E(A14) >= 4/5 | `app:arithA12universal`, `app:arithA13universal`, `app:arithA14fourfifths` |
| A15/A16 | E(A15) >= .8348; E(A16) >= .8703 | `app:arithsummedfeet`, `app:arithgoldenend`; source lines 20065--20099 |
| A3O | E(A3) >= 13/50 in the actual O7+A3 domain | `lem:A3point26`; final dispatch lines 20351--20356 |
| BDtyped | E(B14), E(D14), E(D15) >= .86989+.0487; E(B15) >= .88989+.0487, each after its single paid head construction | `app:arithBD14typedhead`, `app:arithBD15typed` |
| D13 | E(D13) >= .82436 in its full rank-three-spectator domain; the retained remainder jointly costs at least .21 | `app:arithD13rankthree` |
| D14ord | E(D14) >= .84816 beside its actual ordinary rank-two spectator | `cor:D14ordinarydihedral`, lines 3068--3088 |

The costs in BDtyped are one component inequality each. The head gain is
not a second disjoint energy component. A13, A14 and the ordinary dihedral
costs may be used together only in the source-matched ordinary domains;
typed K7 conditions are not transferred to ordinary T2.

## Complete rank-seventeen recipes

Each row below is exactly one group of Table `tab:f17`. A semicolon splits
different retained alternatives, not a choice that may be dropped. An
**X** entry is an explicit named local exclusion input and imposes zero
mass only on its own complete matched branch. An **E** entry lists costs
on disjoint supports, except where a displayed joint support is named.

| Group | Complete component and spectators | Master recipe / conservative total |
|---|---|---|
| 17.1 | A11 + I + (A5 or D5) | E: 13/20 + 5/14 = 141/140; ignore the nonnegative I slot |
| 17.2 | A11 + (N7 or O7) + (A4/B4/D4/F4) | E: 13/20 + 9/100 + 2/7 |
| 17.3 | A12 + A5/D5/(A4 or H4)+A1/T3+T2/2T2+A1 | E: 5/7 + 5/14 = 15/14 |
| 17.4 | A12 + I+T4 / I+T3+A1 / I+2T2 | E: 5/7 plus, respectively, 23/70; 2/7+1/10; 21/100+1/7. I appears in one joint support only |
| 17.5 | A14 + T3 / T2+A1 / I+T2 | E: 4/5 plus, respectively, 1/4; 1/7+1/10; 21/100 |
| 17.6 | A15 + T2 / 2A1 / I+A1 / N7 / O7 | T2: X-A15T2 below; 2A1: E .8348+.2; I+A1: E .8348+.20204; typed: E .9191408+.09 after `app:arithA15typed` |
| 17.7 | A16 + ordinary A1 or actual I | X-A16R1 |
| 17.8 | whole A17 | X-A17, with DATA1/DATA2 and their actual endpoint/domain inclusions |
| 17.9 | B15 + N7 or O7 | E: .88989+.0487+.09 = 1.02859 |
| 17.10 | B16 + actual I | X-B16I |
| 17.11 | whole B17 | E: 1.00184 on the whole original support after its reduction, DATA3 and B17-specific survivor inclusion/consumer |
| 17.12 | D15 + N7 or O7 | E: .86989+.0487+.09 = 1.00859 |
| 17.13 | D16 + ordinary A1 or actual I | ordinary: X-D16R1; actual I: X-D16I |
| 17.14 | whole D17 | E: 1.00184 after its separate genuine-fork reduction, DATA3 and D17-specific survivor inclusion/consumer |

17.3's slash-separated spectator alternatives have rank five; parentheses
retain the original grouping. The complete unambiguous list is
`{A5,D5,(A4+A1),(H4+A1),(T3+T2),(2T2+A1)}`. The same convention in
the remaining table cells denotes unions of whole spectator forests.

## Complete rank-sixteen recipes

| Group | Complete component and spectators | Master recipe / conservative total |
|---|---|---|
| 16.1 | A10 + I + (A5 or D5) | E: 404/625 + 5/14 |
| 16.2 | A10 + (N7 or O7) + (A4/B4/D4/F4) | E: 404/625 + 9/100 + 2/7 |
| 16.3 | A11 + A5/D5/(A4 or H4)+A1/A3+T2/2T2+A1 | E: 13/20 + 5/14 = 141/140; same explicit rank-five forest list as 17.3, with A3 in place of T3 |
| 16.4 | A11 + I + T4 | E: 18/25 on A11+I jointly, plus 2/7 on T4 = 176/175 |
| 16.5 | A11 + O7 + A3 | E: 13/20 + 9/100 + 13/50 = 1 |
| 16.6 | A12 + A4 / H4 / 2T2 | E: 5/7 + 2/7 = 1 |
| 16.7 | A12 + I + T3 | E: 5/7 + 2/7 on I+T3 jointly = 1 |
| 16.8 | A12 + (N7 or O7) + T2 | E: 5/7 + 9/100 + 1/5, using the actual typed propagation |
| 16.9 | A13 + A3 / T2+A1 / I+T2 / O7+A1 | A3: E 53/70+1/4; T2+A1: E 53/70+1/7+1/10=1; I+T2: X-A13I; O7+A1: E 233/300+9/100+2/15=1 after the complete K7 strip proof |
| 16.10 | A14 + T2 / 2A1 / I+A1 / N7 / O7 | T2: E 16/19+3/19=1 after its complete ordinary endpoint/foot proof; 2A1: E 4/5+.2=1; I+A1: E 4/5+.20204; typed: E .91286445198+.09 |
| 16.11 | A15 + ordinary A1 or actual I | X-A15R1 |
| 16.12 | whole A16 | X-A16 |
| 16.13 | B14 + N7 or O7 | E: .86989+.0487+.09 = 1.00859 |
| 16.14 | B15 + ordinary A1 or actual I | X-BD15R1 for ordinary; X-B15I for actual I |
| 16.15 | whole B16 | E: 1.00640742675 after its actual reduction, DATA4 and unique-survivor edge consumer |
| 16.16 | D13 + I+T2 / O7+A1 | E: .82436+.21 = 1.03436, with the remainder owned jointly |
| 16.17 | D14 + ordinary I2(m), m>=4 / N7 / O7 | ordinary: E .84816+3/19; typed: E .86989+.0487+.09 |
| 16.18 | D15 + ordinary A1 or actual I | X-BD15R1 for ordinary; X-D15I for actual I |
| 16.19 | whole D16 | X-D16DYADIC; no whole-D16 census premise |

Every E total is at least one. Equality in 16.5--7, 16.9 and 16.10 is
sufficient: actual total energy is strictly less than one. The LP must not
ask for a uniform positive margin or silently strengthen any local floor.

## Exact local exit interfaces

These are intentionally explicit imports from the existing proof. They are
not conclusions supplied by the coefficient LP or by generic real PSD.
Importing a complete local consumer here avoids expanding its already
proved continuous wall/field subbranches into thousands of scalar rows.
Each exit preserves the complete domain and actual original-vector proof.

| Exit key | Matched branch and complete source theorem |
|---|---|
| X-A15T2 | Whole A15 plus ordinary T2, rank 17; `app:arithA15ordinary`, lines 15066--15458. It splits high T2 energy, silver and golden domains. The displayed total floors include .8348+1/6, 1.000094738, and 1+123/95000. The low-energy native reduction and arithmetic exclusions remain part of this input; there is no blanket T2>1/6 assumption |
| X-A16R1 | Whole A16 plus the same actual rank-one opposite, rank 17; `app:arithA16rankone`, lines 15462--15848. Full nonisolated/isolated/endpoint alternatives end in norm or distinct-wall acuteness contradictions; do not replace this by an invented component floor |
| X-A17 | Whole A17, rank 17; the local `app:wholearithA17`, `app:arithA17datainclusion` and `prop:pair-output` (i),(ii). Actual endpoint isolation, full degree/field reduction and all computed-domain inclusions are retained. The global exact-rank exclusion `cor:no17` is not an input |
| X-B16I | Whole B16 plus actual I, rank 17; `app:arithB16Ishort` through `app:arithB16Ihead`, lines 16031--16954. Full six-role, cubic and quadratic cover. The twelve-row final quadratic table has minimum displayed total 1.00026; rejected field/cycle/square branches remain in the complete local exclusion |
| X-D16R1 | Whole D16 plus ordinary/nonshort A1, rank 17; `app:arithD16ordinary` through `app:arithD16ordinaryfinish`, lines 16957--17530. The ordinary bootstrap and all fork branches precede the final .9084+.1>1; it is not a premise for actual short I |
| X-D16I | Whole D16 plus actual I, rank 17; `app:arithD16Ifoundation` through `app:arithD16Iquadraticfinish`, lines 17536--20044. Full degree-2/3/4 actual-root, both-fork and field-square/dyadic arguments, including unequal pairs. No new computed-output premise |
| X-A13I | Whole A13+I+T2, rank 16; `app:arithA13I`, lines 11952--12133. Both endpoint residuals and every remaining foot are paid jointly. The positive scalar totals are displayed, but no independent A13 floor is substituted for this joint proof |
| X-A15R1 | Whole A15 plus its same actual rank-one opposite, rank 16; `app:arithrankone` and its A15 subcases before the B/D15 ordinary section. Full isolated/nonisolated split, including every h>0 and boundaries |
| X-A16 | Whole A16, rank 16; `app:wholearithA16` through the full actual finite-anchor/time-ratio proof ending at lines 4760--4799. This is a complete manual arithmetic exclusion, not DATA1 or an A17 theorem |
| X-BD15R1 | B15 or D15 plus ordinary/nonshort A1, rank 16; their separate nonisolated reductions and `app:arithBD15ordinary`, lines 5816--6146. The isolated h>=1/2 domain includes equality; actual short I is excluded from this interface |
| X-B15I | Whole B15+actual I, rank 16; `app:arithB15Iwindow` through `app:arithB15Ipurepair`, lines 6150--7141, full 0<h<1/2 domain |
| X-D15I | Whole D15+actual I, rank 16; `app:arithD15Icentre` through `app:arithD15Iforknonisolated`, lines 7144--7890, full 0<h<1/2 domain and both signed forks |
| X-D16DYADIC | Whole D16, rank 16; `prop:d16reduce`, `prop:noD16`, section `sec:d16`: actual quadratic reduction and all odd/inert/ramified/split dyadic cases. This is distinct from D16+I and from every finite-output attachment |

For whole B17/D17, DATA3 is not itself an exit: both necessary survivors
remain until their distinct original-vector consumers supply the numerical
1.00184 floor (`app:arithB17domain`, lines 10041--10123;
`app:arithD17domain`, lines 11083--11164). Likewise DATA4 is not emptiness:
the whole-B16 survivor remains until its 1.00640742675 edge consumer
(`app:arithB16domain`, lines 9400 onward). No surviving algebraic model is
declared a realized polytope.

## The actual master and an explicit dual certificate

Split a table group only when its displayed alternatives use different
recipes. Let b index these finitely many recipes. A symbolic family such as
ordinary I2(m), m>=3, is one recipe only because its imported theorem is
uniform on the whole stated family. Each recipe retains its full source
domain; it is not replaced by an unconstrained printed number.

For each b introduce a mass x_b>=0 and one energy e_bj>=0 for every original
component j of that forest. A joint I+T or A11+I row uses the corresponding
set of these original slots, not a new independent energy. Write its
proved lower bound as c_bp and its support as S_bp. The LP is

    minimize       sum_(b,j) e_bj
    subject to     sum_b x_b = 1,
                   sum_(j in S_bp) e_bj >= c_bp x_b  (every retained cost),
                   x_b = 0                         (every matched X exit),
                   x_b,e_bj >= 0.

Use either both rank tables together or one copy at each exact rank. An
actual polytope maps to one recipe with mass one, its actual component
energies in that recipe and zeros elsewhere. Its objective is less than
one. A cover need not be disjoint: use any matching recipe, or fix a
first-applicable order. Fractional mixtures are permitted; no integrality
conclusion is drawn from them.

The dual can be written with a common eta, nonnegative packing weights
alpha_bp, and free nu_b only on X branches:

    maximize       eta
    subject to     sum_(p:j in S_bp) alpha_bp <= 1       (every energy slot),
                   sum_p c_bp alpha_bp + nu_b >= eta   (every b),
                   alpha_bp >= 0; nu_b=0 on non-X b.

The table gives a rational dual with eta=1: on an E recipe take weight one
on its displayed disjoint cost rows and zero on other retained costs. The
slot capacities are at most one and the displayed total is at least one.
On an X recipe take all cost weights zero and nu_b=1. Thus every feasible
master point has total energy at least one, contradicting its actual
strict budget. If extra overlapping valid costs are retained, their dual
weights can remain zero; they are never added for free.

This is a complete hand LP certificate for the stated imported inputs,
not a claim that an exported matrix or solver run already exists. Its
numbers are explicit, and it needs no numerical optimizer to discover the
dual. An eventual export/replay should check this exact assignment and all
table coverage/support records, without rerunning the successful local
scientific jobs merely to rename them.

## Exact scope of the recovery

This initial LP starts from the **proved pruned terminal covers**. It does
not claim to derive those covers from all spherical forests. Their baseline
and arithmetic proofs remain named inputs, so no changing cost can silently
invalidate the pruning. The next faithful enlargement could start with all
spherical components and use the same valid cost families to derive the
pruning by rank/defect-state flow. That enlargement is not necessary to
reproduce the existing proof at this first explicit scope and is not built
here.

Every branch retains its actual dimension. Append the compact ranks
18--29, then at threshold N retain exactly the branches with dimension at
least N. This is one fixed catalogue with nested threshold feasible sets.
An eta=1 dual at N=16 excludes every admitted rank 16--29; CEILING29 handles
all larger dimensions. Thus the joined recovery does not assume the old
global compact17 theorem as an additional row. Before that join is built,
the original independent compact17 theorem is an alternative explicit
tail dependency with the same one-way logical scope.

This recovers the existing arithmetic 15 conclusion conditional on the
retained source inputs. There is no downward face inheritance or use of
an arithmetic rank-16 theorem to prove the compact branch costs.

Changing the local theorems, weakening their domains or deleting the retained
cover/field/data/lattice inputs is a stronger research problem. The present
LP makes their location explicit rather than claiming that generic Gram
positivity or a solved site-6 coefficient package already proves them.

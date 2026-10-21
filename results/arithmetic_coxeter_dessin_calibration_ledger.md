# Arithmetic Coxeter--dessin calibration ledger

**Status (2026-08-10): exact six-row calibration plus three determinant-line
controls and the universal Gram--Belyi coordinate.**  Every numerical entry below is replayed from a SHA-pinned frozen
artifact.  An unperformed determinant, Heegner, or vector-valued modular-form
calculation is recorded as `open`, never inferred from the existence of a
surface.

The ledger is
[`arithmetic_coxeter_dessin_calibration_ledger.json`](../data/enumerations/arithmetic_coxeter_dessin_calibration_ledger.json),
built by
[`certify_arithmetic_coxeter_dessin_calibration_ledger.py`](../scripts/certify_arithmetic_coxeter_dessin_calibration_ledger.py).

## 1. Exact matrix

Here `F/I/H` means finite vertices, ideal vertices/cusps, and hyperideal
extreme rays.  The arithmetic definition field and the ground field of a
native surface carrier are deliberately separate columns.

| input | exact polytope | arithmetic lane | exact `F/I/H` | carrier direction and status | determinant/Heegner status | honest missing gate |
|---|---:|---|---:|---|---|---|
| Bugaenko | `H8/11`, compact | `Q(sqrt(5))`, totally real compact | `41/0/0` | `P ->` determinant-line arrangement and a noncanonical complete-intersection surface; separately, an existential K3 eigen-nef `surface -> P` carrier | row-specific extension and full norm-two shell proved; `11` facet mirrors, `309` lower-boundary contacts, no interior cut | explicit least safe automorphic correction or a Hilbert-reflective classification with uniform degree bounds |
| toric endpoint | `H6/9`, finite-volume noncompact | `Q`, rational cuspidal | `18/3/0` | automatic smooth-toric `surface -> P` and finite twelve-row `P -> surface -> same P` recognizer; scoped `L3`, not `L4` | not compiled | determinant support and a broad P-only carrier denominator beyond the toric-completable locus |
| Coble Type II | `H10/30`, finite-volume noncompact | arithmetic field `Q`; equations over `Q`, nef return over `C` | `1333/207/0` | strict raw `P ->` B3 plane cluster `->` strong Coble `-> same P`; single-input `L3`, not `L4` | not compiled | make the contraction/incidence/Hilbert grammar exhaustive on a variable P-only class |
| Coble Type I | `H11/26`, finite-volume noncompact | arithmetic field `Q`; a `Q` line--quintic model of a complex Coble surface | `508/42/0` | repaired raw `P ->` line--quintic determinantal locus `->` all 26 curves `-> same P` | not compiled | promote the one-input grammar to a bounded variable family |
| Coble Type VI | `H14/60`, finite-volume noncompact | arithmetic field `Q`; native carrier in characteristic `3` | `27845/1027/0` | raw P recovers the five-boundary package; five plane blowdowns force the conic--four-tangent characteristic-three point; individual-chamber `L3` exact return | not compiled | replace the singled-out template and Kondō exhaustion input by an intrinsic P-only carrier class |
| Coble Type MII | `H17/48`, finite-volume noncompact | arithmetic field `Q`; native carrier in characteristic `3` | `190808/5244/0` | raw P recovers the unique eight-boundary frame, `P1 x P1` grid, and all 24 graph curves; individual-chamber `L3` exact return | not compiled | make the frame/projectivity equations and no-hidden-wall theorem variable-class and intrinsic |

The H11 row is the repaired theorem, not the failed root-seed handoff.  Its
separate
[finite-algebra/Galois calibration](dessin_native_field_effectivity_theorem.md)
certifies the projective rank locus \(a^2-b^2=0\), rejects the nonreduced
\([1:-1]\) branch, and identifies the accepted \([1:1]\) chart point with
reduced coordinate algebra and cyclic-product Vinberg field \(\mathbb Q\).
It does not prove that a variable-class carrier stack is zero-dimensional.
The
hostile replay constructs the ten formerly missing positive-degree cusp-lift
curves and independently recovers `508/42/0`.

The H14 and H17 incidence compilers genuinely start from raw marked P data.
Their last same-P gate still imports Kondō's effective-root exhaustion and the
frozen complete cusp-lift audits.  The finite incidence scheme alone does not
exclude every hidden negative curve.

## 2. What calibrates the general compiler

There are three logically different experiments.

1. **General determinant/Heegner stage.**  Bugaenko is the exact compact
   calibration.  It constructs the sign-correct type-IV arrangement directly
   from P and tests every mirror meeting the closed real chamber.
2. **General cusp stage.**  The toric `H6/9` row is the exact rational
   calibration.  At each of its three cusps the walls split as seven vertical
   and two bounded-penetration walls, and the full Vinberg parabolic
   neighborhood is replayed.
3. **Native divisor-surface stage.**  Toric `H6/9` and Coble
   `H10/30`, `H11/26`, `H14/60`, and `H17/48` provide actual
   `P -> surface -> same P` experiments.  They are stronger geometrically
   than an arbitrary complete-intersection carrier, but are scoped or
   individual inverse theorems, not evidence that every arithmetic P is a
   nef cone.

Thus none of the six principal rows is now *merely* surface-first: each has a
P-first arrangement or native-carrier direction.  The native Bugaenko K3
realization itself remains surface-to-polytope, while its orthogonal
determinant arrangement is P-first.  The two additional K3 determinant
controls below are surface-first/full-root controls, not general native
surface compilers.

## 3. Arrangement, pure product, and compensated product are different

| control | mere Heegner arrangement | pure/exact product | compensated slice-exact product |
|---|---|---|---|
| `E10`, `H9/10`, `9/1/0` | proved | **positive:** `II_(1,9)+<2> ~= 2U+E7(-1)` is complete 2-reflective; multiplicity-one weight `165` product; same complete root chamber | unnecessary for the complete-arrangement return |
| `U+E8(-1)+D4(-1)`, `H13/15`, `42/4/0` | available by determinant transfer | **negative:** the signature `(2,13)` determinant extension admits no exact global 2-reflective product by Wang's complete classification | general rational interior compensation exists, but no row-specific product is archived and boundary behavior remains separate |
| Bugaenko `H8/11` | proved | the standard pure norm-two column is obstructed; the complete globally minimal first compensator shell is also obstructed | a finite complex slice-disjoint correction exists abstractly; an effective Hilbert product is conditional on the stated Hilbert--Weil integral q-structure and has not been made explicit |

For Bugaenko, “pure obstructed” is deliberately limited to the standard
full-group column.  It does not prove that every finite-level reflective form
is impossible.  A compensated divisor can be exact on the marked real slice
while having extra global Heegner support, so it cannot be fed into an exact
reflective-form classification as though it were pure.

The other five principal rows have no row-specific determinant extension,
finite Heegner support, or obstruction matrix yet.  Their Borcherds entries
are therefore `open`, not negative.

## 4. Correct cusp encoding on a surface

Finite volume gives finitely many boundary-at-infinity points: ideal
vertices/cusps.  It does **not** make all cusp-wall curves pass through one
smooth surface point.  After toroidal resolution, the faithful surface datum
must retain the exceptional boundary cluster, its dual complex, the labels of
the meeting mirror branches, and the incidence map to the toroidal cusp
stratum.  This is the boundary decoration used in the ledger.

Ordinary faces whose relative interiors lie in hyperbolic space may still be
recorded by marked interior incidence data.  A face can have an ideal vertex
in its closure without itself lying at infinity.

## 5. The missing classification complexity

One coordinate is now canonical and literal.  The
[Gram--Belyi theorem](gram_belyi_coxeter_dessin_theorem.md) assigns every
marked arithmetic chamber the exact degree

\[
\beta_{\rm GB}(P)=\sum_{i<j}d_{ij}
\]

of its facet-pair Belyi leaves.  Its real-marked SNC surface has a natural
polarization with
$A_P^2=2(\binom f2+\beta_{\rm GB})$, and bounded
\(\beta_{\rm GB}\) gives a finite effective Gram/chamber/field list without
ABSW or Prasad.  This ledger does **not** invent row values: exact minimal
normalized-Belyi degrees for all pair entries have not yet been computed on
the six large examples.  The standard Coxeter-label calibration is now
frozen: \(m=3,4,6,\parallel\) have exact degrees \(5,3,7,2\), \(m=5\)
has degree at least \(29\), and a dotted \(u=2\) packet has exact degree
\(3\).  Therefore
\(\beta_{\rm GB}\ge2\binom f2+3N_3+N_4+27N_5+5N_6\)
before dotted excess.  The frozen `H10` prime-index family already proves that
this coordinate is unbounded on arbitrary finite-index descendants with all
ambient arithmetic data fixed.

The exceptional \(m=5\) lower endpoint has a separate exact search ledger.
Degree \(29\) is still **OPEN**.  The initial stage fixes
\(B(\lambda)=0\) and proves all marked passports with at most three black
vertices impossible.  Its
`3,977` low-vertex eliminations are joined by an exact no-go for all `1,372`
three-black/`WW` records (`784` distinct coefficient systems).  Thus `5,349`
of `42,610,678` Riemann--Hurwitz-compatible marked records are excluded and
`42,605,329` records with at least four black vertices remain.  The `210`
Galois-symmetric quadratic-factor trials also fail exactly.  These counts are
passport records, not dessins or polytopes; see the
[degree-29 reduction](gram_belyi_m5_degree29_reduction.md) and the
[three-black theorem](gram_belyi_m5_degree29_three_black_ww.md).

The four-black ledger is now completely closed by exact coefficient
elimination.  Its `42,994` marked records quotient to `12,752` master systems.
The first pass removes every two-black-endpoint and
one-black/critical-white case.  The second uses a degree-`756` critical
resultant and exact tower-field endpoint gcds to eliminate all `2,560` `WW`
endpoint-status systems / `8,960` records.  The final pass tests both `BW`
and `WB` endpoint values on `28,159` exact critical tower orbits, eliminating
the last `3,458` systems / `10,374` records.  Any degree-`29` witness must
therefore have at least five black vertices.  Including the completed lower
strata, the global remainder is `42,562,335` marked records and remains
**OPEN**; see the
[four-black reduction](gram_belyi_m5_degree29_four_black_reduction.md), the
[four-black `WW` theorem](gram_belyi_m5_degree29_four_black_ww_simple.md), and
the [four-black completion theorem](gram_belyi_m5_degree29_four_black_one_black_simple.md),
with its [hostile audit](gram_belyi_m5_degree29_four_black_one_black_simple_hostile_audit.md).

The five-black denominator is now exact as well: `190,047` raw marked
records quotient to `33,079` binary-status systems, split as
`BB/BW/WB/WW = 11,011/8,450/8,450/5,168`.  Projective modular standard-basis
certificates eliminate all `11,011` `BB` systems / `55,055` records; the
infinity charts and properness implication have passed hostile replay.  The
white-endpoint quotient is endpoint-status independent: `G=DQ+R` forces
every white endpoint to vanish on `Q`, reducing `BW/WB/WW` to
`4,225/4,225/1,292` status-free necessary systems.  Exact factorization `D=29(x-t)^4` closes the
whole `(5,1^24)` row, eliminating `3,876` more records.  The five-black
residual is `9,742` status-free necessary systems (`22,068` in the conservative
binary-status ledger) / `131,116` records and the global degree-`29`
remainder is `42,503,404`.  The next factorized partition `(3^2,1^23)` has a
uniform projective-infinity obstruction; exact bases close its entire
`a>=21` sublocus, another `12` systems / `48` records.  Its residual is
`1,280` systems / `5,120` records, while five-black and global raw remainders
become `131,068` and `42,503,356`.  The fixed `CC` factor `K=x(x-1)` yields
an exact dimension-`24` base algebra for every black record, and its value
equation is always a unit.  Thus all `1,292` `CC` passports are impossible;
after overlap the live remainders are `3,840` in `(3^2,1^23)`, `129,788` in
five-black, and `42,502,076` globally.  The one-critical towers then eliminate
all `1,280` unresolved `CS` rows and the predetermined `41`-row `SC` tranche.
The current exact remainders are `2,519` in `(3^2,1^23)`, `128,467` in
five-black, and `42,500,755` globally.  The broader census is stopped at this
boundary.  Both five-black completion and sharpness remain **OPEN**; see the
[five-black quotient](gram_belyi_m5_degree29_five_black_quotient.md), the
[five-black `BB` theorem](gram_belyi_m5_degree29_five_black_bb.md), and the
[five-black `WW` partition theorem](gram_belyi_m5_degree29_five_black_ww_partition5.md),
with its [hostile audit](gram_belyi_m5_degree29_five_black_ww_partition5_hostile_audit.md),
the [`(3^2,1^23)` high-degree theorem](gram_belyi_m5_degree29_five_black_ww_partition33_high_a.md),
its [hostile audit](gram_belyi_m5_degree29_five_black_ww_partition33_high_a_hostile_audit.md),
the [complete `CC` tower theorem](gram_belyi_m5_degree29_five_black_ww_partition33_cc_tower.md),
its [hostile audit](gram_belyi_m5_degree29_five_black_ww_partition33_cc_tower_hostile_audit.md),
the [`CS` theorem / `SC` tranche](gram_belyi_m5_degree29_five_black_ww_partition33_cs_full_sc_tranche.md),
its [hostile audit](gram_belyi_m5_degree29_five_black_ww_partition33_cs_full_sc_tranche_hostile_audit.md),
and the [surface-bound obstruction theorem](gram_belyi_carrier_dimension_field_bound_obstruction.md).

The larger carrier-classification object remains a **canonically normalized
minimal faithful Coxeter-dessin complexity tuple**

\[
\mathcal T_{\rm dess}(P),
\]

not the single, presently undefined scalar formerly denoted
\(T_{\rm dess}(P)\).  One must first fix the determinant extension, order the
separating principal levels, use the corresponding automorphic polarization,
and then minimize lexicographically over faithful decorated carriers.  The
tuple must include at least:

- the invariant Coxeter-cover degree \(\deg_{\rm Cox}(P)\), together with a
  chosen minimizing maximal reflection overgroup when the minimum has ties;

- the marked Heegner curves and every ordinary face stratum;
- the exceptional boundary clusters and dual complexes at all cusps;
- the sheet-separating arithmetic level; and
- the field degree, embedding dimension, Hilbert polynomial, total
  stratum-zero-cycle degree, curve degrees, and cleanup intersections needed
  to retain the root local system and real-chamber marking.

This normalization/minimization is a theorem target, not yet a computed
invariant on the six rows.  What is missing is a P-computable bound for every
entry of this tuple in terms of

\[
n,\quad [k:\mathbb Q],\quad |\operatorname{disc}L|,\quad
\operatorname{level}(L),\quad \#\operatorname{facets}(P).
\]

Such a bound would put all faithful carriers in finitely many explicit
Hilbert/Hurwitz stacks.  Their emptiness, components, and number-field points
could then be classified and pushed back to give a Coxeter existence,
finiteness, facet, or dimension theorem.  Without this bounded invariant,
the universal decorated surface exists but is too noncanonical and too large
to provide classification leverage.

The degree coordinate is no longer conjectural: the
[Schreier--Coxeter theorem](bounded_degree_schreier_coxeter_dessin_theorem.md)
proves finite classification and a computable congruence/Gram-height bound
when it is bounded.  Its ideal-polygon counterfamily proves that the other
ambient arithmetic coordinates alone cannot bound it.  The open
algebro-geometric leverage theorem is therefore to bound this degree from the
remaining marked Borcherds/Hilbert coordinates.

## 6. Replay and resource bound

```bash
python3 scripts/run_with_rss_limit.py --memory-gib 0.5 --timeout 60 -- \
  python3 scripts/certify_arithmetic_coxeter_dessin_calibration_ledger.py \
  --check-output data/enumerations/arithmetic_coxeter_dessin_calibration_ledger.json

python3 scripts/run_with_rss_limit.py --memory-gib 0.5 --timeout 60 -- \
  python3 -m pytest -q tests/test_arithmetic_coxeter_dessin_calibration_ledger.py
```

The producer also enforces an in-process `0.5 GiB` ceiling.  Its initial
archive build used about `0.027 GiB` peak RSS.

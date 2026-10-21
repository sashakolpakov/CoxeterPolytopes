# TOOT: the `A15_11` `H10/21` chamber has a K3 carrier

## Result

The known `A15_11_1_C02` Coxeter chamber now has an exact positive
algebro-geometric realization.  Its rank-11 chamber lattice (M) has a
unique nontrivial even overclosure (N) preserved by all 21 facet
reflections, with

\[
 [N:M]=4,
 \qquad
 N\cong U\oplus E_8(-2)\oplus\langle-2\rangle.
\]

The certificate gives explicit integer matrices for every link

\[
 M\hookrightarrow N\xrightarrow{\sim}
 U\oplus E_8(-2)\oplus\langle-2\rangle
 \xrightarrow{\sim}NS(X_2)^\iota
 \hookrightarrow NS(X_2).
\]

This produces two deliberately separated K3 statements.

1. A very general (N)-polarized projective K3 surface (Y) has
   (NS(Y)=N).  Hence all 21 reflections act integrally on the **full**
   Néron--Severi lattice of (Y), giving a full-rank mixed
   divisor-reflection subsystem.
2. On Shimada's named singular K3 surface (X_2), the same (N) occurs as
   the actual primitive invariant divisor lattice (F=NS(X_2)^\iota) of a
   rational involution.  All 21 vectors are genuine invariant divisors and
   all 21 reflections act on (F), but only 11 extend integrally to the
   rank-20 (NS(X_2)).

This closes the project's previously open K3/divisor-carrier entry for the
`A15_11` chamber.  It is not a new abstract polytope or arithmetic lattice.

## Exact lattice bridge

The source lattice has signature ((1,10)), determinant (8192), and Smith
form (2^9 4^2).  Exhaustive discriminant-module closure under the 21
reflections finds stable integral and even subgroup orders exactly `[1,4]`.
The unique index-four overclosure is even, has determinant (512), Smith
form (1^2 2^9), and discriminant-form invariants
((r,a,\delta)=(11,9,1)).

The archived matrices are not merely classification labels:

- (P^T G_NP=T), where (T) is the exact standard Gram of
  (U+E_8(-2)+\langle-2\rangle);
- (J^TTJ=G_M), with (|\det J|=4);
- (Q^TG_FQ=T), where (G_F) is the actual fixed Gram computed inside
  Shimada's published rank-20 basis; and
- the composite (QJ) has determinant four and satisfies
  ((QJ)^TG_F(QJ)=G_M).

The fixed basis is independently recomputed as the saturated integral
kernel of the chosen rational involution.  Its Gram inside the published
rank-20 (NS(X_2)) is exactly (G_F).  Thus the named-surface bridge does
not rest only on Nikulin's abstract uniqueness theorem.

## Primitive mirrors and the full-NS boundary

After transport to (F) and primitive normalization, the 21 mirrors are

| type in (F) | count |
|---|---:|
| ((-2,\operatorname{div}_F=1)) | 4 |
| ((-2,\operatorname{div}_F=2)) | 7 |
| ((-4,\operatorname{div}_F=2)) | 10 |

Every reflection is integral, involutive, and Gram-preserving on (F), and
the certificate checks the intertwining identity with the original action
on (M) for every facet.

The primitive embedding (F\subset NS(X_2)) changes divisibility: all 21
vectors have divisibility one in the full rank-20 lattice.  Consequently the
11 square-`-2` reflections remain integral, while the ten square-`-4`
reflections **do not extend** to integral isometries of (NS(X_2)).  They are
therefore reflections of the invariant divisor lattice, not full-NS
reflections and not surface automorphisms.

For a very general (F)-polarized K3 (Y), by contrast, (NS(Y)=F), so all
21 are full-NS integral.  Weyl re-marking can put a strict chamber point in
the ample cone.  The transported point has square `159/2`; its 21 oriented
wall pairings are all positive and have minimum `1`.  K3 Riemann--Roch then
chooses the effective sign of each of
the 11 primitive square ((-2)) classes.  They may be reducible; the ten
square ((-4)) divisors are not asserted to be curve walls.

## Finite volume, twice

The literal Vinberg audit gives:

- 4,159 elliptic rank-9 edge diagrams, each with exactly two endpoints;
- 607 finite vertices;
- 50 cusps;
- full-rank zero-neighbourhoods for every irreducible affine component; and
- zero bad edges and zero bad parabolic neighbourhoods.

An independent GMP `cddlib` double-description replay enumerates all 657
extreme rays of the complete chamber cone: **607 finite**, **50 ideal**, and
**0 hyperideal**.  Of the ideal rays, 18 lie on 11 facets and 32 lie on 12.
The absence of hyperideal rays independently rules out the infinite-volume
failure mode.

Thus the polytope is finite-volume and noncompact.  Its reflection group has
exact index (6720) in the known 13-wall full reflection group for
(U+E_8(-2)+\langle-2\rangle).

## Reusable criterion

The construction isolates a useful positive-search theorem.  Suppose an
even Lorentzian chamber lattice (M) has a finite-index even overclosure
(N) preserved by every facet reflection; (N) embeds primitively in the
K3 lattice; and the chamber passes a literal finite-volume test.  Then a
very general (N)-polarized projective K3 has (NS=N) and realizes the
chamber as a full-rank divisor-reflection subsystem.  Weyl re-marking and
Riemann--Roch orient every primitive square ((-2)) facet effectively.

Promotion to the full nef cone is a separate problem: one must check every
ambient root by the one-sign criterion.  That gate is not passed here, so
this report makes **no claim that the 21 walls form the full nef cone**.
Nor is this Shimada's distinct 37-wall intersection of the theorem-complete
478-wall induced nef chamber with (F).

## Lineage and reproduction

Alexeev supplies the source chamber; Shimada supplies (X_2), its complete
induced chamber, and the genuine rational involution; the known 13-wall
denominator supplies the maximal reflection group.  The exact composed
index-four K3 bridge appears absent from those sources and from the frozen
pre-catch project ledger.  This is bounded evidence only, not a worldwide
priority claim.

Permanent files:

- [certificate](../data/enumerations/a15_11_h10_21_x2_k3_carrier.json)
- [certifier](../scripts/certify_a15_11_h10_21_x2_k3_carrier.py)
- [focused tests](../tests/test_a15_11_h10_21_x2_k3_carrier.py)

Replay from the repository root:

```bash
PYTHONPATH=.:src python3 scripts/certify_a15_11_h10_21_x2_k3_carrier.py --check-archive --replay-cdd
PYTHONPATH=.:src pytest -q tests/test_a15_11_h10_21_x2_k3_carrier.py
```

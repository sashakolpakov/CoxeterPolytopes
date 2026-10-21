# Exact four-port repair after the strengthened three-port scan

The [four-port simplex/prism facet](fresh_sqrt5_joint_four_port_domain_20261017.md) finds **47 violated inequalities at 23 finite vertex forests** in the restricted rank-13 dual after the strengthened three-port repair. Raising vertex masses by the exact per-forest maximum gives

\[
\frac{22302700899796044000000000000000}
{1805838229651556529708934923143}
=12.350331571006466\ldots>12.
\]

The new inequalities improve the diagnostic, but the same fixed face/provenance masses still obstruct target 12. This is not a global optimization, a complete joint local-packing calculation, or a dimension theorem.

## Same masses, stronger simultaneous constraints

The starting vertex masses are exactly those in the frozen [strengthened three-port diagnostic](fresh_sqrt5_joint_triangle_square_ports_20261017.md), whose ratio is \(12.355268075778007\ldots\). Its 2,219 face masses are unchanged: 1,864 simplex/edge rows and 355 actual selected dotted-quad provenances. These are the same masses already certified against the [complete nonnegative local simplex-packing cone](fresh_sqrt5_restricted_flag_local_cones_20261017.md).

For each four-port set \(U\), the new charge has cost three and coefficients

\[
\text{simplex corner: }(|C_S\cap U|-1)_+,
\qquad
\text{selected quad corner: }
\mathbf 1\{\text{at least two occupied square pairs lie in }U\}.
\]

The proof applies to the selected canonical **dotted-prism** aggregates. It is not extended to general non-dotted full-star product aggregates.

The new vertex mass is

\[
\mu'_v=\max\left\{\mu^{(3)}_v,
\max_{U\text{ at }v}\frac{L_4(v,U)}3\right\}.
\]

Thus every old vertex mass increases or stays fixed. All inherited positive-cost constraints, the full local simplex-packing covers and the strengthened three-port cuts remain valid. Face masses are unchanged, so all 94,481 signed constraints of the saved model, including the 77,391 coface-flag transfers, also remain valid. This monotonic inheritance preserves the previously checked model with 82,806 face obligations and 6,558 vertex budgets; this diagnostic does not claim to reconstruct that large matrix again.

## Complete four-port and actual-provenance coverage

| Quantity | Exact count |
|---|---:|
| Positive simplex/edge rows | 1,864 |
| Simplex block orbits | 3,485 |
| Actual selected quad provenances | 355 |
| Actual quad-vertex occurrences | 2,418 |
| Contributing finite forests | 593 |
| Literal four-port subsets | 423,995 = 593 × 715 |
| Four-port automorphism orbits | 90,619 |
| Violated inequalities | 47 |
| Increased vertex masses | 23 |

For each finite forest, all \(\binom{13}{4}=715\) subsets are grouped by the exact fully marked weighted-forest code. The whole literal inequality is averaged over each orbit. Simplex coefficients use the complete literal block orbits reconstructed at these actual masses, including the previously missing 397 simplex rows.

Every selected quad row keeps its one fixed actual provenance from the frozen model map. The checker independently reconstructs its core plus the same finite spectators, obtains the finite vertex forests and occupied pairs from the complete K2,2 union, and compares every omitted-root pair, forest and square footprint with the frozen three-port export. It neither chooses a different representative of a merged row nor replaces prism vertices by all positive-definite two-deletions.

For a simplex block of size \(k\), the coefficients satisfy the exact identity

\[
\sum_{|U|=4}(|C\cap U|-1)_+
=\sum_{j=2}^{\min(k,4)}(j-1)\binom{k}{j}\binom{13-k}{4-j}.
\]

For a prism-corner star with \(e\) leaves, the quad coefficients satisfy

\[
\sum_{|U|=4}\mathbf1\{\text{at least two star pairs in }U\}
=\binom e2(12-e)+\binom e3.
\]

The first term chooses the hub, two leaves and one outside port; the second chooses the hub and three leaves. Both identities are checked for every corresponding block or corner. A one-pair square footprint contributes zero, and all its four-port subsets are accounted for by that zero coefficient.

## A concrete violated orbit

For vertex forest index 218,

\[
5A_1+A_2+A_3+H_3,
\]

the canonical four-port representative \(\{0,1,7,9\}\) has simplex load
\(14672684353/11000000000000\) and quad load
\(3957147/3125000000\). Its required mass exceeds the input three-port vertex mass by exactly

\[
192268781/3000000000000.
\]

The certificate records this orbit and all 46 other violated columns, with their rational coefficients on the common face rows. It also retains every tested orbit load, all exact coefficient-sum identities, reconstructed quad footprints and final repaired masses. Coefficients may be set to zero outside the exported face classes to obtain valid weaker nonnegative charges.

## Artifacts and validation

The [source](../scripts/fresh_sqrt5_joint_four_port_diagnostic_20261017.py) and [saved certificate](fresh_sqrt5_joint_four_port_diagnostic_20261017.json) are a new successor; all frozen predecessors remain unchanged. Generation session 47891 and fresh saved replay session 93338 both terminated with exit 0 and `PASS_EXACT_FOUR_PORT_FACET_AFTER_STRENGTHENED_THREE_PORT_REPAIR`. Replay repeated all reconstruction, orbit enumeration, coefficient identities and rational load checks with no optimizer.

Frontier_literature independently read the source and mathematical normalization: the common actual-provenance map, K2,2 corner usage, all-715 coverage, automorphism averages, sum identities, cost three and monotone repair passed that narrow review. This was a source/proof read, not a duplicate arithmetic replay.

Source SHA256: `ec32211a3a936f066056c0a3152dc105405e66ee06c7125e4d06829b18bb1f1e`.
Certificate SHA256: `e4f7f81c29b12833c14e4485cb534f93d435cef12b306385277626b1044f1fe4`.

No optimizer was run. The result does not exclude stronger five-port or other genuinely joint geometric constraints, and it makes no compact ambient realization or dimension claim.

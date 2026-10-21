# Exact five-port facets after the four-port repair

All fifteen frozen [five-port facets](fresh_sqrt5_joint_five_port_domain_20261017.md) were evaluated on the same restricted rank-13 face masses, starting with the [four-port repaired vertex masses](fresh_sqrt5_joint_four_port_diagnostic_20261017.md). They give **115 violated columns at 18 forests**. The exact simultaneous repair has ratio

\[
\frac{66908102699388132000000000000000}
{5419878674792788039126804769429}
=12.344944732908832\ldots>12,
\]

compared with the input \(12.350331571006466\ldots\). This is a complete exact scan of these fifteen facets, with no optimizer. It is not a complete joint packing cone, a new global optimum, or a dimension theorem.

## Preserved support and constraints

The 2,219 face masses are unchanged: 1,864 maximal-simplex/edge rows and 355 fixed **actual** selected dotted-prism provenances. No merged row is represented by a different source, and no simplex mass is restored or altered.

The scanner reconstructs each finite quad-corner forest and its occupied square pairs, and compares them with the frozen four-port export. As before, quad vertices come from the complete K2,2 union, not all positive-definite two-root deletions. General non-dotted product aggregates are outside the selected scope.

At each vertex the repaired mass is the maximum of its old mass and every new facet load divided by that facet's cost. Therefore all earlier constraints remain valid: inherited positive-cost columns, all 94,481 signed columns including the coface flags, the complete local simplex-only packing covers, and the strengthened three-port and four-port charges. This follows by the checked unchanged face masses and monotone vertex-mass increase; the large inherited matrix is not rebuilt in this diagnostic.

## Complete domain and exact-size histograms

| Quantity | Exact count |
|---|---:|
| Contributing finite forests | 593 |
| Simplex/edge rows | 1,864 |
| Simplex block orbits | 3,485 |
| Actual selected dotted-quad provenances | 355 |
| Quad-vertex occurrences | 2,418 |
| Literal five-port subsets | 763,191 = 593 × 1,287 |
| Five-port automorphism orbits | 143,573 |
| Facet–orbit inequalities | 2,153,595 |
| Violated columns | 115 |
| Increased vertex masses | 18 |

Every \(\binom{13}{5}=1287\) subset is included at every contributing forest. The whole inequality is averaged over the exact marked weighted-forest automorphism orbit. The stored load vector has coordinates

\[
(s_2,s_3,s_4,s_5,q_1,q_2,q_3,q_4).
\]

These are **exact intersection-size histograms**, not threshold indicators. This is necessary because some facet coefficients decrease between successive simplex sizes. For example, facet 1 has simplex coefficients \((1,0,1,2)\).

For each simplex block of size \(k\), every histogram bin satisfies

\[
\sum_{|U|=5}\mathbf1\{|C\cap U|=j\}
=\binom kj\binom{13-k}{5-j},\qquad j=2,3,4,5.
\]

For a selected quad-corner star with \(e\) leaves,

\[
\sum_{|U|=5}\mathbf1\{\text{exactly }j\text{ occupied square pairs in }U\}
=\binom ej\binom{12-e}{4-j},\qquad j=1,2,3,4.
\]

A positive quad bin requires the hub in \(U\), together with exactly \(j\) leaves and \(4-j\) other ports. The scanner explicitly covers all such subsets; hub-absent and leaf-free subsets have zero coefficients. Both identities are checked separately for every block or corner. All fifteen facet coefficient vectors are then dotted with the common exact histogram, including their nonmonotone entries.

## Which facets add constraints

Facet indices refer to the frozen five-port certificate. The counts are before taking a per-vertex maximum; a forest may violate several facets.

| Facet | Violated columns | Forests |
|---:|---:|---:|
| 0 | 0 | 0 |
| 1 | 19 | 15 |
| 2 | 0 | 0 |
| 3 | 19 | 4 |
| 4 | 18 | 6 |
| 5 | 6 | 5 |
| 6 | 1 | 1 |
| 7 | 1 | 1 |
| 8 | 5 | 2 |
| 9 | 8 | 3 |
| 10 | 20 | 5 |
| 11 | 1 | 1 |
| 12 | 1 | 1 |
| 13 | 13 | 3 |
| 14 | 3 | 1 |

The largest individual required-mass increase occurs at forest 435,
\(4A_1+2A_2+H_3+I_2(5)\), for ports \(\{0,1,2,3,11\}\) and facet 5,

\[
s_2+7s_3+9s_4+10s_5+q_1+q_2+q_3+q_4\le10.
\]

It raises the required mass above the four-port baseline by exactly
\(100845889/825000000000\). All 115 violated columns are exported with their rational coefficients on the common face rows. Coefficients outside these represented classes may be set to zero as valid weaker nonnegative payments.

## Verification and scope

The [source](../scripts/fresh_sqrt5_joint_five_port_diagnostic_20261017.py) and [certificate](fresh_sqrt5_joint_five_port_diagnostic_20261017.json) use a fresh successor namespace. The necessary-domain geometry and source received a separate [independent narrow review](fresh_sqrt5_joint_five_port_domain_20261017_review.md). All frozen inputs remain unchanged.

Generation session 79111 and fresh saved replay session 49378 both terminated with exit 0 and `PASS_EXACT_FIFTEEN_FIVE_PORT_FACETS_AFTER_FOUR_PORT_REPAIR`. Replay reconstructed every feature, histogram, identity, facet load, exported violated column and repaired vertex mass, without an optimizer.

Source SHA256: `c54e91602ee8d8d24ee927481b72a29a89cc27c797a744090d50fc1d037e5acc`.
Certificate SHA256: `e42997ac1e8804adf8ead5a8118f656dcb5e23391b82e66af7fb97cda692bb1e`.

The repaired ratio remains above 12. No inference is made about stronger joint configurations, additional signed incidence relations, missing geometric face classes, or ambient realizability.

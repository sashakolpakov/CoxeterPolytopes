# GitHub proof catalogue and complete addenda: dimension bounds 57/17/15


The arXiv article is the concise consolidated presentation of all three
bounds. Complete analytic case proofs, computational proofs, certificates
and reproduction records remain in this existing supplement. No additional
manuscript family is created.

The complete 30 October source and PDF are also preserved byte for byte at
[the repository supplement](../../manuscripts/addenda/dimension_bounds_20261030/).
The immutable links below identify the exact proof version used by the paper.

All mathematical statements are taken from the existing dimension-bound
paper family. The complete source cited by the concise edition is frozen at
[`36eeaddf646eff25d4990c9ab730a20ada21f80f`](https://github.com/sashakolpakov/CoxeterPolytopes/tree/36eeaddf646eff25d4990c9ab730a20ada21f80f).

The [complete manuscript](https://github.com/sashakolpakov/CoxeterPolytopes/blob/36eeaddf646eff25d4990c9ab730a20ada21f80f/arxiv/hyperbolic_coxeter_dimension_bounds/main.tex)
has Parts I--III for the main exposition and Parts IV--VI for the complete
case proofs. Its
[551-page PDF](https://github.com/sashakolpakov/CoxeterPolytopes/blob/36eeaddf646eff25d4990c9ab730a20ada21f80f/manuscripts/build/arxiv_hyperbolic_coxeter_dimension_bounds_20261028_current_20261030.pdf)
is the corresponding reading copy. The proof appendices start on pages
72, 107 and 372 for finite volume, compactness and arithmetic, respectively.

| Concise source label | Complete source labels and proof home |
|---|---|
| `lem:positive` | `compact:lem:inverse`, `compact:lem:projection`; finite/affine feasibility also in the finite-volume local interface |
| `lem:payment` | `fv:lem:payment-comparison` |
| `prop:local` | `fv:thm:local`, `fv:prop:finite`; `fv:app:routing`, `fv:app:weights`, `fv:app:finite-contexts`; dimension-specific final rule |
| `lem:shape` | `fv:lem:shape` |
| `prop:surplus` | `fv:lem:concentration`, `fv:prop:terminal`, `fv:lem:low8`, with the actual factor envelope `fv:eq:envelope` |
| `prop:fvtail` | `fv:app:tail`, `fv:app:narrow-tail`, including the complete terminal-tail lemma |
| `lem:pointpressure`, `prop:shadows`, compact original vector and references | `compact:eq:original`, `compact:eq:allpoint`, `compact:eq:square`, `compact:eq:payment` |
| `prop:compactcover` | `compact:sec:forests`, all listed row consumers, `compact:thm:tail`, `compact:app:tail:completeconjunction` |
| `prop:d18` | `compact:sec:d18-setup`, `compact:sec:composition`, `compact:lem:localA18`; original references, head/fork inventories, and all six terminal fork pairs |
| `lem:absorption` | `arith:app:arithsummedfeet`, `arith:eq:arithAsummed`; the concise lemma isolates its analytic combination |
| `prop:arithreduce` | `arith:prop:f17`, `arith:prop:old17`, `arith:prop:a17reduce`, `arith:prop:f16`, `arith:prop:onlyd16`; full dispatches in Part VI |
| Four finite output propositions | `arith:prop:pair-output`, `arith:sec:algorithms`, and their separate actual-domain inclusions and geometric consumers |
| `prop:d16reduce` | `arith:prop:d16reduce`, `arith:lem:D16initial`, `arith:eq:reference`; full rank-sixteen endpoint and field reduction |
| `prop:dyadic` and its two geometric finishes | `arith:sec:odd-inert`, `arith:sec:ramified`, `arith:sec:split`, `arith:prop:noD16`; the two alternatives collect the explicit terminal cases already printed there |
| Final arithmetic composition | `arith:app:arithfinalimplication`, with the independent `compact:thm:main` |

The concise text retains the main deductions and supplies proof references
for omitted case analyses. Those references are essential mathematical
dependencies. The two terminal dyadic alternatives are necessary residual
possibilities; they are excluded geometrically and are not realized examples.

## Catalogue statements

The identifiers below are used in the concise paper. They name statements
proved in the fixed complete source, not new executions or new examples.
The full checking domains, finite arrays and continuous case proofs are
those at the indicated source labels. Necessary algebraic survivors and
forest lists are not asserted to be realized polytopes.

| Item | Mathematical statement and complete proof location |
|---|---|
| **FV1** | The single intrinsic marked-corner rule in ambient dimension 58 pays at least 2 on every triangular two-face and at least 1 on every quadrilateral. Coverage includes all extra supporting walls and nonsimple ideal vertices. Proof: `fv:thm:local`, `fv:app:routing`, `fv:app:weights`, and the dimension-specific guarded rule. |
| **FV2** | The maximum formal charge over the complete rank 58 spherical-forest domain, including all contextual and guarded terms, is 377. This is a necessary-domain maximum, not a realized forest. Proof: `fv:prop:finite`, `fv:app:finite-contexts`, “The new complete finite maximum.” |
| **FV3** | Over the complete 2552 terminal cusp profiles of total rank 57, the minimum normalized envelope margin is 87/10240 and the minimum absolute margin is 36. Convex-order concentration then gives actual ideal surplus at least 36 for every cusp. Proof: `fv:prop:terminal`, `fv:lem:concentration`, `fv:lem:low8`, `fv:eq:envelope`. |
| **FV4** | No finite-volume Coxeter polytope exists in dimensions 59–163 by their stated dimension-specific rules; none exists from 164 onward by the independent narrow-facet/two-budget argument. Proof: `fv:app:tail`, `fv:app:narrow-tail`. The latter includes the full-budget lemma before first failure; no dimension heredity is assumed. |
| **C1** | At an original farthest vertex, the stated ordinary, defect, chain and prefix energy estimates give the complete necessary rank 18 cover. Every retained row except whole D18 is excluded by its actual-rank consumer. Proof: `compact:sec:components`, `compact:sec:forests`, and Parts II/V's complete component and proper-row case proofs. |
| **C2** | WholeD18 is excluded by the complete head/fork classification and both signed assignments. The separate actual-rank exclusions 19–29, followed by Vinberg's published ceiling, exclude all larger compact dimensions. Proof: `compact:sec:d18-setup` through `compact:sec:composition`, `compact:thm:tail`, `compact:app:tail:completeconjunction`. |
| **A0** | The arithmetic chain inventory, actual golden-field reduction, complete helper packing and simultaneous absorption of all isolated feet give the spectator-uniform component estimates, including e(A15)>.8348, e(A16)>.8703 and e(A17)>.9032. All pressures refer to one original vector. Proof: `arith:app:arithmeticchain`, `arith:app:arithgoldenfield`, `arith:app:arithgoldenpacking`, `arith:app:arithsummedfeet`, `arith:app:arithgoldenend`. Smaller chain and conditioned component inputs have their separate full proofs in Part VI. |
| **A1** | The complete quadratic A17 endpoint domain is empty: 25 discriminants below 88, 53 strip models and 159 ordered pairs. The first rejections are 106 full-field-square failures, 52 long-cycle failures and one two-cycle failure. Exact domain: `arith:eq:a17strip`, `arith:eq:a17cycles`, `arith:prop:pair-output`(i). |
| **A2** | Both cubic A17 finite-anchor domains are empty. The 2432 coefficient triples give nine strip models; all fail the full field-square or long-cycle tests. Exact domain and full integer rings: `arith:prop:pair-output`(ii), `arith:sec:algorithms`, `arith:tab:cubic`. |
| **A3** | The shared quadratic B/D17 algebraic domain has precisely two necessary pairs, in Q(sqrt13) and Q(sqrt17), as stated below. Exact domain: `arith:prop:pair-output`(iii). Whole B17 and whole D17 have **distinct actual geometric inclusions and consumers** in Part VI. |
| **A4** | The complete quadratic whole B16 domain has only U=5+2sqrt6, V=10+4sqrt6. Among 273 pairs, 168 fail half-product integrality and 104 fail full field-square membership. Exact domain: `arith:prop:pair-output`(iv). The sole necessary pair fails the actual geometric consumer in `arith:app:arithB16domain`. |
| **A5** | The two fresh arithmetic covers and their complete local consumers exclude rank 17 and reduce rank 16 to whole D16; its actual head/fork reduction gives the quadratic triple in the paper. Proof: `arith:prop:f17`, `arith:prop:f16`, the complete dispatches, `arith:prop:d16reduce`, `arith:lem:D16initial`. |
| **A6** | The necessary quadratic whole D16 triple is impossible. The full manual proof exhausts odd norm and every inert, ramified and split even-norm branch before the final actual-vector contradictions. Proof: `arith:sec:odd-inert`, `arith:sec:ramified`, `arith:sec:split`, `arith:prop:noD16`. There is no whole D16 census premise. |

### The four finite arithmetic domains

For **A1**, put $Q(X)=X(17X-36)$. The endpoint generators are
in the same actual quadratic field, with
$36/17<U,V<3969/425$ and $0<U',V'<36/17$.
Require $R=Q(U)Q(V)$ to be a square in that field,
$g=(UV-\sqrt R)/18$ and $L=g^2/(UV)$ integral, the physical
positive square root, $g\le0$, and either $g=L=0$ or
$0<L'<4$. The physical own and mutual parameters retain the full
native/dotted tests specified in the source; parallel is excluded
by actual compactness, not silently omitted.

For **A2**, the anchor is the largest root of
$X^3-5X^2+6X-1$ or $X^3-6X^2+9X-3$. Their full integer rings
are used. The complete boxes for $V=a+bU+cU^2$ are
$[-4,4]\times[-9,9]\times[-3,3]$ and
$[-6,6]\times[-9,9]\times[-2,2]$. They come from the embedding
strips and Lagrange interpolation before enumeration. The same actual
time-line, cycle and native tests are applied. Both signed norm branches
are retained in each full field-square decision.

For **A3**, the full strips are
$4<U<11$, $0<U'<2$, $2<W<676/17$, $0<W'<4/17$,
with both generating the actual quadratic field. Set
$g=UW-\sqrt{U(U-2)W(17W-4)}$ and $L=g^2/(2UW)$.
The square is a full field square; $g,L$ are integral,
$0<L'<4$, and the physical mutual parameter is native or dotted.
The two survivors are

| Field | U | W |
|---|---|---|
| Q(sqrt13) | (5+sqrt13)/2 | 11+3sqrt13 |
| Q(sqrt17) | (7+sqrt17)/2 | (21+5sqrt17)/2 |

For **A4**, the exact strips are
$2<U<12$, $0<U'<2$, $4<V<64$, $0<V'<1/4$, and $UV/2$ integral.
The source additionally tests the physical positive field root
$k^2=(U-2)V/[4U(4V-1)]$, integral
$G=4Uk(4V-1)-UV\ge0$ and $L=G^2/(2UV)$, the identity
$L+G=((15U-32)/2)V-2U+4$, the zero branch $G=L=L'=0$
or $0<L'<4$, and the actual-field native/dotted alphabet for U and nonzero L.
A square rational norm is never accepted as sufficient for a field square.

The full algorithms and their soundness, complete-ring coefficient bounds,
termination and exact signs are printed in the complete source. Their
computational proofs remain here in the repository rather than being
repeated in the concise article.

### Verification records and availability

The finite57 current-rule conjunction records 259 actual scoped units,
all 64 selected shards, 28,028 maps, 12,518 nonunit applications and 83 prefixes:
[`finite57_current_rule_release_completion_20261025.json`](../../results/finite57_current_rule_release_completion_20261025.json).
It remains conditional on the explicit inherited geometric and finite-data
inputs; it is not the combined 146-job replay or full transitive regeneration.
The proof maps below bind the other statements to their precise existing
producer, replay and selected conditional Lean records. No mathematical or
Lean job was rerun to prepare this concise article. The
[public-data guide](../../repro/PUBLIC_DATA.md) distinguishes available
repository inputs from the separate 936-file/33-archive bulk collection;
no bulk asset from that collection has been uploaded.


The reproduction entry points at the same version are:

- [Finite volume 57](https://github.com/sashakolpakov/CoxeterPolytopes/blob/36eeaddf646eff25d4990c9ab730a20ada21f80f/docs/finite_volume57_reproduce.md).
- [Compact 17](https://github.com/sashakolpakov/CoxeterPolytopes/blob/36eeaddf646eff25d4990c9ab730a20ada21f80f/docs/compact17_consolidated_proof_map_20261025.md).
- [Compact arithmetic 15](https://github.com/sashakolpakov/CoxeterPolytopes/blob/36eeaddf646eff25d4990c9ab730a20ada21f80f/docs/compact_arithmetic_consolidation_proof_map_20261025.md).
- [Golden-lattice primary-source authentication](https://github.com/sashakolpakov/CoxeterPolytopes/blob/36eeaddf646eff25d4990c9ab730a20ada21f80f/docs/bugaenko_golden_lattice_cutoff_authentication_20261026.md).
- [Reproduction package](https://github.com/sashakolpakov/CoxeterPolytopes/blob/36eeaddf646eff25d4990c9ab730a20ada21f80f/repro/README.md)
  and [input availability](https://github.com/sashakolpakov/CoxeterPolytopes/blob/36eeaddf646eff25d4990c9ab730a20ada21f80f/repro/PUBLIC_DATA.md).

The concise edition adds no new mathematical executions and does not alter
the source bindings or scopes of those historical records.

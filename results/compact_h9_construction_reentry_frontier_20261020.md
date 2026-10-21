# Compact H9 construction reentry: the prime-seven H3+A6 lane

20 October 2026. Read-only source/metadata research, with the paper calculations below. No root enumeration, scientific module import, checker execution, Lean, TeX, installation or download ran. This selects one construction obligation; it reports neither an example nor a new nonreflectivity theorem.

**Recommendation:** retain the literal first-type lattice `L=<-phi>+H3+A6` over `O=Z[phi]`, and replace its old global-height search by one complete outgoing-edge problem. The existing 45-wall prefix has infinite volume. The inspected current records do **not** establish nonreflectivity of its full reflection group. Its mathematical viability means that the known exclusions below do not decide it, not that compact completion is likely or literature-wide novelty has been established.

## Closed sources remain closed

These are literal-lattice/family conclusions. Their stated scope does not silently extend to every invariant lattice, similar form or commensurability class. Earlier route notes calling the same hosts open are historical.

| Host | Current obstruction and authoritative source |
|---|---|
| `<-sqrt3>+I9` | Full scalar facet, complete edges and chamber period of trace `5+3sqrt3`; [conclusion](compact_h9_sqrt3_nonreflectivity_conclusion_20261020.md), SHA256 prefix `00244030f08b`. |
| `<-(3+sqrt13)/2>+I9` | Four actual edges and period trace `268eta+80`; [conclusion](compact_h9_sqrt13_nonreflectivity_conclusion_20261020.md), `1da6e1ce9c04`. |
| Golden E8 binary norm31 / norm59 / k4 norm319 | Full chamber periods `28+42phi`, `40+55phi`, `158+245phi`; [31 proof](compact_h9_golden31_unit_bisector_preparation_20261020.md) `7259c05e994f`, [59 conclusion](compact_h9_golden59_nonreflectivity_conclusion_20261020.md) `0d54c15665f0`, [k4 proof](compact_h9_golden_k4_transverse_edge_and_period_candidate_20261020.md) `7d12748c606e`. The k4 closure is reviewed written proof; executable reproduction remains pending. Golden31 is distinct from the older silver norm31 host. |
| Unit-coupled golden E8 family `L_k` | Closed literal members `k=-1,0,1,2,3,4` and the entire range `k>=2269`; the other finite-range members are not exhausted. The [full-base proof](golden_e8_binary_family_full_base_facet_20261020.md) and [collar/cutoff](golden_e8_binary_family_discriminant_collar_and_finite_cutoff_20261020.md), `841a4b7f9e75`, do not authorize a host-by-host census. |
| Golden `<-phi>+I_n+D4`, n>=5 | Full actual germs and period trace `10+11phi`, then full coordinate facets; [proof](golden_dyadic_r5_chamber_period_candidate_20261020.md), `e2aeff93b42f`. |
| Odd/ramified silver diagonal families | `<-(1+sqrt2)>+I_m`, m>=7, and `<-sqrt2>+I_n`, n>=8, excluded by the exact published cutoff plus full-facet descent; [composition](ramified_silver_uniform_published_facet_closure_20261020.md), `510a6d31d2c3`. |
| Ramified cubic-seven diagonal family | `<-t(2-t)>+I_n`, n>=3, excluded by full facets and incompatible chamber involutions; [proof](cubic7_ramified_diagonal_uniform_noncompactness_20261020.md), `90d097fb04ff`. The separate unimodular cubic-seven primary-classification gate remains unresolved, so it is not a justified new construction lane. |
| Biquadratic diagonal family | `<-(phi+sqrt2-1)>+I_n`, n>=8, excluded by exchanging two full-face operations; [proof](compact_biquadratic_unit_order_exchange_obstruction_20261020.md), `043f1f1f6605`. |
| Cyclic-quartic F4+F4+1 hosts | Both coefficients `T-1` and `(T-1)(2+T)`, `T=2cos(pi/8)`, already have full-chamber infinite-order stabilizers: [first](fresh_h9_quartic_next_frontier_20261013.md) `667bfa9757cf`, [second](fresh_h9_new_host_next_20261013.md) `419bd39fd962`. Their older bounded-prefix reports must not reopen them. |

Also avoid the positive-coordinate shortcut: the [golden diagonal barrier](golden_determinant3_conditional_closure_and_diagonal_barrier_20261020.md), `e0c30a9345ad`, proves a full inherited diagonal facet for every nonunit totally positive coefficient. Its noncompactness conclusion explicitly retains the exact published diagonal-cutoff premise. Rational noncompact constructions and dimension-preserving doubling do not supply compact H9.

## The surviving literal candidate and its distinction

Use the saved H3 Gram `[[2,-phi,0],[-phi,2,-1],[0,-1,2]]` and the ordinary A6 Cartan matrix. Their positive minors are `2,3-phi,2phi^-2` and `det(A6)=7`. Thus

`det(L)=-14phi^-1`, with signatures `(9,1)` physically and `(10,0)` at `phi'=1-phi`. Conjugate positive definiteness rules out every nonzero K-isotropic vector. With the usual admissible arithmetic-lattice theorem, a cofinite integral reflection group here would therefore have a **compact** chamber. A finite prefix escape proves neither full-group nonreflectivity nor an ideal completion.

The inspected report establishes the full initial germ `H3+A6` at `p=(1;0;0)` and the safe complete primitive norm alphabet `{1,2,4,7,14,28}`, modulo unit squares. It records 20 complete time controls through height `H=phi*x^2/q<=30` and 45 accepted walls. The independent audit checks the zero census, controls, roots, reflection integrality, angles and escape; its full positive-height spatial census is source-reviewed, **not independently rerun**. That existing evidence boundary remains explicit.

There is a useful new separation from the *entire specified unit-coupled* E8 family, not only its first closed members. The prime7 is inert in O. Modulo7, `phi^8=-1`, so phi has order16, whereas `-4=3` has order6. Hence `4+phi^(2k+1)` is never divisible by7. The E8 binary determinant `-phi*(4+phi^(2k+1))` has valuation0 there, while `det(L)` has valuation1. In rank10 scalar similarity changes determinant by a square, so these spaces are not K-similar. This does **not** distinguish L from every possible E8-plus-binary form; an arbitrary binary block can contain the prime7.

## A first finite edge follows from the saved census

The following is a paper derivation using that inherited complete census, not a fresh execution. Release the middle H3 simple wall. Write

`w=(1+2phi, 2+2phi, 1+phi)`, so `H3*w=(0,1,0)` and `w^2=2phi^2`.

Saved zero-based wall10 is `r=(phi;w;0)`, of norm1. Its finite contact along `(1;t*w;0)` is

`P=(1;w/2;0)`, with `P^2=-1/(2phi)`.

For any positive-time mirror `m=(x;z)` of norm q meeting `0<=t<=1/2`, spatial Cauchy gives `phi^2*x^2 <= (phi^2/2)*(q+phi*x^2)`, hence `x^2<=phi^2*q` and `H(m)<=phi^3`. The archived complete H30 census has no positive wall below this height and exactly two at it: wall10 and wall9 `s=(phi;0;phi*(1,1,1,1,1,1))`. Wall9 is strictly negative on the entire segment. Consequently r is the unique first new simple wall. The endpoint has full germ `A6+3A1`, with the three A1 roots the retained H3 ends (norm2) and r (norm1). This uses the complete initial germ and the bounded Vinberg census, not just a compatible root list. Independent review should precede promoting this derived endpoint to a newly certified record.

## Precise next edge, and the present bottleneck

From P release an **end A6 wall**, retaining `A5+3A1`. Let Omega be the corresponding A6 end fundamental weight, `Omega^2=6/7`. The outgoing ray is `P+s*(0;0;Omega)` and reaches the light cone at

`s_*^2=7/(12phi)`.

Reading the saved triples `(time,H3-second-coordinate,A6-first-coordinate)` shows that every positive old collision is at least `s=phi^-1`; this is beyond the light point, since `phi^-2>7/(12phi)` is equivalent to `phi<12/7`. Thus none of the 45 existing walls gives a finite upper for this edge. Replaying that prefix cannot close the obligation.

The next mathematical task is the **complete finite-extension profile problem for A5+3A1**, followed by the corresponding binary norm equations with integral reflection constraints. In particular do not omit the isolated prime7 norm classes simply because the saved prefix contains none. An isolated-profile normal has exactly the form

`m=(2x; x*w; 7z*Omega)`, `x,z in O`,

`42z^2-2phi^-1*x^2=q`, `q in {1,2,4,7,14,28}`,

with `q | 2x`, `q | 14z`, and primitivity still required. These conditions follow directly from the full Gram covector. Connected profiles must also be included; solving only this displayed equation is insufficient. A finite explicit upper would allow an exhaustive two-embedding interval certificate with no global-height increase. Without such an upper, a bounded failed search would leave this edge unresolved.

[Allcock's original edgewalking paper](https://arxiv.org/html/2110.03784v1), Introduction and Sections8–9, motivates this finite-corner/binary-plane approach but expressly proves its algorithm for lattices over Z. It is not an off-the-shelf correctness theorem over Z[phi]. A new interval census needs its own complete coefficient bounds. [Allcock2006, Section2](https://msp.org/gt/2006/10-2/gt-v10-n2-p05-p.pdf) supplies the exact silver cutoff already used above; it supplies no H9 compact example. A targeted primary-source search found no exact H3+A6 classification; that limited retrieval is not proof of literature-wide openness.

**Handoff:** first independently check the displayed finite-edge reduction and preflight the literal prime7 host against any exact published classification. Then prepare one finite-profile/interval obligation. Do not rerun the old H30 search or raise its cutoff by default. Any eventual scientific batch must use the existing serial bounded runner, unchanged admission gates,59s/1GiB caps and before/after source pins; no job is proposed as already resource-certified here. This lane may end in another obstruction. At present it is one specific unclosed first-type construction with an explicit arithmetic bottleneck, not a new example.

## Exact inherited bindings

| Input | Bytes | SHA256 |
|---|---:|---|
| [H3+A6 report](fresh_h3_blocks_h9_report_20261013.md) | 8797 | `4545875776e23edf97c0a4a690d8eb0589bf76d3600fe6b7f5daf65081910b19` |
| [H30 certificate](fresh_h3_blocks_h9_a6_phi_h30_20261013.json) | 41893 | `ee0370bc5306d01d0f6954ae0c4759e1ead4194c34190e189ec06d88428b8d5c` |
| [Independent audit](fresh_h3_blocks_h9_independent_audit_20261013.json) | 6113 | `d4acf35f13083ee44419c9e0bc7edd6cda2503d53ca93a7fbc07a38a5ff83863` |
| [Producer source](../scripts/fresh_h3_blocks_h9_a6_search.py) | 6013 | `fcfe23a3beadf58242df05f0c91c0f1c93462b7247bedf2606f717609b4898b7` |
| [Independent verifier source](../scripts/fresh_h3_blocks_h9_independent_verify.py) | 12078 | `06124f9147cb6c70c43928033bc8c57227d80865f1deb3711098976fe2ddc8fd` |

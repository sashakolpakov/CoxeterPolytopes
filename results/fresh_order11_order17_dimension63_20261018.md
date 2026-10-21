# Dimension at most63 for finite-volume hyperbolic Coxeter polytopes

Every finite-sided, finite-volume real hyperbolic Coxeter polytope has
**dimension at most63**. Arithmeticity, compactness, simple cusps and an
extra bound on the number of facets are not hypotheses.

The [source](../scripts/fresh_order11_order17_dimension63_20261018.py)
and [certificate](fresh_order11_order17_dimension63_20261018.json)
combine the new local rules and exact finite/ideal bounds below. The
[completed ceiling64](fresh_deep_finite8_dimension64_20261018.md) already
excludes all dimensions at least65. It remains to exclude dimension64.
This is a complete internally checked computer-assisted proof, with no
identified unresolved mathematical step. It has not been externally
refereed or formalized in a proof assistant.

## A simultaneous complete face-payment theorem

Start from the [complete finite464 payment rule](fresh_ambient64_finite_closure_20261018.md).
It includes the depth11 zero rule for deep finite distance8 pairs, the
finite D16/D17 exchange giving costs116/127, and the affine cycle11
distance5 halves. Its triangles receive at least2 and quadrilaterals at
least1. Add the following compatible changes.

The [order11 tree theorem](fresh_order11_bcd_half_pair_transfer_20261018.md)
halves the two terminal/fork pairs of affine B10, the endpoint pair of
affine C10, and the four opposite-fork pairs of affine D10. Add a finite
half between a whole A1 and the doubled terminal of a separate whole B10.
The ideal family totals are now A99/2, B54, C109/2, D53. The largest is109/2.

The [order17 theorem](fresh_order17_e6_half_coupon_transfer_20261018.md)
lowers four distance10 pairs in the whole affine D16 factor from3/4 to1/2.
Upgrade the finite coupon between an E6 length2-arm endpoint and a D11
short leaf from1/4 to1/2, requiring both components whole. The factor's
cost falls from119 to118. Its companion family totals remain A85/B111/C108.
There are four eligible finite marked pairs per E6/D11 component pair,
so the finite component-pair cost rises from1 to2.

All other finite and ideal weights stay as in the baseline. The new
finite coupons only increase payments. The following actual-corner case
split proves simultaneous quadrilateral coverage and makes the dependency
order explicit.

1. If a cycle11 distance5 corner is present, the independent
   [cycle11 theorem](fresh_cycle11_half_pair_transfer_20261018.md) gives
   its original half plus a distinct retained half, or a finite short unit.
   Its finite outputs have distance at most6; its other ideal outputs have
   distance at most5. None is changed by the new tree or order17 rules.
2. Otherwise, if a changed order11 tree corner is present, apply the
   complete coupled tree theorem at that corner. Its immediate affine
   alternatives contribute a distinct half directly, without recursively
   demanding a full unit from another changed pair. The other alternatives
   are finite whole components of order at most14, whole A9/A18 endpoint
   pairs at distance8, whole D18 leaf9, whole E8 path6, or the stated typed
   finite half coupons. These all survive the depth11 zero and finite
   D16/D17 exchange rules. No alternative uses an order17 distance10 pair.
3. Otherwise, if a changed order17 D corner is present, use its complete
   inherited geometric classification. The64 proper completions remain
   impossible unless the short8 alternative pays. A finite short8 pair
   retains its unit unless it is a deep distance8 zero; in that event the
   constructive depth11 corridor gives a short7 alternative. Its finite
   and ideal replacements avoid the order17 distance10 orbit. The old
   cycle6/7/8 and radius-five tree replacements retain their checked
   depth11 interfaces. In particular every pair at distance at most8 in
   the18 current tree targets of orders12 through17 is a literal unit.
   For cycle arc length8, complementary lengths8 through11 have protected
   nondeep finite outputs, including the sealed cycle16 case; lengths
   at least12 have the strict depth11 cycle-interface gap7/48. Larger
   tree alternatives use retained finite short paths, whole-leaf9 roles
   or separate coupons. If any resulting finite8 output is deep, its
   independent corridor supplies the protected short7 route just described.
   No new order11 corner changes in this branch. Thus the
   short-path route is still sufficient. In each of the four whole-affine
   cases, the19-root radical seals a distinct whole E6/D11 corner; its new
   half plus the original ideal half equals1. No finite same-component
   D16/D17 unit is used.
4. If none of these new ideal corner types occurs, every ideal weight
   agrees with the complete baseline rule, and the only finite changes
   are nonnegative coupon increments. Its old quadrilateral coverage
   therefore remains valid.

The first case also handles a cycle11 corner that was already halved in
the finite464 baseline; applying its standalone theorem again is harmless
and avoids hiding a changed short-path payment inside another proof.
The cases are exhaustive and do not assume that a removed unit survives.

Triangle coverage is inherited from the finite7/ideal5 theorem and the
cycle11 triangle exclusion. The changed tree distances are8,9,10, and the
changed order17 distance is10. Finite paths at most7, ideal paths at most5,
protected ideal6/7 pairs and the small whole affine triangle factors
retain the required payments. The cycle11 theorem excludes its affected
corner on any triangle. Thus every triangle receives at least2 under
the full simultaneous rule.

## The enlarged finite budget is still464

The [complete finite coupon certificate](fresh_e6_b10_half_coupon_forest_budget_20261018.md)
includes both additions: total E6/D11 component-pair cost2, and A1/B10
component-pair cost1/2. It retains all previous coupons, deep11 zeros,
and D16/D17 costs116/127. The new B10 coupon remains dominated by replacing
that finite component with D10; the E6 representative is retained separately.

All finite types and all rank64 forests are covered by753 large-component
vectors,348 coupon environments,38,025 type/partner domination comparisons,
597,168 Bellman upper inequalities and22,272 attaining predecessors.
Every global alternative has a reconstructed attaining forest and literal
coupon recount. The exact maximum is464, attained formally by4D16.
This formal attainment is not an existence assertion about an actual vertex.

The finite curvature coefficient is also464. A finite vertex figure has
H(t)=1+t+...+t^63, V=64 and G=992, giving

    (5G-63V)/2=464.

No positive finite slack is claimed or required.

## Every ideal vertex has positive surplus

The [independent ideal-margin calculation](fresh_ambient64_ideal_margin_closure_20261018.md)
uses the proved order11 envelope109/22 and order17 envelope118/17, with
every other current factor envelope unchanged. The inherited all-rank
concentration theorem reduces total affine rank63 to2,419 profiles:
arbitrary multiplicities of9 through17, at most one1 through8, and at
most one rank at least18. The lower concentration boundary is9 after the
rank10 cost change; the older1,513-profile domain would be insufficient.

For H(t)=product_i(1+t+...+t^r_i), V=H(1) and
G=sum_j[t^j]H(t)*min(j,63-j), every normalized margin

    (5G-63V)/(2V)-sum_i c(r_i)

is positive. Its unique minimum is2385/146432 at(10,11,12,15,15).
The minimum absolute margin is240/23, attained by(63). Both normalized
margin and V decrease under the permitted concentration moves, so the
positive terminal absolute minimum bounds every actual cusp. This covers
all cusp factor counts and ranks, without an arithmetic or simplicity
assumption.

## The strict global contradiction

The completed corner-payment theorem bounds the global defect D5 above by

    D5 <= 464*v_f + sum_(ideal v) V_v*sum_i c(r_vi).

The inherited shape-curvature inequality, with its quantitative two-extrema
strengthening, gives

    D5 >= 464*v_f + sum_(ideal v) (5G_v-63V_v)/2 +5120.

Here5120=5*floor(64²/4). The two actual global extrema and the positive
coefficients of every product-of-simplices vertex figure supply this
constant; there is no new geometric hypothesis. Subtracting the two
bounds and using the absolute ideal surplus gives

    0 >= (240/23)*v_i +5120,

which is impossible even when v_i=0. Dimension64 is therefore excluded.
The completed ceiling64 supplies the entire higher-dimensional tail,
so the resulting unrestricted ceiling is63.

Compact arithmetic remains at most21 and arbitrary compact at most29;
their separate sharper arguments are unchanged. No new finite-volume
example follows from this nonexistence theorem.

## Verification scope

Final generation753c2e and separate fresh replay0bdd94 both terminated
exit0/PASS. All34 direct bindings, Python syntax and nine local report
links pass. The binder reconstructs the five newly changed family interfaces,
checks all short7/short8 roles, including all18 inherited small-tree
short8 targets, and finite coupon amounts, and verifies the
finite464, ideal240/23, constant5120 and inherited higher-tail inputs.
The new ideal source independently recomputes all2,419 scalar rows with
a different enumeration and polynomial implementation. The enlarged
finite budget has exact upper and attainment certificates.

Two separate full source/report mathematical composition reviews passed.
They reread the relevant actual-corner, triangle, deep11 cycle/tree,
finite-budget, factor-fiber, shape and concentration interfaces. The final
review explicitly strengthened the distance8 branch to check all18 current
small tree targets and to state the complete cycle8 and radius-five tree
alternatives; final generation and replay include that change. A separate
saved-only check also confirmed those18 tree pair sets. Neither reviewer
repeated the2,419-row scalar computation, Bellman arithmetic or historical
geometric censuses. The second full ideal-scalar review arrived after its
local report froze and is recorded here without modifying that report.

Geometric classification, whole-component sealing, actual factor fibers,
shape curvature and all-rank concentration remain written theorem inputs,
with their stated internal review scopes. Neither source hashes nor a
numerical solver substitute for those arguments. Historical geometric
censuses are not independently re-enumerated by this composition binder.
No external mathematical refereeing or proof-assistant verification is
claimed. Source, certificate and report are final/frozen.

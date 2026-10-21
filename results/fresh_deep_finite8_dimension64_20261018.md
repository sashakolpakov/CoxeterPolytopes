# Unrestricted finite-volume dimension at most64

**Every finite-sided, finite-volume real hyperbolic Coxeter polytope has
dimension at most64.** The new argument excludes dimension65; the completed
[ceiling65](fresh_ideal5_weighted_dimension65_20261017.md) excludes every
higher dimension. There is no arithmeticity, compactness, simple-cusp or
cusp-factor hypothesis, and no facet-count restriction beyond finite sidedness.

The [source](../scripts/fresh_deep_finite8_dimension64_20261018.py) and
[certificate](fresh_deep_finite8_dimension64_20261018.json) bind the local
and counting theorems and check the current small-factor interfaces.
The written local-to-global composition is given here. Generation, fresh
replay and two separate full internal composition reviews have passed.
This is a complete internally checked computer-assisted proof; the exact
verification scopes are recorded below.

## The single change to the current rule

Retain the complete current face-payment rule, including all ideal weights
and all finite cross-component coupons. Its latest finite refinement is the
[rank16/17 transfer](fresh_finite_d16_d17_long_pair_half_transfer_20261017.md),
and its ideal envelope is the completed
[dimension65 cusp-margin theorem](fresh_ambient65_ideal_margin_closure_20261017.md).

Make exactly one further change. A unit at marked distance8 in a whole
finite A/B/D component becomes a half if both marks have distance at least17
from every endpoint, branch and endpoint of a nonsimple bond of that whole
component. All other weights retain their values. This condition is intrinsic
to the labelled whole component; it makes no choice of an embedding or core.
Only components of rank at least43 can be affected. In particular every
finite short7 unit, every whole small-component payment and every existing
cross-component coupon is retained.

## Triangles and quadrilaterals

The [finite7 triangle theorem](fresh_finite7_ideal5_triangle_payments_20261017.md)
requires only finite distances at most7, ideal distances at most5 or
feature-protected distances6/7, and whole ideal factors of order at most9.
All these payments remain literal units in the current rule and under the
new change. Thus every triangle still receives at least2.

For a quadrilateral, divide into two cases according to its actual corners.
If none has a changed finite pair, every corner weight is exactly its old
value. The completed old local theorem therefore supplies total weight at
least1. This inherits the actual numerical weights, without requiring that
each intermediate branch of its original proof be re-run under new rules.

If at least one actual corner has a changed pair, apply the
[universal deep-corridor theorem](fresh_universal_deep_finite8_corridor_20261018.md)
at that corner. It supplies a path of length at most7 joining the two
opposite black pairs, with only supporting white roots inside the path.
Such a path joins adjacent sides of the actual quadrilateral. Consequently
it is a marked path in another actual corner, with all its supporting roots
present. The statement includes low and high minimal cores, disjoint and
overlapping cores, all native labels and arbitrary spectators. In particular
it does not require the old no-short8 assumption.

If this corner is finite, its marked distance is at most7 and its unchanged
unit pays. If it is ideal, the actual factor-fiber theorem puts the two
marks and the path in the same varying affine factor. The following complete
division treats the possible ideal factors.

* In a factor of order at most11, every pair is still a unit. This includes
  every exceptional affine factor.
* For classical tree factors of orders12 through17, every pair at distance
  at most7 is a literal unit in the current rule. The binder reconstructs
  all18 labelled target graphs and their exact current pair weights, checks
  their total costs against the current envelope, and checks this inclusion.
* For a larger classical tree factor, distance at most5 or distance6/7
  within5 of a feature is a retained elementary ideal unit. In the remaining
  feature-free case the completed
  [radius-five tree transfer](fresh_finite8_ideal5_local_replacement_20261017.md)
  supplies a finite path of length at most7, an ideal path with a persistent
  feature one step from a mark, or the sealed whole finite A7 equality.
  Each output is unchanged. None is a deep finite-distance8 pair.
* For a cycle, distances at most5 are units. Missing short7 pairs have
  distance6 or7. The completed cycle12/13/14 transfers and the universal
  larger-cycle transfer cover these pairs. Their compatibility with the
  present finite halves is precisely the
  [deep finite8 cycle interface](fresh_deep_finite8_cycle_interface_20261017.md).
  The small-cycle distance8 outputs retain a whole-component leaf at
  distance5 or6 from a mark and cannot be deep. Sealed small-component and
  separate-coupon outputs also retain their values. In the sole newly
  residual universal family, an8-edge arc and complementary arc at least18,
  a forced actual spectator neighbor gives the strict determinant gap29/144.
  Thus that family cannot lose the required total payment.

The binder checks all six cycle targets of orders12 through17 as well as
the18 tree targets. For the order13 C factor it uses the73-pair set explicitly
adopted in the completed
[dimension71 proof](fresh_ideal5_weighted_dimension71_20261017.md), whose five
omitted pairs have distance at least10. The earlier support diagnostic alone
is not treated as a local theorem. All other current small-factor sets are
taken from their completed local successors. Their exact costs are:

|Factor order|A|B|C|D|
|---:|---:|---:|---:|---:|
|12|60|64|63|62|
|13|65|73|73|74|
|14|70|84|84|82|
|15|75|94|91|94|
|16|80|107|103|106|
|17|85|111|108|119|

These cases exhaust the affine classification. Each either supplies an
unchanged full unit or a complete retained transfer totaling at least1.
The original changed half need not be added to this payment; no unit is
counted twice. Thus every quadrilateral still receives at least1 under the
simultaneous current ideal rule and new finite rule. There is no obligation
to repeat the later selected-map censuses: faces with unchanged actual
weights inherit their old coverage, and every face with a changed actual
weight has just been handled by the universal short7 alternative.

## The dimension65 contradiction

The [exact conditional finite budget](fresh_depth17_rank65_conditional_budget_20261018.md)
now applies, because its local hypothesis has been proved above. Its maximum
over all rank65 finite Coxeter forests is479, attained by D17+3D16. If a
forest has a component of rank at least18, that component participates in
no current cross-coupon. Strip it off and use the old unrestricted maximum
for the remainder; all48 resulting bounds are at most953/2. The all-small
branch attains479. This argument covers arbitrarily many components and
does not assume the remainder has only small components.

All ideal weights are unchanged, so the complete dimension65 cusp-margin
theorem applies without a new optimization or partition census. It includes
arbitrary products of Euclidean simplex factors in a nonsimple cusp. For
each actual cusp let C_v be its current corner-payment cost and kappa_v its
shape-curvature contribution. Its exact all-profile consequence is

    kappa_v - C_v >=445/23.

Let v_f and v_i count finite and ideal vertices, and let D5 be the sum of
5 minus the number of sides over all two-faces. Nonnegative payments and
the proved triangle/quadrilateral coverage give

    D5 <=2T+Q <=479*v_f + sum_(ideal v) C_v.

The inherited strengthened shape inequality at dimension65 gives

    D5 >=480*v_f + sum_(ideal v) kappa_v +5280.

Combining them yields the impossible inequality

    0 >= v_f +(445/23)*v_i +5280.

Hence dimension65 is impossible. The completed ceiling65 already excludes
all dimensions at least66, proving the ceiling64. In particular the result
is not a dimension65-only obstruction with an unproved higher-dimensional
tail. The separate compact arithmetic ceiling21 and arbitrary compact
ceiling29 are unchanged. No new finite-volume example is claimed.

## Verification scope

The universal corridor and cycle-interface packets have completed exact
generation, fresh replay and separate internal mathematical reads. Their
precise scopes, classification inputs and independent-checker limitations
are recorded in their reports. The new finite budget has likewise completed
generation, replay and a separate full mathematical read. The1,650 cusp
comparisons and their complete reduction domain are inherited from the
frozen ideal-margin theorem and its independent reconstruction. The existing
ceiling65 and its independent reconstruction supply the higher tail.

The new binder checks direct source bindings, all24 current small classical
factor interfaces, the local scalar interfaces and the final strict
inequality. It imports no project arithmetic module. It does not re-enumerate
the historical geometric censuses or re-run their finite-budget optimization.
The present report supplies the mathematical case split connecting these
checked inputs. Final generation eda92d and separate fresh replay102fa2
both terminated exit0/PASS on all24 interfaces. Metadata check5b609a passed
all88 direct bindings, Python syntax and all11 local report links; the
conditional budget's10 bindings, syntax and three local links also passed.

Two other agents completed separate full source/report composition reads.
They checked the actual-corner split, full-star short-path interpretation,
all small-factor provenance including the completed C13 adoption, unchanged
triangle inputs, larger-tree outputs, complete small/universal-cycle
interfaces and whole-component coupon protection. They also read the
conditional finite budget, checked its applicability to the simultaneous
rule, and checked the unchanged absolute ideal surplus, extremum constant,
strict contradiction and higher-dimensional tail. Neither found an
unresolved step. Neither repeated the24-row arithmetic, historical
geometric censuses, Bellman calculation or1,650 ideal scalar comparisons.
These are mathematical/source reviews with the inherited computational
checks explicitly retained, not independent implementations of the whole
proof. Source, certificate and report are final/frozen.

All reviews are internal. External mathematical refereeing and proof-assistant
formalization are absent.

    python3 scripts/fresh_deep_finite8_dimension64_20261018.py
    python3 scripts/fresh_deep_finite8_dimension64_20261018.py --verify

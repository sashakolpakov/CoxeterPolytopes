# Exact monotone repair after the raw23/raw25 incidence constraint

The new [forced-quadrilateral inequality](fresh_sqrt5_core23_25_forced_quad_20261017.md)
can be satisfied by raising only the mass of vertex type `A4+3H3`:

\[
 \mu_{5572}:0\longmapsto\frac{3363761}{3000000000}.
\]

Every face mass and every other vertex mass stays unchanged. The repaired
partial-model obstruction ratio is **12.115434365676277**, down from
12.128700376798244 and still strictly above the target12. This is a
constructive exact repair of one witness, not an optimum or a dimension
bound.

The [source](../scripts/fresh_sqrt5_core23_25_incidence_repair_20261017.py)
and [certificate](fresh_sqrt5_core23_25_incidence_repair_20261017.json)
bind the frozen gate, its [internal review](fresh_sqrt5_core23_25_forced_quad_20261017_review.md),
and the preceding [local](fresh_sqrt5_a3_d5_cap_local_cover_update_20261017.md)
and [joint](fresh_sqrt5_a3_d5_cap_joint_ports_20261017.md) certificates.

The new inequality is

\[
 2\alpha_{26410}+\alpha_{26880}-2\alpha_{580}\le3\mu_{5572}.
\]

Its left side is `3363761/1000000000`; the displayed increment makes it
an equality. The objective numerator remains `3842394855727/309375000000`.
The normalizer is the old sum of vertex masses plus the increment. The
complete exact rational ratio is saved in the JSON, and its comparison
with12 is performed with fractions.

All preceding certificates remain valid without another large scan:

- Every inherited zero-cost signed constraint and all2056 cap checks
  depend only on the unchanged face masses, so their loads are identical.
- All4169 complete local simplex-block marginals are unchanged. The saved
  3712 literal disjoint-packing mixture terms are retained, and their
  required masses have weakly larger vertex allowances.
- The2,711,432 inequalities in the17 joint3/4/5-port families have
  unchanged left sides and nonnegative multiples of vertex masses on
  their right sides.
- The older positive-cost columns likewise have unchanged loads and
  weakly larger nonnegative-cost budgets.

The script checks that the preceding local and joint packages have
identical face and vertex masses, that the target forest belongs to the
existing6558-forest budget universe, and that5572 is the only changed
vertex coordinate. It checks the new equality, the unchanged objective,
the exact sum of masses, and the resulting ratio. The preservation of the
large inherited domains follows from the explicit sign argument above;
the program does not claim to have rerun their arithmetic enumerations.

This repair includes one new geometric incidence constraint. The broader
face/incidence coverage remains incomplete. No optimizer, new global
solve, or dimension theorem is part of this package.

**FINAL/FROZEN:** generation and fresh saved replay passed. All direct
input hashes, local report links, and source syntax were checked. These
small exact checks are distinct from the inherited large computations.

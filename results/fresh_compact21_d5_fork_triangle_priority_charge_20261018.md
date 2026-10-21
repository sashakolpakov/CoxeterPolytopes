# All-face corollary: one unit saved at every E8-triangle-active D5+2E8 vertex

18 October 2026. **A genuine conditional vertex-charge improvement.**
The restriction to quadrilaterals in the subtraction rule of the
[eleven-pattern packet](fresh_compact21_d5_fork_quad_triangle_priority_20261018.md)
can be removed. At every compact arithmetic Q(sqrt5) vertex of type
D5+2E8, subtract its D5 fork-leaf basic unit whenever at least one
designated nonadjacent E8 triangle is present. All existing small-face
inequalities remain valid. The total charge at that vertex decreases
by exactly one. No uniform dimension-bound improvement is claimed.

The required face inequality is sum of corner weights at least 5-k
for a two-face with k vertices. In particular a triangle requires two,
a quadrilateral requires one, and every face with k>=5 needs only
nonnegative payment. This is the existing small-face criterion, stated
in the final weighting paragraph of
[the reviewed compact weighting](fresh_compact_global_bound_20261015.md).

Here is the complete argument. The full D5 triangle domain has no
triangle on the fork pair (3,4), by the
[reviewed complete D5 analysis](review_compact21_d5_two_e8_triangle_coupling_20261018.md).
Thus no triangle corner is modified. In particular every old triangle
bonus, including both H4 bonuses and any optional I2(5)/H3 bonuses,
is exactly unchanged; none is silently replaced by a different route.
For quadrilaterals the eleven-pattern theorem proves total remaining
payment at least one, including when several fork units are subtracted
on the same face. For faces with at least five sides, subtraction of
one basic unit leaves every corner nonnegative, so their required
inequality still holds. All other face weights and every nonnegative
bonus are retained. This proves the extension.

If w_0 is the old charge total at an actual D5+2E8 vertex v, I_E(v)
indicates the presence of at least one designated E8 triangle, and
t_E(v) is the number of those triangles, then the new total is exactly

  w_1(v)=w_0(v)-I_E(v).

The complete E8 triangle-domain theorem gives 0<=t_E<=8, so also

  w_1(v)<=w_0(v)-t_E(v)/8.

This is an actual coupled triangle/quad charge inequality. An additional
weight 1/8 may therefore be placed on every such E8 triangle corner
without increasing the old vertex upper bound: the fork subtraction
funds at least that total. That optional redistribution only increases
triangle payments, and cannot damage a face inequality.

It does not follow that this redistribution changes a particular old
dual. The saved E7/E8-context model has unpaid quadrilateral rows and
an unchanged separately proved triangle baseline. A new positive
triangle payment does not automatically occur in those old rows.
Likewise, an active/inactive refinement may allow the old dual to put
all its relevant mass in the inactive state. A quantitative successor
must connect this state to the actual face rows or prove a complementary
saving when no E8 triangle is present. No old one-row diagnostic is
being interpreted as an obstruction to richer face-dependent rules.

This corollary uses the prior exact certificates and the written
eleven-pattern proof; it needs no new enumeration. It extends that
packet's deliberately conservative quad-only definition, while
preserving its exact source and certificate. Independent review is
pending. No Lean process, optimizer or CPU search was run for this
corollary, and no existing certificate or frozen report was edited.

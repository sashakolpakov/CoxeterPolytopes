# Principal cores extend the deep-corridor alternative to164 corners

The [source](../scripts/fresh_large_finite_corridor_principal_obstruction_20261017.py)
extends the [single-block tail lemma](fresh_large_finite_corridor_tail_obstruction_20261017.md)
to corners with old spectator components. At164 corners in134 of the
[251 remaining supporting contexts](fresh_distance7_direct_contact_features_20261017.md),
a whole finite A/B/D component containing the two marks at distance8,
with both marks at feature distance at least20, forces an opposite-black
path of length at most7 somewhere in the full actual supporting star.

This is again a short-path alternative, not unconditional geometric
impossibility. It excludes the deep corridor only in the branch where the
full supporting star has no such short path. No pair weight is changed.

## Drop spectator roots from the tested Gram

Let P be the displayed connected simple path containing both marks, and
let b and d be the two omitted black roots. The tested principal core is

    P together with b and d.

Every other displayed corner component is dropped from this matrix.
The source selects only cases where this principal core is in physical
norm2 gauge and has exactly one negative eigenvalue and no zero eigenvalue.
Its actual root span is therefore a nondegenerate Lorentz subspace, and
its orthogonal complement in the ambient Lorentz space is positive definite.

An actual white root extending an endpoint of P may simultaneously join
an old spectator path. The immediate first root cannot itself be an old
spectator, since its existing contact to P would put it in the same displayed
block. Later tail roots may be old spectator whites. Neither joining a
spectator block nor this later occurrence is excluded or assumed away.
Each tested tail root is outside the
reduced principal core, and its products with dropped roots do not enter
the necessary Gram condition being tested. This is the reason the argument
does not require the old single-block guard.

In the assumed deep actual corridor, the only contacts from its first
outward tail root to P are its simple contact with the path endpoint.
No later root on that outward tail meets P. The omitted roots b,d are not
in the actual corner, so no tail vertex can be one of them. Thus the
restricted contact profiles are complete even when the full supporting
root contacts dropped old components. The finite A/B/D tree fixes the
tail's internal unit bonds and the orthogonality of the opposite tails.

## Apply the same all-label tail constraints

Both marks having feature distance at least20 supplies an induced A20
white arm outward from each mark. All its vertices are supporting whites,
whether old or newly introduced. Each omitted black belongs to a strictly
positive definite actual edge with those whites. The preceding lemma's
positive-inverse argument therefore applies without change:

* At outward distances3 through18, no omitted-black contact is possible.
* At distance2, any omitted-black contact is simple.
* Contact to both omitted blacks by one white already gives a short2
  opposite-black alternative.

The selected displayed endpoint distances are1 through7, so continuing
to outward distance17 requires only one possibly nonzero restricted contact
vector at the first root of each tail. With principal inverse J and first
contact vectors x,y, the same necessary two-by-two Schur matrix is

    [ (m_left+1)/m_left − x^T J x          −x^T J y             ]
    [             −x^T J y       (m_right+1)/m_right − y^T J y ].

All remaining tail roots are orthogonal to the reduced principal core.
Their possibly nonzero contacts to dropped spectator roots impose further
conditions, but cannot repair a failure of this necessary positive-
semidefinite Gram matrix. The checker treats every permitted unit entry
profile and either exhibits a full-star short7 alternative or a strictly
negative principal Schur minor.

## Complete domain and remaining gaps

All1,004 corners of the251 residual contexts are classified. The164 selected
corners include all64 from the preceding packet and100 additional corners.
Their old actual corner diagrams have one block in64 cases, two in80, and
three in20. Four selected cases have a radical in the full old core while
their reduced principal core is nonsingular Lorentzian; the proof uses only
the latter.

|First applicable classification|Corners|
|---|---:|
|Marks disconnected in the displayed corner|522|
|Marked component is not a simple path|50|
|Displayed marked distance is not8|72|
|A mark is an endpoint of the displayed path|188|
|Reduced principal core is not nonsingular Lorentzian|8|
|Eligible reduced principal core|164|

The540 complete joint contact profiles split into247 negative Schur
diagonals,78 negative determinants with nonnegative diagonals, and215
short7 alternatives. Shortcut profiles are not claimed to fail the Gram
test. Every selected corner has only these outcomes.

The precise excluded-domain inventory and every principal matrix, inverse,
contact profile and Schur outcome are saved in the
[certificate](fresh_large_finite_corridor_principal_obstruction_20261017.json).
In particular,188 distance8 simple-path corners have an endpoint mark,
outside the present positive endpoint-distance hypothesis. Another8 have
principal signatures outside the nonsingular Lorentz case. These are
explicit unresolved obligations, not implied impossible configurations.

The same larger gaps remain as in the single-block packet: disconnected
marks can join through extra supporting whites; earlier structural outputs
still require actual retained pair maps; and low/disjoint and later
universal payment interfaces are not settled here. A short-path alternative
for one corner type does not resolve its entire parent context. This
package provides no new finite budget or dimension bound. The general
ceiling remains65.

## Verification

Generation78056 terminated exit0/PASS (terminal chunkb3be9a). Separate fresh
saved replay3773 terminated exit0/PASS (terminal chunk06508e).
This checker reconstructs native principal cores and all applicable exact
tail comparisons. It inherits the complete supporting census and the A20
all-label lemma; neither the parent census nor a global finite budget is
rerun.

Root and the construction agent each completed a full source/report
mathematical read. They checked the reduced nondegenerate Lorentz span,
the necessary residual PSD condition despite arbitrary dropped-root
contacts, the treatment of old spectator whites later in the tails,
the inherited A20 all-label argument, inclusion of all64 earlier cases,
and the explicit full-star short7 alternative. The clarification about
the immediate first root versus later old spectator roots above records
their common prose correction. Neither review repeated the164-case/540-
profile arithmetic or the inherited supporting census.

    python3 scripts/fresh_large_finite_corridor_principal_obstruction_20261017.py
    python3 scripts/fresh_large_finite_corridor_principal_obstruction_20261017.py --verify

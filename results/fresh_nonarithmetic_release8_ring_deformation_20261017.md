# Continuous ring closure and the cost of releasing an old wall

The missing central Coxeter angle can be reached numerically by changing
the ring geometry. The complete deletion ledger prevents this from being
an improvement: the two valid wall selections have18 or15 open ends,
against15 in the frozen control. No exact chamber or finite-volume example
is established.

The [source](../scripts/fresh_nonarithmetic_release8_ring_deformation_20261017.py)
and [saved result](fresh_nonarithmetic_release8_ring_deformation_20261017.json)
continue the closest [two-center relaxation](fresh_nonarithmetic_two_central_caps_20261017.md),
parent stage25. That relaxation had mutual central product approximately
-0.000638488 and eleven open ends, but the central angle was not Coxeter.
This successor changes a geometric parameter, rather than increasing the
old finite angle alphabet.

## A free scaffold contact

Keep the scaffold wall12's native products with7 and52 fixed at0 and
-sqrt(2), and let its product with old8 be a real parameter t. The norm-two
condition leaves two branches; choose the one through the saved wall.
Rebuild L,P,Q,R,S,T from their original three-contact recipes, retaining
their A2 profiles and prescribed finite labels. The same continuous branch
selection at each step reproduces the original geometry at t=0 to within
2.4e-11 in the numerical coordinates.

The numerical mutual product S.T changes sign between t=0 and t=0.0001.
A scalar root solve gives

`t = 0.000012189912078513614`,

with S.T approximately -4.55e-13. Every other new-new product passes the
numerical native Coxeter/dotted test, as do all retained old-wall products
in each of the two selections below. The inherited common center remains
strictly inside every new half-space.

This is a floating-point root and feasibility diagnostic. Neither the
existence and uniqueness of this root nor its complete Gram matrix has
been certified by interval arithmetic. Finite values are compared using
the inherited numerical tolerances. The root is not promoted to an exact
right angle by rounding.

## Full wall losses change the conclusion

Retaining both old8 and the scaffold12 is invalid: their positive product
t is obtuse. The apparently attractive eleven-end score of that selection
cannot be accepted. The two selections passing all numerical native gates
are instead:

| Wall selection | Old finite-support certificates lost | S6 finite vertices | S6 exterior vertices | S6 open ends |
|---|---:|---:|---:|---:|
| Delete8,47,53; replace12 and addL,P,R,Q,S,T |119|30|10|18|
| Delete12,47,53; retain8 and addL,P,R,Q,S,T |4|33|7|15|

Both use67 listed walls. The loss counts exhaust the frozen198 ambient
finite-vertex supports: a support is lost when it contains a removed or
changed old wall. These are exact incidence losses of old certificates,
not a new ambient finite-vertex census. In particular deleting8 is much
more costly than its seven-sided S6 face alone might suggest.

The intrinsic S6 scores exhaust projective triple intersections numerically.
They include the new left-hand gap created by the deletion, rather than
only the central disk. The second selection provides a numerical closed
central disk but does not reduce the complete open-end count.

## Attempts to fill the resulting left cavity

Retain8 and discard the scaffold12. The new cavity boundary is
7,8,52,L. Using all seven inherited A2 profiles and the unchanged finite
alphabet `2..12,16,24,32`:

* A single replacement framed on7,8,52 has13 candidates satisfying every
  retained old-wall gate. All13 fail the native product withL.
* A two-cap disk uses frames7,8,L and8,52,L. Their complete accepted pools
  have2 and38 candidates. All76 mutual products are tested:53 fail the
  finite-angle/dotted test and23 are dotted. None supplies the finite
  central join needed by this disk.
* The first cap of a three-cap variant is framed on7,L,52 and must make
  both L-end corners finite. Its24 real finite-endpoint branches yield
  no candidate satisfying every retained native gate.

These are bounded searches, not all-label or all-parameter obstructions.
The construction direction remains open, but simply tuning the almost-right
central angle does not by itself improve the whole chamber boundary.

Generation and a fresh saved replay terminate with PASS for the stated
numerical diagnostics and exact old-support loss counts. There was one
scalar root solve and no global optimizer. Frozen predecessors and the
64-wall,198-finite-vertex,15-open-end exact control remain unchanged.

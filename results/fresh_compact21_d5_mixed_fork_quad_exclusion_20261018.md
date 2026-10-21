# The D5 priority subtraction needs no fixed half-units

18 October 2026. **Strengthening needed for the actual golden-field
model interface.** The P13-P23 alternative in the eleven-pattern
fork-quad theorem is arithmetically impossible over Q(sqrt5).
Consequently the active-state D5 fork-unit subtraction preserves
a full basic unit on every affected quad, without assuming any
fixed D5 long-short half-units. This applies directly to the old
golden-field baseline, which has only distance-at-most-two and H4
weights. No original frozen proof or certificate is edited here.

The earlier eleven-pattern theorem was conservatively stated for
a weight rule also containing the D5 half-units; its P13-P23 case
was paid by those units. Those half-units appear in the separate
cubic-field rule but are not part of the fixed baseline G.base
in scripts/fresh_sqrt5_forest_sensitive_weights_20261015_review.py,
which the old golden-field model imports. The strengthening below
removes that extra assumption instead of silently adding a new cost.

## Complete exclusion of the remaining pattern

Keep D5 edges 01,12,23,24, W={0,1,2}, old black leaves A=3,B=4,
and cyclic side order A,B,c,d. Exchange c,d and A,B if needed so
c contacts white 2 by label 3 and d contacts white 1 by label 3.
Their other fixed products are zero, including c dot d=0 and the
products with adjacent old black sides. This is P23-P13. Write
the free opposite products as A dot c=x and B dot d=y.

All roots can be normalized to square two over the field along
their fixed label-three paths. The resulting simple Cartan cycles
show x,y in O=Z[phi]. Opposite indefiniteness, justified by acute
projection and compactness as in the previous packets, gives

  x<-1/2,    y<(1-sqrt5)/2.

Eliminating W gives

  [5/4, -3/4, x-3/4, -1/2;
   -3/4,5/4, -3/4,   y-1/2;
   x-3/4,-3/4,5/4,   -1/2;
   -1/2,y-1/2,-1/2,  1].

Its determinant is (2-x)/2 times

  (2x+1)(1+y-y^2)+3y-5.

Corank one therefore gives that polynomial equal to zero, or

  x=(y-2)^2/(2(1+y-y^2)).                            (1)

Both opposite pairs cannot be dotted or parallel. If y<=-2,
the denominator in (1) is negative and

  x+2=(8-3y^2)/(2(1+y-y^2))>0,

so x>-2. At least one pair is thus finite. The complete finite
product list in Q(sqrt5), together with the strict bounds above,
leaves -1 or -phi for that finite opposite product. There is no
unspecified high-label or opposite-distance cutoff.

For finite y=-1, (1) gives x=-9/2; for finite y=-phi it gives
x=-5phi/4. Neither is an algebraic integer in O. For finite
x=-1, the equation is y^2+2y-6=0, with discriminant 28, which
is neither a rational square nor five times one. For finite
x=-phi the equation is

  sqrt5 y^2+(3-sqrt5)y-(5+sqrt5)=0.

Its discriminant 34+14sqrt5 has rational norm 176, not a rational
square, so it is not a square in the field. Neither equation
has a solution y in Q(sqrt5). These four alternatives exhaust
all cases and prove the claimed exclusion.

## Stronger payment statement and model interface

In the complete eleven-pattern list, this removes the only case
which previously used two fixed half-units. The other seven
excluded patterns stay excluded; each of the three surviving
arithmetic patterns forces E8-triangle absence at every fork
corner. Therefore, after subtracting the fork basic unit whenever
an E8 triangle is present, every affected quad retains another
basic unit. No D5 half-unit is needed. This includes simultaneous
subtractions at several corners.

The all-face extension remains unchanged: D5 fork faces cannot
be triangles, so every old H4 and optional triangle bonus route
is unchanged; faces with at least five sides retain nonnegative
weights. Thus the active-state one-unit saving is valid for the
actual old golden-field baseline. This supplies the missing
interface assumption in the two-row successor proposed in the
extreme-E8 consumer packet. No optimizer or dimension theorem
is claimed.

The [checker](../scripts/fresh_compact21_d5_mixed_fork_quad_exclusion_20261018.py)
and [certificate](fresh_compact21_d5_mixed_fork_quad_exclusion_20261018.json)
verify the exact determinant factor, both nonintegral substitutions
and both field-nonsquare discriminants. Generation passed in
`4f543d`, peak RSS 27,213,824 bytes; a separate replay passed in
`ddf86b`, peak RSS 27,590,656 bytes. Both held the shared lock,
used a 59-second alarm and checked RSS below one GiB. No search,
optimizer, temporary output or Lean process was used. Independent
review of this strengthening is pending.

```sh
python3 -B scripts/fresh_compact21_d5_mixed_fork_quad_exclusion_20261018.py --verify
```

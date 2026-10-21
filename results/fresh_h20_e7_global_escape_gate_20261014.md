# Every added E7-face facet retains the global escape

**Exact scope:** retain the current 87 walls. Any valid new Coxeter wall that
produces an actual facet on the selected E7 face leaves the saved global
spacelike escape inside its halfspace. Consequently a finite collection
consisting only of such additions cannot complete the H20 chamber.

This does not exclude arbitrary completions. A completion could add walls
whose intersections do not form facets of this E7 face, or a construction
could replace old walls. The [complete H13 source](fresh_h13_a11a2_completion_20261014.md)
is therefore being tested in a branch that replaces conflicting old walls.

The [exact certificate](fresh_h20_e7_global_escape_gate_20261014.json) and
[checker](../scripts/fresh_h20_e7_global_escape_gate_20261014.py) use recursive
arithmetic in `Q(sqrt2,sqrt79,sqrt11)`. There are no numerical tolerances,
angle cutoffs, or height bounds in this exclusion.

## Reduction to two lifts

Use the [certified E7 face](fresh_h20_e7_face_20261014.md), its normalized
`A11+A2` simple normals, fundamental weights `lambda_i`, and future timelike
basis vector `T`, of norm `-2`. Write `omega` for the E7 fundamental weight
at old wall 33, of norm `3/2`. The complete necessary contact gate gives,
after normalizing a proposed ambient wall to norm 2,

`r = -sum d_i lambda_i + t T - ell omega`,

where `ell` is either 0 or 1. Its intrinsic projected norm is respectively
2 or `1/2`. Acute Coxeter compatibility with the thirteen face walls gives
`d_i in {0,1} union [sqrt2,infinity)`. This relaxation includes every finite
integer angle label, parallel pair, and divergent pair.

The coefficient `t` is positive for every valid addition. The old finite
corner consists of all twenty `E7+A11+A2` walls. Its positive-definite Gram
inverse has nonnegative entries. Any proposed strict future center has
coordinates `p=uT-sum b_i Omega_i`, with `u>0` and all `b_i>0`; acuteness
gives `r=tT-sum c_i Omega_i`, with all `c_i>=0`. Thus
`r.p=-2tu+c^T G_corner^-1 b` could not be negative if `t<=0`.
The checker verifies the full corner Gram inverse, orthogonality, and
future orientation. Preserving the old apex is a consequence, not an
additional restriction.

## The right lift

Let `E` be the saved global escape, and set
`a_i=-lambda_i.E`, `B=-T.E`, `w=-omega.E`, all positive.
The pairing is `r.E=a.d-Bt+ell*w`.

For `ell=0`, a nonnegative rational combination of the inequalities from
old walls 30, 75, and 76 bounds `r.E` by a linear form whose fourteen
coefficients in `(d,t)` are strictly negative. The multipliers and every
exact residual are saved. Hence every acute right lift retains `E`, even
before imposing its norm or discrete angle restrictions.

## The nonright lift

Suppose `ell=1` and `r.E>0`. Then `t<(a.d+w)/B`. Substitute this upper bound
in the old-wall 16 and 76 inequalities, whose timelike coefficients have
the required sign. The resulting necessary inequalities are linear in `d`.

1. Their combination with rational weights `194/1000` and `843/1000` has
   positive coefficients and forces `d_12<1`. The Coxeter gap gives `d_12=0`.
2. With `d_12=0`, the wall-16 inequality forces `d_0>1`. The next Coxeter
   gap gives `d_0>=sqrt2`.
3. Put `J=diag(C_A11,C_A2)^-1`. The function
   `f(d)=d^T J d-2(a.d+w)^2/B^2` is strictly convex: the checker proves
   `B^2-2a^T J^-1 a>0`.

On the relaxed domain `d>=0`, `d_12=0`, `d_0>=sqrt2`, its minimum is attained
at `d_0=sqrt2`, `d_11=3*a_11*C/D`, with all other entries zero, where
`C=a_0*sqrt2+w` and `D=B^2-3*a_11^2`. Exact nonnegative gradient conditions
certify this minimum, which equals `11/6-2C^2/D` and is strictly greater
than `1/2` (approximately `0.50381449`). But `r.E>0` and `t>0` imply that
the intrinsic norm exceeds `f(d)`, contradicting its required value `1/2`.

This proves that the nonright lift also retains `E`. The original 87
halfspaces and any finite collection of these additions still contain the
spacelike escape; a small perturbation toward their strict future center
makes the inequalities strict while preserving positive norm. The chamber
therefore retains infinite hyperbolic volume.

Replay: `python3 scripts/fresh_h20_e7_global_escape_gate_20261014.py --verify`.

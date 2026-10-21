# Golden `A1^24` uniform ray: exact first cutter

The hyperideal extreme ray

\[
u=(1;(\varphi/4)^{24}),\qquad (u,u)=3+\varphi,
\]

survives the first two wall packets, but it does **not** persist globally.
The first admissible control at which it is cut is
`x=phi^3=1+2phi`.  Choose an octad `O`, put `C=O^c`, choose `p in O`, and set

\[
d_p=\varphi,\qquad d_i={\varphi^2\over2}\ (i\in C),
\qquad d_i=0\quad(i\in O\setminus\{p\}).
\]

Then `r=(x;d)` is a primitive norm-two lattice root.  Both parity supports of
`2d` are the weight-16 Golay word `C`, while the time coordinate is a unit.
The norm equation is

\[
\sum_i d_i^2=\varphi^2+16{\varphi^4\over4}
=9+13\varphi=1+\varphi x^2.
\]

An exact two-place enumeration shows that the only admissible controls up to
this height are

```text
1, phi, phi^2, 2phi, phi^3.
```

At `phi^2`, a complete unordered coefficient census has 43 square patterns,
of which three satisfy the Golay support condition.  Two acute orbits are
accepted: 6,072 roots of type `(phi^2,phi,phi)` and 12,144 of type
`(phi^2,(phi/2)^8)` with the half coordinates on an octad.  The third orbit,
of size 6,072, pairs by `+phi` with its earlier octad wall and is rejected.
At `2phi`, the sole Golay-feasible type is `(phi^3)`; all 24 assignments pair
by `+2phi` with an `x=1` wall, so that control contributes no walls.

For the displayed cutter, the exact maximum pairings with all earlier wall
packets are:

| packet | maximum pairing |
|---|---:|
| height zero | `0` |
| `x=1` | `-2-4phi` |
| `x=phi` | `-2-2phi` |
| accepted integral `x=phi^2` orbit | `-6-8phi` |
| accepted octad `x=phi^2` orbit | `-4-4phi` |

Thus the root survives every earlier wall.  It cuts `u` strictly because

\[
\sum_i d_i-4x=4+\varphi>0,
\qquad
(u,r)={\varphi\over2}(4+\varphi)
=\frac12+\frac52\varphi>0.
\]

The complete multiplicity census at `phi^3` has 383 square patterns, 33 with
Golay-feasible parity supports, and exactly 16 with positive uniform excess.
The displayed root is the particularly simple pattern
`phi^1 + (phi^2/2)^16`.

This closes only the fate of the proposed uniform ray.  It does not enumerate
the full accepted `phi^3` packet or complete the later Vinberg chamber.

- Audit: `scripts/audit_golden_a1_uniform_ray_cutter.py`
- Tests: `tests/test_golden_a1_uniform_ray_cutter.py`
- Archive: `data/enumerations/golden_a1_uniform_ray_cutter.json`

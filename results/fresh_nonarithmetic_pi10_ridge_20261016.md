# The new pi/10 ridges and a stronger escape filter

Every nonzero acute-compatible normal perpendicular to old wall33 pairs **strictly negatively** with the renewed escape of the [affine-free52-wall chamber](fresh_nonarithmetic_affine_filtered_20261016.md). This is an all-real linear-cone obstruction, independent of any norm choice, field, finite-angle alphabet or finite-vertex requirement. In particular, walls perpendicular to the three new pi/10 ridges cannot cut that escape.

The [exact certificate](fresh_nonarithmetic_pi10_ridge_20261016_gate.json) is generated and replayed by [this checker](../scripts/fresh_nonarithmetic_pi10_ridge_20261016_gate.py). It uses the ten exact normal columns

```
M = [n0,n10,n11,n13,n19,n27,n30,n33,n49,n51].
```

Rational outward interval Gaussian elimination certifies that every pivot is nonzero. Thus `c=M^-1 E` is an exactly defined algebraic vector and `E=M*c` is an exact identity. The checker proves `c33<0` and that all nine other coefficients are strictly positive. The interval residual is also checked as a sanity check; the exact identity follows from the invertible linear system, not from a small numerical residual.

If a normal `r` satisfies every retained acuteness inequality and `(r,n33)=0`, then

```
(r,E) = sum_(i != 33) c_i (r,n_i) <= 0.
```

Equality would force orthogonality to all ten independent basis normals, hence `r=0`. The proposed added wall need not preserve the saved center. This theorem only uses its pairwise acuteness and perpendicularity to33.

This does not exclude finite-volume completion of the52-wall chamber. It gives a necessary distinction between tasks: new walls perpendicular to the pi/10 ridges may be needed to create their required finite vertices, but at least one escape-cutting wall must have a nonright contact to33. It may meet33 at a nonright finite angle, be parallel, or be dotted. A finite-volume completion could use both kinds of wall.

The accompanying [bounded numerical endpoint scout](fresh_nonarithmetic_pi10_ridge_20261016_scout.json) was the diagnostic leading to the identity. For the ridge `[33,49]`, the six old common-orthogonal roots are `[0,1,2,3,5,6]`, of type `H4+A2`. It tests normals right to the I2(10) pair, with the seven possible one-root spherical profiles on this H4+A2 support: all right, or a single label3,4,5 to one A2 endpoint. H4 cannot acquire a new connected node in a spherical extension. Both conic branches are tested against each of the44 outside anchors at labels2 through20. All52 old contacts and the renewed escape are scored numerically.

There are11,704 real branches:1,588 fail the saved-center orientation,9,639 have an obtuse old contact, and477 are acute but noncutting. No candidate survives to exact cap construction. This finite pilot is not an all-label rejection or a complete search of normals perpendicular to the ridge. The separately proved linear-cone gate is stronger and has no such cutoff.

The numerical scout uses midpoint input and80-digit arithmetic, with an explicit numerical-only scope. The exact gate uses the frozen algebraic normal recipes and rational outward intervals. Both scripts passed generation and saved `--verify` replay, source-hash checks, syntax and links. Independent review, if supplied, is separate from these producer replays.

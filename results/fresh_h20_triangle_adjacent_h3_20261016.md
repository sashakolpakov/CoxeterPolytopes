# Numerical adjacent H3 faces of the compact triangle

The [source](../scripts/fresh_h20_triangle_adjacent_h3_20261016.py) and
[saved diagnostic](fresh_h20_triangle_adjacent_h3_20261016.json) reconstruct
the residual Lorentz5 geometry of the frozen
[153-wall chamber](fresh_h20_pi8_triangle_face_20261016.md). They restrict
all153 old inequalities to each of the two supports S16+150 and S16+151.
Generation and a fresh tolerance replay pass.

The new compact triangle is their common boundary face. This diagnostic
uses floating point halfspace intersections, a 1e-6 band around the unit
sphere, and deterministic Sobol samples. It is not an exact complete
cone census or a volume certificate.

| Retained additional support | Vertices | Finite | Near sphere | Exterior | Maximum squared radius |
|---|---:|---:|---:|---:|---:|
| 150 | 9 | 3 | 1 | 5 | 1.063870706 |
| 151 | 14 | 3 | 0 | 11 | 1.007803605 |

The affine three-dimensional Euclidean volumes in their respective
orthonormal charts are approximately .00321607 and .000585848. The
sampled portions outside the unit ball are .000922244 and .000139093.
These are chart diagnostics, not hyperbolic volumes.

The projection of the remaining E6 onto each support complement has
squared norm approximately -.07367654 and +.000626488, respectively.
A linear-program test asks whether that projection is a nonnegative
combination of the16 known residual finite vertices, allowing an
arbitrary multiple of the additional support normal. Both numerical
programs report infeasible. Thus this particular finite-vertex cone
test supplies no obstruction to an E6 cutter right to either support.
It does not prove that a Coxeter cutter exists, nor exclude other
positive-cone certificates.

The whole residual face has a much larger exterior region in the wall142
direction. That is addressed by the separate exact
[dominant-axis cap](fresh_h20_dominant142_cap_20261016.md). The small adjacent
H3 faces remain concrete later closure targets; no exhaustive partner
search or finite-volume claim is made here.

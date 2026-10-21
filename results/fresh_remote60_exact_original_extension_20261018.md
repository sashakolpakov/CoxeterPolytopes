# Exact remote cap on the original 63-wall chamber

This packet constructs an exact norm-two wall U75 on the original H9 chamber,
retaining **all 63 original facets and all 199 original finite vertices**.
The resulting chamber has 64 facets and exactly **201 finite vertices**.
In the three-dimensional hyperbolic section orthogonal to the fixed six
white roots, the whole old60 face is now a compact triangle with angles
π/9, π/2 and π/3. A saved strict spacelike escape survives the new wall, so
this particular 64-wall chamber has **infinite volume**.

This is an exact extension of the original scaffold, independent of the
later numerical J/L/X/N/R/S/T construction. It does not establish compactness
of the full ambient facet60, a finite-volume example, or a dimension bound.

The producer is
[fresh_remote60_exact_original_extension_20261018.py](../scripts/fresh_remote60_exact_original_extension_20261018.py),
the certificate is
[the matching JSON](fresh_remote60_exact_original_extension_20261018.json),
and the standalone saved checker is
[check_remote60_exact_original_extension_saved_20261018.py](../scripts/check_remote60_exact_original_extension_saved_20261018.py).

## Exact input and the changed construction

The physical input is the frozen
[63-wall, 199-finite-vertex original chamber](fresh_nonarithmetic_release49_62_replacement_20261017.md).
Its norm-two normals, Lorentz metric, strict common center h, strict facet
witnesses and complete finite census are inherited. The six fixed white
walls are W={0,1,2,3,5,6}. Their Gram matrix is positive definite; W-perp
has Lorentz signature (3,1). No original wall is deleted in this packet.

The numerical
[remote bridge discovery](fresh_remote60_compact_face_bridge_20261018.md)
suggested U with white profile zero and frame

```
frame       [60,58,59]
labels      [ 2, 2, 3]
products    [ 0, 0,-1].
```

Only this one recipe is promoted here. Its numerical coordinates, N/S/T
normals, 420 numerical joins and numerical cone enumeration are not premises
of the exact native-contact proof.

The original 58–59 bond has label9. Thus this calculation is not claimed to
take place solely over Q(sqrt5). The exact construction uses the original
physical roots, cos(π/9), the white projection, rational operations and one
positive square root. Directed rational bounds enclose those quantities.

Most original normal coordinates are reconstructed from the saved facet
witness identity `p_i=h−(h·r_i)r_i/2`. For original wall66, the saved direct
normal enclosure is used: reconstruction through its coarse center product
would be unnecessarily wide. The two resulting formulas for each U product,
the Schur formula and direct physical-coordinate pairing, are intersected.
Both enclose the same exact product; no value is selected from a numerical
approximation. In particular the final U–66 enclosure lies strictly below−2
(approximately [-5.35046,-4.86425]).

Original Coxeter products isolating zero or a single finite label are
intersected with that exact identity. The neighbor-label inequalities are
checked using directed enclosures; floating point only proposes a label to
check. Dotted products retain their interval enclosures. This is the same
original-scaffold identity principle as the frozen
[physical continuous-N certificate](fresh_oblique_low_s_continuous_n_interval_20261018.md).
It is not applied to an unknown new product.

## Defining U and certifying both sheets

Write bars for projection off W. Let F=(bar(r60),r58,r59), G its three-by-three
Gram matrix, b=(0,0,-1), and c=(bar(h)·F_i). The source checks every positive
Gaussian pivot of G. Its exact shape is

```
[(3−sqrt5)/3,           0,              0]
[           0,           2, −2cos(π/9)]
[           0, −2cos(π/9),              2].
```

Define

```
u0 = sum_i (b G^-1)_i F_i
h0 = bar(h) − sum_i (c G^-1)_i F_i
q  = h0²
rho = 2 − b G^-1 b^T
U  = u0 + sqrt(rho/q) h0.
```

The interval proof gives q<0 and rho<0. Thus the positive square root is real,
U²=2 algebraically, U is orthogonal to W, and its three prescribed products
are exact. The complementary subspace is one-dimensional and timelike, so
these equations have exactly the displayed root and the opposite sign of
the square root. The displayed root has h·U<0. The other has h·U>0 and is
rejected by the common-center orientation. The direct coordinate norm
encloses2 as a consistency check, not as the proof of the norm identity.

All 63 original-wall products are checked:

* eight are exactly zero: W and walls60,58;
* wall59 has product exactly−1;
* the other54 have product strictly below−2.

This includes walls30,44,63,64,65, which had been deleted in the newer
numerical construction chain. There is therefore no 19/12/1 deletion debt
here. The old63-wall chamber is the starting point itself.

## All 64 facets and all old finite vertices persist

For each old root r_i, its old witness
`p_i=h−(h·r_i)r_i/2` has strictly negative product with U: h·U<0,
h·r_i<0 and U·r_i≤0. Its old-wall properties are inherited. Conversely
`p_U=h−(h·U)U/2` is timelike, lies on U and has strictly negative products
with every old root. Thus U is an actual new facet, and none of the old
facets becomes redundant. The saved certificate checks all these strict
products and the new witness square directly.

For any original finite vertex with finite support I, let G_I be its
positive-definite Coxeter Gram matrix. A positive-definite matrix with
nonpositive off-diagonal entries has entrywise nonnegative inverse.
Consequently `lambda=G_I^-1(h·r_I)` is entrywise nonpositive. The finite
vertex ray

```
v=h−sum_i lambda_i r_i
```

has `U·v=h·U−sum_i lambda_i(U·r_i)<0`. This proves persistence of every
one of the original199 finite vertices, using their inherited complete
finite-support census. There is no replay of 199 large Gram inversions.

This argument also explains why the present extension differs from the
earlier sequence of wall replacements: there is no lost original finite
support to compensate.

## Complete new finite census and the section triangle

A new finite vertex must contain U. Its other active walls cannot be dotted
to U, so they belong to the nine-element set W∪{60,58,59}.

Finite vertices are simple here. Indeed their active roots lie in a positive
Euclidean space and have nonpositive off-diagonal products. If their Gram
matrix were singular, an irreducible singular component would have a
nonnegative kernel vector. The corresponding nonnegative root relation
contradicts the strict negative products with h. The active roots are
therefore independent, and a vertex in H9 has exactly nine of them.

There are only nine possible eight-old-wall supports to accompany U.
Seven contain both58 and59. Their principal Gram matrix on {58,59,U} has
strictly negative determinant, so they cannot be finite. The two remaining
supports are

```
W + {60,58,U}
W + {60,59,U}.
```

Both are checked positive definite and give timelike future vertex rays
satisfying every other old-wall inequality strictly. They are new because
all original finite vertices lie strictly inside U. The complete finite
census is therefore199+2=201.

The source also checks the old vertex W+{60,58,59}. Projecting h off each
of the three section supports computes all three actual vertices, their
negative squares and every other wall product. On the old60 plane in W-perp,
the walls58,59,U bound the hyperbolic triangle with angles π/9,π/2,π/3.
Their sum is17π/18<π. All three finite vertices satisfy all remaining
walls, so those walls do not cut this triangle. Hence the **whole old60
face of this H3 section** is compact. This is not a statement about the
eight-dimensional ambient facet60.

## The new chamber still has infinite volume

The original certificate's strict spacelike escape, supported near the
old projected region {7,30,47}, has positive square, future orientation
and strictly negative products with all63 old walls. Its product with U
is also strictly negative, with enclosure approximately
[-0.654475745,-0.654475708]. It remains a strict spacelike point of the
homogeneous chamber cone.

Combining it with the strict timelike center crosses the future light cone
inside all wall inequalities. This gives a relatively open patch at
infinity, hence infinite hyperbolic volume. The compact old60 section
triangle removes one obstruction, but does not close the full chamber.

## Verification and limits

Final producer generation: session85882, terminal64d674, log83641f,
exit0/PASS. A preceding native-only generation passed; the final version
adds the complete nine-support argument and inherited strict escape.
Separate fresh replay session61750, terminal5e8c11, and standalone checker
session87019, terminal325403, both completed exit0/PASS (log446feb).

The standalone checker has no project imports. It reconstructs all original
physical roots, original native isolations, white and frame projections,
the exact two-sheet formula, both physical product expressions, facet
witnesses, all three section vertices, all nine new-support trials and the
surviving escape. Its rational interval implementation follows the same
directed arithmetic design as the prior saved checker; it is a separate
saved verification path, not an unrelated arithmetic algorithm. It inherits
the original scaffold's exact root definitions, signature, native theorem
and complete199-vertex census. No historical construction or numerical
N/S/T cone enumeration is replayed.

The order16 agent completed a full producer/report/standalone-helper
mathematical read (3b4420/33f8d5/3f4082), including the original63-wall
report and physical recipe/facet/census source, inherited interval
primitives and all-label theta enclosure. The review checked original-only
native isolation, the two expressions for the same exact U and its unique
center-facing sheet, all64 strict facet witnesses, persistence of199 finite
vertices, the finite active-root independence argument and complete nine
new supports, all three section vertices, and the strict spacelike escape.
The verdict was PASS with no correction. The reviewer did not execute the
producer/helper, replay interval arithmetic, reconstruct the historical
scaffold/census, or rerun numerical cones. Its earlier conceptual check of
the facet and finite-vertex retention argument is separate from this full
read.

Root also completed a full producer/report/helper mathematical read
(a1c187/a21e4d/f2aa94), checking original-only native identity isolation,
white/frame positivity, the timelike-complement two-sheet formula,
Schur/direct-product intersections, strict old/new facet witnesses,
finite-Gram inverse positivity and the complete nine-support census,
the actual section triangle and surviving strict escape. Its helper read
included directed rational rounding, square roots, theta enclosures and
Gaussian inversion as a separate verification path. The verdict was PASS
with no correction. Root did not execute or replay interval arithmetic,
or reprove the historical scaffold/census. All producer, fresh replay,
saved-checker and source/report reviews are complete; this packet is frozen.

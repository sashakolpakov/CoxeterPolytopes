# Review of the full triple-H3 tetrahedral face

20 October 2026. Verdict: **PASS**, with the scope below. This is a
source-only mathematical review; no worker import/body, computation,
generation, replay, TeX, or Lean job was run.

Reviewed in full:

* `results/golden_dyadic_triple_h3_compact_face_20261020.md`,
  SHA256 `665a9861ae70634b20a80fa93007276fb607e254e84a5c3d637d337adfee46a8`;
* the coordinates in `results/golden_dyadic_sixteen_wall_caps_20261020.md`,
  SHA256 `c288b5de619471772c8381d98afa3996ac8851b4d2d146ad520d6474ee8ba92b`;
* the previously reviewed exact integral splitting and norm4-H3 seal in
  `results/golden_dyadic_d4_fourth_vertex_and_sealed_h3_face_20261020.md`,
  SHA256 `d6b148954b09e259d0e0598cc4f754904ba1977b373a3f491196ef87f93f9501`.

The full actual M6 chamber and its corrected first-edge/reversal are
inherited from the frozen four-edge/M6 packet, rather than re-executed.
Their geometric interface was reviewed previously. I co-derived the
sixteen-wall escape and reviewed its fixed checker; that is disclosed
support, not a claim of independent discovery of the present face.
I independently checked the new seal, scalar bounds and full-face step.

## The actual corner and the nonparabolic seal

At M the six active M6 facets are exactly
`(u,beta1,beta2)` and `(w,beta5,beta4)`. Both have norm2 H3 Gram,
determinant `2phi^-2`, and are orthogonal. The other native M6 products
are negative as stated. Thus the inherited whole M6 identification
makes this an actual ambient triple-H3 vertex. Directly,
`M^2=21-13phi=-phi^-7`; `T0=phi^4 M` is integral with square `-phi`.

The unit norm of T0 splits `M6=O*T0+N6` integrally; N6 is positive
definite and unimodular. Each E or F root module is saturated because
its determinant ideal is the prime ideal `(2)`: the square of a
nontrivial integral index cannot divide it. This uses the actual
integral module, not only its K-space or determinant class.

A proper mirror crossing the E+F perpendicular space gives a finite
rank7 reflection group at a timelike point of intersection. A finite
irreducible component containing these actual 5-angle reflections can
only be H3 or H4; the crystallographic finite groups have no pair of
root mirrors at angle pi/5. The two orthogonal rank3 groups cannot be
contained together in rank4 H4. Therefore a new nonorthogonal component
would have to be H4 containing E and orthogonal to F (or vice versa).
This argument does not assume E is parabolic in the H4.

Every H4 reflection is conjugate within that integral finite group to
one of E's reflections. Its primitive norm can therefore be chosen
exactly2. Its simple-root O-module N is unimodular and contains E:
the norm2 H4 reflection matrices preserve their O-root module.
For x in N, the integral projections `R5 -> M6 -> N6` land in
`N6 intersection F^perp=N6 intersection K E=E`. This is precisely the
orthogonal projection to E, so E splits integrally from N. The
nonunit determinant of E makes that impossible in unimodular N.
The same argument treats the other component. Hence the full seal
does cover nonparabolic H3 inclusions, not just a list of simple
one-node attachments.

Consequently the exact perpendicular lattice is `O*T0+P3` and splits
integrally from R5. Every intrinsic reflection extends over its
orthogonal complement, and every crossing ambient mirror is
orthogonal to E and F. This proves the full arrangement identity.
The basis `(a3,s,wP)` is integral with determinant `4phi^4`, equal in
ideal to P3; it is consequently a full basis. Its vectors orthogonal
to s and r form precisely `O*wP`.

## The complete first-wall interval

For `n=a*T0+p`, the finite upper b gives `L/a>=phi+2`, with
`L=p.wP>0`. Positive Cauchy in P3 gives `a^2<=q/2`.
The other embedding gives `a'^2<=q*phi`. These apply to any actual
first wall before b and do not assume finite volume or compactness.

The trace form is `2A^2+2AB+3B^2` for `a=A+Bphi`. For q1 its trace
is below3, so only0 and +/-1 are possible, all excluded by
`0<a^2<=1/2`. For q4, the imported primitive norm4 condition gives
`n.T0 in2O`, hence `a in2O`; a/2 has exactly the excluded q1 bounds.
For q2, trace below5 and `0<a<=1` leave exactly1 and `phi^-1`.

The finite-contact statement for q2 against retained norm4 roots also
checks in the full field. A finite native label requires
`4cos^2(pi/m) in Q(sqrt5)`, allowing m=2,3,4,5,6,10. Labels3 and5
would require a K-square in the nonsquare class2; label6 requires
square6; label10 requires square `2(phi+2)` with rational norm20.
Only the right angle and the 4-bond survive. A 4-bond to either end
of the retained 5-bond gives a negative principal 4,5 chain; adding
the second nonpositive contact cannot repair positive definiteness.
Thus every remaining first root is orthogonal to both s and r.

Write `p=m*wP`, now with m in O by the full basis. At
`a=phi^-1`, `m^2=phi^2/(phi+2)` has field norm1/5, impossible even
as a K-square (and in particular impossible integrally). At a1,
`m^2=1`, and L positive chooses m1. Thus b is the unique first
simple wall; no simultaneous alternative simple wall was suppressed.

## Full compact face and limits

The Gram of `(b,a3,s,r)` has determinant `-16/phi` and all four
rank3 principal subdiagrams positive definite: H3, B3,
`I2(5)+A1` and `I2(4)+A1`. It is the compact [4,3,5] tetrahedron.
All four normals are now actual. For each triple, adjoining the six
white E+F roots gives a positive definite actual rank9 wall set.
Positive acute projection therefore gives its actual finite vertex,
feasible for every other wall. Their convex hull is the tetrahedron.
The full face is contained in that tetrahedron by its four actual
inequalities and contains that hull by convexity, so equality follows.
This is geometric full-face recognition, not recognition from a Gram
or a few vertices alone.

No error or omitted nonparabolic branch was found. The result is a
full compact three-dimensional face of the literal open R5 host.
It gives neither a compact H9 example nor a full-lattice
nonreflectivity certificate, and does not validate further cap
completion without a new global argument.

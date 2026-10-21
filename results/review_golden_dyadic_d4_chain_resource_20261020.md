# Outside-author audit: the golden dyadic D4 edge chain and sealed face

**PASS, written mathematics only.** I read every line of the following
four sources and independently rederived the norm restrictions, scalar
identities, coefficient strips, actual endpoint arguments and full
perpendicular-arrangement proof. No mathematical correction is requested.

|Source under `results/`|SHA256|
|---|---|
|`golden_dyadic_d4_extension_first_edge_20261020.md`|`3fc21952a3d119db4370ca46e39aa607666eaea7eaa54d835c3163e7a281035b`|
|`golden_dyadic_d4_extension_second_edge_20261020.md`|`a894a30bb3b17e63a434049195931e8d6bd51a160b44d85141a820a08dd2e72e`|
|`golden_dyadic_d4_extension_third_edge_20261020.md`|`11ccfde06a5ed9db05fafa803eb7a5fa8b399cf44234b8398f60fe487f350c3b`|
|`golden_dyadic_d4_fourth_vertex_and_sealed_h3_face_20261020.md`|`d6b148954b09e259d0e0598cc4f754904ba1977b373a3f491196ef87f93f9501`|

The actual-wall projection interface was separately fully reviewed at
SHA256 `3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84`.
I also read the optional opposition obstruction,
`golden_dyadic_d4_corner_opposition_obstruction_20261020.md`, SHA256
`09e6bf6fb3156a3103ef35d1512d8eb281cee23d1acd1efc206c4fe828880961`.
The four-source chain does not require that opposition map to be integral.

## Primitive pool and initial full chamber

For a primitive reflector of norm q, reflection integrality gives
`2v/q in R_n^*`. The dual exponent 2 then gives `4v/q in R_n`.
Primitivity over the golden PID forces `q|4`. Both real norms are
positive, and the totally positive units are squares, giving precisely
the necessary normalized pool `1,2,4`. At q4, all pairings with the
lattice are even; this gives the stated even time/I coordinates and
`D4` component in `2D4^*`, equivalently four coordinates congruent
modulo 2. This condition is used, rather than silently treating q1
and q4 as interchangeable squareclasses.

At zero time, expanding a golden vector as `a+phi b` gives exactly
`a^2+b^2=q` and `2(a,b)+b^2=0`. The minimum norms of D4 and
`2D4^*` exclude mixed zero roots and all nonzero b in the stated
small norms. The q4 I-coordinate vectors are imprimitive. Thus the
full zero roots are the B_n system and the two 24-vector D4/dual-D4
sets constituting F4, with the displayed simple system.

The F4-perpendicular arrangement assertion concerns the full group,
not only parabolic simple extensions. A mirror crossing that subspace
gives a finite point stabilizer on the span of F4 and the extra root.
Finite reflection-group classification leaves only F4 plus an
orthogonal A1: no irreducible rank-five finite reflection group
contains this F4 reflection subsystem. The diagonal norm1/2 roots
extend, while perpendicular q4 roots are imprimitive. This justifies
the full initial diagonal-face identity. Standard discreteness/local
finiteness is legitimate here: the conjugate form is positive definite,
so bounded physical isometries have bounded conjugate matrices and only
finitely many integral matrices.

## The first two edges

The first ray, r and its finite endpoint have the stated products and
norms. The complete finite-extension possibilities are the isolated
rows and B_(n+1), B4, F4 rows. In particular a q1 attachment to the
norm4 tail cannot replace the q4 F4 row: its required odd pairing
violates the old norm4 root's even divisibility.

I checked both embeddings of `y^2-phi x^2=S` and the complete S1/S2
strips. Their only square reconstructions fail the displayed parity
conditions. For q4, dividing x by 2 gives exactly the smaller stated
strips; the isolated winner is r and the F4 row is empty. This proves
first simplicity in the full arrangement, rather than only existence
of an integral upper mirror.

For the second ray I obtain
`Q.w=-U`, `t=U/Z` and
`Z^2-3phi U^2=3q-3epsilon-2b^2-3c^2/4` with the stated glue.
All 21 literal finite profiles are represented. The length and pairing
conditions correctly remove the apparent H4 bridges. I checked each
complete U-list, the nonsquare norms 45/99, the inert-prime-3 case,
the q4 divisibility exclusions, and the reconstruction of a4.

Independently, this second census is unnecessary for firstness: a4
is already an actual facet, and it completes the retained actual walls
to a positive-definite full-rank corner. Positive acute projection
produces the feasible finite endpoint. The segment from Q lies in the
one-dimensional retained face, so no earlier cutting wall can intervene.
This shortcut is unavailable for a merely integral, unproved upper mirror.

## Third and fourth endpoints

For the I2(5)-orthogonal third-edge rows, direct norm expansion gives
`V^2-2 beta a^2=6(q-epsilon)-3b^2`. I verified both scalar bounds,
every listed coefficient strip and the nonsquare norms 61/124. The
q4 survivor is 2c and is imprimitive. The remaining roots are c and
`v=c+e1` at the same finite parameter.

For each H3 row, completing the square gives exactly
`beta V^2-2U^2=24(2phi-3)`. Collision by the upper requires
`U^2<=28-16phi<4`, whereas its native contacts give `U>=2`.
Thus neither missing H3 continuation can enter the interval.

The distinction between a root and an actual simple wall is essential
and correctly handled: `c=v+beta1+...+beta_n`. The retained actual
B_n walls and v give the simple B_(n+1) germ; c is a positive
nonsimple short root there, not an extra isolated facet. The finite
inventory consequently identifies the stated full third endpoint.

The fourth endpoint again uses an already actual facet, a3. Its Gram
with the retained walls is the positive-definite
`B_(n+1)+H3`, so actual projection proves firstness. Direct substitution
gives `V^2=8-5phi=-phi^-5`, D4 sum -2, and the claimed whole corner.

## The nonparabolic H4 obstruction and the FULL sealed arrangement

The H4 argument is valid for arbitrary finite overgroups of the present
H3, not just simple extensions in one chosen diagram. All reflections
in H4 are conjugate. Since the conjugating reflections act integrally,
they carry a primitive norm4 H3 normal to primitive vectors of the
**same norm 4** on every H4 reflection line. This exact norm/primitivity
argument is stronger than merely observing a common K-squareclass.

For H4 simple roots the Gram is `2C`, with `det C=phi^-4` a unit.
Their even pairings imply that projection of any lattice vector has
coefficients `C^-1 k in O^4`. Their O-span would therefore split
orthogonally and integrally from R_n. Its determinant ideal is `(2)^4`,
but R_n has ideal `(2)^2`; an integral complementary lattice cannot
have the required negative dyadic determinant valuation. Thus H4 is
impossible.

Any ambient mirror crossing the H3-perpendicular hyperbolic space,
unless already orthogonal to H3, would generate such a finite rank-four
overgroup at a timelike point. Hence it is excluded. This is the full
mirror seal needed for the face statement.

Finally `phi^3 V`, c and the old I_n vectors form an integral orthogonal
basis with unimodular Gram `<-phi>+I_(n+1)`. Pairings recover integral
coefficients for every ambient vector in this K-space, proving equality
of lattices and an integral orthogonal split. Intrinsic norm1/2 roots
extend; an ambient perpendicular q4 root lies in twice this unimodular
summand and is imprimitive. A disjoint ambient mirror has constant sign
on the perpendicular space, fixed by the known feasible endpoint.
It introduces no hidden face inequality. This proves the claimed FULL
diagonal chamber identity, rather than an identification by determinant.

## Scope and remaining obligations

The result is an actual finite edge chain and a complete sealed
dimension-(n+1) diagonal face. For R_5 this face has dimension six;
it does not settle the volume, reflectivity or compactness of R_5.
The cited compactness of particular lower-dimensional diagonal hosts
is an external input, not needed to prove the arrangement identity.
No worldwide novelty or integral equivalence to another host follows.

No source was altered. This audit involved written algebra and small
hand coefficient checks only: no import, enumeration program, checker,
scientific job, build, formal proof or literature classification run.

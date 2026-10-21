# Root review: four actual edges of the golden D4 extension

20 October 2026. **PASS, complete written mathematical review.** No
scientific job, enumeration, build or Lean execution. Four actual edges
and a full-face identity do not establish compactness.

Read in full and independently checked, under `results/`:

|Source|SHA256|
|---|---|
|`golden_dyadic_d4_extension_first_edge_20261020.md`|`3fc21952a3d119db4370ca46e39aa607666eaea7eaa54d835c3163e7a281035b`|
|`golden_dyadic_d4_extension_second_edge_20261020.md`|`a894a30bb3b17e63a434049195931e8d6bd51a160b44d85141a820a08dd2e72e`|
|`golden_dyadic_d4_extension_third_edge_20261020.md`|`11ccfde06a5ed9db05fafa803eb7a5fa8b399cf44234b8398f60fe487f350c3b`|
|`golden_dyadic_d4_corner_opposition_obstruction_20261020.md`|`09e6bf6fb3156a3103ef35d1512d8eb281cee23d1acd1efc206c4fe828880961`|
|`golden_dyadic_d4_fourth_vertex_and_sealed_h3_face_20261020.md`|`d6b148954b09e259d0e0598cc4f754904ba1977b373a3f491196ef87f93f9501`|

## Norms, full initial arrangement, and firstness

For `R_n=<-phi>+I_n+D4(O_phi)`, the dual exponent 2 and primitivity
give `(q)|(4)`. Inertness of 2 and totally positive unit normalization
leave precisely the necessary classes 1,2,4. A norm4 root pairs the
lattice into 2O; this forces its time/I coordinates into 2O and its
D4 part into `2D4*`. The two-coefficient zero-height norm identity
excludes hidden mixed or golden zero roots, leaving the full `B_n+F4`
arrangement. The D4 parity condition is essential throughout.

The first-edge census includes isolated walls and every finite B/F
extension of the whole retained components. Its two-embedding strips
are complete coefficient bounds, not arbitrary height cutoffs. The
norm2 candidate with time `phi-1` fails the D4 parity condition. The
only surviving primitive first root is
`r=(2;0;0,0,0,-2phi)`, of norm4, giving `B_n+B3+A1`.

The second edge is valid. More economically, its upper is already an
ACTUAL facet a4, and its Gram with all retained actual walls is
positive definite. Positive acute projection therefore makes its
intersection feasible against every other actual wall. This proves
firstness directly. The longer saved 21-profile census independently
gives the same result and its reconstruction/parity conditions check.

For the third edge I checked all 11 profiles and rederived the scalar
equations and strips. In the I2-orthogonal rows the rational norm
obstructions 61 and 124, inertness of 2, and the norm4 time-divisibility
test leave just the stated c and `v=c+e1`. In the two H3 rows,
`beta V^2-2U^2=24(2phi-3)` and collision before `beta/2` give
`U^2<=28-16phi<4`, contradicting `U>=2`. Crucially c is a positive
nonsimple root: `c=v+sum beta_i`. The actual simple extension is v,
which enlarges B_n to B_(n+1). Treating c as an isolated new facet
would give the wrong corner.

At the fourth edge, the known actual wall a3 completes the retained
walls to `B_(n+1)+H3`. Positive acute projection proves firstness
without a new census. The stated integral vector V is orthogonal to
this whole corner and has norm `8-5phi=-phi^-5`. These calculations
also fix the sign and outgoing branch.

## Full H3 seal and exact perpendicular lattice

The H4 obstruction is integral, not a determinant comparison between
unrelated spaces. For hypothetical primitive norm4 H4 simples, the
Gram is `2C` with `det C=5-3phi` a unit. Every lattice pairing is
even, so projection has coefficients `C^-1 k` in O. This splits the
H4 span integrally. Its determinant ideal `(2)^4` cannot be an
orthogonal summand of the integral R_n, whose determinant ideal is
`(2)^2`. An integral complement would have a fractional determinant
ideal. The contradiction is sound.

A mirror crossing the H3-perpendicular hyperbolic space generates,
with H3, a finite isotropy group of rank at most4. A nonorthogonal
crossing would require H4 by finite reflection classification; the
H3 cannot lie in the other rank4 types. In H4 all reflections are
conjugate, so primitive normals can be normalized to norm4 as in the
given H3. The obstruction therefore applies even if the given H3
was not chosen as a parabolic subsystem. This establishes a seal for
the entire integral arrangement.

The vectors `V0=phi^3 V`, c and the original n coordinate units are
integral, orthogonal and have unit norms `-phi,1,...,1`. Pairing with
this unimodular basis recovers integral coefficients for every
ambient integral vector in its K-span. It consequently gives the
EXACT perpendicular lattice and an integral orthogonal splitting,
not just a finite-index diagonal sublattice. All its intrinsic norm1/2
reflectors extend. An orthogonal norm4 reflector would lie in twice
this unimodular summand and fail primitivity. Together with the seal,
this proves the full diagonal face identity.

## Opposition and scope

The complete simultaneous corner reversal, including both component
swaps, sends e0 to `e0+Delta/sqrt(5)`. Its D4 coordinate
`1/sqrt(5)` is not integral. The displayed full columns follow from
the two commuting reflection formulas; this particular period
construction is correctly rejected.

For n=5 the sealed face has diagonal dimension6. Its known compactness
does not prove compactness of R_5. No infinite period, complete global
wall list, new finite-volume example, commensurability-class exclusion
or dimension ceiling has been obtained. The initial F4 face and the
H3 face are proved with full arrangement identities; the single r
facet cannot be identified with a standard diagonal chamber just by
its metric. Keep this research as one construction packet, with no
additional primary TeX or Lean expansion yet.

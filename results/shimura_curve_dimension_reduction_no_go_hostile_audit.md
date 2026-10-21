# Hostile audit: higher-dimensional genus-zero curve reduction

**Audited owner bytes:**
`ae08232dc4c5ec9963f950ac63c160e34918acfcba913caee6ebf29e4ac704f5`

**Verdict:** **PASS, with the advertised scoped no-go.**  The report proves a
conditional arithmetic Coxeter-face descent and gives two valid controls
against using genus-zero Fuchsian curves alone as an all-dimensional return
mechanism.  It does not rule out richer decorated curves or other AG
carriers.

## Conditional positive theorem

Bogachev--Kolpakov, Theorem 1.4, says exactly that a Coxeter face of
dimension at least two in a quasi-arithmetic Coxeter polytope is
quasi-arithmetic with the same ground field.  The report separately assumes
that the marked two-face is arithmetic, so it does not upgrade
quasi-arithmeticity silently.  Doubling a finite-area Coxeter polygon gives a
coarse sphere with the stated stack orders and cusps.  A finite map from that
coarse \(\mathbf P^1\) to a maximal arithmetic Fuchsian quotient forces the
target coarse genus to be zero.  Maximal arithmetic Fuchsian groups are
congruence, so the Long--Maclachlan--Reid theorem applies.  The last
finiteness clause correctly retains reflection maximality of the face; a
finite list of overgroups alone would not bound their finite-index reflection
subgroups.

## The (H_6) control

1. The norm-two matrix is Felikson's compact tetrahedron (H_6): a
   four-cycle with one edge labelled (5) and three labelled (3).  The
   selected determinant signs give signature ((3,1)); the conjugate
   Sylvester chain is positive; all four rank-three principal packets are
   positive definite.  The integral cycle field is
   \(\mathbf Q(\sqrt5)\), so Vinberg arithmeticity is correctly applied.

2. Felikson's complete bounded-tetrahedron decomposition classification
   lists no nontrivial decomposition of (H_6).  A proper reflection
   overgroup would tile (H_6) by its smaller chamber and hence give such a
   decomposition.  This is a legitimate literature input for absolute
   reflection maximality, not a reflection-finiteness theorem.

3. The intrinsic facet calculation uses the correct Schur complement, not
   node deletion.  Independent evaluation confirms the four triples of
   angle orders lie in the displayed consecutive-integer intervals.  Every
   intrinsic facet therefore has non-Coxeter angles.  All proper principal
   diagrams and all vertex links are spherical, while the full rank-four
   packet is minimal indefinite.

4. The mirror-centralizer topology survives hostile inspection.  The odd
   diagram is (C_4), so all reflections are conjugate and
   \(\Gamma_\Omega\cong\pi_1(C_4)\cong\mathbf Z\).  More is needed than this
   Betti-number statement, and the report supplies it: Allcock's space (X)
   is the tetrahedral boundary (S^2) with the two closed, disjoint,
   opposite even ridges removed.  It is an open annulus.  The lifted tiles
   form an infinite strip; adding the two mirror-wall ends and quotienting by
   the deck translation gives a bordered mirrored annulus.  Odd-edge
   transports are products of two ambient reflections and preserve the
   selected normal half-line, so they preserve orientation on the mirror.
   The effective orientation subgroup therefore gives the double across the
   two mirror boundary circles, whose underlying closed coarse surface has
   
   \[
   g=2\cdot0+2-1=1.
   \]
   
   Compactness is also consistent: the mirror is (k)-rational and its
   ternary restriction is anisotropic because the ambient arithmetic form is,
   so the arithmetic Fuchsian stabilizer is cocompact.  This is a torus with
   possible orbifold points, not a cylinder whose cusp compactification is a
   sphere.  A larger arithmetic overgroup might have genus zero, as the report
   explicitly allows; genus zero is not forced by the centralizer.

## Same curve, different ambient dimension

The Schur complement of the first two walls of ([5,3,3,3]) normalizes to
the triangle ((2,3,10)).  Its nontrivial cycle
\((2\cos(\pi/10))^2=(5+\sqrt5)/2\) is integral and its nonidentity conjugate
Gram is positive definite, so it is arithmetic over \(\mathbf Q(\sqrt5)\).
The Makarov straight (5)-prism based on ([5,3,3,3,3]) has an isometric
([5,3,3,3]) base, and hence literally contains the same marked triangle.

The cited minimum-covolume inputs do imply absolute reflection maximality.
The rotation subgroups of the (4)-simplex and the Makarov prism are the
unique compact orientable arithmetic minima in dimensions four and five.
Any proper discrete reflection overgroup has finite index, is commensurable
with the arithmetic group and therefore arithmetic, and has a proper
orientation-preserving overgroup of smaller covolume, a contradiction.

## Replay boundary

The owner script and test have hashes
`f421ebb84c20df45bd595675f8b86c4ef27a5c9ac7c5bdda5f88f649f7395595`
and
`c937ea0d32973cf290de4dad9a6830a8c4151a378b8f3669d1a8d82c7a237540`;
all tests pass.  The code exactly checks Gram, inertia, Schur complements,
and angle intervals.  Its variables declaring an annulus and two boundary
components are a regression encoding, not an independent proof of that
topology; the proof is the Allcock/tetrahedral-boundary argument above.
Likewise absolute maximality is, correctly, a cited classification/minimum
input rather than a symbolic output.

No substantive defect remains, and the conclusion is no stronger than the
controls justify: genus-zero arithmetic curves alone neither arise
automatically nor bound transverse rank and Gram gluing.

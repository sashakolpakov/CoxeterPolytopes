# Independent hostile audit: the bi-anticanonical section-algebra gate

**Verdict (2026-08-11): PASS after making the exact-return hypotheses
literal.**  The two ten-point blowups really have the same marked lattice,
proper-point proximity data, labelled-exceptional contraction fan, and formal
jet descriptors but opposite physical answers for \(-2K\).  Proposition 4.1
is an exact intrinsic criterion, and the repaired Corollary 4.2 invokes
reflection maximality only after every hypothesis of the strong-Coble return
theorem has been supplied.

## 1. Same-frame and Bézout audit

For ten proper distinct points of \(\mathbf P^2\),

\[
 K=-3H+\sum E_i,\qquad F=-K,\qquad B=-2K
\]

give \(K^2=F^2=-1\), \(B^2=-4\), \(K.B=2\), and \(p_a(B)=0\).
The exceptional curves are ten disjoint labelled \((-1)\)-curves.  Hence
every one of their \(2^{10}\) labelled subsets contracts, and every one of
the \(10!\) labelled complete orders occurs.  This is correctly called the
**labelled-exceptional** fan; it is not claimed to enumerate every Mori
contraction or every negative curve.

On the smooth-cubic branch, a second cubic through the ten marked points
would contribute at least ten to an intersection of total degree nine.
A sextic double at all ten points would contribute at least twenty to its
intersection with the cubic, whose total is eighteen.  It must therefore be
the cubic plus another cubic through the ten points, hence the square of the
unique cubic.  Thus

\[
                         h^0(-K)=h^0(-2K)=1,
                         \qquad |-2K|=\{2F\}.
\]

On the ten-node irreducible-sextic branch, a cubic through every node would
give local intersection at least \(10\cdot2=20>18\).  A second sextic
double at every node would give at least \(10\cdot4=40>36\).  Bézout and
irreducibility therefore give

\[
                         h^0(-K)=0,\qquad h^0(-2K)=1.
\]

The frozen Winger packet supplies an exact ten-node example, so the positive
branch does not rest on generic-projection language.  Blowing up an ordinary
node separates its two transverse branches; the strict transform is
integral and smooth.  Its arithmetic genus is zero, hence it is rational.

## 2. Proposition 4.1 audit

Let \(K_S^2=-1\) and suppose a reduced irreducible
\(B\in|-2K_S|\) exists.  Any \(D\in|-2K_S|\) has

\[
                              D.B=B^2=-4.
\]

Negative intersection forces the irreducible curve \(B\) to be a component
of \(D\).  The residual divisor is effective and linearly trivial; on a
projective integral surface it is zero.  Thus \(|-2K|=\{B\}\).  If
\(A\in|-K|\), uniqueness would give \(2A=B\), contradicting reducedness.

Conversely, a geometrically integral generator of the one-dimensional
\(-2K\) piece is a reduced irreducible curve.  Adjunction gives arithmetic
genus zero.  For an integral curve over \(\mathbf C\),

\[
                  p_a=g(\widetilde B)+\sum_p\delta_p,
\]

so both the normalization genus and every singularity defect vanish.  The
curve is smooth rational.  The report uses “divisor-prime” precisely for
this geometric-integrality condition; it does not assume finite generation
of the Cox ring.

## 3. Repaired Corollary 4.2 audit

The corollary now explicitly retains the required Coxeter-side data:
finite-sided finite volume, distinct irredundant negative facets defining
\(P\), adapted-host reflection integrality, and no proper discrete
reflection overgroup.  It also assumes a nonempty pairwise-orthogonal
\((-4)\) boundary packet with \(\sum B_i=-2K\), together with the declared
\((-1)\) and \((-2)\) adjunction rows.

Prime sections for all these classes on one and the same marked surface give
integral curves.  The adjunction rows make them smooth rational.  Distinct
effective boundary curves have nonnegative local intersections, while their
numerical intersections are zero, so they are disjoint.  Their reduced sum
lies in \(|-2K|\).  For any other \(D\in|-2K|\),

\[
                        D.B_i=(-2K).B_i=-4
\]

forces every \(B_i\) to be a component.  The residual effective linearly
trivial divisor is zero.  This componentwise uniqueness shows that the
reduced sum is the only member and it cannot be
twice an anticanonical divisor.  This constructs the strong-Coble branch,
after which the cited maximal-return theorem legitimately supplies the full
nef cone and excludes hidden negative walls.

Necessity is immediate because an exact physical facet is an effective
integral curve.  Thus the repaired “if and only if” is carrier-relative and
sound; it neither constructs the marked surface nor proves nonemptiness of
the simultaneous prime-section locus.

## 4. Exact boundary

The theorem does not claim the two calibration surfaces have the same nef
cone or raw Coxeter chamber.  It does not replace coordinate-valued
evaluation matrices by formal descriptors, prove finite Cox-ring generation,
or establish global Property 4.  Its no-go is exactly that lattice,
proximity, labelled-exceptional contractions, degrees, multiplicities, and
formal matrix shapes do not determine divisor primality.

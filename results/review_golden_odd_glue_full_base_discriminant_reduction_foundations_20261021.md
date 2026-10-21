# Full outside-author review: good-prime odd golden glue

21 October2026. **PASS at the explicit written-input scope; no correction
requested.** Reviewed in full:
`results/golden_odd_glue_full_base_and_discriminant_reduction_20261021.md`,
SHA256 `b107f0de4f76da942100e0c82867cd6d4eeee40d833304fee67cebd84d3e41ce`
(248 lines). This is mathematical source review, not a mathematical
program, worker import, enumeration, test, compiler, Lean or TeX run.
No candidate source, paper, guide or accepted bound was changed.

I did not author this source. I authored earlier literal norm79-module
and edge-contract work and contributed to inherited actual-foot/arrangement
foundations. Those shared contributions are disclosed; the new original-
coordinate extension and every implication below were independently
reconstructed. The named analytic and calibration inputs remain inputs,
not newly executed or externally refereed results in this review.

## Exact input boundary

All six direct bindings match:

| Input | SHA256 |
|---|---|
| results/bugaenko_facet_binary_extension_and_collar_audit_20261020.md | `774e87f7abb33680a06879a17a5be9f06e2218bf12a95023bcff5c7046ba1227` |
| results/golden_e8_binary_family_full_base_facet_20261020.md | `ea5d96e331ac91cb9efc457db2ddaa24bac3341e9e6781038150765895675afd` |
| results/golden_e8_binary_family_discriminant_collar_and_finite_cutoff_20261020.md | `841a4b7f9e75cdd49b6384920439d998cc401f00bbf75aed92444892c851d425` |
| results/bugaenko_facet_full_reflection_volume_bridge_addendum_20261020.md | `e61a5531051aac247382b9e1ec2f4536bb7a9b89fba3352febf7fe1bf01cf6a3` |
| results/review_bugaenko_full_reflection_volume_bridge_batch_guard_20261020.md | `56be6455426bb30d66076ba7e197aedc176c622d7c7f735f56988a09e2968761` |
| data/known_examples/bugaenko_compact_h8.json | `c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c` |

I read the full relevant sources, including the general rank-ten volume
bridge and its stated primary-source analytic hypotheses. The complete
Bugaenko calibration, golden PID/unit facts, finite intersection-angle
classification and analytic inequalities are explicit inherited inputs.
The new theorem concerns the LITERAL even module E8(O)+B_c with odd
cross pairing1 and good-prime determinant. It is not a transport theorem
for arbitrary invariant lattices or an entire K-form class.

## Admissibility and every primitive norm class (lines11–50)

For B_c=[[-2phi,1],[1,2c]], det B_c=-D. The hypotheses D>0,D'<0
give physical signature(1,1) and conjugate positive definiteness because
-2phi'>0. Thus E8+B_c is admissible of signature(9,1). Also b=phi*D
is totally positive. Its nonunit proof checks: b=phi mod4, whereas
positive units phi^(2j) have the three residues1,1+phi,2+3phi, using
phi^6=1 mod4. None equals phi.

Evenness gives primitive reflector norm2d; reflection integrality and
a primitive-coordinate Bezout combination give Gm/d integral. The
adjugate then gives d|det G=-D, equivalently d|b. Totally positive
units being squares permits rescaling ideal-associate full-b norms
to exactly2b. These arguments keep all divisor exponents; there is no
prime/squarefree determinant assumption or hidden passage to an
index-two coordinate frame.

## Full-modulus isolation and actual full base (lines52–113)

I expanded in the ORIGINAL X,Y basis. The root r=X+2phi*Y is primitive,
r²=2b and Gr=(0,D), so Gr/b=(0,phi^(-1)) is integral. Its orthogonal
lattice is exactly y0, namely OX+E8 with X²=-2phi.

For norm2b, reflection integrality forces m_E=bZ and
y=2phi*x+b*nu. Division of its exact norm by2b gives

    1=x²+D*x*nu+c*b*nu²+(b/2)Z².

Since D=b/phi and E8 is even, x²=1 mod the ENTIRE ideal(b). For two
such roots the cross-product divided by b is2x1x2 mod b. Its square
therefore equals4 mod b, including repeated prime factors. The complete
finite arrangement-angle list, including non-native arrangement angles,
has A-4 a unit at every prime over(b), since those primes avoid30.
This excludes all finite intersections, including A0. Parallel distinct
mirrors would yield a nonzero K-null vector, excluded by conjugate
definiteness. Distinct norm2b mirrors are consequently ultraparallel.

For a proper divisor norm, (r,m)=b*eta and its nonzero finite invariant
is (b/d)*eta². A prime of b/d contradicts the same finite-angle list.
Thus every mirror meeting r is orthogonal. Orthogonality forces y0;
the product with Y then gives d|x, and E8 unimodularity gives d|m_E.
Primitivity makes d a unit, so these are exactly the norm2 mirrors of
M0. Conversely every such M0 norm2 root reflects L_c integrally.

This identifies the FULL base arrangement with the calibrated Bugaenko
arrangement, rather than merely constructing an acute finite seed.
The finite norm alphabet and two-embedding coordinate bounds give local
finiteness. A generic interior base point has an adjacent chamber on
each side; no nonorthogonal mirror cuts any part of the base. The whole
compact F is therefore an actual facet, with full E8+r vertex germ.

The reversal extension is also on the original module. S0X-X in2M0
makes SY=Y-(S0X-X)/(2phi) integral. It fixes r, is isometric and
involutive, and preserves the inward side and actual full facet.
The full adjacent chamber is thus preserved. This proves one genuine
chamber involution, not an automatic infinite-order period.

## All-divisor collar and the conditional volume cutoff (lines115–163)

For general norm2d write e=b/d and y=2phi*x+d*nu, m_E=dZ. Direct
norm expansion gives the printed equation(4). With gamma=y/phi,
A=e*gamma², subtracting four times that equation gives

    A-4=d*(phi^(-1)*nu²-2Z²),

because b/phi²-4c=phi^(-1). This cancellation is valid for arbitrary
c, not just the old unit-coupled family. Gamma is integral since phi
is a unit. Conjugate positive definiteness and nonproportionality give
0<A'<4 strictly; a nonzero algebraic product cannot have zero conjugate.
The two nonzero integral norms consequently give A>N(e)/4 and
A-4>N(d)/4. Therefore

    256*sinh(ell)²*cosh(ell)²>N(b),
    ell>(1/2)asinh(sqrt(N(b))/8).

This is the distance from the entire base hyperplane for EVERY divisor
norm, including mixed roots. The inward normal flow remains inside all
walls: orthogonal ones retain sign and all others lie beyond that global
distance. Its injectivity and cosh(t)^8 Jacobian give the stated collar
volume. No finite selected-wall prefix is used for containment.

The E8 ball-sector calibration supplies Vol(F)>2^(-38). The elementary
integration is correctly normalized: exp(8rho)>N²/256 and
integral cosh(t)^8 dt>(N²-256)/2^19. Multiplication gives the lower
bound (N²-256)/2^57 for N>=16; for smaller N its negative right side
is handled separately by positive volume. There is no invalid
multiplication of a negative lower bound by Vol(F).

The analytic upper used is the pinned GENERAL full-integral-reflection-
chamber theorem, valid for any admissible rank-ten O-lattice when its
full reflection chamber has finite volume. Its proof retains normality,
the finite chamber stabilizer and its congruence ambient group. It does
not assume the reflection group itself is congruence or maximal.
The bounds U<3^3851 and2^57<3^38 yield
N²<2^57U+256<3^3890, hence N<3^1945. This part is conditional on
that named analytic interface; this review does not replace its primary
analytic proofs. No assertion about arbitrary reflection subgroups is
made.

## Finite marked reduction and the norm239 start (lines165–228)

I checked the proposed integral basis change directly. Y_new=X+phi³Y
has (X,Y_new)=1, its half-square is phi^6*c+phi², and D_new=phi^6*D.
Its determinant is the unit phi³, so the inverse is integral. Iterating
in the correct direction puts1<D<phi^6 and terminates. Equality at1
would contradict D'<0. Then0<c<phi² and
c'=phi*(1+N/D)/4<phi*(1+3^1945)/4. The formulas for integer coefficients
of c=A+Bphi give a finite rectangle. This is a complete necessary
normalized-host algorithm with an enormous bound, not an executed or
practical reflective-host classification. Bad primes/even glue remain
outside its scope.

For c4 the norm239, embedding signs and coprimality check. The value
lies strictly between the old unit-coupled determinant norms131 and319,
so it is not an isometric re-presentation of any literal member of that
strictly increasing sequence. The odd prime239 valuation separates its
K-squareclass from the specifically listed closed golden modules.
Neither fact establishes separation from every later K-similarity class
or a literature-wide novelty claim.

The existing uniform mixed-E8 collar774 Section5 is GLOBAL on H(B),
so no mixed mirror can cut the binary edge axis; this is stronger than
a finite-width claim around one base point. For prime norm239 the
primitive alphabet is exactly2 and2b up to square units. I rederived
both original-module rows from

    2phi*m²=D*y²-(2phi*x-y)².

For norm2, put y=-v and 2phi*x-y=u>0. Then x=(u-v)/(2phi),
u²=Dv²-4phi, parity u-v in2O and collision u/(Dv) follow. For norm2b,
put y=-phi*z and 2phi*x-y=b*w; then x=(Dw-z)/2,
z²-Dw²=4 and collision w/z follow. Since D=1 mod2, both coordinate
integrality and reflection-covector integrality reduce exactly to
z-w in2O. The first covector is -u (=-b*w in the second row).
The conjugate norm equations give precisely the stated necessary
v' and w' strips.

After a genuine finite upper t*<1/sqrt(D), these collisions bound the
physical controls while the conjugate strips bound their other
embeddings, giving a complete finite firstness box. Until an upper is
provided the source rightly asserts no solution, minimum, firstness,
next actual corner or completed finite-volume chamber. Neither a formal
Gram matrix nor a bounded absence is promoted to an example.

## Verdict

All new original-coordinate, full-base, collar, reduction and norm239
profile steps check at their declared inputs. The source correctly
distinguishes the extension of an already known ineffective finite
reduction from a first finiteness theorem. It supplies no new H9
polytope, general odd-glue reflectivity decision, closed norm239 host,
whole form-class exclusion, ambient ceiling or runtime/Lean credit.

# Root review of the growing base collar and finite family cutoff

20 October 2026. Complete written review **PASS**, relative to the
explicit external analytic inputs in the reviewed volume bridge.
Source `golden_e8_binary_family_discriminant_collar_and_finite_cutoff_20261020.md`,
SHA256 `841a4b7f9e75cdd49b6384920439d998cc401f00bbf75aed92444892c851d425`.
No scientific execution, factorization, search, Lean or TeX was used.

I independently reconstructed the whole argument, rather than relying
on the finite cutoff as a numerical report. For every primitive norm2d
root, the unit binary basis, E8 unimodularity and coprimality with2
give m=lambda*r+d*nu*c+d*z and the displayed equation(4). It is valid
for all divisor exponents. With e=b/d and gamma=2lambda-d*nu, direct
expansion gives

    e*gamma^2-4=d*((b-4)*nu^2-2*z^2).

Full-base isolation makes any nonorthogonal distinct mirror
ultraparallel to r. Thus A=e*gamma^2=4cosh(ell)^2>4, while positive
definiteness at the conjugate embedding gives 0<A'<4. Equality in
conjugate Cauchy would imply K-proportional normals and the same
mirror, which has been excluded. The two nonzero integral norms
therefore give A>N(e)/4 and A-4>N(d)/4. Their product is exactly
the strict inequality 256*sinh(ell)^2*cosh(ell)^2>N(b).
Consequently ell>(1/2)asinh(sqrt(N(b))/8). This argument includes
mixed E8 coordinates and all proper divisors, without squarefreeness.

I checked the geometric containment for the full arrangement. An
orthogonal mirror keeps its inequality under the normal flow, apart
from a positive cosh factor. A nonorthogonal mirror misses the entire
base hyperplane and has the same oriented sign on its future sheet;
its global distance exceeds the stated collar width, so it cannot
cut any of these normal segments. The normal exponential map is
injective, with volume factor cosh(rho)^8. This proves a full-base
volume bound, not merely one on the binary axis or a finite prefix.

At the E8 vertex, the three remaining calibrated walls have time
coordinates 1,phi,phi^4 and norm2. The radius1/2 E8 sector therefore
lies in the base. The hyperbolic ball has volume strictly larger
than its Euclidean counterpart, and |W(E8)|<2^30 gives Vol(F)>2^-38.
The exponential integration bounds are correctly normalized:

    integral_0^D cosh(rho)^8 d rho
       >(exp(8D)-1)/2048 >(N(b)^2-256)/2^19.

Thus Vol(C)>(N(b)^2-256)/2^57. For the small negative right-hand
side the source correctly uses positive volume, rather than
multiplying a negative lower bound by an independent positive one.

The full-chamber upper is used with its actual scope. If W is cofinite,
the admissible anisotropic arithmetic input makes C compact. FULL W
is normal in the full integral congruence group Gamma. Its chamber
stabilizer A is finite, and the torsion-free level3 kernel injects A
into GL10(F9), giving |A|<3^200. The bridge retains the m^19 factor
from the quotient conformal volume and the passage back to C. It
does not assume W itself is congruence or A is trivial.

I also checked the original author-hosted Agol--Belolipetsky--Storm--Whyte
PDF on 20 October2026: Theorems3.1 and4.2 give the quoted spherical
quotient growth and spectral/conformal inequality; Theorem5.9 gives
lambda1>=(2n-3)/4 for the congruence quadratic-form group, hence15/4
in dimension9. The explicit hyperplane-count constant and application
to Gamma are in the pinned bridge, whose normal-group argument is
consistent with Section6 of that paper. Reference:
https://homepages.math.uic.edu/~kwhyte/ahfinite.pdf

Finally I checked every coarse integer estimate in the source. They
give U<3^3851, then N(b)^2<3^3890 and N(b)<3^1945. For all k>=-1,
N(b)=15+4(P-P^-1)>4P, and27<phi^7. Hence
2k+1<13615/3 and k<6806/3, which for integral k is k<=2268.
This bound is deliberately coarse. It excludes every literal member
k>=2269, leaving the earlier finite range open apart from separately
closed members. It is neither a practical enumeration authorization
nor a field-wide, commensurability-class or dimension exclusion.

The growing collar was discussed with the author before write-up;
the complete estimates and volume comparison were reconstructed here.
The external analytic and arithmetic inputs remain explicit premises,
not new internal lemmas presented as published facts.

# Outside review: uniform golden collar and finite cutoff

20 October 2026. Independent mathematical/source reading only. No
scientific execution, root census, factorization, Lean or TeX job.

**PASS in its stated literal-family scope.** Full source reviewed:
`golden_e8_binary_family_discriminant_collar_and_finite_cutoff_20261020.md`,
SHA256 `841a4b7f9e75cdd49b6384920439d998cc401f00bbf75aed92444892c851d425`.
It gives the full-base normal collar for every k>=-1 and, conditional
on the explicitly imported analytic upper, excludes finite covolume
for k>=2269. It does not close the remaining finite interval, classify
all marked binary extensions, or improve a dimension bound.

I independently rederived the divisor-normal calculation. In the unit
binary basis(r,c), every primitive normal has the necessary form
`m=lambda*r+d*nu*c+d*z`, including repeated prime exponents in d|b.
The norm equation gives, for gamma=2lambda-d*nu,

    A=(b/d)*gamma^2=4cosh(ell)^2,
    A-4=d*(P*nu^2-2z^2).

This identity retains the full mixed E8 contribution; setting z=0 is
not used. For a distinct nonorthogonal mirror, the previously proved
full-base isolation gives ultraparallelism. The conjugate positive
definite Gram gives STRICT 0<A'<4. Consequently integer norms give
`A>N(b/d)/4` and `A-4>N(d)/4`, independently of squarefreeness and
factorization. Their product is exactly
`256*sinh(ell)^2*cosh(ell)^2>N(b)`. This yields the claimed strict
distance lower bound `ell>(1/2)asinh(sqrt(N(b))/8)`.

The Fermi containment is genuine full-chamber geometry. An ultraparallel
mirror has constant strict sign on the future base hyperplane; its
sign is fixed by the actual base facet. It cannot cross any normal
segment shorter than its global distance from that hyperplane.
Orthogonal mirrors retain the original inequality under normal flow,
and the base wall itself is satisfied on the chosen inward side.
Thus the whole one-sided collar of the full compact base lies in C,
even without assuming C has finite volume. Normal coordinates are
injective and have the dimension-nine volume element cosh(rho)^8.

I checked the calibrated base fixture, SHA256
`c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c`.
Its three non-E8 walls have time coordinates phi,1,phi^4 and norm2.
Their distances from the E8 vertex are at least asinh(sqrt(phi))>1/2.
The radius1/2 E8 Weyl sector therefore supplies the stated weak bound
Vol(F)>2^-38. Integrating cosh(rho)^8>exp(8rho)/256 gives
Vol(C)>(N(b)^2-256)/2^57. The source correctly treats negative right
sides separately instead of multiplying a negative bound by Vol(F)'s
lower bound.

The analytic upper remains the pinned imported theorem, not a new
analytic result of this review. I reread the complete bridge
`bugaenko_facet_full_reflection_volume_bridge_addendum_20261020.md`,
SHA256 `e61a5531051aac247382b9e1ec2f4536bb7a9b89fba3352febf7fe1bf01cf6a3`,
and my prior primary-source review, SHA256
`56be6455426bb30d66076ba7e197aedc176c622d7c7f735f56988a09e2968761`.
Its full integral congruence group and normal full reflection subgroup
apply to these free admissible rank-ten lattices. The finite chamber
stabilizer factor is retained to power19, giving the displayed U with
3^3800; neither trivial stabilizer nor congruence of W is assumed.

Finally, all cutoff conversions check with strict inequalities:
U<3^3851, 2^57<3^38, N(b)^2<3^3890, and N(b)<3^1945.
The exact expression N(b)=15+4(P-P^-1)>4P holds throughout k>=-1.
Combining 27<phi^7 with P=phi^(2k+1) gives k<6806/3, hence
k<=2268 for every finite-covolume member. No factorization, finite
search, first-edge premise, or claim about the other members is hidden
in this conclusion. No source correction is requested.

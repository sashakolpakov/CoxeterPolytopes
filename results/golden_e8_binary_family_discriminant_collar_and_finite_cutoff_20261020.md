# A growing full-base collar and a finite golden-family cutoff

20 October 2026. Paper proof, awaiting source-bound independent review.
No scientific execution, factorization, root census, Lean or TeX job.
The accepted full-base and first-edge sources and all frozen files are
unchanged. The collar was independently reconstructed by root before
this write-up; the deliberately coarse numerical cutoff below uses the
previously reviewed full-reflection volume bridge.

## Statement and inputs

Put O=Z[phi], phi^2=phi+1, and, for k>=-1,

    P=phi^(2k+1), b=4+P,
    L_k=E8(O) orthogonal-sum [[-2phi,phi^(k+1)],
                             [phi^(k+1),2]].

Let C be its FULL integral reflection chamber having the actual
Bugaenko base facet F with outward normal r=(phi^k,2;0), r^2=2b.
Write N=N(b)>0. Then C contains the entire one-sided normal collar
of F of width

    D=(1/2) asinh(sqrt(N)/8).                           (1)

Consequently, whether or not C has finite volume,

    Vol(C) >= Vol(F) integral_0^D cosh(rho)^8 d rho
           > (N^2-256)/2^57.                           (2)

If C has finite covolume, the reviewed fixed-field spectral/conformal
volume bound implies

    k<=2268.                                          (3)

Thus every k>=2269 is nonreflective. This is a literal-family reduction,
not exclusion of every remaining member, a commensurability-class
classification, an enumeration proposal, or a new dimension bound.
It does not require the separately proved first binary edge.

The algebraic and full-facet premises are in
`golden_e8_binary_family_full_base_facet_20261020.md`, SHA256
`ea5d96e331ac91cb9efc457db2ddaa24bac3341e9e6781038150765895675afd`.
The uniform FULL-chamber volume upper, including the chamber-stabilizer
factor and its external analytic inputs, is in
`bugaenko_facet_full_reflection_volume_bridge_addendum_20261020.md`,
SHA256 `e61a5531051aac247382b9e1ec2f4536bb7a9b89fba3352febf7fe1bf01cf6a3`.
Its outside review is
`review_bugaenko_full_reflection_volume_bridge_batch_guard_20261020.md`,
SHA256 `56be6455426bb30d66076ba7e197aedc176c622d7c7f735f56988a09e2968761`.

## 1. Every nonorthogonal mirror moves away from the whole base

Set c=(0,-1;0). The vectors r,c form an O-basis of the binary summand,
with Gram [[2b,-b],[-b,2]]. A primitive mirror normal has square2d,
where d is totally positive, d|b, and d is coprime to2. Write e=b/d.
Unimodularity of E8 and reflection integrality force its E8 part to be
d*z, z in E8(O). Integrality against c forces its c coefficient to be
d*nu. Thus every primitive mirror has the form

    m=lambda*r+d*nu*c+d*z,
    e*lambda^2-b*lambda*nu+d*nu^2+d*z^2/2=1.           (4)

These are necessary equations for ALL divisor norm classes, with their
full prime exponents. No squarefreeness is assumed.

Suppose m is not proportional or orthogonal to r. Full-base isolation
says that their hyperplanes are ultraparallel. With gamma=2lambda-d*nu
and ell their global hyperplane distance, put

    A=e*gamma^2=4cosh(ell)^2.

The conjugate positive-definite Gram gives 0<A'<4. Equation(4) gives
the useful second divisibility, valid also for mixed mirrors:

    A-4=d*xi,      xi=P*nu^2-2z^2 in O.               (5)

Since gamma is a nonzero algebraic integer,
N(A)=N(e)*N(gamma)^2>=N(e). Since xi>0 and xi'<0,
|N(A-4)|=N(d)*|N(xi)|>=N(d). The strict conjugate bounds therefore give

    A>N(e)/4,          A-4>N(d)/4.

Writing t=sinh(ell)^2 and multiplying yields

    256*t*(t+1)>N(d)*N(e)=N.

Equivalently 64sinh(2ell)^2>N, so ell>D in(1).
This is a distance bound from the ENTIRE base hyperplane, independent
of d, nu, z and any chosen finite prefix.

## 2. The full normal collar is actually inside C

Let R=r/sqrt(2b) and take x in F, x^2=-1, future. The inward normal
flow is

    x_rho=cosh(rho)*x-sinh(rho)*R,    0<=rho<=D.

Every mirror orthogonal to r has at x_rho the same inequality as at x,
multiplied by cosh(rho). Every other mirror misses the whole base and
is negative on its future hyperplane, since it is negative on F.
Its global distance from that hyperplane is greater than D by Section1,
so it cannot meet this normal segment. The oriented base inequality
itself is inward. Hence all x_rho lie in C. This argument uses the
full arrangement, and does not infer containment from a bounded seed.

The normal exponential map is injective and its H9 volume element is
cosh(rho)^8 d rho dVol_F, proving the first inequality in(2).

For an explicit, very weak lower bound on Vol(F), use its E8 vertex.
The three non-E8 walls have norm2 and time coordinates 1,phi,phi^4 in
the <-2phi>+E8 calibration. Their distances from that vertex are at
least asinh(sqrt(phi))>1/2. Thus F contains the radius1/2 hyperbolic
ball sector of one E8 Weyl chamber. Since |W(E8)|=696729600<2^30,

    Vol(F) > [pi^4/(24*2^8)]/|W(E8)| >2^-38.          (6)

The strict ball comparison uses sinh(rho)>rho; pi^4/24>1 suffices.
The fixture is `data/known_examples/bugaenko_compact_h8.json`, SHA256
`c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c`.

Finally exp(2D)>2sinh(2D)=sqrt(N)/4, so exp(8D)>N^2/256.
The elementary bound cosh(rho)^8>exp(8rho)/2^8 gives

    integral_0^D cosh(rho)^8 d rho
       > [exp(8D)-1]/2048 > (N^2-256)/2^19.

Together with(6) this proves(2). If its right side is negative the
claim is automatic from positive volume; multiplying the negative
bound by(6) is unnecessary. For N>=16 the displayed bounds apply
directly. This explicitly covers the small-N endpoint as well.

## 3. The genuine full-chamber upper gives k<=2268

Assume the full reflection group W has finite covolume. The admissible
nonphysical definite embedding makes the lattice anisotropic, so C is
compact. For Gamma=O^+(L_k), normality of FULL W gives
Gamma=W semidirect-product A, with A the finite chamber stabilizer.
The pinned volume bridge proves |A|<3^200, using reduction modulo3,
and applies the published congruence spectral/conformal inequalities to
Gamma. Its actual bound is

    Vol(C)<U,
    U=(12/5)^(9/2)*20*50^9*3^3800*pi^5/12.            (7)

The power3800 includes m^19 with m=|A|; no assumption A=1 or that W
itself is congruence is made. The analytic input is the precise one
already isolated in the bridge: Agol--Belolipetsky--Storm--Whyte,
Theorems3.1,4.2,5.9 and their normal-reflection application in Section6.

For a simple integer cutoff use

    (12/5)^(9/2)<3^5, 20<3^3, 50^9<3^36,
    pi^5/12<4^5<3^7,
    U<3^3851,         2^57<3^38.

Combining(2),(7) gives

    N^2<2^57 U+256<3^3889+256<3^3890,
    N<3^1945.                                        (8)

Also N=15+4(P-P^-1)>4P for every k>=-1. Since
27<phi^7=13phi+8, equation(8) implies

    phi^(2k+1)=P<3^1945<phi^(13615/3),
    k<6806/3,        hence k<=2268.

No number-field or reflection search is hidden in this conversion.
The range remains too broad for automatic host-by-host work. In
particular a useful cutoff near50 would require much stronger control
of the actual chamber-stabilizer factor or a sharper geometric/analytic
upper; it is not supplied by this argument. The improvement is a
polynomial determinant lower bound and an explicit finite integer
family reduction, rather than the earlier logarithmic collar bound.

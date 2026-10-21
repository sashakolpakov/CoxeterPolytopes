# Odd golden glue: a full base and a polynomial discriminant bound

21 October 2026. New written supporting candidate, pending outside-author
review and root acceptance. No mathematical program, root enumeration,
Lean or TeX process was run. This extends a geometric reduction to an
entire normalized binary family; it does not establish a new polytope or
ambient dimension bound. The six paper editions remain unchanged.

## 1. The family and its precise hypotheses

Put O=Z[phi], phi^2=phi+1, and let E be the even unimodular E8 lattice
over O in the repository's native basis. For c in O set

    B_c=[[-2phi,1],[1,2c]], D=1+4phi*c, b=phi*D,
    L_c=E perpendicular-sum B_c.

The binary basis is X,Y. Assume D>0 and D'<0 at the two embeddings,
and that the ideal(b) is coprime to30. These are the complete additional
hypotheses. Then b is totally positive and a nonunit, L_c is admissible
of physical signature(9,1), and its other embedding is positive definite.
The determinant of B_c is -D. No primality or squarefreeness is assumed.
For the nonunit assertion, b=phi modulo4, whereas every totally positive
golden unit is phi^(2j), whose residue modulo4 is one of1,1+phi,2+3phi.
Here phi^6=1 modulo4; none of those three residues equals phi.

Use the pinned golden PID/unit, complete finite arrangement-angle and
Bugaenko calibration inputs listed below. In particular a primitive
integral reflection normal has square2d, with d totally positive and
d dividing b; each totally positive unit is a square of an O-unit.
The norm-divisor assertion follows directly: evenness gives the square2d,
reflection integrality gives Gm/d integral, and the adjugate and primitive
coordinates give d dividing det(G)=-D. A finite set of ideal divisors
and the unit normalization therefore covers every primitive mirror.

Let C be a full integral reflection chamber adjacent to the calibrated
base F. We prove that this full base exists for every c above, and that,
with N=N(b),

    Vol(C) > (N^2-256)/2^57.                         (1)

Thus the inherited uniform full-chamber volume theorem implies

    finite volume of C => N < 3^1945.                (2)

The earlier polynomial lower bound was proved for the unit-coupled
one-parameter family. The general marked-binary volume bridge already
gave an ineffective finite reduction for arbitrary glue. The extension
here is the full-base theorem and polynomial lower bound for ALL good-
discriminant odd glue; it is not a first finiteness theorem or a practical
enumeration of the enormous rectangle allowed by(2).

## 2. A full-modulus root calculation in the original basis

Define r=X+2phi*Y. Then

    r^2=2b, Gr=(0,D;0), Gr/b=(0,phi^-1;0).

Hence r is a primitive integral reflector. Its orthogonal lattice on
y=0 is the literal M0=OX+E with X^2=-2phi.

We first show that distinct primitive norm2b mirrors are mutually
ultraparallel. Write such a normal m=(x,y;bZ). The E part is divisible
by b by unimodularity and reflection integrality. Its first binary
covector gives y-2phi*x=b*nu with nu in O. Expanding its norm in the
ORIGINAL X,Y basis gives

    1=x^2+D*x*nu+c*b*nu^2+(b/2)*Z^2,
    x^2=1 modulo b.                                 (3)

Here Z^2/2 is integral because E is even. For two such normals,

    (m1,m2)/b = 2*x1*x2 modulo b.

Indeed the remaining terms are D(x1*nu2+x2*nu1), 2cb*nu1*nu2
and b(Z1,Z2), each divisible by b. Their squared angular invariant
A=(m1,m2)^2/b^2 consequently satisfies A=4 modulo b. Every possible
finite arrangement invariant and its difference from4 has prime support
only over2,3,5, except the invariant zero itself. In particular every
possible A-4 is a unit at all primes dividing b. This contradicts the
congruence even at A=0. Parallel distinct hyperplanes would give a
nonzero K-null vector, impossible by conjugate definiteness. This proves
the assertion without replacing the index-two (X,r) frame by a basis.

For a primitive norm2d root with b/d a nonunit, reflection in r gives
(r,m)=b*eta with eta in O. A nonorthogonal finite intersection would
have invariant A=(b/d)*eta^2, divisible by a prime of b/d. Every nonzero
finite arrangement invariant is a unit at that prime, a contradiction.
When (d)=(b), unit rescaling reduces to the norm2b case just proved.
Thus every mirror meeting r is orthogonal to it.

If m is orthogonal to r, then (m,r)=D*y=0 and y=0. Its product with
Y is x, so reflection integrality gives d|x. Unimodularity gives d|m_E.
Primitivity now forces d to be a unit. After unit rescaling m^2=2.
Conversely every norm2 root of M0 is an integral reflector of L_c.
The full intrinsic arrangement on r is therefore exactly the calibrated
Bugaenko arrangement. No omitted norm class can subdivide F.

Local finiteness follows from the finite norm alphabet and bounded root
coordinates at both embeddings for mirrors meeting a compact set. At a
generic interior point of F there is an adjacent chamber on either side
of r. Every nonorthogonal mirror misses the whole base hyperplane, and
the orthogonal ones cut out precisely F. Thus F is an ACTUAL whole facet
of that chamber. In particular its E8 vertex has full actual germ E8+r.

The calibrated base reversal S0 also extends. Its known property
S0X-X in2M0 makes

    S|M0=S0, SY=Y-(S0X-X)/(2phi)

integral in the original basis. It fixes r, is an involutive isometry
and preserves the chosen inward side. The actual full-facet sector
therefore proves S(C)=C. This supplies one involution, not a second
infinite-order symmetry or a nonreflectivity conclusion.

## 3. A growing collar for every divisor norm

Let m be primitive of square2d and not orthogonal or proportional to r.
Write e=b/d. Reflection integrality gives

    m=(x,2phi*x+d*nu;dZ),
    e*x^2+D*x*nu+c*d*nu^2+(d/2)*Z^2=1.              (4)

With gamma=y/phi=2x+(d/phi)nu, its invariant is A=e*gamma^2.
Subtract four times(4) to obtain the second divisibility

    A-4=d*xi, xi=phi^-1*nu^2-2Z^2 in O.             (5)

The cancellation uses b/phi^2-4c=phi^-1; it does not use c being a
unit or lying in the earlier family. The physical mirrors are
ultraparallel by Section2, so A=4cosh(ell)^2>4. At the definite
embedding Cauchy gives 0<A'<4 strictly. A nonzero product cannot vanish
at the other embedding; equality in Cauchy would make the K-vectors
proportional, already excluded.

The nonzero integral norms in A=e*gamma^2 and(5) yield

    A>N(e)/4, A-4>N(d)/4.

Their product gives 256*sinh(ell)^2*cosh(ell)^2>N, equivalently
64*sinh(2ell)^2>N. Every nonorthogonal mirror stays at global distance
strictly greater than rho=(1/2)asinh(sqrt(N)/8) from the entire base.

For x in F on the future unit hyperboloid, the full inward flow
cosh(t)x-sinh(t)r/sqrt(2b), 0<=t<=rho, stays in C. Orthogonal wall
inequalities retain their signs; no other mirror can cross that segment
before its global distance from the base. This checks ALL arrangement
walls, rather than a selected finite seed. Injectivity of normal flow
and its cosh(t)^8 volume element give

    Vol(C)>=Vol(F)*integral_0^rho cosh(t)^8 dt.

The pinned calibration supplies Vol(F)>2^-38 from its radius1/2 E8
ball sector. Also exp(8rho)>N^2/256, so the integral exceeds
(N^2-256)/2^19. For N>=16 this proves(1) directly; for N<16 its
right side is negative and positive chamber volume suffices.

The general full-reflection volume bridge applies to EVERY admissible
rank-ten O-lattice, not merely the earlier unit-coupled family. It gives
Vol(C)<U<3^3851 when C has finite volume. Since2^57<3^38,
N^2<2^57 U+256<3^3890 and(2) follows. Normality, finite stabilizer
and congruence spectral inputs have their explicit proof home in that
bridge. No claim that an arbitrary reflection subgroup is congruence
or satisfies the same volume upper bound is made.

## 4. The finite reduction is an algorithm, not an executed census

The old marked-binary basis change can be written particularly simply:
Y_new=X+phi^3 Y keeps (X,Y_new)=1 and sends

    c -> phi^6*c+phi^2, D -> phi^6*D.                (6)

Its inverse is integral. Iterate(6) or its inverse until1<D<phi^6.
This halts because phi^6>1. Equality D=1 at any step is impossible:
it would make the algebraic number D rational1, contradicting D'<0.
The normalized c is then physically in(0,phi^2), and

    0<c' = phi*(1+N/D)/4 < phi*(1+3^1945)/4.

For c=A+Bphi, the relations B=(c-c')/sqrt5 and A=c-phi*B bound both
integer coefficients explicitly. Enumerating this finite rectangle,
then imposing the embedding, discriminant and coprimality tests, is a
terminating complete necessary host-list algorithm at the stated input
boundary. It neither decides reflectivity nor supplies actual examples.
The range is enormous; no such enumeration is proposed or executed.
Bad-prime discriminants and the even-glue case are outside this theorem.

## 5. A separate small test module

Take c=4. Then D=1+16phi, b=16+17phi, N(b)=239. The latter is prime:
239 is not divisible by2,3,5,7,11 or13, the only primes below its
square root. The signs are D>0 and D'=17-16phi<0. Thus this literal
module has the complete actual Bugaenko base and the inherited chamber
involution just proved. This is a starting geometry, not a finite-volume
completion. Its discriminant norm is absent from the strictly increasing
old unit-coupled sequence: the adjacent values131 and319 straddle239.
It is not a re-presentation of any member of that sequence. The odd
valuation at the prime of norm239 also distinguishes its K-squareclass
from the specifically closed norm11,19,31,55,59,79,131,319,811,2099
modules. No statement about all later forms or literature novelty follows.

The next actual edge retains E8 and releases r at p0=X. Its axis is
p(t)=X-t*r, 0<t<1/sqrt(D). The uniform mixed-E8 collar implies that
no mixed root meets this binary hyperbolic line. Equivalently, finite
first-corner classification retaining E8 permits only an independent
binary simple opposite. The prime alphabet therefore leaves exactly
two necessary binary profiles:

    norm2: m=((u-v)/(2phi),-v;0), u,v>0,
           u^2=D*v^2-4phi, u-v in2O,
           collision t=u/(D*v);
    norm2b: m=((D*w-z)/2,-phi*z;0), w,z>0,
           z^2-D*w^2=4, z-w in2O,
           collision t=w/z.                         (7)

These are ORIGINAL-module conditions. In the second row integrality
of the two reflection covectors is exactly the displayed parity;
in the first row norm2 makes reflection integrality automatic. The
first binary covector is -u, and the norm identity is
2phi*m^2=D*y^2-(2phi*x-y)^2, which derives both rows. Necessary
conjugate strips are v'^2<=(-4phi')/(-D') and w'^2<=4/(-D').

No solution, upper collision or firstness is asserted. An exploratory
bounded search can earn only a mirror witness or a stated finite-box
absence, never all-height absence. Once a finite upper is supplied,
collision comparison and these conjugate strips give a complete finite
firstness domain; all normals, covectors, minima and ties must then be
checked. A whole-chamber compactness or nonreflectivity proof is still
required. This construction lane remains open.

## Bound inputs and external facts

All paths below are under results/ unless stated otherwise.

| Input | SHA256 | Use |
|---|---|---|
| bugaenko_facet_binary_extension_and_collar_audit_20261020.md | `774e87f7abb33680a06879a17a5be9f06e2218bf12a95023bcff5c7046ba1227` | Golden ring facts, normalized odd glue, mixed-E8 collar, actual calibrated reversal |
| golden_e8_binary_family_full_base_facet_20261020.md | `ea5d96e331ac91cb9efc457db2ddaa24bac3341e9e6781038150765895675afd` | Complete finite arrangement-angle list and its good-prime norm support; no transfer of its unit-basis argument |
| golden_e8_binary_family_discriminant_collar_and_finite_cutoff_20261020.md | `841a4b7f9e75cdd49b6384920439d998cc401f00bbf75aed92444892c851d425` | Calibrated base-volume lower bound and elementary collar integration only |
| bugaenko_facet_full_reflection_volume_bridge_addendum_20261020.md | `e61a5531051aac247382b9e1ec2f4536bb7a9b89fba3352febf7fe1bf01cf6a3` | Uniform upper for full chambers of arbitrary admissible rank-ten O-lattices and marked-basis reduction |
| review_bugaenko_full_reflection_volume_bridge_batch_guard_20261020.md | `56be6455426bb30d66076ba7e197aedc176c622d7c7f735f56988a09e2968761` | Independent review of that precise general analytic interface |
| data/known_examples/bugaenko_compact_h8.json | `c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c` | Native E8 and actual Bugaenko calibration |

Finite Coxeter/arrangement-angle classification, PID/units, the complete
calibration and the named analytic inequalities remain explicit external
or inherited inputs. The original-coordinate full-modulus calculation,
all-divisor collar extension and norm239 starting geometry above are new
arguments, still awaiting independent full review. No runtime or Lean
credit, global ceiling change or new example is claimed.

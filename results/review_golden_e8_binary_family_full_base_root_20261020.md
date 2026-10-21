# Root review: full bases throughout the unit-coupled golden family

20 October 2026. Complete written review **PASS**, binding
`golden_e8_binary_family_full_base_facet_20261020.md`, SHA256
`ea5d96e331ac91cb9efc457db2ddaa24bac3341e9e6781038150765895675afd`.
This is a structural construction theorem, not uniform nonreflectivity.
No scientific job, root enumeration, Lean or TeX was run for this review.

I read the entire source and reconstructed each algebraic step. For
n=2k+1>=-1, b=4+phi^n has conjugate 4-phi^(-n)>0 and determinant
-phi*b; the conjugate binary form is positive definite. The three
residue calculations show gcd((b),(30))=1 for every k, including k=-1.
The argument does not use a conjecture about prime or squarefree terms
in this sequence.

The primitive-norm argument is valid over the explicitly imported
golden PID. Evenness gives q=2d; reflection integrality and a Bezout
combination of primitive coordinates give d|det(G), hence d|b.
If b/d is a unit, it is totally positive and therefore a square of
an O-unit; the normal can indeed be normalized to q=2b. This covers
all repeated divisor exponents rather than assuming two norm classes.

I reconstructed the finite-angle list. A finite-order product has
double-angle cosine in Q(phi); the rational cases and the conductor5/10
cases give exactly the eight displayed x=4cos^2(theta). The listed
norms of x and x-4 are correct. In particular the nonnative 2pi/5
and3pi/10 cases are retained. Every nonzero x and every x-4 in that
list is a unit at primes of b. The positive conjugate form excludes
nonzero K-null vectors, so distinct parallel rational mirrors cannot
be an omitted case.

For a full-b normal v and a proper-divisor normal m, integrality of
the v reflection gives (v,m)=b*gamma. The resulting
x=(b/d)*gamma^2 contradicts any prime divisor of b/d at a nonzero
finite angle. Orthogonality is still allowed. No valuation-one or
squarefree assumption enters.

The new whole-modulus argument also checks. The binary basis (r,c)
has unit determinant. A norm2b reflector has E8 coordinates b*z;
integrality against c and coprimality with2 force its c coefficient
to be b*nu. Dividing its norm equation by2b gives

    1=lambda^2-b*lambda*nu+b*nu^2+b*z^2/2.

Evenness makes the last coefficient integral. Thus lambda^2=1 mod b,
also in a nonreduced quotient. For two full-b normals, pairing divided
by b is 2lambda1lambda2 modulo b, so its square is4 modulo b.
The x-4 norm list excludes every finite intersection. Together with
anisotropy, this proves that distinct full-b mirrors are ultraparallel.

An orthogonal reflector has y=0. Its pairing with Y is the unit
a times x; together with unimodularity of E8 this gives d|x,z.
Primitivity excludes every nonunit d. Conversely each intrinsic
norm2 reflector extends integrally to the whole lattice. The pinned
Bugaenko calibration therefore identifies the complete facet, not
only a root subset. One can also check the possible intrinsic dyadic
norm directly: q=4 and z=2w would imply phi*x^2=0 mod2, hence x=0
mod2 in F4 and the normal would not be primitive. This is consistent
with the imported norm2-only calibration.

Finally S0X-X in2M makes the displayed extension integral for every
unit u. Fixing r gives exactly SY=Y-u(S0X-X)/2. The extension is an
isometric future involution and preserves the chosen inward facet
germ. Local finiteness at a generic facet interior point then makes
it a symmetry of the full adjacent chamber. This step uses the full
facet result, not merely an isometry of the lattice.

The conclusion holds for every k>=-1 with the stated calibrated
and classical inputs. It supplies one actual chamber involution.
Proper-divisor roots may change transverse edges; a second independent
involution and uniform nonreflectivity are unproved. There is no new
example, dimension ceiling or commensurability-class exclusion.

The author and I discussed the initial squarefree angle-isolation
route before this note. The whole-modulus proof was supplied by the
author, and every step above was independently reconstructed here.

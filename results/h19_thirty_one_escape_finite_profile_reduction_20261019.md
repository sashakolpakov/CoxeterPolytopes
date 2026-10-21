# Current H19 edge: a finite contact domain and a reduced Pell equation

19 October 2026. **New paper reduction, awaiting independent review.**
No scientific program, import, root search, matrix calculation program,
census, Lean or build ran. Certificate/source reading and file metadata
inspection were the only operations. No new admissible mirror or
finite-volume completion is claimed. All pinned sources remain unchanged.

## 1. The actual thirty-one-wall escape

Use the fixed native form

    G=[[-2phi,phi],[phi,2]]+C_E8,
    phi²=phi+1, b=4+phi, D=1+5phi=phi*b,
    u=10+17phi, t=12+19phi, 7t²-Du²=1.

The current certificate is
`h19_thirty_one_wall_seed_20261019.json`, SHA256
`5a3b25aa139793f304d535f56d6033ad103d88cffe07c20278bc7b42c2712d57`;
its fixed checker is SHA256
`e81cd38a68c5262acc6421c16f8d8114e2480f5f453edde503a32de30710c76e`.
The coordinator's generation and separate replay PASS are recorded in
`h19_thirty_one_wall_execution_root_20261019.md`; they are not new runs
by this note's author.

The new escape starts at vertex3=p_s and frees old3, not the binary
walls freed by the preceding K/L extensions. In native coordinates,
with omega_j=C_E8^-1 e_j,

    V=(1,-phi/2,(Du/(14t))*omega1), V²=-D/(14t²),
    delta=(0,0;omega2-(10/7)omega1),
    delta²=10/7, V.delta=0,
    E=V+h*delta, h=(-42+27sqrt(5))/220.

The first hit is q. The certificate gives
E²=6227/4840-691sqrt(5)/1210>0 and all31 products nonpositive.
The current common white system is

    {old0,old1,old4,old5,old6,old7,old8,s}
     = A4(old4--old5--old6--old7)+4A1.

The four isolated roots are old0 (binary), old1 (E8 index0),
old8 (E8 index7), and s. Every one has norm two. None of the
prior A6-edge profile lists is a complete list for this new edge.

## 2. Complete finite short endpoint profiles:63 before congruences

Let R be a new canonical norm-two integral reflector which closes this
edge at a finite point. Every product with a white root occurs at that
new finite corner. Since products lie in K=Q(sqrt5), the only nonright
finite short/short contacts are -1 and -phi, labels3 and5. This follows
from the degree of2cos(pi/m): if it lies in K, the finite degree-two
possibilities are m=3,4,5,6, and the values at4 and6 lie in the
different quadratic fields Q(sqrt2),Q(sqrt3). Thus no label cutoff is
assumed.

Write its nonnegative contact magnitudes on W as a profile. The finite
Coxeter classification gives the following complete own list:

* If R touches the whole A4, it contacts exactly one of its nodes,
  with label3. At either endpoint it may also contact zero, one,
  or two of the four isolated roots, all with label3. The connected
  components are A5,A6,D7 respectively. At either internal A4 node
  it may contact zero or one isolate, giving D5 or E6. Two A4
  contacts create a cycle; three isolates or a double fork are not
  finite. A5-or-larger connected type cannot contain a5-bond.
  This gives2*(1+4+6)+2*(1+4)=32 profiles.
* If R does not touch A4, it may touch zero isolates (one profile),
  one isolate with label3 or5 (eight profiles), two isolates with
  labels(3,3),(3,5),(5,3) (eighteen profiles), or three isolates
  all with label3 (four profiles). Four3-contacts are affine D4;
  two5-contacts or a5-contact with three leaves are not finite.
  This gives1+8+18+4=31 profiles.

Total:63. These are necessary finite-corner profiles, not63 admitted
reflectors. Integrality, norm, all31 native pairs, orientation and the
actual endpoint must still be checked for each candidate.

## 3. Exact profile Pell equation independent of large coefficients

Write R=(x,y,z). Let beta=-R.old0, epsilon=-R.s. Let eta_j be
the contact on old(j+1) for j in J={0,3,4,5,6,7}. All these
entries belong to {0,1,phi}, as specified by Section2. Write

    z=T*omega1+U*omega2+sum_(j in J) eta_j*omega_j,
    R1=sum_j Cinv[1,j]*eta_j,
    R2=sum_j Cinv[2,j]*eta_j,
    Q0=sum_ij eta_i*Cinv[i,j]*eta_j,
    B=7R2-10R1.

Here Cinv=C_E8^-1 has integral entries. The two remaining white
equations are

    2y=beta-phi*x,
    14t*T+20t*U+t*R1=Du*x-epsilon.                (1)

All x,y,z,T,U are integral over O=Z[phi] for an integral native
root. In particular t divides Du*x-epsilon. Multiplying by u and
using Du²=7t²-1 gives t divides x+u*epsilon. Thus

    n=(x+u*epsilon)/t lies in O.                  (2)

Put A=10U+B. The norm equation R²=2kappa is equivalent to

    A²-5D*n² = B²+140kappa-35(beta²+epsilon²)
                       -70Q0+5R1².             (3)

Here kappa=1 for a short root. To derive it, completing the square
in the E8 block gives

    z²=((Du*x-epsilon)/t)²/14
          +(10/7)U²+(2/7)B*U+Q0-R1²/14.

The binary square is beta²/2-D*x²/2. Multiply their sum by14t²,
use D*u²=7t²-1, and complete x to x+u*epsilon. Completing U to
A=10U+B gives exactly(3). No profile-dependent determinant or
floating-point square test is required.

Conversely a field solution A,n of(3) reconstructs a NORM candidate by

    x=t*n-u*epsilon,       U=(A-B)/10,
    T=(Du*n-7t*epsilon-20U-R1)/14,
    y=(beta-phi*x)/2.                             (4)

Both signs of A must be retained before checking T,U,y integrality.
The integrality conditions in(4) and x>0 remain necessary; a field
square alone is not a root. The finite profile list does not impose
orthogonality on any unlisted current wall.

## 4. The pure short pencil is obstructed at every height

For the zero profile beta=epsilon=eta=0. Equation(1) and the first
equation give x=2v and t divides Du*v. Norms are N(t)=11,
N(u)=-19,N(D)=-19, so t and Du are coprime. Hence v=t*n0,
x=2t*n0. The remaining equations reduce to

    7T+10U=Du*n0,
    5U²-D*n0²=7.                                 (5)

Modulo the ramified prime(sqrt5), the residue field is F5,
phi=3 and D=1. Equation(5) would force -n0²=2, hence n0²=3
in F5. This is impossible. Thus no integral norm-two reflector
orthogonal to all eight current whites exists at ANY height.
The short finite profile domain therefore reduces63 to62 before
any bounded search. This excludes just this orthogonal pencil,
not the nonorthogonal finite profiles or other completions.

## 5. The long endpoint pencil remains genuinely separate

The reviewed H19 norm alphabet is {2,2b}. A canonical long root of
norm2b cannot have a finite nonright bond with a canonical short root:
the cyclic-square field test permits only m=3,4,5,6,10, while the
rational norms of its required squared raw products are respectively
19,76,19,171,95, none rational squares. Thus a long root closing
the present edge at a finite vertex must be orthogonal to ALL
eight white short roots.

The reflector integrality condition forces its E8 covector divisible
by b. Write z=b*(T0*omega1+U0*omega2). As above x=2t*n0 and
y=-phi*t*n0, but now the equations are

    5b*U0²-phi*n0²=7,
    T0=(phi*u*n0-10U0)/7 lies in O.               (6)

This is not excluded by the short mod5 proof: modulo(sqrt5) it
requires n0²=1. At the prime(b), phi=15 in F19, it requires
n0²=16. Both have local solutions. No all-height long-pencil
obstruction is asserted.

One small exact field solution illustrates why norm equations alone
are insufficient. Set n0=1+6phi, U0=-phi². Then(6) has
T0=27+43phi integral. But the native product with RETAINED old3
is -b*U0=b*phi²>0. It is not an admissible appended mirror. Changing
only the sign to U0=+phi² gives T0=(169+281phi)/7, not integral.
Thus neither sign supplies the required next wall. These are direct
paper substitutions, not a bounded enumeration or a root-height
completeness assertion.

## 6. Actual endpoint test and bounded-search interface

For every short-profile norm candidate, the exact products on the
current edge are

    (V,R)=-D*(x+u*epsilon)/(14t²),
    (delta,R)=(10U+8eta3+6eta4+4eta5+2eta6+5eta7)/7.

A candidate with positive slope cuts the edge at

    tau=-(V,R)/(delta,R).

It must have0<tau<h, a future timelike endpoint V+tau*delta,
all31 old inequalities nonpositive, exactly the intended finite
active diagram, and all new/old pairs native and acute. These are
geometric acceptance conditions; they do not follow from(3).
The endpoint is timelike precisely when tau²<D/(20t²).

An authorized bounded worker can enumerate the62 profiles, integral n
in a declared physical/conjugate box, and test the field square in(3).
Both A signs must be reconstructed by(4). A separate bounded branch
can test(6). All current31 roots must be recovered from their pinned
native covectors/coordinates, not from an obsolete28-wall list. Any
survivor must carry its full31-pair/endpoint ledger. Preservation of
the compact cap facets and prior named finite vertices must be checked
separately before reusing their counts. No enumeration box is chosen
or run here; the coordinator owns resource gates and scheduling.

The source-only achievement is a complete current-edge finite contact
domain, an all-height obstruction for its short orthogonal profile,
and the explicit common reduced Pell identity. It is not a new
compatible mirror, root census, lattice nonreflectivity theorem or
finite-volume H9 example.

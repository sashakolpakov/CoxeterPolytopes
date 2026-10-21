# Golden isolated-five B stars: the invariant-lattice interface

20 October 2026. Paper proof only; no scientific execution. This is a
supporting exclusion of two contact positions, not a dimension bound.

Put phi=(1+sqrt5)/2, O=Z[phi], L=O^(n+1), and

    f(t,x)=-phi*t^2+sum x_i^2,   n>=8.

Write O(M) for K-linear isometries preserving a full O-lattice M.
The following lattice assertion is independent of a reflectivity cutoff.

## Exact lattice assertion

Let H contain all signed spatial permutations and the reflection in
either norm-two vector

    c_1=(phi^-1;phi,0,...,0),
    c_5=(phi^2;phi,phi,phi,phi,phi,0,...,0).

For EVERY H-invariant full O-lattice Lambda,

    O(Lambda)=O(L).                                      (1)

Here is the complete proof. The ring O is a PID: its quadratic
Minkowski bound sqrt5/2 is less than2. Let A be the O-algebra generated
by H. Signed coordinate flips give 2E_ii in A for i>0; multiplication
by spatial permutations gives ALL spatial 2E_ij. For a conjugate cap
c, its time coordinate and every nonzero spatial coordinate are units.
The norm-two formula R=I-c B(c,.) therefore lets us subtract the
known spatial entries from (2E_ii)R and R(2E_ii), obtaining 2E_i0
and 2E_0i. Permutations give these for every i. Finally R(2E_i0),
after subtracting the known spatial-time entries, gives 2E_00.
Thus 2 End_O(L) is contained in A. If I is the coordinate-content
ideal of Lambda, then 2IL is contained in Lambda, which is contained
in IL. Principal-ideal scaling reduces to

    2L <= Lambda <= L.

Now O/2O=F4, with phi^2=phi+1 and phi^3=1. Both cap types reduce to
the same formula, with |T|=1 or5 respectively:

    (R_T-I)y=(phi^2*t+sum_(i in T) x_i)
                  *(e_0+phi^2*sum_(i in T)e_i).

All T of that size occur. The common kernel of these functionals is
D=<(phi,1,...,1)>: for |T|=5, differences of two T's force every x_i
equal, and then the displayed functional gives their value; |T|=1
is immediate. The orbit span of the displayed image vectors is the
hyperplane E={phi^2*t+sum x_i=0}=D^perp. For |T|=5, differences give
the whole spatial augmentation hyperplane, and one image vector adds
a time direction. We use n>5 here.

If an invariant subspace U is not contained in D, one functional is
nonzero on U, so its image vector belongs to U. Permutations then
give E<=U. Consequently U is exactly 0,D,E, or V=L/2L. There are no
unlisted lifts: Lambda is the full inverse image of U. Besides L and
2L, the possibilities are

    Lambda_D=2L+O*d,  d=(phi,1,...,1),
    Lambda_E={y in L:B(d,y) in 2O}=2 Lambda_D^*.

It remains to identify their orthogonal groups, rather than assume
that commensurable lattices have equal groups. The exact norm is

    f(d)=n-1-2phi.

Its dyadic valuation is0 for n even and exactly1 for n odd: in the
second case (n-1)/2-phi is nonzero in F4. Therefore

    {y in Lambda_D/2 : f(y) in O}=L.                    (2)

Indeed write y=z+b*d/2. If b is not divisible by2, the last term in
f(y)=f(z)+b B(z,d)+b^2 f(d)/4 has negative dyadic valuation, while
the first two are integral. This proves O(Lambda_D)<=O(L).
Conversely d mod2 is the unique characteristic vector satisfying
B(d,y)=sqrt(f(y)) over F4, so O(L) fixes it and preserves Lambda_D.
Thus equality holds. Duality gives the same group for Lambda_E;
L and2L are immediate. This proves (1), including every lift.

## Arithmetic-star consequence and the precise external input

Suppose a compact arithmetic Coxeter polytope has a WHOLE B_n vertex
(path 1,...,n, terminal4 at n-1,n), and an actual outside wall has
only a label5 contact, at a LONG node k. Normalize old long roots to
norm2 and the short root to norm1. The cap can have norm2 and spatial
part phi*(e_1+...+e_k); its time coefficient square is

    a_k=k*phi^2-2.

This full-rank Lorentz star has cyclic field Q(phi). Any larger
arithmetic ground field would have a nonphysical embedding fixing
this entire indefinite Gram matrix, contradicting conjugate
positivity. Hence the actual ground field is exactly Q(phi).
The other embedding requires k*(phi')^2<2, so k<=5.
For k=1, a_1=phi^-1; for k=5, a_5=phi^5. Rescaling the time axis
puts these two stars in the displayed L, with caps c_1,c_5.
Arithmeticity supplies a full invariant O-lattice for the entire
reflection group Gamma (e.g. the sum over finitely many coset
representatives of a lattice preserved by an arithmetic subgroup).
Assertion (1) therefore gives Gamma<=O(L) in these exact coordinates.

**External cutoff G:** the FULL integral reflection group of this
literal L has no compact chamber in dimension n>=8. With G, neither
position k=1 nor k=5 can occur: a full O(L)-arrangement chamber
refines the Gamma chamber and is contained in the assumed compact
polytope. This uses actual group containment, not commensurability.

The literal cutoff is stated in Kellerhals, *Hyperbolic Orbifolds of
Minimal Volume* (2014), section4.1.1, printedp474, equation(4.3) and
its following sentence, attributing it to Bugaenko1984:
https://homeweb.unifr.ch/kellerha/pub/CMFT.pdf . The primary Vinberg
1985 article, Introduction printedp38 and section5.6 printedp64,
discusses the same unimodular golden family and cites Bugaenko[71]:
https://www.mathnet.ru/php/getFT.phtml?jrnid=rm&option_lang=eng&paperid=2140&what=fullteng .
The original Bugaenko1984 theorem was not retrieved in this audit;
G remains an explicitly named imported input, not an own proof.

The unresolved positions are k=2,3,4, with time square coefficients
2phi,3phi+1,2phi^3 respectively. Positions2/4 have the determinant-two
class and position3 the norm-five class; (1) does not transfer them
to L. Dotted caps, other contacts, other finite types, and whole B_n
energy are also not settled. No new uniform ceiling is claimed.

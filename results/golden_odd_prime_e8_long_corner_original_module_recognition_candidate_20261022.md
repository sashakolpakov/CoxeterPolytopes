# Golden odd-prime modules: an E8 plus long corner determines an integral map

22 October2026. New complete manual supporting candidate by root, pending
full other-author mathematical review. No scientific worker, enumeration,
symbolic computation, test, Lean or TeX ran. No existing proof, paper,
receipt, guide or bound is changed. This is a module-recognition lemma for
the construction lane, not a claim that a new corner exists.

The first outside read found a missing parity hypothesis: the binary
diagonal m must be EVEN for the quadratic reduction modulo2 below.
The statement is corrected accordingly, retaining both intended hosts
m8 and m10. The previous8247-byte candidate (SHA256
5923505e18476dd329a9279e5577131dedd33c279c180bc1bcd26b34e136434c)
is byte-verified in
`source_archives/golden_e8_corner_before_even_binary_scope_20261022.md.gz`.
No theorem for odd m is claimed.

## 1. Statement and exact boundary

Put O=Z[phi], phi²=phi+1, K=Q(phi). Let m>=2 be an EVEN integer,

    D=1+2m phi, b=phi D=2m+(2m+1)phi,
    L=E8(O) perpendicular-sum [[-2phi,1],[1,m]].

Use the original basis (X,Y,e0,...,e7). The E8 Gram has diagonal2
and products-1 on02,13,23,34,45,56,67. Assume N(b)=4m²-2m-1
is a rational prime. This includes the literal m8/norm239 and
m10/norm379 hosts; it is not an exhaustion of arithmetic forms.
The original isolated long root is r0=X+2phi Y, with r0²=2b and
r0.X=0. Fix the future cone containing X.

**Module theorem.** Suppose actual original vectors e'_0,...,e'_7 have
the same ordered E8 Gram, and w is an original primitive integral
reflection root of square2b perpendicular to them. There is a UNIQUE
future integral isometry G of the ORIGINAL L satisfying

    G(e_i)=e'_i, G(r0)=w.

The vector N=G(X) is the unique future vector of norm-2phi on the
one-dimensional perpendicular of e'_0,...,e'_7,w. In particular it is
original integral and

    G(Y)=(w-N)/(2phi) is original integral.

The word actual is not needed for the module conclusion: integral vectors,
the stated Gram, and an integral long reflection suffice. Actual full
chamber germs are needed only for the geometric corollary in Section5.

The original-module convention and calibrated E8 data are the same as
`golden_odd_glue_full_base_and_discriminant_reduction_20261021.md`
(SHA256 b107f0de4f76da942100e0c82867cd6d4eeee40d833304fee67cebd84d3e41ce)
and `data/known_examples/bugaenko_compact_h8.json`
(SHA256 c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c).
Paths are under results/ except the displayed data path. The proof below
does not assume a reflection chamber is finite-volume, a rational
E8-matching map is integral, or an active subspace is an integral summand.

## 2. Elementary ring and integral splitting facts

The ring O is norm-Euclidean. Indeed round the two rational coefficients
of any element of K in the basis1,phi to nearest integers. The remainder
a+bphi has |a|,|b|<=1/2 and

    N(a+bphi)=a²+ab-b²,
    -5/16 <= (a+b/2)²-5b²/4
            =N(a+bphi)
            =5a²/4-(b-a/2)² <=5/16.

Multiplicativity gives Euclidean division with strictly smaller absolute
norm. Thus O is a principal ideal domain. Also O/(2) is F4, since
T²+T+1 is irreducible over F2. The odd prime ideal (b) is distinct
from (2). The displayed L is even: v² belongs to2O for every v in L.

Let S be the O-span of the e'_i. Its E8 Gram has determinant1, so
its inverse is integral. For every v in L the orthogonal projection
onto S has coefficients in O, because all eight products (v,e'_i)
are integral. Hence

    L=S perpendicular-sum M, M=S-perpendicular intersect L.

This is an ACTUAL integral direct sum. M is free of rank2 over the
principal ideal domain O. Its physical signature is(1,1), and its
other real embedding is positive definite. The latter follows also
directly from the displayed original binary Gram: -2phi'>0 and
-D'>0 for m>=2.

The primitive vector w in M extends to an O-basis (w,z) of M. With
the ordered e'_i this gives a basis of L. Let delta be its determinant
relative to the original ordered basis; delta is a unit of O. Since
the original full determinant is-D and det(E8)=1, the binary Gram
in the basis (w,z) has determinant-D*delta². This records the unit
square; its value is not silently set to1.

## 3. The long covector and an integral vector of the correct time norm

Integrality of the actual reflection in w implies

    (w,L) is contained in bO.

For completeness, its original matrix is I-w*(G_full*w)^t/b. At the
prime ideal (b), some coordinate of w is a unit: otherwise w²=2b
would have valuation at least2 instead of1. The corresponding row
of this integral reflection forces every covector entry into bO.
At all other primes division by b causes no issue. This proves
whole-module divisibility, not just divisibility on the active vectors.

Write (w,z)=b eta and z²=c, eta,c in O. The determinant identity is

    b*(2c-b eta²)=-D delta²,
    2c-b eta²=-delta²/phi.                      (1)

Consequently the integral vector

    N0=(phi/delta)*(2z-eta w)

is perpendicular to w, with

    N0²=(phi²/delta²)*(4c-2b eta²)=-2phi.

It is also perpendicular to S. Choose N to be N0 or-N0 according
to the prescribed future cone. This proves existence of an ORIGINAL
integral future vector of precisely the required norm. The perpendicular
line is timelike; any two vectors on it with that norm differ by sign,
so the future choice is unique. Changing the auxiliary basis z changes
no final N.

## 4. The dyadic glue is forced, not guessed

Reduce M modulo2. Its determinant-D delta² is a unit modulo2, since
D=1+2m phi. The bilinear form there is alternating and nondegenerate:
L is even and the rank is2. Both w and N have nonzero reductions,
since their quadratic values v²/2 are respectively b and-phi, units
modulo2. Their orthogonality therefore implies

    N_bar=lambda*w_bar for some lambda in F4-star.

The quadratic function q(v)=v²/2 is well-defined modulo2: replacing
v by v+2u changes q by2(v,u)+2u². Its two values satisfy

    q(N)=-phi=phi=b=q(w) modulo2.

Homogeneity gives lambda²=1, hence lambda=1 in the field F4.
Therefore N-w belongs to2M. The vector

    Y'=(w-N)/(2phi)

is ORIGINAL integral. Its three binary products are exactly

    N²=-2phi, (N,Y')=1,
    Y'²=(b-phi)/(2phi²)=m.                      (2)

Map X to N, Y to Y', and each e_i to e'_i. Equation(2) and the
orthogonal E8 block prove that this K-linear map preserves the full
original form. All ten images are original integral. Its determinant
squares to1 by the Gram identity, hence is+1 or-1, a unit of O;
the adjugate proves its inverse integral too. Thus G(L)=L, not
merely containment in a same-rank submodule. G preserves the future
cone because it sends X to future N, and G(r0)=N+2phi Y'=w.
The unique future N and the fixed ordered E8 frame prove uniqueness.

This argument explains why matching a finite E8+long diagram can be
upgraded here to an ORIGINAL-module map. Both determinant square and
dyadic glue were checked. The same conclusion must not be asserted
for a different norm alphabet, invariant lattice or field without
reproving those facts.

## 5. Actual full-germ corollary and its strict limit

Suppose the original full outward germ (-e_i,r0) and another FULL
outward germ (-e'_i,w) belong to the SAME chamber C of the arrangement
of all integral reflection mirrors of L. The map G above preserves
that entire arrangement, because it and its inverse preserve L and
the form. It maps the first actual outward local chamber to the second.
Chambers of the arrangement have disjoint interiors, so G(C)=C.
This conclusion needs FULL actual germs, not merely a finite subdiagram
or a feasible point against a chosen finite wall list.

There is NO automatic infinite-order conclusion. If a separate exact
argument proves G has infinite order, then its powers form an infinite
subgroup of the chamber stabilizer. Let W be the group generated by
all integral reflections. Its fundamental-chamber action is simply
transitive on chambers, so W meets that stabilizer trivially. Since
W is normal in the integral orthogonal group, the powers of G give
infinitely many different cosets, and W has infinite index. This is
the usual nonreflectivity obstruction. Finite order of G, a new actual
germ, or any covolume statement is not established by this lemma.

The module theorem applies PROSPECTIVELY if the norm379 fourth-edge
computation eventually produces a full E8+norm2b germ. That computation
has not run here and no such output is assumed. The lemma may simplify
its later period map; it does not itself identify a fourth wall, rule
out the host, construct a polytope or change any dimension bound.

# A biquadratic admissible unit form and its complete first B8 edge

20 October 2026. Own paper preflight and complete finite interval.
No worker, root-prefix search, Lean or TeX execution; no compact example.

Put `s=sqrt2`, `phi=(1+sqrt5)/2`, `u=phi-1`, and

    K=Q(s,phi), O=Z[s,phi], a=u+s,
    N9=<-a>+I9, d=2-s, X=2+s.

## 1. Exact arithmetic preflight

The coprime quadratic discriminants8 and5 give the full product
integral basis (1,phi,s,s*phi) and field discriminant1600.
Minkowski's ideal bound is `(4!/4^4)*sqrt(1600)=15/4`. There is
no prime ideal of norm2: modulo the ramified s, phi has irreducible
polynomial `z^2+z+1`, so the residue field is F4. Nor is there a
prime ideal of norm3: the Frobenius acts nontrivially on both
quadratic factors, giving residue degree2 and prime norm9.
Consequently every ideal class contains the unit ideal; O is a PID.

The four embeddings and the coefficient signs are explicitly

|embedding|(s,phi)|a|
|---|---|---|
|0|(+s,+phi)|u+s>0|
|1|(-s,+phi)|u-s<0|
|2|(+s,1-phi)|s-phi<0|
|3|(-s,1-phi)|-s-phi<0|

Thus N9 is admissible and indefinite at EXACTLY one real place.
Moreover `Norm_(K/Q(phi))(a)=u^2-2=-phi`, so Norm(a)=-1.
The units -1,alpha=1+s,phi,a have independent sign vectors
`1111,0101,0011,0111`. By the unit rank3 theorem these exhaust
unit classes modulo squares, so every totally positive unit is a
square. These standard elementary number-field inputs are stated
explicitly rather than replaced by a bounded unit search.

The sole dyadic prime is (s), with `(2)=(s)^2`. Primitive reflection
integrality in this unimodular lattice gives q dividing2. Total
positivity and the preceding unit signatures leave exactly the
necessary normalized norm pool

    q=1,2,d.

All three norms divide2, so an integral vector of any displayed
norm is reflective. The value d is the valuation-one representative
of the same unit-square class as X (`X=alpha^2*d`).

For any coordinate x=A+Bphi+s(C+Dphi),

    Tr(x^2)=U(A,B)+2U(C,D),
    U(A,B)=4A^2+4AB+6B^2.                       (1)

The smallest nonzero trace is4, attained only by +/-1. At time0,
norm1 roots are therefore the coordinate unit vectors. Norm2 has
trace8, giving only two coordinate +/-1 entries or a single +/-s;
the latter is imprimitive and gives an already present unit mirror.
Norm d also has trace8, but the same possibilities all have rational
norm2, so none has norm d. Thus the FULL zero-height arrangement is
B9, with actual simples

    beta_i=-e_i+e_(i+1), i1,...,8; beta9=-e9.

The initial point is p=(1;0^9). Local finiteness is the usual
bounded-integral-matrix argument using all four embeddings; the
three nonphysical forms are positive definite.

## 2. Finite upper and complete endpoint profiles

Release beta1, retaining the whole B8. Since `phi^2-a=d`,

    r=(1;phi,0^8), r^2=d

is an integral primitive reflector. It hits the coordinate ray
`P(t)=(1;t,0^8)` at `t=a/phi`, strictly before light sqrt(a).
The endpoint square is `-a*d/phi^2<0`. Its initial halfspace is
correctly oriented and r is orthogonal to every retained root.

A finite simple extension of the whole B8 is either an isolated
root of norm1,2,d or the unique B9 long-end extension, whose norm
is2 and whose retained spatial coordinate is1. This is the complete
rank-nine finite classification: no exceptional finite diagram has
the required rank and an old4-bond, and multiple contacts create a
cycle. In the connected row that spectator contributes square1.

All FOUR profiles consequently have the binary equation

    y^2-a*x^2=S,  S=1,2,d,1 respectively.

An actual first normal has x>0, and comparison with the finite upper
gives `a*x/y<=a/phi`. With the conjugate positive-definite bounds:

    x_0^2<=S_0/d_0,
    x_i^2<=S_i/(-a_i), i=1,2,3.                 (2)

No mirror is omitted by a chosen height, norm-square assumption or
already-selected cap. These bounds include the entire closed
interval up to r.

## 3. Complete four-embedding coefficient strips

The identity

    a^-1=(u-1)+s*u

gives `Tr(1/a)=-6` and `Tr(d/a)=-8`. Summing (2), using a>2
and 1/d<2, gives

|S|trace upper|complete even-trace cutoff|
|---|---|---|
|1|`6+1/d+1/a<9`|Tr(x^2)<=8|
|2|`12+2/d+2/a<18`|Tr(x^2)<=16|
|d|`9+d/a<10`|Tr(x^2)<=8|

The first trace cutoff in (1) gives only

    x positive in {1,phi,u,s}.

For S1, phi and s violate the physical bound, whereas1 and u
violate embedding3. For Sd, the physical bound x<=1 leaves1,u;
u violates embedding3, and1 satisfies all four bounds. Thus S1
has no possible x and Sd has exactly x1.

Here is the FULL S2 enumeration. Write x=z+s*w, z,w in Z[phi].
The trace16 cutoff allows only the following families; every sign
is retained before the four bounds are applied:

|w|all possible z|
|---|---|
|0|`0,+/-1,+/-2,+/-phi,+/-u,+/-phi^2,+/-(phi-2)`|
|+/-1|`0,+/-1,+/-phi,+/-u`|
|+/-phi,+/-u|`0,+/-1`|

Indeed U can only be0,4,6,14,16 below17, and 2U can only
be0,8,12. This proves completeness of the families without a
rectangular coefficient enumeration.

For a transparent final finite filter, first retain x>0 and the
physical and embedding3 inequalities in (2):

|family|only survivors of these two tests|
|---|---|
|w0|phi|
|w=+/-1|phi-s, s-u|
|w=+/-phi,+/-u|s*phi-1|

All omitted entries fail one of those tests. The sole close physical
test has the exact identity

    (1+s*u)^2-(2+s)=(s-1)(2phi-3)>0.

For the first surviving entry, embedding1 fails since
`phi^2*(s-u)=s*phi^2-phi>2` (already s>7/5 and phi>8/5 suffice).
For the last three entries, their embedding1 absolute values are
respectively phi+s, a, and s*phi+1, all exceeding2. But
`s-u>3/4`, so the embedding1 bound is below sqrt(8/3)<2.
These elementary comparisons use, for example,
`8/5<phi<13/8` and `7/5<s<10/7`. Every remaining S2 entry is
therefore excluded. The unused embedding2 bound can only shrink
these lists; no selection depends on ignoring it.

This proves the exact profile result

    S1: no x; S2: no x; Sd: x=1 only.

At Sd,x1, the norm equation gives y^2=a+d=phi^2, and y>0
selects y=phi. It is exactly r. No connected or simultaneous
alternative remains. Hence r is uniquely the first SIMPLE wall of
the FULL integral arrangement, and the endpoint has whole actual
finite germ B8+A1(r).

## 4. Distinct architecture and the next boundary

The known old conductor-sixteen quartic host has field
`Q(2cos(pi/8))`, a cyclic quartic field, whereas this field is
biquadratic with discriminant1600. The repository's earlier
Q(sqrt2,sqrt5) finite-diagram/high-dimension field calculations do
not themselves exclude a compact dimension-nine lattice over this
field. A bounded read-only check found no prior literal <-a>+I9
closure. This is not a worldwide novelty or classification claim.

Admissibility itself also rules out scalar descent of this form to
a proper subfield: every indefinite real place of a descended form
would extend to at least two indefinite places here; multiplying by
a nonzero scalar cannot turn an indefinite real form definite.
This does not classify all quartic hosts.

Do not reuse the silver facet identification. The natural binary
orthogonal vector (phi;a) has norm -a*d and its determinant with
r is d, a nonunit. Thus this is not an integral orthogonal splitting
as <-a>+a unit line, nor an identification with a closed silver
lattice. Only a complete first edge is proved. No compact cone,
nonreflectivity, finite-volume conclusion or new example is claimed.

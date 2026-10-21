# A distinct quadratic13 host: complete first two B-type edges

20 October 2026. **Own paper derivation, not executed.** This gives a
literal admissible lattice, its complete primitive norm alphabet, and
two complete actual finite edges. It does not give a compact chamber,
finite volume, reflectivity or nonreflectivity. No scientific worker,
matrix program, Lean or TeX job was run.

Put s=sqrt(13), omega=(1+s)/2, eta=omega+1=(3+s)/2. Then

    omega^2=omega+3,   eta^2=3eta+1,   eta*eta'=-1.

Consider the O=Z[omega] lattice

    L=<-eta> orthogonal-sum I9.

Its two signatures are (9,1) and (10,0). Positive definiteness at the
conjugate place proves anisotropy over K=Q(s). This is an odd unimodular
diagonal host with B-type finite corners, rather than the E8-plus-binary
golden construction. Coordinates below are (time; nine spatial entries).

## 1. Complete norm alphabet by elementary arithmetic

O is norm-Euclidean. Round both coefficients of a+b*omega to integers;
the remaining a,b have absolute value at most1/2. Their norm is
a^2+ab-3b^2=(a+b/2)^2-13b^2/4, so it is at least -13/16 and at most
1/2. Thus its absolute value is less than1. In particular O is a PID.

The positive units are exactly eta^Z. For a unit t with 1<t<eta, its
conjugate has absolute value 1/t<1. If t=a+b*omega, then
0<b=(t-t')/s<(eta+1)/s<2, whence b=1. The inequalities give a=-1 or0;
both omega-1 and omega have norm -3, not +/-1. Thus no such unit exists.
Multiplication by powers of eta reduces every positive unit to this
interval. Since N(eta)=-1, the totally positive units are precisely
eta^(2Z), hence squares of units.

The prime2 is inert: omega^2-omega-3 becomes X^2+X+1 modulo2. For a
primitive root r, unimodularity makes both its coordinate ideal and its
covector ideal O. The integral reflection condition

    2*r_i*(Gr)_j/(r,r) in O   for every i,j

and two Bezout identities imply (r,r) divides2. Its positive value at
both embeddings therefore is a totally positive unit times1 or2. Unit
rescaling gives the complete alphabet

    (r,r)=1 or2.

Conversely every integral vector of either norm is an integral reflector.
No other primitive norm class has been omitted. These arguments apply
in every spatial rank, although only rank9 is used here.

Finite fundamental bonds have only these possibilities: norm2/norm2
label3, norm1/norm2 label4, or orthogonality. Indeed an angle pi/m has
4cos^2(pi/m) in K, which restricts m to2,3,4,6 because the other real
quadratic cyclotomic fields are not K. The norm ratios and the integral
product rule exclude label6, mixed label3, and same-short label3.
This concerns fundamental finite angles, not arbitrary non-simple
mirror angles.

## 2. The full initial B9 corner

At p=(1;0^9), an incident root has time coordinate0. A spatial coordinate
z=a+b*omega satisfies |z|,|z'|<=sqrt(2), so

    |b|=|z-z'|/s<=2sqrt(2)/sqrt(13)<1.

Hence all such coordinates are integers. The incident roots of norm1
and2 are exactly +/-e_i and +/-e_i+/-e_j. Choose the actual chamber germ
with outward simple normals

    alpha_i=-e_i+e_(i+1)  (i=1,...,8),
    alpha_9=-e_9.

It is B9. The integral reflection arrangement is locally finite because
L is admissible, so this finite germ determines an actual chamber C.

Release alpha_1, retaining B8 on coordinates2,...,9. The edge is
p+tau*e_1. A finite one-root extension of this white B8 is either an
isolated norm1 or norm2 root, or a norm2 B9 extension at its long end.
No other finite component contains B8; in particular adding a short
root at the other end would give two double bonds, not a finite diagram.
The complete normal profiles are therefore

    r=(x;z,0^8),       norm1 or2,
    r=(x;z,1,0^7),     norm2.

Their equations are z^2-eta*x^2=R, with R=1,2,1 respectively, and
their collision parameter is tau=eta*x/z, x,z>0.

The integral norm2 mirror

    v=(1;omega,0^8)

has norm omega^2-eta=2 and gives a finite upper tau_v=eta/omega<sqrt(eta).
Its product with alpha_1 is -omega<-2; all other alpha_i are orthogonal.
It is therefore a valid upper before the light endpoint, without first
assuming it simple.

Collision by this upper implies x^2<=R/2. The conjugate equation gives
x'^2<=R/(-eta')=R*eta. Thus for every profile

    0<x<=1,   |x'|<=sqrt(2eta)=sqrt(3+s)<s-1.

The last strict inequality is equivalent to 11>3s, since121>117.
Writing x=a+b*omega, we obtain |b|=|x-x'|/s<1. Hence x is an integer.
R=1 has no positive integer x with x^2<=1/2. R=2 forces x=1, then
z^2=eta+2=omega^2 and positive contact forces z=omega. The bounded
interval therefore has exactly v: it is the actual first simple wall.

The endpoint has apex

    q=(omega;eta,0^8),       q^2=-2eta,

and full finite germ B8+A1(v). The two apex norm classes differ by2,
which is not a square in K; there is no K-isometry carrying the ray of
p onto the ray of q. An immediate same-germ translation is unavailable.

## 3. A second finite upper and its complete six profiles

From q release alpha_2. Retain B7 on coordinates3,...,9 and the isolated
norm2 root v. The edge is q+tau*e_2. The explicit integral norm1 root

    u=(eta;eta+1,eta+1,0^7)

satisfies u^2=1, (u,v)=-1 and u perpendicular to B7. Also
(u,alpha_1)=0, (u,alpha_2)=-(eta+1), with no other old-wall contact.
It gives the finite upper

    tau_u=eta^2/(eta+1)<sqrt(2eta).

It is only an upper at this stage; the root that is simple at this
same contact is determined below.

For any finite new profile write

    r=(X/2)q-(ell/2)v+(Z/2)e_2+t*e_3,
    X=omega*r_time-r_spatial1 in O.

Here t=0 or1, and ell=-(r,v)=0 or1. The complete finite inventory is:

| Diagram after adjoining r | norm r | ell | t | R |
|---|---:|---:|---:|---:|
| isolated short | 1 | 0 | 0 | 4 |
| isolated long | 2 | 0 | 0 | 8 |
| B2 on v | 1 | 1 | 0 | 2 |
| A2 on v | 2 | 1 | 0 | 6 |
| B8 extension | 2 | 0 | 1 | 4 |
| B9 bridge | 2 | 1 | 1 | 2 |

The only connected finite diagrams containing B7 are B8 and B9; to
connect also the retained v the latter must have two simple bonds.
All other possibilities are the displayed isolated A1, A2 or B2 part.
Thus this is a complete list, not an imposed contact sparsity assumption.

In every case

    Z^2-2eta*X^2=R,    tau=2eta*X/Z,
    X^2<=R*eta^2/2,    X'^2<=R*eta/2.

The last two inequalities follow from the proved upper and conjugate
positivity. Binary integrality gives ell congruent omega*X modulo2.
For ell=0 this is X=2Y; for ell=1 both coefficients of X in 1,omega
are odd, since omega^-1=omega+1 modulo2.

The complete finite domains and square tests now fit on one table.
No unit-height loop or floating comparison is involved.

| R, ell | Complete domain | Reduced square test and outcome |
|---|---|---|
| 4,0 | Y=X/2=1 | (Z/2)^2=2eta+1; norm3, impossible |
| 8,0 | Y=1,omega,eta | (Z/2)^2 respectively2eta+2, 2(eta+1)^2, 4(eta+1)^2; first norm12 impossible, second requires sqrt2 outside K, last gives Z=4(eta+1) |
| 6,1 | X=eta-2,eta,eta+2 | (Z/2)^2 respectively eta+1,5eta+3,13eta+5; norms3,29,51, all impossible |
| 2,1 | X=eta | Z=2(eta+1) |

For completeness, derive these domains directly. At R4, Y<=eta/sqrt2
and |Y'|<=sqrt(eta/2). If Y=a+b*omega then b=0 or1 by the broader
bounds Y<5/2, |Y'|<3/2. For b=1 the physical bound forces a<=0,
whereas omega'^2>eta/2 (again11>3s), excluding it. For b=0 only Y=1
remains. At R8, Y<=eta and |Y'|<=sqrt(eta) give b=0 or1; their integer
constants give precisely1,omega,eta. At R2 the same bounds on X and
the odd parity give only X=eta. At R6, X<6, |X'|<7/2 give b=0,1,2;
odd parity gives b=1 and a=-1,1,3, all of which meet the exact bounds.
These elementary strips prove completeness of the table.

The three surviving roots, all at tau_u, are

    u                         (norm1),
    d=u+e_3                   (norm2, the B9 bridge),
    w=(2eta+1;3eta+1,2eta+2,0^7)  (norm2, isolated from B7+v).

Let A=alpha_3+...+alpha_9=-e_3. Then

    u=d+A,          w=2d+v+2A.

Thus u and w are positive combinations of d and the retained outward
normals. The actual simple wall at the simultaneous contact is d,
not u. Its products with v and alpha_3 are -1, with alpha_2 are -eta,
and with all other old normals are0. Consequently d is acute/native
to the entire ten-wall seed; no cap-product filter was used in the
census. The full endpoint is a B9 vertex with outward chain

    v, d, alpha_3,...,alpha_9

and future apex

    H=(3eta;4eta+1,3eta+1,0^7),      H^2=-eta.

All old inequalities hold there. The roots alpha_1 and alpha_2 are
strict; the displayed nine new-corner roots vanish and have positive
definite B9 Gram. Local finiteness and the finite upper justify the
first-contact passage, without assuming C has finite volume.

## 4. A forced finite symmetry, not an infinite-order exclusion

The full old and new B9 germs give a K-isometry G with

    Gp=H,
    Ge_1=(-(eta+1);-2eta,-(eta+1),0^7),
    Ge_2=(-eta;-(eta+1),-(eta+1),0^7),
    Ge_i=e_i  (i=3,...,9).

Indeed it maps alpha_1 to v, alpha_2 to d and fixes alpha_3,...,alpha_9.
The finite-corner Gram and equal apex norms prove full Gram preservation.
Its columns are integral. Direct multiplication using eta^2=3eta+1
gives G^2=I; the inverse is integral and time orientation follows from
Gp=H. Since it permutes the entire integral mirror arrangement and maps
two full actual inward germs of C, it preserves C. Its trace is6,
not an infinite-order trace. Hence this is a genuine chamber involution,
with no nonreflectivity consequence. It also preserves the eleven-wall
seed alpha_1,...,alpha_9,v,d by swapping the first two pairs.

This finite return is the key difference from the recently closed
golden59 host. Further edges, compactness, finite volume and the index
of the full reflection subgroup remain unresolved. Eleven compatible
walls are not a new compact9 example, nor an inference about the
eventual number of walls. No worker is proposed before independent
review of this complete paper reduction.

## 5. Bounded primary-source screen

The repository screen found no literal quadratic13 diagonal-host
exclusion. This is not a complete search of every historical source.
The following retrieved sources do not exclude the displayed host:

* [McLeod's2013 thesis](https://etheses.durham.ac.uk/7743/) states that
  its diagonal classification is over Q. Rational forms -m+I_n do not
  cover the present eta coefficient over Q(sqrt13).
* [Bogachev's2020 research paper](https://nvbogachev.netlify.app/Bogachev-GeomArithmCoxeter-2020.pdf),
  p13, lists the existing unimodular classifications over sqrt2, sqrt5
  and cos(2pi/7). No quadratic13 classification is stated there.
* [Belolipetsky's2016 survey](https://arxiv.org/pdf/1506.03111),
  Proposition5.2 and Table1, pp15-16, retains field discriminant13
  for dimensions8 and9 under its EXTRA congruence assumption. This is
  neither an exclusion nor a proof that the present form is reflective.

The original full text of Bugaenko's1990 unimodular paper was not
retrieved in this screen; its title alone is not evidence for a
quadratic13 classification. Higher-dimensional field lists have not
been imported at dimension9. The own algebra above, rather than an
assertion of novelty or an incomplete literature search, is the reason
to retain this as a concrete next construction route.

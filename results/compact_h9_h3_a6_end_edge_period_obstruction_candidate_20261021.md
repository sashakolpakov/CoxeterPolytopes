# The literal golden H3+A6 host: a forced end-edge period obstruction

21 October 2026. **Paper candidate, pending root and independent full review.**
No enumeration, Pell search, scientific import, checker, numerical job,
Lean, TeX or canonical-file edit ran. The earlier construction frontier
is preserved. This argument concerns only

`L=<-phi> orthogonal-sum H3 orthogonal-sum A6`, `O=Z[phi]`, `phi^2=phi+1`.

**Claim.** The full integral reflection chamber is noncompact, hence its
full reflection group is noncofinite by the explicit anisotropic
arithmetic-lattice input. The mechanism is a forced integral chamber
period, not a proof that the selected edge has no finite endpoint. No
similarity-class, commensurability-class or field-wide conclusion follows.

## 1. Bound starting corner and hypotheses

Use the native H3 Gram `[[2,-phi,0],[-phi,2,-1],[0,-1,2]]`, and A6
simple-root basis `E1,...,E6` with diagonal2 and adjacent product-1.
The initial outward normals are the negatives of these bases. In
time/H3/A6 coordinates put

`w=(1+2phi,2+2phi,1+phi)`, `P=(1;w/2;0)`, `r=(phi;w;0)`.

Then `H3*w=(0,1,0)`, `w^2=2phi^2`, `P^2=-1/(2phi)`, and `r^2=1`.
The previously derived first middle-H3 edge gives P as an actual full
corner with outward germ `a=-h1`, `b=-h3`, `r`, and `-E1,...,-E6`.
The first three are mutually orthogonal of norms2,2,1 and perpendicular
to P. This derivation is source-bound in the
[reentry frontier](compact_h9_construction_reentry_frontier_20261020.md):
spatial Cauchy bounds every intervening wall by height `phi^3<30`, and
the inherited complete Vinberg census has only walls9/10 at that height,
of which only wall10 meets the segment. Its positive-height spatial
exhaustion is source-reviewed, not independently rerun in this task.

The complete primitive norm pool, after unit rescaling, is
`q in {1,2,4,7,14,28}`. The bound follows from `14G^-1` integral,
primitivity and `q|28`; the golden PID and totally-positive-unit-square
facts give these representatives. These are inherited exact inputs.

Suppose the full chamber C is compact. Release `-E1` at P, retaining
`A5+3A1`. Its next actual endpoint is finite. Let m be its arriving
outward simple normal, primitive of norm q in that pool. Its products
with the eight retained normals are finite native Coxeter products;
the full replacement Gram is positive definite. Let Omega be the A6
end weight, `Omega=(6,5,4,3,2,1)/7` in the E basis, `Omega^2=6/7`.
The outgoing ray is `P+t*Omega`, with `t>0`.

## 2. Coordinates and finite contact alphabet

Write `(m,a)=-l`, `(m,b)=-n`, `(m,r)=-p`, with `l,n,p>=0`.
The orthogonal decomposition of its time/H3 part is

`m_H = T P -(l/2)a -(n/2)b -p r`, `T=m_time+phi*p in O`.

Consequently

`m_H^2=-T^2/(2phi)+(l^2+n^2)/2+p^2`.                 (1)

The finite-angle/field/integrality test gives this complete table:

| New q | Product magnitude with a retained norm2 root | With retained norm1 root r |
|---|---|---|
| 1 | 0 or1 (label4) | 0 |
| 2 | 0,1,phi (labels3,5) | 0 or1 (label4) |
| 4 | 0 or2 (label4) | 0 |
| 7,14,28 | 0 | 0 |

Indeed the native angle possibilities over Q(phi) have labels
2,3,4,5,6,10; taking the required field square roots excludes label6
and label10 here, and excludes every nonorthogonal prime7 contact.
For q4 against norm1, the otherwise possible products1,phi violate
`2(m,L)/q subset O`. These are actual finite facet angles, not an
alphabet assertion about arbitrary intersecting arrangement mirrors.
There is at most one contact on retained A5, since two would create a
cycle in its positive finite replacement.

## 3. No endpoint orthogonal to A5

If m has no A5 contact, its A6 part is `7z*Omega`, z in O, and (1)
gives

`T^2-84phi*z^2=-phi*R`, `R=2q-l^2-n^2-2p^2`.        (2)

At7, phi is a nonsquare in F49: its norm is -1, a nonsquare in F7.
Multiplication by any nonzero rational residue preserves this status.
For q2 the positive finite A1-extension list is: isolated; single3;
single5; two3; 3+5; and 3 together with the norm1 label4.
The norm1-only label4 row is also included. Their respective R are
`4,3,3-phi,2,2-phi,1,2`. The rational R rows are excluded modulo7;
`R=2-phi` is excluded because `N(-phi*R)=-1` modulo7 is nonsquare.
For `R=3-phi`, equation(2) has right side `1-2phi=-sqrt5`.
Modulo the prime `pi=sqrt5`, one has `phi=3` and `84phi=2` in F5.
Thus `T^2-2z^2=0` forces both T,z divisible by pi. The left side of(2)
then has valuation at least2, but its right side has valuation1.

The q1/q4 single-label4 rows have R=1 and R=4, respectively, and are likewise
excluded modulo7. Two label4 contacts already give an affine rank3
principal, so there are no further q1/q4 rows. For prime7 norms all
contacts vanish. Integrality then gives `T=2x`, and

`q=42z^2-2phi^-1*x^2`, `q|2x`, `q|14z`.

Odd q7 is impossible modulo2. For q28 the divisibilities make x and z
even, contradicting primitivity. For q14 put x=7X and divide by14:
`3z^2-7phi^-1*X^2=1`. Modulo3 this gives `X^2=-phi` in F9, impossible:
`N_F9/F3(-phi)=-1`. The isolated q1/q4 rows are already excluded by
the same rational-R test. This exhausts the no-A5-contact alternative.

## 4. Contacts on A5: a complete short filter

Number its five nodes by j=1,...,5, corresponding to old E_(j+1).
Let d be the product magnitude at that unique contact. A6 integrality
is exactly

`m_A6=7z*Omega+d*s_j`,
`s_j=(-j,1-j,...,-1,0,...,0)` (j initial nonzero entries),

since its free E1 covector plus `(j+1)d` is divisible by7. Hence

`m_A6^2=42z^2-14jdz+j(j+1)d^2`.

Set `Y=6z-jd` and `sigma=l^2+n^2+2p^2`. Combining with(1) gives

`3T^2-7phi*Y^2=phi*[3sigma+j(6-j)d^2-6q]`.           (3)

The full finite-classification list, followed by its rational bracket
on the right of(3), is:

| Finite profile | q,d | j | sigma | Bracket |
|---|---|---|---:|---|
| B6, no A1 contact | 1,1 or4,2 | 1,5 | 0 | -1 or-4 |
| A6/D6/E6, no A1 contact | 2,1 | 1,2,3,4,5 | 0 | -7,-4,-3,-4,-7 |
| A7/E7, one norm2 A1 contact3 | 2,1 | 1,2,4,5 | 1 | -4,-1,-1,-4 |
| B7 (norm1 contact4), or D8 (two norm2 contacts3) | 2,1 | 1,5 | 2 | -1 |

This is exhaustive: a finite component containing A5 has rank6--9,
so only A/B/D/E6/E7/E8 can occur. A weighted bond must be terminal
in a B path. A new contact to an additional A1 makes the A5 bond
internal, excluding q1/q4. For q2, one added norm2 A1 gives A7 or
E7; its central attachment gives affine E6. Two norm2 A1s give D8
only at an A5 end; an internal attachment already contains the
affine two-fork. The norm1 A1 gives only terminal4 in B7. More
contacts give a weighted branch or degree4 and are not finite.
In particular an E8 diagram is not omitted: the present attachment
shapes at rank8 are the displayed D8 or a nonpositive two-fork.

Since phi times a nonzero rational is nonsquare modulo7, equation(3)
eliminates every row except q2, d1, sigma0, j=1 or5. These survivors
also have `7|T`. They are orthogonal to a,b,r, so `T=2x` by original
H3 integrality, and `7|x`.

For j1 their equation is `x^2=7phi*z*(3z-1)`. The outgoing direction
forces z>1/3; conjugate positivity forces `0<z'<1/3`. Thus z is
totally positive. Since `(z,3z-1)=1`, unique ideal factorization
gives `z=epsilon*s^2` or `7epsilon*s^2`; every non7 prime has even
valuation. Total positivity of z makes epsilon totally positive,
hence a unit square. Therefore z=s^2 or7s^2 after absorbing it.
The equation makes s divide x. Cancel s^2 and reduce modulo3:
`(x/s)^2=-phi` in F9, since 7 and49 are1 modulo3. This is impossible.
Cancellation also handles arbitrary 3-adic content; there is no
unproved restriction to solutions coprime to3.

Only j5 remains. Put u=z-1. Its equation and orientation are

`x^2=7phi*u*(3u+1)`, `x>0`, `u>0`, `-1/3<u'<0`, `7|x`.          (4)

Here the far-end orientation is essential: replacing z by1-z changes
the physical sign. The previous totally-positive-unit argument does
NOT exclude this row. It is retained and gives the period below.

## 5. The surviving endpoint forces an integral chamber period

Let `H=E1+...+E6` and `Y=1+6u`. The arriving normal is exactly

`m=2xP+H+7u*Omega`.

Define B to fix a,b,r and every vector in A5=span(E2,...,E6), and set

`B(P)=Y P+(x/phi)Omega`,
`B(Omega)=Y Omega+(12x/7)P`.

Equation(4) gives `Y^2-12x^2/(7phi)=1`, so B is a future Lorentz
isometry. Both B and its inverse preserve L. Explicitly

`B(E1)-E1=7u*Omega+2xP` is integral;

for every original time/H3 lattice vector v its P coefficient is
`-2phi*(v,P) in O`, and

`B(v)-v=[-2phi*(v,P)]*[6uP+(x/phi)Omega]`

is integral because `6uP` is integral and `7|x`. The inverse has the
same formulas with x negated. There is no rational-isometry-to-integral
group substitution here.

Let R be the ordinary integral A6 Coxeter rotation, fixing time/H3,
with `R(E_i)=E_(i+1)` for i<6 and `R(E6)=-H`. Then S=BR fixes a,b,r
and maps the complete old outward germ as follows:

`-E1,...,-E6  -->  -E2,...,-E6,m`.

Also `S(P)=B(P)` is exactly the finite edge endpoint ray. Since m is
the actual first arriving wall, both complete oriented germs belong
to the same full chamber C. Local finiteness of the full integral
mirror arrangement therefore gives `S(C)=C`.

For transparency the trace is independent of unknown unit size.
On A6, `tr(R)=-1` and `(Omega,R(Omega))/Omega^2=-1/6`.
On the time/H3 block R is identity. The boost's off-diagonal terms
have zero contribution to `tr(BR)`, so

`tr(S)=3+(Y-1)*(1-1/6)=3+5u`.

Discreteness here is elementary: over the other embedding every integral
isometry lies in the compact orthogonal group of the positive definite
conjugate Gram. A bounded set of physical matrices therefore bounds all
entries at both embeddings. Only finitely many O-valued matrices satisfy
these bounds, proving discreteness of the physical `O^+(L)`.

Write u=A+Bphi in O. The two signs in(4) imply B>=1 and
`u=u'+B*sqrt5>sqrt5-1/3>5/3`. Thus `tr(S)>34/3>10` and S has
infinite order. An infinite subgroup of the discrete integral
isometry group cannot stabilize a compact full-dimensional chamber.
This contradicts the assumed compactness of C.

## 6. Scope and exact bindings

This proves noncompactness if the inherited actual first-corner input
and the displayed finite/profile arguments are accepted. The additional
standard arithmetic input is explicit: admissibility and conjugate
positive definiteness make `O^+(L)` a uniform arithmetic lattice.
A cofinite reflection subgroup would have finite index and a compact
chamber, contradicting the argument. Thus the literal full reflection
group is noncofinite. No alternate lattice in this K-space is excluded.
No exact first far-end mirror, fundamental Pell unit or fresh replay
has been computed or claimed. Earlier literal-host closures are untouched.

| Inherited input | SHA256 |
|---|---|
| [H3+A6 source report](fresh_h3_blocks_h9_report_20261013.md) | `4545875776e23edf97c0a4a690d8eb0589bf76d3600fe6b7f5daf65081910b19` |
| [45-wall H30 certificate](fresh_h3_blocks_h9_a6_phi_h30_20261013.json) | `ee0370bc5306d01d0f6954ae0c4759e1ead4194c34190e189ec06d88428b8d5c` |
| [Independent audit](fresh_h3_blocks_h9_independent_audit_20261013.json) | `d4acf35f13083ee44419c9e0bc7edd6cda2503d53ca93a7fbc07a38a5ff83863` |
| [Producer source](../scripts/fresh_h3_blocks_h9_a6_search.py) | `fcfe23a3beadf58242df05f0c91c0f1c93462b7247bedf2606f717609b4898b7` |
| [Independent verifier](../scripts/fresh_h3_blocks_h9_independent_verify.py) | `06124f9147cb6c70c43928033bc8c57227d80865f1deb3711098976fe2ddc8fd` |
| [Reentry frontier](compact_h9_construction_reentry_frontier_20261020.md) | `668a075c6001f0860d05fb848a8c2357de59ac133d286f8dbbf599c603a49df6` |

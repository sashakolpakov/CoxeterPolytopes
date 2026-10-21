# A positive-part opposition bound: candidate exclusion of mixed dimension twelve

19 October 2026. **New paper proof candidate, pending independent review
and exact coefficient replay.** No scientific process, optimizer,
enumeration, matrix calculation, Lean or TeX job ran in this lane.
Earlier reviewed sources are unchanged. No dimension table or manuscript
is changed by this preparation note.

The proposed universal upper bound at dimension12 is

    2B <= (1420561/12500) f = 113.64488 f.             (1)

Here f is the number of FINITE vertices and B the number of canonical
P3 three-faces, not polygonal bad triangles. The established polygon/
three-face coupling lower bound is

    2B > (184460/1619) f.                            (2)

Since

    184460*12500 - 1420561*1619 = 5861741 > 0,

the two bounds contradict each other whenever f>0. Conditional on
acceptance of the actual-object and moment interfaces below and the
new polynomial proof, this excludes mixed dimension12. The previously
established general ceiling12 then gives a mixed ceiling11. This is
NOT yet a promotion or a claim to exclude dimension11.

## 1. Actual opposition objects and the full formal row

Use the reviewed actual P3 opposition geometry from
`right_angled_mixed_opposition_moments_research_20261018.md`.
At a finite vertex q of an actual finite-volume right-angled12-polytope,
let q^2=-1 and n_1,...,n_12 be its orthonormal outward incident normals.
A graph triangle S in its bad-triangle graph determines an actual
canonical three-face P3 and its other finite vertex

    q_S=2q-sum_(i in S) n_i.

Distinct triples at q give distinct opposite vertices. Each P3 has
exactly two finite vertices, so the actual finite opposition graph
has B undirected edges and degree sum2B. It is finite. Its invariant

    h=3q-sum_i n_i,     h^2=3,     h.q=-3

is the same at both endpoints of every opposition edge, including
vertices with incomplete bad-triangle graphs. On EACH finite opposition
component choose one fixed nonzero future null vector Z perpendicular
to h. It exists because h is spacelike. For isolated vertices make
the same harmless componentwise choice. Put

    E(q)=-Z.q>0,   p_i=-Z.n_i/E(q),   g(q)=E(q)^(-2).

The equations Z.h=0 and Z^2=0 give sum p_i=3 and sum p_i^2=1.
There is no sign assumption on p_i. Put x_i=p_i-1/4. Then
sum x_i=0 and sum x_i^2=1/4.

For ALL220 formal triples S, including those that are not actual
oppositions, q_S=2q-sum_S n_i is future unit timelike: its square
is -1 and q.q_S=-2. Consequently its formal energy ratio is

    r_S=E(q_S)/E(q)=2-sum_S p_i=5/4+d_S,
    d_S=-sum_S x_i,     -3/4<=d_S<=3/4.             (3)

The last bound is zero-sum Cauchy: the three-subset indicator minus
its mean has squared norm9/4, and the x squared norm is1/4.
Thus 1/2<=r_S<=2 at EVERY finite vertex. A formal q_S is used only
to define these220 positive row numbers; it is not presumed feasible
as a vertex of the original chamber.

## 2. Pointwise moments and their signs

Uniform averages over the220 formal triples at ONE vertex satisfy

    average d=0,
    average d^2=9/176,
    average d^3=-(27/220)*sum_i x_i^3,
    average d^4=3/440-(3/220)*sum_i x_i^4 <=3/440.

These are the reviewed exact subset-moment identities, valid for
every vector satisfying the two displayed x equations. Moreover
abs(sum x_i^3)<=sum |x_i|^3<=(sum x_i^2)^(3/2)=1/8, so

    -average d^3 <=27/1760.                       (4)

No average over vertices, geometric height maximization, restricted
three-face spectrum, finite-edge hypothesis or graph census occurs
in these identities. The cubic coefficient of the polynomial below
is NEGATIVE, so (4) is the needed upper bound with the correct sign.
The quadratic and quartic coefficients are positive, so their stated
upper bounds are used in the ordinary direction. The linear average
is exactly zero.

## 3. The proposed positive-part polynomial majorant

For -3/4<=d<=3/4 define

    P(d)=19/125 -(979/1000)d +(768/625)d^2
                   -(971/3125)d^3 +(53186/9375)d^4.

The new scalar claim is the complete interval inequality

    P(d) >= max(0, (5/4+d)^(-2)-1/2).              (5)

Both parts are necessary: nonnegativity permits inclusion of absent
formal triples later, while the reciprocal part controls every
actual opposition. Here is a rational polynomial proof on the full
domain, rather than point sampling.

### 3.1. Nonnegativity

For d<=0 every term of P has nonnegative sign and its constant
is positive. For d>=0,

    P''(d)-12/5
       =36/625 -(5826/3125)d +(212744/3125)d^2
       >=36/625-2d+60d^2 >0.

The last quadratic has minimum36/625-1/60>0. At d0=1/4,

    P(d0)=1627/1200000 >13/10000,
    P'(d0)=-10243/150000,       |P'(d0)|<7/100.

Taylor's convex lower bound with P''>=12/5 therefore gives,
for every d>=0,

    P(d)>=13/10000-(7/100)^2/(24/5)
          =67/240000>0.

This proves the first half of(5).

### 3.2. Reciprocal majorization at d>=0

Write r=5/4+d and mu=5/4. The exact reciprocal-square identity is

    r^(-2)=mu^(-2)-2d/mu^3+3d^2/mu^4-4d^3/mu^5
                     +d^4*(4r+mu)/(mu^5*r^2).

On r in[1/2,2], the factor(4r+mu)/r^2 is at most13, its value
at r=1/2. Thus the earlier Taylor majorant for r^(-2)-1/2 is

    P0(d)=7/50 -(128/125)d +(768/625)d^2
                  -(4096/3125)d^3 +(13312/3125)d^4.

The new polynomial differs by

    P-P0 =3/250 +(9/200)d +d^3 +(106/75)d^4.

For d>=0 this is nonnegative, proving reciprocal majorization
on that half of the interval.

### 3.3. Reciprocal majorization at d<=0: one explicit degree-six polynomial

Set t=-d in[0,3/4] and r=5/4-t. Direct multiplication gives

    N(t)=r^2*(P(-t)-r^(-2)+1/2)
        =3/160 -(321/3200)t +(249/2000)t^2 -(643/400)t^3
             +(27949/3000)t^4 -(130052/9375)t^5
             +(53186/9375)t^6.                    (6)

It remains to prove N>=0 on this one compact rational interval.
The following elementary rational certificate supplies that proof.

For a polynomial Q of degree m on[a,b], write
Q(a+(b-a)u)=sum_k c_k u^k and set

    beta_i=sum_(k<=i) c_k*C(i,k)/C(m,k),  0<=i<=m.

Then Q(a+(b-a)u)=sum_i beta_i*C(m,i)u^i(1-u)^(m-i).
The Bernstein basis is nonnegative and sums to one on[0,1].
Thus bounding each beta_i bounds Q on the WHOLE interval.
The displayed coarse rational bounds can each be checked directly
from(6) and this formula; no roots, numerical approximation or
unrecorded interval subdivision is needed.

|Polynomial and interval|Bounds on all its successive Bernstein coefficients|
|---|---|
|N' on[0,1/5], degree5|less than respectively -1/10,-9/100,-99/1000,-98/1000,-79/1000,-9/200|
|N' on[3/10,1/2], degree5|greater than respectively 3/50,11/100,17/100,21/100,23/100,1/5|
|N'' on[1/2,3/4], degree4|less than respectively -3/5,-19/10,-37/10,-29/5,-36/5|

All comparisons are strict; they were hand-expanded as rational
arithmetic during this preparation. Exact replay remains a separate
requested verification layer, not an executed fact here.

The first row makes N decreasing on[0,1/5]. Direct substitution
gives N(1/5)>1/1000, hence N>0 there.

For the central interval[1/5,3/10], the quadratic N'''' is convex
and satisfies

    N''''(1/5)<-27,       N''''(3/10)<-91,
    N'''(3/10)>9/10,     N''(1/5)>4/5.

It follows successively that N''''<0, N'''>0, and N''>1/2
throughout that interval. At its middle point1/4,

    N(1/4)=679/1200000 >1/2000,
    N'(1/4)=3061/600000,       |N'(1/4)|<1/100.

The convex quadratic lower bound consequently gives
N>=1/2000-(1/100)^2>0 throughout the central interval.

The second Bernstein row makes N increasing on[3/10,1/2], so
positivity propagates from the already covered central endpoint.
The third row makes N concave on[1/2,3/4]. Both endpoints are
positive, the first by the preceding interval and the second by

    N(3/4)=57/64000>0.

A concave function is at least its chord, so N>0 there too.
These four intervals cover[0,3/4], proving(6) nonnegative.
As r>0, this proves the remaining half of(5).

## 4. Row deletion and ordinary finite-vertex counting

Let t(q) be the degree of q in the actual P3 opposition graph.
For each actual neighbor with ratio r, (5) gives

    r^(-2)-1 <= P(r-5/4)-1/2.

Sum over actual triples and then include all absent formal triples
using P>=0. At every finite vertex this gives

    sum_actual (g(q_S)/g(q)-1)
          <=sum_all220 P(d_S) -t(q)/2.             (7)

The pointwise moments of Section2 and the signs of the coefficients
give the uniform bound

    average_all220 P
      <=19/125 +(768/625)*(9/176)
              +(971/3125)*(27/1760)
              +(53186/9375)*(3/440)
       =1420561/5500000.                          (8)

Sum(7) over ALL finite vertices, allowing the componentwise choices
of Z specified in Section1. On each actual undirected edge the
left side contributes r^(-2)+r^2-2>=0, because the reverse ratio
is reciprocal. There is no use of an unweighted height-maximum
bound, no assumption that g is globally constant, and no replacement
of a maximum by an ordinary average. It follows directly that

    0 <=220*(1420561/5500000)*f -(1/2)*sum_q t(q),
    sum_q t(q)=2B <=(1420561/12500)*f.

This is exactly(1). The constant1/2 in(5) is an algebraic
deletion allowance, not a new move weight. The actual move set
is the old undirected P3 graph; no quad reversals or carrier
multiplicities are added or changed.

## 5. The opposing strict lower bound and exact scope

The reviewed source
`right_angled_mixed_polygon_three_face_coupling_research_20261018.md`
proves at dimension12

    101B-5130a2+69080f+22049280c >=0.

Its established strict weighted row a2>15f+4246c gives

    101B>7870f-267300c.

Independently the canonical P3 sum gap and that same weighted row
give B>50f+880c. Multiply this last inequality by267300 and
the preceding one by880. Their cusp terms cancel, yielding

    (101*880+267300)B >(7870*880+50*267300)f,
    356180 B >20290600 f,
    2B >(184460/1619) f.

Both multipliers are positive; strictness is retained. Here c
counts ambient ideal vertices, while f counts only finite ones.
The same B is the canonical P3 count used in the actual opposition
graph. Thus no two-dimensional bad-triangle count T, total-vertex
count f+c or noncanonical three-face count is substituted into
the comparison.

A hypothetical mixed12-polytope has f>0 by definition. Purely
compact/ideal cases already have separate accepted ceilings below12;
this candidate does not reprove those. If all stated reviewed
geometric/count interfaces and the new exact polynomial certificate
are accepted, (1)--(2) exclude dimension12. Dimension11 remains
outside this numerical contradiction; no example, novelty claim
or formal geometry verification is supplied here.

## 6. Dependency and verification boundary

The actual-opposition/moment source, reciprocal-square identity source,
polygon coupling source and weighted-triangle source are the unchanged
reviewed18October packets with their matching stems under results/.
Their source hashes should be bound by the companion exact certificate.
The new finite arithmetic obligation is only (6), its stated Bernstein
and central-interval bounds, the coefficients in(8), and the final
cross-multiplication. It requires no optimizer, geometry enumeration,
giant JSON input, matrix calculation or historical census. A standard-
library Fraction checker can verify these identities in a tiny bounded
root-scheduled batch after independent source/mathematical review.

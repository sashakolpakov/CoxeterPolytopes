# Complete five-dimensional corners: a sharp twenty-facet gap and a corner-density bound

Date: 2026-10-18. Mathematical research packet; no Lean verification and no change
to the universal dimension bounds. This packet extends, without modifying, the
reviewed weighted-triangle and low-facet-rigidity packets of the same date.

## Statements and hypotheses

Let F be a finite-volume right-angled hyperbolic five-polytope. At a finite
vertex q, label its five orthogonal outward unit normals n_1,...,n_5, with
q.q=-1. A *bad triangle* has one finite and two ideal vertices. The corner
graph joins i,j exactly when the two-dimensional face through those edge
directions is a bad triangle. A *complete corner* means that this graph is K5.

The following statements hold.

1. A complete corner contains a canonical P5 seed: the convex hull of the 26
   vertices in the sixteen-wall model below lies in F.
2. If F has a complete corner, then either F is canonical P5, or F has at least
   twenty facets. The bound twenty is attained by the reflected double of P5.
3. If F is not P5, and L and V denote its number of complete finite corners and
   its number of finite vertices, then **11L <= V**. The same twenty-facet
   double attains equality, with L=2 and V=22.

In particular, a K5 corner alone does not force P5. These are actual local
geometric restrictions; they do not exclude mixed dimension twelve.

## The normal gap and the initial fifteen walls

Use the Lorentz form with negative time coordinate. Distinct outward unit
facet normals of an acute-angled polytope have product at most zero. For a
right-angled polytope the product is either zero or at most -1. For
self-containment, if two normals u,v had g=u.v in (-1,0), their positive
definite Gram matrix G would have a nonnegative inverse. For an interior
future timelike point x, put b=(x.u,x.v), whose entries are strictly negative.
Its orthogonal projection p=x-[u v]G^{-1}b is x plus positive multiples of
u,v. It is timelike because p.p=x.x-b^T G^{-1}b<0, is future because p.x<0,
and is strict for every other wall by acuteness. Thus the two facets meet
along an actual ridge at a nonright angle, a contradiction. This is also the
normal-gap fact used in the reviewed rigidity packet.

Every bad triangle ij has opposite outward unit normal

    m_ij = q - n_i - n_j.

Thus a complete corner supplies the five n_i and the ten distinct m_ij.
Three-subset corners give actual P3 faces with opposite finite vertex

    q_S = 2q - sum_{i in S} n_i,       |S|=3.

Four-subset corners give actual P4 faces, hence actual cusps

    v_i = 2q - sum_{j != i} n_j.

These facts, including feasibility against every remaining wall, were proved
in the reviewed weighted-triangle packet. They also follow directly from
writing an additional wall as r=a q-sum c_i n_i. Here a>0, c_i>=0 and
c_i+c_j<=a, so all the relevant finite points have nonpositive wall product.
For the null vertices, a positive sum of known outward normals has
nonpositive product with every other outward normal.

## A canonical seed is always contained

In Klein coordinates x=q-sum y_i n_i, the selected sixteen inequalities are

    y_i >= 0,        y_i+y_j <= 1,        sum y_i <= 2.

The last wall has normal rho=2q-sum n_i; it need not yet be a facet of F.
The vertices of this canonical P5 model are:

| Coordinates | Count | Type |
| --- | ---: | --- |
| origin | 1 | finite |
| one coordinate 1 | 5 | ideal |
| three coordinates 1/2 | 10 | finite |
| four coordinates 1/2 | 5 | ideal |
| one coordinate 2/3 and the other four 1/3 | 5 | finite |

The first four classes are feasible against every wall of F by the preceding
corner constructions. For an additional wall r=a q-sum c_i n_i, the five
remaining finite vectors are z_i=3q-2n_i-sum_{j!=i}n_j. It suffices to show

    2c_i + sum_{j!=i} c_j <= 3a.

If c_i>=a/2, sum the four inequalities c_i+c_j<=a to get a left side at most
4a-2c_i<=3a. If c_i<=a/2, partition the other four indices into two pairs;
their sum is at most 2a, giving the same conclusion. The original five n_i
also contain these points. Thus every vertex of the seed is in F, and
convexity proves seed containment. This proof does not assume that rho is
an actual wall.

## Cusp widths prove the twenty-facet gap

At v_i, the six walls m_ab with a,b different from i form three opposite
pairs, each pair summing to v_i. Along with n_i these are seven of the eight
facets of the rectangular four-dimensional cusp box. Its unique remaining
wall, opposite n_i, has the form

    r_i = lambda_i v_i - n_i,          lambda_i > 0.

The coefficient of n_i is -1 because its image in the Euclidean cusp link
is the opposite unit normal; the remaining difference is a multiple of the
null vector v_i. Evaluation at q fixes its positive sign. For j different
from i,

    r_i . m_ij = 1-lambda_i.

The normal gap therefore gives lambda_i=1 or lambda_i>=2. For i different
from j, direct expansion gives

    r_i . r_j = lambda_i+lambda_j-lambda_i lambda_j
              = 1-(lambda_i-1)(lambda_j-1).

If any lambda_i=1, the last product is 1. Distinct acute outward normals
cannot have product 1, so all these r_i coincide. Comparing coefficients
then forces every lambda_j=1, and the common wall is rho. Since the seed is
contained in F, while the sixteen actual walls imply F is contained in the
seed, F=P5.

Otherwise every lambda_i>=2. The five r_i are pairwise distinct: equality of
two would first force equal lambdas by their q coefficients, and comparison
of the two exceptional n coefficients would force that common lambda to be
1. None is one of the original fifteen walls (its q coefficient is at least
4). Hence F has at least twenty facets. A further genuine width restriction
is that (lambda_i-1)(lambda_j-1) is 1 or at least 2; the value 1 requires
both lambdas to equal 2.

## Incoming opposite corners give the sharp finite-vertex density

Assume F is not P5. Every complete corner q has ten distinct P3 opposite
finite vertices q_S, indexed by triples S. At q_S the five incident normals
are the three m_ij with i,j in S and the two untouched normals n_l, l not
in S. They form an orthonormal basis. If the untouched pair were a bad
triangle at q_S, its opposite wall would be

    q_S - sum_{l not in S} n_l = rho.

The preceding seed argument would make F=P5. Thus q_S is not complete, and
the complementary pair is a missing edge of its corner graph.

Conversely, suppose a finite vertex p receives a P3 opposite from some
complete corner. Label the corresponding triangle at p by S={1,2,3}, so the
complete opposite corner is p_S=2p-n_1-n_2-n_3. Completeness at p_S forces,
for each k in S and l not in S, the wall

    p_S - m_ij - n_l = p-n_k-n_l,      {i,j,k}=S.

These force all six cross edges in the corner graph at p. To check that an
actual wall of this form indeed gives the bad triangle, its two required
null endpoints are

    p-n_k = (p-n_k-n_l)+n_l,
    p-n_l = (p-n_k-n_l)+n_k.

They are future null, feasible against every other wall by acuteness, and
span the two ideal vertices of the selected three-wall two-face. The
triangle S already supplies its three internal edges. The complementary
pair is missing by the previous paragraph. Thus the recipient graph is
exactly K5 with that single complementary edge removed. Its unique missing
edge determines S, and S determines the opposite point. A recipient
therefore receives at most one complete corner.

Count directed opposite incidences. Their ten sources at each complete
corner land at distinct noncomplete finite vertices, each receiving at most
one. Hence 10L<=V-L, or 11L<=V. A weaker graph-only observation gives at most
four candidate incoming triangles: their complementary pairs are edges of
the complement graph that are vertex covers, and these pairs intersect
pairwise, so form a star or triangle. The geometric cross-edge condition is
what improves four to one.

## Sharp model and exact finite checks

Reflect P5 across its wall rho. The other ten m_ij are orthogonal to rho and
remain the same exterior walls; the reflected five n_i are

    n_i+2rho = 2v_i-n_i.

The double has precisely these twenty walls: the original fifteen and the
five reflected n_i. All distinct normal products are 0,-1,-2. Gluing across
a right-angled facet gives a convex right-angled double; the following
exact enumeration also checks its entire Klein body. It has 37 vertices,
22 finite and 15 ideal, and exactly two complete corners. Thus the
twenty-facet gap and the density fraction 1/11 are both sharp.

The companion checker uses rational arithmetic, with no optimizer:

* It reuses the previously checked 26-vertex seed.
* It solves all 15,504 five-wall choices for the double, checks the complete
  vertex set against the union of the two reflected seeds, checks the
  ball condition and all twenty walls, and checks the finite corner graphs.
* It checks all 1,024 labelled five-vertex graphs. The graph-only candidate
  maximum is four; the strengthened geometric maximum is one, attained
  precisely by the ten labelled graphs K5 minus one edge.

The finite checks support the models and elementary graph arithmetic. The
universal cusp-box, seed-containment and opposite-incidence arguments above
remain written mathematical proofs, not machine-checked universal claims.

## Global consequences and a surviving aggregate assignment

For an n-dimensional polytope, let D count canonical P5 faces, and let J
count noncanonical five-faces containing at least one complete corner. The
twenty-facet gap, compared with the generic seventeen-facet gap, gives

    D >= 17 a_5 - (n-4) a_4 + 3J.

Let K be the total number of K5 subgraphs over all finite corner graphs.
Each belongs to a unique five-face. A canonical P5 has sixteen complete
corners; every other five-face has complete-corner count at most one
eleventh of its finite vertices. Since the total finite-vertex/five-face
incidence count is binomial(n,5) f,

    11K <= binomial(n,5) f + 160D.

This excludes the earlier frozen aggregate witness with D=3000. It does not
exclude its repair obtained by changing only D to 4000. The checker replays
every assertion of that frozen packet, then checks every D-dependent
inequality with the new D and the new bound. All pass. In that repair

    K=64000=16D=8C,       J=0,

so all complete corners have been assigned to canonical five-faces. The new
bound has strictly positive slack

    (binomial(12,5) f + 160D)/11 - K
       = 78172078997921 / 68750000000.

Consequently the new genuine geometric restrictions still leave this
aggregate route feasible. They sharpen the remaining question to geometric
compatibility and propagation of canonical five-faces at the dense finite
corners. The repair is an incidence assignment, not a realizable polytope.
No universal bound on the number of complete corners equal to two is
claimed; two is the count in the exhibited double only.

## Reproduction and scope

Run `python3 scripts/check_right_angled_mixed_complete_five_corner_research_20261018.py`
to replay the saved JSON, or use `--write` to regenerate it before freezing.
The checker is tiny exact arithmetic, not a proof process. It runs no Lean,
linear optimizer, background worker or large numerical search. This packet
does not modify any frozen earlier evidence. It claims neither a new
dimension ceiling nor a new example beyond the explicit reflected double.

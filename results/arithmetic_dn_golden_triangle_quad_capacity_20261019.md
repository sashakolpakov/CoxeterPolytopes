# A shared capacity for golden triangle and quadrilateral exceptions

19 October 2026. **New written capacity theorem for actual wall
profiles, awaiting independent review.** No computation or build was
run. No full refund, global example or dimension improvement is claimed.

## 1. The two profiles

Let a compact hyperbolic Coxeter polytope have a whole D_n finite
vertex, with fork a, whiskers b,e and long arm a-l1-...-lL.
Here n is 19 or 20, N=n-2 and L=N-1. Unit outward roots are used.
Put

    x=(1+sqrt(5))/4,     4*x^2=2*x+1.

There are two different kinds of potentially deficient face at the
fork/whisker angle (a,b), under a proposal deleting its baseline unit:

* A golden singleton triangle has third wall t_b, with
  b.t_b=-x and every other old contact zero.
* The golden opposed-fork quadrilateral has a new opposite fork
  c_b, with a.c_b=-x, e.c_b=lL.c_b=-1/2, and every other old
  contact zero. Its second new wall is the d from the
  [explicit golden quadrilateral](arithmetic_dn_golden_opposed_fork_quad_survivor_20261019.md).

The analogues at (a,e) are denoted t_e and c_e and exchange b,e.
These are definitions of exact profiles. We do not assume that they
are globally realizable, or that they exhaust every possible deficient
face. The statements below are conditional on the corresponding walls
being actual walls of one compact polytope.

**Capacity theorem.** At most one of these two angles can carry any
of the listed golden exceptional profiles. More explicitly, the pairs
(t_b,t_e), (c_b,c_e), (c_b,t_e), and (t_b,c_e) are all incompatible
with acuteness.

The first two exclusions have simple exact products; the mixed
exclusion needs a further calculation given below. This theorem is
independent of any adjustable weight coefficients.

## 2. Actual-wall projection onto a full finite vertex

Let G be the old positive-definite unit D_n Gram matrix. For an actual
new wall r, let v_r be its row of contacts against the old roots, and
write

    q_r=v_r*G^-1*v_r^T.

The positive old span has a one-dimensional timelike perpendicular.
The old finite vertex lies strictly inside every new wall: it already
has its n incident facets and the compact chamber is simple. Thus
q_r>1, and the timelike component of each actual new wall has the
same orientation relative to that vertex. For two such unit roots,

    r.s = v_r*G^-1*v_s^T
            -sqrt[(q_r-1)*(q_s-1)].              (1)

This sign is essential; the two timelike components cannot be chosen
independently. The formula follows by separating the positive old
projection from the common future timelike perpendicular line.

For the D_n inverse entries we need

    G^-1_aa=2*N,
    G^-1_ab=G^-1_ae=N,
    G^-1_bb=G^-1_ee=n/2,
    G^-1_be=N/2,
    G^-1_(a,lL)=2,
    G^-1_(b,lL)=G^-1_(e,lL)=1,
    G^-1_(lL,lL)=2.

They follow directly from the ordinary D_n fundamental weights and
the factor-two change for unit roots.

## 3. Two triangles or two quadrilaterals

For two contact rows exchanged by the b,e automorphism, with all
other contacts equal, their projection costs coincide. If the contact
rows differ by h*(e_b-e_e), their squared difference cost is

    h^2*(n/2+n/2-2*N/2)=2*h^2.

Consequently their cross projection is q-h^2. Equation (1) gives

    r.s=1-h^2.                                  (2)

For t_b,t_e, h=x<1, so their product is

    1-x^2=(5-sqrt(5))/8>0.

For c_b,c_e, h=1/2, so their product is 3/4>0.
Both contradict the nonpositive scalar product of distinct acute
outward normals. This also recovers the golden quadrilateral capacity
from its original note.

## 4. A triangle at one angle and a quadrilateral at the other

It remains to exclude c_b together with t_e; the other mixed case
is its b,e interchange. Define

    Q=q_(c_b)-1,
    T=q_(t_e)-1,
    C=v_(c_b)*G^-1*v_(t_e)^T.

Using the entries in Section 2, then 4*x^2=2*x+1, gives

    Q=2*N*x^2+(N+2)*x+(N+2)/8
      =2*(N+1)*x+(5*N+2)/8,

    T=(N+2)*x^2/2-1
      =(N+2)*x/4+(N-6)/8,

    C=N*x^2+(N+4)*x/4
      =(3*N+4)*x/4+N/4.

All three numbers are positive at N=17,18. In particular these
profiles have nonzero timelike components as required for (1).
Expanding and reducing the one remaining x^2 gives the exact identity

    C^2-Q*T=[22*(N+1)*x+7*N+3]/16>0.             (3)

For transparency, before the final reduction the numerator of
64*(C^2-Q*T) is

    4*N^2*x^2+(-2*N^2+88*N+88)*x-N^2+28*N+12.

The relation 4*x^2=2*x+1 cancels both quadratic-in-N terms and
leaves 88*(N+1)*x+28*N+12, proving (3).
Since C>0, equation (3) implies C>sqrt(Q*T). The forced actual
product in (1) is strictly positive. Thus c_b and t_e cannot be
distinct actual outward walls. Their prescribed contact rows differ,
so they cannot coincide either. This proves the mixed exclusion and
the capacity theorem.

## 5. Field separation from the other two explicit profiles

The silver opposed-fork quadrilateral has
x=beta=(1+sqrt(2))/2 and forces ground field Q(sqrt(2)), as shown
in its [separate note](arithmetic_dn_silver_opposed_fork_quad_survivor_20261019.md).
It therefore cannot occur in the same arithmetic polytope as either
of the golden profiles, which force Q(sqrt(5)).

For the golden singleton triangle, this last field assertion follows
without guessing the ground field. Its two-cycle puts x^2 in K,
hence Q(sqrt(5)) is a subfield. The old D_n plus t_b determinant
ratio is 1-n*x^2/2<0 physically, and is nonzero. Any nonidentity
embedding of a larger totally real K fixing Q(sqrt(5)) would preserve
that negative determinant, contrary to arithmetic conjugate positive
definiteness. Thus K=Q(sqrt(5)). Its nonidentity projection cost is
n*(3-sqrt(5))/16, lying strictly between 0 and 1 at n=19,20;
the triangle is not eliminated by that local arithmetic test.

The distinct triangle with a.c=-1/2, b.c=-sqrt(3)/2 and c perpendicular
to W likewise forces Q(sqrt(3)). Its whole I2(6) angle has the
additional unit under the already proved high-label field rule;
the bare distance-two baseline alone must not be substituted for
that rule. In particular it is not a demonstrated counterexample to
the current complete compact20 payment.

## 6. Saved golden cross-role audit

The current rooted golden source registry is
`results/fresh_sqrt5_dimension21_a4_source_orbits_20261016_cap2.json`.
Its only `A1xA18` and `A1xA19` orbit records have marked root 1.
Their literal simple-edge near-end path is

    marked 1 -- 5 -- 6 -- 4 (endpoint).

The other neighbor of 1 is 11. Thus this saved mark lies three edges
from the nearer endpoint. In both exceptional faces here, the needed
large path starts at e and has the marked fork next to it:

    endpoint e -- marked a (or c) -- remaining long chain.

Its nearer-end distance is ONE. These rooted orbits differ; a weight
on the saved three-edge role cannot be transferred to the required
second node. The literal record prefixes and labels were read directly;
no graph reconstruction program was executed. The independent angle
router in `scripts/check_compact21_full_face_primal_independent_20261018.py`
encodes the whole target component and the actual mark, chooses an
`A1x` orbit key, and uses zero when that key is absent from its
coefficient map. The mismatch is therefore substantive, not a choice
of the path's displayed direction.

This verifies the absent saved role. It does not change the old
compact20 proof: that proof did not delete the D_n unit and its
original triangles and quads were already paid. It shows why the
unchanged saved golden cross registry cannot automatically pay this
new proposed deletion.

## 7. What a balanced activation argument would still need

The capacity theorem means that charging both old fork/whisker angles
for these exact exceptional profiles would overcount. If a complete
consumer classification showed that these are the only unpaid
exceptions to a fork-unit deletion, one could retain the unit on an
exceptional face and attempt to refund the other unit at each whole
D_n vertex. At most one golden profile, or at most one silver
profile, could obstruct the two angles; the field separation prevents
mixing the golden and silver families.

That is a possible use of the proved capacity, not yet a valid
activation theorem. It still needs the complete triangle and quad
case lists, actual bonus roles, and simultaneous retention for every
selected face. In particular the actual cross destination in the
A1+A_(n-1) corner is the second node of its long path, one edge from
the nearer endpoint. An orbit for a different node cannot be used
as its payment.

Even a justified saving of one unit per whole D20 vertex would only
reduce its bare baseline from 38 to 37. A uniform bound of 35 needs
further balanced savings and all other rank-20 forests checked. No
such uniform bound, global realization, or new computation is
claimed by this note.

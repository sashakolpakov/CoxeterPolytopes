# Integral opposition frames in dimension ten

20 October2026. Source-only paper lemma; no scientific calculation,
census, optimizer, Lean or build ran. This supplies an exact model of
one opposition component, not a dimension-nine theorem or an arithmetic
classification of the ambient polytope.

Let q be a finite vertex of a right-angled10-polytope, q²=-1, with
orthonormal outward normals n1,...,n10. Put

    h=3q-sum_i ni,  h²=1,  h.q=-3.

The actual P3 lemma and full-wall feasibility argument are inherited
from the dimension10 manuscript. Across its opposition indexed by a
three-set S, q becomes qS=2q-sum_S ni, the three corresponding
normals become m_jk=q-nj-nk, and the other normals are unchanged.

## Reflection and lattice

Set rS=q-sum_S ni. Then rS²=2 and h.rS=0. Reflection
x -> x-(x.rS)rS sends q to qS, each ni in S to the m indexed
by the other two members of S, and every other ni to itself.
It preserves the integral lattice

    L=Zq + sum_i Zni,  Gram diag(-1,1,...,1),

and fixes h. Thus every finite vertex and incident normal of this
opposition component has integral coordinates in the initial frame.
This is a constructed local lattice; facets outside the component
need not belong to it. It does not make P arithmetic.

The vector h is primitive of norm1 and characteristic: for x in L,
x² is congruent to x.h modulo2. Consequently

    L=Zh orthogonal-direct-sum L0,  L0=L intersect h-perp,

and L0 is even unimodular of signature(9,1). An explicit basis is

    a_i=ni-n_(i+1), i=1,...,9,
    r=q-n1-n2-n3.

Indeed x=tq+sum zi ni lies in L0 exactly when sum zi=-3t;
subtracting tr leaves an integer zero-sum spatial vector, generated
by the a_i. These ten roots have square2. The a_i form the A9
chain, while r meets only a3, with product-1. This is the explicit
T_(2,3,7), usually called E10, Gram. No classification theorem is
needed for the basis or the displayed reflection identities.

## Ideal nine-simplex frame and external walls

For every incident normal n in this component one has h.n=-1.
Thus w=n+h is a primitive future-null vector in L0. At one vertex,

    wi.wj=-1 (i distinct j),
    sum wi=3v,  v=q+3h,  v²=-10.

Future orientation follows from q.wi=-3. Their ten ideal endpoints
therefore form a regular ideal nine-simplex frame in the hyperbolic
nine-space h-perp. Under an opposition its three null vectors become
v-wj-wk, while the other seven stay fixed.

Every halfspace of an incident component facet contains the entire
hyperbolic nine-space h-perp strictly, because for future timelike
x perpendicular to h, x.n=x.(n+h)<0. A finite set of only such
facets cannot bound finite volume: choose a future-null point in
h-perp different from their finitely many w-rays. Every wall is
strict there, giving an open ideal-boundary patch and hence infinite
volume. Thus a finite-volume P needs further facets outside this
component-normal set. This necessary boundary fact gives no uniform
number of missing oppositions by itself.

In fact some actual facet r must satisfy h.r>-1. Otherwise put
t=h.r<=-1 and u=r-th. Then u lies in h-perp, u²=1-t²<=0
and u.(q+3h)=r.q+3t<0, so u is nonzero future causal. Its
halfspace again contains h-perp; choosing a null ray there different
from all finitely many possible null u-rays makes every wall strict.
That gives the same open ideal patch and infinite volume. This proves
existence of a boundary wall beyond the opposition-normal level, but
does not count its contacts or give a uniform component-boundary loss.

## A discrete height from one component facet

Fix Z=wj for one incident normal in the initial frame. At any finite
vertex of the component, put F=-Z.q and ti=-Z.ni. All are integers.
Since each ni+h is future-null, ti>=0. The null/frame equations give

    sum ti=3F,  sum ti²=F².

Therefore F²>=sum ti=3F, so F>=3. At F=3 every ti is0 or1,
with exactly nine1s and one0. Conversely ti=0 means Z is parallel
to ni+h. Both are primitive integral future-null vectors, hence equal;
then F=-Z.q=3. Thus height3 is exactly incidence with the chosen
component facet, and every ti is positive at heights greater than3.

At height4, writing ui=ti-1 gives sum ui=2 and sum ui²=2.
Thus precisely two ti are2 and eight are1. Every formal triple move
has height F'=2F-sum_S ti. At height3 there are84 formal moves
of height3 and36 of height4. At height4 the formal moves have
heights3,4,5 in counts8,56,56 respectively. These are exact
formal-frame counts; missing actual oppositions are not supplied by
these numbers.

This offers a discrete boundary/transport formulation. A finite actual
component has a maximum height, but excluding its outward moves at
that maximum alone gives no uniform global loss: the size of that
level relative to all vertices has not been bounded. No claim about
all component boundaries, universal spectral radius, compactness or
dimension reduction follows from this note.

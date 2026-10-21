# Selective distance 10 and controlled supporting extensions

**Complete internally checked computer-assisted proof.** Primary generation
and a fresh saved-output replay pass. A separate mathematical/source review
and independent exact-arithmetic reconstruction also pass; the linked review
records its final source bindings and saved replay. No unresolved mathematical
step is identified in the argument below.

The argument below gives dimension at most **116** for every finite-sided,
finite-volume real hyperbolic Coxeter polytope. It requires no arithmeticity,
compactness, simple-cusp, affine-factor-rank, or facet-count assumption.
Its compact-facet consequence is ambient dimension at most **87**.

The [primary source](../scripts/fresh_selective_distance10_dimension116_20261017.py)
and [exact certificate](fresh_selective_distance10_dimension116_20261017.json)
record the complete finite split. A [separate review](fresh_selective_distance10_dimension116_20261017_review.md)
reconstructs the arithmetic and checks the geometric argument. The preceding
[dimension-126 proof](fresh_selective_distance11_dimension126_20261017.md)
supplies the inherited incidence and all-rank estimates and the excluded
higher-dimensional tail.

## Actual face setup

For an actual quadrilateral, write W for its supporting white roots and
U,V for the two opposite pairs of black side roots. Normalize every root
to square 2. W and W together with any one black root are positive
definite; W with an adjacent black pair is positive semidefinite. The
opposite sides do not meet, even in the closed quadrilateral, so W+U and
W+V are indefinite. Choose minimal indefinite cores L and M. Each
contains its two black roots, and deleting either black root gives a
positive-definite diagram.

Assume there is no path of length at most 10 from U to V with white
internal nodes. The complete order-6 through order-10 catalogue has
384 marked contexts, with maximum marked radius 7. Cores of order at
most 4 have radius at most 2. Order-5 cores have radius at most 3 and
require a separate treatment in the disjoint case below.

An external supporting white contacting node i of a black-deletion
Gram A has squared projection at least (A inverse)ii. This uses the
entrywise nonnegative spherical inverse, nonpositive contacts, and
minimum nonzero Coxeter magnitude 1. If either black-deletion inverse
diagonal is at least 2, positive definiteness of W plus the remaining
black root is impossible. Equality 2 is excluded, since this is an
edge-support test. Of the 527 high-order ports at marked distance at
least 3, exactly 161 pass.

Whenever a required positive-definite diagram fails, the producer
records a nonnegative nonzero vector with nonpositive square. Symmetric
elimination of a Z-matrix yields such a vector at its first nonpositive
pivot; the saved vector and square are checked exactly. Its sign makes
the obstruction valid after arbitrary additional nonpositive contacts
or larger bond magnitudes are restored. Label-4 diagrams are rationalized
by a positive diagonal congruence with squared scales 1 or 2. Inverse
diagonals are returned to the norm-two gauge before comparison with 2.

## Disjoint high-order cores, including multiple cross edges

Disjoint indefinite cores cannot be orthogonal in an ambient Lorentzian
space. Each cross bond joins white ports with marked distances summing
to at least 10. Each distance is therefore at least 3.

The 161 allowed high-order ports give 1,478 ordered candidate joins.
In 363, the simple cross bond already gives a nonnegative vector of
nonpositive square in W. The other 1,115 possible edges occur in 1,007
ordered marked-context pairs: 905 admit one edge, 96 admit two, and six
admit three. Thus one cannot assume uniqueness of the cross bond at
this distance threshold.

Every possible cross bond is nevertheless simple. Each endpoint has a
three-node white path. Coefficients (1,2/3,1/3) on that path have square
at most 4/3, since native bond magnitudes are at least 1. The two paths
together have square at most 8/3-2 beta for cross magnitude beta;
additional acute contacts only decrease this value. Positive definiteness
forces beta<4/3, leaving only the simple Coxeter magnitude 1. This
excludes all larger finite labels, parallel contacts, and dotted contacts
without a label cutoff.

All 35 underlying single-edge Grams have exactly two negative directions.
For the context pairs admitting multiple edges, the source tests every
subset of size at least two. All 120 such choices have a nonnegative
nonpositive-square witness in W plus one black root. Thus every nonempty
subset is impossible, and no unspecified cross bond can repair the
single-edge inertia obstruction.

## Complete small-order supplement

The radius inequality alone does not exclude an order-5 core here:
radius 3 plus radius 7 permits a connecting edge at distance 11. This
supplement is essential. Orders at most 4 are excluded because
2+7+1<=10, and two small cores are also excluded.

Number the black roots 0,1 and a white root at marked distance 3 by 4.
The other whites must form the path 4-3-2. Neither black can contact
3 or 4. Connectivity leaves exactly three tree attachments:

    0-1-2-3-4;
    1-0-2-3-4;
    0-2, 1-2, 2-3, 3-4.

The only additional possibility would contain the triangle 0,1,2.
Three nonright Coxeter contacts give a nonpositive square on the
positive all-ones vector, with equality only for the affine simple
triangle. In the equality case its attached white 3 makes the proper
four-node diagram indefinite. Thus the triangle is impossible.

All edges in these five-node trees have finite labels at most 5. Indeed,
each edge belongs to a connected proper four-node principal diagram,
which must be positive semidefinite. An edge of magnitude at least
sqrt(3), together with an adjacent nonright edge, has a three-node Gram
determinant at most zero. Strict inequality contradicts semidefiniteness;
equality is the affine [6,3] path with positive kernel. Its connected
fourth-node extension again contradicts semidefiniteness. Parallel and
dotted contacts are excluded by the same argument. Hence enumerating
labels 3,4,5 is exhaustive, rather than a search cutoff.

The 3 times 3^4=243 marked tree recipes split exactly as follows:

| Outcome | Recipes |
|---|---:|
| Indefinite proper principal diagram | 209 |
| Full diagram positive semidefinite | 14 |
| A black deletion is not positive definite | 4 |
| Far-port inverse test fails | 12 |
| Surviving far ports | 4 |

For the semidefinite tests, the one- and two-node minors are positive
from the allowed bond alphabet; every larger proper principal minor and
the full determinant are checked exactly. The four survivors are the
two black orientations of the path [5,3,3,3], and the fork with respective
edge labels (3,4,3,3) or (4,3,3,3). Their far-port inverse diagonals are
(4/5,3/4), (3/4,4/5), (1,4/5), and (4/5,1).

After the high-order inverse test, the only radius-7 ports are L25/node 9
with black pairs {1,2}, {1,3}, {2,3}. In each small/high pair the only
possible cross edge is 4-to-9. The same white-tail argument makes it
simple. Each of the four resulting 15-root Grams has exact inertia
(13,2,0), excluding all 24 ordered marked joins. The two labels 5 in
this small supplement are handled exactly over the real algebraic field,
not by rounding their cosines.

## Complete overlap split and fixed core Grams

Every shared white root has marked distances summing to at least 11,
so both distances are at least 4. Small-order cores therefore cannot
occur in an overlap. Each marked high-order context has at most four
eligible nodes. Enumerating injective partial identifications gives
22,881 candidates and 4,862 Gram-consistent overlaps.

External-neighbor inverse tests exclude 4,296, and explicit nonnegative
vectors in W plus one black root exclude 254. The remaining 312 contexts
are exactly the following; mixed counts include both orientations.

| Core intersection | Contexts | Core order | Extension mechanism |
|---|---:|---:|---|
| L10/node 5 = L25/node 9 | 6 | 15 | Actual radical relation |
| L18/node 7 = L25/node 9 | 6 | 17 | Actual radical relation |
| L22/node 4 = L25/node 9 | 18 | 18 | At most one zero-residual root |
| L22/node 8 = L25/node 9 | 54 | 18 | Actual radical relation |
| L22 terminal edge = reversed L25 terminal edge | 18 | 17 | Residual 1/2; at most one attached A1 |
| L25/L25 shared node 9 | 114 | 19 | Actual radical relation |
| L25/L25 reversed terminal edge | 60 | 18 | No additional core contact |
| L25/L25 reversed terminal three-node path | 27 | 17 | At most one zero-residual root |
| L25/L25 reversed terminal four-node path | 9 | 16 | Residual 0 or 2/3; at most an A2 tail |

Every unspecified bond between the exclusive core parts is tested.
Exactly 21,393 create a path of length at most 10. The other 168 have
nonnegative nonpositive-square witnesses in W plus one black root.
Stronger or additional nonpositive bonds cannot remove these obstructions.
Thus the actual Gram is fixed in all surviving cases before the extension
arguments are applied.

## Preferred corners and actual radical relations

For each surviving context, the certificate checks all four adjacent
corners and chooses one whose marked roots lie in a connected component
of order at most 17. Every corner is decomposed into finite or affine
components. A new white can contact a finite component only at normalized
inverse-diagonal entries at most 2; it is orthogonal to an affine component
by the positive kernel. The two omitted black roots are included separately
in the possible contact domain.

In the radical cases the Gram has one negative direction and one null
direction. The radical gives an actual root relation: its linear
combination is orthogonal to a timelike vector in the core span and has
square zero, so it vanishes in the positive-definite perpendicular space.
For rationalized Grams, the positive diagonal scaling changes coefficients
but not their signs or zeros.

The selected contact domains contain no zero-radical coordinate. The
certificate assigns the two radical signs to the two marked pairs and
checks that any pair of opposite-sign contacts gives a marked path of
length at most 10 through the new root. Acuteness and the actual relation
force both signs whenever there is any contact. Hence every additional
supporting white is orthogonal to these cores. Their marked components
are isolated with order at most 17.

## Nonsingular cores and the residual norm budget

All other surviving cores are nonsingular Lorentzian. In each preferred
corner, the allowed contact domain splits into two end sets, each within
distance 4 of its own marked pair. A root contacting both ends gives a
path of length at most 10. Within either end, the full inverse-Gram block
is strictly positive.

For a norm-two new root with core contact vector b, write its orthogonal
decomposition as e=e_core+w. The perpendicular space is positive definite,
so

    b^T G^{-1} b <= 2,     w^2 = 2-b^T G^{-1}b.

Exact inverse entries show that two nonzero contacts at one end exceed
the budget. Every permitted diagonal is greater than 1, so the next
Coxeter magnitude sqrt(2) also exceeds the budget. Consequently a
nonorthogonal supporting white must have a sole contact -1 at a listed
node i. Its projection square h_i=(G inverse)ii is one of

    2, 3/2, 4/3,

with residual square rho_i=0,1/2,2/3, respectively.

There is at most one such core-contacting white in total. For each pair
of candidate nodes i,j, the exact certificate checks

    ((G inverse)ij)^2 > rho_i rho_j.

If the inverse entry is positive, Cauchy's inequality makes the actual
product positive, violating acuteness. If negative, the actual product
is strictly negative and therefore creates an edge between the two
new whites. The certificate checks that this edge gives a marked path
of length at most 10. The same test with i=j excludes two distinct
roots at a single candidate node. In the zero-residual case the root
is uniquely determined in the core span as well.

It remains to control whites attached indirectly through that one root.
All other supporting whites are core-orthogonal. Their Gram E is a
positive-definite principal submatrix of W. If c is the extra root's
contact vector with these whites, projection gives

    c^T E^{-1} c <= rho.

Each contacted connected component costs at least 1/2. A connected
component of order at least two costs at least 2/3, and one of order at
least three costs at least 3/4. These bounds follow by projecting onto
a connected two- or three-node principal subdiagram through a contacted
node; the explicit inverse inequalities are proved in the
[supporting-extension lemma](fresh_distance10_extension_20261017.md).
Larger diagrams and additional contacts cannot lower projection cost.

Thus rho=0 permits no attachment, rho=1/2 permits only one A1, and
rho=2/3 permits only one A1 or A2. This bounds the entire connected
attachment, with no ambient-rank cutoff. The core-contacting root plus
this tail adds at most one, two, or three nodes, respectively. Every
selected marked component still has order at most 17, as checked
context by context. This explicitly handles supporting extensions with
positive Euclidean residual; unchanged isolation is not assumed.

At an ideal corner, the inherited factor-fiber theorem places the entire
varying affine factor inside W plus the two side roots. The bounded
marked component is therefore the actual component relevant to payment.
Omitted spectator cusp walls cannot enlarge it. The argument includes
nonsimple cusps.

## Face payments and dimension inequalities

Pay 1 on a corner pair if its distance is at most 10 or its connected
component has order at most 17, and 0 otherwise. The classification
above pays every quadrilateral at least 1. The inherited triangle
argument gives two paths of length at most 8, so every triangle receives
at least 2. Nonnegative payments suffice for larger polygons. Every
k-gonal two-face therefore receives at least 5-k.

The inherited all-rank graph counts give finite vertex cost at most 10n.
At an ideal vertex with affine factor ranks r_i and link vertex count
V=product(r_i+1), factor-fiber multiplicity gives

    C(v)/V <= sum min(10,r_i/2).

The shape inequality has normalized ideal curvature

    kappa/V = 3d/4 - (5/2) E|sum X_i|,    d=n-1,

where X_i is centered uniform on r_i+1 consecutive integers. The proved
convex-order transfer applies below and above cap 20, preserving cost
and not increasing curvature. It terminates by increasing the sum of
squared ranks. Reduced profiles have copies of 20, at most one part
1 through 19, and at most one part at least 21.

For d=116 through 125, all 1,015 reduced profiles have positive
curvature-minus-cost. Integer convolution and the separate moment formula
agree exactly. The smallest margin is **145234/4084101**, at ranks
(20,20,20,20,20,16). The finite slack beta_n-10n is positive for n=117
through 126, and the global shape inequality has its strict positive
term 5(n-1). These dimensions are excluded. The completed 126 theorem
excludes all higher dimensions, giving 116.

At ambient dimension 116, the profile (20,20,20,20,20,15) has margin
-673055/16336404. The same scalar weights do not prove 115; that local
profile is not a finite-volume realization.

For a compact facet of dimension d=n-1, finite payments give
beta_d<10(d+1). For d=87+2t and 88+2t, eight times the opposite
difference has coefficient lists (12,176,4) and (96,180,4), positive
for all t>=0. Hence d<=86 and the ambient dimension is at most 87.

Published classifications, actual-face and cusp incidence, the all-rank
distance counts, convex-order transfer, and shape inequality are explicit
inherited mathematical inputs. The arithmetic certificates do not replace
the written geometric proof. This package makes no finite-volume
construction, external-refereeing, proof-assistant, or publication-priority
claim.

Primary final generation (session 75668) and fresh verification (session 7019)
both exited zero with status
`PASS_SELECTIVE_DISTANCE10_CONTROLLED_EXTENSIONS_DIMENSION116`.
The verifier reconstructs the complete classification and exact inequalities:

```sh
python3 scripts/fresh_selective_distance10_dimension116_20261017.py --verify
```

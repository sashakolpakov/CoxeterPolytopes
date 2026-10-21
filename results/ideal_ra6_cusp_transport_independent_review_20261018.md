# Independent review of the original4425 cusp-transport bound

**PASS.** This is a separate mathematical and source review of
[the original proof](ideal_ra6_cusp_transport_research_20261018.md) and
[its checker](../scripts/fresh_ideal_ra6_cusp_transport_20261018.py).
The original proof's universal conclusion
`V >= ceil(3645F/187)`, hence `V >=4425` for `F>=227`, is correct under
the explicitly inherited ideal-five carrier inputs. No original file
was changed. The stronger reciprocal estimate found during this review
is deliberately recorded in a separate packet.

## Actual triangle walls and endpoint normalization

The centered box frame has orthogonal positive e_i, future null v,z,
v.z=-1, and n_i^±=±e_i+(w_i/2)v with unit normals. In the two-face
obtained by omitting n_i^epsilon_i, the other m-1 normals support its
plane. Its sides through v are n_i^+ and n_i^-; their third side has
an actual ambient facet normal c_i, orthogonal to that support. All three
polygon vertices are ideal, so c_i has product-1 with each n_i^±.
Those equations force its z coefficient2/w_i and its e_i coefficient0.
Orthogonality to the other support normals fixes every remaining e_j
coefficient; square one fixes the v coefficient. This derives formula(1)
for arbitrary positive widths, not merely the88 saved model frames.

For a good edge all these actual normals exist simultaneously. Direct
products give c_i.c_j=0 for i!=j and n_i+c_i=v_epsilon/w_i. The displayed
v_epsilon is future null and is distinct from v because v.v_epsilon=-2.
Its expression as a positive sum of two actual normals proves feasibility
against every other wall; against either summand its product is zero.
It is therefore the second ideal endpoint of the actual supporting edge,
not an unproved formal incidence. The2m displayed incident normals are
distinct and fill the rectangular cusp link at that endpoint. Opposite
pair sums are v_epsilon/w_i, so the new widths are precisely1/w_i in
this chosen endpoint normalization. In particular cubicity is preserved
in both directions along good edges.

## The normal-product gap and the noncubic count

The product gap used in the proof has an actual-face justification. If
two outward unit normals have product in(-1,0), their positive definite
two-root Gram is an M-matrix with nonnegative inverse. Project a strict
timelike interior vector off their span. The coefficient vector is
negative, and all other outward-normal products are nonpositive, so the
projection remains strictly inside every other halfspace. It is timelike
and lies on an actual nonright ridge, contrary to right angles. Products
of distinct outward normals are nonpositive by acuteness. Thus only0 or
values at most-1 occur.

If two good sign vertices differ in direction k, the normals c_i before
and after reversal and the common c_k are all actual and distinct. They
are orthogonal to the n-3 unchanged supporting normals. If the first two
had product zero, these would be exactly n mutually orthogonal unit
normals. The same interior projection gives an actual finite vertex,
contrary to the entirely ideal hypothesis. Thus
`1-2w_k^2/w_i^2 <= -1`, and w_k is maximal. This counts the normals
correctly: three triangle walls plus n-3 support normals, not n-1.

A noncubic cusp has a strictly nonmaximal direction. The perfect matching
in that cube direction pairs all32 sign vertices, and at most one endpoint
of each matching edge is good. Hence g<=16. A link vertex is bad exactly
when at least one incident link edge is not triangular. Counting missing
edge endpoints gives2(80-t)>=32-g, so t<=64+g/2<=72. Here t counts
triangular two-faces at this cusp, equivalently triangular link edges.

## Canonical cubic frame, invariant and unweighted graph sum

At a cubic cusp scaling all widths to1 makes v canonical: it is the sum
of either opposite pair of actual unit normals. The span of the e_i is
unchanged by signed permutations, so its complementary z and
h=z-v/8 are canonical as well. Across a good edge the normalized pair
sums are exactly v_epsilon. The proposed new frame satisfies
e_i'=n_i-v_epsilon/2 and z'=h+v_epsilon/8, with z'^2=0,
z'.v_epsilon=-1 and z'.e_i'=0. Thus h is the same actual vector at both
ends, including around cycles in a good-edge component.

Since h is spacelike of square1/4, its orthogonal complement has future
null vectors. Any such Z gives E=-Z.v>0: Z cannot be proportional to v,
as Z.h=0 while v.h=-1. This avoids a hidden zero-energy vertex. Writing
Z=(E/8)v+Ez+E sum p_i e_i gives sum p_i^2=1/4 and the exact neighbor
energy ratio3/2-sum epsilon_i p_i. Every formal endpoint is future null
with h-product-1, so its ratio is positive even when it is not an actual
good neighbor. Adding those missing formal ratios is therefore legitimate.

The second/fourth moment identities, fourth-moment upper13/80 and
mu^2-d^2>=1 give the stated mean reciprocal bound71/90. The reciprocal
identity and coefficients are correct. With f=1/E, each actual oriented
good edge contributes f(u)/f(v). The sum at a vertex is bounded by the
32-sign sum. Each undirected edge contributes x+1/x>=2 to the global
sum. Thus the result bounds the **unweighted** total good degree,
not merely a weighted degree. Componentwise choices of Z cause no
problem because the inequality is summed within each finite component.

## All face-incidence factors

The5-cube link has32 vertices and80 edges. Good degree is the count g
of good ambient edges at the cusp. Hence cubic triangular corner counts
sum to at most64V_c+(1136/45)V_c/2=(3448/45)V_c. Noncubic cusps contribute
at most72 each, which is smaller. Every ideal triangular two-face has
three cusps, so the total is **3T**, with T counting triangles. This
factor3 is present in the original report and the final bound.

Independently, each cusp has degree32, giving a_1=16V. Summing the
ideal-five identity over all facets gives3a_3=30V+6a_4-6F because a
three-face has exactly three containing facets and a four-face has two.
Euler then gives a_2=25V+a_4-F. Polygon/cusp incidences equal80V, and
every nontriangle has at least four vertices, so T>=4a_2-80V. Combining
this with3T<=(3448/45)V yields a_4-F<=187V/135. Finally the inherited
56-facet minimum in ideal dimension five gives2a_4>=56F; therefore
27F<=187V/135. The rational ratio3645/187 and ceiling4425 atF=227
are correct. No facet incidence or triangle corner is counted twice.

## Verification and scope

The saved checker replay passed in chunk5ce0fa in0.022s, peak RSS26,755,072
bytes, under the shared lock and59s alarm. It checked88 exact frame
models and the stated exact moment/coefficient consequences; those
models are correctly described as identity checks, not a classification
of possible widths. The independent review above supplies the universal
geometric reasoning. It inherits the already reviewed56 and227 bounds
from the linked earlier ideal-six packet; it does not independently
reprove that earlier carrier chain or external literature. No Lean,
new example or dimension exclusion is asserted. No process remains.

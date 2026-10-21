# Dimension ten assembly: the remaining analytic contract is explicit

20 October 2026. **Assembly candidate, not yet a promoted ceiling.**
The ten-direction rational-envelope certificate and full independent
composition review are still required. The following proof identifies
every connection between the separate packets. It assumes no average
bound from one distinguished corner or from sampled frames.

## Statement and count input

Assume the structured subset estimate below for every subset size
k=0,...,11. Then no mixed finite-volume right-angled hyperbolic
eleven-polytope exists. Together with the already reviewed mixed
ceiling eleven, this gives a mixed ceiling ten. The compact and entirely
ideal cases retain their smaller published ceilings four and six.

For a hypothetical mixed eleven-polytope use the actual counts
f>0 finite vertices, B canonical P3 three-faces and

    E=sum_(three-faces other than P3)(facets+cusps-10)>=0.

The complete reviewed ordinary count theorem is

    3722B>159195f+2794E.                          (1)

Its source is `right_angled_mixed_eleven_three_face_corner_balance_20261020.md`;
the full root review, including its one-cusp and R7 inputs, is
`review_mixed_eleven_one_cusp_and_corner_balance_root_20261020.md`.
The all-ideal-edge error theorem used both there and below gives

    H=sum_(finite vertices q)[C(k(q),3)-t(q)]<=2E, (2)

where k(q) counts incident edges with ideal other endpoint, and t(q)
counts actual P3 three-faces through q. Thus sum_q t(q)=2B. Equation
(2) counts actual noncanonical all-ideal-edge corners; no graph
completion is used to assert that every ideal triple is P3.

## One frame on each actual opposition component

At q write n_1,...,n_11 for its orthonormal outward normals. An actual
P3 triple S has opposite finite vertex

    q'=2q-sum_(i in S)n_i.

The new incident normals are n_i for i outside S and

    m_i=q-sum_(j in S, j!=i)n_j       for i in S.

They are actual facet normals by the reviewed P3 attachment/feasibility
proof. Hence the vector h=3q-sum_i n_i is unchanged across the
undirected opposition edge q--q'. Its square is11-9=2>0.
For each connected component choose ONE nonzero future null vector
Z perpendicular to h, possible because h is spacelike in signature
(11,1). The same choice applies to every finite vertex in that
component, including both ends of every opposition edge. Components
consisting of an isolated vertex cause no difficulty.

Set F(q)=-<Z,q>>0 and p_i=-<Z,n_i>/F(q). The null equation and
<Z,h>=0 give exactly

    sum_i p_i=3, sum_i p_i^2=1.

For an actual triple its ratio is

    r_S=F(q')/F(q)=2-sum_(i in S)p_i>0.           (3)

The reversed actual edge has ratio r_S^-1. All formal triples in
this frame also have positive ratio by the full-frame Cauchy range,
but formal triples need not represent actual chamber points.

## The precise analytic input and actual-subset deletion

For any eleven real coordinates satisfying the two frame equations,
and any k-element coordinate subset I, the required contract is

    sum_(S subset I, |S|=3) (r_S^-2-37/100)_+
             <=756/25+(27/100)C(k,3).            (4)

At q let I be its actual ideal-ending directions. Its t(q) actual
P3 triples form a subset of the C(k,3) formal triples on I. For
each actual triple,

    r_S^-2-1 <=(r_S^-2-37/100)_+-63/100.

The positive terms can be enlarged to all formal triples on I.
Consequently (4) implies the pointwise inequality

    sum_actual(r_S^-2-1)
      <=(9/25)(84-t(q))+(27/100)(C(k,3)-t(q)).    (5)

This deletion uses only nonnegativity. It neither replaces the actual
triangle graph by a complete graph nor substitutes a local maximum
for the ordinary global average.

## Reversal, excess and the contradiction

Summing the left side of (5) over all finite vertices pairs every
actual undirected opposition edge in its two directions. By (3), its
contribution is

    r^-2+r^2-2=(r-r^-1)^2>=0.

Thus (2) and sum t=2B give

    0<=(9/25)(84f-2B)+(27/100)H
      <=(9/25)(84f-2B)+(27/50)E.

Equivalently

    2B<=84f+(3/2)E,
    4B<=168f+3E.                                (6)

Multiplying (1) by two and (6) by1861 gives respectively

    7444B>318390f+5588E,
    7444B<=312648f+5583E.

They would require5742f+5E<0, impossible because f>0 and E>=0.
This is a contradiction from ordinary counts and one valid null vector
per component, with no extremal-height substitution.

## Verification boundary

The small subset sizes k<=6 have a direct elementary proof in
`right_angled_mixed_k10_subset_moment_reduction_20261020.md` Section7.
The k7,8,11 contract is proved by written moment arguments and a replayed
rational certificate in
`right_angled_mixed_eleven_structured_small_and_full_subsets_20261020.md`.
The k9 contract is proved by an analytic equality neighborhood and a
replayed complete rational interval certificate in
`right_angled_mixed_eleven_nine_direction_subset_20261020.md`; an
independent global analytic proof is being prepared.

The k10 moment and rational-envelope reductions are written in
`right_angled_mixed_k10_subset_moment_reduction_20261020.md` and
`right_angled_mixed_k10_rational_envelopes_20261020.md`. Their polynomial
positivity obligations must pass exact checking and the complete
argument must receive independent written review before the current
ceiling is changed. Lean and the final TeX manuscript follow that
mathematical verification.

# Independent review of the mixed dimension-ten assembly

20 October2026. **PASS: no mathematical flaw found in the complete written
assembly, subject to its named published and inherited geometric inputs.**
I read the complete assembly, ordinary-count sources, one-cusp proof,
small-three-face and attachment interfaces, all scalar subset proofs used
below, and all three exact workers. I independently rederived the critical
geometry, moments, clipping corrections, ordinary count elimination and
final inequality. No scientific program, certificate replay, optimizer,
graph census, Lean or TeX build was run in this review. Reading small
existing replay ledgers and streaming source hashes are not new proof runs.

The assembly source's statements that the k10 certificate and k9 global
proof were still pending are historical source status. The current saved
workers/replay ledgers and written proofs discharge those obligations as
described below. This review does not edit that source or any frozen file.

## 1. Source-bound scope

|Read source|SHA256|
|---|---|
|`right_angled_mixed_dimension10_assembly_candidate_20261020.md`|`70ef7de40e9fbd62dbe763d64d258c40ab777a42d7550e1608074e4c29aa0326`|
|`right_angled_mixed_eleven_three_face_corner_balance_20261020.md`|`73b8073843171d4bf11ad0fc11829b0855df36b2a78307fc024cd6d26dd450df`|
|`right_angled_one_cusp_ten_facet_exclusion_research_20261020.md`|`19590da4197ee04940df141a492ff89d69aea44b49bf87fe3c63699c252b56d4`|
|`right_angled_mixed_nu11_all_ideal_corner_exclusion_research_20261019.md`|`fc047b78e9a8fc9eaf283bd5d689b69d72a8bbf0976a4765600ffbe8625955b0`|
|`right_angled_mixed_nu10_edge_transition_research_20261018.md`|`0a4282db244068b3a94941819369e51bdd621d37537fa334eb621848b40e6d31`|
|`right_angled_mixed_weighted_triangles_research_20261018.md`|`2e9fac6d7c82493e8e7e1e0172235cc3b910f5f3b05d7ded0d7056bdae217a60`|
|`right_angled_mixed_polygon_three_face_coupling_research_20261018.md`|`f2949232322ec28ba5c2d10a97a30a0a47a519c3dc6061e88dd95c4cf0b31a41`|
|`right_angled_mixed_eleven_structured_small_and_full_subsets_20261020.md`|`04d3b8fa433365210600234709894d777803e7db8591cbb32876a0796e87a1ce`|
|`right_angled_mixed_eleven_nine_direction_global_proof_20261020.md`|`bf288985ade5acf93e13faafb25544b0642b38f13e422ad38a3f1223e9706d15`|
|`right_angled_mixed_k10_subset_moment_reduction_20261020.md`|`8c2051b3be0f692e5097eedd81cf494f15313037d36b0f6625ca78cd01aabd34`|
|`right_angled_mixed_k10_rational_envelopes_20261020.md`|`7fdcad68b9b4ae113c482e9e17266a6f5e9aea0d21409eb7db332ec915627e45`|

These paths are in `results/`. The older coefficient-two carrier-successor
proof was also read, but the direct H budget below makes its D3 recognition
and its sharper nu12/13 classification unnecessary for this assembly.

The conclusion is absence of a mixed finite-volume right-angled
eleven-polytope. Face heredity, together with the published compact
ceiling4 and entirely ideal ceiling6, gives the universal tail and hence
the mixed ceiling10. It does not require the earlier internal ceiling11
as an additional theorem. The geometry and real-variable arguments remain
written proofs, not consequences of the exact scalar certificates or a
fully formalized Lean theorem. No external priority or novelty audit is
claimed here.

## 2. Actual P3 opposition, one null vector and reversal

The P3 attachment interface is genuinely geometric. At a finite vertex
q, write its unit outward incident normals n_i orthonormally, q^2=-1.
An actual bad triangle in directions i,j has actual opposite wall
m_ij=q-n_i-n_j. Three pairwise bad triangles give three mutually
orthogonal such walls and the future unit timelike point

    q'=2q-n_i-n_j-n_k.

Every other wall r=a*q-sum c_l*n_l has a>0,c_l>=0. Acuteness against
the three actual m walls gives the three pair bounds c_i+c_j<=a,
hence c_i+c_j+c_k<=3a/2 and r.q'<=-a/2<0. Thus q' is an ACTUAL
finite vertex. The original ideal neighbors and the two finite vertices
give precisely the canonical P3 three-face. No complete triangle of an
abstract graph has been mistaken for an arbitrary realizable face.

At q', the new three incident normals are

    m_i=q-sum_(j in S, j!=i)n_j,

while all outside normals are unchanged. They are actual facet normals.
Their sum on S is3q-2sum_S n, so

    3q'-sum new normals=3q-sum old normals=h.

This is literal equality of ambient vectors, including outward signs.
Reversing the same new triple gives2q'-sum_i m_i=q. Different old
triples have different opposite points, by independence of the old
orthonormal normals; even a multigraph interpretation would preserve
the paired directed-edge sum. Each P3 has exactly two finite vertices,
so sum_q t(q)=2B without an orientation factor being lost.

In dimension11, h^2=2>0. Its perpendicular is Lorentzian and has a
nonzero future null vector Z. Choose ONE Z for EACH actual finite
opposition component; do not choose it again at neighboring vertices.
Then F(q)=-Z.q>0 everywhere in that component. Writing
Z=F(q)*(q-sum p_i*n_i) gives

    sum p_i=3, sum p_i^2=1,
    F(q')/F(q)=2-sum_S p_i=r_S.

The reversed ratio is exactly r_S^-1 because the SAME Z is used.
There is no equal-height, maximal-height or vertex-weighting premise.
Isolated components also admit this choice and simply have no edge terms.

Centering all eleven p_i at3/11 gives squared deviation2/11. The centered
indicator of a triple has norm squared24/11, so ALL formal triples have

    (13-4sqrt3)/11<=r_S<=(13+4sqrt3)/11,
    11/20<r_S<20/11.

In particular the analytic denominators are positive even on formal
triples that are not chamber vertices. Negative individual p_i are allowed.

## 3. The exact ordinary three-face budgets

For a three-face write F for facets,s for cusps,p for finite vertices,
nu=F+s, and A_j for its finite vertices with j ideal-ending edges.
Euler and the two link degrees give p+4s=2nu-4. The reviewed small
types have

    P3: A3=2;
    Q8: A1=8;
    R7: A0=1,A2=3.

I reread the complete Q8/R7 arguments. Zero excess at Q8 excludes bad
triangles by the attachment lemma;32 polygon corners on8 facets force
all to be(3finite,1ideal) quads, giving one ideal-ending edge at each
finite vertex. At R7, excess1 and at least4 triangles force one all-ideal
triangle, three bad triangles and three such quads. Its finite-edge graph
is a three-leaf star, so its center's A0 contribution MUST accompany the
three A2 contributions. All three small types satisfy A2-3A0=0.

For E=sum_nonP3(nu-10), define

    D=sum_faces(A2-3A0),
    H=sum_nonP3 A3.

The direct proof H<=2E avoids D3 recognition entirely. If a non-P3
face has A3>0, its all-ideal-edge corner gives three distinct polygon
excess units: each incident nontriangle contributes1, and each incident
bad triangle gives a DIFFERENT positive-excess partner across its ideal
edge. The attachment lemma prevents a bad partner, and the normals
q-n_i-n_j distinguish those partners and exclude incidence at q. Thus
F+2s>=15. At nu<=10 this contradicts s>=15-nu and4s<=2nu-5.
At nu11 it forces s4,p2, so A3<=2. At nu12, A3>=5 would require
p>=5 and4s>=3A3>=15, hence4s>=16 and p+4s>=21>20. Thus A3<=4.
For nu>=13, p>=A3 and4s>=3A3 give
A3<=(nu-2)/2<=2(nu-10). Summation proves H<=2E. No planarity,
classification of D3, global finite-edge restriction or full graph
completion is needed for this error budget.

The second bound D>=-30E uses only the small types and the one-cusp
lemma. At nu11, F+2s>=12 gives F<=10, and the one-cusp exclusion
removes F10,s1; thus p<=10. At nu>=12, p<=2nu-4<=10(nu-10).
Since A2-3A0>=-3p, summation gives D>=-30E.

I checked the one-cusp proof's nontrivial boundary step. The four cusp
quadrilaterals have eight DISTINCT finite vertices because opposite cusp
planes meet only at the cusp and adjacent ones share only their cusp edge.
The four successive pentagon labels cannot repeat at opposite positions:
two repeated-facet intersections would close a triangular pentagon boundary.
They therefore give four distinct pentagon adjacencies. The two remaining
pentagons each meet all five others, giving nine further distinct edges.
That contradicts the total eleven pentagon-pentagon edges. This establishes
F>=11 with one cusp, without using the stronger published classification.

Consequently the exceptional incidence count X, omitting P3,Q8, satisfies

    X<=18E+4R,
    R=sum_nonP3(8-F).

The pointwise budget difference is20nu-144-8F. At nu10,R7 it is0;
at nu11 the one-cusp exclusion leaves F<=9; at nu>=12,F<=nu suffices.
R is signed and is never treated as nonnegative.

## 4. Balanced ordinary count and strict elimination

At each actual finite vertex let k be its ideal-ending edge count,
m=11-k,t its actual P3 degree, and put

    U=C(m,3), V=C(k,2)*m, W=C(k,3).

Every triple with zero or two ideal directions is exceptional, hence
X(q)>=U+V. Exactly D(q)=V-3U and H(q)=W-t. In particular W-t
counts NONCANONICAL actual all-ideal-edge three-face corners; it does
not presume every formal ideal triple is P3. The twelve exact values
of40U+28V-16(165-W), for k0 through11, are

    3960,2160,972,288,0,0,180,432,648,720,540,0.

I independently recomputed these integer values by hand. They imply
31X(q)>=16(165-t)+3D(q)-16H(q). Ordinary summation, using sum t=2B,
D>=-30E and H<=2E, gives

    31X>=16(165f-2B)-122E,
    680E+124R>=2640f-32B.

I also independently checked all exact elimination coefficients. Incidence
and three-face Euler give

    E=B-10a3+9a2+11520c,
    R=8a3-9a2-2B,
    a3=(9/2)a2-(165/4)f-5760c.

Here11520 is the number of2-faces of the ten-cube cusp link. Thus

    R=-5B/4-165f/8+5760c-3E/4.

The strict weighted polygon row a2>(25/2)f+1930c is the dimension11
case of the source-bound Dehn-Sommerville/weighted packet. Its strict
sign uses the published82-facet bound for finite-volume RA eight-faces,
the absence of compact RA eight-faces, and elementary ten-face facet
counting. The printed lambda=2123/1280 gives exactly1930, not a rounded
coefficient. The same identities yield the STRICT cusp substitution

    360c<B-75f/2-E.

Substituting R first gives
587E-123B-10395f/2+714240c>=0. Since714240/360=1984>0, the strict
cusp substitution puts its left side strictly below
1861B-159195f/2-1397E. Hence

    3722B>159195f+2794E.                         (L)

The signs, strictness and use of ordinary counts are correct. There
is no transport-weighted quantity hidden in(L).

## 5. Every scalar subset contract is complete

The required statement, for ANY k-coordinate subset I in ANY full frame,
is

    sum_(T subset I,|T|=3)(r_T^-2-37/100)_+
       <=756/25+(27/100)C(k,3), 0<=k<=11.       (C)

I independently checked all written reductions used to establish it.

For k<=5, the full-frame lower ratio6/11 bounds at most ten summands
by673/225 each, totaling less than756/25. At k6 complementary triples
have ratio sum>18/11 and each ratio>6/11. Convexity and monotonicity
of the positive-part reciprocal bound each pair by12461/3600; ten
pairs are strictly below891/25, the exact k6 target.

For k7,8, Cauchy on the omitted coordinates gives the exact concave
variance envelope, and the moment coefficients in the source follow from
the one-, two- and three-index inclusion probabilities. The negative
quartic coefficient is used with the correct Cauchy LOWER bound on the
fourth coordinate moment. The rational interval worker lowers the mean,
raises the variance and square roots, and uses a justified minimum of
global, variance or zero clipping corrections. Its S intervals cover
the complete real domain; a skipped interval has negative maximum
variance and contains no frame. No sampled direction is substituted
for a full interval.

At k11 the exact mean13/11 and variance2/11 give the joint fourth
moment bound15V^2/121-T3^2/5. Replacing its common remainder factor
by5 is valid because the WHOLE bracket is nonnegative; equivalently
it bounds the nonnegative actual fourth moment. The remaining quadratic
in T3 increases up to the sharp allowed endpoint bound1/15. The
variance clipping correction uses the positive gap5/11. This yields
the exact rational target2493/5500.

For k9 I prefer the independent global analytic proof bf288 above.
I rederived its moments V/4,-3T3/28,5V^2/28-3T4/28, its cubic
bound and the global clipping inequality(c-r^-2)_+<=3(r-mu)^4.
The clipping is corrected on ALL ratios, not only the positive-part
support. The reciprocal remainder coefficient is at most1280/121.
After V<=x(4-11x)/2, x=delta/3, the cleared difference is x^2R(x).
I independently expanded its polynomial and checked

    R(x)>=233/24-(1643/2904)*16=635/968>0,

because(4-11x)^2(1+x)^5 decreases on[0,4/11]. This proves the
stronger mean bound63/100-delta/18. Its equality frame is exactly
nine selected1/3 entries and two omitted zeros; no unproved optimizer
classification is used.

For k10 I checked both moment and rational-envelope sources completely.
The cubic/quartic expression multiplying the common remainder K has
lower coefficient259/2700>0, so replacing the ENTIRE K by7 is valid.
The negative-v branch maximizes at the allowed cubic endpoint and uses
the square-root tangent only in the positive cubic contribution. The
nonnegative-v branch completes the concave square over a containing
real domain. The zero-correction and variance-correction domains cover
the complete outside-coordinate interval. Their denominator signs are
positive. All three cleared polynomials agree with the printed formulas.
Points with negative formal V outside the true domain are used only as
polynomial test points, not as moment frames. The Bernstein conversion
is exact and positive coefficients establish every closed subinterval.

In particular none of these arguments assumes nonnegative p_i, uniform
selected entries, a single local optimum or a maximum-height average.

## 6. Actual-subset deletion, reversal and the final strict gap

For the actual ideal-ending set I(q), all t(q) ACTUAL P3 triples lie
among its C(k,3) formal triples. For every actual triple,

    r^-2-1 <=(r^-2-37/100)_+-63/100.

Enlarging only the nonnegative positive-part terms to all formal
triples, then applying(C), gives exactly

    sum_actual(r^-2-1)
       <=(9/25)(84-t)+(27/100)(C(k,3)-t).         (Ulocal)

This is the key safe deletion step. Formal noncanonical triples are
charged to the explicit H error; they are not declared actual moves.

Sum(Ulocal) UNWEIGHTED over all actual finite vertices. Each undirected
opposition edge uses one common Z at both ends, so its two contributions
are r^-2+r^2-2=(r-r^-1)^2>=0. Thus

    0<=(9/25)(84f-2B)+(27/100)H
      <=(9/25)(84f-2B)+(27/50)E,
    4B<=168f+3E.                                (U)

This does not replace an average by the value of a favorable or
height-maximizing corner. Every actual finite vertex participates.
There is no missing graph-reversal factor, orientation count or
weight-normalization hypothesis.

Twice(L) has left coefficient7444B and right318390f+5588E.
Multiplying(U) by1861 gives right312648f+5583E. Their strict
incompatibility is exactly5742f+5E<0. A mixed polytope has f>0,
and E>=0 by the complete small-three-face gap. Contradiction.

## 7. Saved exact checks, separate from this written review

|Worker|SHA256|
|---|---|
|`check_mixed_eleven_balanced_and_subset_rows_20261020.py`|`e403456a9f880af672407173e51b18e4e3b03b3cb474cdb1030ce54ba381f830`|
|`check_mixed_eleven_nine_direction_global_20261020.py`|`357070546d69aaa634c01292fb15c174f6d85552e87cb42be6adbf4d2cdb0b62`|
|`check_mixed_k10_rational_envelopes_20261020.py`|`d0520c374cb3e8689dedec78284dc50ea8d7d9c182b47d71ab2e81430bcf99fb`|

These paths are in `scripts/`. I read all three workers without importing
or executing them. The first reconstructs all twelve balanced integer
rows and every k7/8 rational interval. The second checks the global k9
polynomial identities and rational comparisons. The third constructs all
three k10 polynomials with Fraction arithmetic, then all840 full closed
intervals and12550 Bernstein coefficients. No numerical positivity or
stored solver answer is trusted by their read-only replay paths.

The existing separately generated certificates have hashes

* `mixed_eleven_balanced_and_subset_rows_20261020.json`:
  `eb0b79f20b856e13ec4ec78ba57f147328eb07ad08b846d8c4d87baf2134e257`;
* `mixed_eleven_nine_direction_global_20261020.json`:
  `bc2b9a390237f1de48922de0c091dcfb7140fd37e964935fe6962cc828846090`;
* `right_angled_mixed_k10_rational_envelopes_20261020.json`:
  `a04a84b1d997b1479d5478408e277684228afd1b1d65e67114bf747c7f4c6163`.

The small inspected replay ledgers all record PASS_BOUNDED_RESEARCH_CHECK:

* `mixed_eleven_balanced_subset_replay_20261020.json`, hash
  `ce89fa94118d7e572a0e0aade9a931fcc869b0ef7e3ae61d15e90dbd25c2eb25`;
* `mixed_eleven_nine_direction_global_replay_20261020.json`, hash
  `23cd3bfadc0ab32c9ee1667b21e2e5aaef7bd049f4ab729c08bc160c12b90c02`;
* `mixed_eleven_k10_envelope_replay_20261020.json`, hash
  `f2a2f71c2b328bc274780668634cb379e45be60cb1e47406fb66144dfa6f7066`.

Those PASS statements are attributed to the saved executions. They do
not mean this reviewer ran them, nor that a scalar checker proves actual
P3 attachment or the one-cusp geometry. The written review above supplies
the missing source-to-geometry and source-to-count composition check.

## 8. Universal tail and final scope

The argument rules out mixed dimension11. Compact dimension11 is already
excluded by the published compact ceiling4, and entirely ideal dimension11
by the published ideal ceiling6. Hence there is no finite-volume RA11
polytope of any vertex type. Any hypothetical higher-dimensional one has
a finite-volume right-angled11-face, by the standard face-heredity input,
which would contradict that conclusion. This gives the higher-dimensional
tail without importing the earlier mixed11 proof and avoids a circular
ceiling assumption.

Within the stated published and inherited geometric inputs, the complete
written dimension-ten proof and its exact analytic obligations PASS this
independent review. It remains an internal mathematical proof, with geometry
and real-valued transport not fully Lean-formalized. Manuscript editing,
build scheduling, catalogue updates and user-facing promotion belong to
the coordinator; no such mutation or scientific job was made here.

# Root review: metric facets and good-prime tree arguments

20 October 2026. Independent paper review; no scientific computation,
Lean compilation or TeX build was needed. These are restricted structural
statements, not improvements of the main dimension ceilings.

## Metric Coxeter facets: PASS

Reviewed source:
`coxeter_metric_facet_angle_rigidity_20261020.md`, SHA256
`f688e520f3ff7ac43225aaacd89413c03f5dd3d2d008d7e4dcdcc48326ac83fb`.
The independent resource-lane review is
`review_coxeter_metric_facet_angle_rigidity_resource_20261020.md`, SHA256
`84061c99b780b2522bab20b794a22328e08e147e2220a18f2627f6543941e6dc`.

I rederived the projected-normal formula and all eight marked A3/B3/H3
values. Cyclotomic degree excludes every integer Coxeter angle among
these values except the B3 value cos²(pi/4)=1/2. Thus a nonright intrinsic
angle other than pi/4 forces both parent walls perpendicular to the
metric facet. The compact arithmetic [4,3,5] simplex realizes the pi/4
exception: its displayed Gram matrix has the stated leading minors,
spherical proper diagrams and integral reflection coefficients.

If every metric facet is Coxeter, every finite vertex component has rank
at most two. Conversely such components make every metric face Coxeter:
white-root elimination cannot couple distinct components. At dimension
nine, a vertex therefore contributes at most four nonright two-face
corners. For a compact k-gon the number of nonright corners is at least
5-k: two for k=3, one for k=4, and the inequality is automatic for k>=5.
Summing gives B>=5f2-36f0, whereas B<=4f0.

For clarity, the last identity can be checked without a script. If
h0,...,h9 is the dual simplicial h-vector, then
f0=sum(h_i), f2=sum(binomial(9-i,2)*h_i). Pairing the terms with indices
i and 9-i gives

    5f2-40f0 = 100h0+60h1+30h2+10h3 > 0.

Dehn--Sommerville symmetry, h_i>=0 and h0=1 supply the contradiction.
Higher dimensions descend to an actual nine-face. The ceiling eight is
therefore proved for this restricted class; its sharpness and published
priority have not been established. Arbitrary compact dimension29 and
arithmetic compact dimension20 remain unchanged.

The application to a metric Bugaenko eight-facet is also valid. Every
node of its eleven-wall diagram touches a 3- or 5-edge, so every parent
side wall is perpendicular to it. That application was already present
in earlier repository arguments; the new note supplies a complete
field-free local classification and its sharp exception.

## Good-prime trees: PASS with a limited conclusion

Reviewed sources:

* `general_coxeter_good_prime_two_face_trees_20261020.md`, SHA256
  `fc1d24889f35fd73aa557392e81ec51486e60da61fcc37b1f470879e55673b6b`;
* `compact_coxeter_good_prime_global_tree_action_20261020.md`, SHA256
  `3574174f4f6dd1d7eb995264dad16ce1189e8f006372bcad1e3121e859449d53`.

The rank-one projector proof is valid at odd residue characteristic.
Contractivity propagates one norm value across each connected white
component; the integral inverse Gram supplies the reverse inequality.
Unimodularity and self-duality fix that value to zero. This identifies
the entire white fixed set with the complement's building, rather than
only embedding some fixed norms. A codimension-two white set leaves a
three-dimensional complement, hence a tree or a point.

Finite edge and vertex stabilizers supply the claimed consecutive
subtree intersections. Chordality and the subtree Helly property give
the bounded-diagonal triangulation. The trace formula then proves the
stated obstruction to two bad opposite pairs. The explicit Lambert
quadrilateral passes the real compactness and local-valuation checks;
one bad opposite pair is possible already on a quadrilateral. Appending
white roots is only a local model, not a finite-volume completion.

For the global successor, I checked the cited norm/lattice-chain model
in Cornut, Section3, printed pages6--8. Its orthogonal-space construction
applies to the nondegenerate spaces used here. The constant unimodular
white lattice in every radius interval 0<=t<1/2 ensures that each fixed
tree uses edges of the same ambient chain complex. This is the essential
extra point beyond abstract CAT(0) convexity. The equivariant map from
the tessellation's two-skeleton to the building's one-skeleton lifts
coherently to the universal covering tree because that two-skeleton is
simply connected. A fixed point upstairs would give one in the building.
All white systems must be good at the SAME place.

This conclusion is not an arithmeticity or dimension theorem. In
particular no fixed-point property for these Coxeter groups has been
proved. A Coxeter diagram with an infinite pair i,j already has the
visual splitting

    W = W_(S\{i}) *_(W_(S\{i,j})) W_(S\{j}).

Consequently the mere existence of a nontrivial tree action cannot supply
the missing dimension estimate. A controlled separator or another
geometric/global restriction is still needed. The final open-status
paragraph of the earlier local note is superseded only as to existence
of the conditional global tree action, not as to any quantitative gain.
Park this approach unless such additional input is found; no new
manuscript or Lean expansion is justified by the current statements.

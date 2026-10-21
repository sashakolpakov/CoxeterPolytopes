# All direct A1×A1 couplings of the fixed H18 source escape

Every coherent direct A1×A1 coupling of the fixed40-facet H18 hybrid has infinite volume. This covers every finite inducing-contact label: the initially unbounded single-contact label families are reduced by exact source-anchor constraints, not by a search cutoff. A single exact spacelike point survives the entire remaining168-normal profile bank.

The [producer](../scripts/fresh_h20_a1a1_h18_coupling_20261015.py) and [certificate](fresh_h20_a1a1_h18_coupling_20261015.json) pass generation and saved replay. They reconstruct all1600 source Gram entries, the24-node forcing component, the literal source anchors and the common escape. There is no assignment enumeration, full ambient vertex census or finite-volume example.

This is a new direct-coupling gate. The earlier [arbitrary-height bare lift](fresh_h18_general_height_lift_20261014.md) assumed zero in the convex hull of specified height sets and therefore did not directly cover this support-wall construction. The [A2](fresh_h20_a2_h18_coupling_20261014.md) and [B2](fresh_h20_b2_h18_coupling_20261014.md) results also have different inverse-Gram constraints. All earlier artifacts remain unchanged.

## Profiles and the complete pair formula

Take two orthogonal norm2 support roots. Their fundamental weights have squares1/2 and mutual product0. Normalize every inducing ambient wall to norm2. For a source root n of square q, an inducing wall has the form

`r=√(h/q)n−a0ω0−a1ω1`,

where each nonzero `a_c=2cos(π/m_c)`, and

`h=2−(a0²+a1²)/2>0`.

A zero contact is a right angle. Positive projected norm rules out parallel or dotted support contacts. It initially permits a single arbitrary finite label `m≥3`. If both contacts are nonzero, positivity permits only

`(3,3),(3,4),(4,3),(3,5),(5,3)`.

Indeed if both labels are at least4 the projection is nonpositive; if one is3, the other must be below6. The ordinary profile has both contacts zero and h=2.

For any two source facets the complete metric rule is

`(r_i,r_j)=√(h_i h_j/(q_i q_j))G_ij+(a_i0 a_j0+a_i1 a_j1)/2`.

Thus orthogonal source facets can both be special when their nonzero contact colors are disjoint. The A2 argument that any orthogonal special pair is forbidden would be incorrect here.

## Three orthogonal facets force the24-node component ordinary

On a source π/3 edge an ordinary/special pair has product `−√(h/2)∈(−1,0)`. No acute Coxeter angle between norm2 roots has such a product. Therefore ordinary/nonordinary status agrees at both ends of every source π/3 edge.

The source has one connected π/3 component with24 vertices and27 edges. It contains three pairwise orthogonal source roots, indices `0,2,4`. If any component root were special, all24 would be special. These three roots would then require three pairwise disjoint nonempty contact supports in a two-color set, which is impossible. All24 component roots are consequently ordinary. The certificate checks the component through an explicit spanning tree and all three orthogonal pairings.

The16 remaining facets consist of the three core roots `17,23,24`, eleven other norm-one roots and two norm-two roots. Each core root has a forced ordinary norm2 anchor with old product−1. Every other remaining root has such an anchor with old product−2. Their literal indices and products are saved.

## Anchor constraints remove the entire label tail

For a single contact of label m, `h=2sin²(π/m)`. The squared ambient product with the chosen ordinary anchor is h for a core root, 4h for another norm-one root, and 2h for a norm-two root. These are products squared, not squared cosines; the allowed finite values begin at1,2,`φ²`,3.

For a core root, m=3 gives `1<h<2`, between the π/3 and π/4 values; m≥5 gives `0<h<1`. Hence only the single label4 remains. Among double contacts only `(3,3)` remains, with h=1.

For another norm-one root, single labels3 and4 give dotted and parallel anchor relations, while label6 gives π/4. Label5 gives squared product `5−√5`, strictly between `φ²` and3, hence between π/5 and π/6. Labels7 and8 give squared product strictly between1 and2. To close every larger label, put `u=√(7/8)` and `f(c)=4c³−3c−1/2`. Since

`f(u)=(u−1)/2<0`

and f is increasing above1/2, its root `cos(π/9)` is larger than u. Therefore `2sin²(π/9)<1/4`. Monotonicity gives `0<4h<1` for every m≥9, again forbidden. No maximum label has been assumed.

For a norm-two root, single labels3,4,6 give π/6, π/4, π/3 respectively. Label5 gives a squared product strictly between1 and2. Every m≥7 has `0<2h<1`. The double profiles `(3,5),(5,3)` give `h=(3−√5)/4`, so their anchor product has magnitude below1 in all three root classes; they are forbidden.

The complete remaining necessary profile lists are therefore:

| Source class | Profiles, with0 denoting a right contact | Count per facet |
|---|---|---:|
| Forced24 | `(0,0)` | 1 |
| Core3 | `(0,0),(4,0),(0,4),(3,3)` | 4 |
| Other13 | `(0,0),(3,0),(0,3),(4,0),(0,4),(6,0),(0,6),(3,3),(3,4),(4,3)` | 10 |

Including the two support walls gives `24+3·4+13·10+2=168`. Other pair constraints may reject further choices; the escape covers this larger relaxed bank.

## A uniform strict escape

The literal source point x inherited from the B2 experiment is rechecked against the raw source Gram. It has

`x²=(3857−2800√2)/529<0`

and normalized gap at least1 on every potentially special source facet: `(x,n_i)≤−√q_i`. Its pairings with the other facets are nonpositive. The sum h of all40 source normals is future timelike and strict against every source wall.

Set

`E=(2x+h/100000,−ω0−ω1)`.

In the norm2 simple-support basis the two transverse coordinates are `−1/2,−1/2`. Their square is1 and both support pairings are−1. For each remaining special profile, the certificate verifies

`4h_profile>(a0+a1)²/4`.

The smallest margin belongs to `(3,4)` or `(4,3)` and equals `(5−2√2)/4>0`. Consequently

`(E,r)<−2√h_profile+(a0+a1)/2<0`.

All ordinary profiles are strict because of the source-center perturbation. The exact square is

`E²=19945826604799/661250000000−(87506613029/4132812500)√2>0`.

The point `C=(100h,−ω0−ω1)` is a common strict timelike point for the same168-normal bank, and `(E,C)<0`. The segment from C to E crosses the future null cone while all wall inequalities remain strict. Every coherent direct assignment therefore has an open ideal-boundary region and infinite volume.

The theorem concerns the two retained support walls and facets inducing the fixed40 source facets. Additional walls with new source projections or nonpositive face projections, changed retained support walls and different source metrics remain outside it. In particular this does not exclude subsequent cap constructions with genuinely new source directions; it shows that changing only the direct A1×A1 inducing profiles cannot provide the missing finite-volume H20 construction.

# Affine-tree quadrilaterals with short outer end blocks

There is a further all-rank paid branch for the proposed mixed finite/ideal 3n rule. In an all-simple finite-volume Coxeter polytope of dimension n>=14, let an affine B/C/D corner have central-path marked roots a,b at distance P>=4. Write the diagram as

`F_r — a — (P−1 white path nodes) — b — F_s`,

where each end is B_r (r>=1) or D_r (r>=3). Here B1 means A1 joined to a by a bond4; D3 means A3 joined at its middle node. If **r+P>=9 and s+P>=9**, the opposite marked roots c,d are adjacent. The quadrilateral is therefore paid.

This removes the large-outer-end restriction from the [earlier corridor lemma](fresh_simple_cusp_affine_tree_corridor_20261017.md). It leaves only 13 central-position types, up to exchanging the two ends, together with fork-leaf and fork-center positions. It is not a full quadrilateral theorem or a dimension bound.

The [checker](../scripts/fresh_simple_cusp_affine_tree_short_ends_20261017.py) and [data](fresh_simple_cusp_affine_tree_short_ends_20261017.json) check the closed inverse formulas, the five distinct join costs, and literal positive kernels for the two exceptional affine exclusions. All-rank completeness comes from the arguments below, rather than the range of the confirming matrix calculations.

## Complete extension classes

Deleting a leaves the short outer F_r and an inner B/D component of rank s+P>=9. Every new connected component contains whole old components. An extension touching only the outer component must be finite: an affine component plus the untouched inner component cannot be the connected affine diagram of a simple ideal vertex. Its contact vector may otherwise be arbitrary, including unbounded labels for an A1 outer component. An extension touching only the inner component is its unique finite long-end extension by a bond3. Extending it affinely would again leave a disconnected corner. The isolated case is also allowed.

If c touches both components, its connected diagram has all n nodes. The complete finite/affine classification now gives a classical diagram: exceptional diagrams have order at most9. The large retained B/D component determines its long port. At the outer end the possibilities are precisely:

* B1: bond3 for a finite classical extension or bond4 for an affine extension;
* B2: bond3 at either endpoint;
* B_r, r>=3: bond3 at its long port;
* D3: bond3 at either path endpoint for a finite extension, or at the middle for an affine extension;
* D4: bond3 at any of its three leaves;
* D_r, r>=5: bond3 at its long port.

The contact to the inner component is always bond3 at its long port. This list follows by deleting a node from classical paths and single/double forks; there is no remaining exceptional or cycle case with the large B/D component. A proposed type incompatible with the other end may be retained in the table as a harmless relaxation. In particular c.b=0 because P>=4. The same reasoning gives d.a=0. It suffices to contradict c.d=0, since a nonzero nonpositive product is an adjacent pair under the rule.

## Exact end identities and arbitrary outer-only contacts

Normalize roots to norm2. Let J=G(F_r)^−1 and let v be the positive contact vector from a to F_r. Except for B1 it is the unit long-port vector; for B1 it is sqrt(2) times that vector. In every case v^T J v=1. The cusp coefficients on F_r are t_j=(Jv)_j.

For B_r, r>=2, the diagonal of J is 1,2,...,r−1,r/2, and the t_j are 1,...,1,1/sqrt(2). For B1 the diagonal is1/2 and t_1=1/sqrt(2). For D_r the diagonal is 1,...,r−2,r/4,r/4, and t_j are 1,...,1,1/2,1/2. The checker displays and multiplies the full closed inverse matrices in 78 examples; sparse multiplication of the same formulas proves them for every rank.

For an outer-only contact vector w>=0 put q=w^T Jw and t=v^T Jw. Its finite corner requires q<2. Every nonzero coordinate of w is at least1, and all t_j>=1/2, so t>=1/2. Cauchy–Schwarz in J gives q>=t^2. If x=−a.c>=1, then

`(2−q)/(x+t)^2 <= (2−t^2)/(1+t)^2 <= 1`.

If x=0, the same bound follows from q+t^2>=2. With two or more contacts, t>=1 and Cauchy–Schwarz already proves this. With one contact of magnitude at least sqrt(2), the displayed diagonals prove it as well. The only remaining possibilities are one bond3 at B1, the non-port endpoint of B2, or a fork leaf of D3,D4,D5,D6.

For B1/B2 that extra leaf and the first old path nodes give the five-node affine F4 path with bonds(3,4,3,3). For D3..D6 the added leaf extends one fork arm to length2; the other fork arm has length1, and the old long direction supplies length5. This is the nine-node affine E8 with arms(1,2,5). The full dimension n>=14 supplies the five simple steps before the other decoration. These are induced diagrams: the new root has precisely its one outer contact and x=0. Both are proper affine subdiagrams, forbidden by the [all-simple affine exclusion](fresh_simple_cusp_affine_exclusion_20261017.md). The checker verifies their positive kernels and positive definite one-deletions. Thus every outer-only profile contributes at most1, without an angle or rank cutoff.

## Heights of isolated and inner-only roots

Give the old cusp kernel coefficient1 at a,b and on the intervening path. For an isolated c its height is h=−c.nu=x. The value x=0 would make c orthogonal to the full simple cusp, hence proportional to its null vector. The value x=1 supplies a proper affine B/D diagram: retain the whole outer end, a, and the two simple leaves c and the first inner path node. This diagram has order r+3<n. The value x=sqrt(2) gives the four-node affine star with bonds(4,3,3), unless the outer bond itself is4, when the two bond4 contacts give the three-node affine C2. Hence the conservative bound h>=sqrt(2) is valid and 2/h^2<=1.

For the inner-only bond3 profile, h=x+1. At x=0 the same affine B/D construction uses the outer end extended by a, with the first inner node as its branch point, and c and the second inner node as its leaves. Its order is r+4<n. At x=1, a, the first inner node and c form an affine A2. Therefore h>=1+sqrt(2), giving `2/h+1/h^2<=1`. These exclusions work with the entire outer block and do not require it to have large rank.

## Join costs and the corank contradiction

For a join, write q=w^T Jw and t=v^T Jw for its outer contact. Its height is h=x+1+t, so h>=1+t. The complete join list has only the following five distinct scalar data:

| Outer profile | q | t | bound on 2/h+(1−q)/h^2 |
|---|---:|---:|---:|
| B1, finite bond3 |1/2|1/sqrt(2)|less than13/9|
| original affine port |1|1|1|
| B2, exchanged endpoint |1|1/sqrt(2)|less than13/9|
| D3, path endpoint |3/4|1/2|13/9|
| D4, other leaf |1|1/2|4/3|

All entries are checked exactly in Q(sqrt(2)). No restriction on the opposite contact x beyond x>=0 is needed here.

Eliminate all white roots W. The a,b block is `(1/P)[[1,−1],[−1,1]]`, because each eliminated end costs exactly1. For c use moment l=0 for isolated/outer-only and l=1 for inner-only/join. Its residual norm is `C=c0+l^2/P`, with c0 respectively 2,2−q,1,1−q; its cross-coupling to b is l/P. The analogous quantities for d are m,d0. Since the two possible inner contacts lie at opposite ends of the intervening path, their only projection cross term is lm/P. The [formally checked determinant identity](fresh_simple_cusp_affine_a_quadrilateral_20261017.md) consequently reduces to

`det(S)/(h^2 k^2) = 1 − [g_c+g_d]/P`,

where `g_c=2l/h+c0/h^2`, and similarly for d. The preceding exhaustive classes give each g<=13/9. Therefore

`det(S)/(h^2 k^2) >= 1−26/(9P) >= 5/18 > 0`.

The actual quadrilateral star has corank1, so det(S)=0. This contradiction proves c.d<0.

## Remaining positions and verification scope

Failure of an inner-rank condition means, after exchanging ends, r+P<=8. Since P>=4, B ends give ten pairs (r,P), and D ends give three. The data save all 13 pairs explicitly; the other end may have arbitrary rank. Fork leaves and fork centers are separate unclosed positions. This is an inventory of remaining positions, not an assertion that they fail.

Generation and a fresh saved replay both pass all 78 inverse-formula checks, five exact join costs and the two positive-kernel witnesses. A separate internal mathematical/source review by the root agent passes the complete join list, the small-affine exclusions, arbitrary outer-contact bounds and the Schur estimate; that review does not duplicate the 78 matrix computations. The all-finite-corner transfer and the complete mixed quadrilateral theorem remain separate obligations. No previously frozen proof or dimension bound is changed.

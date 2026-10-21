# Review of the 76-wall changed-support replacement

**PASS.** The [frozen construction](fresh_h20_contact_closure_20261014_replacement.json) coherently replaces eight norm-418 mates by eight norm-two caps over `Q(sqrt(79))`. It retains the other 68 facets, creates eight actual `A12 + E8` finite vertices, and joins them to eight retained `D13 + D7` vertices by eight actual finite edges. It cuts the previous escape but still has infinite volume, witnessed by a new exact strict spacelike point.

The [separate checker](../scripts/fresh_h20_contact_closure_20261014_replacement_review.py) and [hash manifest](fresh_h20_contact_closure_20261014_replacement_review.json) pass generation and saved replay. The primary [replacement checker](../scripts/fresh_h20_contact_closure_20261014_replacement.py) also passes a fresh `--verify`. The review binds the source and certificates, not the [construction report](fresh_h20_contact_closure_20261014.md).

The separate implementation imports no producer arithmetic. It uses nested quadratic arithmetic in `Q(sqrt(2))(sqrt(79))`. Its exact sign test repeatedly squares opposite-sign terms, instead of using the primary checker's fixed-point radical enclosures. The radicand 79 is nonsquare in `Q(sqrt(2))`: a square `(a+b sqrt(2))^2=79` with rational `a,b` would force `ab=0` and then a rational square equal to 79 or `79/2`.

It independently checks:

- The four rank-20 linear recipes, their specified contacts, the eight resulting normals, and the norm equation `54 t^2 - 214 t + 35=0`, with `t=(107+11 sqrt(79))/54` and the recorded sign variants.
- Every one of the 572 new pairs: 144 right, 16 at pi/3, and 412 strictly dotted. Validity of the old 68 Coxeter relations remains an explicitly checked source dependency.
- The strict timelike center, all 16 specified rank-20 positive-definite vertex Grams, and all 1,216 vertex inequalities, with exact zero sets and strict nonincident signs.
- The complete finite-neighbor list for each new cap, all 160 choices of nineteen neighbors, and the unique spherical vertex for each new cap. This is a complete finite-vertex census on the new facets, not a whole-chamber corner census.
- The eight common rank-19 supports defining finite edges, the cuts of all four old-escape sign choices, the three preserved cusp-ray signs, and the final strict escape of norm `269344645863312077/40000000000000000`.

The geometric conclusions follow directly from these checks. For any facet normal `r`, projecting the strict center `H` to `r`-perpendicular gives `H-(H,r)r/r^2`. It is timelike, future oriented, and strict at every other wall because all cross pairings are nonpositive. Thus every retained and new inequality defines an actual facet. At each saved vertex, positive definiteness and the exact strict nonincident signs give an actual finite vertex; the Dynkin-type identification is additionally checked in the primary source from its literal rational Gram trees. Two such vertices sharing nineteen independent positive-definite incident normals are the two endpoints of their feasible one-dimensional face, so the eight saved segments are actual edges.

Every added wall is strict on each of the three inherited cusp rays, while the deleted mates were strict there as well. Their existing bounded links therefore remain unchanged; this review inherits the earlier link-completeness proofs. The complete per-new-cap census uses the standard fact that a finite Coxeter vertex has a simplicial spherical link. It does not exclude additional old-only vertices or classify ideal endpoints elsewhere.

Finally, the strict future spacelike point and the strict timelike center give an open feasible patch at infinity. This proves infinite volume independently of any incomplete global frontier census. The inherited conditional non-quasi-arithmetic statement is outside this arithmetic replay; no claim of a finite-covolume lattice, global completion, or exhaustion of alternative replacements is made.

```sh
python scripts/fresh_h20_contact_closure_20261014_replacement_review.py --verify
```

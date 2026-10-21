# An exact endpoint obstruction for one quartic cusp seed

The fixed quartic section `local_record=0`, `u=1/2` cannot be completed to a finite-volume Coxeter polytope by retaining its 22 cusp walls and adding projected **integral ambient normals of square 2 or 4** with nonpositive pairings to those walls. This covers arbitrarily high caps in that specified class. It is not a claim about arbitrary real normals, other ambient norm classes, other sections, or replacing seed walls.

The [exact gate](fresh_h20_quartic_orbit_cap_gate_20261014.json) leaves three necessary cap families, so it does not say that every individual cap is impossible. Its decisive conclusion is that two of the three triangle-edge types cannot acquire a finite or ideal opposite endpoint. They account for 168 of the seed cusp's 252 edge germs. A finite-volume completion must close those edges as well.

This is distinct from the original-chamber section obstruction: the present calculation permits new ambient mirrors and discards every original inequality except the 22 retained cusp walls. The earlier crystallographic-angle theorem for original-chamber sections is not used.

## Fixed section and its retained cusp

The source is the literal rational Borcherds H21 chamber and its cusp4. The active affine factors have ranks 11, 3, and 6, with 12, 4, and 7 walls respectively. The rank-three factor is `[7,23,50,200]`. Record0 of the [complete local census](fresh_h20_rank3_omitted_sections_20261014.md) omits wall7 and projects the other three to an affine A2 triangle.

Use the exact local notation

`e=e0-(L/2)D`, `e0=k23*r23+k50*r50+k200*r200`, `q=e²=e0²>0`.

Here D is the integral future null ray, the k coefficients are positive in the designated quartic field, and `k200=L`. The checker reconstructs the source Gram, these identities, and the positive kernel of the projected triangle. If `rbar` denotes projection to the section, then

`k23*rbar23+k50*rbar50+k200*rbar200=(L/2)D`.

The other two affine factors are unchanged and have positive weighted sums D and D/2. Each one-node deletion has positive-definite Gram, checked exactly. The 22 retained projected roots have Gram rank 19 and actual ambient span rank 20: the unprojected retained span has rank 21 and contains the spacelike section normal, so orthogonal projection lowers its dimension by one. The three positive affine components span the entire 19-dimensional horosphere and bound its Euclidean link. The section has the claimed cusp with all 22 link facets actual.

There are `12*3*7=252` edge germs, one for each choice of a deleted wall from each affine component. Each germ retains two of the triangle walls and seventeen other cusp walls.

## Why new compatible walls are strict at this cusp

Let an additional projected facet normal be lbar, with positive square, and write `h=(l,D)=(lbar,D)`. Assume its pairings with every retained seed facet are nonpositive. The positive affine relations force `h<=0`.

If h were zero, each positive relation would force every individual seed pairing to vanish. Then lbar would be orthogonal to the full 20-dimensional incident span in the 21-dimensional section space. Its perpendicular is the null line spanned by D, contradicting the positive norm of a facet normal. Hence every admissible additional wall has `h<0`.

Thus added acute walls cannot subdivide or remove the retained cusp link. For a finite collection they are all strict sufficiently near the cusp, and all 252 original edge germs persist. This conclusion does not assume that the original saved center remains inside the enlarged chamber.

## Complete three-root endpoint test

Suppose a new cap ends one such edge at a finite or ideal point. Its projected normal and the two retained triangle normals lie in the perpendicular of that point, so their Gram is positive semidefinite. This necessary condition includes ideal endpoints and nonsimple endpoints with several new walls.

The two triangle normals meet at label3. After normalization, write the other two off-diagonal entries as `-x,-y`, with x,y nonnegative. The determinant condition is

`3/4-x²-y²-xy>=0`.

If both contacts are nonright Coxeter angles, then `x,y>=1/2`, so equality forces `x=y=1/2`. This is precisely the affine A2 triangle with both cap labels3. Otherwise one contact is right and the other has `y²<=3/4`, allowing exactly labels `2,3,4,5,6`. A parallel contact is impossible in this three-root positive-semidefinite Gram. This elementary argument is complete and uses no field-degree angle cutoff.

Let `d=l²` be 2 or 4, `g_i=(l,r_i)` the integral old pairings, `a_i=(e,r_i)`, and `T=(e,l)` in the quartic field. Projection gives

`lbar²=d-T²/q`, and `(lbar,rbar_i)=g_i-a_i*T/q`.

The old three-root Gram is the projected positive-semidefinite Gram plus a positive rank-one matrix, and is therefore positive semidefinite too. Consequently the integer old pairings satisfy `g_i²<=d*r_i²`, together with its determinant inequality. There are exactly 120 possibilities over the three seed pairs and the two ambient norms.

For a right contact, `T=q*g_i/a_i`. The other contact is checked against all five allowed labels. Label5 is included by the exact polynomial for its squared cosine, with the larger-root inequality selecting the correct conjugate. No unsupported assumption about which cyclotomic subfields lie in the quartic field is made.

For two label3 contacts, T must be a common root of two quadratic angle equations over the quartic field. The checker computes their exact polynomial gcd. All nonconstant gcds in this census are linear; there is no unresolved quadratic factor or untested root in the field. Positive projected cap norm, nonpositive pairings, both exact angle equations, and the full three-root determinant are then checked. Fifteen literal angle-pass records remain; repeated right-contact choices intentionally produce duplicate records.

## The rational identity and its singular cases

The ambient cap also satisfies

`T=k23*g23+k50*g50+k200*g200-(L/2)h`.

For each candidate T and fixed pair of old integer contacts, equality of the four rational field coefficients constrains the third old contact and h. A full-rank two-column system is solved exactly and all four coefficient equations are checked. The resulting third contact and h must be integral, with `h<0`.

The two columns are sometimes proportional because `k200=L`. The checker does **not** discard these rank-one systems. It checks consistency and reduces the equation to an integer congruence, retaining every family that has integral solutions with h negative. This is essential: discarding the singular cases would falsely rule out all caps.

After removing repeated records, the three necessary survivor families are:

| Ambient norm | Seed pair | Contacts and field value | Remaining relation |
|---:|---|---|---|
| 2 | 23,50 | both right; `g23=g50=T=0` | `h=2*g200<0` |
| 4 | 23,50 | both right; `g23=g50=T=0` | `h=2*g200<0` |
| 4 | 23,50 | both label3; `g23=0`, `g50=-2`, `T=a200` | `h=2*g200-8<0` |

These are necessary algebraic families, not asserted realizations of integral lattice vectors or globally compatible facets. For each other seed pair, every angle-pass record forces `h=0`. Such a projected spacelike cap cannot retain the full acute seed cusp, as proved above.

## Finite-volume consequence and limits

The edge types retaining `{23,200}` or `{50,200}` have no possible ending cap in the specified ambient norm class. Each occurs with twelve choices from the first affine factor and seven from the third, giving `2*12*7=168` unclosable edge germs.

If a finite-volume polytope retained this cusp, each actual edge issuing from it would have a finite or ideal opposite endpoint. At least one added facet would be incident at that endpoint; its Gram with the retained triangle pair would satisfy the necessary test just exhausted. The seed inequalities alone have no second finite vertex and only a deficient opposite ideal endpoint along each such edge. No new facet can be incident there either, since the same positive-semidefinite test applies. Equivalently, with finitely many added walls all strict at that opposite null ray, one may continue slightly beyond it and perturb inward to obtain a strict spacelike escape. Hence finite volume is impossible within the stated retained-wall and norm-class scope.

The [producer](../scripts/fresh_h20_quartic_orbit_cap_gate_20261014.py) uses exact rational quotient-field arithmetic from the frozen global-section checker. The certificate binds this arithmetic dependency, the local algebraic record, the complete corner walk, and the raw source. There is no coefficient-height bound or numerical decision in the gate.

```sh
python scripts/fresh_h20_quartic_orbit_cap_gate_20261014.py --verify
```

The [independent review and replay](fresh_h20_quartic_orbit_cap_gate_20261014_review.md) passed using separate SymPy quotient-field arithmetic, including every old positive-semidefinite case, all accepted angle records, singular rational systems, and the seed component kernels.

The earlier [one-reflection orbit scout](fresh_h20_quartic_orbit_completion_20261014.json) remains numerical and concerns many section parameters. This certificate itself concerns only record0 at u=1/2. The subsequent [rational-cap theorem](fresh_h20_quartic_rational_cap_gate_20261014.md) now excludes every positive real shift of this same spatial direction for arbitrary rational ambient cap vectors, without a norm bound. Genuinely nonrational ambient directions, other spatial directions, and coherent deletion of seed facets remain outside that stronger theorem.


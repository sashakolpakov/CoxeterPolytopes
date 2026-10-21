# Independent review of the quadratic-five triangle restrictions

**PASS.** The [standalone checker](../scripts/fresh_arithmetic_quadratic_bound_20261015_review.py) and [saved certificate](fresh_arithmetic_quadratic_bound_20261015_review.json) independently verify the arithmetic in the [primary package](fresh_arithmetic_quadratic_bound_20261015.md). The source and proof review finds no gap in the stated ambient triangle restriction. The quadrilateral conclusion concerns exactly 23 saved recipes; it is not an exclusion of all quadrilateral completions and does not give a dimension bound.

The checker imports neither the producer nor a symbolic algebra library. It uses a separate `Fraction` implementation in the bases \(1,\sqrt5\) and \(1,\sqrt2\) over \(\mathbb Q(\sqrt5)\), exact Gaussian inversion, exact signs, and explicit square-root membership equations. It binds the frozen primary source/data and the earlier literal quadrilateral source/data. All primary source-manifest entries are checked.

## Independently reconstructed arithmetic

| Test | Independent result |
|---|---:|
| E7 contact rows | 2,187 |
| E8 contact rows | 6,561 |
| Lorentzian rows with positive conjugate | 242 E7; 286 E8 |
| Rows with at least two finite replacements | 5 E7; 3 E8 |
| Possible marked triangle pairs | 10 in each component |
| Literal quadrilateral failures | 20 outside the field; 3 nonintegral |

Every accepted contact vector, projection cost, conjugate cost, finite-deletion inequality and marked-pair index is compared with the primary data. The independently inverted literal E7/E8 matrices also give the minimum inverse diagonals \(3/2\) and \(2\). The review verifies the complete subcritical lists and the pure-label conjugate-cost exclusions.

The finite-contact alphabet \(\{2,3,4,5,6,10\}\) is an analytic all-label reduction. The real cyclotomic subfield has degree at most two; the totient condition allows only primes 2, 3 and 5 with bounded exponents and yields \(2,3,4,5,6,8,10,12\). Labels 8 and 12 have the wrong quadratic fields. Products of any two nonzero contacts close a cycle through the E tree. Exact squareclass checks then leave mixed magnitudes \(1,\phi\), or a uniform magnitude \(\sqrt2,\sqrt3,\sqrt{2+\phi}\). The latter three fail conjugate positivity. No numerical label cutoff is used.

## Orthogonal components: the necessary rescue argument

It would be incorrect to infer the ambient assertion from the isolated Lorentzian census alone. A new root can have the subcritical E7 profile consisting of a single \(\pi/3\) contact at node 5, of cost \(3/2\), and acquire the missing negative direction from another orthogonal component. The checker independently reconstructs a genuine example: an additional A1 with a \(\pi/5\) contact gives total cost

\[
\frac32+\frac{\phi^2}{2}=\frac{9+\sqrt5}{4}>2,
\qquad \frac{9-\sqrt5}{4}<2.
\]

Its finite E7 replacements are precisely nodes 3, 4 and 5.

The primary proof correctly controls every such rescue. Write the spectator projection cost as \(c\). Physical indefiniteness and conjugate positivity require \(c>1/2>c^\sigma\), so \(c\) is nonrational. In a finite corner retaining node 5, the component containing the new root and any contacted spectator is a finite tree. If it contained only labels 3 and 4, the relevant determinant ratios would lie both in \(\mathbb Q(\sqrt2)\) and the cyclic ground field \(\mathbb Q(\sqrt5)\); the spectator cost would be rational. Thus this finite component contains a label 5. A finite connected component with a label 5 and at least three nodes is H3 or H4, hence has rank at most four. It already contains the new root and at least one spectator node, leaving at most two nodes in the retained E7 component containing node 5.

The independently calculated sizes of that E7 component after deletion of nodes 0 through 6 are \(6,5,3,2,1,0,6\). Both marked nodes must therefore be among 3, 4 and 5, with deletion of node 5 included separately. Their three pairs already occur in the ten-pair census. E8 has no nonzero subcritical profile. A zero E-contact profile cannot make the full star indefinite while either marked replacement is finite; an affine principal star is excluded by compactness. These observations cover arbitrary orthogonal spherical spectators, without asserting that spectators can never rescue an isolated profile.

## The 23 quadrilateral recipes

For each literal recipe the checker reconstructs a closed cycle through the saved dotted edge, using the explicit old E-tree path and the new finite contacts. Additional chords do not invalidate a cyclic-product test. Magnitudes can replace signed products here because field membership and algebraic integrality are unchanged by sign.

The saved dotted magnitude is \(a+b\sqrt D\), with \(a,b,D\in\mathbb Q(\sqrt2,\sqrt5)\). The independent checker solves the square-root membership equations over this biquadratic field. A nonzero radical coefficient with nonsquare \(D\) cannot belong to its quadratic subfield. Otherwise exact reduction tests the \(\sqrt2\) coefficient and then the ring \(\mathbb Z[(1+\sqrt5)/2]\). This uses no producer field-membership routine.

The three field-valued failures give primitive irreducible integer polynomials

\[
9X^2-84X+16,\qquad 5X^2-25X+11,\qquad 4X^2-36X+1,
\]

for E7 marks 03, E7 marks 56 and E8 marks 03 respectively. None is monic. All other 20 tested cyclic products lie outside \(\mathbb Q(\sqrt5)\). The earlier geometric realizability of these recipes is inherited from its frozen reviewed package; the present checker independently replays their arithmetic obstructions, not the earlier full face geometry.

## Validation and scope

Generation and fresh saved replay pass:

```sh
python3 scripts/fresh_arithmetic_quadratic_bound_20261015_review.py
python3 scripts/fresh_arithmetic_quadratic_bound_20261015_review.py --verify
```

The result is a complete necessary triangle marked-pair restriction for the literal E7/E8 components in a compact arithmetic Coxeter polytope with ground field \(\mathbb Q(\sqrt5)\), including orthogonal spectators, and a bounded rejection of 23 specific quadrilateral recipes. No simultaneous star completion, exhaustive quadrilateral census, global chamber, improved dimension bound, or priority claim is established.

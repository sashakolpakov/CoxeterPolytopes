# Rational caps cannot complete an entire quartic cusp family

For **local direction record0 and every real shift `u>0`**, no finite collection of projections of rational ambient normals can complete the retained 22 cusp walls to a finite-volume Coxeter H20 chamber. The ambient norms and rational coefficients are unrestricted. This strengthens the [integral norm-two/norm-four gate](fresh_h20_quartic_orbit_cap_gate_20261014.md) from one shift to the full positive real pencil; it is not a statement about every quartic spatial direction or arbitrary H20 constructions.

The [exact producer](../scripts/fresh_h20_quartic_rational_cap_gate_20261014.py) and [certificate](fresh_h20_quartic_rational_cap_gate_20261014.json) use quartic field arithmetic and rational polynomial gcds. They do not enumerate a bounded mirror bank. The original-chamber finite-angle theorem is not used: these caps may replace all original nonincident inequalities, and need not be ambient reflection roots.

The [independent arithmetic and proof review](fresh_h20_quartic_rational_cap_gate_20261014_review.md) passes both generation and saved replay. It reconstructs all 21 rationality gcds and twelve records using a separate algebraic-field implementation, checks the quadratic-subfield exclusion by the reverse field-degree test, and verifies the three identities forcing a rational shift.

## Necessary data at an endpoint

Use `u=1/2` as an algebraic anchor. Its normal is `e=e0−LD/2`, with `e0=Σ k_i r_i`, `e²=q>0`, and the three retained triangle walls indexed by `23,50,200`. Write `a_i=e·r_i` and `n_i=q(r_i²)−a_i²`. The projected triangle is affine A2, with all three labels3, and its positive kernel satisfies `k_i² n_i=k_j² n_j`.

For any rational ambient cap vector l, put

```text
T=e·l,   g_i=r_i·l,   h=D·l,   d=l².
```

Here `g_i,h,d` are rational and T belongs to the quartic field K. A spacelike projected cap requires `q d−T²>0`, so `d>0`. A cap ending an outgoing cusp edge has `h<0` in the outward orientation. In particular `h≠0` is a necessary condition independent of orientation. The normal identity gives

```text
T = k_23 g_23 + k_50 g_50 + k_200 g_200 − L h/2.
```

Consider either two triangle walls contained in that edge. At a finite or ideal opposite endpoint their projected Gram matrix together with the cap is positive semidefinite. With the old triangle pair already having label3, the only possible cap-contact patterns are `(2,m)` or `(m,2)` for `m∈{2,3,4,5,6}`, or `(3,3)`. This follows directly from the normalized three-by-three determinant and includes the ideal cases. A parallel pair cannot occur in this block. Thus no assumption about a larger angle alphabet is necessary.

Label5 would require `cos²(pi/5)=(3+sqrt5)/8` to belong to K. Exact factorization of `x²−5` over the irreducible quartic field leaves it irreducible, excluding that label here.

## Removing the ambient norm bound

If the first contact is right, then `T=q g_i/a_i`; every `a_i` is nonzero. Rational projective scaling reduces `g_i` to 0 or1, with both orientations retained. If `g_i=0`, then `T=0`; the other pairing is either0 or can be normalized to1. These cases are checked separately, including the family with both contacts right and arbitrary positive norm.

For `g_i=1` and other pairing `g_j=r∈Q`, a nonright contact with squared cosine c forces

```text
d = T²/q + (q r−a_j T)²/(c n_j q).
```

For `c=1/4,1/2,3/4`, this is a quadratic polynomial in r over K. A rational cap requires d rational. In the basis `1,alpha,alpha²,alpha³`, its three nonconstant coefficients are rational polynomials of degree at most two. The checker takes their gcd and exhausts every rational root. No case has all three coefficient polynomials zero, so there is no omitted continuous branch. The second right-contact case is handled directly rather than dividing by c=0.

For the `(3,3)` pattern, the two projected negative pairings have the prescribed length ratio, giving

```text
T = q (k_i g_i−k_j g_j)/(k_i a_i−k_j a_j).
```

The denominator is checked nonzero for all three pairs. Normalize the first nonzero rational pairing to1 and again impose rationality of the resulting quadratic expression for d. Both-zero pairings force T=0 and zero projected norm, so are impossible. This covers the ideal affine A2 endpoint pattern as well as the right-contact finite or ideal patterns.

Finally, the four rational coefficients of the normal identity form a linear system for the unused triangle pairing and h. Rank-deficient systems are retained with their rational free parameters. A system forcing `h=0`, or having no solution, cannot supply a closing cap. The implementation deliberately retains possible orientations and conditions that may fail later; this enlarges the necessary candidate set and is safe for exclusion.

## The anchor result

The 21 rationality gcd calculations yield twelve recorded necessary cases. Only three records permit nonzero h, and all concern the pair `{23,50}`; two of those are duplicate orientations of the both-right family. **No rational cap can end an edge whose retained triangle pair is `{23,200}` or `{50,200}`.** Individual compatible caps are not all excluded.

## From the anchor to every positive real shift

Now set `e(u)=e0−uLD` for arbitrary real `u>0`. The values q, a_i and every seed Gram entry stay unchanged, since D is null and orthogonal to every incident root. The positive triangle relation becomes `Σ k_i rbar_i=uLD`; the other two positive relations are unchanged. The same 22 walls therefore retain the bounded full-dimensional cusp link throughout the positive real pencil. Its 252 edge germs are still grouped into three classes of84 by the retained triangle pair.

The angle census also remains complete before assuming u algebraic. At an endpoint, the three-root determinant forces either a right contact or the `(3,3)` pattern. The right-contact formula forces `T∈K`; the `(3,3)` linear formula does the same, with every denominator checked nonzero. Thus the field and rational-norm arguments above apply for every real u, including the exclusion of label5.

For the pair `{23,50}`, the entire angle census, **before** imposing the anchor identity, contains just these two projective types:

- Both contacts right: `g23=g50=T=0`. Since `k200=L`, the normal identity gives `u*h=g200`.
- Both labels3: normalize `g23=0,g50=1`. The checked identity `k50−T=2L` gives `u*h=g200+2`; without normalization this is `u*h=g200+2*g50`.

All quantities on the right and h are rational, and an ending cap requires `h≠0`. Thus **any cap closing an edge of the `{23,50}` class forces u rational**. For irrational u, those84 edge germs are unclosable.

For nonzero rational u, set `h*=2u*h`. This is still a nonzero rational number. The identity `T=Σ k_i g_i−uLh` becomes exactly the anchor identity `T=Σ k_i g_i−Lh*/2`. The anchor's rational linear systems already permit every rational h*, and exclude the other two triangle pairs. Hence for rational `u>0`, those168 edge germs are unclosable.

These two cases exhaust the positive real pencil. The certificate records all three literal `{23,50}` angle records, including the duplicated both-right orientation, and independently checks the rational constants in `u*h=g200+constant`.

## Finite-volume consequence and scope

The [integral cap report](fresh_h20_quartic_orbit_cap_gate_20261014.md) supplies the exact seed geometry and the positive-kernel argument that an additional spacelike acute facet cannot have zero cusp pairing. The rational certificate binds its frozen certificate and checks its source hashes. That geometric argument remains valid for every positive u: all added facets are strict near D, so all252 edge germs persist.

A finite-volume Coxeter chamber would require every such edge to end at a finite or ideal vertex, with some added facet supplying the opposite endpoint. The necessary three-wall block forbids that cap for84 germs when u is irrational and168 when u is rational. The seed alone cannot provide a full ideal endpoint link. Therefore no allowed completion has finite volume at any positive real shift.

This excludes projections of rational ambient cap directions, even if their norms are arbitrary and they are not ambient reflection roots. It does not exclude genuinely nonrational ambient caps, deletion or alteration of the retained cusp walls, other quartic spatial directions, or other constructions. It does not assert that projections of rational ambient vectors always generate an arithmetic group.

```sh
python3 scripts/fresh_h20_quartic_rational_cap_gate_20261014.py --verify
```

# Deleting a wall before lifting two H8 normals

No new finite-volume H9 example results. This construction changes the seed, rather than selecting more caps on the preceding single-endpoint lift. It gives a valid ten-wall Coxeter cone in H9 and **117 distinct additional mirrors**. Seventy-one tested compatible configurations have 40–51 facets and exact infinite-volume certificates. More strongly, one strict spacelike witness survives all 117 caps simultaneously, excluding finite-volume completion from this entire finite dictionary. The unrestricted seed remains open.

## The changed construction

Use the eleven norm-two normals of the archived compact [Bugaenko H8 chamber](../data/known_examples/bugaenko_compact_h8.json), with its original zero-based numbering. Put φ=(1+√5)/2. Delete wall 3 **before** lifting into one additional positive coordinate. Give walls 0 and 9 opposite heights with squares

\[
h_0^2=2\phi,\qquad h_9^2=\frac{2+8\phi}{11}.
\]

Every other retained height is zero. Normalize each lifted normal back to square two. Exactly these pairings change:

| Original pair | New normalized pairing | Relation |
|---|---:|---|
| 0–1 | −1 | π/3 |
| 0–9 | −1 | π/3 |
| 9–10 | −√(3+2φ) | strictly disjoint |

The old 0–1 relation was π/5; 0–9 was orthogonal. The old 9–10 pairing was −2φ. Deletion removes 2–3, 3–4 and 3–9; in particular wall 9 loses its old label-three neighbor. No further retained pairing is changed by the lift. These facts are checked by the exact squared-cosine identities, including the signs of the two heights.

Write q=3+2φ. Rescale the four normalized normals with original indices 0,1,2,9 by √q. All Gram entries then lie in Q(φ): these four diagonal entries are 2q, their nonzero internal pairings are −q, and the 9–10 pairing is −q. The other six diagonal entries are two. The resulting ten-by-ten Gram matrix has exact signature (9,1). Its weighted Gram differs from the previous released single-endpoint seeds already at the squared dotted pairing.

This avoids the earlier [all-eleven-directions lift obstruction](fresh_h8_lift_h9_20261013.md): that theorem retained original wall 3, while this construction deletes it first and changes the height of wall 9 as well as wall 0. It is not a claim to lift the complete compact H8 polytope as a facet.

The other embedding of Q(φ) also gives signature (9,1). In particular this Gram is not an admissible arithmetic quadratic host with positive-definite conjugate. The experiment does not assume arithmeticity and makes no new commensurability claim.

## Exact caps and the remaining gap

Let G be the rescaled seed Gram and I=G⁻¹. A cap with seed pairings p has seed-basis coefficients c=Ip and square pᵀIp. The first bounded dictionary tests squares 2 and 2q, at most five nonzero pairings, and small Q(φ) pairing magnitudes adapted to the two seed norm types. Numerical norm matches are only proposals; exact arithmetic checks all accepted norms, angles, and a common strict timelike center.

This gives 62 caps. Thirty-one deterministic maximal-clique proposals give compatible 40–47-facet configurations, all escaping. In fact a common exact witness survives all 62 caps, independently of their mutual compatibility.

The second construction makes one seed pairing free, keeps at most four other pairings nonzero in the declared dictionary, and solves the norm equation

\[
I_{uu}x^2+2(If)_u x+f^TIf-N=0,
\qquad N\in\{2,2q\}.
\]

The free pairing must be strictly dotted. Up to five numerically scored proposals per free index and norm type are replayed exactly. The 55 accepted caps all cut the first common witness. Thus none is a mirror from the initial dictionary. Exact normalized pairing comparisons also find no aliases among these 55. The initial 62 are distinct: equal-norm coordinates are distinct, and the two norm classes cannot define the same mirror over Q(φ), since Norm(q)=11 is nonsquare.

Forty further compatible configurations containing the new caps have exact spacelike escapes. They use up to 51 facets. Every selected pair is verified using its actual algebraic pairing and an exact finite-angle polynomial; no volume conclusion comes from the preliminary numerical compatibility graph.

The final common witness, archived as `all_dictionary_escape`, has rational coordinates y with

\[
y_i>0,\qquad y^TIy>0,\qquad y^TIh<0,
\]

where h is the fixed positive timelike center in seed-pairing coordinates. Moreover c·y>0 for **every one of the 117 caps**. These strict inequalities place an open patch of future ideal directions inside their intersection. Hence every subset of this finite cap dictionary has infinite volume, whether or not it passes the Coxeter-pair gate. No bound on future caps or impossibility theorem for the seed is asserted.

For any tested pairwise acute configuration, all listed walls are actual facets: projecting the strict timelike center orthogonally to one wall keeps its norm negative and all other wall inequalities strict. The retained seed has full rank ten.

## Reproduction and independent replay

The [seed and bounded-cap producer](../scripts/fresh_h8_twoheight_lift_20261013.py) and [free-distance producer](../scripts/fresh_h8_twoheight_free_caps_20261013.py) archive their candidates and exact witnesses. The separate [standalone replay](../scripts/fresh_h8_twoheight_verify_20261013.py) imports no producer, floating-point library, optimizer, or SAT solver. It uses earlier independently implemented Fraction arithmetic, with separate reconstruction logic.

It rebuilds the Gram from the original H8 data, verifies both signatures, every accepted cap and selected pair, all 71 escapes, both common witnesses, and distinction of the new mirrors. Mutations of a seed pairing and substitution of a timelike escape are rejected. Another agent independently reviewed the lift and reran the final verifier successfully. The [audit JSON](fresh_h8_twoheight_independent_audit_20261013.json) records source hashes and exact counts.

```sh
PYTHONPATH=. python3 scripts/fresh_h8_twoheight_lift_20261013.py --support 5
PYTHONPATH=. python3 scripts/fresh_h8_twoheight_free_caps_20261013.py
PYTHONPATH=. python3 scripts/fresh_h8_twoheight_verify_20261013.py --verify
```

The two principal archives are [the 62-cap seed bank](../data/enumerations/fresh_h8_twoheight_lift_20261013.json) and [the 55 targeted caps and common escape](../data/enumerations/fresh_h8_twoheight_free_caps_20261013.json).

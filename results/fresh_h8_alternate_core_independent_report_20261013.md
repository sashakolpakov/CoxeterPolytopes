# Independent finite-volume obstruction for the 493-cap H9 bank

**No selection from the saved 493-cap bank can complete the ten fixed seed walls to a finite-volume Coxeter polyhedron.** This is a finite-dictionary theorem. The numerical cap proposal process was not exhaustive, so unrestricted H9 completions remain open.

The bank consists of 348 golden-field caps, one polar cap, and 144 coupled caps over quadratic extensions of the golden field. All use the same ten independent seed normals and contain a common strict timelike interior. The [independent cover certificate](fresh_h8_alternate_core_cover_independent_20261013.json) records five exact spacelike witnesses, every one of their cutters in the full bank, and a small Boolean proof that no compatible selection cuts them all.

## The positive construction retained by the audit

The first new coupled cap is an exactly realized additional mirror, even though the bank does not close. With `φ=(1+√5)/2`, write

```
δ=49−24φ,   N(δ)=649,
x=−11+7φ+(φ−1)√δ.
```

It satisfies `(φ/2)x²+(−7+4φ)x−16+8φ=0` and `x>2`. Its only nonzero pairings with old seed walls are `−x` at original wall 0 and `−(2+φ)` at original wall 9 (seed indices 0 and 8). Both are dotted pairs. The cap has norm two and pairs by −1 with the polar of norm `(6−2φ)/5`, giving dihedral angle π/10. The norm 649 proves δ is not a square in Q(φ), so the coefficient field has degree four over Q.

The [first-cap replay](fresh_h8_alternate_first_cap_independent_audit_20261013.json) checks all 66 wall pairs of this twelve-wall seed, its Gram rank ten and inertia `(9,1,2)`, a common strict timelike interior, and its cut of the previous escape. It also gives a strict relative-interior point on every facet. These are obtained by projecting the common timelike interior onto each wall; all other pairings stay strictly negative because wall normals are pairwise nonpositive. Thus this is a real twelve-facet Coxeter polyhedron in H9. No finite-volume or compactness assertion follows from this seed certificate.

## Complete angle audit on the necessary core

Five witness clauses, numbered `[1,21,23,30,40]` in the archived clause concatenation, suffice. They have respectively 29, 12, 12, 13 and 14 cutters. Their union has only **54 nodes**, so all 1431 pairings among those nodes were independently classified. Any compatible selection hitting the five clauses can be reduced to nodes in this union; vertices elsewhere in the bank cannot rescue an incompatible core selection.

A pairing belongs to `Q(φ,√d,√e)`, of degree 2, 4 or 8 over Q. Therefore a finite Coxeter label m must satisfy `φ_Euler(m)/2 | 8`. The elementary prime-factor inequality `φ_Euler(m)²/m≥1/2` bounds m by 512 when `φ_Euler(m)≤16`. Exhausting that range gives the complete safe pool

```
2,3,4,5,6,8,10,12,15,16,17,20,24,30,32,34,40,48,60.
```

The audit constructs the polynomial of `4cos²(π/m)` from the exact cyclotomic resultant, and isolates its **largest real root** by rational intervals. Exact compositum arithmetic then tests membership. Isolating the largest root prevents a forbidden Galois-conjugate angle from being accepted merely because it satisfies the same polynomial. No floating-point prefilter participates in these decisions. SymPy supplies exact cyclotomic polynomial algebra and rational root isolation; pairings, extension arithmetic and signs use Fraction arithmetic.

The [complete pair certificate](fresh_h8_alternate_core_angles_independent_20261013.json) gives:

| Pair type | Count |
|---|---:|
| Positive pairing | 782 |
| Non-Coxeter finite angle | 408 |
| Valid finite Coxeter angle | 64 |
| Dotted | 175 |
| Parallel | 2 |

There are **no missing compact-compatible edges** in the saved graph on this core. The two omitted parallel pairs are `(20,114)` and `(21,114)`. Both are added in the finite-volume proof below.

## Complete cutters and exact contradiction

For each of the five witnesses, the second checker independently recomputes its positive norm, strict seed inequalities and negative pairing with the common timelike center. It then recomputes the cutter sign against **all 493 caps**. The resulting complete lists agree with the saved lists; no cutter was excluded because of a preferred angle or intermediate selection.

A finite-volume selection must cut every spacelike witness. Otherwise that witness, perturbed slightly toward the common strict timelike interior, gives a strict spacelike point and an open patch at infinity. The actual Coxeter facets must also be pairwise compatible. Each original seed wall remains a genuine facet under arbitrary selections: its projected interior point has strictly negative pairings with every cap, since every cap pairs nonpositively with that seed wall. Redundant optional caps can be discarded without changing this argument.

A standard-library DPLL search finds the five covering clauses incompatible with the complete core graph. The certificate stores all branches of its **17-state proof**, and replays every branch independently of the search. It remains UNSAT after **both parallel edges are included**. This establishes a finite-volume obstruction, rather than only a compactness obstruction, for the entire saved bank and all its selections.

```sh
PYTHONPATH=. python3 scripts/fresh_h8_alternate_first_cap_review_20261013.py
PYTHONPATH=. python3 scripts/fresh_h8_alternate_core_angles_independent_20261013.py --verify
PYTHONPATH=. python3 scripts/fresh_h8_alternate_core_cover_independent_20261013.py --verify
```

The [cap-producing agent's separate Fraction audit](../data/enumerations/fresh_h8_alternate_verify_20261013.json) independently verifies all cap certificates, the full set of 41 cutter clauses, and the parallel-augmented five-clause contradiction under a different branch order. The two checks agree. This does not prove nonreflectivity of an ambient arithmetic lattice or exclude caps outside the finite bank.

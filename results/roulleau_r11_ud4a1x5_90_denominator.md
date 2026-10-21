# Exact closure of Roulleau's 90-wall rank-11 stock

Roulleau's published `U+D4(-1)+A1(-1)^5` nef chamber is now an exact
machine stock.  It is a known finite-volume noncompact Coxeter polytope in
`H10`, not a new abstract example, and it cannot produce any nontrivial fold
in dimension at least 10.

The source formulas have a rigid combinatorial interpretation.  The first 20
walls contain eight matched pairs.  The other 70 walls are canonically all
4-subsets of an 8-set: `A_J` is the subset `J` not containing the eighth
point, while `B_J` is its four-element complement containing that point.
Permuting the eight points and complementing every 4-subset gives
`S8 x C2`.  Conversely, weighted degree profiles fix `A2` and `A4`, allow
only the common `A1/A3` layer swap, and reduce the 16-wall cell to eight
intrinsic matched pairs.  Every one of the 70 remaining walls has a unique
incidence signature against that cell.  Therefore the weighted-wall group has
order at most `2*8!`; the explicit integral generators attain the bound:

`Aut_weighted(Gram) = S8 x C2`, of order `80,640`.

The exact coordinate reconstruction verifies:

- the basis `A1,...,A11` has signature `(1,10)`, determinant `128`, and Smith
  invariants `1^4 2^7`, agreeing with `U+D4(-1)+A1(-1)^5`;
- all 90 coordinate rows are distinct integral roots of square `-2`, and the
  full wall Gram has rank 11 and inertia `(1,10,79)`;
- the published coordinates of `A12`, all eight fiber-partner relations, all
  70 closed `A_J,B_J` formulas, and their intersection formulas hold exactly;
- the ample class `D18` has square 18 and the stated wall degrees, while the
  degree-two class has the stated degrees.

There is also a literal Vinberg noncompactness certificate.  The primitive
isotropic fiber `F=A2+A4` is nonnegative on all 90 walls.  Its zero
neighbourhood consists of nine disjoint affine `A1` components,
`(A2,A4)` and `(A5,A13),...,(A12,A20)`, with total parabolic rank 9, the full
cusp rank in `H10`.  Finite covolume is certified separately by the
algebro-geometric arithmetic route: Roulleau's complete 90-wall nef chamber is
a Weyl chamber and the general K3 has finite automorphism group, so its Weyl
group has finite index in the arithmetic orthogonal group.  No claim of a
literal enumeration of every 90-wall elliptic edge is made.

Finally, the walls span the rank-11 Néron--Severi lattice, so the full diagram
action on it is faithful.  Every nontrivial subgroup contains a nonidentity
linear action and hence fixes rank at most 10.  Thus every nontrivial fold is
at most `H9`; the entire fold denominator relevant to dimensions above 9 is
closed without enumerating the subgroup lattice of a group of order 80,640.
Roulleau's `Aut(X)=C2 x C2` and tag `tn` are not silently identified with the
much larger weighted-diagram group.

Replayable files:

- [`roulleau_r11_ud4a1x5_90_denominator.json`](../data/enumerations/roulleau_r11_ud4a1x5_90_denominator.json)
- [`certify_roulleau_r11_ud4a1x5_90_denominator.py`](../scripts/certify_roulleau_r11_ud4a1x5_90_denominator.py)
- [`test_roulleau_r11_ud4a1x5_90_denominator.py`](../tests/test_roulleau_r11_ud4a1x5_90_denominator.py)

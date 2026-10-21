# Scoped review of the 76-wall finite-edge bridges

**PASS.** Source/proof review and a fresh [primary](../scripts/fresh_h20_contact_closure_20261014.py) saved replay confirm the scope of the [76-wall certificate](fresh_h20_contact_closure_20261014.json): sixteen specified actual finite vertices, eight actual finite edges, and an exact strict spacelike escape. The old sixty walls are retained; the previous optional thirty-two caps are replaced.

The [separate Fraction checker](../scripts/fresh_h20_contact_closure_20261014_review.py) imports no producer. It independently rebuilds the PD18 support's three-dimensional rational nullspace and its Gram, then verifies the literal norm-two cap, norm-418 orthogonal mate and anchor pairings −1 and −37. It also checks the sixteen distinct incident sets, eight distinct edge supports, and exact coverage of the eight claimed positive-dual support records. Its [hash-bound result](fresh_h20_contact_closure_20261014_review.json), generation and saved replay pass. This small separate arithmetic check does not duplicate all 1,080 new pair gates.

The primary replay passes 296 right, eight π/3 and 776 strictly dotted new pairs. The common timelike center is strict for all walls; the explicitly checked projections to individual walls establish actual facets. At each proposed vertex the physical rank-20 Gram has positive elimination pivots, its projected point is future timelike, and all nonincident inequalities are strict. Thus these are actual simple finite vertices, not merely positive subdiagrams.

Each stated edge has two distinct such endpoints whose incident sets intersect in precisely the same nineteen positive-definite walls. Their future cone segment satisfies every retained inequality; all other walls are strict in its interior. The two distinct remaining facet equalities terminate opposite directions of that one-dimensional face. This proves a finite edge and avoids inferring closure from a finite anchor angle alone. No complete vertex or edge census is claimed here; the separate frontier extension is outside this review.

The height flips permute the old walls, and all transported normals, norms, pairings, facet points and vertices are tested directly. All sixteen new walls are strict at the three previously repaired cusp rays. Their old bounded links therefore persist. The old four escape images are cut, while the previous 92-wall strict escape remains strictly feasible for the new 76-wall chamber and retains positive norm. The volume is consequently still infinite.

The positive-dual proof is correctly restricted to the eight **height-19-unreflected** support instances. There the escape functional, restricted to the PD18 complement, is a positive combination of three retained old-wall functionals. Every acute perpendicular normal therefore pairs nonpositively with it. Equality would force all three old-wall pairings to vanish; together with the eighteen support equalities, the checked full rank forces the normal to be zero. Every nonzero such normal is thus strictly negative on the escape. Neither the source nor this review extends that conclusion to the reflected support instances, where an acuteness-only counterexample exists.

The inherited old signature, old-old Coxeter gates, bounded-link proofs and conditional non-quasi-arithmetic criterion are not independently re-enumerated here. The finite-edge construction does not establish finite volume, a dimension record, or nonexistence of further completion.

```sh
python scripts/fresh_h20_contact_closure_20261014.py --verify
python scripts/fresh_h20_contact_closure_20261014_review.py --verify
```

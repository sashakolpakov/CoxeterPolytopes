# The A5 source can receive at most one H4 restoration

The previously proved H4 transfer can be strengthened: **all four source
capacities are one**. Only the A5+2E8 entry changes, from W+2t to W+t.
All frozen earlier packets remain unchanged and valid with their conservative
capacity two.

Number the source A5 path 0--1--2--3--4. The complete earlier H4 triangle domain
shows that its potential restored angles are (0,3) and (1,4). A restoration at
(0,3) requires the actual neighboring H4 wall x, whose only contact with these
five roots is a label-five bond to root 0. The symmetric restoration requires
y, with a label-five bond only to root 4. Both walls are orthogonal to the two
whole E8 spectators. These contact statements follow directly by reversing
the deficient H4 pattern with contacts label three at H4 root 3 and at its
whole A1 spectator. The original omitted H4 root has zero contact with the new
third wall, so no contact with source root 3 (respectively 1) is missing.

Normalize the source roots and x,y to square two. For the A5 Gram matrix J,
its inverse has endpoint diagonals 5/6 and opposite endpoint entry 1/6. Write
x=p_x+sqrt(H)e and y=p_y+sqrt(H)e, where e is the future unit timelike vector
of the source vertex. The 21 independent finite-vertex normals span the whole
positive spatial hyperplane; its orthogonal complement is one-dimensional
and timelike. Thus no extra positive projection direction was omitted.
Halfspaces containing this finite vertex strictly have
the same positive height sign. With phi=(1+sqrt(5))/2,

    H = 5 phi^2/6 - 2 = (-9+5 sqrt(5))/12 > 0,
    p_x.p_y = phi^2/6.

Consequently the two actual wall normals would have product

    x.y = phi^2/6 - H = (3-sqrt(5))/3 > 0.

This contradicts acuteness. The two required neighboring walls cannot coexist,
and therefore the two restored source angles cannot both occur. This is a
restriction on actual source incidence, stronger than counting two abstract
angle orbits. It uses no additional arithmetic integrality assumption beyond
the inherited complete H4 domain. All source-angle restoration and small-face
proofs in that packet remain valid; only the A5 vertex allowance decreases.

The constant-size checker verifies the inverse entries, heights and strictly
positive forced product exactly. Generation and a fresh replay can be run as
separate tiny batches:

```sh
python3 -B scripts/fresh_compact21_a5_h4_source_capacity_20261018.py
python3 -B scripts/fresh_compact21_a5_h4_source_capacity_20261018.py --verify
```

Generation passed in chunk `317deb`, and fresh separate-process replay passed
in `5d7879`; both were tiny shared-lock batches and left no temporary output.
The displayed algebra is exact. No Lean proof,
full-model feasible weight vector, or new dimension bound is claimed.

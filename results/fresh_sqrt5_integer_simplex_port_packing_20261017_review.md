# Independent source review of integer simplex-port packing

**PASS; FINAL/FROZEN.** This is an independent mathematical and source
review of the [primary package](fresh_sqrt5_integer_simplex_port_packing_20261017.md).
It is not a second arithmetic replay of its 62,353 new columns. The primary
generation and fresh saved replay passed separately. No correction was
identified.

The reviewed immutable files are:

| File | SHA-256 |
|---|---|
| [Checker](../scripts/fresh_sqrt5_integer_simplex_port_packing_20261017.py) | `4d08168bf847420e48c1ba3c3e9f653f9ce1f5443ec7d03c4311f937d1dd9686` |
| [Certificate](fresh_sqrt5_integer_simplex_port_packing_20261017.json) | `5d25a38ff2f57d50b6081642c8e8996065ac336190e679236608599449d72c75` |

## Incidence proof

The [compact disjoint-port theorem](fresh_simplex_disjoint_ports_20261017.md)
provides pairwise disjoint port blocks C for the distinct maximal simplex
faces at a vertex. For every port subset U and integer t at least two,

\[
 \sum_C\left\lfloor {|C\cap U|\over t}\right\rfloor
 \le\left\lfloor {\sum_C|C\cap U|\over t}\right\rfloor
 \le\left\lfloor {|U|\over t}\right\rfloor.
\]

The first step is the elementary superadditivity of the floor of a sum;
the second uses disjointness. When U has three elements and t=2, the block
coefficient is exactly its number of contained pairs minus twice its
complete-triple indicator. For intersection sizes 0,1,2,3 this gives
0,0,1,1, respectively. The resulting coefficient is nonnegative. For
the full thirteen-port set, t=2,3,4 gives capacities 6,4,3 and simplex
corner coefficients `floor(d/t)`. The compact simplex-size bound d<=4
makes higher t vacuous. Maximal edges contribute zero to these new charges.

Prism rows have zero in this new charge family. This does not remove any
older prism charge or shared-base inequality. A forced prism base is
assigned to its unique maximal simplex block as in the inherited shared
capacity proof; it is not counted again as every possible subsimplex.

## Source and normalization checks

The source reconstructs literal thirteen-node forests, complete pair and
triple orbits, and their incidence matrix. For an orbit O of triples and
an orbit p of pairs, the raw incidence divided by the size of p is the
number h(O,p) of triples containing a fixed pair. Constancy follows from
the full forest automorphisms, and the program checks integrality and
`sum_p |p| h(O,p) = 3 |O|`. Averaging does not assume an ambient chamber
symmetry.

For a simplex row, the implemented coefficient is its pair-incidence
sum weighted by h, minus twice its complete-triple feature. I checked the
inherited feature definition: on simplex rows it does count every triple
inside the simplex port block. The source checks nonnegativity and the
corner total `C(d,2)(13-d)+C(d,3)`, including the newly positive triangle
contributions. Full-port corner multiplicities are recovered from the
one-port total divided by d, and the complete simplex has d+1 corners.

The vertex correction takes the maximum required mass over the new
columns, while retaining every face mass and actual provenance mass.
The subsequent exact inequalities, including inherited mixed and signed
columns, are visibly included in the checker. This source review does
not independently evaluate their entire exported matrix.

## Scope

The certificate's ratio is greater than 12 and obstructs only its stated
partial weighting model. It neither optimizes all new variables nor
proves an ambient dimension ceiling. Its feature export covers positive
dual-support rows; zero-load rows require reconstruction before a full
primal solve. A restriction of a nonnegative charge to specified face
classes is a valid weaker charge, but must be labelled as such.

Finally, the proof uses disjoint compact simplex port blocks. In a
noncompact chamber, the finite-ended-port theorem allows different blocks
to share ideal-ended ports. The present integer bound cannot be applied
to all such ports without an additional argument. The primary package
keeps the compact hypothesis explicit.

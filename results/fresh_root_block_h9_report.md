# Fresh E7+A2 positive blocks in H9

These searches change the positive root block, rather than raising the
height bound on an E8 host. They do **not** produce a new polytope. Both
completed height-30 cones have independently verified spacelike escapes.
The full reflection groups remain undecided: later mirrors may cut off
these escapes.

## Hosts and exact scope

Work over `O=Z[phi]`, with `phi^2=phi+1`, on

```
L_alpha = <-alpha> + E7 + A2,   alpha in {phi, 2phi}.
```

The E7 diagram has edges `(1,2),(2,3),(3,4),(4,5),(5,6),(3,7)`;
A2 has edge `(8,9)`. All spatial diagonal entries are two, all indicated
off-diagonal entries minus one. The spatial determinant is six, and
`det L_alpha=-6alpha`. The exact signatures are `(9,1)` at the chosen
embedding and `(10,0)` at the conjugate embedding. Thus this is a search
for compact arithmetic H9 chambers if the reflection chamber closes.

These block presentations were absent from the earlier repository search
lanes. This is not a claim that their reflection groups are new. In
particular, the `2phi` host has determinant `-12phi`, whose square class
equals that of the previously screened `<-3phi>+E8+<1>` host. Integral
lattice equivalence and reflection-group equivalence with that host have
not been established or excluded. The exact determinants have different
ideals, so the displayed integral lattices themselves are not isometric;
an overlattice relation or a shared group remains possible.

## The complete height-zero chamber is E7+G2

Using only E7+A2 simple walls would miss a zero-height reflection.
On the A2 coordinates, the vector `(1,-1)` has norm six and integral
reflection coefficients. It extends the Weyl group of A2 to the full
hexagon symmetry group, of Coxeter type G2.

The nine initial walls are `-e_1,...,-e_8,e_8-e_9`. The last two have
norms two and six, with pairing minus three, hence angle `pi/6`.
In Dynkin weights `w=R*y`, their inequalities are the ordinary E7
dominance inequalities and `w_9 >= w_8 >= 0`. The omitted wall `-e_9`
is redundant; its pairing with `e_8-e_9` is positive three, so it must
not be retained as another simple wall.

Completeness follows without a numerical census. For a norm-two vector
`y=a+phi*b` in the positive block,

```
Tr((y,y)) = 2(a+b/2,a+b/2) + (5/2)(b,b).
```

The rational E7+A2 lattice has minimum two. If `b` is nonzero, the right
side is at least five, while a norm-two vector has trace four. Therefore
all norm-two vectors are rational roots. These roots span the block,
so every O-linear isometry of it has rational, hence integral, entries
and belongs to the ordinary integral automorphism group. Its factors
are `W(E7)` and `W(G2)`; the E7 diagram has no nontrivial symmetry.
The regression suite independently closes the simple-root reflection
orbits and counts 126 E7 roots and six A2 roots.

## Complete norm alphabet and enumeration

The spatial discriminant module has exponent six. In both hosts,
`6*B^-1` is integral. If an integral root `r` is primitive and has norm
`q`, reflection integrality gives `q | 2*B*r`; multiplication by
`6*B^-1` and primitivity then give `q | 12`.
The ring is a PID, the primes two and three are inert, and every totally
positive unit is a square. Up to unit rescaling the complete candidate
alphabets are therefore

```
alpha=phi:   1,2,3,4,6,12
alpha=2phi:  2,4,6,12.
```

The second lattice is even, excluding odd norms. This is a candidate
alphabet, not an assertion that each class actually occurs.

For each norm, the search enumerates all time coordinates satisfying the
exact height bound `h=alpha*x^2/q <= 30`, conjugate positivity, and time
reflection integrality. It enumerates dominant Dynkin weights with exact
Schur-complement lower bounds for `2*w^T R^-1*w`, verifies reconstructed
coordinates are integral, and enforces every spatial reflection
coefficient. Previous wall inequalities prune the recursion using the
nonnegative entries of `R^-1`. Each equal-height packet is checked for
pairwise acuteness. Parallel nonprimitive representatives are removed;
the smallest available norm is retained. Checkpoints explicitly distinguish
completed bounds from unfinished searches.

## Results

| alpha | Complete height | Control records | Recursion nodes | Walls | Exact spacelike norm |
|---|---:|---:|---:|---:|---|
| phi | 30 | 24 | 36,480 | 29 | `(44-27phi)/5 > 0` |
| 2phi | 30 | 18 | 26,245 | 30 | `(1514-774phi)/1083 > 0` |

The first wall list has 25 roots of norm two, one of norm six, two of
norm one, and one of norm three. The second has 27 roots of norm two,
one of norm six, and two of norm four. Searches took about seven seconds
each on the working machine. No increased cutoff was attempted.

Each artifact includes an exact strict timelike interior point and a
spacelike vector satisfying all wall inequalities. Numerical linear
programming only proposes a support; exact field arithmetic verifies
the resulting ray. The standalone verifier imports no project arithmetic
and replays the root norms, integral reflection coefficients, chronological
heights, pairwise acuteness, strict interior, and escape inequalities
using Python `Fraction`. A strict interior and a feasible spacelike vector
imply an open set of ideal directions in the cone, hence infinite volume
for this finite prefix. No full extreme-ray census was needed.

This independently checks the saved certificates. It does not repeat the
enumeration to independently prove its exhaustion, and it does not prove
that either full lattice reflection group has infinite covolume.

## Artifacts and reproduction

- [Search](../scripts/fresh_root_block_h9_search.py)
- [Geometry](../scripts/fresh_root_block_h9_geometry.py)
- [Standalone verifier](../scripts/fresh_root_block_h9_verify.py)
- [phi artifact](fresh_root_block_h9_phi_h30.json)
- [2phi artifact](fresh_root_block_h9_2phi_h30.json)
- [Structural and certificate regressions](../tests/test_fresh_root_block_h9.py)

```sh
PYTHONPATH=src:. python3 scripts/fresh_root_block_h9_search.py \
  --alpha 0,1 --height 30 --output results/fresh_root_block_h9_phi_h30.json
PYTHONPATH=src:. python3 scripts/fresh_root_block_h9_geometry.py \
  results/fresh_root_block_h9_phi_h30.json --escape-only
python3 scripts/fresh_root_block_h9_verify.py results/fresh_root_block_h9_phi_h30.json
python3 scripts/fresh_root_block_h9_verify.py results/fresh_root_block_h9_2phi_h30.json
PYTHONPATH=src:. python3 -m pytest -q tests/test_fresh_root_block_h9.py
```

Use `--alpha 0,2` and the second artifact path to reproduce the other host.

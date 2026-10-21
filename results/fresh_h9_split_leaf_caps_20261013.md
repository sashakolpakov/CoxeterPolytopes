# Breaking leaf symmetry before completing an H9 chamber

This search changes the seed itself, rather than retaining the preceding
18-wall symmetric seed. A first cap dictionary produces an exact 22-wall
Coxeter configuration in H9; explicit caps with dotted leaf pairings then
allow **30-wall configurations**. Exact strict spacelike escapes prove
infinite volume for every compatible selection from the resulting 228-cap
bank. Other caps remain available for further construction.

## A new two-block seed

Let `phi^2=phi+1`, `phi>0`, and put

```
a = 1/2 + phi,
b = -a/(12a+1) = -5/22 + phi/11,
epsilon = (1,1,1,1,-1,-1,-1,-1).
```

There are sixteen square-two normals `u_i,v_i`, with

```
(u_i,u_j) = (v_i,v_j) = 2 delta_ij,
(u_i,v_j) = -2 delta_ij - a - b epsilon_i epsilon_j.
```

Every cross-species pair is strictly dotted, because `a+b>2`.
The parameters distinguish two groups of four leaf pairs. This breaks the
old permutation symmetry and increases the leaf span from rank nine to
rank ten without adding center walls.

Set `s_i=u_i+v_i`, `d_i=u_i-v_i`, `S=sum s_i`, and
`E=sum epsilon_i s_i`. The two spans are orthogonal. Their Grams are

```
(s_i,s_j) = -2a - 2b epsilon_i epsilon_j,
(d_i,d_j) = 8 delta_ij + 2a + 2b epsilon_i epsilon_j.
```

The sum span has squares `S^2=-128a`, `E^2=-128b` and six zero modes.
The difference span has six eigenvalues eight and the remaining two
`8+16a`, `8+16b`. Since `a>0` and `-1/2<b<0`, the full seed has inertia
`(9,1,6 zero)`. Moreover `S` is a strict timelike point, pairing `-8a`
with every wall.

These statements describe actual half-spaces. Projection of a strict
interior point onto any square-two wall preserves all other strict
inequalities when the wall pairings are nonpositive; hence every retained
wall is a genuine facet.

## Constructing and exhausting finite leaf caps

Write the negative leaf pairings of a cap as `x_i,y_i`. The six kernel
relations force `x_i+y_i` to be constant within each four-block. Put

```
t = (T_plus+T_minus)/2,  h = (T_plus-T_minus)/2,
Delta_i = x_i-y_i,
D = sum Delta_i, F = sum epsilon_i Delta_i.
```

Because the seed spans the full Lorentz space, these pairings determine the
normal uniquely. Its square is

```
-t^2/(2a) - h^2/(2b) + sum Delta_i^2/8
-a D^2/[32(1+2a)] - b F^2/[32(1+2b)].
```

A cap with a balanced `{0,-1}` pattern on one four-block and all zero
pairings on the other has norm two precisely when `1/a+1/b=-12`.
That is the constructive reason for the chosen parameter `b`.

For leaf angles restricted to `pi/2, pi/3, pi/5`, the pairings belong to
`{0,-1,-phi}`. Each four-block has exactly 51 ordered assignments: three
constant pair types and three unequal pair types with sixteen orientations
each. Exact enumeration of all **2,601** two-block assignments gives
**192** square-two normals. This is complete for all real normals with
these leaf pairings, not merely a bounded coefficient search.

All 192 contain the same strict center. Their mutual pairings lie in
`Q(phi)`, whose finite square-two Coxeter values are exactly `0,-1,-phi`.
The exact compatibility graph has **4,986 edges** and **8,478 maximal
cliques**, of maximum size **six**. One such six-cap choice gives the
22-wall configuration. The census is larger and its compatible selections
are larger than the preceding symmetric seed's 70-cap dictionary.

## Exact failure of every compatible subset of this dictionary

For each sign vector `eta` with two plus signs in each four-block, define

```
Y_eta = (1/(2a)+1/1024) S + sum eta_i d_i.
```

There are 36 such vectors. Their squares are
`64-128a(1/(2a)+1/1024)^2>0`; their pairings with all sixteen seed
walls are strictly negative, and their pairing with `S` is negative.

For every one of the 8,478 maximal compatible cap selections, the exact
checker exhibits a `Y_eta` whose pairings with all selected caps are also
strictly negative. Thus every compatible subset of this **literal
192-cap dictionary** has a strict spacelike escape and infinite volume.
This does not exclude other finite leaf angles or dotted leaf pairings.

The [certificate](../data/enumerations/fresh_h9_split_leaf_caps_20261013_exact.json)
contains the full Gram, all normals, every compatible edge, every maximal
clique and its strict escape. The
[exact constructor and checker](../scripts/fresh_h9_split_leaf_caps_20261013_exact.py)
uses Fraction arithmetic over `Q(phi)` and imports no numerical search
code. It checks all normal equations and pairings directly and enumerates
maximal cliques by a complete pivoted recursion.

```sh
python3 -m scripts.fresh_h9_split_leaf_caps_20261013_exact
python3 -m scripts.fresh_h9_split_leaf_caps_20261013_exact --verify
```

## Broader numerical construction scout

A [separate scout](../scripts/fresh_h9_split_leaf_caps_20261013.py) fixes
`a` at `phi+1/2, 9/4, 5/2, 3, 4` and solves the quadratic norm equations
for `b` from 324 assignment orbits. It finds 434 parameter points satisfying
`-1/2<b<0` and `a+b>2`. These yield 15,062 caps counted across parameters
and 12,346 tested greedy compatible selections. Every tested selection
has a numerical spacelike vertex; none is a finite-volume candidate.

Those numerical counts are not counts of distinct polytopes. Greedy
selection is not a complete clique census, and numerical vertex enumeration
is not an exact exclusion of all 434 parameter points. The exact result
above concerns only the one specified seed and its complete finite-leaf
192-cap dictionary. The [numerical archive](../data/enumerations/fresh_h9_split_leaf_caps_20261013_scout.json)
keeps these evidentiary levels separate.

## Explicit dotted caps leave the finite-leaf restriction

Set `d=phi(1+sqrt(phi))>2`, choosing the positive square root. For each
of the same 36 balanced sign vectors, put

```
r_eta = (d+1) S/(16a) + (d-1) sum eta_i d_i/8.
```

Its negative leaf pairings are `x_i=1` or `d` according to the sign,
and `y_i=1+d-x_i`. Thus every leaf relation is either pi/3 or strictly
dotted. The normal has square two because

```
(d-1)^2 - (d+1)^2/(2a) = 2.
```

These **36 new exact caps** lie outside the finite-leaf dictionary. Each
cuts its corresponding earlier balanced escape. They are not claimed to
exhaust possible dotted-leaf caps.

Combining these caps with the 192 previous caps gives an exact graph with
**228 vertices, 8,460 compatible edges, and 156,990 maximal cliques**.
The maximum compatible selection has **fourteen caps**, yielding a
30-wall configuration. The exact archive identifies a literal such
selection.

The field is `L=Q(phi,sqrt(phi))`. If a pairing has form
`u+v sqrt(phi)` with `v` nonzero in `Q(phi)`, its conjugate under the
negative embedding of phi is nonreal. It therefore cannot equal a finite
Coxeter value `-2cos(pi/m)`, all of whose conjugates are real. This proves
an exact finite-angle rejection for every label, rather than trusting a
finite numerical angle list. Pairings with `v=0` use the earlier golden
field classification; dotted and parallel cases use exact signs.

A constructive search attempts the compatible selections and obtains
**twelve rational-coefficient strict spacelike escapes**. Every maximal
clique has at least one of these escapes. The verification reconstructs all
228 normals, all pair classifications, the entire maximal-clique census,
and the coverage by the twelve escapes. No numerical linear program is
used during verification. The discovery program used linear programming
only to propose rational vectors, then checked every sign exactly.

The [combined-bank certificate](../data/enumerations/fresh_h9_split_leaf_caps_20261013_dotted_exact.json)
and [exact checker](../scripts/fresh_h9_split_leaf_caps_20261013_dotted_exact.py)
therefore exclude finite-volume completion from this fixed 228-cap bank.
They do not exclude other caps or completions of the sixteen-wall seed.

```sh
python3 -m scripts.fresh_h9_split_leaf_caps_20261013_dotted_exact --verify
```

Two independent replays confirm the result. The
[Gram-inverse audit](fresh_h9_split_leaf_caps_independent_20261013.md)
reconstructs the finite-leaf census by solving on a ten-normal basis and
uses different arithmetic and clique enumeration. A separate
[Fraction and NetworkX audit](fresh_h9_split_leaf_caps_20261013_dotted_independent.md)
rechecks all 25,878 pairs in the combined bank, reproduces its clique digest,
and verifies the selected fourteen-cap configuration and all escape coverage.
Both reject corrupted geometric data.

## A larger dotted-cap scout

The [larger numerical construction](../scripts/fresh_h9_split_leaf_caps_20261013_dotted.py)
lets one or both four-blocks have an unknown common dotted pairing,
solving 765 quadratic norm templates. It proposes 2,848 additional normals,
for a total of 3,040 with the finite-leaf bank. A numerical pair graph and
120 seeded greedy searches prioritize cutting the known balanced escapes.
All 120 tested selections, containing 95 to 137 walls including the seed,
still have numerical spacelike vertices. The best tested maximum Klein
vertex square is approximately 4.04636, well outside the unit ball.

These scout results alone are numerical proposals and numerical failures.
The exact 36-cap construction above is an independent algebraic extraction
from this scout. Its broader
[numerical archive](../data/enumerations/fresh_h9_split_leaf_caps_20261013_dotted_scout.json)
remains preserved as the discovery record.

In a [later exact audit](fresh_h9_split_leaf_3040_independent_20261013.md),
the 3,040 profiles were reconstructed and six repeated-root duplicates
identified, leaving 3,034 distinct mirrors. The subsequent
[endpoint-cover proof](fresh_h9_split_leaf_endpoint_cover_20261013.md)
excludes every finite-volume compatible completion from this entire bank.
It uses exact endpoint conditions and a universal spacelike escape,
independently of the numerical failures above. It does not exclude caps
outside the bank.

# Four exact edge-endpoint clauses for the split seed

Every finite-volume acute-compatible polytope obtained by adding finitely many facets to this sixteen-wall seed must meet four necessary conditions: at each end of its two orthogonal eight-wall edges, some added cap must have **nonnegative** pairing. Zero pairing must be included.

For the archived 3,040 cap entries, the exact replay checks all **12,160 endpoint pairings**. Each necessary list has fourteen entries: twelve positive pairings and two zero pairings. These exact lists agree with the numerical scout. The clauses alone are satisfiable; a compatibility and escape argument is separate.

## Exact endpoints

The seed has orthogonal square-two groups `u_0,...,u_7` and `v_0,...,v_7`, with

`(u_i,v_j)=−2 delta_ij−a−b epsilon_i epsilon_j`,

where `a=1/2+phi`, `b=−5/22+phi/11`, and `epsilon=(1,1,1,1,−1,−1,−1,−1)`.

Write `S=sum(u_i+v_i)`, `E=sum epsilon_i(u_i+v_i)`, `U=sum u_i`, and `Uepsilon=sum epsilon_i u_i`. The orthogonal complement of the eight `u_i` is Lorentzian of dimension two. It has perpendicular vectors

`P_U=S+4a U` and `R_U=E+4b Uepsilon`,

with

`P_U²=−128a(1+2a)=−384−640phi<0`,

`R_U²=−128b(1+2b)=(1664−384phi)/121>0`.

Its two future null endpoints are

`P_U ± sqrt(45+74phi) R_U`.

Swapping `u` and `v` gives the other edge. Each endpoint is orthogonal to exactly its designated eight seed normals and has strictly negative pairing with the other eight. Its center pairing equals `P_U²<0`, so its future orientation is explicit.

The rank/signature statement also follows directly from the sum and difference decomposition: the sum block has one negative and one positive eigenvalue and six zero modes; the difference block is positive definite because `a>0` and `−1/2<b<0`. Thus the full seed has inertia `(9,1,6 zero)` and these rank-eight orthogonal complements have signature `(1,1)`.

## Why nonnegative pairing is necessary

Let `z=P_U+sqrt(d) R_U` be one endpoint. If every selected cap has strictly negative pairing with `z`, move slightly outward along

`z_eta=P_U+(sqrt(d)+eta) R_U`, with `eta>0`.

Its norm is `R_U²(2 eta sqrt(d)+eta²)>0`. The eight `u_i` pairings remain zero, and every other seed and selected-cap pairing remains strictly negative for sufficiently small `eta`. Perturb slightly toward a common strict timelike interior to make the eight zero inequalities strict as well, while retaining positive norm and future orientation. This gives a strict spacelike escape and proves infinite volume.

The other three endpoints have the identical argument, using the corresponding outward sign. Consequently, finite volume requires at least one cap with pairing **greater than or equal to zero** at every endpoint. If a cap has zero pairing, this argument no longer permits the outward perturbation while preserving that inequality; such a cap may help form an ideal endpoint. Omitting the zero cases would therefore be an unjustified stronger condition.

These are actual edges in any full-dimensional acute-compatible completion. Use that completion's own strict timelike interior `p` and project it orthogonally onto its eight retained orthogonal walls. Its norm becomes more negative. Its pairing with every other wall becomes strictly more negative or stays strict because all normal pairings are nonpositive. This leaves an open segment in the Lorentzian two-dimensional complement.

No cap can have nonnegative pairing with both ends of the same edge: that projected strict interior is a positive linear combination of the two future endpoint rays, and its pairing with the cap is strictly negative. A positive pairing at one end is possible; the literal bank contains twelve such caps at each endpoint. This excludes only simultaneous nonnegative pairings at both ends of one edge.

## Exact lists

Indices are zero-based entries in the authoritative exact-normal archive, which has 3,040 entries representing 3,034 distinct mirrors.

| Edge end | Nonnegative cap entries | Zero entries |
|---|---|---|
| `U−` | 6,7,8,9,10,11,195,218,219,220,221,231,263,1767 | 263,1767 |
| `U+` | 0,1,2,3,4,5,194,214,215,216,217,230,262,1766 | 262,1766 |
| `V−` | 6,7,8,9,10,11,193,197,210,211,212,213,233,297 | 233,297 |
| `V+` | 0,1,2,3,4,5,192,196,206,207,208,209,232,296 | 232,296 |

The `U±` convention is reversed relative to the numerical eigensolver's choice of spatial basis; the actual clauses agree.

The [exact archive](../data/enumerations/fresh_h9_split_leaf_edge_endpoints_20261013.json) includes every endpoint coefficient, seed pairing, all cap signs, and both source hashes. The [standard-library replay](../scripts/fresh_h9_split_leaf_edge_endpoints_20261013.py) imports no search or numerical producer. It independently checks every cap's norm, common strict center, nonpositive seed pairings, and literal normal-to-covector identity.

For cap field `Q(phi,sqrt(Delta))`, pairing with an endpoint lies in the compositum with `sqrt(45+74phi)`. The replay computes its sign by two successive exact squared-norm comparisons. It does not assume these radicals are independent; zero and dependent radicals have explicit cancellation checks.

```sh
python scripts/fresh_h9_split_leaf_edge_endpoints_20261013.py --verify
```

This endpoint certificate asserts four necessary clauses for the literal bank. It does not by itself assert a Boolean contradiction or rule out caps outside the bank. The geometric endpoint condition itself applies to arbitrary full-dimensional acute-compatible polytopes with finitely many facets retaining this seed.

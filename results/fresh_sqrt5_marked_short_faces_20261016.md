# A short-face restriction on both maximal marked vertex types

Let `H` count vertices of the two marked baseline38 types on a selected
minimality14 facet of a hypothetical compact21-dimensional arithmetic Coxeter
polytope over `Q(sqrt5)`, and let `L` count its other vertices. The exact local
analysis below gives the new conditional incidence restriction

**`34 H < 155 L`, hence `H/(H+L) < 155/189`.**

The inherited baseline average is consequently less than `7148/189`, about
37.8201. This improves its previous maximum38 but remains above the required
35; it does **not** exclude dimension21.

Artifacts: [source](../scripts/fresh_sqrt5_marked_short_faces_20261016.py) and
[literal certificate](fresh_sqrt5_marked_short_faces_20261016.json).
Inputs are the frozen [complete marked transition test](fresh_sqrt5_marked_maximum_transitions_20261016.md)
and its recursive dependencies. Both maximal marked types are included:
`A1+D20` with the isolated `A1` marked, and `D21` with its long-arm endpoint
marked. The calculation includes all narrow cross transitions between them.

## All narrow parameters, with the new pairing free

Every attachment and every admissible framing product from the complete
minimality14 transition catalogue is used. Across380 port pairs, there are
3,028 pairs of edge choices. The old finite21-root matrix is `A`; the two new
root contact vectors are `r,s`. Eliminating `A` gives

`[[a,x-c],[x-c,d]]`, where `a=2-r^t A^-1 r`,
`d=2-s^t A^-1 s`, and `c=r^t A^-1 s`.

The saved exact calculations give `a,d<0`, `a',d'>0`, and `c>0`. The new-new
product `x` is free. Ambient dimension22 requires
`D(x)=ad-(x-c)^2=0`; a finite fourth vertex requires `-2<x<=0`.
For3,010 choices, `D(-2)>0` and `D(0)>0`, so strict concavity excludes
**every real** finite fourth-corner product. The remaining18 choices are
tested exactly against the complete finite-label catalogue.

For each permitted finite label, `x^2` belongs to `K=Q(sqrt5)`. Since `c!=0`,
the corank equation itself forces `x` into `K`. The products for labels4,6,10
are not in `K`: their squares are2,3,`2+phi`, respectively. Rational-square and
norm checks establish this exactly. Only products0,`-1`,`-phi` remain.
These give54 exact tests, with exactly one finite quadrilateral closure.

The unique closure starts at `A1+D20`, using the last two roots18,19 of the
long arm. Both returns copy the old attachments, both framing products are
`-2phi`, and the new-new product is `-1`. The fourth vertex is again the same
marked `A1+D20` type. The certificate saves all four finite future rays,
their wall inequalities, the root dependence, physical signature `(21,1,1)`
and conjugate signature `(22,0,1)`; the final entry counts the one zero
eigenvalue. It is a genuine local quadrilateral. No quadrilateral closure
exists from the maximal marked `D21` type with both neighboring vertices
maximal.

These statements allow **any finite fourth vertex**, not just a maximal one.
In particular, a quadrilateral with three maximal vertices must have four.

## Triangles and the shared end-corner budget

An all-maximal triangle is impossible. Its two continuations from a maximal
vertex would use the same new wall; its contact at either released root would
have to be both0/`-1` (a simply laced endpoint contact) and `-phi` or dotted
(an admissible maximal-edge framing product). Those possibilities are
disjoint.

A triangle with exactly two maximal vertices must use a maximal-to-maximal
edge whose framing product is `-phi`; a dotted product cannot occur in the
third finite vertex. Exhausting those choices and all possible second
released roots leaves three directed cases:

| Maximal source | Maximal neighbor | Third vertex | Marked baseline at third vertex |
| --- | --- | --- | --- |
| `A1+D20` | `A1+D20` | `A1+D18+I2(5)` | 35 |
| `A1+D20` | `D21` | `D18+H3` | 35 |
| `D21` | `A1+D20` | `D18+H3` | 35 |

All three have exact finite future endpoint rays. At a maximal `A1+D20`
vertex, these triangles use precisely the same end corner `{18,19}` as the
unique quadrilateral. At a maximal `D21` vertex, they use its single end
corner `{17,18}`; there is no maximal quadrilateral there. Thus **each maximal
vertex has only one corner available to either kind of face**.

## Counting on the selected20-dimensional facet

Let `T_2` count triangles with exactly two maximal vertices, and `Q_4` count
all-maximal quadrilaterals. A simple polytope has one2-face for each pair of
incident edges. The shared end-corner restriction gives

`2 T_2 + 4 Q_4 <= H`.

Every other triangle has at least two lower vertices; its contribution2 to
the positive face excess can be charged at most1 per lower-vertex corner.
Every other quadrilateral also has at least two lower vertices and contributes
only1. There are `binom(20,2)=190` corners per lower vertex. Therefore, writing
`T,Q` for all triangles and quadrilaterals,

`2 T + Q <= H + 190 L`.

The independently frozen h-vector calculation gives `f_2>45 f_0` on a
simple20-polytope. Also `sum k N_k=190 f_0`, where `N_k` counts its `k`-gonal
2-faces. With `f_0=H+L`,

`35(H+L) < 5 f_2 - 190(H+L) <= 2T+Q <= H+190L`.

This proves `34H<155L`. The exact marked baseline census gives cost38 on the
two maximal types and at most37 elsewhere, proving the stated strict average
bound `37+155/189=7148/189`.

## Verification and scope

Generation and a fresh saved-output replay both pass using exact `Fraction`
arithmetic in `Q(phi)`, without numerical optimization. Recursive input hashes
are checked. Replay command:

```sh
python3 scripts/fresh_sqrt5_marked_short_faces_20261016.py --verify
```

The local catalogue is exhaustive for the two specified maximal marked
types and their narrow transitions. The incidence conclusion is conditional
on the selected narrow facet of a compact arithmetic21-polytope; it does not
transfer to arbitrary whole-polytope frequencies. The saved triangles and
quadrilateral are local faces, with no finite-volume completion or new
polytope claimed. The remaining gap to35 is explicit. No new LP was run.

# Graph completion removes the two formal-row witnesses

19 October 2026. Paper-only successor to
`right_angled_mixed_eleven_formal_row_barrier_research_20261019.md`.
No script, optimizer, graph census, Lean or build ran. This proves two
explicit graph-dependent inequalities, not a universal graph inequality
or a new dimension bound.

For a graph `G` on eleven labeled directions and an admissible formal
frame `sum p=3`, `sum p^2=1`, write

    F(G,p)=sum_(S a triangle of G) [(2-sum_S p)^(-2)-16/25].

Thus the deletion constant is `c=9/25`, with exponent two. A uniform
bound `F(G,p)<=C*c`, valid for every actual corner graph and every
admissible row, would imply `2B<=C*f` by the same actual reciprocal-edge
pairing as in the reviewed dimension-twelve proof. This observation
does not assert such a uniform bound. The old two-frame obstruction
used arbitrary subsets of triples, not necessarily graph triangles.

## 1. Frame A: an exact completion penalty

Take two high coordinates `6/11` and nine low coordinates `7/33`.
The weights for triples containing two, one and zero high coordinates
are respectively

    w2=9/25+560/529,
    w1=9/25-67/1156,
    w0=9/25-104/225=-23/225.

The first two are positive. Adding an edge incident with a high
vertex creates only these positive triangle types. Therefore every
graph can be enlarged without decreasing `F` to `K2 join H`, where
`H` is a graph on the nine low vertices. If `e,t` are its edges and
triangles, then the exact completed row is

    F=9*w2+2*e*w1-(23/225)*t.                  (1)

For an edge `ij` of `H`, its common-neighbor count is at least
`d_i+d_j-9`. Summing over edges and applying square-sum Cauchy gives

    3t>=sum_i d_i^2-9e>=4e^2/9-9e,
    t>=4e^2/27-3e.                            (2)

Use the elementary rational comparisons

    w2<71/50,    2*w1<303/500,
    303/500+23/75=1369/1500<137/150.

Equations (1)--(2) imply, for every possible real value of `e`,

    F <639/50+(137/150)*e-(92/6075)*e^2
      <=639/50+506763/36800
       =977067/36800
       <666/25=74*(9/25).                     (3)

The middle inequality is the maximum of a concave quadratic; it
does not use a graph census. The final strict comparison has positive
cleared difference `666*36800-977067*25=82125`.

Consequently **every graph at Frame A has `F/c<74`**. In particular,
the scalar selection of all 81 triples involving a high coordinate
cannot be used: its 72 one-high triples force all low--low edges,
and hence also the 84 low triples. The completed-graph argument above
allows arbitrary omissions and proves a stronger quantitative statement
than merely pointing out that inconsistency.

## 2. Frame B: the scalar selection is an actual clique

Take one distinguished coordinate `(3-2sqrt(5))/11` and ten equal
coordinates `(15+sqrt(5))/55`. The respective ratios for triples
containing or avoiding the distinguished vertex are

    rH=(65+8sqrt(5))/55,
    rL=(65-3sqrt(5))/55.

Since `sqrt(5)>2`, `rH>81/55>5/4`, so every distinguished-vertex
triangle has negative weight. Also `1<rL<5/4`, so every other
triangle has positive weight. It follows immediately that `K10`
plus the isolated distinguished vertex maximizes the row:

    F(G,p)<=120*(rL^(-2)-16/25).                (4)

Using `sqrt(5)<9/4` gives `rL>233/220`. Hence

    rL^(-2)<48400/54289<223/250,

where the second comparison has cleared positive difference
`223*54289-250*48400=6447`. Substitution into (4) gives

    F(G,p)<120*(223/250-16/25)
           =120*(63/250)=84*(9/25).            (5)

Thus **every graph at Frame B has `F/c<84`**. Unlike Frame A,
the selected 120 triples are a graph triangle family; their exact
row is nevertheless below the current lower degree threshold.

## 3. What remains

The actual-subset/reciprocal summation is unchanged. What remains is
a universal graph-dependent row estimate (or an adequate global
replacement), at all admissible eleven-dimensional frames, below the
strict incidence threshold `117825/1379`. The two exact inequalities
above remove the old formal witnesses but do not cover other frames,
other graphs, or the compatibility of rows across a whole opposition
component.

For orientation, `K9` with coordinates `1/3` on its nine vertices
and zero on the other two has `F=84*c`. Therefore any such universal
bound has `C>=84`. This is a formal graph-row example, not a claimed
hyperbolic polytope or a finite actual opposition component. It leaves
a possible analytic interval between 84 and `117825/1379`; no bound
in that interval is established here.

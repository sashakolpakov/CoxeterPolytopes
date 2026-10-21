# Exact triangle-density minimum of the current linear relaxation

Date: 2026-10-18. A bounded numerical minimization followed by exact primal
and dual replay gives

    min (2B/f) = 6620/59 = 112.2033898305...

for the stated dimension-twelve aggregate linear relaxation. This is not a
realizability result or a dimension bound. All geometric strict rows were
relaxed to nonstrict rows. The nonlinear Macaulay conditions, actual facet
gluing/cusp widths, and a quantitative closed-cluster defect condition are
not included.

## Linear model and universal polygon cone

Normalize f=1. Variables are c, a6,...,a11, D, 77 nonnegative weights on
complete multipartite twelve-vertex clique vectors, and four nonnegative
minimal-polygon counts. The total is 89 variables, 76 nonnegative rows,
and three equalities. DS determines a0,...,a5; the equalities normalize f,
sum the corner weights to one, and count the non-T polygons.

The clique-vector hull reduction is proved in the separate frozen
opposition-moments packet. T,2B,5C,K are their weighted edge, triangle,
four-clique and five-clique counts. The model includes the previous facet
and facet-plus-cusp minima; the stated NK, weighted density and vertex
ratio rows; canonical P3/P4 facet and sum gaps, the P5 facet gap, canonical
face/finite/ideal/coface caps; K>=16D; 11K<=792f+160D; and the subsequent
universal reciprocal-square cap 2B<=516992f/3125.

Every noncounted hyperbolic polygon dominates one of the minimal
(finite,ideal) pairs (0,3),(2,2),(3,1),(5,0). The four polygon variables
sum to a2-T; their finite and ideal corner sums are at most 66f-T and
11264c-2T. Their two nonnegative residuals allow arbitrary extra corners.
This universally valid cone is used instead of assuming that the four
polygon classes in an earlier particular witness exhaust all polygons.

## Exact certificate

The JSON stores a rational feasible primal assignment. Its nonzero entries
are

    c=9/1298,       a6=248/59,       D=3201/236,
    w_(1,1,1,1,1,1,2,2,2)=662/1121,       w_(12)=459/1121,
    polygon_(3,1)=4140/1121,       polygon_(5,0)=3972/1121.

All a7,...,a11 vanish at this relaxed optimum. The (12) partition is the
empty twelve-vertex graph; the other partition has six singleton parts and
three doubled parts. Each of its complete nine-cliques selects all six
singletons and one vertex from each doubled part. Nothing asserts that
these weighted graphs can be glued as corners of an actual polytope.

The exact dual has only three nonzero nonnegative row multipliers:

| Nonnegative row | Multiplier |
| --- | ---: |
| weighted two-face gap | 200/11 |
| ratio-upper gap (1298/9)c-f | 720/59 |
| P3 facet-plus-cusp gap B-10a3+10a2+28160c | 2 |

After the DS identity a3=5a2-55f-14080c, their sum is exactly

    2B-(6620/59)f.

The saved primal attains equality and the exact replay checks every row,
every nonnegative variable, every dual coefficient/reduced cost, and
equality of primal and dual objectives. This proves the stated linear
minimum without trusting the numerical solver's residuals.

There is also a short direct proof of the universal strict lower bound.
The canonical P3 sum gap and weighted inequality give

    B >=40a2-550f-168960c >50f+880c.

Since c>9f/1298, we obtain 2B>6620f/59. A universal upper cap reaching
6620/59 would therefore exclude dimension twelve by this route; the current
165.43744 cap does not. Additional geometric restrictions can also change
the relaxed minimum, so this is a target for this model, not a necessary
target for every possible proof.

## Further exact consequences

The stronger lower bound and t(G)<=64+10k4(G) imply

    C > (1422/1475)f

in dimension twelve. A K4-free corner graph has at most 64 triangles, and
any twelve-vertex graph at most 220. Thus more than 237/767 of the finite
vertices have a K4 in their corner graph.

In dimension eleven the analogous P3 sum gap is

    B >=10a3-9a2-11520c
       =36a2-(825/2)f-69120c
       >(75/2)f+360c.

The reviewed weighted row a2>(25/2)f+(21230/11)c and c>3f/236 therefore
give 2B>4965f/59. Using the integer K4-free triangle cap 49, the edge
deletion coefficient nine and maximum 165 triangles gives

    C > (2074/2655)f,
    fraction of K4-bearing finite corners >1037/3422.

These are necessary occurrence bounds, not a dimension reduction. The
earlier frozen core packet is unchanged.

## Resource use and reproduction

The successful bounded recovery run used the shared exclusive proof lock,
one foreground optimizer with one thread, a 59-second wall-time/RSS guard,
and resource gates. Peak RSS was 81 MiB and total work took a few seconds.
A preceding diagnostic run failed in exact dual construction; it supplied
no claimed optimum. The coordinating agent explicitly authorized this
single recovery run after immediate snapshot saving and exact free-dual
handling were implemented. No parameter sweep or optimizer farm was run.

Run
`python3 scripts/fresh_right_angled_mixed_triangle_minimum_research_20261018.py`
for optimizer-free exact replay of the final certificate. Fresh replay
c9f826 passed. The optional `--recover` redoes rational reconstruction from
the compact numeric snapshot without optimization. `--optimize` is a new
guarded optimizer run and requires a new task justification; it was not
used for the final replay. No Lean compilation was performed.

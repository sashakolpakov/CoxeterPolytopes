# Intrinsic payment rule and the local-certificate obligation for dimension57

Staging appendix,20 October2026. This text defines the **current ambient58
rule**, rather than a succession of improvements. It is intended for integration
with the geometric and algorithmic proof sections. No computation, import,
replay or build was run while preparing it. The routing assertions below are
separated from the arithmetic definition: defining the weights and reproducing
their sums does not by itself prove the geometric payment theorem.

## 1. Objects on which the rule is evaluated

At a finite vertex let `F` be its full spherical Coxeter forest, with two
distinct marked nodes `a,b`. Components in this appendix always mean **whole
components of F**. Define

    w_F(a,b) = s_T(a,b;F),             if a,b belong to the same component T;
    w_F(a,b) = j_{S,T}(a,b)+k_F(a,b),  otherwise.

Here `s` is Sections2–3, `j` is Section4 and `k` is Section5. An omitted
cross-role has weight zero. An absent context witness contributes zero;
additional witnesses do not multiply a Boolean contribution. Consequently
this is one function of the full marked forest, not a choice of the most
convenient historical local rule.

Number `A_r` by the path `0,...,r-1`, `B_r` by the same path with the last
bond doubled, and `D_r` by spine `0,...,r-3` with leaves `r-2,r-1` at its
last vertex. The finite names `B_r,C_r` describe the same labelled Coxeter
graph here. Use the standard labelled exceptional diagrams. The two long-arm
ends of `E6` are its two leaves on arms of length two. An `A_r` endpoint is
in `E_r={0,r-1}`; its kth-from-an-end set is
`P_k(A_r)={k-1,r-k}`, with repetitions removed. A `D_r` long endpoint is0,
its fork center is `r-3`, and its short leaves form `L_r={r-2,r-1}`.

The complete finite alphabet is `A_r`, `B_r (r>=2)`, `D_r (r>=4)`,
`E6,E7,E8,F4,H3,H4`, and `I2(m),m>=5`. The last family is **not truncated**.
Its ordinary same-component pair has weight one. In cross tables it is a
small rank-two component whenever a row is stated for every small component;
it is not `A2` in rows expressly naming `A2`.

At an ideal vertex apply Section6 separately to each whole irreducible affine
factor of its Euclidean link. The weight is assigned to the varying factor
of the marked two-face. Finite cross-role and context tables are not applied
to a cusp. The factor-fiber counting lemma explains the resulting division
by the order of that affine factor.

All tables below are invariant under the automorphism group of the labelled
component. Thus any canonical identification gives the same number.

## 2. Same-component finite rule in all large ranks

For a finite classical component `T` of rank `r>=18`, let `d` denote its graph
distance. A feature is a vertex of degree other than two, or an endpoint of
a nonsimple bond. Put `delta(v)=distance(v,features)`. Start with weight one
on pairs satisfying

    d(a,b)<=8, or d(a,b)=9 and at least one mark is a leaf of T,

and zero otherwise. Delete a started unit when

    d(a,b) in {6,7,8} and min(delta(a),delta(b))>=6.

For `D_r,r>=20`, change each distance-nine pair containing a short fork leaf
from one to one half. These are the only further large-rank ordinary changes.
There is one contextual `D28` exception in Section5. The ordinary component
sums, a useful independent identification of this definition, are

    b(A_r)=8r-34-(r-18)_+-(r-19)_+-(r-20)_+,
    b(B_r)=8r-34-(r-19)_+-(r-20)_+-(r-21)_+,
    b(D_r)=8r-32-1_{r>=20}-(r-20)_+-(r-21)_+-(r-22)_+.

Here `(x)_+=max(x,0)`. These formulas are sums of marked-pair weights, not
substitutes for their marked locations. The distance-six, seven and eight
deletions use the independent shorter-path routes of Section8; they cannot
be justified by the fact that their total cost has this formula.

## 3. Complete finite tables in ranks at most17

Every pair in every component of rank at most11 has weight one. In ranks
12–17 start with one on every unordered pair and apply exactly the following
ordinary replacements. In a row `Z`, `H`, `Q`, `O` mean respectively weights
zero, one half, one quarter and one eighth. Other explicitly given weights
have their literal value. Components not mentioned retain all their units.
The contextual restorations are in Section5, not in the ordinary rows.

```text
D12 H: (0,10),(0,11).

A13 Z: (0,12).
B13 Z: (0,12).
D13: no ordinary replacement.

A14 Z: (0,12),(1,13),(0,13).
B14 Z: (0,13),(1,13),(2,13),(3,13).
D14 Z: (0,12),(0,13); H: (0,11).

A15 Z: (0,12),(0,14),(2,14),(0,13),(1,14).
B15 Z: (0,12),(0,13),(0,14),(1,13),(1,14),(2,14),(3,14),(4,14).
D15 Z: (0,11),(1,13),(1,14);
    H: (0,12),(0,13),(0,14),(1,12).

A16 Z: (0,12),(0,14),(0,15),(1,14),(1,15),(3,15),(1,12),(3,14);
    H: (0,13),(1,13),(2,14),(2,15);
    O: (0,11),(4,15).
B16 Z: (0,11),(0,12),(0,13),(0,14),(0,15),
       (1,12),(1,13),(1,14),(1,15),(2,13),(2,14),(2,15),
       (3,14),(3,15),(4,15),(5,15).
D16 Z: (0,10),(0,11),(0,12),(0,13),(1,11),(2,12),(2,14),(2,15);
    H: (0,14),(0,15),(1,12),(1,13),(1,14),(1,15),(2,13).

A17 Z: (0,10),(0,11),(0,12),(0,14),(0,15),
       (1,10),(1,11),(1,13),(1,14),(1,16),
       (2,12),(2,14),(2,15),(2,16),(3,15),(4,14),(4,16),
       (5,15),(5,16),(6,15),(6,16);
    H: (0,16),(2,13),(3,14);
    weight7/8: (1,12),(4,15).
D17 Z: (0,10),(0,11),(0,12),(0,13),
       (1,10),(1,11),(1,12),(1,13),(1,14),(1,15),(1,16),
       (2,12),(2,14),(3,15),(3,16),(4,14),(4,15),(4,16);
    H: (0,14),(2,13),(2,15),(2,16),(3,13),(3,14);
    Q: (0,15),(0,16).
```

For `B17` replace the default-one convention by the following positive set:
distance at most eight, distance nine at a whole-component leaf, and the
three additional pairs `(2,11),(3,12),(4,13)`, all with weight one. All
other `B17` pairs have weight zero.

The ordinary sums are

|rank|A|B|D|
|---:|---:|---:|---:|
|12|66|66|65|
|13|77|77|78|
|14|88|87|177/2|
|15|100|97|100|
|16|433/4|104|217/2|
|17|453/4|105|227/2|

These rows also make clear that an envelope representative in the forest
algorithm need not be the actual finite type present geometrically.

## 4. Complete ordinary cross-component marked roles

The tables are symmetric in the two marked components. A Cartesian product
of mark sets means **every** marked pair in that product has the displayed
weight. First apply the broad table, then replace the weights at any roles
listed in the exception table. Exception entries are replacements, not
increments. Unlisted roles have weight zero. There is no conflict among
exception rows having different final values.

In the broad table `S` is any finite irreducible component of rank at most5
and `all(S)` is its whole node set.

|small component and marks|large component and marks|weight|
|---|---|---:|
|`S,all(S)`|`D10,L10`|1/2|
|`S,all(S)`|`D11,D12,D13,D14`, center and short leaves|1/2|
|`S,all(S)`,rankS<=2|`D15`,center and short leaves|1/2|
|`S,all(S)`,3<=rankS<=5|`D15,L15`|1/2|
|`A1`|`D16`,center and short leaves|1/2|
|`S,all(S)`|`A13` middle6 or `A15` middle7|1/2|
|`S,all(S)`|`A17`,marks5,11|1/2|

The remaining ordinary assignments, including replacements in this table,
are as follows. `P_k` always refers to the `A` component named in its column.

|small component and marks|large component and marks|weight|
|---|---|---:|
|`A1` or `A2`,all marks|`D11,L11`|1|
|`A4,E4`|`D11,L11`|1|
|`A5,P2`|`D11,L11`|1|
|`E6`,long ends|`D11,L11`|1|
|`D5`,long end0|`D11,L11`|1|
|`A1`|`D12,L12`|1|
|`D5,L5`|`D12`,center and short leaves|1|
|`A2`,all marks|`D13,L13`|1|
|`A2`,all marks|`D13`,center10|23/64|
|`A1` or `A4,E4`|`D13`,mark9|1/2|
|`A1` or `A2`,all marks|`D14`,mark10|1/2|
|`A1`|`D14`,mark9|1/16|
|`A1`|`D15`,mark11|1/16|
|`A2`,all marks, or `A3`,middle1|`D10`,center7|1|
|`A1` or `A2`,all marks, or `A3,E3` or `A4,E4`|`B9`,doubled terminal8|1/2|
|`A1`|`B10`,doubled terminal9|1|
|`A2`,all marks|`B10`,doubled terminal9|1/2|
|`A1`|`B11`,mark9|1/2|
|`A1`|`B12`,mark9|1/2|
|`A_r,E_r`,1<=r<=8|`D9,L9`|5/32|
|`A_r,E_r`,1<=r<=4|`A9`,middle4|1/2|
|`A1`|`A11`,middle5|1|
|`A2`,all marks|`A11`,middle5|1/2|
|`A1` or `A2`,all marks, or `A4,E4`|`A11,E11`|1/2|
|`E6`,long ends, or `D5`,long end0|`A11,P2`|1|
|`A1`|`A12`,marks5,6|1|
|`A1`|`A12,E12`|1/4|
|`A4,E4`|`A12,E12` or `A12,P3`|1|
|`D5`,long end0|`A12,P2`|1|
|`D5`,long end0|`A12,E12`|3/8|
|`D5,L5`|`A12,P3`|1|
|`A1`|`A13`,middle6|41/64|
|`A2`,all marks|`A13`,middle6|23/64|
|`A1`|`A13,E13`|1/2|
|`A1`|`A13,P2`|1/4|
|`A2`,all marks|`A13,E13`|1|
|`A4,E4`|`A13,P2`|1|
|`A4`,internal marks1,2|`A13,P2`|1/2|
|`A2`,all marks|`A14,P2` or `A14,P3`|1/2|
|`A1`|`A14,P3` or `A14,P4`|1/2|
|`A1`|`A15,P3`|1/2|
|`A1`|`A16,P4`|1/2|
|`A1`|`A16,E16`|1/4|

For example the ordinary `D5` long-end/`D12` fork-center weight is one half,
by the broad row; it is not one of the short-leaf upgrades to one. This
distinction is needed in the order17 replacement proof. Likewise no ordinary
`A1/B8` payment is present. The fractional order11 scout containing new
`A2/B8,A3/B8` payments is not part of this ambient58 baseline.
The D5-long/D11-short-leaf row is the order16 unit transfer; its two
upgrades make the full D5/D11 component-pair charge17/2, not15/2.

## 5. Context is part of the local rule

All witnesses below are other, unmarked whole components at that corner.
The notation `A1 present` is a Boolean test, not a multiplicity factor.
The following same-component entries replace the ordinary entry:

|marked component/pair|guard|new value|
|---|---|---:|
|`A13,(0,12)`|A1 present|1|
|`A14,(0,13)`|A1 present|1/2|
|`A15,(0,13),(1,14)`|A1 present|1/4|
|`D15,(1,13),(1,14)`|A1 present|1/4|
|`A16,(1,12),(3,14)`|A1 present|5/8|
|`D16,(0,14),(0,15)`|A1 present|1|
|`D14,(0,12),(0,13)`|A1 or A2 or A3 present|1/2|
|`D13,(0,11),(0,12)`|at least three A15/D15 components|1/2|
|`D28,(0,12)`|at least two A15/D15 components|ordinary value+1/2|

There are also the following cross-component additions. They add to the
ordinary cross rule; all these ordinary roles are zero.

|marked pair|guard|added value|
|---|---|---:|
|A4 endpoint/A12 second node|A1 present|1/8|
|A2 mark/A4 endpoint|A1 and A10 present|1/4|
|A2 mark/D5 long endpoint|A1 and A9 present|1/2|
|A1/D5 long endpoint|A11 present OR both A9,A2 present|3/8|
|A1/D5 long endpoint|A9 and A2 present|5/8|
|two distinct A2 components,any marks|A10 and at least two A1 present|9/32|
|A12 endpoint/D16 long endpoint|at least two A15/D15 present|1/2|

The two A1/D5 entries deliberately add: its weight is one when A9 and A2
are present, even if A11 is also present. There is no additional3/8 for a
second witness. The marked A2 components in the sixth row are distinct.

For comparison with the forest recurrence, put `a_i=#A_i`, `d_i=#D_i`,
`k=a15+d15`. The total difference between the full rule and the ordinary
same/cross sums, in scale64, is exactly

```text
36 a2(a2-1) [a1>=2,a10>0]
+64 a2 d5 [a1>0,a9>0]
+24 a1 d5 [a11>0 OR (a9>0,a2>0)]
+40 a1 d5 [a9>0,a2>0]
+32(a15+d15)[a1>0]
+64 d14[a1+a2+a3>0]
+[a1>0](80a16+32a14+64d16+64a13
          +32a4 a12+64a2 a4[a10>0])
-64d13[k>=3]+32d28[k>=2]+64a12 d16[k>=2].
```

Each coefficient follows by counting the explicit marked pairs. For
instance there are four A2/A2 mark pairs per unordered component pair, so
`4*(9/32)*binom(a2,2)*64=36a2(a2-1)`. The A1/A16 endpoint-quarter cost
is ordinary and therefore is not included in this polynomial. This identity
is what permits the same intrinsic function to feed both geometric routing
and the finite-forest dynamic program.

## 6. Complete ideal rule

Here `m` is the **order** (number of nodes), so the affine rank is `m-1`.
All pairs in an irreducible affine diagram of order at most9 have weight
one. This includes all exceptional affine diagrams. Classical diagrams
of orders10–17 use the following tables; classical orders at least18 use
the formula at the end of this section.

Canonical affine numberings are:

* `A[m]`: the cycle on0,...,m-1;
* `B[m]`: spine0,...,m-3, double bond(m-4,m-3), leavesm-2,m-1 at0;
* `C[m]`: path0,...,m-1 with both terminal bonds doubled;
* `D[m]`: spine0,...,m-5, leavesm-4,m-3 at0 and m-2,m-1 at m-5.

For `A[10]` start with all units and halve its five antipodal pairs. For
`A[11]` halve its eleven distance-five pairs. For every `A[m],m>=12` keep
exactly the pairs of cycle distance at most five at unit weight.

For `B[m],C[m],D[m]` in the next table start with all units; the displayed
zeros and fractions are exhaustive. `X×Y` denotes all cross pairs, and
fractions given after a semicolon replace the default unit.

```text
m10:
 B H:(7,8),(7,9).
 C Z:(0,9).
 D weight27/32: {6,7}x{8,9}.

m11:
 B Z:(8,9),(8,10).
 C Z:(0,9),(0,10),(1,10).
 D H:{7,8}x{9,10}.

m12:
 B Z:(0,9),(9,10),(9,11); H:(8,10),(8,11).
 C Z:(0,10),(0,11),(1,10),(1,11).
 D Z:{8,9}x{10,11}.

m13:
 B Z:(0,10),(9,11),(9,12),(10,11),(10,12);
   H:(1,10),(8,11),(8,12).
 C Z:(0,10),(0,11),(0,12),(1,11),(1,12),(2,10),(2,12).
 D Z:(0,8),{9,10}x{11,12};
   H:(0,11),(0,12),(8,9),(8,10).

m14:
 B Z:(0,10),(0,11),(1,10),(1,11),(2,11),
     (9,12),(9,13),(10,12),(10,13),(11,12),(11,13).
 C Z:(0,10),(0,11),(0,12),(0,13),(1,11),(1,12),(1,13),
     (2,11),(2,12),(2,13),(3,13).
 D Z:(0,9),(0,12),(0,13),(9,10),(9,11),{10,11}x{12,13};
   H:(1,12),(1,13),(8,10),(8,11).

m15:
 B Z:(0,11),(0,12),(1,11),(1,12),(2,12),(3,12),
     {9,10,11,12}x{13,14}.
 C Z:all path pairs of distance>=10 except (2,12).
 D Z:(0,10),(0,13),(0,14),(1,13),(1,14),
     (9,11),(9,12),(10,11),(10,12),{11,12}x{13,14};
   weight15/16:(2,13),(2,14),(8,11),(8,12).

m16:
 B Z:all pairs of graph distance>=10.
 C Z:all path pairs of distance>=10 except (1,11),(4,14).
 D Z:(0,11),(0,14),(0,15),(2,14),(2,15),(9,12),(9,13),
     (11,12),(11,13),{12,13}x{14,15};
   H:(0,9),(0,10),(1,11),(2,11);
   weight7/16:(1,14),(1,15),(10,12),(10,13).

m17:
 B Z:all pairs of graph distance>=10.
 C Z:all path pairs of distance>=10.
 D Z:(0,10),(0,12),(0,15),(0,16),(1,11),(2,12),
     (3,15),(3,16),(9,13),(9,14),(12,13),(12,14),
     {13,14}x{15,16};
   H:(0,11),(1,12),(1,15),(1,16),(11,13),(11,14);
   weight1/8:(2,15),(2,16),(10,13),(10,14).
```

For every classical affine diagram of order `m>=18`, features are branching
vertices (degree at least three) and endpoints of nonsimple bonds. With
`delta` the distance to that set, keep a unit precisely when

    d(a,b)<=5, or 6<=d(a,b)<=9 and min(delta(a),delta(b))<=5.

All other pairs have weight zero. A simply-laced cycle has no features, so
only its distance-at-most-five pairs remain. This convention differs from
the finite feature convention: a simple affine endpoint alone is not a
feature. There are no ideal context bonuses.

The exact family maxima divided by order give the following envelope for
affine rank `r=m-1`:

```text
r1..8: r/2;
r9..16: 71/16,53/11,31/6,11/2,40/7,367/60,411/64,227/34;
r17..23:119/18,128/19,137/20,146/21,155/22,163/23,85/12;
r>=24:(5r+56)/(r+1).
```

These are the `c(r)` used by the terminal cusp algorithm. In particular,
orders10–17 have maxima355/8,53,62,143/2,80,367/4,411/4,227/2.
The ideal rule above, not a larger interim table from a local predecessor,
is the function whose maximum is taken.

## 7. What an exhaustive local certificate must establish

The geometric input is an **actual** triangular or quadrilateral two-face.
The white roots support that face; adjacent black roots mark its vertices.
For a quadrilateral write the opposite pairs `{a,c}`,`{b,d}`. Every white
plus one-black edge diagram is positive definite; the corner diagrams are
positive semidefinite; each white plus an opposite black pair is indefinite;
the full Gram has at most one negative direction. None of these assertions
is a combinatorial property of an arbitrary abstract four-marked graph.

For reference the common geometric section proves: positive inverse
projection gives actual finite vertices; a positive affine kernel gives
an actual cusp factor; a radical of a timelike-spanning core maps to zero;
nonnegative inverse/residual bounds control a first new supporting contact;
and a full affine factor splits from the positive spectator forest. Those
lemmas are the semantic meaning of the certificate fields below.

**Local certificate theorem (finite domain plus analytic residual clauses).**
Fix a nonnegative intrinsic weight function on whole marked components.
Suppose a certificate has the following data and verified properties.

1. An explicitly generated list contains every minimal opposite-side core
   with its two finite black deletions, every allowed labelled join and
   every injective partial identification of two such cores. Core orders,
   radii, finite attachment alphabets and all prunings are proved before
   generation. A symbolic clause covers each remaining infinite family.
2. For each retained core a finite record specifies the displayed actual
   corner, the marked pair, all compatible **induced** embeddings into
   a candidate whole finite/affine component, including failed searches,
   and any possible first outside supporting contact. Extension records
   use one of the geometric lemmas just stated, or explicitly produce
   a further independently paid actual corner.
3. Every terminal payment names the marked nodes, their complete whole
   components, every context witness and the required lower weight from
   Sections2–6. The sum at distinct actual corners is at least two for
   a triangle and at least one for a quadrilateral.
4. Replacement references are well founded: their dependency graph is
   acyclic, or a displayed integer measure strictly decreases. Every
   retained old output is checked at its actual marked role. A decreased
   output is either impossible under its guard or replaced by distinct
   actual corners whose weights suffice without referring back to it.

Then every actual face has the stated payment. To prove this, choose its
two minimal opposite-side cores. Completeness in1 puts their union in a
generated record, unless an analytic clause disposes of it. Reintroduce
the remaining supporting forest by its first contact with the displayed
core. The radical/inverse/attached-forest alternatives in2 either exclude
that first contact, preserve the asserted whole component, or pass to a
strictly earlier paid route. Induct on this measure (and on the number
of contacts before a first-entry seal). On a terminal record use3.
For a replacement use4 and sum the weights at the stated distinct actual
corners. Isomorphism invariance of the intrinsic rule removes canonical
labelling choices. Thus no abstract embedding, saved status flag, or sum
of component costs is being substituted for a geometric payment.

In an implementation the search-domain generator, the exact certificate
verifier and the mathematical completeness proof are separate objects.
A verifier that checks only successful records cannot establish1 or2.
Nor does a source digest prove that a caller's routing list was exhaustive.

## 8. Current geometric route families and their retention requirements

The following organization is mathematical, not a history of decreasing
weights. It specifies the families that a complete flattening of the current
certificate must retain.

**A. Short paths and triangles.** Resolve an adjacent marked path before
opening the two-opposite-core domain. Finite paths of length at most five
retain units. Deep finite distance-six, seven and eight outputs are reduced
by the dedicated shorter-path/actual-support arguments, with their protected
whole-small-factor alternatives; their large-rank deleted units are not
silently reused. An ideal short path invokes the whole-factor cycle or tree
route. Elementary triangles use their explicitly retained same-component
and cross-role pairs, and still have sum at least two. A terminal triangle
pair list must be checked against the full current tables, including zeros.

**B. No path of length at most eight: the core domain.** The complete
higher-order minimal-indefinite catalogue has26 diagrams and768 ordered
finite-deletion pairs. Orders four and five have respectively32 and14
canonical low-core types. All allowed black markings, ports, joins,
identifications and inherited nonedges are part of the domain; it is not
enough to embed the successful records from a prior search.

An order-two core has no white. An order-three core has one white and Gram

    [[2,-x,-u],[-x,2,-v],[-u,-v,2]].

The finite black-black label is unbounded: `x=2cos(pi/m)` can approach two;
the parallel value two is included. A dotted black-black pair is already an
indefinite proper order-two core and is handled there. The actual-edge and
port conditions bound the white-incident labels, not `m`. The only surviving
high port for this low/high join is node9 of coreL25, where the full inverse
diagonal is zero. Eliminating the high core therefore gives identically
zero Schur correction to this small Gram, for **every real** `x,u,v` for
which it is indefinite and every admissible cross column. It supplies a
second negative direction. This exact identity is the unbounded-I2 clause;
replacing it by a sampled finite list of labels would leave a gap.

For orders four/five the derived finite label bounds are six/five. In
disjoint low/high joins actual inverse bounds give cross-bond square
`<7/3`, leaving only simple/double bonds. The domain has726 single-edge
labelled joins:668 have two negative directions and58 violate white positive
definiteness. All321 compatible multi-position subsets fail the simple-bond
actual-edge test by nonnegative vectors, so increasing any native contact
or adding further such contacts cannot repair them. Overlaps include every
eligible size, including two shared nodes. The surviving low/high route
produces a whole corner component of order at most12 after arbitrary
supporting roots are restored.

For high/high overlap the generator considers **every** injective partial
identification preserving labelled edges and nonedges and the marked-distance
condition, and tests all remaining cross contacts. There are82197 consistent
identifications and2614 retained contexts on154 cores. Its six semantic
methods comprise1329 radical seals,823 positive-residual seals,259 bounded
component cases,120 isolated-inverse cases,74 protected feature/leaf cases
and9 adaptive cases. Their counts sum to2614, but the extension lemmas,
not this sum alone, prove completeness with arbitrary supporting forests.
High/high disjoint cores have their independent two-negative obstruction.

**C. Whole affine factors.** A displayed affine factor has a strictly
positive kernel. At an actual cusp all further supporting rows annihilate
that kernel; since their products are nonpositive, each is orthogonal to
every factor root. A nonzero coefficient of a full-core radical then
restricts the other old-core contacts. A first additional white contacting
both relevant black sides produces an earlier short path; otherwise the
whole affine factor and the typed finite alternatives stay whole. This
argument applies to arbitrarily many attached whites, not only the saved
one-white trials.

The small-factor branch has1153 eligible native patterns. Its induced-map
searches into whole affine orders10–17 must be complete at every target
order; the intermediate target partitions into bounded smaller factors,
protected/adaptive cases, displayed-too-large components and eligible
patterns. Large factors use the formula in Section6. Cycles use uniform
path/tail bounds; branching and doubled-end trees use positive-kernel
Schur profiles and the nonnegative inverse bounds. The local proof must
keep the actual distinct finite coupons from Section4 or5, not merely the
ideal-family total or a corner-independent representative envelope.

**D. Sparse finite roles and context guards.** Whole finite ranks12–17
are classified using the complete induced embedding domain, then all
arbitrary-support seals. The sparse tables in Section3 are not justified
by absence from selected maps alone: universal ideal/cycle/protected outputs
may consume a role absent from that particular search. Every chosen output
from A–C is retained unless the corresponding physical star has a proved
replacement. The following two final interfaces illustrate required strength.

* The zero `B13,(0,12)` has two genuine old affine-B singleton consumers,
  contexts166 and1744. Their full-core first-contact seals preserve a
  **distinct** whole A1/D12-short-leaf unit. Other universal whole-B13
  outputs include `(0,10)`, `(0,11)` and `(1,12)`, which remain units.
  In particular332 positive D12 outputs (280 radical/saturated and52 other)
  require retention of `(1,12)`; its absence in the801-pattern finite
  search does not license deleting it. Current216 literal configurations,
  564 later finite-transfer corners and24 D14 variants are retention
  interfaces, not a new enumeration of every actual face.
* For the rank58 D13 guard, exactly three unmarked whole rank15 components
  remain. Its two decreased leaf payments are repaired by the full physical
  classification, not by ranking formal forests. In the quiet affine-D12
  branch at most one rank15 component is active; two inactive components
  remain whole witnesses. The remaining black attachment gives either
  a distinct D28 pair `(0,12)` or an A12-end/D16-long-end pair, each paying
  the new guarded half. The old quiet D13 half plus that distinct half is
  one. The other finite and ideal consumers retain their unchanged roles;
  no replacement points back to a decreased guarded D13 pair.

An automorphism orbit in these tests means a marked orbit in the full
labelled whole component, not in an accidentally truncated displayed path.
The A1/A16 endpoint quarters, the A1/D16 fork halves and the guarded
same-component restorations are different resources. A proof that retains
one of them does not automatically retain the others.

## 9. Recovered definitions versus the outstanding exposition obligation

Sections1–6 are a finite, intrinsic, evaluable definition, with analytic
formulas in all unbounded ranks and no label cutoff on `I2(m)`. Their totals
feed the existing Bellman certificate exactly through the context polynomial
in Section5 and the cusp `c(r)` table in Section6. The current finite and
ideal scalar conclusions are respectively charge at most377 and absolute
cusp surplus at least36. Those conclusions also require the separate
complete counting certificates; they are not proved by listing the rule.

Sections7–8 give a precise route-certificate theorem and the actual family
partition needed to use it. The current local reports provide the physical
extension lemmas and affected chosen-output checks in a distributed atlas.
This staging text does not reconstruct or require a new flattened
terminal-record run. The appropriate organization of the existing proof
is the distributed family theorem in Section12, with its finite domains,
exact receipt predicates and well-founded retention order. A declaration
that each stored `status` is PASS is not its proof. The generic certificate
theorem explains the receipt predicates; the individually recorded exact
finite propositions and the analytic extension proofs supply them.

The higher-dimensional tail59–97 uses earlier local variants of these
tables. It is not legitimate to substitute the stronger final ambient58
context rule into their old counting rows. A faithful tail specification
must give each row's marked-rule deltas or its exact finite generator,
then its own routing and charge proposition. The broad distance/feature
tail from98 onward and the narrow-facet tail from164 onward have separate
definitions and proofs. Their scalar envelopes alone do not reconstruct
their local search domains.

## 10. Source bindings for table recovery and route expansion

These bindings identify mathematical data used to recover the definition;
they are not a substitute for Sections1–8.

* `fresh_d17_a16_coupled_ambient58_frontier_20261018.md`, SHA256
  `6f510a745724c64197d4b202068d6e29d2f3ebe381821db0c5f67db07fadec97`:
  ordinary/base and full context composition before the final guard.
* `fresh_d17_leaf_a16_context_coupled_transfer_20261018.md`, SHA256
  `19ed13114b8dddcf5a03febea1180209fa53f52b360ca0ec8c58f2653f4c0738`:
  the A17 table, A16/D16/A14 contextual restorations and distinct-coupon
  compatibility. Its JSON `local_proof.A17_weighted_pairs` is a full table.
* `fresh_b13_endpoint_retention_20261018.md`, SHA256
  `fedf5747f39c26e006084a7af623402205d5114841e31a941d29ef5bd23a384d`:
  actual endpoint consumers and positive output retention.
* `fresh_d13_three_rank15_context_transfer_20261018.md`, SHA256
  `91e93c3f7ff5eb4e64b435eea2892bf20af08711aed131110f40e898a9322014`:
  final guards, full quiet-cusp classification and377 closure.
* `fresh_order10_low8_ambient58_frontier_20261018.md`, SHA256
  `77149e2f7e6d6a96f720bbca5bbef5ba7647ca1039be6386bac2c9f0216cf983`:
  the actual order10 cross roles and final cusp envelope/terminal closure.

For precise finite table recovery the relevant JSON mathematical fields are:
`fresh_finite_rank14_existing_payment_transfer.rules`,
`fresh_a15_d15_quarter_a16_zero_center_retention.weighted_tables`,
`fresh_d16_a1_long_leaf_d17_zero_ambient58_frontier.local_proof`,
`fresh_d17_a9_a2_context_zero_ambient58_frontier.local_proof.weighted_pairs`,
and the final A17 field just named. Apply the explicit Section3/5 deltas,
including the B13 zero, rather than the predecessor costs.

For affine tables the fields are `weighted_pairs` for order13D and
order14B/C/D; `rules[].weighted_pairs` for order11B/C, order12B/C,
order13B/C, order15B/D and order16C/D; `records[].weighted_pairs` for
order16B; and `local_proof.weighted_pairs` for final order17D. The order13D
array has only74 entries: its four absent opposite-fork-leaf pairs are zero.
By contrast most later arrays explicitly store every pair. This is why a
parser-independent mathematical default convention is essential.

The route expansion is in the low-core, high-disjoint, high-overlap,
finite8/ideal5 combined, depth6 finite6/7/8, whole-cycle/tree, rank12–17
finite and order10–17 affine local theorems named by Sections7–8. These
must be cited as mathematical lemmas with their domain and support clauses;
their scripts or historical filenames alone are not the lemmas.

## 11. A precise generation and receipt specification for the routing appendix

This section makes the finite algorithm in Section7 concrete enough to be
stated as a mathematical algorithm. It does not assert that a new flattened
final-weight receipt file has been computed. The previously recorded local
calculations are distinguished from that possible reorganization below.

### 11.1 Core generation and the infinite-label clause

A labelled Gram has diagonal2 and off-diagonal entries zero or
`-2cos(pi/m)`. A *minimal indefinite core* has indefinite full Gram and
positive-semidefinite every proper principal submatrix. Its graph is
connected: otherwise a connected component is already an indefinite proper
submatrix. For an opposite-side core additionally mark two nodes such that
both corresponding deletions are positive definite.

Generate the high cores as follows. List the connected affine diagrams
of ranks4–8: A/B/C/D in each rank, F4 in rank4, and E6,E7,E8 in their
respective ranks. There are24. To a source on `l` nodes attach one new
node by **each** nonzero vector in `{0,1,sqrt(2)}^l`. Retain the extension
if every one-node deletion is a finite or affine Coxeter forest and its
full Gram is indefinite. Quotient by complete labelled graph isomorphism;
refinement signatures may accelerate this test but are not acceptance
criteria. Retain all ordered finite-deletion markings of every result.

This generator has146262 candidate attachments,84 retained labelled
attachments and26 isomorphism classes, of orders6–10 with counts12,3,4,4,3.
There are768 ordered markings, or384 if the opposite black pair is
unordered. These are recorded results of the high-core census.

The generator is exhaustive under the published Lanner order-five bound
and minimal-indefinite order-ten bound. Indeed a core of order at least6
is not compact Lanner, so has a singular positive-semidefinite proper
principal submatrix and hence a connected affine component A. Connectivity
provides a node adjacent to A. Positivity of A's kernel makes their union
indefinite; minimality makes that union the whole core. Thus A has
codimension one. An edge label at least5 could be enlarged to a connected
five-node proper subdiagram, whereas no connected finite or affine
five-node diagram has such a label. A parallel pair with any adjacent
third node is already indefinite by its positive kernel; a dotted pair
is indefinite alone. Consequently the attachment alphabet above is derived,
not imposed as a numerical search cutoff. Testing all one-node deletions
tests all proper principal submatrices, because positive semidefiniteness
passes to principal submatrices.

Generate low cores of orders4,5 from all connected labelled graphs on
their node sets, using labels2–6 at order4 and2–5 at order5, accepting by
the same minimal-indefinite criterion, and retain every finite-deletion
marking. To prove these smaller alphabets exhaustive, enlarge any edge
to a connected proper three-node or four-node graph and apply the complete
finite/affine list there. The retained canonical counts are32 and14.
The order-two dotted family and order-three Gram in Section8B are analytic
families, not finite rows in this generator. In particular the black-black
parameter of order3 is never replaced by an upper bound on its finite
label. Its zero-Schur-correction exclusion works throughout the full
physical parameter interval and at the parallel endpoint.

### 11.2 Every union, including every nonedge and failed identification

For each marked core `C` put `delta_C(v)=min(d(v,a),d(v,c))`. In the
high-core catalogue the largest such value is seven. Suppose no adjacent
marked white-interior path has length at most eight. For two opposite-side
cores, identify only white nodes. A shared pair `(v,w)` must satisfy
`delta_C(v)+delta_D(w)>=9`. Generate every injective partial matching
of every positive size satisfying this inequality. Accept a matching only
if the full labelled relations on its matched nodes agree, **including
nonedges**. Black nodes are never identified, and a black node cannot
become a white node through a matching.

After forming this induced union, any still free cross bond must join an
unmatched left white to an unmatched right white. If its end distances sum
to less than eight it gives an earlier short path. Otherwise test all
allowed endpoint pairs. For disjoint cores test every subset of cross-bond
positions, not merely each single bond. A magnitude bound follows from
the actual white and black-deletion inverse tests. For example a bond
of magnitude `b` between two positive white blocks requires

    b^2 (G_left,white^-1)_vv (G_right,white^-1)_ww < 1.

Its two black-deletion tests also require the normalized single-end
projection to be less than2. These inequalities remain necessary in
the presence of further contacts, because the inverse of a positive
Coxeter Gram is entrywise nonnegative.

At this stage a pruning receipt is one of:

* a nonnegative nonzero vector with nonpositive square in a matrix required
  positive definite;
* a nonnegative vector with **strictly negative** square in a corner matrix
  only required positive semidefinite;
* a proved exact negative index of at least two;
* a displayed adjacent marked path of length at most eight.

A nonnegative-vector receipt at the smallest possible bond is valid for
every larger bond and every additional bond on its tested support, since
those changes only decrease its square. Negative index alone is not used
as a monotonicity assertion when further bonds are added. The domain
enumerates multiple bonds separately where that monotonic argument is
unavailable.

For disjoint high cores the saved domain has12135 ordered candidate joins,
9276 single-join labelled cases and10990 multiple-position subsets.
Of the latter,10987 have a simple-bond nonnegative edge witness. The
remaining three two-position subsets have six allowed label assignments;
all have an actual-edge obstruction or negative index two. Of the
single-position types262 have an exact zero inverse diagonal in one
full core, giving a zero Schur correction and hence two negative directions
for every contact magnitude. The other four endpoint types are checked
exactly. With no cross bond the cores already give two negative directions.
This proves the disjoint pruning, rather than merely observing that its
successful-record list is empty.

For high overlaps the saved induced matching domain has82197 consistent
matchings. Inverse/edge prunings reject75947 and3636, leaving2614.
For each retained union every free cross bond is either a short path or
has a nonnegative actual-edge witness at its minimum magnitude; there
are154502 shortcut tests and5528 such edge witnesses. Thus all omitted
cross relations are justified, rather than set to zero by convention.
Low/high joins and identifications use the same generation rules, their
smaller radius bounds, and their separate analytic/finite receipts.

### 11.3 Extension receipts and finite completion domains

The following receipt meanings suffice for arbitrary additional supporting
whites. All vectors and matrices are in the physical normalized Gram;
a positive diagonal rational congruence is permitted if its squared
scales are recorded and used in each norm test.

* **Radical.** Give a nonzero full-core kernel vector `z`. The core contains
  a timelike vector, so its radical vector is zero as an ambient vector.
  Every new row therefore annihilates `z`. Record every possible contact
  position. If all have nonzero `z_i`, a nonempty nonnegative contact vector
  must meet both sign sets. List the resulting adjacent marked paths;
  if all have length at most eight, either that independent route applies
  or no further supporting root contacts the core. The first root of
  an indirect attaching chain is covered by the same statement.
* **Inverse/residual.** Give the inverse `J` of the nonsingular Lorentz core
  and partition possible contact positions into two positive end blocks.
  Contacts meeting both ends give the displayed short path. Otherwise
  nonnegative inverse entries at one end exclude every contact whose
  projection square exceeds2, and every multiple-contact pattern whose
  projection already exceeds2 at native minimum magnitudes. For an
  allowed unit contact at i put `rho_i=2-J_ii`. New residual vectors live
  in a positive-definite orthogonal complement. For two first contacts
  i,j, record `J_ij`, `rho_i rho_j`, and the path bound. If
  `J_ij^2>rho_i rho_j`, Cauchy either excludes the contact pair or forces
  their nonzero bond and hence the recorded short path, depending on
  the sign of `J_ij`. Norm-zero residuals are actual zero vectors, so
  coincidence and forced-product tests are valid there.
* **Attached forest.** If every root of a new positive forest contacts the
  old core only at one unit position p, write u for its unit-or-zero
  contact column and r for `J_pp>1`. Lorentz signature requires
  `G_Y-r uu^T>=0`, hence `u^T G_Y^-1u<=1/r`. Every connected component
  of Y has at least one attachment, so inverse positivity and the
  all-rank bound `(G_Z^-1)_ii>=rank(Z)/(rank(Z)+1)` give
  `u^T G_Y^-1u>=M/(M+1)` for total rank M. Therefore
  `M<=1/(r-1)`. This bounds the entire indirect attached forest.
  When two old positions remain, their exact residual Gram and inverse
  provide the corresponding joint receipt before applying this bound.
* **Protected marked role.** Record the fixed incident neighbors and
  nonsimple bonds at each protected feature or leaf. Inverse tests
  forbid a first new contact there; a radical can provide the same
  prohibition. The receipt must protect a leaf in the **whole** component.
  A leaf of an induced subpath is insufficient. The resulting distance-nine
  finite leaf or radius-five ideal feature is then looked up in the
  current tables, including any fractional fork-leaf exception.
* **Affine factor.** Give a positive kernel and the selected actual corner.
  Every supporting-white row of its positive-semidefinite corner
  annihilates the kernel, so nonpositive contacts are all zero. This
  preserves the whole factor and splits the positive spectator forest.
  Any further full-core radical is handled as in the first item.

Each check is an equality or inequality over an explicitly specified exact
algebraic field. An unbounded contact ray is written as zero or a native
minimum plus a nonnegative parameter. A polynomial nonnegativity receipt
specifies the substitution and its entire coefficient list; nonnegative
coefficients with a strictly positive constant prove the inequality on
the whole ray. A sampled list of parameter values is not such a receipt.

The high-overlap receipts yield a selected whole-component bound at most17
except for the protected/adaptive outputs. For each bounded component H,
generate every whole classical target of every allowed rank up to that
bound, every exceptional target at its rank, and every **induced** labelled
embedding `H -> T`. Search failures are part of this quantified domain.
The marked images, rather than just `T`'s total cost, select the table entry.
For affine targets do the same by order. The four irreducible targets in
each order14–17 give16 target searches on1153 eligible patterns; their
whole-order bounds have multiplicities313,354,346,140. An actual bounded
completion cannot escape these searches: the old H is an induced diagram
in its full whole component, and the rank/order bound applies to that
component, not just its displayed subset.

An unbounded component with an `I2(m),m>=7` edge cannot acquire a third
node in a finite or affine corner, by the finite/affine classification.
The `m=6` affine G2 possibility has order3 and a full unit rule. Together
with the analytic order-three low-core clause, this prevents an implicit
large-label cutoff from entering the target step. The ordinary rank-two
cross weight is independent of m in every broad row of Section4.

### 11.4 Triangle receipts, including all extensions

For a triangle take a minimal indefinite core containing its three black
roots. A corner path avoiding the third black persists in the whole actual
corner, and a displayed branch/nonsimple bond persists there. A core of
order at most6 has a Steiner tree through the marks with at least two
marked paths avoiding the third mark, each of length at most five. This
handles all labels in the small orders without enumeration.

For orders7–10 generate **every** marked triple in every high core. The
counts are `3*binom(7,3)=105`, `4*binom(8,3)=224`,
`4*binom(9,3)=336`, and `3*binom(10,3)=360`, totaling1025.
Delete the third black when computing each corner path and feature.
The elementary test retains two payments except at four triples of the
simple10-node E10 core, with edges

    0-1,0-2,2-3,0-4,4-5,5-6,6-7,7-8,8-9.

The exceptions and distinct alternate pairs are

|triple|alternate pair|length|fixed simple leaf and neighbor|
|---|---|---:|---|
|0,1,9|0,9|6|3 and2|
|0,2,9|0,9|6|1 and0|
|1,2,8|1,8|6|1 and0|
|1,2,9|1,9|7|1 and0|

In the first, second and fourth rows, the alternate corner is protected
by the other corner's affine E8 subdiagram0,...,8. Its positive kernel
is `(6,3,4,2,5,4,3,2,1)`. Kernel orthogonality forbids all extra supporting
contacts to those nodes. In the third row use finite E8 on0,...,7; its
inverse is entrywise nonnegative with diagonal
`(30,8,14,4,20,12,6,2)`. A new native contact at0,...,6 has projection
greater than2, so only a contact at7 is possible. This still protects
the indicated leaf and its degree-two neighbor.

No classical irreducible affine diagram has a simple leaf with a degree-two
neighbor. Hence each alternate marked component is finite or exceptional
affine of order at most9. Its alternate unit persists under arbitrary
additional support. For the elementary pairs, all54 distance-six cases
have a persistent feature within five steps, as do all distance-seven
elementary cases. Thus the depth-six finite deletions remove none of
these two selected triangle payments. The three alternate length-six
paths have protected leaves within two, one and zero steps of a mark;
the fourth length-seven path has its protected leaf. Their sparse
finite-rank modifications occur at larger marked distances and do not
alter these entries. All ideal pairs used here have distance at most5,
or the stated protected short features/exceptional whole factors; no
changed affine10–17 entry is spent as a missing elementary unit.

This gives a complete mathematical triangle proof once the finite1025-row
elementary path assertion is supplied as its recorded exact proposition.
Its generator is the marked-triple loop just specified, not a list of
previously successful triples. The kernel/inverse arguments prove the
supporting-white scope of the four exceptional receipts.

### 11.5 Final-weight receipts and the limit of the present extraction

A final receipt consists of a physical Gram template, the supporting-white
indices, its actual marked corners, extension certificates from11.3,
and a rational lower contribution at each chosen corner. Reconstruct the
whole finite/affine roles and guards, apply Sections2–6, and check that
each stated lower contribution is at most the resulting actual weight.
Distinct corner indices must remain distinct under the template's
identifications. Sum these contributions and demand at least one (two
for triangle receipts). A printed scalar family cost cannot be used
as any individual corner contribution.

The final local scripts perform this role retention in several separate
domains, including the216 physical literal stars,564 later finite corners,
24 D14 variants,332 positive D12 outputs and analytic universal outputs.
For the final negative D13 guard, the selected completion domain is801
patterns and862 maps with33 affected contexts/76 maps. Its exact physical
division is five impossible proper completions (forced projection6,12,30
or102), four radical-sealed cores1985,2378,2575,2581, twelve cases with a
whole small orthogonal component, and twelve attached-forest cases.
In the latter the old small rank is1,2,3 and the inverse square is4/3,3/2,2;
11.3 bounds new ranks by3,2,1, making the whole component rank at most4.
Each possible extra small factor contradicts the exact rank58 guard
`D13+three rank15`. This is an arbitrary-support proof of the guard, not
just76 failed attempts at adding a missing root.

Universal quiet affine-D12 outputs require a separate receipt. Its marks
are opposite fork leaves of the whole order13 factor and its spectator
forest splits orthogonally. A native adjacent-black bond or a shared
spectator contact gives an earlier length-one/two payment. Otherwise the
two nonnegative spectator contact vectors have disjoint support. Complete
Schur profiles give end costs at most8, nonextreme costs at most4. In the
nonextreme equality8+4 branch the spectator vectors vanish. The three
nonextreme cases give respectively a retained A1/D12 unit, an independently
handled whole affine-B13 corner, or a distinct whole D13 unit with a
small component forcing its new guard to fail. In the two-extreme branch
at most one of the three rank15 components is active; two remain whole
unmarked witnesses. The remaining attachment yields the D28 or A12/D16
guarded half of Section8D. Thus the quiet half plus its distinct new half
is one. These are analytic family receipts, not records to be inferred
from the finite801-pattern search.

This specification supplies a well-defined candidate flat certificate:
take the complete domain11.1–11.2, the exhaustive target embeddings11.3,
the short/cycle/tree analytic clauses, and expand each nonterminal reference
until it ends in a final receipt evaluated as above. A verifier must either
check acyclicity of the explicit reference graph or verify a strict
decrease measure for each analytic reference. A route may not be marked
accepted merely because the corresponding predecessor script exists.

**Recorded versus unperformed.** The local reports record their own
complete searches and exact receipt checks; no such recorded run has
been repeated here. Section12 uses those distributed finite propositions,
not a new single-file flattened certificate or a new claimed computation.
The generic flat verifier is an optional future audit organization. The
domain counts alone are not the finite propositions: their exact receipt
identities and inequalities must be retained as well.

## 12. Distributed routing-retention theorem

The following is the mathematical composition to retain in a standalone
presentation. In its finite propositions, *all* means the domain generated
in11.1–11.3, not the support on which a previous search succeeded. The
Gram lists, marked images, positive kernels, inverse entries, coefficient
lists and chosen corner contributions are finite supplementary tables;
their row domains and checking predicates have just been defined. Native
root magnitudes and unbounded parameters are handled analytically where
stated. No proposition here asserts that a finite list enumerates every
actual supporting forest.

**Distributed finite propositions.** The recorded local computations have
the following forms.

1. The core domain and pruning receipts have the results in11.1–11.2.
   Every one of the2614 surviving high-overlap rows admits a receipt in
   one of the six classes counted in Section8B. Every low/high surviving
   row has a component bound at most12 or an independently paid short
   path. Every protected/adaptive output is a whole finite marked leaf
   at distance9, a retained feature-near ideal pair, or a direct
   multi-corner receipt.
2. For each selected bounded component, the full induced target searches
   are those of11.3. A final nonunit entry in Sections3 or6 has a
   proper-completion exclusion, a whole-core seal, a classified additional
   white, or a residual/attached-forest receipt. After that receipt the
   sum of the **chosen current** contributions is at least one. Records
   referring to a previous contribution are retained only through the
   intrinsic role test and the family clauses below.
3. The triangle generator has1025 rows and the four exceptions in11.4.
   Every other row has two persistent elementary pairs; their explicit
   path/feature data pass the final-weight lookup, including the54
   distance-six cases. The four exceptions have the protected alternate
   unit proved there.
4. The current literal and later-output retention domains in11.5 retain
   their positive contributions at the stated actual marked roles.
   The universal outputs not stored in those fixed Gram lists have
   precisely the family descriptions below. Their proof is by the
   same analytic profile and arbitrary-support arguments, not by the
   finite row counts.

These propositions are checked by finite exact loops and exact symbolic
coefficient identities. Termination follows from the explicit core orders,
finite labelling alphabets, finite partial matchings, finite target lists,
and the derived whole-component bounds. Algebraic sign tests use the
physical embedding; no approximate eigenvalue acceptance is allowed.
The required conclusions are the receipt predicates, not just the printed
word PASS. The existing sources record these separate computations;
this appendix does not claim another execution of them.

### 12.1 Independent short receipts and the order of use

First discharge any actual path of length at most five, allowing the
cycle10 antipodal and cycle11 distance-five halves. Their independent
whole-face arguments give a full payment through a second distinct half,
shorter paths, or whole finite alternatives. The cycle10 mixed equality
has only whole A7 or D6 active spectators, yielding whole A9 endpoints
or E8 leaves at the opposite corner. Its other outputs are the explicit
small/A9 halves. The cycle11 length-six single-contact alternatives have
a fixed whole leaf within four steps of the old mark. Thus none uses a
deep finite6 zero or a decreased sparse finite long pair.

An independent zero-C proof is also a complete payment of the **face**,
not a unit assigned to its zero corner. The order11 and order12/13
near-terminal C proofs give a strict determinant gap unless an actual
length-one/two path already pays. The C10 endpoint proof has no finite
payment output. Their use may therefore precede a B proof or finite
selected proof that happens to encounter such an actual corner.

Next use the deep finite6, finite7 and finite8 theorems in that order.
The finite6 theorem's no-short-five branch has its complete E10 forced
white alternatives and exact Schur/equality receipts; its positive outputs
are shorter paths, protected whole marked leaves, small whole affine
factors or the A1/A12 middle-node unit. The finite7 theorem either uses
that previous route or its whole A8/D9/affine-D9/affine-E8 marked leaf.
Its cycle boundary outputs at distance7/8 have feature depth at most5.
The finite8 theorem's no-short-seven branch is excluded by its forced
E10 enlargement and strictly positive determinant polynomials; otherwise
it uses one of those earlier routes or a finite8 output with feature
depth at most5. These are not the deep finite8 entries being deleted.
This is a strict dependency order6→7→8. An apparent reference to a
finite path of the same length is harmless only if its fixed feature
or whole leaf has been proved; the final table lookup then already
retains it and no new recursive call is made.

The universal large-cycle and affine-tree corridor proofs use exactly
those protected short outputs, whole finite leaf9 outputs, exceptional
small factors and the typed cross roles in Section4. A large finite D
short-leaf9 half is not mistaken for a unit: its independent whole-star
fork proof supplies the distinct small/D10 half. A path outside all
finite/ideal feature neighborhoods has a simple corridor; the exact
inverse and Schur profile bounds either force an earlier short path or
one of the stated equality components. Thus corridor length is treated
by the all-rank inverse formulas, not a finite list of sampled lengths.

### 12.2 Small affine trees: closed profile systems

For the affected affine factors of orders10–17 the positive kernel first
splits off the spectator forest. If adjacent black roots are joined, or
both contact a common spectator node, a distinct path of length one/two
pays. Otherwise the two nonnegative spectator columns have disjoint
supports. The complete inverse-diagonal tests classify the native profiles
on each end; opposite-mark products remain zero or at least one with no
upper cutoff. The physical four-root determinant and the spectator
difference energy then give either a strict contradiction, an earlier
independent face payment, or an equality receipt. At equality Gram-metric
Cauchy determines which whole spectator components are active. This is
why a named finite alternative stays whole even if the *other*, omitted
black root has arbitrary additional spectator contacts.

The terminal roles of these closed systems are:

* Order10 B uses the A1/B9 and A3-end/B9 halves, the existing A2/B9 and
  A4-end/B9 halves, distinct affine B10 halves, or a current D10 fork
  payment. Its C alternative is handled independently. Order10 D uses
  its current27/32 pair and a distinct B10 or D10 payment where those
  appear. Its small finite outputs are A-end/D9 `5/32` roles. The genuine
  whole A17 endpoint output needs only5/32 and now has1/2, so remains paid.
* Order11 C zeros are independent. Order11 B zeros give a short path,
  whole finite B11/B12 unit, upgraded A1/B10 unit, or two distinct halves:
  A2/B10 plus A1/D10-short-leaf or affine-D10-opposite-leaf. The D11 half
  proof may previously have preferred a B11 half. If that actual B zero
  corner occurs, its complete B proof just described already pays the
  face; otherwise no B debit is spent. This is an actual-corner split,
  not a circular invocation of an obsolete B half.
* Order12/13 B near-terminal profiles give an earlier BC11 zero face,
  a distinct affine-D10 half, a whole B11–B13 unit at its specified
  endpoint/near-end role, or the A1/B11/B12 halves at canonical node9.
  Their C proofs are independent strict-gap/short-path proofs. The
  B13 long endpoint `(0,12)` is an exception handled by the full seal
  below; `(0,10),(0,11),(1,12)` retain units.
* The affine D12–D14 profile systems have whole A/D alternatives or
  their specified small/fork halves. Their extreme mixed spectator
  blocks are finite all-rank inverse-spectrum problems; nonunit native
  contacts strictly exceed the relevant energy threshold. The ordinary
  new weights are precisely the A1/A11,A1/A12,A2/D11,A2/D13,
  A4/D11,A5/D11,E6/D11,A1/D13 and A1/D14 roles in Section4.
  A mere induced small subgraph cannot activate one of these coupons.
* Orders15–17 use their full-core radical or saturated inverse receipts
  to protect the explicit A4/A11,A4/A12,A4/A13,D5/A11,D5/A12,D5/D12,
  A2/A13,A2/A14 and A2/D14 roles; the order15/16 near-fork sixteenths
  are A1/D14 at9 and A1/D15 at11. Their ranks and marked positions
  distinguish them from a rank3–5/D15-center payment, which is zero.
  The order17 center system uses the whole A2/D5 guarded half and
  a distinct D5-long/D12-center half. It does not consume a restored
  same-component D13 unit or one of the lowered D15 centers.

Every immediate affine alternative in these systems either contributes
its **current** fraction as a different corner of the same receipt, or
is one of the independent earlier zero/cycle proofs. In the first case
there is no recursive theorem call. In the second case the entire face
is already paid and the later computation stops. This is the well-founded
interpretation of their otherwise mutually suggestive diagram changes.

### 12.3 Joint finite sparse receipts and the guards

After the independent short/ideal routes, apply the complete finite-target
embedding and extension receipts. All nonunit same-component roles are
exactly Sections2–3 and5; the algorithm does not choose a rule by its
old aggregate cost. The following simultaneous substitutions are receipts
on **all four actual corners**, so their arrows are not recursive calls
from one decreased finite weight to another.

* At a D14 long-leaf half, the full proper-completion test excludes its
  missing native leaf when required. The whole cases have a radical
  seal or one uniquely classified saturated extra white. Their distinct
  other half is D15 `(1,12)`, D16 `(2,13)`, or the broad A4/D13-center
  half. Each remains one half in the current table.
* The four joint D14/A15/D15/A16 variants have contributions
  `1/2,1/4,1/4,0`. Where the quarters or the D14 half need context,
  their full physical seal proves the required whole A1/A2/A3 witness;
  an absent guard is a separately excluded or independently paid
  branch. An induced small spectator is not counted as a witness.
* The A16 guarded5/8 receipt has final contributions
  `0,1/8,1/4,5/8`. Its distinct guarded A2/A4 quarter has actual whole
  A1 and A10 witnesses. The D17/A16 coupled receipt has final
  contributions `0,1/4,1/4,1/2`. In each case all positive corners are
  distinct and their complete whole roles are reconstructed by the
  seal, so a quarter does not stand for a second use of the half.
* The A13,A14,A15,D15,A16,D16 restorations in Section5 are Boolean
  witness statements. The intrinsic full-component proofs show either
  that a chosen old positive role has its whole small witness and
  hence the restored value, or that a different retained receipt pays.
  The simultaneous check uses the actual witness, not the number of
  witnesses or the forest representative that will later bound cost.
* Rank3–5/D15-center roles can be zero because none is a chosen positive
  contribution in the complete fork/cycle/selected/later-output theorem.
  The original fork profiles have a **whole A1**, sealed by the affine
  kernel/full radical, not a small path free to grow to ranks3–5. The
  same distinction retains A1/D16 fork halves while permitting the
  other small/D16 roles to be absent.

Two final consumers need their own universal clauses. A B13 endpoint
unit was used in actual affine-B singleton models166 and1744, even
though the finite801-pattern selected maps do not use it. Their
radical/first-contact seals give a distinct whole A1/D12-short-leaf
unit, which is independent of the endpoint deletion. All other chosen
whole-B13 roles listed above retain their unit. Finally the guarded
D13 fork role uses the exact rank-capacity and quiet-affine-D12 proof
in11.5: if the guard is false its old unit remains; if true its quiet
half plus the guarded D28 or A12/D16 half is one, unless an earlier
short/BC13 route already paid. These positive guards use two inactive
whole rank15 components. They do not assume a changed D13 unit to
finance the replacement.

### 12.4 Composition and scope

For an actual triangle use11.4. For an actual quadrilateral first test
the independent short/zero/cycle clauses. A call to one of those clauses
ends the face proof. The deep6→7→8 order and protected-current outputs
justify every use of a shorter path. If no such path pays, choose the
two minimal cores and use the complete union domain. A pruned union is
impossible or gives an earlier path; an unpruned union has a whole bound,
protected role or adaptive receipt. Classify its whole target and apply
the small-affine closed system or joint finite receipt. Every positive
terminal role has its current value by12.2–12.3, including the two
unstored universal consumers just separated. Its distinct actual corners
sum to one. Additional supporting roots are covered by first contact,
radical orthogonality and the attached-forest bounds, not by an implicit
cap on the number of whites.

This proves the distributed routing-retention theorem from the explicitly
stated finite propositions and analytic clauses. Its dependency order is
semantic: independent actual-face payment, shorter-path/protected receipt,
whole-factor profile system, joint finite receipt. A joint receipt is a
single terminal inequality, not an arrow back into another decreased
role. The original source sequence is one way of having verified these
finite propositions and their chosen-role interfaces; its dates are not
part of the mathematical rule or proof.

For publication, the finite supplementary tables must remain available
at their specified domains, and their recorded verification results must
be cited next to the propositions. This text has recovered the rule and
proof architecture, not independently rerun or recursively rechecked every
coefficient of that distributed supplement. That limitation concerns the
audit performed while preparing this staging text, not a new requirement
for a flat checker, a new claimed run, or a substitute for the analytic
support proofs.

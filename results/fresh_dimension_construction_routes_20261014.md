# A two-coordinate lift of the nonarithmetic H18 hybrid

The literature comparison led to a concrete new experiment: lift the known Vinberg H18 hybrid through **two transverse coordinates**, instead of doubling it or lifting the old H19 chamber through one coordinate. The result is an exact **60-facet Coxeter chamber in H20**, with three explicitly bounded repaired cusp links. A small exact spacelike escape proves that it still has **infinite volume**. This is a reproducible new starting geometry within the project, not a new finite-volume dimension or a worldwide novelty claim.

The raw construction admits 8,192 labelled axis assignments modulo exchanging the two axes. Every raw assignment has an explicit escape. A finite coloring calculation chooses the unique assignment minimizing a specific cusp-repair obligation; its three deficient links are then repaired coherently. The global escape test is performed after those repairs, rather than inferring volume from the number of closed links.

Artifacts: [exact script](../scripts/fresh_dimension_construction_routes_20261014.py) and [frozen field/vector certificate](fresh_dimension_construction_routes_20261014.json). Run:

```sh
python3 scripts/fresh_dimension_construction_routes_20261014.py --verify
```

## What the primary constructions actually provide

Vinberg's H18 construction glues the chambers of `I_(18,1)` and `II_(17,1)+<2>` across their common orthogonal H17 facet. The gluing theorem adds their volumes in the **same ambient dimension**; its arithmetic parents supply a nonarithmetic hybrid. It does not construct an H19 or H20 polytope from the H18 output. The explicit parent normals appear in §7. [Vinberg, *Non-arithmetic hyperbolic reflection groups in higher dimensions*, §§3, 5, 7](https://sfb701.math.uni-bielefeld.de/preprints/sfb14047.pdf).

Allcock's doubling construction likewise takes the union of a chamber and its reflection across a wall whose meeting angles have even denominators. Lemma 1.3 and Theorem 4.1 preserve the ambient dimension. They are useful once a chamber in the desired dimension exists, but do not provide the missing dimension-raising step. [Allcock, *Infinitely many hyperbolic Coxeter groups through dimension 19*](https://arxiv.org/pdf/0903.0138).

The 2024 Makarov-garland result produces infinitely many commensurability classes in H4 and H5 by gluing compact blocks. The July 2026 paper constructs noncompact families in dimensions 4, 5, 6, 7 and 9; its gluing operation is explicitly dimension-preserving. Neither supplies an automatic compact H9 extension or a finite-volume H20 extension. [Bogachev–Douba–Raimbault, §§3–4](https://arxiv.org/html/2309.07691), [Guo–Ma–Zang–Zheng, §3.3](https://arxiv.org/html/2607.14715v1).

These are construction-scope observations, not a claim that no undiscovered dimension-raising method exists.

## Comparison with the existing searches

The repository already reconstructed Vinberg's 40-wall hybrid and checked its relevant lower-dimensional faces. It also exhausted direct seam twists and first parent doubles. The [H18 gluing/doubling report](fresh_gluing_h18_twists_doubles_20261013.md) records why those changes did not repair the targeted face angles.

The [H19-to-H20 lift](fresh_dimension_bridge_20261013.md) retains the projected directions of a different 36-wall source and uses one added positive coordinate. Its complete one-coordinate height classification is not a classification of the present two-coordinate lifts. The two [fixed-direction H8-to-H9 lifts](fresh_h8_lift_h9_20261013.md) are also a different construction and have their own positive-null-vector compactness obstructions.

The experiment below uses the published H18 hybrid's **norm-one** walls, keeps its nonrational cross-parent pairings, and changes its metric rank by two. It does not search an integral reflective lattice in a forbidden arithmetic dimension.

## Exact signed-pair construction

Let the source roots be `r_i`, using the labels and Gram matrix in the [archived hybrid](../data/enumerations/vinberg_h18_mixture_facets_2026_10_13.json). There are 26 roots of square two and 14 of square one. Adjoin orthogonal unit vectors `e_0,e_1` to the old space of signature `(18,1)`.

Keep each norm-two root unchanged. Assign each norm-one root a color `a_i` in `{0,1}`, and replace it by the pair

```text
R_i,+ = r_i + e_(a_i),      R_i,- = r_i - e_(a_i).
```

Both new roots have square two and their mutual pairing is zero. Their average recovers the old root. There are 54 walls. If both colors occur, their span is the full space of signature `(20,1)`; the constant-color case has an unused transverse direction and is immediately infinite-volume.

Every color assignment satisfies the complete Coxeter pair gate. Between different norm-one source roots, every old pairing is one of `-1,-3,-5,-7,-9`. Added height pairings are `-1,0,1`. An old parallel pair therefore becomes right, label three or parallel; every more negative pair remains parallel or dotted. The script also directly checks every mixed short/long entry, including the hybrid's irrational entries. Those pairs are right, label three, label four, parallel or dotted. All old long/long relations are unchanged.

The old strict timelike point `h=sum r_i`, extended by zero heights, stays strict. A complete nonpositive pair gate then makes every root an actual facet: its orthogonal projection of `h` is timelike and strict on every other wall.

## The global test that prevents a false positive

The three source roots `L18,L25,L26` have square one and pairwise pairing minus one. Every two-color assignment has a monochromatic pair `i,j` among them. Put

```text
delta = r_i+r_j.
```

It is a feasible future null ray. If the common color is zero, its unused height is `e_1`; exchange axes in the other case. The script verifies, with every other short-root color allowed,

```text
X = delta + h/100000 + e_1,
X² > 0,      (X,h)<0,      (X,R)<0 for every lifted wall R.
```

Thus **all raw axis assignments are infinite-volume**, without testing their global geometry one by one. This does not exclude later walls with new projected directions.

There are 25 old parallel pairs among the 14 norm-one roots. Enumerating the 8,192 colorings modulo global color swap gives a unique minimum of **three** monochromatic pairs: assign `L18,L25,L26` to one axis and the other eleven short roots to the other. This is a minimum only for these 25 pair constraints; it is not a complete inherited-cusp census or a global optimality theorem.

## Six coherent repairs, followed by another global test

For each of those three monochromatic pairs, add

```text
C_delta,+ = 2delta + e_missing,
C_delta,- = 2delta - e_missing.
```

These roots have square one and form a parallel pair. The complete 60-wall pair census is:

| Relation | Count |
| --- | ---: |
| Right | 971 |
| Label three | 127 |
| Label four | 4 |
| Parallel | 79 |
| Dotted | 589 |
| Total | 1,770 |

All 60 walls are actual facets by the common-center projection argument. At each of the three repaired null rays the script checks all 60 signs, exactly 24 incident walls, normal span 20 and horospherical Gram rank 19. The five orthogonal affine components have ranks `8,8,1,1,1`. It computes a strictly positive kernel for every component and verifies that its weighted normal sum is a positive multiple of `delta`. Hence the link has no nonzero recession direction and is bounded. Rank alone is not used as a boundedness test.

Nevertheless the final chamber escapes. In the first nineteen source roots used as an old basis, followed by the two added coordinates, set

```text
v = (3,6,10,9,8,7,6,5,4,3,2,1,0,-1,-2,-3,-4,-5,5,5,1).
```

It has square 25 and is feasible on all 60 walls. The strict point

```text
E = v + h/100
```

has square `28281/1250-(36/125)sqrt(2)>0`, negative pairing with `h`, and strictly negative pairing with every wall. It certifies an open ideal-boundary patch and infinite volume. An initial numerical linear program suggested the direction; the final checker uses only the displayed literal exact point, with no numerical acceptance predicate.

No complete vertex, edge or cusp census of the 60-wall chamber is asserted.

## Arithmetic gate before further work

The unchanged norm-two pair `L1,R10` has pairing `-2-2sqrt(2)`, so its Cartan two-cycle is `12+8sqrt(2)`. The cycle field therefore contains `Q(sqrt(2))`. The unchanged pair `R23,R24` has rational Gram matrix

```text
[[2,-6],[-6,2]],      determinant -32.
```

This principal block remains indefinite under the conjugation of the cycle field. In any hypothetical quasi-arithmetic finite-volume completion retaining these walls, total reality would extend that conjugation to the full Gram field, contradicting Vinberg's positive-semidefiniteness condition. Thus such a completion would be **non-quasi-arithmetic**, in particular nonarithmetic. This is conditional on finite volume; the current group is not being called a nonarithmetic lattice. See the criterion stated in [Guo–Ma–Zang–Zheng, Theorem-Definition 3.1](https://arxiv.org/html/2607.14715v1#S3.SS1).

## An actionable global closure test

This lift offers a lower-dimensional global test rather than another unstructured cap bank. Write a point as `(x,z_0,z_1)`. The signed pairs force the old projection `x` into the source chamber cone. Let `f_a(x)` be the minimum of the positive linear quantities `-(x,r_i)` assigned to axis `a`, together with `-2(x,delta)` for repairs on that axis. The full new inequalities are exactly

```text
x in the old H18 chamber cone,
|z_0| <= f_0(x),      |z_1| <= f_1(x).
```

Consequently a spacelike direction exists precisely when some old feasible `x` satisfies

```text
x² + f_0(x)² + f_1(x)² > 0.
```

For the finite-volume direction of this test, partition the old cone by the finitely many choices attaining the two minima. On each cell both functions are linear. Its lifted box is generated by the four signed lifts of each extreme base ray: interval boxes add under Minkowski sums, so decomposing a base point also decomposes every allowed fiber point. It is therefore sufficient and necessary to check the Lorentz norms of those finitely many lifted cell rays. Future timelike rays give finite vertices; future null rays must be included, not discarded. If every ray is future causal, the finite polyhedral chamber lies in the closed Lorentz cone and has only finite or ideal vertices, giving finite volume. This is a closure criterion, not a claim that the current failed model satisfies it.

The next bounded experiment should optimize or replace the two minimum systems so that they cut the literal `E` and then apply this cell-ray test. Merely changing the two-color assignment cannot remove the raw parallel-triangle obstruction, and merely counting the repaired links cannot establish volume. General non-axis height vectors, changed projected directions and different coupled repair profiles remain outside the exact family tested here.

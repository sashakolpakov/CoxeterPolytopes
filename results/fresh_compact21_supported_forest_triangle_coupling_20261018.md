# Seven nonadjacent triangles at a supported compact arithmetic rank21 vertex

18 October 2026. **New local incidence theorem; no improved dimension bound.**
Let P be a compact arithmetic hyperbolic Coxeter polytope of dimension21,
and suppose a vertex has diagram A1+D5+E7+E8. Among the42 incident two-faces
whose marks are nonadjacent roots in one of D5,E7,E8, at most **seven** are
triangles. Thus at least35 of these faces are nontriangular.

The complete necessary seven-triangle case uses the E7 third wall with
contact5 at root5, and both E8 third walls with contact5 at roots0 and6.
It has no designated D5 triangle. This equality description concerns the
checked necessary domain; it does not prove global realization.

The [producer](../scripts/fresh_compact21_supported_forest_triangle_coupling_20261018.py)
and [certificate](fresh_compact21_supported_forest_triangle_coupling_20261018.json)
classify the complete actual-corner triangle domains with **whole D5 and
A1 spectators included initially**. The program does not infer their
orthogonality from a selected old witness. It combines this with the
separate complete
[D5+A1 triangle theorem](fresh_compact21_d5_a1_triangle_coupling_20261018.md).
This targets the largest supported forest in the saved exact dimension21
weight obstruction, while preserving every old face-payment rule.

## Actual geometry and the finite domain

Use root square two and number E7/E8 as the path0-1-...-(n-2), with leaf
n-1 attached to2. Number D5 as0-1-2 with leaves3,4 at2. A nonadjacent pair
a,b specifies a two-face at the old finite vertex. If it is triangular,
its third wall c must make both old21-root corners obtained by replacing
a or b with c positive definite.

For D5 marks, the whole E7 and E8 are orthogonal to c by the preceding
D5 theorem. For E7 marks, the whole E8 has no connected finite extension
and is orthogonal to c. The full E7+D5+A1 is retained in the enumeration.
For E8 marks, suppose c contacts the untouched E7. Finite classification
forces E7+c to be E8, and maximality of that finite E8 forces c orthogonal
to every other old root retained in either new corner. Their union is the
whole old forest. Then the old21 roots together with c form the positive
rank22 forest (E7+c)+A1+D5+E8, impossible in signature(21,1). Hence this
whole E7 is orthogonal too; the enumeration retains E8+D5+A1.

Every contact appears in a finite corner. A label at least six at any
old D/E root would force an isolated I2(m), incompatible with its retained
old neighbors. At a marked root the other mark is nonadjacent, so an old
neighbor remains in the opposite corner. At a white root the two omitted
singletons cannot both contain all its old neighbors.

There is a separate argument for the whole isolated A1. If its contact
with c has label at least six, the isolated-I2 classification forces c
orthogonal to every other root retained in each finite corner, hence to
their union. The remaining exceptional spectator is already orthogonal.
The resulting height is (2cos(pi/m))^2/2-2<0, impossible for an actual new
wall containing the finite vertex strictly. This proves the complete
label domain two through five, including the A1 spectator, with no high
label sampling or unstated cutoff.

The producer enumerates the first finite corner as a finite tree. The
new root has degree at most three, meets each old connected component at
most once, and at degree three has only label-three bonds. Conversely
these support choices followed by the exact finite classification include
every first-corner possibility. It then restores the other marked root
with each label two through five and checks the second corner. These
are complete actual duties, rather than a small finite witness bank.

|Marked type|Old roots retained|Nonadjacent marked pairs|First-corner trials|Accepted first corners|Second-corner trials|Finite-eligible occurrences|
|---|---|---:|---:|---:|---:|---:|
|E7|E7+D5+A1|15|7,497|574|2,296|379|
|E8|E8+D5+A1|21|12,492|842|3,368|470|

The existing D5 theorem contributes the six D5 nonadjacent pairs, giving
6+15+21=42. The isolated A1 supplies no internal pair.

## Field, height and compactness checks

The completed internal compact21
[field reduction](fresh_compact21_single_field_20261016.md), pinned by
the producer, gives actual ground field k=Q(sqrt5). This is not the
three-field theorem valid only from dimension22.

Every simple cycle in an old component plus c consists of c and the
unique old path between two contacts. All old bonds have label three,
so the endpoint contact product decides its membership in k. Whole D5,
marked exceptional and isolated A1 components are treated separately;
their contacts are all included in the projection height. Write
c=p+sqrt(H)e at the common future unit vector e. The strict old vertex
inequality requires H=p squared minus two>0. Once cycle membership
passes, H belongs to k. The old Gram is rational positive definite under
the other real embedding, so conjugate positive semidefiniteness requires
H-conjugate at most zero by the exact Schur complement.

Every remaining candidate is also tested against all connected affine
principal subsets containing c, with caching only of identical contacts.
Subsets not containing c are finite already; disconnected affine subsets
have a connected affine component. Finite subsets may be skipped, as may
connected subsets containing a label-five bond, which no affine Coxeter
diagram possesses. Compactness excludes the affine alternatives through
their positive null kernels. No remaining necessary arithmetic root has
such an obstruction.

The exact filters give:

|Type|Nonpositive H|Cycle outside k|Conjugate not PSD|Remaining occurrences|
|---|---:|---:|---:|---:|
|E7|195|5|174|5|
|E8|231|4|231|4|

Finite/affine Coxeter classification, compact acute simplicity and the
Vinberg arithmetic criterion are explicit theoretical inputs. The field
reduction is an inherited internal theorem. None is newly established
by summing the finite records.

## Complete surviving third-wall types

Write s=sqrt5 and phi=(1+s)/2. The D5 list has the old five types and the
new A1-contact type proved in the preceding packet. The exceptional
domain, with spectators included, leaves exactly the following:

|Name|Exceptional contact magnitudes|Whole A1 contact|Squared height H|Nonadjacent triangle marks|
|---|---|---|---|---|
|E0|E7 root0: phi|0|1+s|(0,2)|
|E5|E7 root5: phi|0|(1+3s)/4|(2,4),(2,5),(3,5)|
|E*|E7 root5: 1|phi|(1+s)/4|(3,5)|
|F0|E8 root0: phi|0|4+2s|(0,2)|
|F6|E8 root6: phi|0|1+s|(3,5),(3,6),(4,6)|

All these survivors have zero contact with the whole D5. This is an
output of its complete inclusion in the two finite corners, not an input
restriction. No extra E8 type survives its A1 spectator. The E7 type E*
is genuinely new to the necessary domain: its height is not obtained by
pretending that the A1 is orthogonal.

Contact vectors plus positive height determine a unique outward
orientation of each root. There are eleven candidate third walls in
total: six D5, three E7 and two E8. If a candidate wall occurs, its checked
finite corners are actual vertices: projection of a strict interior
timelike vector onto the orthogonal complement of a positive principal
wall span preserves the other strict inequalities, since the inverse
of a positive Gram with nonpositive off-diagonal entries is nonnegative.
The two corresponding face edges and their finite endpoints close a
triangle by convexity. Thus identical candidate walls appearing in
several marked domains really account for the listed triangle bundles.

## Coupling the eleven candidates

For any two candidates set p=p_i dot p_j and Q=H_i H_j. Their actual
mutual product is p-sqrt(Q). Their old projections are orthogonal across
different old components, except that the two A1-contact types share that
one projection. Within one component its exact inverse is used.

All55 unordered pairs are checked. Eighteen have p>0 and p^2>Q, so their
actual product is positive and forbidden by acuteness. Twenty-three
have p=0 and Q nonintegral in Z[phi], violating the length-two cyclic
integrality condition. Two have p>0 and p^2=Q and therefore are
orthogonal necessary pairs. The remaining twelve have p=0 and Q an
algebraic integer with either an explicit allowed finite bond or Q>=4.
Those twelve are only necessary possibilities, not assertions that
the walls simultaneously extend to a polytope.

For a+bs in k, integrality means that 2b and a-b are integers. This is the
literal test used. For the surviving finite-bond cases, the code checks
the known values1,2,phi squared,3,(5+s)/2, which are respectively the
squared root products for labels3,4,5,6,10. It does not assume arbitrary
integral squares below4 are native.

The existing theorem already makes the six D5 types pairwise
incompatible. Within E7, E5 is incompatible with E0 and E*, while E0 and
E* are orthogonal. Thus E7 contributes at most three listed nonadjacent
triangles. The two E8 types are orthogonal and together contribute four.
If a D5 triangle occurs, all E7 types except E0 are forbidden; the total
then is at most1+1+4=6. With no D5 triangle the total is at most3+4=7.

For example, E* has an A1 projection in common with the new D5 type C*.
Their projection product is (3+s)/4 and the square root of their height
product is1/2. Their actual product is (1+s)/4>0, excluding coexistence.
The other D5/E5 and D5/E* exclusions have zero projection product and
nonintegral height product. The first one cannot be inferred from
different-component orthogonality: its shared A1 must be included.

The code independently enumerates all2^11=2048 root subsets consistent
with these necessary pair gates, counts the union of their marked
triangle faces, and finds maximum seven with the unique necessary
maximizing subset {E5,F0,F6}. This verifies the last finite combinatorial
step. No assertion of global realization of that subset is made.

## What this does and does not change

The new result constrains actual incidence at a repeated supported
rank21 forest. It also describes the complete necessary equality case,
which may be useful when classifying adjacent quadrilaterals or edge
endpoints. It gives no reduction in an old per-vertex charge allowance
by itself. The established global rule needs at least two total payment
on each triangle and at least one on each quadrilateral; a nontriangle
may still be a charged quadrilateral. No weight LP was solved, no new
charging column was inserted, and no old unit-payment consumer changed.
The compact arithmetic ceiling remains21, not20.

Generation completed exit0/PASS in terminal chunk `d95b4b`. The process
used about29.7MiB peak RSS (31,113,216 bytes), the shared exclusive proof
flock, a59-second wall alarm, and a periodic one-GiB RSS guard. It used
standard-library exact arithmetic, no optimizer, NumPy, Lean or parallel
worker. No temporary proof outputs or bytecode cache were produced.
The certificate is418,007 bytes of saved exact data plus a small Python program,
not an unbounded catalogue.

```sh
python3 -B scripts/fresh_compact21_supported_forest_triangle_coupling_20261018.py --verify
```

Separate-process complete saved replay passed in chunk `cff80c`, exit0,
in0.515 seconds, with measured peak RSS30,900,224 bytes. Independent
mathematical review is pending at this writing. The full geometric proof has not been encoded in Lean; no Lean
or external human verification is claimed.

Producer SHA-256:
`9637283f5f1c83c732f3b3ca84c0df3b855cdd564eb8bdd04c0c9d1992f60e9b`.
Certificate SHA-256:
`d663dc0650a6e2681b5884119727b41b2ca0770f7b0c871187f6194eda77b5c6`.

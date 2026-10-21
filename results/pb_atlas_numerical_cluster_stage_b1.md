# PB-Atlas numerical Stage B1: the arrays-only cluster compiler

## Result

The second bounded PB-Atlas layer is implemented.  Starting from **one**
frozen Stage-A endpoint program, with no surface name, Coxeter label, stored
geometric point, or classified incidence pattern, the compiler now produces:

1. the chronological reverse-blowup program;
2. the exact Enriques proximity DAG, including free and satellite centres;
3. the endpoint degree/bidegree and multiplicity row of every supplied wall;
4. every proximity inequality and excess;
5. exact intersection, genus, Euler-characteristic and virtual-dimension data;
6. the class matrices for `K`, `-K`, `-2K`, and the weighted boundary;
7. the standard-section columns and local-jet rows of the `-K` and `-2K`
   evaluation matrices.

The compiler is
[`scripts/pb_atlas_cluster_stage_b1.py`](../scripts/pb_atlas_cluster_stage_b1.py).
Its schema, hostile fixtures, certifier, frozen archive, and tests are:

- [`schemas/pb-atlas-numerical-stage-b1.schema.json`](../schemas/pb-atlas-numerical-stage-b1.schema.json);
- [`data/fixtures/pb_atlas_numerical_stage_b1_programs.json`](../data/fixtures/pb_atlas_numerical_stage_b1_programs.json);
- [`scripts/freeze_pb_atlas_stage_b1_fixtures.py`](../scripts/freeze_pb_atlas_stage_b1_fixtures.py);
- [`scripts/certify_pb_atlas_numerical_stage_b1.py`](../scripts/certify_pb_atlas_numerical_stage_b1.py);
- [`data/enumerations/pb_atlas_numerical_cluster_stage_b1.json`](../data/enumerations/pb_atlas_numerical_cluster_stage_b1.json);
- [`tests/test_pb_atlas_numerical_cluster_stage_b1.py`](../tests/test_pb_atlas_numerical_cluster_stage_b1.py).

The frozen Stage-A core was **not edited**.  The Stage-B1 archive checks its
SHA-256 against the Stage-A archive before compiling anything.

## 1. Exact reverse-blowup extraction

Let the Stage-A contraction order be

\[
                 C_1,C_2,\ldots,C_N,
\]

and let `E_1,...,E_N` be the orthogonal total-exceptional basis emitted in
that same blowdown order.  Reversing the birational map gives chronological
blowup centres

\[
                 p_1,p_2,\ldots,p_N
                 \quad\leftrightarrow\quad
                 E_N,E_{N-1},\ldots,E_1.
\]

After reversing both rows and columns, the contracted strict curves give an
upper-triangular matrix `P` with diagonal one.  Its only other permitted
entries are `-1`.  Define

\[
 Q_{cp}=1
 \quad\Longleftrightarrow\quad
 p_c\text{ is proximate to }p_p .
\]

The compiler proves coefficientwise

\[
                 P=I-Q^{\mathsf T},
 \qquad
                 (C_i.C_j)=-PP^{\mathsf T}.                 \tag{1.1}
\]

Every parent precedes its child, every centre has zero, one, or two
exceptional parents, and a claimed two-parent satellite is accepted only if
the two parent strict exceptionals intersect once immediately before that
blowup.  Thus a three-parent centre and a two-parent centre on disjoint
exceptionals are rejected, not silently coerced into a free cluster.

The invariant object is the DAG/partial order.  A raw sequential witness is
still retained so the next layer can build charts.  Independent centres may
commute; the archive checks an explicit centre/wall bijection between two
transported witnesses, verifies `Q` and every divisor coefficient under that
bijection, and obtains the same relabelling-invariant cluster profile.
Source-wall reversal and the unoriented `F_0` ruling swap are tested
separately.

## 2. Divisor and proximity rows

For a plane endpoint the chronological coordinates are written

\[
             D=dH-\sum_i m_iE_i .
\]

For a ruled endpoint they are

\[
             D=aS+bF-\sum_i m_iE_i,
 \qquad S^2=-e,\quad S.F=1,\quad F^2=0 .
\]

Every supplied wall receives the following exact row:

- `d`, or the section/fibre bidegree `(a,b)`;
- all chronological multiplicities `m_i`;
- `D^2` and `K.D`;
- arithmetic genus `1+(D^2+K.D)/2`;
- `chi(O(D))=1+(D^2-K.D)/2`;
- Riemann--Roch virtual projective dimension `(D^2-K.D)/2`;
- the endpoint section count, when the standard endpoint class is effective;
- the ordinary cluster-condition length `sum_i binom(m_i+1,2)`, when all
  multiplicities are nonnegative;
- for every centre `p`, the proximity excess

  \[
                   m_p-\sum_{q\succ p}m_q.                 \tag{2.1}
  \]

Negative multiplicities, negative excess, special cohomology, and dependent
conditions are recorded rather than erased.  In particular, expected
dimension is not promoted to `h^0`.

## 3. The `-K` and `-2K` matrices

The class matrix uses columns

```text
P2:  H | E_p1 ... E_pN
F_e: S | F | E_p1 ... E_pN
```

and freezes the four rows `K`, `-K`, `-2K`, and weighted boundary sum.  The
last two agree coefficientwise.

The evaluation-matrix template is also complete at the numerical level.
For `P2`, its columns are all homogeneous exponent triples of degree `3` or
`6`.  For `F_e`, the columns are the standard pushforward summands

\[
 H^0(F_e,aS+bF)
   =\bigoplus_{i=0}^{a}H^0(\mathbf P^1,\mathcal O(b-ie)),   \tag{3.1}
\]

with nonnegative summands only.  At each centre the `-K` matrix has the
value row; the `-2K` matrix has value and the two first-jet rows.  Hence the
row counts are `N` and `3N`, respectively.  Every column and row has an
explicit integer index and a formal entry token

```text
jet(center, local_multiindex)(endpoint_section_column).
```

The polynomial replacing that token depends on universal proper/free/
satellite chart coordinates.  Generating those entries is Stage B2.  B1
therefore asserts neither matrix rank nor nonemptiness.

## 4. Four frozen calibrations

Every one of the `282` Stage-A terminal orbit representatives compiles.  The
calibration aliases in this table live outside the recognizer.

| calibration | orbit programs | labelled endpoints | proper centres in representatives | free centres in representatives | proximity edges |
|---|---:|---:|---:|---:|---:|
| `H10/30` | 197 | 3,864 | 748 | 1,057 | 1,057 |
| `H11/26` | 69 | 440 | 224 | 475 | 475 |
| `H14/60` | 15 | 560 | 128 | 68 | 68 |
| `H17/48` | 1 | 1 | 16 | 0 | 0 |
| **total** | **282** | **4,865** | **1,116** | **1,600** | **1,600** |

The labelled endpoint histogram is unchanged from Stage A:

```text
F_0  1029      F_1  2586      F_2  564
F_3    44      F_4    12      P2   630
```

All four calibrations are free-only at the exceptional-proximity level.  No
satellite is inferred from them.  That absence is why the separate hostile
satellite fixture is essential rather than decorative.

Across the representatives the compiler freezes `8,652` wall-divisor rows.
The representative-centre total is `2,716`; weighting by terminal orbit size
gives `20,104` proper and `26,998` free centre occurrences.

## 5. Hostile and mutation tests

The synthetic arrays contain:

1. a free forest with two commuting proper roots;
2. the same forest with the independent order transported;
3. a genuine satellite `p3` proximate to both `p1` and `p2`;
4. a ruled `F_2` free cluster;
5. a forbidden three-parent centre;
6. a forbidden satellite whose two exceptional parents are disjoint.

They verify the exact satellite matrix

\[
 Q=\begin{pmatrix}0&0&0\\1&0&0\\1&1&0\end{pmatrix},
 \qquad
 P=\begin{pmatrix}1&-1&-1\\0&1&-1\\0&0&1\end{pmatrix},
\]

including the parent-intersection check.  The fixtures also include one
proximity-consistent divisor and an explicit negative-excess control.

Mutation coverage includes:

- reversal of the complete source wall array in each of the four real cases;
- exchange of the two `F_0` rulings;
- exchange of commuting reverse-blowup centres with all divisor columns
  transported;
- rejection of an added surface/type label;
- rejection of malformed strict/total exceptional data.

The test command is

```bash
python3 -m pytest -q tests/test_pb_atlas_numerical_cluster_stage_b1.py
```

and currently gives `8 passed`.  It recompiles and schema-validates all 282
real packets, checks every frozen packet/divisor-table digest, and replays the
complete archive.

## 6. Exact scope boundary and next gate

This result closes only **PB-Atlas numerical Stage B1**.  It does not prove:

- a universal-coordinate chart or a polynomial evaluation entry;
- incidence-stack nonemptiness or independence of conditions;
- collision, Jacobian, irreducibility, Hilbert, or avoidance bad factors;
- reduced/disjoint bi-anticanonical boundary geometry;
- hidden-wall exclusion, full-nef generation, same-P return, or L4.

The smallest next layer is not another numerical census.  Stage B2 should
consume the frozen `center_records`, `Q`, divisor rows, and evaluation row/
column indices; construct standard `P2/F_e` gauges and inductive free/
satellite blowup charts; replace each formal jet token by its polynomial;
then generate the finite incidence and avoidance factors needed to prove a
surviving strong-Coble/full-nef carrier and run the literal same-P replay.
The separate maximal-carrier theorem can handle the ensuing finiteness step;
no additional uniform facet/Hilbert-degree bound is part of B1.

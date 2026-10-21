# H9/14 record 16: complete solver-independent exclusion

## Theorem and scope

The second-layer graph shadow of record 16,

```text
M???????????o?o?_
```

cannot be the nerve of a compact hyperbolic Coxeter 9-polytope.  More
precisely, there is no simplicial 8-sphere on these fourteen vertices that
satisfies the published compact-Coxeter missing-face reductions, all of its
polytopal link constraints, and the Dehn--Sommerville schedules.

This is a theorem about the **compact** H9/14 search.  It does not assert
finite or infinite volume for any noncompact Coxeter chamber, and it does not
classify Coxeter labels or Lorentzian Gram matrices.

The proof is independent of the archived Z3 run, which timed out with status
`unknown`.  It uses blocker duality, every ridge equation, an analytic graph
bound, exact small-graph enumeration, and a second independently written
hostile replay.  The resulting graph-shadow partition is

\[
122\text{ independently impossible}
+7\text{ explicit ordinary-polytopal}
+13\text{ theorem-undecided}.
\]

## 1. Shadow, schedules, and local-link support

Put

\[
U=\{0,1\},\qquad I=\{2,\ldots,11\},\qquad P=\{12,13\},
\qquad W=U\sqcup I.
\]

The missing-edge graph has the five edges

\[
\{0,12\},\{1,12\},\{0,13\},\{1,13\},\{12,13\}.
\]

The edge `\{0,1\}` is the deleted edge of this `K4-e` on `U\cup P`; the
other ten vertices `I` are isolated.  Its inclusion-minimal vertex covers
are

\[
P,\qquad U\cup\{12\},\qquad U\cup\{13\}.
\tag{1.1}
\]

Let `a` be the number of missing four-faces.  The exact schedules are:

| `a` | missing-five offset | facets | summed vertex-link `g4` |
|---:|---:|---:|---:|
| 0 | 8  | 110 | 42 |
| 1 | 14 | 108 | 36 |
| 2 | 20 | 106 | 30 |
| 3 | 26 | 104 | 24 |
| 4 | 32 | 102 | 18 |
| 5 | 38 | 100 | 12 |
| 6 | 44 | 98  | 6  |

Every one of the fourteen vertex links has strict `g4>=1`; hence the last
two rows, whose sums are below 14, are impossible.  The graph-determined
baseline link vectors are

\[
g_3=(4,4,6^{10},1,1),\qquad
g_4=(5,5,7^{10},1,1).
\tag{1.2}
\]

The exact local-link intervals at vertices 12 and 13 allow no missing
four-face through either vertex.  A four-set of isolates would itself be a
missing four-face in both special vertex links and would spend their sole
baseline `g4`; therefore every missing four-face lies on `W` and meets
`U`.  The same strict special-link bound forces every additional minimal
five-face to lie on `W`.  The mechanical support replay checks all 735
independent four-sets and all 1,212 independent five-sets; its relevant
counts are 495 four-sets on `W`, 210 isolate-only four-sets, 792 five-sets on
`W`, and 420 five-sets meeting `P`.

## 2. Blockers and the complete ridge equations

Let `H` be the clutter on `W` consisting of the missing four-faces and the
additional minimal five-faces, and let `Z` be the members of `H` avoiding
`U`.  Because a facet complement is an inclusion-minimal transversal of the
complete minimal-nonface clutter, and no transversal may have order at most
four,

\[
\tau(H)\ge3,\qquad \tau(Z)\ge2.
\tag{2.1}
\]

Using the three graph covers (1.1), every five-element facet complement is
of exactly one of two forms:

\[
P\cup T \quad (y_T=1),\qquad
U\cup\{12\}\cup S,\ U\cup\{13\}\cup S \quad (z_S=1),
\tag{2.2}
\]

where `T` is a three-transversal of `H` and `S` is a two-transversal of `Z`
contained in `I`.  The two `z` complements occur as a pair.  Condition
`tau(Z)>=2` also gives `y_{01i}=0` for every isolate `i`.  Thus, in the row
with `a` missing four-faces,

\[
|y|+2|z|=110-2a.
\tag{2.3}
\]

The complement of a ridge is a six-set, and the number of selected
five-complements inside it is zero or two.  Enumerating all
`binom(14,6)=3003` six-sets gives exactly the following nonempty equations:

\[
\sum_{T\in\binom Q3}y_T
+2\,\mathbf 1_{Q=\{0,1,i,j\}}z_{ij}\in\{0,2\}
\quad(Q\in\tbinom W4),
\tag{R1}
\]

\[
z_{ij}+z_{ik}+z_{jk}\in\{0,2\}
\quad(i,j,k\in I).
\tag{R2}
\]

There are 735 nonempty ridge groups: 495 of type (R1), including 45 with a
doubled `z`, and 240 of type (R2).  The replay checks all 1,980 `y`
incidences and all 810 `z` incidences, so these are not selected necessary
equations but the entire ridge system.

Equation (R2) forces the `z` graph on the ten isolates to be a complete
bipartite cut

\[
z=K_{r,10-r},\qquad 0\le r\le5.
\tag{2.4}
\]

An exhaustive seed replay reconstructs all 512 labelled solutions, with
edge-count histogram

\[
0:1,\quad9:10,\quad16:45,\quad21:120,\quad24:210,\quad25:126.
\]

## 3. The `X` graph and the analytic pruning

Define a graph `X` on the isolates by `x_ij=y_{0ij}`.  Equation (R1) gives

\[
y_{1ij}=x_{ij},\qquad
y_{ijk}=x_{ij}\mathbin\oplus x_{ik}\mathbin\oplus x_{jk}.
\tag{3.1}
\]

Its exact part says that `X` is triangle-free and induced-`2K2`-free, and
that `x_ij=0` on every `z` edge.  Edges on both sides of the cut would form
an induced `2K2`, so all edges of `X` lie on one side, of size at most
`s=10-r`.

If `e=|E(X)|` and the nonzero degrees are `d_i`, direct counting in (3.1)
gives

\[
|y|=12e-\sum_i d_i^2\le 9s.
\tag{3.2}
\]

The inequality is Cauchy--Schwarz.  Combining it with (2.3) eliminates all
cuts for `a=0,1,2`; the `a=0` details are separately frozen in the branch
certificate.  Only five cases remain:

| `a` | cut `r+(10-r)` | required `|y|` | result of exact `X` replay |
|---:|---:|---:|---|
| 3 | 2+8 | 72 | no graph |
| 3 | 3+7 | 62 | `FCxv?`, degrees `(3,3,3,3,3,3,2)`, no A4 candidate |
| 4 | 2+8 | 70 | `G?o~f_`, degrees `(4,4,3,3,3,3,3,3)`, no A4 candidate |
| 4 | 3+7 | 60 | `F?zv_`, `F?~v_`, only two A4 candidates each |
| 4 | 4+6 | 54 | `EFz_ = K3,3`, exactly four A4 candidates |

For the four discarded rows, the number of possible four-blockers is
strictly below the scheduled `a`.  The last graph has degree sequence
`(3,3,3,3,3,3)` and is structurally verified to be `K3,3`.  Its four
four-blockers are

\[
\{0,2,3,4\},\ \{0,5,6,7\},\
\{1,2,3,4\},\ \{1,5,6,7\}.
\tag{3.3}
\]

The small-graph input is an archived exact `nauty geng 2.9.3` enumeration
of all unlabeled triangle-free graphs:

| order | command | graphs | retained after independent induced-`2K2` test |
|---:|---|---:|---:|
| 6 | `geng -q -t 6` | 38 | 22 |
| 7 | `geng -q -t 7` | 107 | 41 |
| 8 | `geng -q -t 8` | 410 | 82 |

The downstream verifier decodes every graph6 word, rechecks
triangle-freeness, performs an exact pairwise nonisomorphism test, tests the
induced-`2K2` condition, evaluates (3.1) directly, and constructs the A4
blockers.  Thus `geng` is trusted only for exhaustive generation, not for
the retained properties.

The raw graph6-stream SHA-256 digests in orders 6, 7, and 8 are respectively
`f0da23b1ee4c292a4e77422a38b5fe2b3774997a541249c9a3e93f5a3f4e8ae5`,
`c21188cdee83a78eb19f73a67d22fe81fbb8349d52e99d1373e6322110baea15`,
and `25ef89391c8de758e689783e1b930edd36e5fe8c2c4963e4e9821f6b44eee829`.

## 4. The final `K3,3` case

The four A4s (3.3) force a shadow of 30 five-element nonfaces.  The `a=4`
schedule requires offset 32, so exactly two further minimal A5s are needed.
The executable enumeration considers 762 candidate minimal A5s and all
289,941 unordered pairs.  Exactly 126 pairs satisfy `tau(Z)>=2`; in every
one, both A5s avoid `U`, are disjoint, and partition the ten isolates.
Consequently their two-row blocker is always `K5,5`, with 25 `z` edges.
But the complete ridge equations fixed `z=K4,6`, with 24 edges.  This is a
contradiction.

There is also a shorter, stronger independent check.  Construct the 102
facet complements directly from the ridge-forced `K4,6` and `X=K3,3`, then
dualize this full facet family over all `2^14` subsets.  Its complete list of
minimal transversals has size histogram

\[
2:5,\qquad4:4,\qquad5:2,\qquad6:2,\qquad7:2.
\tag{4.1}
\]

The order-two and order-four entries are precisely the five prescribed
missing edges and the four A4s (3.3).  The order-five entries are

\[
\{8,9,10,11,12\},\qquad\{8,9,10,11,13\},
\]

which violate the already proved special-link support restriction.  The
four entries of orders six and seven independently violate the compact
Coxeter/Lanner bound that every minimal nonface has order at most five.
Thus the final facet system is impossible even without relying on the
two-A5 schedule argument.

This exhausts `a=0,...,6` and proves the theorem.

## Trust and reproducibility

Frozen artifacts:

- full manifest:
  `data/enumerations/h9_14_record16_certificate.json`;
- canonical full-certificate SHA-256:
  `ac7f386a36ffe8a0ddbeab74c0fdd60913eb891365735588ee7a111fc8b7c8d5`;
- full-manifest file SHA-256:
  `2c4377089da7487c1d0694d5abd9759e6ac1ca5dde75a66b209fc46f459b8888`;
- positive-branch proof-module SHA-256:
  `6fd9d9842db6b5e907a9e5f65a4f5c04dd9b0485b4d71a9682167075993abe19`;
- hostile-replay source SHA-256:
  `2961fe9223689d744e2a2ab8a55e92459bb3520f21aef5859a0ac0442c1f3139`;
- hostile summary SHA-256:
  `bfd07dd04a18632fa7a8ca1c868e48fed60d8b944e594282b83427241633c3d5`;
- dependency-closure source SHA-256s (`local_link.py`, `missing_edges.py`,
  `h9_14_record92.py`, `nerve_support.py`, and the chained A0 certifier):
  `8f8169edab8b87aff682785667c2a8eab5c52b92a7c68cef130aeec4d3ae807d`,
  `231e85fa98e27e864f997cd4ebe59136707f6aa8dd2d8136b80cff543dfbfbb8`,
  `6b165e2a9e3ed120949c10784e42786364a7197014bc57da6b2b0a1cb30e1e29`,
  `aad8fb25c2f9a56e45121fa1d1cbbdedfdff8ba1f67e379be74eed73940ad2a2`,
  and `9b032b4afce9a76bf67bd06ad52374b2fcaf8c28dcc74f1333bae7c23c92f9c2`;
- X-archive file / internal SHA-256:
  `1d6915efe97c2ec5fc81e0022d415c1073e17bbed0e012038a1171132630a1ce` /
  `3a3af277fb725b44b136990e7096073dbfd404a3291d0f8bba9a7c725fda7479`;
- X-archive generator source SHA-256:
  `d9b6457e483e24e0fe78ab89e7f10a5d8f802d9536dcfd5731d4bb83f8f7171f`.

Replay with

```text
PYTHONPATH=src:. python3 scripts/certify_h9_14_record16.py --verify-only
PYTHONPATH=src:. python3 scripts/independent_h9_14_record16_audit.py
PYTHONPATH=src:. pytest -q tests/test_h9_14_record16_a0_certificate.py tests/test_h9_14_record16_certificate.py
```

The archived bounded Z3 result is explicitly marked
`trusted_for_proof: false`; its timeout contributes no logical step.

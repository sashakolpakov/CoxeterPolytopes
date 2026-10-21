# Solver-independent exclusion of H9/14 record 46

> **Current status.** The later record-37 and record-34 certificates and
> record-35 rational multiwedge give the current partition (121+6+15); see
> `results/h9_14_facet_completion.md`.

## Result

**Record 46 is impossible.** Its graph6 string is

```text
M?????????O?o?o?_
```

and its missing-edge graph has edge set

\[
\{(0,11),(0,12),(1,12),(0,13),(1,13),(12,13)\}.
\]

The proof combines strict local `g`-vector inequalities with an exact
symmetry-compressed enumeration of four missing five-faces and blocker
duality. It uses no SAT/SMT backend.

**Historical stage snapshot.** Together with the previous certificates, this
record-46 result changed the theorem-facing partition of the 142 graph
shadows to

\[
119\text{ independently impossible}
\; +\;5\text{ explicit ordinary-polytopal selections}
\; +\;18\text{ theorem-undecided}.
\]

This excludes one graph shadow, not all compact H9/14 nerves.

## 1. Strict link sums leave only `a=0,1`

Let `a` be the number of missing four-faces. The four
Dehn--Sommerville schedules give

| `a` | sum of vertex-link `g4` |
|---:|---:|
| 0 | 24 |
| 1 | 18 |
| 2 | 12 |
| 3 | 6 |

Every one of the fourteen vertex links must have strictly positive `g4`.
The sum must therefore be at least 14, excluding `a=2,3` immediately.

Before selecting missing four- or five-faces, the exact vertex-link entries
are

\[
g_3(v)=g_4(v)=(1,2,4,4,4,4,4,4,4,4,4,3,1,1)_v.
\tag{1.1}
\]

## 2. The `a=1` branch is locally impossible

There are 690 independent four-sets. Applying every exact vertex-, edge-,
and triple-link incidence interval leaves 294 possibilities in three types:

| active part of A4 | multiplicity | forced independent A5 shadow | vertices with nonpositive `g4` |
|---|---:|---:|---|
| empty | 126 | 10 | `0,2,3,4,5,12,13` for a representative |
| `{1}` | 84 | 8 | `0,1` |
| `{11}` | 84 | 9 | `11,12,13` |

The displayed vertices are computed after inserting the A4 and all of its
forced missing-five extensions, before adding any further minimal A5.
Additional missing five-faces can only decrease the relevant `g4` entries.
Thus every locally admissible A4 already contradicts strict positivity, and
`a=1` is impossible.

## 3. Exact four-A5 census for `a=0`

The remaining schedule has exactly four distinct minimal A5s and requires

\[
f_8=90.
\tag{3.1}
\]

Equation (1.1) and the exact local intervals force every A5 to avoid

\[
\{0,12,13\}.
\]

Hence its support lies in

\[
\{1,11\}\sqcup I,
\qquad I=\{2,3,\ldots,10\}.
\]

The nine vertices in `I` are isolated in the missing-edge graph and are
fully symmetric. The pair-link interval for `{1,11}` has A5 incidence upper
bound zero, so no selected A5 contains both distinguished vertices. There
are 378 individually allowed independent A5s before the simultaneous local
constraints are imposed.

The finite compiler represents a family of four A5s by four incidence rows.
The columns for vertices `1` and `11` remain distinguished. Each isolate has
a four-bit column recording the rows in which it occurs, and the nine
isolate columns are stored only through their 16 pattern multiplicities.
This is an exact quotient by the full isolate permutation group; a further
canonicalization quotients the four row permutations.

For every state the audit checks:

- every row has size five and the four rows are distinct;
- every row is independent in the missing-edge graph;
- all exact A5-incidence intervals in every vertex, edge, and triple link;
- absence of transversals of sizes one through four;
- the number of five-element transversals of the complete clutter consisting
  of the missing edges and the four A5s.

The absence of smaller transversals makes every five-transversal minimal.
By blocker duality, their complements are exactly the candidate
eight-dimensional facets, so their number is `f8`.

The exhaustive census is

| stage | states |
|---|---:|
| active/isolate row-count joins | 55,458 |
| locally valid before distinctness | 55,458 |
| four distinct rows | 53,902 |
| canonical local states | 2,523 |
| no transversal of size at most four | 330 |

Across the final 330 classes, the five-transversal count ranges from 48 to
83. Its maximum is 83, attained by three classes. No class has the required
count 90. This excludes `a=0` and completes the proof.

## 4. Reproduction and trust ledger

Run

```sh
PYTHONPATH=src:. python3 scripts/certify_h9_14_record46.py --verify-only
PYTHONPATH=src:. python3 scripts/certify_h9_14_record46.py --verify-only --recompute
PYTHONPATH=src:. python3 -m pytest -q tests/test_h9_14_record46_certificate.py
```

The proof object is
`data/enumerations/h9_14_record46_certificate.json`. Its frozen canonical
payload digest is

```text
5f28acdd31ae4c043e799c4dfb62c63dd4b22f95af7032c6b8af83a5fa04b3f3
```

The full, uncompressed census has canonical digest

```text
2cbafdc1d17c6690aa1aa92b2c22d943245f23053c203ffc43b64070097a6af0
```

The proof trusts Burcroff's no-missing-triples theorem for the compact
H9/14 problem, the upstream strict-link and Macaulay inequalities, the
Lanner order-five bound, Dehn--Sommerville, and blocker duality. The archived
Z3 computation ended in a bounded timeout; it is recorded only as untrusted
corroboration and contributes no logical step. The optional `--recompute`
path replays the same exact census implementation; the separately reported
hostile audit independently checked its state quotient, branches, and counts.

# Solver-independent exclusion of H9/14 record 74

> **Current status.** The later record-46, record-37, and record-34 certificates and the
> multiwedge
> realizations of records 43, 44, 47, and 35 give the current partition
> \(121+6+15\); see `results/h9_14_facet_completion.md`.

## Result

**Record 74 is impossible.**  Its graph6 string is

```text
M?????????O?_@_@_
```

and its missing-edge graph is

\[
H=K_4[A]\sqcup10K_1,
\qquad A=\{0,11,12,13\},\quad I=\{1,\ldots,10\}.
\]

The proof uses strict local `g`-vectors and blocker duality.  It needs no
SAT/SMT backend and no enumeration of missing-five families.

**Historical stage snapshot.** Together with the previous certificates, this
changed the theorem-facing partition of the 142 graph shadows to

\[
118\text{ independently impossible}
\; +\;2\text{ explicit ordinary-polytopal selections}
\; +\;22\text{ theorem-undecided}.
\]

This excludes one graph shadow, not all compact H9/14 nerves.

## 1. Strict link sums leave `a=0,1`

Let `a` be the number of missing four-faces.  The four
Dehn--Sommerville schedules are

| `a` | missing independent five-sets | facets `f8` | sum of vertex-link `g4` |
|---:|---:|---:|---:|
| 0 | 4  | 90 | 24 |
| 1 | 10 | 88 | 18 |
| 2 | 16 | 86 | 12 |
| 3 | 22 | 84 | 6  |

Strict positivity of `g4` in all fourteen vertex links requires the last
column to be at least 14.  Hence `a=2,3` are impossible.

## 2. Every missing four-face is locally impossible

Before selecting missing four- or five-faces, the exact vertex-link values
are

\[
g_3(v)=g_4(v)=
\begin{cases}
1,&v\in A,\\
4,&v\in I.
\end{cases}
\tag{2.1}
\]

An independent set contains at most one vertex of the `K4`.  Thus the 690
independent four-sets form exactly two `S4 times S10` orbits:

| A4 type | multiplicity | independent five-shadow |
|---|---:|---:|
| four isolates | `C(10,4)=210` | 10 |
| one active plus three isolates | `4 C(10,3)=480` | 7 |

If an A4 `Q` contains an active vertex `v`, then the actual local entry is

\[
g_3(v)=1-1=0,
\]

contradicting strict positivity.  Equivalently, the exact active-vertex A4
incidence interval is `[0,0]`.

It remains only to consider `Q subset I`.  Every `Q union {v}`, `v in A`, is
an independent five-extension forced to be missing.  Therefore at every
active vertex

\[
g_4(v)=1-1=0.
\]

The four vertices of `Q` also have `g4=4+6-10=0`.  Thus this second orbit is
locally impossible as well, before any additional minimal A5 is chosen.

There is also a separate blocker check.  The isolate A4 has ten independent
five-extensions, exactly the scheduled number of missing five-sets for
`a=1`, so there is no additional minimal A5.  Any minimal three-cover `C` of
the `K4`, together with any `q in Q`, is then a transversal of order four.
This independently contradicts purity.

Hence the whole `a=1` branch is impossible.

## 3. The `a=0` facet count is divisible by four

The remaining schedule has four distinct minimal A5s.  If an A5 contains an
active vertex `v`, equation (2.1) gives

\[
g_4(v)=1-1=0.
\]

Consequently all four A5s lie in the ten-set `I`.  The audit checks all 1,092
independent five-sets directly: 840 meet `A` and fail this strict gate, while
the `C(10,5)=252` isolate-only sets are the complete locally allowed list.

The inclusion-minimal vertex covers of `K4` are exactly

\[
C_i=A\setminus\{i\},\qquad i\in A.
\tag{3.1}
\]

Let `G` denote the family of four isolate-only A5s.  If their total
intersection contains `q`, then `C_i union {q}` is a transversal of order
four.  Purity therefore forces

\[
\bigcap_{F\in G}F=\varnothing.
\tag{3.2}
\]

Under (3.2), the five-element transversals have an exact disjoint
classification.

1. A transversal with exactly three active vertices is a unique `C_i`
   together with a two-subset of `I` hitting every member of `G`.  Every
   hitting isolate pair occurs with all four choices of `i`.
2. A five-set with all four active vertices has one isolate `q`; it hits `G`
   precisely when `q` lies in the intersection (3.2), so no such set
   survives.

The machine audit independently reconstructs all 190 five-element graph
covers: 180 have type three-active-plus-two-isolates, with uniform
multiplicity four for each of the 45 isolate pairs, and 10 have type
four-active-plus-one-isolate.  Hence, writing `tau2(G)` for the number of
hitting isolate pairs,

\[
f_8=4\,\tau_2(G).
\tag{3.3}
\]

But the schedule requires `f8=90`, and

\[
90\equiv2\pmod4.
\]

This contradiction excludes `a=0` and completes the proof.

## 4. Reproduction and trust ledger

Run

```sh
PYTHONPATH=src python3 scripts/certify_h9_14_record74.py --verify-only
PYTHONPATH=src python3 -m pytest -q tests/test_h9_14_record74_certificate.py
```

The proof object is
`data/enumerations/h9_14_record74_certificate.json`.  Its canonical payload
digest is

```text
0a8812c4dbade45cabe3c544c968d265581e4d3f0fdd6b48673ff7847148ceb4
```

The audit reconstructs the missing graph, its four minimal covers, every
schedule and summed link value, all 690 independent A4s, both exact A4
orbits, all 1,092 independent A5s, and both classes of five-element graph
covers.  It trusts the upstream strict-link theorem, the Lanner order-five
bound, Dehn--Sommerville, and blocker duality.

The archived Z3 run ended in a bounded timeout.  That record is stored only
as untrusted corroboration and contributes no logical step.

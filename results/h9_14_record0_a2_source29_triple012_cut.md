# Record 0, compact `a=2`, source 29: the triple-012 cut

## Outcome

Source orbit 29 remains **OPEN**.  This note proves neither a compact
incidence exclusion nor a positive Coxeter realization.  It does replace one
1,600-variable undifferentiated search by six exact symmetry-normalized
branches, and it freezes an exact scheduled-shadow family only one strict
triple unit from feasibility.

## Analytic strict-link cut

The source has missing four-sets `0123,0124` and stabilizer

\[
S_{\{0,1\}}\times S_{\{3,4\}}\times S_{\{5,\ldots,12\}},
\qquad |G|=2\cdot2\cdot8!=161{,}280.
\]

Its 1,600 candidate additional A5s form 22 `G`-orbits.  For the allowed
triple `T=012`, the pre-addition link vector is `(g1,g2,g3)=(1,1,1)`.
If `d_T` selected A5s contain `T`, strictness gives

\[
g_2(T)=1-d_T\ge1,
\]

so `d_T=0`.  This deletes the complete 28-element orbit with signature
`(2,1,0,2,0)`.

Every allowed six-set containing `T` is `012F`, where `F` is a three-set in
the free eight-block.  After `d_T=0`, exactly three candidate A5s can destroy
it:

```text
01F, 02F, 12F.
```

Writing `m_T` for the number of destroyed such six-sets, the last strict
entry is `g3(T)=1-m_T`, hence `m_T<=1`.  Therefore all selected A5s in the two
killer orbits (sizes 56 and 112) lie over at most one free triple.

There are now six exhaustive branches:

- zero support: all 168 killer candidates are absent, leaving 1,404
  undecided candidates;
- active support: normalize the unique free triple to `567`; its three
  candidates have masks `227,229,230`.  Their nonempty subsets have five
  orbits under the remaining `S2` swapping `0,1`.  Fixing the three fibre
  bits leaves the same 1,404 undecided candidates and remaining selection
  counts `20,20,19,19,18`.

The support stabilizer has order `2*2*3!*5!=2,880`.

## One-defect exact shadow seed

The frozen 21-A5 family realizes the exact missing-shadow schedule

```text
(M6,M7,M8) = (147,420,629).
```

Every vertex and edge strict link passes.  Exactly one triple fails:

```text
face 123: final g=(2,1,-1), selected-A5 incidence 2.
```

Thus its total strict penalty is one.  It is not a compact incidence
witness: its blocker histogram is `4:9, 5:99, 6:5`.

An exhaustive canonical neighbourhood replay on the **post-triple-012-cut
ground** finds:

| replacements from seed | post-cut exact-shadow families | strict-link families |
|---:|---:|---:|
| 1 | 596 | 0 |
| 2 | 515,354 | 0 |

At distance two, 167,656 families fail a vertex/edge/triple-degree cap and
347,698 reach the final triple-`g3` gate and fail there.  The smallest
remaining penalty is still one.  Consequently any strict-link scheduled
shadow seed differs in at least three selected A5s from the frozen family.
As a domain cross-check, the uncut radius-one ground has 599 exact-shadow
families: the 596 displayed above plus exactly three using the already-proved
forbidden `d_012>0` orbit.  Those three cannot be strict.  No all-ground count
is being mislabeled as a post-cut count, or conversely.

## Exact next step

The subsequent
[`triple-123 fibre-budget theorem`](h9_14_record0_a2_source29_triple123_cut.md)
refines these six branches into 144 exact outside-pair configurations with
killer-fibre budgets `4,2,1,0`.  Solve those bounded supports with the exact
order-6/7/8 shadows and all strict links.  Equivalently, search only repairs
of at least three A5 replacements from the frozen near seed.  A survivor must
then pass blocker purity, all 176 facets, ridge degree `0/2`, Coxeter labels,
and the Lorentzian Gram test.

The archive is
[`h9_14_record0_a2_source29_triple012_cut.json`](../data/enumerations/h9_14_record0_a2_source29_triple012_cut.json),
and the certifier is
[`certify_h9_14_record0_a2_source29_triple012_cut.py`](../scripts/certify_h9_14_record0_a2_source29_triple012_cut.py).
The default certificate audit is small.  The optional full neighbourhood
replay is guarded by an explicit 1 GiB peak-RSS ceiling:

```bash
PYTHONPATH=src:. python3 scripts/certify_h9_14_record0_a2_source29_triple012_cut.py --replay-radius-two
```

No timeout, unknown solver status, or bounded search disposition is used as
mathematical evidence.

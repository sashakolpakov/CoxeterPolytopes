# Bounded global lineage of `X2-RAT-C2-H10-37`

## Verdict

As of 2026-08-07, no published exact match for the literal 37-facet
diagram has been found.  This is **not** a first example in dimension 10,
not a new arithmetic lattice, and not a new maximal reflection group.  The
defensible priority statement is narrower:

> `X2-RAT-C2-H10-37` is a newly extracted genuine
> K3/rational-quotient divisor realization and a bounded literal-diagram
> priority candidate.  Worldwide priority of the abstract 37-facet
> polytope remains open.

The lineage separates as follows.

| layer | status |
|---|---|
| hyperbolic dimension $10$ | previously populated |
| lattice $U+E_8(-2)+\langle-2\rangle=(11,9,1)$ | known |
| full root-reflection group | known 13-wall finite-volume group |
| 37-wall group | proved finite-index subgroup, exact index $241{,}920$ |
| literal 37-node relation diagram | unmatched in the bounded audits below |
| extraction from a rational involution of Shimada's $X_2$ | not present in the audited primary source |

Thus the new catch does not fill a formerly empty dimension or arithmetic
class.  Its possible novelty is the larger fundamental chamber for a
finite-index reflection subgroup and the algebro-geometric route producing
it.

## Published denominator, but not the 37-wall numerator

[Shimada's primary paper](https://arxiv.org/abs/1412.6904), *The
automorphism groups of certain singular K3 surfaces and an Enriques
surface*, supplies the complete induced chamber $D^{(0)}$ for $X_2$,
the group

\[
\operatorname{Aut}(X_2,a_2)\cong \operatorname{PGL}_2(\mathbb F_9),
\]

and its 45 symplectic and 36 rational involutions (Proposition 1.6).  It
does not take an $X_2$ rational-involution fixed section or state a 37-wall hyperbolic
polytope.  The paper's analogous projected construction is instead the
Enriques quotient $Z_0=X_0/\langle\varepsilon\rangle$: Section 10 defines
$D_{Z}^{(0)}=\mathcal P(Z_0)\cap D^{(0)}$ and explicitly obtains 40 wall
vectors in orbits of sizes 30 and 10.  Consequently, Shimada is the
published algebraic denominator for the new extraction, not a publication
of the extracted `H10/37` chamber.

The arithmetic denominator is likewise known.  [Alexeev--Engel](https://arxiv.org/abs/2208.10383)
compute the full reflection diagrams for the 2-elementary lattices on the
$g=1$ line and prove finite covolume in the relevant cases; their
$(11,9,1)$ diagram has 13 walls.  [Alexeev's companion
classification](https://arxiv.org/abs/2209.09110) gives the broader
2-elementary reflective-lattice setting.  Independently, the maintained
[Scharlau--Kirschmer Lorentzian census](https://www.math.rwth-aachen.de/~Markus.Kirschmer/lorentz/index.html)
places the scaled-dual/sign-convention version at `List10`, entry 6.  The
project's exact Gram identification and Euler-characteristic quotient give

\[
W_{37}\subset W_r(11,9,1),\qquad
[W_r(11,9,1):W_{37}]=241{,}920.
\]

This closes the lattice and maximal-group lineage.  General finite-index
reflection-subgroup results such as
[Felikson--Tumarkin (2003)](https://arxiv.org/abs/math/0305093) and
[Felikson--Tumarkin (2007)](https://arxiv.org/abs/0705.0426) provide
structural constraints, but the audited papers do not list this particular
37-generator subgroup or chamber.

## Strengthened CoxIter exclusion

The earlier check used only the declared dimension and facet count.  The
new audit scans **every 37-node record**, ignoring its declared dimension,
in the two lanes of the
[pinned CoxIter repository](https://github.com/rgugl/CoxIter/commit/b800be48d8240fbabd2d1c9542e61054584ad7fe):

| lane | all records | 37-node records | declared dimensions |
|---|---:|---:|---|
| ordinary | 1,214 | 2 | both 18 |
| index-two doubling archive | 2,959 | 43 | 42 in dimension 3; 1 in dimension 13 |
| total | 4,173 | 45 | none declared dimension 10 |

More importantly, all 45 records are exact nonmatches before any reliance
on their headers.  The complete off-diagonal CoxIter relation-label
histogram of the candidate is

\[
\{0^{121},1^{55},2^{405},3^{15},4^{70}\},
\]

where labels $0,1,2,3,4$ mean parallel, generic ultraparallel,
$m=2,m=3,m=4$, respectively.  Neither ordinary record nor any of the 11
histogram classes represented by the 43 doubling records has this
histogram.  Since a relation-diagram isomorphism preserves the multiset of
all 666 off-diagonal labels, every one of the 45 records is rigorously
excluded.

CoxIter label 1 deliberately forgets the finer metric decorating an
explicit dotted ultraparallel edge.  That does not weaken the negative
result: equality of the finer diagrams would still force equality of these
coarser label histograms.

The frozen corpus identifiers are:

```text
CoxIter commit                    b800be48d8240fbabd2d1c9542e61054584ad7fe
ordinary archive SHA-256         ae04648b291c1e5f3d0e47d666c846dc805b244469a20486bca59674b4d171c3
doubling archive SHA-256         3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036
all 45 records, canonical SHA-256 b271caf9b12ec3a4a2a5038a43645cd555cd09ce347379e56de1f2b24e0935fa
```

Exact searches using the unusual index, Euler characteristic, and face
census also produced no relevant match.  Those searches are discovery
screens only.  Neither they nor the CoxIter exclusion establish worldwide
absence: an older diagram can be missing from CoxIter, encoded under
different invariants, or present in literature outside the audited citation
neighbourhood.

## Frozen conclusion

The correct ledger wording is therefore:

> Known $H^{10}$ arithmetic class and known 13-wall full group; newly
> extracted genuine K3/rational-quotient 37-wall finite-index chamber.  All
> 45 37-node records in the pinned CoxIter ordinary and doubling lanes are
> exact relation-histogram nonmatches.  Worldwide abstract priority remains
> open.

This is materially stronger than “no `37 10` header,” while making no
worldwide novelty claim.

## Reproducibility

- [bounded-lineage certificate](../data/enumerations/shimada_aut3singk3_x2_h10_37_global_lineage.json)
- [memory-light certifier](../scripts/certify_shimada_aut3singk3_x2_h10_37_global_lineage.py)
- [regression test](../tests/test_shimada_aut3singk3_x2_h10_37_global_lineage.py)
- [complete `H10/37` theorem and construction certificate](shimada_aut3singk3_x2_rational_involution_h10_37.md)

Replay with

```bash
python3 scripts/certify_shimada_aut3singk3_x2_h10_37_global_lineage.py \
  --check-archive
PYTHONPATH=src:. pytest -q \
  tests/test_shimada_aut3singk3_x2_h10_37_global_lineage.py
```

Frozen SHA-256 values:

```text
composite H10/37 certificate  12914344403b7d3d6b69a2866c4150f0a0420a64914e4bc8584fa372799599a7
bounded-lineage certificate  52b4800b07608804b7834aad65a4671b3255ef918c07779497742f11091a984d
bounded-lineage certifier    128be8c89c6ad1fcfb75b7294dcff83acdf7c80809b1203cc7ce885439d7de98
bounded-lineage test         4eca0a996994c6083905a3d30c179cfdf4c68dc75befd22c2c194a7ec9b94711
```

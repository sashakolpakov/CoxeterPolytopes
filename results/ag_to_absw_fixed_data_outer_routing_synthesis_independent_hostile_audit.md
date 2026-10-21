# Independent hostile audit: AG-to-ABSW fixed-data closure and outer routing

## Verdict

**PASS, in the repaired conditional and classification-importing scope.**
The six terminal synthesis documents now agree with the three frozen
owner/audit theorem pairs.  The downstream statement is exact:

> A finite set of native pairs \((k,m=n+1)\), together with maximality among
> discrete reflection groups or among arithmetic reflection groups, gives
> finitely many unmarked normalized maximal chambers in dimension at least
> two, using Nikulin's 1982 fixed-data reflective-lattice theorem.

This is a conditional ABSW-type closure.  It is not an AG-only proof, it does
not construct the missing exhaustive carrier compiler, and it does not prove
global ABSW.

The two outer controls also have their exact terminal quantifiers.  `X3` is
**one counterexample** to a universal assertion over the class of all
finite-volume exact-module projected-root inputs; it neither obstructs every
input nor settles a maximal-only conjecture.  The Borcherds `H21` chamber, by
contrast, is absolutely reflection-maximal and empties the **current complex**
rational/weak-log/full-trace/projected-root stack, while retaining a positive
characteristic-two mixed curve/Cremona carrier.

No owner theorem or synthesis byte was edited by this audit.

## 1. Frozen inputs

The independent replay pins all fifteen inputs individually.

### Six synthesis documents

| key | artifact | SHA-256 |
|---|---|---|
| `goal` | [`GOAL.md`](../GOAL.md) | `0b1736db46f8872ddf3868159c273620642a328ac94bba0488e806ae83f80b4d` |
| `plan` | [`RESEARCH_PLAN.md`](../RESEARCH_PLAN.md) | `d4533215663c8398fad35de9d3a764b8f6dcf1e6fdea103eda8a5fc34515f983` |
| `completion` | [canonical completion audit](canonical_chamber_four_property_completion_audit.md) | `f1c56a0ade1b37a498b759262b3309c97f96a54ae00900d11d4035a490a05e36` |
| `inverse` | [inverse carrier ledger](inverse_carrier_transfer_ledger.md) | `ec572e1ba7f0936f3e1ed1eb2c0837780b58b807c2af5009988a498c8e46b2f9` |
| `nikulin` | [Nikulin-to-AG ledger](nikulin_to_ag_translation_ledger.md) | `f016006bbe9424a757b949feb17ad73cd39c9f7303448e1265912159eb858237` |
| `synthesis` | [fixed-data/outer-routing synthesis](ag_to_absw_fixed_data_closure_and_outer_routing_obstruction.md) | `9c7fdf14ab62629ae1a47d04a6e7ede91a8605e86b936305d2f2b957a6ff77b8` |

### Owner synthesis certificate package

| artifact | SHA-256 |
|---|---|
| [certifier](../scripts/certify_ag_to_absw_fixed_data_outer_routing_synthesis.py) | `cb78b46f171f6fcc07c99e6139a4198d0f628be628b79e3317471b8179a57487` |
| [archive](../data/enumerations/ag_to_absw_fixed_data_outer_routing_synthesis.json) | `d0adf354e185f601d2e96ed94ef8225742cb07afa6d6790626f0550932ee5911` |
| [tests](../tests/test_ag_to_absw_fixed_data_outer_routing_synthesis.py) | `ff941226e19d85b50d881669840b15cd973e2b7dba2f3dc315f8239c0813162d` |

### Frozen theorem/audit evidence

| endpoint | owner SHA-256 | independent-audit SHA-256 |
|---|---|---|
| [finite-pair theorem](projected_k3_maximal_abstract_chamber_finiteness.md) / [audit](projected_k3_maximal_abstract_chamber_finiteness_independent_hostile_audit.md) | `59dd458acf14d18021f6c5583096da13f39c83115fc136d6ba91bcc850ce3603` | `edc663c0e7f35880123149d81d832812d5f28837c837df16eaa0f1daf725916a` |
| [`X3` theorem](projected_k3_universality_rank_saturation_no_go.md) / [audit](projected_k3_universality_rank_saturation_no_go_independent_hostile_audit.md) | `7ed5176c0c0a36e036c0e12dfa6acf1ee15d90d5ec7143ca5df59ebc3b7143d2` | `89296f04bac6b0b4ed68d9b86936a08ede515ba5bf9c847abeaadbc628926071` |
| [`H21` theorem](borcherds_h21_outer_routing_counterexample.md) / [audit](borcherds_h21_outer_routing_counterexample_independent_hostile_audit.md) | `fa899977f418d6fa93f7d5270e067c28a648c5e55f2636cb2b362e381534821d` | `8d705a619e9c1ffff3ad69765e90bb1a1a2a5d8713015ecca79382e8a1ccbc86` |

The audit reads these files directly and does not import the owner producer.
It separately cross-checks the owner's version-2 archive against the live
six-document and six-evidence hashes.

## 2. Finite native pairs really do close the maximal chamber count

For one finite-sided finite-volume arithmetic chamber, choose a preserved
native \(\mathcal O_k\)-lattice \(L\).  Its side group and complete native
reflection group satisfy

\[
                 W(P)\subseteq R(L)\subseteq O^+(L).
\]

Both endpoint groups have finite covolume.  Hence \(W(P)\), and therefore
\(R(L)\), has finite index in \(O^+(L)\); \(L\) is reflective in the sense
used by Nikulin.  For fixed totally real \(k\) and fixed \(m\ge3\), Nikulin
1982, Theorem 5.2.1, leaves finitely many such lattices **up to similarity**.
Similarity is the right equivalence because common scaling cancels from

\[
             x-\frac{2h(x,r)}{h(r,r)}r
\]

and from normalized Gram data.

Reflection maximality is the indispensable final hypothesis.  It gives
\(W(P)=R(L)\); then a finite generic gallery proves that the complete locally
finite mirror arrangement has one chamber orbit.  Without maximality,
finite-index reflection descendants on one fixed host can supply infinitely
many chambers.  The terminal synthesis documents now retain maximality in
every finite-pair summary.

This reasoning also identifies exactly where classification enters.  The AG
rank/field bounds may produce a finite pair set, but the passage from fixed
\((k,m)\) to finitely many native similarity hosts is Nikulin's theorem.  The
blanket claim that Nikulin is no longer a dependency is therefore false for
this downstream route and is absent from the terminal goal.

## 3. Field and ambient-host scopes

The finite-field theorem is not triggered by an arbitrary eigenslice or a
bare Neron--Severi embedding.  Its complex-K3 input must retain:

1. exact finite-volume return in dimension at least two;
2. an admissible full restriction-of-scalars trace;
3. the finite projective order module and its exact reflection action; and
4. the complex-K3 Picard embedding required by the declared branch.

The terminal documents now keep these gates even in their local summaries.
The audit has explicit regressions for the previously dangerous shorthand
“exact finite-volume full-trace K3 carrier,” which omits admissibility and the
order-module data.

Two different statements must also stay separate:

- The rank-six projected \(\Delta(2,4,6)\) construction proves **existence for
  one packet**.
- Every exact projected complex-K3 return, whenever it exists, obeys the K3
  rank ceiling \(n\le19\).

The first statement is not a class-wide existence theorem; the second is a
carrier-relative rank implication.

Finally, ambient projected Picard-host finiteness is false in general.  The
golden scaling/glue family changes ambient discriminant and glue while fixing
the native similarity class and normalized chamber.  This does not obstruct
the abstract maximal-chamber count: after finite native pairs **and reflection
maximality** are known, only the native similarity host matters.  Ambient
projected-host data remain relevant to exact geometric return; they are
irrelevant only to this downstream abstract count.

## 4. The `X3` quantifier

The logical content is

\[
 \neg C(X_3)\quad\Longrightarrow\quad
 \neg\bigl(\forall P\in\mathcal F_{\mathrm{vol}},\ C(P)\bigr),
\]

where \(C(P)\) is exact-module projected-root complex-K3 return.  It is **not**

\[
                 \forall P\in\mathcal F_{\mathrm{vol}},\ \neg C(P).
\]

At rank twenty the projector complement vanishes.  The primitive norm-(-2)
and divisibility-three norm-(-6) facet lines, under one rational common
similarity, would require

\[
                     ca^2=1,\qquad3cb^2=1,
                     \qquad (a/b)^2=3,
\]

impossible over \(\mathbf Q\).  That certifies one counterexample to the
universal all-finite-volume proposal.  Absolute reflection maximality of
`X3` is not proved, so the example does not settle a maximal-only proposal.
Its weaker mixed-divisor complex-K3 realization remains positive.

The independent replay verifies that \(3\) is not a rational square and
scans every synthesis document for both the positive one-counterexample
wording and the rejected pointwise-overclaim patterns.

## 5. The `H21` endpoint

The Borcherds chamber has a different and stronger maximality status.  The
frozen hostile audit verifies that the 42 short roots integrally span
\(L=U\oplus D_{20}(-1)\), have connected odd graph, and induce
\(GL(2,\mathbf F_2)\) on \(A_L\cong(\mathbf Z/2)^2\).  Any lattice stabilized
by a discrete reflection overgroup lies between \(L\) and \(L^\vee\); the only
stable choices are the endpoints, and they have the same full reflection
group.  Thus the side group is absolutely reflection-maximal.

Its current complex return stack is empty for three independent reasons:

- rank \(22>20\) excludes full-trace and projected-root complex K3 return;
- the direct rational branch is empty; and
- exhaustive rational physical square classes and adjunction reduce
  `RLogComp` to a \(PG(2,4)\) point-line incidence that cannot occur over
  \(\mathbf C\).

For the embedded Fano subplane, the replay obtains determinant \(-2\), which
is nonzero in characteristic zero but zero modulo two.  This is consistent
with, and does not erase, the positive characteristic-two supersingular K3
carrier: its 42 short roots are curve walls and its 168 long reflections are
Cremona actions.  Accordingly `H21` disproves exhaustiveness of the **current
complex** stack, not existence of every possible AG carrier.

## 6. Six-document interface audit

Each document is checked on its own, rather than only through one aggregate
concatenation.

| document | local links checked | terminal interface |
|---|---:|---|
| `GOAL.md` | 502 | PASS; Nikulin dependency exception, finite-pair maximality, repaired `X3`, `H21`, and exact field gates |
| `RESEARCH_PLAN.md` | 347 | PASS; one-packet/category-rank split and all terminal scopes |
| canonical completion | 56 | PASS; classification import, maximality, current-complex obstruction, and all six evidence links |
| inverse ledger | 300 | PASS; projected-host finiteness is false/unneeded, not an open gate |
| Nikulin ledger | 69 | PASS; exact logical arrows and hostile-scope table |
| synthesis theorem | 15 | PASS; conditional closure and outer-routing obstruction |

All 1,289 local Markdown links resolve inside the repository.  The known
stripped-TeX token scan has zero hits in every document.  The completion
report does not repeat the characteristic-two construction verbatim, but it
links both frozen `H21` owner/audit reports and makes only the narrower current
complex-stack claim; the other five synthesis documents state the positive
escape explicitly.

## 7. Exact terminal scope

| assertion | status |
|---|---|
| finite native pairs + reflection maximality imply finite normalized maximal chambers | **proved**, importing Nikulin 1982 |
| same conclusion from AG alone | **not proved** |
| global ABSW from the current carrier stack | **not proved** |
| projected ambient-host finiteness | **false in general and unnecessary for the conditional abstract count** |
| exact projected \(\Delta(2,4,6)\) return | **proved for one packet** |
| \(n\le19\) for every exact projected complex-K3 return | **proved carrier-relative rank implication** |
| universal exact-module projected-root return over all finite-volume inputs | **false by one `X3` counterexample** |
| maximal-only projected-root universality | **unsettled by `X3`** |
| `H21` absolute reflection maximality | **proved** |
| current complex return stack on `H21` | **empty** |
| characteristic-two mixed curve/Cremona return on `H21` | **positive** |
| class-wide mixed divisor/action compiler | **not proved** |

## 8. Independent replay

- [audit script](../scripts/audit_ag_to_absw_fixed_data_outer_routing_synthesis.py)
- [stable archive](../data/enumerations/ag_to_absw_fixed_data_outer_routing_synthesis_independent_hostile_audit.json)
- [focused tests](../tests/test_ag_to_absw_fixed_data_outer_routing_synthesis_independent_hostile_audit.py)

Run

```text
python3 -m py_compile scripts/audit_ag_to_absw_fixed_data_outer_routing_synthesis.py
python3 -m scripts.audit_ag_to_absw_fixed_data_outer_routing_synthesis --check-output
python3 -m pytest -q tests/test_ag_to_absw_fixed_data_outer_routing_synthesis_independent_hostile_audit.py
```

The archive records live input hashes and the audit payload digest, but does
not contain its own SHA-256 pin.

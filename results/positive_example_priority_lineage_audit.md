# Positive-example priority audit: what is new, and what was already known

Status date: 2026-08-07.

## Verdict

The positive fishing results are genuine finite-volume noncompact hyperbolic
Coxeter polytopes, but the word **new** has to be attached to the correct
layer.

- None of the present `H10`--`H14` examples is first in its dimension.
  Allcock proved that there are infinitely many finite-covolume Coxeter
  groups in every dimension through 19.  The `H19` and `H21` ledger entries
  are themselves classical reconstructions, not new catches.
- No row currently proves a new rational arithmetic commensurability class.
  Most literal chambers are finite-index reflection-subgroup chambers inside
  an explicitly identified classical full group.  Three rows still have an
  open *integral maximal-reflection* lineage, but their rational classes are
  already known.
- Some literal weighted chambers may be unrecorded.  That is the live
  priority question.  For every such row the present evidence is bounded;
  there is no worldwide-priority theorem.
- Several algebraic/divisor realizations are genuinely new extractions from
  the maintained source stocks even where the abstract chamber is old.

The sharp current wording is therefore: **new construction in this project,
often a new algebraic realization, sometimes a bounded literal-diagram
candidate; never a first-in-dimension claim, and presently never a proved new
arithmetic class.**

[Allcock's primary theorem](https://arxiv.org/abs/0903.0138) gives infinitely
many finite-covolume examples in every `H^n`, `n <= 19`.  The
[Felikson--Tumarkin finite-index theorem](https://arxiv.org/abs/math/0305093)
and its [Coxeter-group extension](https://arxiv.org/abs/0705.0426) constrain
the number of facets of a finite-index reflection subgroup, but do not
classify the possible subgroup chambers.  They therefore do not close the
literal priority question for the large chambers in this ledger.

## Exact negative-control theorem for the pinned CoxIter lanes

The two pinned CoxIter lanes are **not** a denominator for published
hyperbolic Coxeter polytopes.

The exact scan covers all 1,214 ordinary records and all 2,959 records in the
index-two doubling archive at commit
`b800be48d8240fbabd2d1c9542e61054584ad7fe`.  The 56 positive ledger rows use
43 distinct `(dimension, facet-count)` headers.  Of those 43 headers, 30 are
absent from both lanes; equivalently, 35 of the 56 rows have no matching
header.

Most decisively, the following six **prior-known controls** also have zero
matching headers in both lanes:

| control | header | independent prior lineage |
|---|---:|---|
| `R11-D4A1x5-NEF-H10-90` | `H10/90` | explicitly published in Roulleau's K3 atlas |
| `K3-D4-H12` | `H12/17` | full determinant-six chamber, Scharlau--Kirschmer `List12` entry 6 |
| `SS2-2C-H14-98` | `H14/98` | exact doubled-dual copy of Alexeev `(15,9,1)` |
| `D16-D5-H14` | `H14/60` | exact entrywise copy of Alexeev `(15,7,1)` |
| `X3-H19` | `H19/36` | Vinberg's most algebraic K3 chamber |
| `BOR-D4-H21-210` | `H21/210` | Borcherds' classical Lorentzian-lattice chamber |

The relevant primary or maintained sources are
[Roulleau's atlas](https://arxiv.org/abs/2003.08985),
[the Scharlau--Kirschmer census](https://www.math.rwth-aachen.de/~Markus.Kirschmer/lorentz/index.html),
[Alexeev's 2-elementary classification](https://arxiv.org/abs/2209.09110),
[Vinberg's K3 paper](https://doi.org/10.1007/BF01456933), and
[Borcherds' Lorentzian-lattice paper](https://doi.org/10.1016/0021-8693(87)90245-6).
Each control's exact local identification is linked in the JSON certificate.

This proves the scope limitation without relying on a vague warning:

> Zero CoxIter headers do not imply that a chamber is new.  Even an exact
> normalized-Gram nonmatch against every same-header CoxIter record proves
> only nonoccurrence in those frozen lanes.

The Scharlau--Kirschmer archive has a different denominator: it lists
reflective Lorentzian lattice genera and face/cusp counts for their **full**
reflection chambers.  Its own documentation states a conditional
completeness qualification in dimensions above five.  In any event, it is
not a census of all finite-index reflection-subgroup chambers.  Thus neither
CoxIter nor Scharlau--Kirschmer can by itself settle worldwide priority for a
37-, 50-, 123-, or 196-wall subgroup chamber.

## Evidence grades used from now on

The ledger's lineage phrases have the following exact meanings.

| grade | evidence | conclusion allowed |
|---|---|---|
| `K` | explicit facet permutation and entrywise normalized-Gram equality, or an exhaustive all-root proof identifying the full chamber | same abstract chamber |
| `L` | integral lattice isometry or exhaustive reflection-stable overclosure, plus containment in an identified full reflection group | lattice/full-group lineage closed; literal subgroup chamber may still be new |
| `B2` | every same-header candidate in a pinned corpus is an exact relation/Gram nonmatch | literal diagram absent from that bounded corpus |
| `B1` | no pinned record has the candidate `(dimension, facets)` header | header absent from that bounded corpus only |
| `S` | bibliographic or web search found no match | discovery aid only; no mathematical absence claim |

`B1`, `B2`, and `S` never imply worldwide priority.  Facet-count equality,
face-vector equality, determinant equality, or commensurability also never
imply chamber equality; grade `K` requires an explicit exact identification.

## The `X2-RAT-C2-H10-37` lineage, layer by layer

The rational-involution chamber is the strongest current algebro-geometric
example, but not a new dimension or arithmetic class.

| layer | exact status |
|---|---|
| dimension | `H10` was already populated, indeed infinitely, before this project |
| fixed divisor lattice | the known even 2-elementary lattice `(11,9,1)=U+E8(-2)+<-2>` |
| maximal reflection group | the known 13-wall `W_r(11,9,1)` chamber of Alexeev--Engel / Scharlau--Kirschmer |
| subgroup chamber | a 37-wall finite-index reflection subgroup, exact index `241,920` |
| pinned literal-diagram audit | all 45 37-node CoxIter records, in every declared dimension, are excluded by the complete relation-label histogram |
| worldwide literal priority | open |
| algebraic realization | newly extracted here from a genuine rational involution of Shimada's singular K3 surface `X_2` |

[Shimada's primary paper](https://arxiv.org/abs/1412.6904) supplies the
complete 478-wall chamber, the genuine group
`Aut(X_2,a_2)=PGL_2(F_9)`, and the rational involutions.  It does not state
the `X_2` 37-wall fixed section.  Its analogous quotient construction is the
40-wall Enriques chamber obtained from `X_0` in Section 10.  The arithmetic
denominator is supplied by
[Alexeev--Engel](https://arxiv.org/abs/2208.10383), not by a new lattice found
here.

The strengthened X2 comparison is much better than a missing `37 10`
header: it selects every 37-node record regardless of declared dimension and
excludes all 45 by an isomorphism invariant.  But the negative-control theorem
above shows why even this remains bounded.  An older diagram could be absent
from CoxIter or appear only in an unindexed paper.

## What remains genuinely open in the current positive ledger

Three candidates still lack an identified maximal **integral** reflection
chamber:

- `A15-C2-H10-34`;
- `JK-H12`; and
- `A20-C2a-H12-19`.

Their rational commensurability classes are already known, so closing these
three rows can identify a lattice and subgroup index but cannot, without new
evidence, create a new rational arithmetic class.

Two current abstract chambers have multiple exact divisor carriers inside
the project:

- the 27-wall `H10` chamber occurs as `D16-D4-C7-H10-27`,
  `A15_7-C2-H10-27`, and `GHESS-C2-H10-27`;
- the 34-wall `H10` chamber occurs as `A15-C2-H10-34` and
  `SS2-V4-H10-34`.

These are algebraic-realization multiplicities, not additional abstract
polytopes.  Conversely, the large finite-index chambers such as
`F3-C2-H10-196`, `F3-V4-H10-123`, `SS2-2A-H13-50`, `JK-H12`, and
`X2-RAT-C2-H10-37` remain legitimate bounded literal-priority candidates.
The ledger records the exact comparison denominator for each one.

## Reproducibility

- [negative-control certificate](../data/enumerations/positive_example_lineage_header_controls.json)
- [memory-light certifier](../scripts/certify_positive_example_lineage_header_controls.py)
- [regression test](../tests/test_positive_example_lineage_header_controls.py)
- [X2 strengthened bounded audit](shimada_aut3singk3_x2_h10_37_global_lineage.md)
- [complete positive ledger](higher_dimensional_polytope_lineage_ledger.md)

Replay with

```bash
python3 scripts/certify_positive_example_lineage_header_controls.py \
  --check-archive
PYTHONPATH=src:. pytest -q \
  tests/test_positive_example_lineage_header_controls.py
```

The replay is a corpus/header theorem and provenance check.  It is not, and
is not advertised as, a worldwide literature exhaustion.

Frozen SHA-256 values:

```text
ordinary CoxIter archive  ae04648b291c1e5f3d0e47d666c846dc805b244469a20486bca59674b4d171c3
doubling CoxIter archive  3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036
negative-control script   ab71e62c7e74147d9d0c36918cd0146c316c205a92404015d7ef6e8ede2b7678
negative-control JSON     58e00b546f50a9f3a3bcd4a21e16414aad6bb0d65f471475b6bfc199a57421aa
negative-control test     1840d9a60013648ada5e0a46c2d6086989914dfc7390c9191d352fc663f4ce07
```

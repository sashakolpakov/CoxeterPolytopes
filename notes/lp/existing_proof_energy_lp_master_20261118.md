# An LP master for recovering the existing compact energy proofs

This is a hand design, authored by `/root/pressure_review`. It aggregates
proved local inequalities and their complete applicability domains. It is
not a new geometric relaxation that attempts to reprove those inequalities.
No master table has yet been generated or solved; no new recovery or bound
is claimed. No scientific module was imported or executed for this note.

## Inputs and the one-vector rule

For an actual compact polytope choose one strict interior point `z` and one
original farthest vertex `p`. The accepted geometry writes

    z = E(p-x),    x in p-perp,    ||x||^2 = 1-E^(-2) < 1.

The incident spherical irreducible components give orthogonal subspaces.
For a component slot `j` put `e_j=||proj_j x||^2`. Thus `e_j>=0` and
`sum_j e_j<1`. Every energy in a master branch refers to this SAME `p,z,x`.
An adjacent-vertex vector has no additive energy account here unless a
proved local theorem has already transferred its conclusion to these slots.
Auxiliary points and original-wall checks belong to those retained local
theorems, not to a new assertion of feasibility by the master.

The primitive inputs to the master are:

* the accepted original-vector geometry and spherical classification;
* a complete necessary cover at each actual rank under consideration,
  including all actual defect and terminal alternatives;
* named local theorems, with their original ranks, component types, actual
  wall/mark hypotheses, spectator quantifiers and full parameter ranges;
* any explicit external ceiling used to make the remaining rank range finite.

These inputs must be retained genuinely, with their source statements and
proof authority. Their conclusions may be used directly. Reusing a proved
theorem does not require rerunning its thousands of semantic applications.
A new data-to-theorem applicability join still needs review; a source label
alone does not establish that join. Historical execution and conditional
inputs retain their recorded scopes.

In particular, neither `n<=17` nor a proposition excluding an entire rank is
a primitive row. The purpose is to recover those global conclusions from
local inputs. A theorem excluding a precisely specified local terminal
configuration may instead be an explicit contradiction leaf, as below.

## Branch data and admissible energy cuts

For each actual rank `n`, let `B_n` be a finite complete collection of branch
labels. Every hypothetical actual configuration at rank `n` must map to at
least one label. Labels retain the component ranks/types and every hypothesis
needed to apply their rows. Infinite native parameters are allowed only when
the local theorem covers their entire retained range. A finite sample of
`I_2(m)` labels is not coverage of the unbounded ordinary family.

A branch has a fixed finite set of energy slots and rational rows

    A_b e >= c_b.                                      (1)

The usual row is `sum_(j in S) e_j >= c`, obtained by safely weakening a
proved strict lower bound on exactly that support. Rational lower bounds for
irrational expressions require a proved inequality; decimal rounding alone
is insufficient. A theorem about a proper coordinate subspace of a component
gives a lower bound for that whole component, but two overlapping subspace
bounds cannot be added as though their subspaces were orthogonal.

More generally one may refine the slots into proved orthogonal subspaces,
or use another proved rational linear energy cut. Every coefficient must
refer to that fixed energy account. An unsquared pressure, a coefficient
price, or a bound on a differently optimized auxiliary vector is not an
independent energy cost. The local reference theorem must first pay all its
helpers, endpoints and heads and prove the asserted simultaneous bound.

Conditional estimates must keep their hypotheses. For example, a theorem
proved only under a component budget or one terminal pattern does not become
an unconditional floor. One can use a complete proved branch split, retaining
the applicable guard and its boundary treatment, or use an already proved
unconditional consequence. The master does not silently replace strict
guards by weak guards and then apply a theorem outside its stated domain.

The resulting branch polyhedron is an outer set of *energy consequences*,
not an assertion that its points are realizable diagrams or polytopes.
Feasibility of this master can therefore mean that more of the existing
local proof must be encoded. It does not refute the accepted proof.

## Shared defects and overlapping bounds

Use the accepted compatibility theorem: at most one actual `I/J` profile or
one typed `N_7/O_7` exception, with no coexistence. `J` occupies two old
`A_1` factors. `N_7` is the stated old-`A_2`/new-isolated-7 profile and `O_7`
the stated old-`I_2(7)`/new-isolated-3 profile; a bare `I_2(7)` is not `O_7`.
All occupied ranks count in the same rank total. An actual exceptional
component does not simultaneously receive the ordinary cost it replaces.

Joint bounds remain rows of (1). If both

    e_I + e_A >= a,    e_I + e_B >= b

are valid, retain both. Adding their right sides as a cost for
`I+A+B` charges `e_I` twice. The LP may use their legitimate combined effect
by choosing the shared `e_I`, or its dual may combine them with coefficients
whose total charge to every energy slot is at most one. Several alternative
floors on one component similarly give their maximum through simultaneous
rows, not their sum. A single disjoint partition into funded supports is a
sufficient certificate, but is not the only valid dual combination.

A generic spectator formula, such as a rank rate less one defect allowance,
is permitted only with a theorem covering every type and defect to which it
is applied. Rank alone does not justify dropping the actual defect identity.

## A normalized, ordinary rational feasibility LP

Fix a target threshold `N` and an explicit external integer ceiling `U`.
Keep all branch labels in

    B_[N,U] = union_(n=N,...,U) B_n.

For every nonterminal branch introduce a mass `tau_b>=0` and scaled energies
`E_bj>=0`. Put `delta_b=tau_b-sum_j E_bj` and impose

    A_b E_b >= c_b tau_b,
    delta_b >= 0,
    sum_b delta_b = 1.                                 (2)

This is a rational LP with no strict inequality. It is a normalized
convex-cone lift of the union of the branch energy polyhedra with total
energy strictly below one. It has no fractional ambiguity about a component
type: each positive branch mass has its own complete type/branch data.

To prove the equivalence, an energy point `e` in one branch with
`s=1-sum e>0` gives `tau_b=1/s`, `E_b=e/s`, all other masses/energies zero.
Conversely, (2) has at least one branch with `delta_b>0`. For that branch
`tau_b>0`, and `e=E_b/tau_b` satisfies (1), `e>=0` and `sum e<1`.
Thus (2) is feasible exactly when at least one listed branch has an energy
point satisfying these rational consequences below the unit budget.
This is an exact statement about the listed polyhedra, not about geometry.

An actual polytope supplies one covered branch and its original energies,
so it supplies a feasible point of (2). An exact rational infeasibility
certificate for (2) therefore excludes all actual ranks `N,...,U`, at the
retained local and cover premises. No whole-rank exclusion was inserted.

Merely replacing complete branches by fractional type counts, rank balance
and a few selector inequalities may enlarge the master. It remains sound
for exclusion if every actual configuration has a lift, but it need not
recover the proof. Any claimed exact compression must preserve the relevant
branch energy projection, rather than assume integrality of its selectors.

## The transparent capacity dual

A particularly simple exact witness gives, for every nonterminal branch,
a rational vector `alpha_b>=0` such that

    A_b^T alpha_b <= 1  (coordinatewise),
    c_b^T alpha_b >= 1.                                (3)

Indeed, using `E_b>=0`,

    sum_j E_bj >= alpha_b^T A_b E_b
                 >= (alpha_b^T c_b) tau_b >= tau_b.

Hence `delta_b<=0` for every branch, contradicting (2). The coordinatewise
capacity condition is the explicit no-double-counting rule. It applies to
each original energy slot, including one shared defect slot. Inconsistent
linear branch rows may also be closed by their ordinary exact Farkas
certificate; no numerical solver status substitutes for a checked witness.

For these rational branch inequalities, standard finite-dimensional LP
duality/Farkas alternatives supply a complete way to check the final linear
step. The proof of validity of the geometric rows is the retained local
proof, not a consequence of the numerical master solve.

## Non-energy terminal leaves and proof scope

Some accepted local arguments end with a feasible null vector excluded by
compactness, a forbidden actual native contact, or a final inner-product
contradiction. They need not assert any component floor. A complete branch
matching such a theorem can carry the explicit row `tau_b=0`; its energy
capacity then gives `E_b=0`. Record the local theorem and its full hypotheses
as the authority for this row. Do not manufacture an energy cost from it.

This yields a legitimate LP aggregation of a retained theorem proof tree.
Its exact description must distinguish numerical energy witnesses from
these local contradiction leaves. It does not claim that linear programming
has independently reproved the leaves' geometry or field theory. Importing
the final `n<=17`, `n<=15` or an entire exact-rank exclusion as such a leaf
would defeat the stated proof-recovery task. Local whole-component terminal
theorems may be retained only under their full local domain and with their
role disclosed, rather than hidden behind a rank-level citation.

## Threshold monotonicity and the separate external ceiling

Use one frozen collection of actual-rank branches, rows and theorem joins.
Then increasing `N` only removes branches from (2): feasibility at `N+1`
implies feasibility at `N`, and exclusion at `N` implies exclusion at every
larger threshold through `U`. The two statements follow by zero extension
or deletion of branch masses; no geometric rank deletion is involved.
Changing cost recipes or recompiling a weaker basis at a new threshold
would require a separate containment argument.

For the general compact proof the permitted external ceiling is the named
published `n<=29`. Recovering exclusion in each actual rank `18,...,29`,
or equivalently one complete master (2) on their union, gives `n<=17`.
The external theorem handles `n>=30`. Exact-rank-18 exclusion alone does
not exclude higher dimensions: neither faces nor deleted diagrams are
assumed to inherit the problem's geometric properties.

For compact arithmetic recovery, the already recovered *independent*
general compact theorem can supply `U=17`. Complete arithmetic branches at
actual ranks 16 and 17, with their own local field/native and computed-input
premises, would then recover `n<=15`. Until that compact dependency has
actually been recovered by the master, arithmetic recovery remains
conditional on it; no reciprocal use of arithmetic15 to prove compact17
is permitted. A fully independent arithmetic run can instead retain the
external29 ceiling and cover every required intervening rank. Existing
arithmetic computed-output and golden-lattice premises are not erased by
putting their consequences into an LP.

## Uniformity and source anchors

The assembler is independent of the case: consume complete branch data,
instantiate slots and valid rows, and verify (2)/(3) or an ordinary Farkas
witness. Within each local family use its uniform theorem/recipe over its
entire stated domain. Uniformity does not require replacing different
proved condition families by a single formula or discarding their guards.
This design authorizes no new local enumeration or numerical run.

Source anchors inspected for the general compact design are the accepted
working article's `eq:original`, `lem:pointpressure`, `eq:pointenergy`,
ordinary/defect table and `eq:joint`, `prop:prefix`, `sec:forests`,
`thm:tail`, `sec:verification`, and
`app:tail:prop:exact19cover`. The latter explicitly demonstrates distinct
actual-rank costs and defect accounting. The arithmetic dependency scope
is read from its accepted proof map, not newly reverified here.

| Retained source | SHA256 |
|---|---|
| `manuscripts/revisions/compact_farthest_vertex_bounds_consolidated.tex` (936,002 bytes) | `6222f5c6cee8b30b2c1019d648bc31cde8bfcddf08f031a0355a3e7137c107ac` |
| `docs/compact17_consolidated_proof_map_20261025.md` | `85004df24d7199fecf365d8fb39bc058342da9284142c52ff7113988699ad81e` |
| `manuscripts/revisions/compact_arithmetic_dimension_bound_consolidated.tex` (opaque binding only) | `d0f040b253146f318bafd2391c8809c9d995f53fcd920477f78f42038d424295` |
| `docs/compact_arithmetic_consolidation_proof_map_20261025.md` | `4a5defd68509bf9ada54434be0ac03fc01100c2ce5b11302da2cf9b934db1765` |

This note is a design and hand proof. Actual complete tables, applicability
joins and an exact checked master witness remain to be supplied and reviewed.

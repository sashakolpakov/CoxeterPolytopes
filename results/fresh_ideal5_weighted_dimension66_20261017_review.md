# Independent reconstruction of the dimension66 comparison

The [dimension66 primary proof](fresh_ideal5_weighted_dimension66_20261017.md)
passes independent reconstruction of its complete ideal arithmetic and
the newly changed finite rank17 interfaces. The reviewer imports only the
Python standard library; it imports no primary or project arithmetic.

[Independent source](../scripts/fresh_ideal5_weighted_dimension66_20261017_review.py)
and [saved certificate](fresh_ideal5_weighted_dimension66_20261017_review.json).

## Independent complete ideal domains

The primary enumerates sorted partitions by adjoining one part at a time.
This reviewer instead chooses the optional low part0 through9, the
optional rank12 part, and the optional high part at least18 first. It then
enumerates multiplicities of the repeated ranks10,11,13,14,15,16,17 to fill
the remaining rank. The categories are disjoint, and duplicate profiles
are explicitly forbidden. Thus this is a separate enumeration of every
terminal profile in the proved concentration domain.

For each profile, the reviewer multiplies the uniform polynomials using
a sliding window. The primary uses direct nested integer sums. Every
coefficient is compared, along with the total link vertex count V and
symmetry. The curvature is reconstructed in three mutually checked forms:

    G = sum_j h_j min(j,d-j),
    G = floor(d/2)*V - 2*sum_t a_t binom(floor(d/2)-t+s,s+1),
    kappa/V = 3d/4 - (5/2)*E|sum U_i-d/2|.

Here s is the number of factors and a_t are the coefficients of
product_i(1-z^(r_i+1)). The second identity sums the cumulative bounded
composition formula by the hockey-stick binomial identity. Binomial
terms outside their natural domain are zero. All arithmetic uses integers
and exact fractions, with no floating-point comparisons.

Every primary profile polynomial, G, normalized curvature, factor cost
and margin agrees with the independent reconstruction:

|Ambient dimension|Profiles|Minimum ideal margin|Minimizing ranks|
|---:|---:|---|---|
|67|1,761|851/165888|(11,11,11,11,11,11)|
|68|1,906|411499/3234816|(11,11,11,11,11,12)|

All3,667 margins are strictly positive. The independent source also checks
the duplicate12 cost increment23/1092, the linear low-rank cost through10
and the constant high-rank cost from17. The all-rank convex-order proof
is inherited, as is the actual cusp-factor multiplicity. The changed
finite rule leaves both that ideal envelope and the complete terminal
domain unchanged.

## Finite rank17 interfaces reconstructed separately

The reviewer constructs the canonical finite A17, B17 and D17 labelled
graphs directly. It independently counts their short8 and actual-leaf9
basic pairs and reconstructs the prescribed retained rules of costs
114,105,127. It computes all diagram automorphisms and checks invariance
of the removed sets.

For all140 eligible contexts from the complete2,614-context overlap
certificate, it recovers native labels from the original positively
gauged Gram matrices using4*g_ij^2/(g_ii*g_jj). It then independently
enumerates induced labelled maps using a fixed breadth-first pattern
ordering. This differs from the producer's adaptive next-vertex choice.
Every permitted target image is tried, testing all earlier edges and
nonedges. Necessary degree pruning does not discard a valid embedding.

The complete independently enumerated map sets agree exactly:

|Target|Contexts with maps|Maps|Distinct marked images|
|---|---:|---:|---:|
|A17|15|48|10|
|B17|5|18|3|
|D17|75|168|20|

Every selected marked image is retained. The reviewer also recomputes
each debit from the saved whole-component marked images under the three
independently reconstructed rules, then recomputes all174 literal totals.
Every total remains at least1. This part checks exact payment arithmetic;
the literal core inertias and the whole-component sealing proofs are
inherited from the frozen local theorem and its separate full reviews.
The complete original geometric census is not rerun here.

## Finite budgets and strict composition

The exact finite-forest upper bounds are supplied by the separately
reviewed [finite-budget certificate](fresh_rank17_reduced_forest_frontier_20261017_review.md).
This reviewer checks that its full variant uses A17=114, D17=127 and
B17<=127, satisfied by the actual B17 cost105, and that all76 saved maxima
agree between the finite producer and its independent verifier.

The attaining forests at the two excluded dimensions are checked directly:
construct finite D67 and D68, count every basic paid pair, and obtain
504 and512, respectively. Neither witness contains a rank17 component
or a cross-component coupon. The Bellman upper-bound and type-domination
certificate is not duplicated by this global reviewer.

The independently reconstructed global coefficients are:

|Ambient dimension|Finite maximum|Finite curvature|Finite slack|Global constant|
|---:|---:|---:|---:|---:|
|67|504|1023/2|15/2|330|
|68|512|527|15|335|

Together with the positive ideal margins, these give the two strict
contradictions in the primary proof. The completed ceiling68 supplies
the tail in every dimension at least69. This composition covers arbitrary
finite-volume Coxeter polytopes, including nonsimple cusps and either
absence of finite vertices or absence of ideal vertices.

The next comparison obstruction is also independently checked. At ambient66,
the formal profile(13,13,13,13,13) has margin-1565/10976. Direct counting
in D66 gives496, against finite curvature495, a deficit1. The entire
primary boundary record agrees. Neither obstacle is claimed to be a
realized hyperbolic polytope or an optimality statement.

## Review and verification scope

The reviewer read the complete primary source and mathematical report.
The simultaneous local deletions, unchanged coupons, separation of A17
from D17 in the finite budget, complete concentration, actual cusp
multiplicity, two strict global contradictions and completed higher tail
all pass. The arithmetic comparison has terminated exit0/PASS on every
profile and all stated interfaces.

Generationbd9375 and fresh separate-process saved replay58927 both
terminated exit0/PASS after binding the final primary source, certificate
and report. All25 direct bindings, Python AST and four local report links
pass.

Root read the entire independent source and its final interface changes.
The independent optional-first domain, sliding-window/CDF/absolute-deviation
identities, fixed breadth-first induced-map enumeration, exact retained
rules,174 saved literal debits/totals, direct D67/D68 counts, witness fields
and full ambient66 boundary comparison passed that mathematical/source
review. Root did not independently rerun this reviewer's arithmetic.
Source, certificate and report are final/frozen. No external refereeing,
publication-priority claim, proof-assistant formalization or new
construction is asserted.

    python3 scripts/fresh_ideal5_weighted_dimension66_20261017_review.py
    python3 scripts/fresh_ideal5_weighted_dimension66_20261017_review.py --verify

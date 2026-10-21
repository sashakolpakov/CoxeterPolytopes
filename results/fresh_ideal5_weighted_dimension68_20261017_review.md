# Independent scalar and composition review of the dimension68 proof

The independent reconstruction agrees with all2,060 exact ideal-profile
rows in the [dimension68 proof](fresh_ideal5_weighted_dimension68_20261017.md).
The complete mathematical composition read finds no remaining logical
obligation in that proof. It excludes dimension69, and the completed
dimension69 theorem excludes every higher dimension. There is no
arithmetic, compactness, simple-cusp or cusp-factor restriction.

The [reviewer source](../scripts/fresh_ideal5_weighted_dimension68_20261017_review.py)
imports only the Python standard library. Its
[certificate](fresh_ideal5_weighted_dimension68_20261017_review.json)
records the independent scalar data, pair weights, finite domination
checks and exact input bindings. Final generation/replay status is given
below. This is internal verification, not external refereeing or a
proof-assistant formalization.

## Complete profile reconstruction

The reviewer enumerates optional exceptional parts first, then recurses
over multiplicities of10,11,13,14,15,16,17. This differs from the primary's
nondecreasing-part recursion. The exceptions are at most one low part
between1 and9, one12, and one part at least18. At total rank68 the two
independent domains agree exactly on2,060 profiles, with no duplicates.

For each profile the reviewer multiplies the uniform-factor polynomials
by a sliding window, instead of the primary's nested coefficient loops.
It checks symmetry and total mass V. It then calculates G in two ways:

    G = sum_j h_j min(j,d-j),
    G = floor(d/2)*V - 2*sum_shift a_shift
                           * binom(floor(d/2)-shift+s,s+1).

The second formula is the cumulative bounded-composition identity, with
the coefficients a_shift obtained from the numerator product
`product(1-z^(r+1))`. It also checks the independent absolute-deviation
formula

    kappa/V = 3d/4 - (5/2)*E(abs(sum_i U_i - d/2)).

Every profile's h-polynomial, V, G, curvature, budget and margin agrees
with the primary. All margins are positive; the minimum is3641/13608
at `(11,11,11,11,11,13)`.

The changed concentration argument is essential. Ranks at most10 have
linear cost and concentrate only up to10; rank11 is retained. Reusing
the old cap11 argument would be invalid. The exact gain for
`12+12 -> 11+13` is23/1092, positive. High-rank concentration preserves
the constant cost. Every move increases the bounded integer squared-rank
sum, so iteration terminates. Curvature decreases and cost is preserved
or increased; positivity in the complete terminal domain therefore
implies positivity for every cusp profile of the given total rank.

The unchanged large-affine envelope is checked again from its15 exact
piecewise-linear spine branches and endpoint ratios, yielding163/23.
The all-rank graph-count proof and published affine classification remain
inherited inputs.

## One simultaneous local rule

The reviewer independently constructs the four canonical12-node affine
graphs and their graph distances. It reconstructs all264 pair weights:

|Family|Cost|Changed pairs|
|---|---:|---|
|A|60|Six antipodal distance6 pairs removed|
|B|65|Two terminal/fork-leaf distance10 pairs receive one half|
|C|65|The distance11 endpoint pair removed|
|D|62|Four opposite-fork-leaf distance9 pairs removed|

Every elementary ideal triangle pair still has unit weight. The cycle
has no feature, so its antipodal pairs were not elementary triangle
payments. Exceptional whole-factor triangle alternatives have order at
most9. The sole finite deletion has distance16, outside finite short8.

The complete source and proof reads cover the following interfaces.
The universal cycle12 transfer returns finite units/coupons or ideal
paths of length at most4. The universal D11 transfer returns finite
payments or ideal paths of length at most5. Neither can recurse into
another changed ideal pair. Their finite outputs avoid the deleted
finite A17 endpoint.

For the B rule, the14 affected cases split into eight radical sealed
cores, two already-whole finite B11 components, and four nonsingular
cores admitting zero or one forced supporting root. The reviewer
recounts all16 saved sealed payment variants: an ideal half plus at
least a finite half gives one. It independently recomputes the three
immediate-D13 strict-gap arithmetic identities,1/9 and1/11. The core
projection/contact classification, actual radical argument and whole
component sealing are mathematically reviewed inherited local proofs;
their28 maps and exact core inertias are not re-enumerated here.

The C theorem checks all758 eligible high patterns,66 low closures with
three optional variants, immediate whole B/D factors and all downstream
outputs. Its fifteen successful C contexts/30 maps avoid the endpoint.
Its proof also checks that the finite A17 deletion's only debited legacy
case662 still has two finite D17 units. The complete inherited high/low
classifications and closure proofs remain explicit inputs. This review
does not describe those earlier computations as a new independent census.

The short8-first reduction, universal transfers, complete no-short8 core
split, protected payments and literal later half-weights thus form one
noncircular rule. In particular, an old whole-small-factor phrase is not
silently reused after deleting its pair. The local reports inspect the
actual alternate factors, short paths or finite coupon outcomes.

## The finite budget includes the A11 coupons

The reviewer reconstructs the new coupon interaction totals directly,
without importing the finite-budget producer. Each A11 has internal
cost55, equal to D11. Its new total interaction is1 with A1 or A2 and0
otherwise. Upgraded D11 has respective totals5/2 and3 and dominates
every other partner. All228 finite-type partner rows agree with the
saved domination certificate, including the symbolic dihedral family.

Interactions between two large components vanish. Replacing A11
components successively by D11 therefore preserves rank and never
decreases total charge. The frozen exhaustive Bellman bound remains an
upper bound for the enlarged rule. Its separate no-producer review is
bound here; its full Bellman enumeration is not repeated.

The reviewer independently recounts all76 typed attaining forests with
the new coupons. Their costs agree exactly; none contains A11. Hence
they still attain the inherited upper bounds. In rank69 both
`A1+4D17` and `A1+4A17` have cost544. The A17/D17 distinction is retained:
their internal costs are135 and136, respectively, and their coupon
eligibility differs.

At n69 the finite curvature coefficient is544. All ideal margins are
positive and the global shape inequality has the extra constant340.
Subtracting the payment upper bound from that shape lower bound is
strictly positive even with zero cusps. This excludes69; the frozen
ceiling69 supplies n>=70, proving68.

At ambient68 the checked profile `11^5+12` has positive ideal margin
411499/3234816, but the finite maximum remains544 while its curvature
coefficient is527. The finite deficit17 prevents inferring ceiling67
from this comparison. The reviewer makes no such further claim.

## Verification status and limits

Independent scalar preparation and the full primary-row/interface
comparison both terminated exit0/PASS. After the primary report froze,
final hash-bound generation and a separate fresh saved replay also
terminated exit0/PASS. All45 direct bindings, report links and Python
syntax passed. Source, certificate and report are final/frozen.
The new arithmetic consists of2,060 profile reconstructions,264 pair
weights,228 partner comparisons,76 attaining forests and the explicitly
stated saved local-payment sums. Full mathematical/source review covers
the local composition, concentration and strict global contradiction.

Root separately read the complete reviewer source and mathematical
algorithms, passing the optional-first domain, sliding convolution,
cumulative and absolute-deviation identities, canonical BFS/pair weights,
independent partner comparisons and typed witnesses, and the explicit
inherited-upper-bound scope. That source review did not separately rerun
the2,060-row arithmetic.

Complete all-label core/overlap classifications, arbitrary spectator
geometry, physical supporting-root closures and the exhaustive Bellman
upper-bound certificates retain their frozen validation scopes. They
are not replaced by finite samples or claimed as newly reimplemented
here. The result is an internally checked computer-assisted proof,
with no claim of outside refereeing or Lean-style formalization.

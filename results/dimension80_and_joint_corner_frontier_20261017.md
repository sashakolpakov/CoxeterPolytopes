# Dimension80 and joint actual-corner constraints

**Completed internally checked computer-assisted proof:** every finite-sided,
finite-volume real hyperbolic Coxeter polytope has dimension at most80.
There is no arithmetic, simple-cusp or facet-count restriction. This improves
the [completed86 proof](order16_and_large_core_frontier_20261017.md).
The compact arithmetic21, arbitrary compact29 and compact-facet ambient79
ceilings are unchanged. No new finite-volume construction is claimed.

The [dimension proof](fresh_ideal6_piecewise_dimension80_20261017.md) and
[independent review](fresh_ideal6_piecewise_dimension80_20261017_review.md)
state the full dependency and verification scopes. This is a proof with
no currently identified gap, not merely numerical evidence for a conjecture.
All reviews are internal. External refereeing and proof-assistant formalization
remain outstanding; no publication-priority or current-best literature claim
is made.

## What completes the unrestricted bound

The [preceding checkpoint](ideal6_local_and_large_core_model_frontier_20261017.md)
had proved the stronger local triangle and quadrilateral lemmas but left
63 preferred contexts affected by nine long affine-D14 pairs unresolved.
The new [universal D14 transfer](fresh_affine_d14_long_pair_transfer_20261017.md)
handles all nine pairs directly, including arbitrary supporting roots.
Complete high-rank D extension profiles force one exceptional equality case:
two replacements of height1/2 and a single isolated A1 spectator contact.
The exact72 configurations give finite unit payments or two sealed finite
half-payments. Their whole-component condition is proved, not inferred
from a small core alone.

The [finite budget](fresh_broad_finite_half_coupons_20261017.md) includes these
payments and all inherited coupons, with total at most9n for n<=108.
The [complete local composition](fresh_ideal6_combined_payments_20261017.md)
therefore pays every triangle at least2 and every quadrilateral at least1.
Its ideal-factor budget per link vertex is

    b(r)=r/2 for1<=r<=13;
    b(14)=20/3; b(15)=111/16; b(16)=127/17;
    b(r)=95/13 forr>=17.

The whole-factor embedding census and all-rank spine counts justify the
exceptional and unbounded ranges. The actual cusp-factor multiplicity is
V/(r+1), so nonsimple ideal vertices are included.

Concentration reduces every partition of n-1 to copies of13,16,17, at most
one rank1,...,12, at most one14, at most one15, and at most one rank>=18.
The repeated middle-rank moves terminate because each increases the bounded
integer sum of squared ranks. Exact curvature comparisons are positive on
all5,008 terminal profiles in ambient dimensions81,...,86; the completed86
theorem supplies the entire higher-dimensional tail.

| Ambient dimension | Profiles | Least normalized ideal margin | Finite margin |
|---:|---:|---:|---:|
|81|742|940675/11294304|31|
|82|776|1132393/9143008|41|
|83|813|3418145/9143008|105/2|
|84|852|15382621/27429024|63|
|85|892|221295/326536|75|
|86|933|12970945/15059072|86|

Primary generation and separate saved replay pass for both the local binder
and the dimension comparison. A second agent independently reconstructs
every scalar row using different enumeration and convolution formulas, and
independently derives the12 affine cost branches. Separate full-source and
written mathematical reads pass for the local D14 argument and composition.
These do not duplicate all inherited classifications; those remain named,
hash-bound proof inputs with their own verification records.

At ambient80 the formal profile(1,13,13,13,13,13,13) has margin
-1131729/7529536. The current comparison stops there. This is a limitation
of the payment rule, not an actual chamber or a sharpness result. Further
progress requires cheaper order14 payments or stronger global allocation
and incidence information; a new global method is still needed to approach
the dimensions of known constructions.

## Compact arithmetic lane

The [joint raw87/88/89 package](fresh_raw87_88_89_actual_corners_20261017.md)
reconstructs89 actual-corner classes,712 corners and380 forests over
Q(sqrt5), with19,820 literal orientations. It distinguishes the two nontrivial
dotted squareclasses from the third family whose contact root lies in the
base field. Exact tests exclude2,530,570 of2,567,854 literal pairs; all37,284
unknown pairs remain allowed. Cross-family testing uses common forests and
full native contacts and future-height signs.

The1,002 capacity columns comprise440 clique covers and562 coloring covers.
At the fixed4,410 face masses,113 violated columns require75 vertex-allowance
increases, including21 additional vertices from cross-family constraints.
The repaired exact ratio is

    1342109395844304360/109100031098664527
    =12.301640818329087... >12.

This is a certificate for a partial finite model, not a new dimension theorem.
Every inherited inequality remains valid under the monotone allowance repair,
including signed coefficient columns whose variables are nonnegative. No
optimizer was invoked. The prior solved finite optimum
12.4145168684--12.4145233755 remains a historical result for its earlier model;
the enlarged model has not been solved. Additional actual-corner constraints
or complete implicit coverage are needed before another solve is useful.

Generation and fresh replay pass. Separate source/mathematical review checks
the native two-cycle field tests, actual-corner capacities and monotone repair;
independent rational checks reconstruct the1,002 capacities and fixed-witness
repair conditional on the saved pair-exclusion certificates. These checks
do not repeat the pair geometry or claim complete coverage of all ambient
configurations.

## Construction lane

The [completion gate](fresh_double_flip_completion_gate_20261017.md) excludes
both tested completions of the six saved double-flip prefixes. It covers450
complete label domains and900 norm branches, including unbounded high-label
families as intervals. Necessary finite-corner conditions reduce the possible
WZ and BX labels; all retained native contacts are checked. The certificate
records572 native exclusions and328 exclusions from all allowed finite WZ
labels. Generation, fresh replay and a separate narrow source/mathematical
read pass; an additional rational check validates the saved interval
inequalities without independently re-enclosing the physical Gram data.

This excludes the six fixed prefixes and two specified incidences. Deformed
prefixes remain open. The next construction direction varies R and Z jointly
and permits RZ<=5, with a derived admissible parameter domain. The exact
control still has63 walls,199 finite vertices and10 target-face open ends;
it has certified infinite volume. No new chamber has been promoted.

## Preservation and replay

The [incremental archive checker](../scripts/audit_dimension80_and_joint_corner_frontier_20261017.py)
and [archive certificate](dimension80_and_joint_corner_frontier_20261017_audit.json)
rehash every predecessor source and report before adding the seven new
packages. The archive records its independent rational checks separately
from the geometric proofs, the fresh package replays and the internal
source reviews. Frozen historical reports are preserved, including their
then-current bounds and explicitly conditional results.

    python3 scripts/audit_dimension80_and_joint_corner_frontier_20261017.py
    python3 scripts/audit_dimension80_and_joint_corner_frontier_20261017.py --verify

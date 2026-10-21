# Complete nonzero spacelike-face classification: only the original poles

For each original raw23 cone, consider proposed norm-two Coxeter normals
whose Lorentz projection is nonzero and belongs to the cone generated
by its32 spacelike extreme rays. The complete classification is now:

|Original pair|Possible normal in that domain|
|---|---|
|(1,1)|The original pure pole on ray31 only|
|(1,6)|None|
|(6,6)|The original pure pole on ray6 only|

This needs no extension-field, arithmeticity or finite-angle cutoff.
The surviving normals have zero passive part. Zero Lorentz projection
and directions involving timelike extreme rays remain open. No raw23
pair exclusion, arbitrary-completion theorem, shared capacity or
dimension-bound improvement follows from this restricted-domain result.

The [producer](../scripts/fresh_raw23_spacelike_pole_classification_20261018.py),
[certificate](fresh_raw23_spacelike_pole_classification_20261018.json)
and [saved checker](../scripts/check_raw23_spacelike_poles_saved_20261018.py)
complete the [contact-sum reduction](fresh_raw23_spacelike_contact_sum_20261018.md).
They preserve the [original pole packet](fresh_raw23_union_polar_completion_preflight_20261018.md).

That complete face-lattice theorem reduces every nonzero projection in
the stated cone to rays6 or31 in each case, excluding all1,929 face
interiors of dimensions2..8 and90 ray interiors. The present packet
uses individual contacts to settle the remaining six ray directions.

Let r be a ray of positive norm h, and write a projected normal as
x=t*r with t>0. The known scaled roots have norms2s_j. Put
p_j=−r·r_j>=0. Since a full physical norm-two normal has a positive-
definite passive part, t²h<=2. Its squared physical contact to wall j
therefore has the exact sharp upper

    t²*p_j²/s_j <= 2*p_j²/(h*s_j).

If p_j is nonzero and that upper is strictly below1, the contact
cannot be a Coxeter contact: every nonzero finite angle has square
at least1, and parallel or dotted contacts have square at least4.
This excludes ray6 of(1,1), both rays of(1,6), and ray31 of(6,6).
All strict comparisons are exact in Q(√5); they impose no field
restriction on the proposed normal.

For each of the other two rays, a nonzero contact to a norm-two
known wall has upper exactly1. Coxeter compatibility forces equality.
Thus t²h=2, and the passive part has squared norm zero in a positive-
definite space, so it is zero. The projection scale is uniquely
t=1/p_j, since s_j=1 and the contact itself is nonpositive.

The producer constructs that exact vector, verifies its norm2 and
all known squared contacts, and finds precisely two label3 contacts
and one label5 contact; all others are right. These are the same
original pure poles already present in the frozen packet. The new
result proves their uniqueness throughout the entire nonzero
spacelike-generated cone, rather than only among extreme-ray candidates
or within a fixed field's finite-angle catalogue.

The original after-pole cones still have exterior directions. This
classification does not remove those directions or account for
zero-projection normals. It supplies a larger exact local domain
for future completion arguments, not a complete arithmetic exclusion.

Generationef9a97, separate replay3d4027 and standalone saved checker3aef06
terminated exit0/PASS: four additional strict exclusions and two forced
pure poles. A peer completed a full mathematical producer/report/helper
read, also reading the complete contact-sum dependency and original
polar proof. The review checked the exhaustive nonzero face reduction,
physical contact scaling, strict gap, equality forcing projection norm2
and passive part0, exact pure-pole roles and explicit zero/timelike open
domain. It found no correction and executed no enumeration, exact
arithmetic, native classification or historical-catalogue replay.

The checker imports neither producer; it uses the earlier standalone
Q(√5) arithmetic to rebuild every contact bound, strict comparison,
forced scale and pure pole. The complete contact-sum theorem and
original ray completeness remain bound inputs. Source, certificate
and report are final/frozen after these checks and review. No external
refereeing or proof-assistant verification is asserted.

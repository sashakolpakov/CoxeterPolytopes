# Complete affine-circuit screen of the retained 49-wall H9 chamber

**PASS:** the retained 49-wall chamber has no irreducible affine circuit.
In particular it has no persistent affine cusp orthogonal to one of its
40 forbidden-angle pairs. This removes one possible obstruction to the
present search; it does not prove that a finite-volume completion exists.

The [exact pair inventory](fresh_h9_persistent_cusp_screen_20261016.json)
and [its independent geometry reconstruction](../scripts/fresh_h9_persistent_cusp_screen_20261016.py)
classify all 1,176 pairs as 367 right, 91 finite nonright and 718 dotted.
There are no parallel pairs. The finite census includes 40 edges of label
five or ten. Exact finite equalities are checked in the previously
reviewed nested field; strict dotted inequalities use rational intervals.
No numerical nearest-label decision is used.

The [affine enumerator](../scripts/fresh_h9_persistent_cusp_affine_20261016.py)
and [complete enumeration](fresh_h9_persistent_cusp_affine_20261016.json)
use that pair matrix. Both programs pass generation and fresh saved
replay. They leave the frozen 49-wall construction unchanged.

## Why the search is complete for this obstruction

An irreducible affine Coxeter diagram has only finite labels two, three,
four and six, except for the two-root parallel diagram. Every proper
principal subdiagram is positive definite. These standard facts follow
from the affine classification; see
[Felikson, section 1.2 and Table 1.1](https://www.maths.dur.ac.uk/users/anna.felikson/talks/cl_cox.pdf).

Start at every singleton. Add a root joined by a nonright allowed edge,
requiring every pair in the enlarged set to have an allowed finite or
parallel label. Compute its exact normalized Gram inertia. Continue only
from positive-definite sets, recording a first positive-semidefinite
singular set as an affine circuit. An irreducible affine diagram admits
a connected ordering of its vertices, and each proper prefix is positive
definite, so this procedure cannot miss one. Indefinite sets cannot have
a positive-semidefinite extension.

The matrix calculations use exact arithmetic in Q(sqrt(2),sqrt(3)),
including recursive sign decisions and exact symmetric elimination.
The whole-chamber search visits 630 distinct connected admissible sets;
all 630 are positive definite. The separately enumerated common-orthogonal
domains for all 40 forbidden-angle pairs contain no affine circuit either.

The ambient Lorentz signature bounds positive-definite rank by nine.
Orthogonality to a forbidden-angle pair removes a positive two-plane,
bounding it by seven in those domains. The search permits a final
singular extension of those sizes, so these conservative stopping bounds
do not discard a potential affine circuit.

## How this guides cap selection

A positive affine relation defines a nonzero future null vector nu.
Every additional acute-compatible wall r satisfies r.nu<=0, so it
cannot remove that ideal direction. If an incident pair has angle pi/5
or pi/10, finite volume would force a bounded Euclidean cusp link with
that pair. The product of its two reflections has irrational linear
trace; a cocompact Euclidean group preserves a full translation lattice
and every such trace must be integral. See
[Gelander, Theorem 2.9](https://www.math.utah.edu/pcmi12/lecture_notes/gelander.pdf)
and the repository's [earlier detailed argument](fresh_h9_fixed_core_cusp_obstruction_20261013.md).

Because the old chamber has no affine circuit at all, any newly introduced
one must involve a new wall. For a mutually dotted cap collection, an
affine circuit can contain at most one of its new walls. This makes the
screen local to each cap while still covering larger affine diagrams;
checking only affine triangles would miss possible obstructions.

Absence of this obstruction leaves all other finite-volume obligations
open. The retained chamber still has its independently certified strict
spacelike escape and infinite volume.

```sh
python3 scripts/fresh_h9_persistent_cusp_screen_20261016.py --verify
python3 scripts/fresh_h9_persistent_cusp_affine_20261016.py --verify
```

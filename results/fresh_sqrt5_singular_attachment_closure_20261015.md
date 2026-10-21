# Complete closure of the new singular attachment branches

**PASS for this specified frontier.** The five nonsealed singular starts
from the preceding 51-class stage, and all 16 longer single-A-chain stars,
are paid under the extended forest-sensitive rule for every further unpaid
white attachment. No rank cutoff is used in this closure. Nonsingular
multiple-component constructions and the symbolic dotted families remain
outside the result, so it establishes no global dimension bound.

The [checker](../scripts/fresh_sqrt5_singular_attachment_closure_20261015.py)
and [certificate](fresh_sqrt5_singular_attachment_closure_20261015.json)
bind the frozen [51-class stage](fresh_sqrt5_singular_extension_frontier_20261015.md),
the exact [single-chain completion census](fresh_sqrt5_single_chain_completion_20261015.json),
and the [extended forest rule](fresh_sqrt5_extended_forest_weights_20261015.json).
They preserve all preceding artifacts.

The seven accepted continuations of the five old starts and the 16 stars
with chain length greater than one give 23 literal roots. Full marked
canonicalization reduces them to 17 states, with a checked literal map for
every seed. All white subsets of eligible sizes are searched for one of
the 158 known sealed subcores. This gives 3,054 exact marked-subdiagram
tests; six states reduce to the known catalogue. The certificate checks
that the remaining whites are orthogonal to each displayed sealed subcore.

The other eleven states are newly sealed and already paid. Their complete
next-white domains contain 138 unpaid structural profiles:

| Exact outcome | Profiles |
|---|---:|
| Violates the physical core relation / conjugate kernel | 115 |
| Fails the finite-cycle ground-field condition | 22 |
| Creates a persistent entire-H4 payment | 1 |
| Survives | 0 |

The exact tests are inherited from the separately frozen 51-class checker:
complete finite-label profile enumeration, all four spherical corner
conditions, the short-path payment, affine witnesses, cyclic field,
physical relation and physical/conjugate inertia. This new primary run does
not claim an independent second implementation of those helpers.

The eleven new sealed states have indices 1,2,3,4,5,6,7,13,14,15,16 in this
certificate. Each has total weight exactly one both without and with an
external entire A1 component. Together with the 158 previous classes this
gives a catalogue of 169 sealed marked classes paid by the current rule.
No new pair weights are required at this step.

## Why this is an arbitrary-size closure

In an unpaid completion, an external white node with any nonright contact
to a sealed core would have one of its enumerated profiles. All such
profiles either fail a necessary condition or give a persistent short-path
or H4 payment. Thus every external white is orthogonal to the whole sealed
core. Longer paths cannot evade this first-contact argument. Its marked
corner components are entire; the remaining external forest changes the
rule only through its entire-A1 flag, and both possibilities are checked.

For each of the five initially nonsealed starts there are two cases. If no
external white touches it, its original core components are entire and the
checked weight is already one in either A1 context. Otherwise choose a
touching node. Unless it already gives baseline payment, its profile is one
of the seven enumerated continuations. Every such continuation belongs to
one of the 17 states now reduced to a paid sealed core or proved newly
sealed. This closes all subsequent attachments without assuming a maximum
number of white vertices. All 16 longer-chain stars are handled by the same
sealed-state argument.

The implementation permits an explicit unresolved queue if a future
frontier fails to terminate. In this run the queue actually exhausts after
17 states, before its safety cap of 100, and no new accepted child is
generated. Termination here is a checked finite reduction followed by the
all-size sealing lemma, not a finite-rank extrapolation.

Generation and saved replay pass, as do syntax, source manifests and local
links. The payment conclusion is conditional on the bound extended weight
rule and the inherited arithmetic/compactness hypotheses. Actual global
compact realization of every principal star is not asserted.

```sh
python3 scripts/fresh_sqrt5_singular_attachment_closure_20261015.py
python3 scripts/fresh_sqrt5_singular_attachment_closure_20261015.py --verify
```

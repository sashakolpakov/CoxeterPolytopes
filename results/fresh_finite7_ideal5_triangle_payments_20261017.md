# Finite distance seven suffices for the triangle payments

Every triangle of a finite-sided finite-volume real hyperbolic Coxeter
polytope receives at least two units under the following rule:

* At a finite corner, pay a marked graph distance at most7.
* At an ideal corner, pay distance at most5, or distance6 or7 when a
  mark is within5 steps of a branch vertex or an endpoint of a nonsimple
  bond in the actual varying affine factor.
* Retain the whole-factor ideal payment for factors of order at most9.

This sharpens the finite-distance requirement in the completed
[finite8 triangle theorem](fresh_ideal5_finite8_triangle_payments_20261017.md).
It does not yet allow removing finite distance8 units from the global
rule: quadrilateral coverage and the full changed finite budget remain
separate obligations. The current ambient ceiling stays65.

The [source](../scripts/fresh_finite7_ideal5_triangle_payments_20261017.py)
and [certificate](fresh_finite7_ideal5_triangle_payments_20261017.json)
re-enumerate all1,025 high-core marked triples and verify the inherited
four exceptional interfaces.

## Proof and arbitrary supporting diagrams

Choose an inherited minimal indefinite core containing all three black
side roots. Every other root of the core supports the triangle. A corner
path avoiding the third black persists after all supporting roots are
restored. Its displayed branch or nonsimple bond also persists. At an
ideal corner, the actual factor-fiber theorem places the two marks and
their path in the same varying affine factor.

For core order at most6, the unchanged all-label Steiner-tree argument
gives at least two distinct marked-pair paths avoiding the third mark,
each of length at most5. They satisfy the new rule.

For orders7 through10, the complete inherited core catalogue gives
105,224,336,360 marked triples, respectively. The source repeats every
triple and every actual corner path/feature test, now discarding every
elementary payment whose path has length8. No new exception appears.
Exactly the same four triples of the simply-laced10-node L25 core need
a non-elementary payment:

|Marked triple|Alternate marked pair|Path length|Protection|
|---|---|---:|---|
|(0,1,9)|(0,9)|6|Affine E8 positive kernel|
|(0,2,9)|(0,9)|6|Affine E8 positive kernel|
|(1,2,8)|(1,8)|6|Finite E8 inverse|
|(1,2,9)|(1,9)|7|Affine E8 positive kernel|

Each has one elementary payment, and its alternate pair is distinct
from that payment. The predecessor's exact kernel and inverse argument
protects a simple leaf and its degree-two neighbour in the entire
alternate marked component, allowing arbitrary additional supporting
roots. No classical irreducible affine Coxeter graph has this protected
simple-leaf configuration. Consequently the component is finite or
exceptional affine of order at most9. The displayed finite path has
length at most7; the exceptional affine case has its whole-factor unit.
Both cases retain the second payment.

The proof does not require that a protected leaf be a mark. It does not
infer whole-component finiteness from a truncated diagram. It uses the
same arbitrary-support E8 argument as the predecessor, with the sharper
explicit path lengths checked here.

## Verification scope

Generation61e45a and a separate fresh saved replay332ec3 both terminated
exit0/PASS. The source repeats all1,025 triples, checks exactly the four
old exceptions and reconstructs the E8 kernel/inverse interfaces. The
graph and matrix routines are shared with the predecessor; this is not
an independent arithmetic implementation or a new core classification.

A separate full source and inherited-proof read passed the low-core
all-label argument, new path threshold, four distinct protected
alternatives and arbitrary-support scope. That reviewer did not rerun
the numerical census. The complete new report also passed that review.
All11 direct source bindings, three local report links and Python syntax
pass. Source, certificate and report are final/frozen.
No quadrilateral theorem, changed finite maximum, new ambient bound,
new polytope, external refereeing or formal verification is asserted.

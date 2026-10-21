# Universal finite extensions: an exact conservative diagnostic

The [checker](../scripts/fresh_sqrt5_uniform_face_weights_20261015.py) and
[certificate](fresh_sqrt5_uniform_face_weights_20261015_rank29.json) give exact
optimum **154** for a deliberately conservative rank-29 face-weight system.
The target needed by the intended face-count argument is84. This result does
not exclude dimension29 polytopes or arithmetic-aware weighting methods.

The inputs are all8,298 marked classes of the
[complete quadrilateral union census](fresh_sqrt5_all_quad_unions_20261015.md)
and all76 minimal triangle markings. They give3,090 distinct face inequalities.
There are202 partial marked-pair types and9,946 entire-component marked-pair
types in all92 allowed irreducible finite Coxeter diagrams of rank at most29.
Every connected induced subtree on at most8 vertices is enumerated, yielding
38,038 extension inequalities. The auxiliary weight of each partial marked
pair is required to be no larger than its weight in every such extension.
Disconnected marked pairs receive zero.

An exact rank-partition dynamic program supplies the worst vertex forest.
Linear programming proposes coefficients; Fraction arithmetic checks every
primal inequality, every dual inequality and equality of their objective
values. The saved dual uses18 face constraints,169 extension constraints and
the D29 vertex constraint. Its value is154.

There is also a simple optimal primal rule: weight1 on connected marked pairs
at graph distance at most6 and zero otherwise. The checker independently
evaluates this rule on all face and extension constraints and all rank29
forests. Its maximum is154, attained by D29.

The obstruction comes from imposing extensions without their arithmetic
context. Some cores are sealed against precisely those extensions. This
calculation therefore closes this conservative formulation, while the
[forest-sensitive partial rule](fresh_sqrt5_forest_sensitive_weights_20261015.md)
uses that missing context. The optimum154 is a method diagnostic, not a
dimension theorem. The complete geometric interpretation has not received
a separate independent review.

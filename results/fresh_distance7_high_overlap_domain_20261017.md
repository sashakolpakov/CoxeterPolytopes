# Complete necessary base-overlap domain at distance seven

The [source](../scripts/fresh_distance7_high_overlap_domain_20261017.py)
and [certificate](fresh_distance7_high_overlap_domain_20261017.json)
enumerate the complete base-identification domain for two high minimal
opposite-side cores under absence of a marked path of length at most7.
There are331,593 label-consistent identifications. The necessary inverse
and actual-edge tests leave6,691 base contexts with348 literal core specifications.

Exactly2,614 surviving contexts are the preceding distance8 domain;
the additional4,077 have minimum marked distance8. This identifies the
new geometric work needed to lower the finite path threshold. It does
not supply quadrilateral coverage or a new dimension bound.

## Enumeration completeness

The inherited complete high-core catalogue gives384 marked contexts.
The two chosen marks have finite deletions, and every node is at graph
distance at most7 from their pair. If a left node at distance d and a
right node at distance e are identified, the union has a marked path
of length at most d+e. Hence a context with no such short marked path
requires d+e>=8. Both identified nodes are white: a black node has
distance0 and the other distance is at most7.

For each ordered context pair, the source visits the eligible left
nodes in increasing order. At each node it either omits it from the
identification or maps it to one unused eligible right node. It tests
every edge and nonedge against all prior assignments. Consequently
every nonempty label-consistent partial injection occurs exactly once.
No size cutoff is imposed.

|Identified white nodes|Consistent identifications|
|---:|---:|
|1|129056|
|2|128028|
|3|57949|
|4|14238|
|5|2105|
|6|208|
|7|9|

## Native inverse gates and the new black-neighbour boundary

Deleting either mark from a core leaves a positive-definite Coxeter
diagram. In its positive rational gauge, the physical norm-two inverse
diagonal is the saved inverse diagonal multiplied by the gauge square.
All native nonzero contacts have norm-two magnitude at least1, and the
inverse of a positive-definite Coxeter Gram is entrywise nonnegative.

An unshared neighbour from the other core must contact the shared root.
If that neighbour is white, it lies with the displayed one-black diagram
in an actual edge. That edge is strictly positive definite, so an inverse
projection bound at least2 excludes the configuration.

At the new distance1 ports, this neighbour may instead be one of the
other core's black marks. It then lies in an actual corner, which need
only be positive semidefinite. Equality2 is allowed. The new source
rejects this case only when the inverse bound is strictly greater than2.
It does not copy the predecessor's stronger white-only threshold into
this larger domain.

Remaining base unions undergo an actual shortest-path test and every
strict actual-edge test. The latter uses exact positive diagonal native
congruences where available and otherwise the justified unit-bond lower
comparison. A nonnegative nonpositive-quadratic witness remains an
obstruction when further native cross bonds are added. Thus the rejected
base unions cannot be rescued by unspecified extra cross-core contacts.

|Classification|Contexts|
|---|---:|
|Native inverse rejection|313937|
|Actual edge not positive definite|10965|
|Survives the base checks|6691|

The surviving minimum marked-distance counts are8:4077,9:1665,10:637,
11:198,12:78,13:27,14:9. Every old distance8 survivor occurs literally
in the new set.

## What remains open

The surviving base unions have not had all additional cross-core edges
enumerated or excluded. Their arbitrary supporting-white extensions,
whole marked components and retained finite/ideal pair images are not
yet classified. The low-core and disjoint-core domains also require
their own distance7 analysis. The separate finite7 triangle theorem
does not fill these quadrilateral obligations.

Accordingly this package is a complete necessary-domain diagnostic,
not a finite7 payment theorem. It records every survivor and a
deterministic digest of the full ordered classification. Replaying it
repeats the enumeration and all exact rejection tests; a count or hash
alone is not offered as a proof of geometric completeness.

Generationa3fd36/72929 and a separate saved replay669307/17587 completed
exit0/PASS. A separate full source and inherited-helper read passed the
partial-injection completeness, native normalization, strict-edge versus
PSD-corner distinction and conservative scope. That read did not repeat
the enumeration. The complete report also passed that review. All eight
direct source bindings, two local report links and Python syntax pass.
Source, certificate and report are final/frozen within this stated scope.

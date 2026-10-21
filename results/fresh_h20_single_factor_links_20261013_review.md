# Independent single-factor section review

The [separate SymPy checker](../scripts/fresh_h20_single_factor_links_20261013_review.py)
imports no producer or repository arithmetic helpers. It reconstructs the
literal form and 210 roots, and checks all 156 saved ambient face points
against every wall, including positive definiteness of their defining
Gram matrices. It symbolically recomputes the section and ridge Gram
identities for all 79 leaf-sum planes and checks every rank-two auxiliary
configuration. Its [saved output](fresh_h20_single_factor_links_20261013_review.json)
binds the reviewed producer and packet hashes.

The continuous argument in the [report](fresh_h20_single_factor_links_20261013.md)
was separately reviewed. Positive combinations of the saved timelike
points stay in the same future cone, and their zero walls are exactly the
intersection of endpoint zero sets. The projected rank-two Gram
determinant is positive throughout the sum-plane interval; the positive
off-diagonal entry therefore gives a genuine obtuse ridge.

For C2, the normalized three face points realize the three vertices of
the coefficient triangle. A strictly interior weighted mean exists
precisely in the stated parameter interval. If xy<=0, x+t and y+t have
opposite signs there, including the zero-coefficient cases. If xy>0,
orientation reversal and swapping reduce to x>=y>0; the second pair of
faces gives the stated obtuse ridge throughout -x<t<0. The independence
argument is valid: dependent projected normals would force the section
normal into their original span, contradicting its nonzero evaluations
on face points on both of those walls.

For A2, every mixed-sign coefficient triple has the stated form after
orientation reversal, including one zero coefficient. The two auxiliary
face evaluations have opposite strict signs. The exact formula
q-4t²=4((x-y)²+xy+t(x+y))>0 places the actual squared cosine strictly
between 1/4 and 1/2. Thus no untested angle labels or parameter boundaries
remain in either rank-two argument.

The producer's distinctness tests were reviewed as well: each normal
plane has a unique null line, and at the same null line the incident
pairing support distinguishes the new planes from old difference planes
and from each other. The independent checker does not duplicate those
distinctness counts. The arbitrary higher-rank single-factor problem
and reflection subdivisions of these new sections remain outside scope.

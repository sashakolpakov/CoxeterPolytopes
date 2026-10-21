# Internal review of the twenty-ring central-cap obstruction

Root independently read the
[proof](fresh_nonarithmetic_central_ring_contact_gate_20261017.md) and
[source](../scripts/fresh_nonarithmetic_central_ring_contact_gate_20261017.py).
**PASS for the specified twenty exact recipe rings and one central facet
with four finite endpoints.** This is a conceptual and source review;
it is not a second implementation of the interval engine or another
35,840-expression arithmetic replay.

The low-label reduction is complete after the seven retained-support
profiles: the four ring projections span the residual four-space, so four
prescribed products determine one vector. Its norm equation is therefore
necessary without any numerical acceptance threshold. Every tuple in
{2,3,4,5}^4 is included for each support profile and ring.

For the high-label tail, a bond of magnitude at least sqrt(3) cannot share
an endpoint with another nonright bond in a positive-definite corner.
The stated three-by-three determinant has the correct negative mixed
term. Consequently both endpoints have zero A2 contacts and the adjacent
ring contacts must be right. Adjacent high central bonds are impossible;
the literal eligible lists also exclude an opposite pair. Thus there is
at most one high contact. The other three contacts are exactly the two
right angles and one enumerated label 2 through 5.

Those three products leave a one-dimensional affine solution line. The
script certifies a nonzero positive kernel square and enumerates both
real norm branches, or proves the required residual square negative.
It does not silently drop a zero or unresolved branch: those cause an
assertion failure. A certified product above -sqrt(3) then excludes every
finite high label at once, without an angle-recognition cutoff.

Numerical coordinates select which exact norm branch defines a ring
recipe. Either branch is still an exact vector obtained from the stated
contacts and norm equation. This selection does not certify all ambient
contacts of the ring, and neither the proof nor this review claims that
it does. The result is an obstruction for these recipe objects; it does
not exclude changed ring parameters, the full predecessor search bank,
several new central walls or a different completion topology.

The source uses the frozen outward rational interval engine and prior
exact support reduction as explicit inputs. Its generation and fresh
saved-data replay are reported separately in the producer package.
No improved chamber, finite-volume completion or new dimension follows.

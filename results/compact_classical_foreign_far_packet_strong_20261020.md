# Foreign contacts make a classical far-end packet strong

20 October 2026. Root paper proof, awaiting a separate full review.
No scientific computation, Lean, TeX or accepted-source edit. This is
an interface for extending the existing whole-component arguments;
it is not a dimension bound.

## Statement and conventions

Let p be a farthest vertex from a strict interior point of a compact
Coxeter polytope. Suppose one whole component is B_r or D_r, in the
standard unit-root coordinates. For a long-arm index i put

    omega_i=sqrt(2)(e1+...+ei), q=sqrt(2i), h_i=x.omega_i.

Assume 5<=i<=29 and that the actual opposite c in direction i has,
within this component, precisely these nonzero contacts:

    -c.n_i=gamma>=cos(pi/5),   -c.n_1=1/2.

All other contacts in that component are zero. Contacts in arbitrary
orthogonal finite components are allowed; let Q be their complete
unit inverse-Gram cost. If Q>0, then

    h_i>sqrt(2i)-61/100.                            (1)

Here i is on the long chain, so deleting i leaves the path
A_(i-1) on its left; no old fork or terminal-4 root is among this path.
The theorem does not restrict r to29. In applications to compact
polytopes the published ambient ceiling29 bounds all indices.

## Actual Schur data, including every foreign component

An actual finite replacement deleting i retains all old foreign
components. Every nonzero foreign contact is finite and native, so
its contact square is at least1/4. In a finite unit Gram with
nonpositive off-diagonals, inverse entries are nonnegative and each
diagonal entry is at least1. Consequently Q>=1/4 if Q>0.

The sole contact with the retained A_(i-1) is simple3 at its far
endpoint1. The inverse endpoint diagonal is2(i-1)/i, giving cost
(i-1)/(2i). Thus the actual white residual and outgoing derivative
are

    delta=(i+1)/(2i)-Q>0,
    a_i=(2i*gamma+1)/sqrt(2i).

The actual edge length L_i therefore satisfies

    cosh L_i=(2i*gamma+1)/sqrt(i+1-2iQ),
    cosh L_i/q >= R(i)=(2i*gamma+1)/sqrt(i(i+2)).    (2)

Increasing Q only increases the length. This calculation does not
assume that foreign contacts vanish or that a metric facet is Coxeter.

## Exact finite-range comparison

Write c0=61/100 and define

    T(q)=[1+(1-c0/q)^2]/[q*(1-(1-c0/q)^2)]
        =1/c0-(q-c0)/[q*(2q-c0)].                  (3)

For q>=sqrt(10), T is increasing: the derivative of the subtracted
fraction has numerator -2q^2+4c0*q-c0^2<0. This follows already from
q>3 and c0<1. Since gamma>809/1000, R(i) is strictly above

    R0(i)=[(809/500)i+1]/sqrt(i(i+2)).

The derivative of R0 has numerator ((309/500)i-1), divided by a
positive denominator, so R0 is increasing for i>=5.

The following rational bounds prove R(i)>T(q) on the entire required
range; the q upper bounds are strict except at i8.

| i range | q upper bound | R(i) lower bound | T(q) upper bound |
| --- | ---: | ---: | ---: |
| 5..8 | 4 | 1.536 | 1.530 |
| 9..12 | 5 | 1.563 | 1.550 |
| 13..29 | 7.62 | 1.577 | 1.577 |

For the three R bounds, square the positive quantities at i5,9,13.
The positive cleared differences are, respectively,

    9.09^2-35*(1.536)^2=.052740,
    15.562^2-99*(1.563)^2=.321913,
    22.034^2-195*(1.577)^2=.546001.

For the T bounds, its numerator and denominator before division are

    2q^2-2c0*q+c0^2,     q*(2c0*q-c0^2).

At q4,5,7.62 they are respectively

    (27.4921,18.0316),
    (44.2721,28.6395),
    (107.2045,68.003166).

Multiplying the denominator by the displayed upper bound gives a
strictly larger numerator in every row. All decimals here are exact
rationals; sqrt(24)<5 and sqrt(58)<7.62 follow by squaring.
The elementary native-angle estimate cos(pi/5)>809/1000 follows from
cos(pi/5)=(1+sqrt(5))/4 and 2.236^2<5.

By (2)--(3), cosh L_i>[1+t^2]/[1-t^2], where t=1-c0/q is in(0,1).
The half-tangent identity then gives tanh(L_i/2)>t. The actual
farthest-edge inequality h_i/q>=tanh(L_i/2) proves (1).

## Inputs and scope

The actual-edge convention and farthest comparison are those of
`compact_farthest_vertex_energy_obstruction_20261020.md`. The
whole-B weak-packet notation is in
`compact_whole_b_bulk_weak_packet_and_endpoint_pressure_20261020.md`.
Finite native contacts and nonnegative inverse finite-Gram entries
are the same explicit classical inputs used there.

This proof settles the foreign-contact branch of the far-end3 packet
only. It does not classify all actual opposites, clear foreign contacts
from isolated or near-end packets, or transfer the whole-B23 theorem
on its own. No asymptotic conclusion beyond index29 is claimed.

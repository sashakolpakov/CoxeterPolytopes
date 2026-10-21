# Executed28-wall check and a new exact spacelike escape

19 October 2026. **Generation and separate replay actually passed.**
This is a new execution/consequence note. All pinned root-construction
reports, the28-wall source and its preparation note remain unchanged.
The coordinator ran the scientific jobs; no such job ran in this
authoring lane.

**Conclusion.** The28-wall H19 seed still has infinite hyperbolic
volume. Its28 actual facets, two full compact H8 facets, at least86
finite vertices and nine finite initial edges remain valid. An edge
walk beyond the initial vertex found the exact spacelike point below.
Its polar line contains no lattice reflection direction. This does
not exclude different completing mirrors or prove nonreflectivity.

## Actual execution and bounds

Source `scripts/check_h19_twenty_eight_wall_seed_20261019.py`, SHA256
`f9ccdb1dbb5fe6a810bf456be56b3f15d1b83be69fd0d080931713eae3738600`,
checked all378 wall pairs,28 actual-facet projections, four distinct
finite-vertex witnesses, and all nine initial p-edge endpoints.
The bounded walk stopped on a feasible spacelike endpoint after19
finite vertices had been discovered. Its64-vertex cap was not reached;
its queue did not close. No full vertex or cone census occurred.

Generation/replay took0.544/0.546 seconds, with peak process-tree
RSS56,721,408/57,032,704 bytes. Both started at187.47GiB free disk;
available RAM was6.19/6.18GiB. The wrapper owned the exclusive lock,
59-second/1GiB limits and owner-death control; the child retained
its45-second alarm. These are actual ledger measurements.

|Executed output|SHA256|
|---|---|
|`results/h19_twenty_eight_wall_seed_20261019.json`|`3996be94b07f251685a15650e50321456a4f3ed5548e0ad927d94abba1525b83`|
|`results/h19_twenty_eight_wall_seed_generation_20261019.json`|`1b9ca4b6201d2ad4c41dece5a7ea9b380e9b31f0586ea02c42b0bb3aac1ea281`|
|`results/h19_twenty_eight_wall_seed_replay_20261019.json`|`e10086668c070c00fa33da6f5c2c9ecb74a552bf6b8f95a071e11087da4a7a1b`|

## A compact exact formula for the returned witness

Keep phi²=phi+1, b=4+phi, D=1+5phi and omega_j=C_E8^{-1}e_j.
The walk starts in the affine slice x=1. Its vertex2 is

    V=(1,-phi/2,(phi³/4)*omega0),
    V²=(3-2phi)/4<0.

Its nine active walls are old0,old2,...,old8,z. Freeing old2 gives
the exactly normalized edge direction

    d=(0,0,omega1-(7/4)*omega0).

Here d_x=0, its product with old2 is-1, and its products with the
other eight active normals are0. The first positive wall hit is
old18 at parameter1/5. Thus the returned feasible point is

    E=V+d/5.

The identities omega0²=4, omega1²=14, omega0.omega1=7 give

    d²=7/4,  (V,d)=0,
    E²=(3-2phi)/4+7/100
       =(41-25phi)/50
       =57/100-sqrt(5)/4>0.                        (1)

For positivity it suffices that phi<13/8<41/25. This simple formula
is the same point as the full ten-coordinate executed witness.

## Complete28-wall feasibility ledger

Expanding the certificate's exact products into the native phi basis
gives the following exhaustive table. The products not zero are
strictly negative.

|Wall|Product with E|
|---|---:|
|old0,old3,old4,old5,old6,old7,old8,old18,z|0|
|old1,old12|(1-5phi)/10|
|old2|-1/5|
|old9|(5-9phi)/10|
|old11|-1/2|
|old13|-3phi/10|
|old14|(3phi-5)/5|
|old15|-7(1+phi)/10|
|old16|-(5+phi)/10|
|old17|(phi-4)/10|
|old19|-(5+4phi)/10|
|r|-D/2|
|k|-(17+30phi)/10|
|w|-(12+22phi)/5|
|a|-(1+8phi)/5|
|d-root|-(18+26phi)/5|
|s|-(33+36phi)/10|
|q|-9/10|
|h|-(36+56phi)/5|

The label d-root distinguishes the earlier wall d from the new edge
direction d. No untested wall is omitted from this volume argument.
The future time orientation is fixed by(p,E)=-D<0. The common strict
timelike center c was checked against all28 walls. A small positive
perturbation of E towards c is still spacelike and strict on every
inequality. The segment from c then crosses the future light cone
inside all halfspaces, exposing an open ideal patch. Thus the chamber
has infinite volume; no global edge-walk completeness theorem is
needed for this negative conclusion.

## The new polar line is not a lattice reflection direction

The rational field norm is N(a+bphi)=a²+ab-b². From(1),

    N(E²)=31/2500.

The unchanged H19 primitive norm-alphabet theorem allows canonical
reflection norms2 or2b, whose rational norms are4 and76. The two
ratios are

    N(E²)/N(2)=31/10000,
    N(E²)/N(2b)=31/190000.

Each has odd valuation at31, so neither is a rational square.
Multiplication of E by any scalar in K changes N(E²) by the
rational square N(scalar)². No such scaling can therefore give
either canonical reflection norm. There is no missed real scaling:
E has first coordinate1, so any K-coordinate root parallel to E
has its scalar in K. This excludes the polar line at every height
and unit scaling, conditional on that inherited complete norm alphabet.

A new mirror cutting E need not be polar-parallel. The previous
nonpolar construction already demonstrates why this local norm
obstruction must not be mistaken for whole-lattice nonreflectivity.
The next endpoint white support and its necessary finite profiles
will be analyzed separately. No compact H9 or dimension record
has been obtained.

# Marked narrow transitions do not remove the baseline38 obstruction

The exact max-to-max transition test retains the obstruction at **38**, above
the narrow20-dimensional facet threshold **35**. The marked type
`A1+D20`, with its isolated `A1` as the facet root, has an admissible
minimality14 arithmetic return at every one of its20 edge ports, preserving
the marked port orbit. The separate marked type `D21` does lose five ports.

Artifacts: [source](../scripts/fresh_sqrt5_marked_maximum_transitions_20261016.py)
and [literal certificate](fresh_sqrt5_marked_maximum_transitions_20261016.json).
The two extremizers and the threshold are established in the frozen
[marked-facet diagnostic](fresh_sqrt5_marked_facet_barrier_20261016.md).

## Exhaustive local test

The source uses norm2 finite diagrams and retains the facet root literally in
the common20-root support. At each of the20 ports of each extremizer, it tests
all1,351 subsets of common roots of sizes0 through3 as the new root's3-bond
attachments. These are complete: both candidate endpoint diagrams are simply
laced, with maximum degree3. Exact marked tree isomorphism accepts66
reattachments among54,040 subsets:43 from `A1+D20` and23 from `D21`.

Eliminating the positive common Gram matrix leaves
`[[s,g-q],[g-q,t]]`, with rational positive `s,t`. Every accepted attachment
has `q-sqrt(st)>=0` and `q+sqrt(st)<=2`. Conjugate positive definiteness and
minimality14 therefore require `0<g'<2` and `-14<g<=0`.
Every row saves a cycle containing the framing product `g`, whose other
factors are3-bond products `-1`; the arithmetic cyclic-product condition
thus puts `g` in `O_K=Z[phi]`. This derives, rather than assumes, integrality
of the normalized framing entry.

Writing `g=a+b*phi`, the two embedding bounds give the complete42-element box
`-4<=a<=1`, `-7<=b<=-1`. The program tests every scalar exactly. Any surviving
entry is either dotted or `-phi` (a5-bond); candidates strictly between `-1`
and0 are rejected as non-Coxeter products. Physical Lorentz signature and
conjugate positive definiteness are verified by the exact Schur determinants.
The55 feasible reattachments have109 admissible scalar choices in total;
one explicit pair of future finite endpoint rays is saved per reattachment.
The open edge segment lies in the retained facet.

## What survives

In `A1+D20`, label the marked isolated root0 and shift the canonical `D20`
roots to1 through20. At the fork-center port2, copied attachments
`{1,3,20}` fail, but the reversed long-arm attachment `{1,19,20}` succeeds:

`s=t=1/18`, `q=19/18`, `g=-2-5*phi`, `1<g'<10/9`.

The isolated marked root0 remains isolated throughout. Both endpoints are
the same marked type, and the framing roots occupy the same marked root
orbit. Every other port also has such a narrow same-port return.

In `D21` with its long-arm endpoint19 marked, ports1,2,4,5,6 have no narrow
return to either extremizer. Marking does invalidate many unmarked long-arm
reversals, but that improvement does not remove the other extremizer.

Put vertex mass1 on the marked `A1+D20` type. Giving each saved same-port
return undirected edge mass1/2 supplies incidence1 at every port. This is an
exact balanced finite-state port-flow relaxation. Its central-E8 incidence
inequality remains `0<=0`, and its baseline cost remains38. Consequently
these marked narrow edge constraints, even with that incidence inequality,
cannot force the required upper bound35.

## Verification and limits

Generation and a fresh saved-output replay both pass. Arithmetic is exact
`Fraction` arithmetic in `Q(phi)`; no optimizer or numerical tolerance is
used. Input source manifests are recursively checked. Replay command:

```sh
python3 scripts/fresh_sqrt5_marked_maximum_transitions_20261016.py --verify
```

These are separate local arithmetic edge witnesses. They need not realize
one simultaneous star, one globally compatible narrow facet, or a
finite-volume polytope. No dimension bound or new polytope follows. The test
settles this particular proposed restriction; a useful next restriction must
couple multiple edges or change the baseline/payment family. No new LP or
additional search family was launched.

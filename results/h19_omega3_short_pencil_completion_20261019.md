# All-height obstruction for the omega3 short pencil

19 October 2026. The bounded exact source has now actually passed generation
and a separate replay under the coordinator's single-lock resource wrapper.
The following independent paper argument explains its negative outcome and
removes the height cutoff for this pure pencil. It does not exclude walls
with nonzero contacts to the eight fixed walls of the initial edge.

## Paper theorem

Put O=Z[phi], phi^2=phi+1 and D=1+5phi. There is no solution

    10t^2-Du^2=1,  u,t in O.

Indeed N(D)=1+5-25=-19, so (D) is a prime ideal of index19.
The quotient is F19 under phi maps to15: 15^2-15-1=209=11*19
and1+5*15=76=4*19. Reducing the equation gives10t^2=1, or
t^2=2. But2 is a nonsquare modulo19, since2^9=512=-1 modulo19.
This contradiction holds at every height, with no conjugate bound.

It even rules out a field square t^2=(1+Du^2)/10 for integral u:
the right side is integral at (D), so any t in Q(sqrt5) with this
square is also integral at that prime and has a residue in F19.
The same nonsquare contradiction applies. Thus zero *field-square*
solutions, before the integral-t filter, were to be expected.

For the current lattice the pure omega3 short reflector has the form
r=(2u,-phi*u,t*omega3), omega3^2=20, and its norm-two condition is
exactly this equation. This excludes the entire integral norm-two pure
pencil. The already reviewed ramified-prime argument excludes its
canonical norm2b counterpart. The omega2 pure pencil was independently
excluded in both canonical norm classes by the separate F9 argument.
Canonical norm-alphabet completeness remains an inherited lattice fact.

These conclusions do not exclude a short wall meeting a fixed edge wall
nonorthogonally. They neither prove nonreflectivity nor settle any general
completion/finite-volume problem. The existing26-wall chamber remains
infinite-volume by its separate explicit escape proof.

## Executed bounded evidence

Source `scripts/search_h19_omega3_short_pencil_20261019.py` remains
SHA256 `02750aed919d5d507f6e32816718484e81cf8f5f1570ad783295aa591ae1edbd`.
Its source-only preparation, complete coefficient-box proof, exact square
test and centre construction are recorded separately in
`results/h19_omega3_short_pencil_preparation_20261019.md`.

The actual certificate checked all12,080 coefficient pairs in the
conservative box and found537 values satisfying0<u<=300 and|u'|<=2.
It found zero field-square solutions and hence zero integral candidates
or surviving roots. Consequently the26-wall product/centre checks have
no candidate instances in this run; the paper theorem above, rather
than this vacuous downstream test, supplies the all-height exclusion.

|Evidence|SHA256|
|---|---|
|`results/h19_omega3_short_pencil_20261019.json`|`d8d048b18c408bcf71f80aea358aaf413b27c99793eb8a8b1a2773fec653aa11`|
|`results/h19_omega3_short_pencil_generation_20261019.json`|`b6c2590c434f6529a45dbefe6fc92dbf2068717b4c04f7fc058353175ef7aa0d`|
|`results/h19_omega3_short_pencil_replay_20261019.json`|`0d69c3e8d5451767f155caaf7a423c27023009b5117f449cc28f4f64de199a0f`|

Generation and replay both report `PASS_BOUNDED_RESEARCH_CHECK`, each
in0.258 seconds, with process-tree peaks56,082,432 and55,820,288 bytes.
Start resources were186.67GiB disk and6.84/6.66GiB RAM respectively.
The wrapper owned the single lock,59-second deadline and1GiB RSS cap;
the child added its45-second alarm. No second search, optimizer, Lean
process or build was run by this review lane.

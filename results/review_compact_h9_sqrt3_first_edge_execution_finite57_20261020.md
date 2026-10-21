# Independent interpretation of the generated first diagonal endpoint

20 October 2026. **PASS for the explicit root and local finite corner.**
Written algebra/source inspection only; this reviewer did not execute the
checker or reproduce its finite enumeration.

Evidence inspected:

* `results/compact_h9_sqrt3_diagonal_first_edge_20261020.json`, SHA256
  `10c8afa329dd9328c924c17562ba7fd33500dde2b4686f5e6eccb6f876d05eb3`.
* `results/compact_h9_sqrt3_first_edge_generation_20261020.json`, SHA256
  `67ac03539fece90725aacc2184e27393e1152e9c271c57550f6f2cd8428fbd14`.

The certificate binds worker56149b65 and route16732ee8, already reviewed in
`results/review_compact_h9_sqrt3_diagonal_first_edge_finite57_20261020.md`.
The generation ledger records PASS in0.104 seconds, peak54,902,784 bytes,
with the wrapper owning the lock. Its start gates were180.81GiB free disk
and6.03GiB available memory. The separately attempted replay was reported
pre-launch deferred at4.4GiB available memory: **replay remains pending**;
there is no replay success or replay ledger asserted here.

## Independent exact checks

Write s=sqrt(3). The unique accepted record is isolated profile5, power4,

    a=6230+3597s,       b=2733+1578s,
    r=(a;b,b,b,b,b,b,b,b,b),       q=9+5s.

Direct multiplication gives

    a^2=77628127+44818620s,
    b^2=14939541+8625348s,
    9b^2-s*a^2=9+5s=q.

Both embeddings of q are positive and its field norm is6. The root is
primitive: a common nonunit factor would make the square of its nonunit
integer field norm divide6, which is impossible. Here the norm-Euclidean
ring property supplies the common-factor description.

The two coefficients

    2s*a/q=1223+705s,       2b/q=309+179s

are integral. Thus `2<r,x>/q` is integral for every lattice vector x,
and the root reflection preserves the lattice directly, independently of
trusting the stored norm or a numerical approximation.

The endpoint is

    p=(1;lambda,...,lambda),
    lambda=s*a/(9b)=(339+536s)/2889>0.

It lies on all eight old A8 walls and on r, while its product with the
released wall -e9 is -lambda<0. Its square is exactly

    p^2=-s*q/(9b^2)<0.

It is therefore a future timelike endpoint, strictly before the diagonal
light cone. No near-zero floating-point sign is used. All eight products
of r with the retained simple roots vanish. Their Gram block is the A8
Cartan matrix and the ninth block is q, so the incident Gram determinant
is9q>0 and the full-rank finite corner is A8+A1. The released old wall is
ultraparallel to r: its raw product is -b and b^2>q.

The saved representative `(z0,X0)=(9+6s,6+3s)` independently satisfies
`z0^2-s*X0^2=81+45s=9q`. The stored z=9b and X=3a have the same identity.
The certificate reports138 z-box elements,80 X-box elements and only this
one accepted root across the complete ten-row reduction. Consequently
there is no simultaneous-choice ambiguity in this output.

## Meaning and remaining boundary

The finite enumeration's completeness and minimum interpretation use the
previously reviewed primitive norm inventory, finite A8-extension theorem,
relative-unit strip, and reflection-chamber argument. The explicit algebra
above independently checks the returned object; it is not a second
implementation of that enumeration. Conditional on those written inputs,
the successful generation identifies the first finite endpoint of this
particular initial edge. The rank-nine corner is full dimensional locally;
no larger ambient independent incident root is needed.

This is a local result in the integral reflection arrangement. It neither
closes any other initial edge nor supplies a complete acute root list or a
global extreme-ray certificate. No finite-volume or compact H9 example is
claimed. Separate replay and any next-edge computation remain separate
verification steps.

# Fixed R5 chamber-period verification contract

20 October 2026. SOURCE ONLY. No import of the worker or scientific
execution, generation, replay, Lean or typesetting is claimed here.

The single worker is
`scripts/check_golden_dyadic_r5_chamber_period_20261020.py`; its fixed
certificate is `results/golden_dyadic_r5_chamber_period_20261020.json`.
The final source hash is supplied in the source-review handoff, avoiding
a circular source/preparation hash. No separate coordinate JSON is used.

The mathematical input is
`results/golden_dyadic_r5_chamber_period_candidate_20261020.md`, SHA256
`e2aeff93b42f87673e55399804df28a9eee3713ac2e7cf064c6f885e5bfce671`.
The source binds the unit-roof and actual-N proofs and both complete
written reviews (root74de10e5 and outside-author81bea9b5). It also binds
the original immutable25-file M6/H3 packet: freeze SHA256
`22dd341c2d79651742aa93ebdc8675e8fc55025fdd0cbcb678b128f7e59aa37b`.
Every enclosed file's size and hash is checked. The own source and all
scientific/control input hashes are checked again before a certificate
can be returned. Frozen sources and prior certificates are never changed.

## Exact scope

The worker reconstructs the integral basis
`(V0,e1,...,e6,a3,s,wP)` with Gram
`diag(-phi,I6,[[4,-2,-2],[-2,4,0],[-2,0,phi+2]])`.
It checks equality with the ORIGINAL lattice
`<-phi>+I5+D4(O)` in both directions, using a standard full D4 basis.
It independently projects e1/e5 off the new unimodular H4, verifies
P1/g and their explicit original coordinates, and computes the actual
rank8 positive projection S directly from its Gram.

All six fixed finite continuation profiles are retained. Their B3
white projections, norms and contacts are checked, including the B2
row's separate `c^2=1` and `kappa=0`. The proven complete coefficient
rectangle `m=-4,...,4`, `n=-2,...,1` contains36 pairs per norm class,
108 fixed strip tests altogether. Exact physical/conjugate inequalities
give the three negative-A lists. Every earlier case is rejected by
nonintegrality of b^2 or b^2/4. At a=0 the rational-square test excludes
sqrt2 and sqrt(1/2) over Q(sqrt5); 2g is explicitly imprimitive.
Both primitive simultaneous roots remain: g=Q+beta6. Their tie and
oriented inequalities prove the written nonsimplicity of g; no root
is dropped by a filter chosen to force uniqueness.

The full period D is reconstructed from the p0 dyadic glue, then the
independent V0/u/beta relation and F formula. All100 entries are compared
with the paper's array parsed by a restrictive integer/linear-phi
regular expression, never eval. The worker verifies Gram preservation,
both integral inverse identities, every positive full9-root germ image,
both future finite apices and trace10+11phi. It checks integral
reflection coefficients for the two complete germs and g/Q.

No arbitrary height bound, Pell search, cone census or global root
enumeration occurs. Completeness of the six finite weighted profiles,
full H4 seal, actualness of imported walls, the positive-projection
theorem and chamber-germ implication remain explicit written inputs.
The certificate verifies their fixed arithmetic interfaces, not those
general geometric theorems by program. The resulting obstruction is
ONLY for the literal R5 lattice, not every lattice in its similarity
or commensurability class, and not a new dimension ceiling/example.

## Execution boundary

Root alone may schedule one bounded generation and one separate fresh
replay after source review. Use the existing parent-lock wrapper with
`python -B`: it owns the shared lock, resource gates,1GiB RSS limit,
owner/time/watchdog and59-second outer lifetime. The child requires
`COXETER_CHECK_OWNER=os.getppid()`, owns no inner lock, and sets a
45-second alarm. It has exclusive `--write` creation, fresh exact
`--verify`, a128KiB certificate hard cap and a brief metadata-only
stdout result. Any assertion, resource/timeout failure or changed
source stops before certificate creation. No expanded retry or new
search is authorized by this preparation contract.

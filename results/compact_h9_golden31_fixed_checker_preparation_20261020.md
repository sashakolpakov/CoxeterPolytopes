# Fixed golden31 unit-bisector checker contract

20 October 2026. SOURCE ONLY. No worker import, scientific execution,
generation, replay, Lean or TeX run is claimed.

The single worker is
`scripts/check_compact_h9_golden31_unit_bisector_20261020.py`.
Its fixed certificate is
`results/compact_h9_golden31_unit_bisector_20261020.json`.
The final worker hash is supplied in the source-review handoff, avoiding
a circular worker/preparation pin.

The written proof is
`results/compact_h9_golden31_unit_bisector_preparation_20261020.md`,
SHA256 `7259c05e994fd48d8bb5f70c172e98f9c0120bbfa13fc58f3a6b08174e5f9ea1`.
Root's complete independent mathematical review is
`results/review_compact_h9_golden31_unit_bisector_root_20261020.md`,
SHA256 `4cf76e054d1b8461a46d86fce9421a0b7139b741f6f45938761a6a15a71a80f9`.
An outside-author source review is required before root schedules execution.

## Fixed arithmetic scope

The worker uses only the Python standard library and its own exact
Fraction-pair implementation of a+b*phi. Its sign and conjugation routines
are exact. It reconstructs the native E8 Gram from the seven fixed edges,
compares it to the calibrated fixture, verifies all eleven root products
and their native indexing, and derives the three fundamental weights used
in the proof by exact inversion. No project/scientific module is imported.

It reconstructs B,r,c,g,t and the two old/new finite binary endpoints;
checks the unit binary basis change, complete necessary norm ideal and
finite-angle prime31 seals, the first-edge norm bound, both germ Cauchy
norm seals and the unit size contradiction, and the E7 Gram/rank and
orthogonal complement. It independently builds S from native columns,
W from the E7 longest-action formula, H from point inversion, and checks
all involutions, exact Gram identities, integral adjoint inverses, root
images, future signs and the four trace contractions. The product has
the checked trace28+42phi, not a floating approximation.

There is no loop in height, Pell units, candidate mirrors, or cones.
The loops have only fixed matrix dimensions or the fixed finite-angle
list. The certificate stores the fixed matrices and arithmetic seals.

The golden PID/unit facts, primitive-reflector ideal argument, completeness
of the finite-angle list, prime-isolation geometry, actual full Bugaenko
facet, positive actual chamber signs, local chamber-germ theorem and group
coset/covolume implications remain explicit written inputs. The checker
does not enumerate all mirrors or prove those geometric/classical inputs
by program. Its conclusion concerns ONLY the literal golden31 lattice;
it is neither the older silver31 host nor a uniform infinite-family or
commensurability-class exclusion, new example or dimension ceiling.

## Ownership, resources, and binding

Root alone may schedule one bounded generation and a separate fresh
replay after full source review, through:

    python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock -- scripts/check_compact_h9_golden31_unit_bisector_20261020.py --write
    python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock -- scripts/check_compact_h9_golden31_unit_bisector_20261020.py --verify

These are a contract, not executed commands. The existing wrapper owns
the sole lock, resource gates, owner monitor and59s/1GiB outer limits.
The child requires COXETER_CHECK_OWNER=os.getppid(), owns no inner lock,
sets a45s alarm, requires -B and explicitly rejects disabled assertions.
It hashes all pinned proof, fixture,
inherited mathematical and control files plus itself before checking
and again before returning a certificate. The certificate is capped
at128KiB, --write refuses overwrite, and --verify snapshots the saved
certificate with a bounded read of at most128KiB+1 BEFORE recomputation.
After recomputing, it compares the snapshot with the exact new payload
and with a second bounded saved-certificate read, rejecting oversized or
changed bytes. Stdout is one brief metadata line.
Any failure leaves no completed certificate; no expanded retry/search
is authorized by this contract. Frozen files are never changed.

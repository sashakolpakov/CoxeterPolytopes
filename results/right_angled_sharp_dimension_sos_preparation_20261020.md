# Source preparation: universal tensor-square algebra

20 October 2026. **Prepared, not executed.** Python AST parsing passed;
the checker module was not imported, its body was not called, and no
certificate or execution ledger was generated in this lane.

Worker: `scripts/check_right_angled_sharp_dimension_sos_20261020.py`,
154 lines, SHA256
`bd0ed1c4c1b7e18a219a2a8aca2b178dbe1d559d7e9ec90ad609150aedf8b2c3`.
It pins the reviewed full geometric proof
`results/right_angled_dimension8_descending_height_proof_20261020.md`,
SHA256 `9470fa1fa3a973a12c412fbc1257537655536dba4ffc4e9ac23c320615babdbe`,
and streams both that input and its own source before and after its work.

## Exact checked statements

The ring is `Q[u1,...,u6,v1,...,v6]`, represented by exact monomial
dictionaries with Fraction coefficients. All 36 rank-two tensor entries
and all 216 rank-three entries are constructed directly with Kronecker
deltas, then their products are expanded. The checker independently
reconstructs `su=sum u_i^2`, `sv=sum v_i^2`, `t=sum u_i*v_i` and checks

    <A(u),A(v)>=t^2-(su+sv)/6+1/6,
    <T(u),T(v)>=t^3-(3/8)t(su+sv)+(3/8)t.

For both formal vectors it also verifies the tensor traces and squared
norms with explicit multiples of the unit-ideal generator `su-1` or
`sv-1`. Thus unit normalization is handled by universal polynomial
remainders, not samples or unverified substitution.

The independently expanded scalar polynomial is

    f(t)=(t+1/2)^2(t-1/4)=t^3+3t^2/4-1/16.

The SOS kernel satisfies the full arbitrary-coordinate identity

    kernel-f=(2-su-sv)(3t/8+1/8),

including its explicit representation in the ideal `(su-1,sv-1)`.
Normalized rank-two/rank-three kernels and their positive Gegenbauer
coefficients are checked separately. The final rational scalar checks
are `f(1)=27/16`, the ratio to the constant coefficient is 27, and four
unit vectors with all mutual products minus one half would have squared
sum minus two.

The saved output records compact canonical polynomial summaries and
hashes, the exact coefficients, contraction counts and source hashes.
The program's full symbolic equalities are recomputed on every replay.
It checks no Lorentz-space construction, face inheritance, published
facet minimum, spherical-code census or hyperbolic example.

## Root-owned bounded execution

The standard shared wrapper must own the lock and resource watchdog;
`owns_lock=false`. The exact argv arrays are

    ["python3", "scripts/check_right_angled_sharp_dimension_sos_20261020.py", "--write"]
    ["python3", "scripts/check_right_angled_sharp_dimension_sos_20261020.py", "--verify"]

The optional `--stdout` mode recomputes the same algebra and prints the
bounded certificate, so it too is a job for the coordinator to schedule.
The worker sets a 45-second alarm, exclusively creates its certificate,
refuses an existing output on generation, and caps the result at 256KiB.
It uses only the standard library, with no scientific producer imports,
optimizer, floats, samples or background children. The expected work is
small, but no timing or memory measurement is claimed before execution.

The 40GiB disk/6GiB available-memory start gates and 8/3 abort gates apply,
with one process-tree job at a time, a 59-second parent deadline and a
one-GiB cap. This preparation leaves all earlier proof packets and
unexecuted logarithmic/LP sources unchanged.

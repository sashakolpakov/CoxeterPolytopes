# Fixed exact escape check for the sixteen golden dyadic walls

20 October 2026. Source preparation only. No worker body, generation,
replay, cone enumeration, or new compactness result has run here.
Root alone schedules execution through the shared bounded wrapper.

The fixed inputs are `golden_dyadic_sixteen_wall_caps_20261020.md`,
SHA256 `c288b5de619471772c8381d98afa3996ac8851b4d2d146ad520d6474ee8ba92b`,
and its exact coordinate JSON,
SHA256 `80f18b60003375876f5bae9c122c9466ff1eb1c69fed582d87c7facf50fd5745`.
The JSON uses rational pairs `a+b*phi`, not `a+b*sqrt(5)`.
The worker independently reconstructs all sixteen normals, the original
embedding, strict interior point and spacelike escape, then compares
them with the fixed list. All scientific/control input hashes and its
own hash are bound before and after the body.

## Exact obligations

In the coordinates `(V0;e1,...,e6;a3,s,r)` the form is
`diag(-phi,I6,G3)`, with

    G3 = [[4,-2,0],[-2,4,-2phi],[0,-2phi,4]].

The original lattice is the literal `R5=<-phi>+I5+D4(O)`,
`O=Z[phi]`: the last four coordinates have sum in `2O`.
The worker reconstructs and compares all ten original embedding columns
and their entire Gram matrix. It checks lattice membership of every
normal and checks its reflection on a complete original integral basis
(time, five unit coordinate vectors, four D4 simple roots). Reflection
integrality is consequently not inferred from a determinant or a
K-space similarity. It also verifies the positive definite nine-wall
old Gram and inverse. The original conjugate form `diag(-phi',I9)` is
positive definite; a nonzero K-rational vector cannot be null.

The 120 unordered native pairings are recomputed in both coordinate
systems and compared with the full supplied ledger. Each must be right,
have the exact angle `pi/3`, `pi/4`, or `pi/5`, or be strictly
ultraparallel. Unexpected angles or tangent pairs fail closed. The
explicit strict point q and its sixteen strictly negative products
are independently checked, including

    q^2 = -phi + (383+24phi)/(4*100^2) < 0.

For the old nine normals b_j, write `C=Gram(b)^-1` and
`omega_j=sum_k C[j,k]*b_k`. The corner change is

    X(y)=y0*V0-sum_j yj*omega_j,
    X(y)^2=-phi*y0^2+y_spatial^T*C*y_spatial.

The strict point has corner coordinates `(1,1/100,...,1/100)`.
All sixteen full row vectors in these coordinates are computed, without
any extreme-ray enumeration.

## The single counterexample and its exact support

The fixed witness is

    E_M=(1;phi^-1,phi^-1,phi^-1,phi^-2,phi^-2,phi^-2),
    E_P=-phi^-3*omega1-phi^-5*omega3.

The separate squares are `21-13phi` and `26-16phi=2phi^-6`,
so `E^2=47-29phi>0`. Its original integral multiple is

    2E=(8-2phi;(2phi-2)^2,(4-2phi)^3;
        1-2phi,1-2phi,3-2phi,5-6phi).

All sixteen products are recomputed and must be nonpositive. The exact
active list is `beta1,beta2,beta4,beta5,s,u,w,W,Z`. The worker performs
one pivoted exact rank calculation on these nine full rows and obtains
rank9. Hence E is an actual extreme ray of this pointed corner cone:
the active homogeneous equations cut out its line, and all remaining
inequalities have the recorded correct signs. Extremality is extra
evidence; a feasible spacelike vector alone settles noncompactness.

The strict point and E give the stronger open-boundary obstruction.
For a sufficiently small positive epsilon, `E+epsilon*q` is still
spacelike and is strict for all sixteen walls. Its inner product with
the future V0 is negative. The segment from q to this point therefore
crosses the future light cone in a point strictly inside all the
halfspaces, giving an open ideal-boundary patch. Thus this finite
sixteen-halfspace prefix has infinite volume. This does not exclude
later integral mirrors or the full lattice chamber, and it supplies
no compact dimension-nine example.

## Bounded source-only contract

The single fixed worker is
`scripts/check_golden_dyadic_sixteen_wall_cone_20261020.py`.
It has only `--write` (exclusive creation) and `--verify` (fresh complete
recomputation and comparison). Its certificate is
`results/golden_dyadic_sixteen_wall_cone_20261020.json`, at most128 KiB;
stdout is a short status/count/hash record. There is no DD algorithm,
height search, Pell search, cone census, or automatic larger retry.
The earlier unfinished DD preparation was replaced before any worker
prototype or execution because this fixed witness decides the gate.

The child has a 45 second alarm. The unchanged parent-lock wrapper
owns the shared lock, live owner watchdog, 59 second wall limit and
1 GiB process-tree RSS limit. It runs with `-B`; no inner lock is taken.
Scientific/helper/control sources and the worker are rehashed after
all calculations before a certificate can be written. A correct
root-scheduled command has the form

    python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock \
      --ledger results/NEW_EXCLUSIVE_LEDGER.json -- \
      scripts/check_golden_dyadic_sixteen_wall_cone_20261020.py --write

Replay uses a different new ledger and `--verify`. Actual execution
statuses remain separate from this source preparation. No frozen M6
source or scientific certificate is modified.

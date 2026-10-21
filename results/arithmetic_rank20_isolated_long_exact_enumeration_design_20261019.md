# Complete exact enumeration design for the bounded isolated long-contact rows

19 October 2026. **Source/paper design only. No enumeration has run.**
No worker, scientific process, optimizer, Lean or build was started.
This note supplies the coverage proof and exact interfaces for a future
small scheduled batch; it does not report that any row is empty.

The proved height inputs are
`arithmetic_rank20_isolated_long_dotted_height_bounds_20261019.md`,
SHA256 `9fc0583124cba17afed393bb569f83835e8f2c4be8ce5cfa64839f0369254516`.
Those written bounds still require independent review. This design is
conditional on accepting them, the actual-quad arithmetic interfaces,
and the integral-ring input stated explicitly below.

## 1. Exactly what would be enumerated

A batch is ONE field K, ONE N in{8,...,14}, and ONE rooted contact
type/endpoint. The spectator remains explicitly orthogonal to c,d.
It is not a general facet catalogue, arbitrary-height search, or a
completion of the ten half-refund consumers. The complete target types
and integral variables are:

|type|A|physical interval|integral scalar identity|
|---|---|---|---|
|endpoint3 only|2x|2<=A<24|C*delta=4(N+1)|
|next-to-end3 only|2x|2<=A<16|C*delta=8|
|endpoint4 only|2sqrt2*x|A>0,8<=A^2,A<16|C*delta=8|
|endpoint3 plus e3|2x|2<=A<18|C*delta=N+2|
|endpoint3 plus e4|2x|2<=A<4|C*delta=8|

Near and far roles remain distinct. The two r10 finite E8
extensions were excluded by the separate C-unit proof, rather than
silently dropped from finite classification. The next-to-end3 and
endpoint4 rows at N=8,12 also have an explicit odd-unit exclusion
in the height source. A future manifest must record these proved
omissions as named dependencies, not as rows that an enumeration
purportedly checked.

For the e4 type K is restricted to K16 or K24=Q(sqrt2,sqrt3),
both degree4. All other types initially retain all eight fields.

For clarity, the exact conjugate bounds used to construct the boxes
are listed here. In every row i>0 denotes a nonphysical embedding;
sqrt2_i has its actual sign at that embedding. Each strict interval
may be replaced by a certified closed OUTER rational interval.

|type|linear expression Z|necessary nonphysical inequality|
|---|---|---|
|endpoint3|NA+h, h=1 or N-1|Z_i^2<N+1|
|next-to-end3|NA+h, h=2 or N-2|Z_i^2<4|
|endpoint4|NA+2h, h=1 or N-1|Z_i^2<4|
|endpoint3 plus e3|2NA+h, h=N+2 or3N-2|Z_i^2<2(N+2)|
|endpoint3 plus e4|2NA+2h+N*sqrt2, h=1 or N-1|Z_i^2<4|

Thus, for example, endpoint3 has

    (-h-sqrt(N+1))/N < A_i < (-h+sqrt(N+1))/N.

The e4 intervals must be constructed separately for each embedding,
using certified intervals for that embedding's sqrt2_i. Replacing
all its conjugates by the positive sqrt2 would be unsound.

## 2. Exact field presentations and the external integer-ring input

Use ascending power coordinates in the following generator theta:

|K|theta|monic polynomial f(X)|
|---|---|---|
|Q(sqrt3)|sqrt3|X^2-3|
|Q(sqrt6)|sqrt6|X^2-6|
|K9|2cos(2pi/9)|X^3-3X+1|
|K11|2cos(2pi/11)|X^5+X^4-4X^3-3X^2+3X+1|
|K15|2cos(2pi/15)|X^4-X^3-4X^2+4X+1|
|K16|2cos(2pi/16)|X^4-4X^2+2|
|K20|2cos(2pi/20)|X^4-5X^2+5|
|K24|2cos(2pi/24)|X^4-4X^2+1|

**External integer-ring input:** O_K=Z[theta] for these presentations.
For the two quadratic fields this is the ordinary quadratic integral
basis. For the other six it is the standard maximal-real-cyclotomic
integer-ring theorem. Merely checking irreducibility of f does NOT
prove this integer-ring statement. A proof map must list it as an
external fact or separately provide a valid integral-basis certificate.
It is essential because enumeration of Z[theta] alone would otherwise
miss integral elements in a larger order.

The stated field presentations also use their declared minimal
polynomials and degrees. A later standalone checker must either
name that field input explicitly or certify irreducibility separately;
Sturm counts alone do not prove irreducibility.

The physical embedding is the LARGEST real root of f in each row.
All other real roots are retained, in a fixed certified order. Sturm
counts must show exactly deg(f) distinct real roots and one in each
disjoint rational isolating interval. No floating-point root ordering
can exclude or relabel a field embedding.

In the e4 type encode the ACTUAL sqrt2 element as

    K16: sqrt2=theta^2-2;
    K24: sqrt2=theta^3-3theta.

Reduction modulo f proves its square is2, and its positive physical
sign selects the intended square root. Its sign at every other
embedding is then determined by that embedding, not chosen freely.

## 3. Small exact field and sign primitives

Represent a field element by a length-d vector over Fraction, reduce
products modulo f, and compute inverses by rational polynomial Euclid.
Integrality is EXACTLY integer coordinates in the stated integral basis.
No approximate denominator test is allowed.

For a sign of b(theta_i), first reduce b modulo f. Exact zero is
recognized by the zero coordinate vector. Otherwise f is irreducible
and b has degree<d, so b is nonzero at every theta_i. Refine that
root's rational isolating interval until rational interval Horner
evaluation excludes zero. This terminates and decides the sign
exactly. Strict inequalities that happen to be equality are handled
by the exact-zero branch, with the correct strict/weak convention.

Roots of auxiliary rational polynomials, if used to enclose square
roots in a box construction, must also be isolated and counted
exactly. The final necessary conditions can instead be evaluated as
polynomial field inequalities, avoiding those auxiliary roots.

Existing sources read for possible reuse:

* `fresh_h20_global_fraction_verify_20261013.py`: small Fraction
  polynomial trim/arithmetic/division routines;
* `fresh_h20_triple_verify_20261013.py`: independently written
  positive-pseudo-remainder Sturm chains and root counts.

These are read-only, standard-library routines; their full scientific
drivers must NOT be invoked or imported for side effects. They have
not been pinned as dependencies of a new worker yet. The exact
quadratic model helper in
`check_arithmetic_dn_opposed_fork_models_20261019.py` is insufficient
for the six higher-degree fields; the large radical-tower helper is
also not a suitable replacement for this small quotient-field core.
A future implementation should use a small pinned pure helper or
an attributed copy with its own independent replay, without changing
the existing frozen sources.

## 4. Trace-dual boxes: a proof of complete integral-element coverage

Let beta_j=theta^j,0<=j<d. Compute the rational trace matrix

    S_jk=Tr(theta^(j+k))

using exact Newton recurrences from f, and its inverse over Q.
Define the trace-dual elements

    lambda_j=sum_k (S^(-1))_jk*theta^k.

Check Tr(beta_i*lambda_j)=1[i=j] exactly. For any integral
A=sum_j c_j*beta_j, c_j is the integer

    c_j=Tr(A*lambda_j).

Write A_i=sigma_i(A), with i=0 physical, and T=Tr(A) in Z.
Then A_0=T-sum_(i>0) A_i, so

    c_j=T*lambda_j(theta_0)
        +sum_(i>0) A_i*(lambda_j(theta_i)-lambda_j(theta_0)).     (1)

This is the useful correlated formula; one need not start from a
large symmetric coefficient box around A_0 alone.

For each contact type the height source supplies a finite physical
bound H and exact finite conjugate intervals for EVERY i>0.
Replace any irrational endpoint by a proved OUTER rational interval,
using certified root intervals. Closed outer intervals are harmless:
they may add candidates but never omit an actual element. These
give a finite closed rational interval for T, so enumerate every
integer from its ceiling to its floor. For each such T, rational
interval evaluation of(1) gives closed bounds L_j<=c_j<=U_j.
Enumerate ALL integers in each coefficient interval.

During the enumeration enforce the exact trace identity

    sum_j c_j*Tr(theta^j)=T.

Since Tr(1)=d, c_0 can be recovered from c_1,...,c_(d-1)
and rejected unless the required quotient is integral. This avoids
one enumeration dimension without an unproved lattice reduction.

Every actual integral A is in one of these trace boxes: the trace
identity and(1) are exact, and every interval used was an outer
enclosure of its actual embedding value. Thus the finite coefficient
boxes are an exhaustive consequence of the theorem's hypotheses.

Root intervals may be refined to reduce overestimation BEFORE
enumeration. They must never be narrowed without their root counts.
If the derived boxes exceed a resource budget, the status is
DEFERRED or an explicitly complete trace-range partition is prepared.
An arbitrary coefficient cut is never a replacement for the bounds.

## 5. Reconstructing every denominator and checking candidates

For each integer coefficient vector reconstruct A exactly and first
test the original physical range and all nonphysical proper-block
inequalities. The useful variables and denominators are:

* Endpoint3: h=1 or N-1, Z=NA+h,
  C=[2Z^2-(N+2)(N+1)]/N, delta=4(N+1)/C.
* Next-to-end3: h=2 or N-2, Z=NA+h,
  C=[Z^2-2N-4]/N, delta=8/C.
* Endpoint4: h=1 or N-1, Z=NA+2h,
  C=[Z^2-2N-4]/N, delta=8/C.
* Endpoint3 plus e3: h=N+2 or3N-2, Z=2NA+h,
  C=[Z^2-(N+2)^2]/(4N), delta=(N+2)/C.
* Endpoint3 plus e4: h=1 or N-1,
  Z=2NA+2h+N*sqrt2,
  C=[Z^2-2N-4]/N, delta=8/C.

Each C must first be checked integral by exact coordinates; the
displayed algebra predicts this and an implementation failure must
not be treated as candidate exclusion. C=0 is incompatible with
the nonzero right side and rejects that vector. Check delta has
integer coordinates, put z=delta+4, and test

    delta_0>0,       0<z_i<8/r for every i>0.

Recheck the original corank identity exactly, and the original qc
cost physically>1 and nonphysically<1. These checks are cheap
field arithmetic and exact signs, not a numerical tolerance.

In particular the original quantities can be evaluated using these
polynomial expressions, with r=N+2:

|type|R=(x+kappa)^2/t|qc|
|---|---|---|
|endpoint3|Z^2/[2N(N+1)]|1-(N+1)/(2N)+Z^2/(2N)|
|next-to-end3|Z^2/(8N)|1-2/N+Z^2/(2N)|
|endpoint4|Z^2/(8N)|1-1/N+Z^2/(4N)|
|endpoint3 plus e3|Z^2/[4N(N+2)]|1-(N+2)/(4N)+Z^2/(8N)|
|endpoint3 plus e4|Z^2/(8N)|1-1/(2N)+Z^2/(8N)|

Every row must satisfy (R-r/(4N))*delta=1. The endpoint4 formulas
are deliberately expressed without separately adjoining sqrt2: the
integral variable A=2sqrt2*x and its squared products lie in K even
when sqrt2 itself does not. This differs from the e4 type, where
the actual triangle cycle forces sqrt2 into K.

For endpoint4, additionally check w=4x^2=A^2/2 is integral.
Integral A=2sqrt2*x does NOT imply this automatically. For all
other types w=A^2 is automatically integral. In the e4 type the
additional triangle cycle sqrt2*A is integral as well and can be
checked directly. All native finite long/e contacts are fixed by
the type, not inferred from a rounded Gram matrix.

One may check that z generates K by the rational rank of
1,z,...,z^(d-1). This is necessary because the physically negative
old-plus-d determinant forces field equality. Uniformly testing
z is safe. Do NOT require A alone to generate K in the e4 type:
its original cost also contains sqrt2, so a proper-subfield A may
still generate the full K together with sqrt2. Omitting this
distinction would lose legitimate candidates. The absence proof
could also retain such extra candidates; it must never discard
them using a false A-primitivity premise.

No passing vector is automatically an actual global polytope or
even an actual face. It is an unresolved necessary-condition hit.
Every such hit must be saved and flagged; a nonempty hit set
blocks a claimed exhaustive absence result.

## 6. Compact certificate and genuinely separate replay

A batch certificate should bind field/type/N, the height theorem,
all source hashes, the integer-ring input status, f and its physical
root designation, rational root intervals and Sturm counts, exact
trace/dual matrices, trace ranges, coefficient boxes, enumeration
counts and rejection totals, and every surviving exact coordinate
vector. Include a deterministic stream digest of visited coefficient
vectors and outcome codes; do not save a giant per-rejection log.

The checker must RECONSTRUCT the coefficient boxes and enumerate
their complete integer trace slices. Rechecking only the producer's
saved survivors or claimed counts does not prove absence. Prefer
separate enumeration/control flow and sign/polynomial checks. A
shared pinned pure arithmetic helper must be named explicitly as
a common dependency, not described as independent implementation.
Any unhandled sign, field, role, source mismatch or output overflow
must fail closed and retain a compact unknown-hit record.

Use distinct statuses for design prepared, execution incomplete,
complete with unresolved hits, and complete with zero necessary
hits. Only the last can discharge THAT field/N/type arithmetic
branch, conditional on the geometric inventory and ring input.

## 7. Resource and promotion boundary

No runtime promise is justified before box counts are prepared.
The intended unit is one scheduled field/N/type batch, serial under
the shared lock, at most59seconds and1GiB RSS, with bounded output
and the repository's disk/memory gates. No worker runs in this
design task. A timed-out or resource-deferred batch is incomplete;
it cannot be replaced by an unscheduled concurrent attempt.

After any future complete exact replay, the whole-factor finite
inventory, spectator-contact rows, triangle consumers and
simultaneous half-source charge still require explicit composition.
This design alone neither changes the frontier cost nor lowers the
compact arithmetic dimension ceiling.

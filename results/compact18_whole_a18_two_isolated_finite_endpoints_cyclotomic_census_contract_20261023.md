# Whole A18: isolated finite endpoint pairs and a complete cyclotomic domain

23 October2026. **SOURCE-ONLY PREPARATION, UNREVIEWED and UNEXECUTED.**
There are zero actual batch or combine outputs. The root proposed the
geometric/coupled-cosine reduction; I independently rederived it and supplied
the sharper n<=105 domain and the checker below. No worker import, mathematical
execution, enumeration, test, compiler, Lean, TeX, plan or Git operation ran.

Companion standalone stdlib source:
`scripts/check_compact18_a18_isolated_finite_endpoint_pairs_20261023.py`.
It has sixteen units (m=7,...,14; mutual k=3 or4) and a SEPARATE authenticated
coverage combine. Only root may execute them, serially through the unchanged
59s/1GiB/one-thread research guard, with its parent-held exclusive lock and
40GiB-disk/6GiB-RAM admission. No successful or zero-survivor outcome is assumed.

## 1. Exact local theorem interface, with no arithmetic restriction

Fix an ACTUAL whole-A18 farthest vertex p of a compact hyperbolic Coxeter
18-polytope. Assume BOTH actual opposite facets of its two endpoint edges
are isolated WITHIN the old A18 and have FINITE own contact. This isolation
and finiteness are hypotheses, not proved for all whole-A18 configurations.
No arithmeticity, ground-field bound, arbitrary face inheritance, component
floor, different farthest vector or disconnected root scale is used.

The only inherited mathematical source pin is

| Source | SHA256 | Restricted interface |
|---|---|---|
|results/compact_farthest_vertex_mixed_replacement_inequalities_20261020.md|3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84|Actual compact acute finite vertices and Lorentz cap conventions; no numerical rank transfer|

Standard external structural facts are the acute-wall Gram rule, finite
Coxeter native angles, and the positive-null affine exclusion for compact
chambers. The latter follows directly: a positive null combination of actual
facet normals is future by its strict interior sign, orthogonal to those
facets and nonpositive against EVERY other actual facet. Its ideal ray would
belong to the closed chamber, contradicting compactness.

Use p^2=-1 and UNIT old normals n1,...,n18 with diagonal1 and adjacent
product-1/2. The ACTUAL endpoint caps are unit normals

    cL=AL*p-gL*Omega1, cR=AR*p-gR*Omega18,
    Omega1^2=Omega18^2=36/19, Omega1.Omega18=2/19,
    AL^2=(36/19)*gL^2-1, AR^2=(36/19)*gR^2-1.       (1)

Their positive time coefficients follow from the old actual vertex and
the fixed interior/future convention. These two caps are distinct: each
has positive own contact and zero contact at the other old endpoint.

Own0/3/4 have nonpositive time square (strictly negative in these cases).
Native6 is impossible: on actual c,n1,n2 the native6/3 chain has positive
null combination sqrt3*c+2n1+n2. The reflected endpoint is identical.
Every other finite own label is5 or at least7.

## 2. Complete pair reduction: native5 disappears; mutual3/4 only

Write D=36/19, e=2/19, and g5=cos(pi/5). The actual mutual magnitude is

    M(g,h)=sqrt(D*g^2-1)*sqrt(D*h^2-1)-e*g*h.       (2)

It is strictly increasing in each variable throughout g,h>=g5. Indeed
g/sqrt(D*g^2-1)>1/sqrt(D), and
(D^2-e^2)h^2>D holds because D^2-e^2=68/19 and h^2>=g5^2>9/17.
This proves the positive derivative, or the equivalent inequality after
squaring positive quantities, without a numerical monotonicity guess.

If either own label is5 and the other is finite, then

    M>=M(g5,g5)=(17sqrt5-25)/76>0,
    M<M(g5,1)<sqrt(A5^2)<1/2,
    A5^2=(9sqrt5-11)/38<1/4.

Thus the two DISTINCT actual walls would have positive mutual magnitude
strictly below1/2, which is not a native finite value and is neither
parallel nor dotted. Native5 is excluded in either orientation.

Consequently g=cos(pi/m),h=cos(pi/n), m,n>=7. Relabel the old A18 path so
m<=n; this is a coordinate relabelling, not a chamber symmetry. Equation(2)
is positive, since already M(cos(pi/7),cos(pi/7))>0; for example cos(pi/7)>
7/8 by cos(x)>1-x^2/2 and pi<22/7. Its supremum at g=h=1 is15/19,
which is below cos(pi/5). The native mutual label therefore is EXACTLY
k=3 or4. Orthogonal, parallel and dotted mutual contacts are impossible.

If m,n>=15, monotonicity gives

    M>=(34/19)*cos^2(pi/15)-1>1/sqrt2.

For a rational seal, cos(pi/15)>10783/11025 follows from pi<22/7 and
cos(x)>1-x^2/2. Also sqrt2<99/70, and

    ((34/19)*(10783/11025)^2-1)-99/140
       =1498515515/323324662500>0.

This contradicts the largest allowed mutual value cos(pi/4). Hence every
actual pair has one smaller label m in7,...,14 and n>=m.

## 3. The exact polynomial equation and complete degree cutoff

Put a=2cos(pi/m), b=2cos(pi/n), c=2cos(pi/k). Then c=1 or sqrt2. The
positive-time equation(2) is equivalent to

    2sqrt((9a^2-19)*(9b^2-19))=a*b+19c.

Both sides are positive for m,n>=7. Squaring and dividing by19 gives

    E=17a^2*b^2-36(a^2+b^2)-2c*a*b+76-19c^2=0.    (3)

Thus no wrong sign is introduced by the squaring in this domain. Viewed
as a polynomial in b, its leading coefficient17a^2-36 is strictly positive.
Let d=[Q(a):Q]=phi(2m)/2, dc=[Q(c):Q]=1 or2. For m7,...,14 the complete
d sequence is3,4,3,4,5,4,6,6. Equation(3) implies

    [Q(b):Q]<=2*d*dc, hence phi(2n)<=4*d*dc<=48.     (4)

This uses a containing field Q(a,c); no equality or linear disjointness
of its degree is assumed. The prime-power formula gives

    phi(q)^2/q=product_(p^a||q) p^(a-2)*(p-1)^2>=1/2:

only the factor p=2,a=1 can be below1, where it equals1/2; all odd-prime
factors are at least4/3 and higher exponents only increase them. Thus
phi(2n)^2>=n and the simple universal bound is n<=2304.

There is a useful COMPLETE sharpening: if q is EVEN and phi(q)<=48,
then q<=210. At least five distinct primes have minimum totient480.
Four distinct primes force exactly2,3,5,7, all to exponent1, giving210.
For three primes write q=2*p*r*M with odd p<r and
M=2^(a-1)*p^(b-1)*r^(c-1). Then phi(q)=(p-1)(r-1)M.
Thus ALL higher exponents are retained through the integer M, and
q<=2*p*r*floor(48/((p-1)*(r-1))). No exponent is silently set to1.
The COMPLETE (p,r) possibilities and upper bounds on q are

| p | r | q upper bound |
|---|---|---:|
|3|5,7,11,13,17,19,23|180,168,132,156,102,114,138 respectively|
|5|7,11,13|140,110,130 respectively|

Here p>=7 is impossible because6*10>48. For at most two distinct primes,
q=2^a*p^b gives q/phi(q)=2p/(p-1)<=3, so q<=144; a power of2 is even
smaller. This proves q<=210, NOT an executed list. Hence the final domain is

    m=7,...,14; k=3,4; m<=n<=105;
    phi(2n)<=4*(phi(2m)/2)*(1 if k3 else2).         (5)

No n removed by(5) can be actual by(4). This is the entire unordered,
finite-isolated-endpoint pair domain, conditional only on Section1.

## 4. Exact cyclotomic algorithm; no floating approximation

For one pair let L=lcm(2m,2n,2k), zeta=exp(2pi*i/L). Substitute

    a=zeta^(L/(2m))+zeta^(-L/(2m)),
    b=zeta^(L/(2n))+zeta^(-L/(2n)),
    c=zeta^(L/(2k))+zeta^(-L/(2k)).

Construct E as a sparse integer polynomial modulo X^L-1 FIRST, so every
exponent lies in0,...,L-1. The standard minimal-polynomial fact for a
primitive Lth root is external elementary algebra: E(zeta)=0 iff the
remainder of E on division by Phi_L is zero. This tests the exact physical
positive cosines too; it does not select a different conjugate or a
rational numerical tolerance. All equality/survivor branches are retained.

Phi_L is constructed by the defining exact recursion

    Phi_d=(X^d-1)/product_(e|d,e<d) Phi_e,

for divisors d of L in increasing order. Each monic integer division must
have zero remainder. This is justified inductively by X^d-1's partition
into primitive roots. Each new Phi has degree phi(d), leading coefficient1
and nonzero constant coefficient. No CAS, field table, cached external
polynomial or approximate complex arithmetic is used.

For(5), L=2*lcm(m,n,k)<=11760. Also
phi(lcm(x,y))*phi(gcd(x,y))=phi(x)*phi(y) implies
phi(L)<=phi(2m)*phi(2n)*phi(2k)<=2304. For ONE pair, the cached divisor
polynomials have total degree sum_(d|L)phi(d)=L, then the cache is released.
There is no cross-pair polynomial cache. Dense arrays have length at most
L+1, integer operations terminate, and no enumeration of roots is required.

The integer bit guard is4*11760+128. Cyclotomic coefficients are bounded
by2^degree, since their roots have absolute value1. Exact-division
intermediates are bounded by the full dividend plus partial quotient times
divisor. For the final E division, coefficients of 1/Phi_L through degreeL
are bounded by binomial(phi(L)+L,L)<=2^(phi(L)+L); multiplying the at mostL
input coefficients (absolute value below1024) and the divisor gives a
bound below2^(4L+128). Thus the guard is conservative; crossing it aborts
with NO PASS. The source has degree/order/file/line bounds as well. A
runtime resource PASS is still required; this is not a measured RSS claim.

## 5. Source, output, guard and combine contract

The worker imports stdlib only and never imports repository workers. It
requires its actual --expected-source-sha256 and checks itself plus all
fixed pins before and after each successful unit. In addition to this
contract and the Section1 source, it binds the unchanged bounded runner
and its monitor/signal-safe/base-runner dependencies. Their current hashes
are explicit in the source. Root's launch MUST use --parent-lock, a fresh
ledger, fresh output and the unchanged resource gate. No command is run here.

For each(m,k), the batch emits one header, every n in(5) in increasing order,
and one final summary. Each exact pair record contains labels/degrees/L,
remainder hash/degree/nonzero count, and its first nonzero coefficient as
signed hexadecimal with index. A zero remainder is explicitly SURVIVOR;
survivors are NEVER asserted absent. No huge coefficient array is saved.
The complete output is create-only, bounded and streamed. After it closes,
the worker prints one small stdout seal with its path/bytes/SHA256 and
terminal summary, binding the file to the guard's recorded actual stdout.

The separate combine takes exactly sixteen part paths, output hashes,
ledger paths and ledger hashes. It checks every current pin and file hash,
successful exact-command59s/1GiB/one-thread/parent-lock ledgers, admitted
initial resources, stdout seal, strict duplicate-free JSON, exact ordered
label domain and full counters/summary. No missing, repeated, extra or
stale unit is credited. Files are rehashed after their serial read. It
copies every survivor and emits per-part bindings plus the final counts;
it does NOT repeat the cyclotomic mathematics of sixteen actual successful
batches. Combine PASS means authenticated COMPLETE COVERAGE, not a new
sixteen-batch mathematical execution or an assumed empty result.

Prepared arithmetic fixtures exercise totients and the complete cutoff,
several exact cyclotomic polynomials, exact/nonexact monic division,
cyclic convolution and both zero/nonzero remainder paths. They have NOT
run. Their eventual execution is part of a guarded unit, not current
verification evidence.

Batch CLI shape (root only, UNDER the guard):

    check_compact18_a18_isolated_finite_endpoint_pairs_20261023.py
      --expected-source-sha256 CURRENT_HASH batch --m 7 --k 3
      --output results/NEW_PART.jsonl

For combine repeat --part/--part-sha256/--ledger/--ledger-sha256 sixteen
times and supply a fresh --output. The outside source review and root
preclear must bind actual final code/contract hashes before any launch.
No global whole-A18 exclusion, bound17, example or acceptance follows from
this preparation. Even an eventual empty(5) would exclude only the BOTH
isolated-and-finite endpoint branch until separate geometry covers others.

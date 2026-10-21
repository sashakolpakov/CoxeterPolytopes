# Norm239: a small folded domain for the released-m edge

22 October2026. **Prepared mathematical contract, pending full source and
outside-author review; no checker has been imported or executed.** This
replaces neither frozen input below. It folds a large finite upper into
four short phase intervals and specifies an exact, constant-storage
checker. No first-wall result, full new germ or chamber period is assumed
or reported. All scientific/Lean/TeX launches remain root-owned.

## 1. Inputs and precise target

Use the literal ORIGINAL module L=E8(O) perpendicular-sum
`[[-2phi,1],[1,8]]`, O=Z[phi], phi²=phi+1. Put D=1+16phi,
b=phi D, A=phi², E=7-4phi and F=3-phi=A E. All vector notation,
the original covectors and the retained native E7 are as in these inputs:

| Input | SHA-256 | Boundary |
|---|---|---|
| results/compact_h9_prime239_first_e7_golden_transverse_edge_candidate_20261021.md | 2a207fb82414203e34625d323238abc4c297541b6e34f23c4d28b16b7c3efa15 | Actual full starting germ Q, E7, m, R; pending separate review |
| results/compact_h9_prime239_second_h2_edge_finite_upper_contract_candidate_20261021.md | 1d11ee290abb97a9da0aaab1d86496071de46d2cf74a02489457fb7f369f89b4 | Original-integral finite upper V³m; exact V, its quarter-integral powers and original parity |
| results/compact_h9_prime239_full_base_and_first_binary_edge_candidate_20261021.md | 879a3feba630072a4b127e68a3729c4da3fca4f3c44630d9aae623ddd899de02 | Original binary module/glue, actual m and base wall r0 |
| results/golden_odd_glue_full_base_and_discriminant_reduction_20261021.md | b107f0de4f76da942100e0c82867cd6d4eeee40d833304fee67cebd84d3e41ce | Prime norm239, primitive norm alphabet2/2b, local finiteness and actual base |

The outgoing line is P(z)=Q-zW, W=m+(phi/2)R, with
Q²=-2DE, W²=F/2 and Q.W=0. Here E7 and R stay zero,
whereas P(z).m=-zF/2<0. Its ideal parameter is2sqrt(D)/phi.
The previous source proves that V=SBS fixes E7 and R, is future,
has ORIGINAL entries in(1/4)O, and translates this Lorentz line in
the outgoing direction by a distance l0>0 with

    tr_active V=t=563+914phi, cosh(l0)=t/2.

V^+-2 also has quarter-integral entries by
`V²=tV-I+(2-t)Pi`, and V^+-3 is original-integral. The finite
upper is V³m at the future point V³Q, at distance3l0. None of
V, V², V³ is assumed to preserve the chamber. In particular, an
original-integral image is NOT automatically an actual simple wall.

## 2. Why the small H2 orbit does not itself continue the edge

Let n=-e7. Its products with m,R are(0,-phi³), and Q.n=-D.
Reflections in m,R fix Q. The complete five-element orbit has

| (r.m,r.R) | r.W |
|---|---|
| (0,-phi³) | -phi⁴/2 |
| (-phi⁴,phi³) | -phi⁴/2 |
| (phi⁴,-phi⁴) | phi²/2 |
| (-phi³,phi⁴) | phi²/2 |
| (phi³,0) | phi³ |

These follow successively from `(a,b)->(a+phi b,-b)` and
`(a,b)->(-a,b+phi a)`. Every row still has Q.r=-D.
The only positive collisions are z=2D/phi⁴; they exceed the ideal
parameter because D-phi⁶=-4+8phi>0. All other collisions are
negative. Thus this finite orbit supplies no finite upper or next
facet. It proves nothing about other integral mirrors. The nonintegral
H2 diagram swap in the preceding source remains unavailable.

## 3. Four phases and the complete first-simple-wall profiles

For ANY first new facet at signed line distance l in(0,3l0], choose
a nearest integer j to l/l0. Then j is0,1,2 or3 and
`|l-jl0|<=l0/2`. Put s=V^-j r. It is in(1/4)L, and after
reconstruction r=V^j s it must be an ORIGINAL-integral reflector.
The endpoint phases0 and3 must both be included. Overlap at half
steps is intentional; canonical original vectors remove duplicates.

V fixes E7 and R, so every first-simple-wall profile remains

    s=(Tq-Km)/2+hH+epsilon J, s.R=-kappa,
    norm2: (epsilon,kappa)=(0,0),(0,1),(0,phi),(1,0),
    norm2b: (epsilon,kappa)=(0,0),                     (1)

where H=Omega7, J=Omega6 and T,K,h lie in(1/4)O. Folded h
may have EITHER SIGN. Acuteness to old m,-e7 is imposed only on
the reconstructed original root. These five profiles exhaust spherical
one-node extensions of E7+R: E7+I2(k), or E8+A1 at node6.
The exact O strip1<=kappa<2, |kappa'|<2 gives1 orphi. A norm2b
nonzero finite contact would have angular invariant b eta² of norm
at least239, but both embeddings below4, and is impossible.
These are profiles of a first SIMPLE facet, not all arbitrary mirrors.

The contact and norm equations are exact even on the folded lattice:

    Lh=2A(2h+3epsilon)+2kappa/phi,
    K=DT-Lh,
    M0=4h²+12epsilon h+12epsilon²-4d,
    D(D-1)T²-2D Lh T+Lh²+M0=0,                     (2)

d=1 or b. Thus for each h there are at most two field possibilities

    Z=D[Lh²-(D-1)M0],
    T=(D Lh +/- sqrt(Z))/(D(D-1)).                   (3)

Both signs, a zero square and the rational/pure-sqrt5 subcases must
be retained. Exact Q(phi) square membership is required, not only
a rational-norm-square test. The implementation uses Q(sqrt5): for
Z=a+c sqrt5 with c!=0, test both signs of the rational square root
of a²-5c², then rational-square `(a+q)/2`, setting y=c/(2x).
Verify `(x+y sqrt5)²=Z`. For c=0 test the rational and pure-sqrt5
possibilities separately. This is complete because2xy=c.

The seed ORIGINAL coordinates must satisfy4s in L. Merely placing
T and h in quarter-O does not replace this test. Reconstruction then
tests r in L and every original reflection matrix entry, with exact
norm2 or2b. Norm classes are primitive: a nonunit square cannot
divide2 or2b, since the primes2 and b occur only to exponent1.

## 4. A small proved scalar domain

The centered window has C=cosh(l0/2)=sqrt((t+2)/4)<23, using
phi<81/50 and t<2044. Let U=Q/sqrt(2DE). For a norm2d seed
whose mirror meets this window, Lorentz Cauchy gives

    |s.e7|<=sqrt(2d)[sqrt(C²-1)D/sqrt(2DE)
                           +C sqrt(2+D/(2E))].       (4)

The argument works for both signs of the collision distance and both
signs of h=s.e7. At the conjugate positive form, |h'|<=2sqrt(d').
Since D/(2E)<673/26<26, the two square roots in(4) are respectively
<51/10 and<53/10. Therefore

    d1: |h|<342, |h'|<=2.                            (5)

For the norm2b row we also have h in(b/4)O. Here is a divisibility
proof that does not pretend V integral. Write f=4s in L. Its norm
is32b, so some coordinate f_i is not divisible by the prime b.
The seed reflection is `V^-j s_r V^j`, with matrix in(1/16)O.
Since it also equals `I-f(Gf)^t/(16b)`, every f_i(Gf)_k is
divisible by b. Primality forces every covector(Gf)_k divisible
by b. Its E8 covector at7 is4h, hence b divides4h. This uses
the original reflection of r and both quarter-integral inverse powers;
it is not a primitive-root assertion over an unrelated quarter-module.

Put jh=h/b. Formula(4) now gives

    |jh|<53, |jh'|<43/50.                            (6)

Indeed b>216/5 yields sqrt(2/b)<11/50, and b'>273/50 gives
2/sqrt(b')<43/50. No original norm2b solution is lost by(6).

For d1 write4h=a+c phi. From(5) and sqrt5>11/5,
`|c|<=625`, and the exact conjugate interval
`|a-c(phi-1)|<=8` contains at most17 integers a. This gives at
most1251*17=21267 scalar rows per norm2 profile/phase. For db
write4jh=a+c phi. Equation(6) gives |c|<=97 and
`|a-c(phi-1)|<86/25`, at most7 integers a. There are at most
195*7=1365 scalar rows. All physical and conjugate inequalities
are checked exactly, further shrinking these boxes. The total over
four phases and all five profiles is at most

    4*(4*21267+1365)=345732 scalar rows,              (7)

each with at most two roots from(3). This is a mathematical loop
bound, not a measured runtime or a claimed checker result.

## 5. Exact filtering, firstness and full-germ recognition

For a seed, the signed folded collision is

    zseed=D(2AT-2h-3epsilon)/(K+phi kappa/2).

The denominator must be nonzero. Require

    zseed² <= [4D/A]*(t-2)/(t+2),                    (8)

which is exactly the half-step window, including zero and both signs.
After reconstruction r=V^j s, require ORIGINAL coordinates integral,
norm2d, and every coefficient of `r(Gr)^t/d` integral. Require
Q.r<0, -W.r>0 and0<z=-Q.r/(-W.r)<=zup, with the original
upper from the preceding source. Check acuteness to all explicitly
known actual walls(-e0..7,r0,m,R). All equality comparisons are
exact field comparisons and all minimum ties are retained.

Every actual first simple wall is in one of these20 phase/profile
domains by the folding proof. Conversely every admitted original
reflector is an arrangement mirror. No such positive collision can
occur strictly before the actual first collision while the outgoing
segment stays in C. Thus the least collision across all20 complete
parts equals the actual first collision. Empty individual parts are
allowed; an empty aggregate would contradict the proved finite upper.

For each tied root r, form the nine vectors(-e0..6,R,r). Check a
positive definite Gram by exact elimination and reconstruct every
other tied oriented root as a NONNEGATIVE combination of these nine
vectors. Exact vector equality as well as coefficient signs is checked.
The unique such new ray identifies the actual simple cone: the actual
first facet is in the complete tie list and its cone contains all
oriented local roots. If this cone test is unresolved, the checker
reports it explicitly and does not certify a full germ. No row-order
selection or hardcoded model winner is permitted.

At a finite point the local reflection group is finite. Its rank is at
most9; canonical primitive classes2/2b admit finite mutual labels
only2,3,5. Finite spherical classification bounds its positive roots
by121 (E8+A1 or H4+H4+A1 attains the numerical maximum). Hence
128 minimum ties suffice; exceeding that bound is an ERROR, never
truncation. The worker stores no raw row arrays. It keeps counters,
the exact minimum and its at most128 ties. Each output is new-only.

## 6. Prepared execution interface and limits

Worker: `scripts/check_compact_h9_prime239_second_h2_edge_batches_20261022.py`.
It has no scientific imports and does not import another worker.
Twenty units are phase0..3 times profile orthogonal/simple/golden/e8/
longnorm. Each runs separately under the unchanged59s/1GiB serial
runner and fresh40GiB/6GiB admission; the unchanged8GiB/3GiB abort
floors and live watchdog remain required. No runtime bound is inferred
from(7); a timeout earns no domain PASS. Example queued command:

    python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock \
      --ledger results/prime239_folded_phase0_orthogonal_execution_NEW.json -- \
      scripts/check_compact_h9_prime239_second_h2_edge_batches_20261022.py \
      --phase 0 --profile orthogonal \
      --output results/prime239_folded_phase0_orthogonal_NEW.json

The separate `--combine --part PATH` interface requires exactly all20
phase/profile outputs under the SAME current worker/contract/input pins,
rechecks recorded roots and performs the cone test. It does not substitute
for missing part executions. Root must bind actual runner ledgers before
claiming a complete runtime PASS. Inputs and worker are SHA-checked
before and after each unit/combine; existing destinations are refused.

All candidate sources here remain supporting research. The dyadic
original module, actual starting germ and first-simple-wall boundary
are explicit. This contract supplies neither an unobserved census
outcome nor T(C)=C, finite volume, nonreflectivity, a compact H9
example, a form-class exclusion or a dimension bound.

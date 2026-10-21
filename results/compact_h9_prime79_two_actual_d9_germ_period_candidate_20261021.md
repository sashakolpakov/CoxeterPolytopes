# Norm79: an original-module period between two full D9 germs

21 October 2026. **Manual supporting candidate, pending full outside-author
and root review.** No worker, mathematical import, enumeration, test, Lean
or TeX process was run in preparing this proof. The END0 execution below
is inherited evidence, not a new execution. No paper, guide, checkpoint,
global bound or example was changed. The proposed literal-host conclusion
requires acceptance of the explicitly pinned actual-germ boundary.

I derived the original-coordinate glue proof below. The root supplied the
two equal-norm corner normalizations and the short finite-order obstruction
in Section4; those shared contributions are not independent reviews.

## 1. Exact actual-germ premises

Use the ORIGINAL module

    O=Z[phi], phi²=phi+1, a=3+3phi, b=13+18phi,
    D=18+31phi, L=E8(O) perpendicular-sum [[-2phi,a],[a,2]].

The ordered native basis is `(Bx,By,e0,...,e7)`. Every E8 root has
square2; its edges are02,13,23,34,45,56,67. Put `c=-By`,
`p1=2Bx-a By`, `Q=omega0=(4,5,7,10,8,6,4,2)` and
`W=omega7=(2,3,4,6,5,4,3,2)` in E8 coordinates.

| Input | SHA-256 |
|---|---|
| `results/compact_h9_nonunit_coupling_prime79_full_base_and_first_transverse_20261021.md` | `33f6cb9dd20ba91e7e7974586075309f2f437564f2e8dd978a16594349e48da3` |
| `results/compact_h9_prime79_actual_a9_end_edge_upper_and_complete_contract_candidate_20261021.md` | `615d2915b1f95ffe5dc4d0ec634840b67991a7afb7b620aaca63d498bb10c2d3` |
| `results/compact_h9_prime79_actual_a9_end_interval_checker_spec_20261021.md` | `b21fe0b47693d0767b9b27283bdc7ba78e5e141857577b7ebd153aad7df94d4d` |
| `results/prime79_actual_a9_end0_interval_20261021.jsonl` | `7f35075ed06766a4d78bc328ac9f8b8a819e9a860ec8a32b92d13d2b123319f0` |
| `results/prime79_actual_a9_end0_interval_execution_20261021.json` | `ff8d43f1dd228712db84147e2ce23e8813d286a23f4c4dbc00fd39ca016146a5` |

Section5 of33f6 proves an OLD full actual D9 germ, ordered as the
long-end-to-center chain

    (c,n,-e7,-e6,-e5,-e4,-e3), leaves (-e1,-e2),
    n=(phi/2)p1-c/2+W+phi³Q.

The other two sources prove the complete six-profile END interval and
original-module reconstruction/filtering, including closed first-collision
ordering, every tie and the positive-simple-cone germ test. The terminal
END0 receipt records all810436 raw scalar trials, six integral mirrors
and the unique minimum id4. Its NEW full D9 proposal is

    (-e2,-e3,-e4,-e5,-e6,-e7,s), leaves (c,h),

where s is the observed actual A9 bridge and h is id4. Their ORIGINAL
coordinate vectors are

    s=(16+25phi,-61-99phi;
       162+260phi,259+416phi,324+520phi,486+780phi,
       389+624phi,292+468phi,195+312phi,98+156phi),
    h=(18+30phi,-72-117phi;
       189+309phi,301+493phi,377+617phi,565+925phi,
       452+740phi,339+555phi,226+370phi,113+185phi).

The recorded minimum has controls `A0=1,E0=11+9phi` and collision
`t=(810+1790phi)/139`. Acceptance of the new germ means the whole
firstness/coverage/geometric implication of these sources, not merely
the receipt's PASS label. No period argument below substitutes for it.

## 2. Equal negative vectors and a unique K-isometry

The nine OLD wall roots are orthogonal to the future corner vector

    u=(3+4phi)p1+(Dphi/2)Q, u²=-D.

Indeed `2(3+4phi)²-Dphi²=1`, `p1²=-2D`, `Q²=4` and
`Q.n=2(3+4phi)`. Thus u is future and orthogonal to all nine old
roots, including n. It need not belong to L; this is a rational
orthogonal-complement calculation, not a change of integral module.

Let P and Omega be the NEW interval's recorded corner and released
A9-end dual, so `P²=-10D`, `Omega²=9/10`, `P.Omega=0`.
Set

    E=11+9phi, E²-9D=40,
    v=(E P-10D Omega)/20.

The id4 formula `h=w_B+P/10-(E/9)Omega` shows `v.h=0`;
both P and Omega are orthogonal to the retained eight roots.
Consequently v is orthogonal to the NEW D9, has square
`-D(E²-9D)/40=-D`, and is future since `v.P=-DE/2<0`.

Define T on the seven ordered chain roots by their displayed order,
send `-e1` to c and `-e2` to h, and send u to v. The two positive
D9 Grams agree. Together with u,v these are nondegenerate ten-vector
bases with the same Gram; hence T is a unique K-linear isometry.
It is future preserving because u and v are future. Explicitly,

    TBy=e2, Te1=By, Te2=-h, Te3=-s,
    Te4=e7, Te5=e6, Te6=e5, Te7=e4, Tn=-e3.       (1)

This does not infer integral-module preservation from matching diagrams.
The missing images TBx and Te0 are checked next.

## 3. Original glue: an explicit mod4 proof, not an index-D shortcut

Write `A=6+8phi` and `B=62+98phi`. The Bx/e0 coefficients
of n are `(phi,A)`, and those of u are `(A,B)`. Their determinant is

    phi B-A²=-2.                                      (2)

Let subscripts denote original coordinates and define integral vectors

    N=-e3-n_y TBy-sum(j=1..7,n_j Tej),
    R=2v-2u_y TBy-sum(j=1..7,2u_j Tej).

All2u_j and2v_j are integral. The equations for Tn and Tu give

    TBx=(A R-2B N)/4,  Te0=(2A N-phi R)/4.            (3)

Here2A and2B are divisible by4. It therefore suffices to show
`R in4L`. For reproducibility, the exact vector2v is

    (1106+1786phi,-4338-7017phi;
     11412+18469phi,18234+29507phi,22806+36907phi,
     34200+55345phi,27360+44276phi,20520+33207phi,
     13680+22138phi,6840+11069phi).

Substituting (1) in R, and reducing in `O/4O`, gives

    R = 2v+(2+2phi)e2-(3+phi)By
           +(1+3phi)h+(2+2phi)s
           -(2+2phi)e6-(2+2phi)e4  (mod4L).          (4)

The complete coordinate residue table is below. The middle columns are
the two displayed scalar multiples, not unnormalized wall products.

| coordinate | 2v | (1+3phi)h | (2+2phi)s | R |
|---|---|---|---|---|
| Bx |2+2phi|2phi|2|0|
| By |2+3phi|1|2phi|0|
| e0 |phi|3phi|0|0|
| e1 |2+3phi|3phi|2+2phi|0|
| e2 |2+3phi|3phi|0|0|
| e3 |phi|3phi|0|0|
| e4 |0|0|2+2phi|0|
| e5 |3phi|phi|0|0|
| e6 |2phi|2phi|2+2phi|0|
| e7 |phi|3phi|0|0|

The By,e2,e4,e6 corrections in (4) complete the last column. Thus
(3) is integral on the original basis. Since T is an isometry of a
nondegenerate K-form, `det(T)²=1`, so its determinant is a unit and
its inverse is integral by the adjugate formula. **T(L)=L**, not just
preservation of an index submodule or of the K-space.

For completeness the other leaf assignment does not preserve L.
Swapping Te1 and Te2 changes R by
`2(1-phi)(h+By)` modulo4L. Its e0 coordinate is `2phi mod4`, and multiplication
by the unit phi cannot clear the denominator4 in the Te0 formula.
This map is not an alternative integral period. No independent D9
leg-swap integrality is assumed.

## 4. Full chamber stabilization and infinite order

Original-module preservation implies that T permutes the **full**
integral reflection arrangement: it conjugates every integral reflector
and its inverse, not only the known finite wall set. At the accepted
actual-germ boundary, T sends the full outward OLD simple germ and
its finite corner to the NEW full outward simple germ and its corner.
Their local chamber sectors therefore agree. Two arrangement chambers
whose interiors share a local sector are the same chamber; hence
`T(C)=C` for the literal full chamber C.

There is a short obstruction to T having finite order which avoids a
trace computation or a spectral classification. From (1),

    Tn=-e3, T²n=s, n²=2, n.Tn=0,
    n.s=-80-132phi<-6.                               (5)

The last product follows, for example, from the exact reflection
covector of s in the ORIGINAL basis,

    (-530-859phi,1;0,32+52phi,0,0,0,0,0,1).

If a future Lorentz isometry F has finite order, averaging a future
timelike orbit gives a fixed future unit vector y. Write `z=ay+w`,
`w perpendicular y`, for any `z²=2` with `z.Fz=0`. On the positive
space perpendicular to y, put `W=w²=a²+2`; then `w.Fw=a²`.
Cauchy applied to w and `Fw+F^-1w` gives

    4a^4 <= W(2W+2w.F²w),
    z.F²z >= -2-4a²/(a²+2) > -6.                    (6)

Taking F=T,z=n contradicts (5). Thus T has infinite order.

The admissible arithmetic input for this literal O-module says its
future integral orthogonal group G is discrete and has finite
covolume. Discreteness can also be seen directly: the conjugate form
is positive definite, so a compact physical set of isometries has
bounded matrices at both embeddings; integrality then leaves finitely
many matrices. Let W_ref be the group generated by ALL integral
reflections. It is normal in G. Its action on full arrangement chambers
has trivial chamber stabilizer, whereas distinct powers of T stabilize
C. Therefore `[G:W_ref]` is infinite. The finite-covolume subgroup/index
formula gives infinite covolume for W_ref and infinite volume for C.

**Conditional proved conclusion:** once both full actual-germ premises
above are accepted, this literal norm79 host is nonreflective. The
argument does not exclude arbitrary reflection subgroups, other
invariant lattices, the K-similarity/commensurability class, a field,
or an infinite binary family. It constructs no compact H9 example and
changes no dimension bound. Independent review and root acceptance
are still pending; the six papers stay frozen.

## 5. Prepared exact reproduction and algorithm contract

The source-only checker is
`scripts/check_compact_h9_prime79_two_d9_germ_period_20261021.py`,
12,107 bytes/300 lines, SHA-256
`186da22652389ea7eb88fca17045c0c375d96c96b4b05b6d792466f1869d4944`.
**It is unexecuted and pending independent source review.** This proof
was derived first; the checker is a small exact reproduction of its
two finite linear maps, not a missing mirror census.

Every field element is a pair of rational numbers for `a+bphi`.
Multiplication reduces `phi²=phi+1`; inverse uses the nonzero rational
norm. Exact physical signs use `2a+b+b sqrt5` and a squared rational
comparison when its summands have opposite signs. No approximate
square root, matrix library, symbolic import or worker import is used.

The checker hashes its own source, this proof (required CLI hash),
all five direct inputs above, the unchanged resource wrapper, and the
conflict-free entire nested input map from the inherited END0 header.
All pins are rechecked before output. It reads the bounded16,733-byte
JSONL, verifies its complete END0/all-profile/unique-id4 record, trial
count and matching actual terminal ledger/output digest. It does NOT
rerun that interval or derive geometric completeness from these fields.

It constructs the ORIGINAL10-dimensional Gram, old/new root vectors
and the two equal-norm negative vectors from their formulas. It checks
both ordered D9 Grams, every corner orthogonality, the exact negative
normalizations and future signs. For EACH of the two target-leaf
assignments, it obtains the unique map as

    T = [new simple roots | v] [old simple roots | u]^-1.

Rational-field Gauss elimination has exactly ten pivot stages and
fails if a pivot is absent; therefore it neither assumes nor hides
invertibility. The checker verifies full ambient Gram preservation,
both inverse identities, all100 entries of T and its inverse for
O-integrality, all ten residual mod4 coordinates, the two-step root
orbit and its scalar obstruction (5). It requires the manual outcomes
integral/not-integral for the two assignments, rather than selecting
an assignment by row order. Both complete matrices and inverses,
their traces, negative vectors and residue vectors are retained.

There are two maps, three ten-by-ten inversions, and a fixed finite
number of scalar field operations. Hashing streams64KiB blocks;
storage is bounded by these matrices and the small pinned JSONL.
There is no root enumeration, square search or unknown collision
minimum. Field arithmetic and finite elimination prove termination;
the explicit two leaf permutations prove domain completeness for
this ordered D9-to-D9 isometry problem. Neither gives completeness
of all possible chamber isometries or any new edge census.

Output is a NEW exclusive-created file directly under `results/`;
existing evidence is never overwritten. The worker requires the live
bounded owner's environment and does not launch children. Only root
may admit it under the existing40GiB/6GiB start gate,59s/1GiB limit,
8GiB/3GiB abort floors and shared serial lock. After independent
source review, the exact reproducibility command is the following,
with `PROOF_SHA256` set to this frozen file's reported digest:

```sh
python3 -B scripts/run_bounded_research_check_20261018.py \
  --ledger results/prime79_two_d9_germ_period_execution_20261021.json \
  --parent-lock -- \
  scripts/check_compact_h9_prime79_two_d9_germ_period_20261021.py \
  --verify \
  --expected-source-sha 186da22652389ea7eb88fca17045c0c375d96c96b4b05b6d792466f1869d4944 \
  --proof-sha "$PROOF_SHA256" \
  --out results/prime79_two_d9_germ_period_20261021.json
```

The prospective PASS label is
`PASS_EXACT_TWO_D9_ORIGINAL_MODULE_PERIOD_ALGEBRA`; it means only
these exact checks at the bound inputs. The output explicitly retains
the external actual-germ and arithmetic finite-covolume premises and
does not assert a new example, dimension or form-class exclusion.
No such checker execution has occurred in preparing this candidate.

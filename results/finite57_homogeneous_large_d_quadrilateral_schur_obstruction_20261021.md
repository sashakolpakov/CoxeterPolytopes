# Finite57 continuation: no homogeneous large-D long-end/center quadrilateral

21 October2026. **New manual supporting candidate, pending root and
outside-author review.** No worker, enumeration, scientific or symbolic
import, test, Lean or TeX process was run. No current paper, rule, guide,
receipt or global bound was changed. This is a restricted actual-face
obstruction, not a reduction of the formal finite maximum372.

## 1. Exact result and mathematical boundary

Let F be an ACTUAL quadrilateral two-face of a hyperbolic Coxeter
polytope in dimension n. Let r>=7, k=n-r, and suppose:

1. its full supporting-white diagram is A_(r-4)+2A1+Z, where Z is a
   positive definite forest of rank k;
2. Z is orthogonal to all four black sides, so it is the SAME unchanged
   whole spectator forest at every corner;
3. at each of the four corners the other component is D_r, and the two
   black marks are its long endpoint and its fork center.

**Theorem.** These hypotheses are inconsistent.

Every corner is finite under these hypotheses. There is no restriction
on arithmeticity, labels or types within Z, or on its number of factors.
The full white rank and persistence of Z are essential: an unsealed
subdiagram at one corner alone does not supply them. The theorem uses
actual roots in one Lorentz space, not a purported realization of a
formal Gram matrix. It needs no field, norm-integrality or lattice premise.

Source-bound context, read as written mathematics:

| Source | SHA256 | Role |
|---|---|---|
| `manuscripts/finite_volume57.tex` | `c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0` | Actual face/white-rank convention, norm2 acute roots, current finite weights and shape coefficient |
| `results/finite57_complete_context662_source_composition_20261021.md` | `7e53dc180baf41bce6d5cab815426a3391ca2de40fc2aa0849c4a6b58226664f` | Exact scope of the already financed selected producer; NOT an input to the Schur calculation |
| `results/finite57_current_a17_endpoint_producer_scope_audit_20261021.md` | `c9e2107004641fd3bcab1fc53e1b797aa264b38cf110e6e83f7bc0e65f8f1d71` | Current A17 provenance/372 boundary; NOT a substitute for the proof below |
| `results/fresh_d17_node3_large_d_guard_transfer_20261018.md` | `eac18fe44e8e1d76175adfbf7446a73ff16c3f104bf92d7f4170edec4138ffa4` | Earlier held guarded savings, distinguished from this new long-end/center argument |

The only general geometric facts needed for the theorem are the actual
positive rank n-2 white support, a Lorentz space of dimension n+1 and
norm2 acute Coxeter wall normals. For any two distinct such walls, minus
their product is either0 or at least1: an elliptic pair has its ACTUAL
intersection by the positive projection lemma and hence native angle
pi/m, m>=3; a parallel/ultraparallel pair has magnitude at least2.
No claim that a general face inherits Coxeter angles is used.

## 2. Exact marked white block and all its products

Put m=r-4>=3. Number the A_m white path q1,...,qm and the two isolated
white leaves f1,f2. All roots have square2; consecutive path products
are-1. List the four cyclic black sides as a,b,c,d, with a,c in the
long-end role and b,d in the center role. The hypotheses force

    a.q1=c.q1=-1,
    b.qm=d.qm=b.f1=b.f2=d.f1=d.f2=-1,                 (1)

and all other black/white products in this block are zero. Each
adjacent black product is also zero:

    a.b=b.c=c.d=d.a=0.                               (2)

Indeed in D_r a marked long end has one neighbor, whereas the marked
center has three. Deleting these two marks leaves exactly this path
and two leaves. With b fixed, white+b is D_(r-1), so a new long-end
completion must attach at q1. With a fixed in the long-end role, a
center completion must join qm and both isolated leaves. Thus (1)
is not a relabeling that forgets the white identifications. At r>=7
the two marked nodes are nonadjacent, giving(2).

Write the opposite products as

    a.c=-g, b.d=-h,   g,h>=0,
    g,h each belong to {0} union [1,infinity).       (3)

Their precise finite labels or disjoint distances are unrestricted.
All roots of Z are orthogonal to this entire displayed core, by the
unchanged-spectator hypothesis.

## 3. Four-by-four Schur complement and the native gap

Let B be the white Gram of A_m+2A1. The path inverse has

    (B_A^-1)11=(B_A^-1)mm=m/(m+1),
    (B_A^-1)1m=1/(m+1),

while each isolated leaf has inverse diagonal1/2. Schur elimination
of B gives black diagonal entries

    S_aa=S_cc=(m+2)/(m+1), S_bb=S_dd=1/(m+1).

Every adjacent entry is -1/(m+1), and

    S_ac=-g-m/(m+1),
    S_bd=-h-(2m+1)/(m+1).                            (4)

These entries include BOTH isolated-leaf contributions in S_bd;
omitting them would change the conclusion.

Use the orthonormal opposite-difference and opposite-sum coordinates.
The two difference eigenvalues are2+g and2+h, strictly positive.
The remaining two-by-two sum block is

    T = [ 2/(m+1)-g       -2/(m+1)       ]
        [ -2/(m+1)       -2m/(m+1)-h     ].          (5)

Consequently det S=(2+g)(2+h)det T, where direct expansion gives

    det T = ((m+1)g*h+2m*g-2h-4)/(m+1).             (6)

If g=0, the numerator is -2h-4<0. If g>=1, it is

    h*((m+1)g-2)+2m*g-4 >= 2m-4 >0,

since m>=3. The native gap in(3) exhausts all possibilities, so
S is always nonsingular. In fact the g>=1 branch is even impossible
without the final rank count: both diagonals of T are negative and
its determinant positive, giving two negative directions. No real
root or numerical approximation enters the calculation.

## 4. Why actual dimension forces the missing zero eigenvalue

The core has m+2 white roots and four black roots, hence r+2 roots.
B is invertible and S is nonsingular, so its Gram has rank r+2.
The unchanged orthogonal Z contributes another k independent positive
directions. Thus the full actual Gram has rank r+2+k=n+2.
But ALL these vectors lie in the original Lorentz space of dimension
n+1. Its Gram rank cannot exceed n+1. This contradiction proves
the theorem, including g=0. Equivalently, an actual two-face with
precisely this complete orthogonal white support would require det S=0,
which(6) forbids.

The argument is not a statement that every indefinite core radical
would be a feasible cusp. It avoids that issue: the required radical
is absent. Nor is an independent energy at different corners added.

## 5. What this does and does not change in the ambient57 budget

For r17, the relevant marked pair is(0,14) in the current D17
numbering; its current weight is1/2. With unchanged spectators
Z=D17+D23 this theorem excludes a quadrilateral whose FOUR active
corners all have that same long-end/center role. The conclusion is
about the full carried white block, not an arbitrary four appearances
of D17. A side changing its marked role, a split active corner, an
ideal corner or spectator attachment remains outside this theorem.

The current finite cost calculation is still

    C0(2D17+D23)=2*(227/2)+145=372.

The simple-vertex shape coefficient at ambient57 is364: for a
56-simplex G0=28²=784 and (5*784-56*57)/2=364. Over the fixed
inert D17+D23 charge258.5, the active rank17 allowance would have
to be at most105.5. The new obstruction supplies NO such uniform
allowance. In particular another active A17 still costs453/4=113.25,
giving total371.75 and deficit7.75. Even a hypothetical whole endpoint
half deletion by itself leaves deficit7.25. The earlier held large-D
node3 guard has different marked pairs and its own extra coupons;
it must not be counted here as a newly proved saving.

A complete reallocation must control every other actual corner reached
by EACH changed marked role and the number of debits/credits per actual
vertex. Nonhomogeneity alone supplies neither a cheap neighbor nor a
once-only injection. The existing Context662 allocation applies only
to its chosen actual faces and ideal sources; it is not a universal
363.75 allowance for these other vertices. There is presently no proof
that enough good finite/ideal recipients are forced to absorb the full
eight-unit deficit at each heavy all-D vertex, including competition
with heavy A17 corners.

**New proved scope:** the complete homogeneous long-end/center
quadrilateral obstruction above. **Unresolved scope:** a universal
reallocation across nonhomogeneous corners, the full rank57 finite
maximum and a global56 contradiction. No charge372 improvement or
new ambient bound is claimed, and all six manuscripts remain frozen.

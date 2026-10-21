# The ramified silver diagonal has an odd-silver facet bridge

20 October 2026. Own paper proof, awaiting root and outside review.
No computation, root enumeration, import, Lean or TeX execution. Put
p=√2, α=1+p, O=Z[p], and

    M_n=<-p>+I_n,   q0=2−p=p/α,   X=2+p=pα.

**Candidate literal-family conclusion:** for n≥10 the full integral
reflection arrangement of M_n has no compact chamber. This does not yet
cover arbitrary lattices in the same quadratic space or actual node-two
stars preserving such other lattices.

## 1. Complete norm pool and the exact facet over-lattice

The complete primitive M_n norm pool is {1,2,q0}, up to unit-square
rescaling. This is the pool argument in Section 3.1 of
[compact_h9_odd_silver_diagonal_two_edges_20261020.md](compact_h9_odd_silver_diagonal_two_edges_20261020.md),
SHA256 `85095d203aea614961425bf356cfbad21529cd9f6f01d32b77f73cc5b2ce234c`.
Here is the needed argument explicitly. For a primitive root v=(a;x),
reflection integrality gives q|2p. Total positivity and the units ±α^Z
leave valuations 0,1,2,3 at (p). Valuation 3 would force all xi∈pO;
primitivity then forces a to be a dyadic unit. Its norm contribution
−pa² has valuation 1, whereas every spatial square has valuation at
least 2. It cannot have valuation 3. The three surviving classes can
therefore be represented by 1,q0,2. All these norms divide 2, so an
integral vector of any displayed norm is reflective. Replacing q0 by
X=α²q0 is merely a unit-square change of representative.

Use the integral reflector

    c=(1;1,1,0,...,0),   c²=q0,   2/q0=pα∈O.

In c-perp the exact integral lattice N=M_n∩c-perp satisfies
x1+x2=pa. Set

    v0=(1;p,0,0,...),   u=(0;1/p,−1/p,0,...).

Then N has basis v0,pu,e3,...,en, and the index-p over-lattice

    N~=N+O u

has the integral orthogonal basis

    h=α(v0−u), u, e3,...,en,
    h²=−α, u²=1, h.u=0.

Thus N~ is literally isometric to L_(n−1)=<-α>+I_(n−1). No
identification of the ambient indefinite lattices M and L is asserted.

## 2. Complete crossing-product table

Let v be a primitive ambient reflector of norm q∈{1,2,q0}, and put
t=c.v∈O. If its hyperplane crosses c-perp and is independent of c,
the physical two-wall Gram and the other-embedding form give

    t²<q q0,   (t')²<q' q0'.                            (1)

Writing t=A+Bp, addition gives
2A²+4B²<Tr(q q0). The latter trace is 4,8,12 respectively. Testing
the individual physical and conjugate bounds in this tiny coefficient
set gives the complete table:

| q | all possible t |
|---|---|
| 1 | 0 |
| 2 | 0, ±1, ±α^(-1) |
| q0 | 0, ±α^(-1) |

This table concerns arbitrary intersecting arrangement mirrors. In
particular the q=2, |t|=α^(-1) case is retained even though its acute
angle is 3π/8 rather than a native fundamental angle.

Project v to c-perp: vbar=v−(t/q0)c. For t=0 it is already integral
in N⊂N~, of norm dividing 2, hence defines an integral N~ reflection.
For nonzero t, the complete possibilities, up to signs, are:

| q,t | p vbar | norm of p vbar |
|---|---|---|
| 2, 1 | p v−αc | q0 |
| 2, α^(-1) | p v−c | X |
| q0, α^(-1) | p v−c | q0 |

Each displayed vector belongs to N, and each norm divides 2. Since
N~ is integral, its reflection preserves N~. Therefore every induced
mirror on c-perp belongs to the standard full L_(n−1) arrangement.
Only a SUBSET assertion is needed; extension of all intrinsic mirrors
to the ambient lattice is neither used nor claimed.

## 3. The full-arrangement compactness contradiction

The c mirror bounds a full chamber at any generic point of c-perp
avoiding the other mirrors. Such points exist by local finiteness.
It need not be the first wall of the B_n germ, and no first-wall
assertion is made. Full-arrangement chambers are transitive under
their reflection group. Thus hypothetical compactness of a full
chamber would make this chamber, and its actual c facet F, compact.

Every boundary mirror of F is covered by Section 2. Choose a generic
interior point of F and its full standard N~ chamber. This entire
standard chamber satisfies every selected bounding halfspace of F,
so it is contained in F. The inherited standard odd-silver theorem
says its volume is infinite when n−1≥9:

[odd_silver_stable_integral_isometry_and_uniform_obstruction_20261020.md](odd_silver_stable_integral_isometry_and_uniform_obstruction_20261020.md),
SHA256 `2ad82b9a52c7b81364a070b985932932eab72438922e7b0d717b392e5e4c5ee3`.

That theorem's stable integral identification, strong-approximation
input and old full E8-chamber argument remain written imported inputs.
The subset containment contradicts compactness of F. It does not assume
that an arbitrary metric facet is itself a Coxeter polytope.

This proves the stated compact literal-family obstruction, subject to
the requested full mathematical review. Using the standard arithmetic
lattice theorem, a finite-covolume full reflection subgroup here would
have finite index in the anisotropic arithmetic group O(M_n), and
would also be cocompact; the same argument would therefore give its
nonreflectivity. No execution of a chamber search is involved.

The outstanding node-two obligation is the invariant-lattice interface:
an actual star fixes this K-quadratic space but need not preserve the
literal M_n. Its pair-cap dyadic operators differ from the singleton
and triple operators, so the earlier flag theorem cannot simply be
reused. No new uniform dimension ceiling is claimed.

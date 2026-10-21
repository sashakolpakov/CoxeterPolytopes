# Good-prime two-face trees: a local restriction and its limit

20 October 2026. New written research lemma, awaiting independent review.
No scientific computation, census, optimizer, Lean or build was run. No
dimension bound or arithmeticity theorem follows from this note.

## 1. Precise local hypothesis

Let F be a nonarchimedean local field of odd residue characteristic, with
valuation v, valuation ring O, and a nondegenerate symmetric space (V,b).
Let e_1,...,e_r span a nondegenerate subspace W, and put U=W-perp. Assume:

* each b(e_i,e_i) is an O-unit;
* every b(e_i,e_j) is integral, and every nonzero such entry is a unit;
* the white Gram matrix C=(b(e_i,e_j)) has unit determinant.

The connectivity graph uses the nonzero off-diagonal entries. Let s_i be
the reflection in e_i and H the finite group they generate. Finiteness is
needed for the Coxeter application; the splitting assertion below only
uses invariance under the individual s_i.

Use the standard model of the orthogonal Bruhat--Tits building by
self-dual additive valuation norms. Its points are splittable norms alpha
with alpha(av)=v(a)+alpha(v), alpha(v+w)>=min(alpha(v),alpha(w)), and

    alpha*(x) = inf_y [v(b(x,y))-alpha(y)] = alpha(x).

The building-model input, its norm duality and its apartment description
are explicitly stated in Cornut, *On p-adic norms and quadratic
extensions, II*, Section3, printed pages6--8:
https://webusers.imj-prg.fr/~christophe.cornut/papers/npeq2.pdf .
This is the external geometric input; the white splitting below is proved
directly. The action of the full orthogonal group, including reflections,
is used. There is no claim that a black reflection preserves U.

## 2. Exact white fixed set

**Lemma.** Every self-dual norm invariant under all s_i is exactly

    alpha(w+u)=min(alpha_W(w),alpha_U(u)),
    alpha_W(sum a_i e_i)=min_i v(a_i),

where alpha_U is an arbitrary self-dual norm on U. Consequently Fix(H) is
the embedded building of O(U), with the standard white norm fixed.

**Proof.** The rank-one projector

    P_i=(1-s_i)/2,   P_i x=b(x,e_i)e_i/b(e_i,e_i)

is contractive for an invariant alpha: alpha(P_i x)>=alpha(x). Write
c_i=alpha(e_i). Across a nonzero white bond, applying P_i to e_j and P_j
to e_i gives c_i>=c_j and c_j>=c_i, because both coefficients are units.
Thus c_i is a constant c on each connected component.

On one component, for w=sum a_i e_i, the ultrametric inequality gives
alpha(w)>=min_i v(a_i)+c. Conversely contractivity bounds all
v(b(w,e_i))+c below by alpha(w). The inverse of that component's Gram
matrix is integral, so it bounds min_i v(a_i)+c below by alpha(w).
Equality follows. Every component Gram is unimodular, since C is a block
diagonal integral matrix with unit determinant.

For a general x, the same integral inverse and contractivity show that
the orthogonal component projector P_W is contractive. Hence 1-P_W is
also contractive, and alpha splits orthogonally over W+U. On a white
component the dual norm has constant -c, because its Gram is unimodular.
Self-duality therefore forces c=0. The restriction to U is self-dual.
Conversely these orthogonal sums are self-dual and white-invariant. This
proves the assertion, including the componentwise case. The norm
description identifies the fixed set with the complement building. It is
convex as a fixed set in the ambient CAT(0) building. QED.

In particular, when dim(U)=3, this fixed set is a tree or a point: its
apartments have dimension at most one. This conclusion does not require
the residue characteristic to be prime to |H|. For simply laced raw
norm-two white diagrams, odd primes dividing neither their Cartan
determinants nor any nonzero bond satisfy the hypotheses. For example
A_r requires p not dividing r+1; D_r, E7 and E8 satisfy them at every odd
prime. Primes above2 are deliberately excluded.

## 3. An actual compact two-face

Suppose the representation of a compact hyperbolic Coxeter polytope is
defined over F. Let e_i be the n-2 actual walls through a two-face, and
let r_1,...,r_k be its cyclic side walls. Impose the hypotheses above on
the white walls. Their orthogonal complement has dimension three. Put

    Y=Fix(H),     T_j=Y intersect Fix(s_{r_j}).

Each T_j is a nonempty closed subtree: H together with a side reflection
is the finite stabilizer of the corresponding edge, so has a fixed point
in the ambient building. Consecutive T_j meet, because the white walls
and those two side walls form a finite vertex stabilizer. These statements
use the bounded-orbit fixed-point property of complete CAT(0) buildings.

The intersection graph of finitely many subtrees is chordal and has the
clique Helly property. Here is a direct proof. Choose one point for each
nonempty intersection and take their finite convex hull in Y; intersect
each subtree with this hull. Root the resulting finite tree. A subtree
whose nearest-to-root point is deepest has that point in every subtree
meeting it. Its neighbors are therefore a clique. Remove it and repeat.
This is a perfect elimination ordering, proving both assertions. The
point case is immediate.

It follows that the cyclic k-gon admits a combinatorial triangulation
whose k-3 added diagonals join pairs T_i,T_j with nonempty intersection.
Every triangle of that triangulation has a common fixed point for H and
its three side reflections. In particular:

* a triangular two-face star generates a bounded p-adic subgroup;
* a quadrilateral has at least one bounded opposite-pair extension;
* two opposite black pairs cannot both have nonintegral product trace.

For the last assertion, if d=dim(V) and q_i=b(r_i,r_i), then

    tr(s_{r_i}s_{r_j})=d-4+4 b(r_i,r_j)^2/(q_i q_j).

A bounded matrix subgroup preserves a lattice, so its traces are
integral. Negative valuation of this displayed trace thus precludes the
corresponding fixed-set intersection. This criterion concerns the full
white-plus-black group, not just its Schur complement.

The statement above is made for compact polytopes. Extending it to ideal
vertices requires explicitly checking boundedness of their parabolic
representations at the selected place; real parabolicity alone is not
being silently substituted for the finite-group argument here.

## 4. Sharp local barrier to a dimension argument

The lemma does not force a bad pair to occur only on long faces. For
every odd rational prime p, set

    x=(p^2+1)/(2p), A=(p^2-1)/(2p), c=1/2,
    y=sqrt(p^4-p^2+1)/(p^2-1), T=c x/A.

Thus A^2=x^2-1 and (x^2-1)(y^2-1)=c^2. In Lorentz coordinates with
form -t^2+X^2+Y^2 take four outward unit normals

    n1=(0,1,0),       n2=(0,0,1),
    n3=(A,-x,0),      n4=(T,-c,-y).

Their cyclic products are 0,0,0,-1/2, their opposite products are -x,-y,
and their Gram determinant vanishes. The cyclic angles are three right
angles and pi/3. In the time-one Klein chart their four corner points are

    (0,0), (-A/x,0), (-A/x,-c/(Axy)), (0,-T/y).

Every corner is strictly inside the unit disk. For the third corner the
needed inequality is c^2<A^2 y^2=A^2+c^2; for the fourth it is c^2<1.
Each nonincident wall inequality is strict there. A small point with
both spatial coordinates negative is strict for all four walls. Hence
this is an actual compact Lambert quadrilateral, with no orientation or
formal-Gram ambiguity.

The number field Q(sqrt(p^4-p^2+1)) has an embedding into Q_p, since its
radicand is1 modulo p. At either such local embedding, v_p(x)=-1 and
v_p(y)=0. Thus one opposite product trace is nonintegral and the other is
integral. Appending any orthogonal good finite white root system gives
the same actual local two-face in arbitrarily large ambient dimension,
with genuine finite corner groups. The ambient local chamber is not
claimed to have finite volume. This is precisely the remaining global
issue: a good-prime bad pair can be fully compatible with a quadrilateral.

## 5. Research boundary

The fixed-set lemma is a new local tool proposed for review, not a
uniform improvement of57 or29. A useful global argument would have to
couple these triangulations across different two-faces, control the bad
white primes and primes above2, and use finite volume essentially. In
particular it cannot assume that all black reflections act on one common
three-dimensional complement. A possible tree-action construction from
the whole dual two-skeleton also requires a compatible global graph
model; that step, any resulting splitting theorem and any quantitative
dimension consequence remain unproved.

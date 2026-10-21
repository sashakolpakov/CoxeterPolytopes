# Exact weak full-Picard K3 carrier theorem

## Verdict on the proposed equivalence

The fixed-lattice implication in `GOAL.md` is correct, but the Coxeter-level
wording

> the wall lattice (or a reflection-stable overclosure/scaled dual) embeds
> primitively in the K3 lattice

is not an exhaustive `if and only if` as written.  It becomes exhaustive only
after one fixes, or exhaustively enumerates, the **integral root datum**: the
positive rational scaling chosen on every wall-normal line.  For one fixed
integral root datum, reflection-stable even overclosures are exhaustive.  A
scaled dual is a useful way to find a different root datum, but is not by
itself a theorem that all possible root data have been found.

The safe Coxeter-level statement is therefore existential: a cone has a weak
full-Picard complex-K3 carrier if and only if it has **some** oriented even
integral reflection realization `N` that embeds primitively in the K3
lattice.

## Definitions

Let `P` be a finite-sided hyperbolic Coxeter cone in a real quadratic space
`(V,q)` of signature `(1,n)`, with a specified positive-cone component and
specified positive side of every wall.

An **oriented even integral reflection realization** of `P` consists of:

1. an even nondegenerate lattice `N` of signature `(1,n)`;
2. a similarity
   
   ```text
   phi : (V,q) -> (N_R,( , )),       (phi(x),phi(y)) = c q(x,y), c>0,
   ```
   
   carrying the chosen positive component to a positive component of `N_R`;
3. for every facet, a primitive vector `v_i in N`, with `v_i^2<0`, such that
   `phi` takes that wall to `v_i^perp` and takes its specified positive side
   to `(x,v_i)>=0`; and
4. integrality of every displayed reflection
   
   ```text
   s_i(x)=x-2(x,v_i)/v_i^2 v_i in O(N).
   ```

For primitive `v_i`, the last condition is the exact divisibility test

```text
v_i^2 divides 2 div_N(v_i),
```

where `div_N(v_i)=gcd{(v_i,x):x in N}`.  When the original Coxeter data have a
rational structure, `phi` is required to be rational as well.

A **weak full-Picard K3 carrier** is a projective complex K3 surface `X` and a
marking `NS(X) ~= N` for such a realization.  “Weak” means only that the
mirror vectors are algebraic divisor classes and their reflections are
integral on the full Picard lattice.  It does not mean that the vectors are
effective, irreducible, nef walls, or induced by surface automorphisms.

## Theorem (fixed integral realization)

Let `N` and the oriented vectors `v_i` be an oriented even integral
reflection realization as above, and put `r=rank(N)=n+1`.  Then the following
are equivalent.

1. There is a weak full-Picard complex-K3 carrier `X` with `NS(X) ~= N`.
2. There is a primitive lattice embedding
   
   ```text
   N -> Lambda_K3 := U^3 + E8(-1)^2.
   ```

Necessarily `1 <= r <= 20`.  If these conditions hold, the locus with exactly
`NS(X)=N` has complex dimension `20-r` for each primitive-embedding/period
component.  For `r<20` it is the relevant type-IV period domain with the
extra-Noether--Lefschetz hyperplanes removed.  For `r=20` it is
zero-dimensional and the orthogonal complement is a positive-definite even
rank-two lattice.

### Proof

If `X` is a projective complex K3 surface, `H^2(X,Z)` is the even unimodular
K3 lattice of signature `(3,19)`.  The subgroup
`NS(X)=H^2(X,Z) cap H^(1,1)(X)` is primitive.  Hodge index gives signature
`(1,r-1)`, so a full-Picard carrier supplies the required primitive
embedding.  The remaining reflection and orientation conditions are part of
the definition of the fixed integral realization.

Conversely, fix a primitive embedding `N -> Lambda_K3` and put
`T=N^perp`.  Then

```text
signature(T)=(2,20-r).
```

Choose a period

```text
[omega] in Omega_T = {[omega] in P(T_C): omega^2=0,
                                          (omega,bar(omega))>0}
```

outside every hyperplane imposed by a class of `Lambda_K3` not in `N`.
For `r<20`, these are countably many proper analytic hyperplanes and cannot
exhaust `Omega_T`.  For `r=20`, `T_R` is a positive plane and its period has
no nonzero real orthogonal vector, so the same conclusion is automatic.
Thus

```text
omega^perp cap Lambda_K3 = N.
```

Surjectivity of the K3 period map gives a marked complex K3 surface with this
period, hence with `NS(X)=N`.  Since `N` is hyperbolic it contains a class of
positive square, so `X` is projective.  The integral `(1,1)` classes are
divisor classes, and the already-certified `s_i` act integrally on the full
Picard lattice.  Composing the marking with `-id` if necessary aligns the
chosen positive-cone component.  The similarity and positive wall scalings
make the original oriented Coxeter cone literal.

This proof uses period surjectivity, not strong Torelli: no automorphism is
being constructed.

## Exact primitive-embedding decision

Let

```text
A_N = N^vee/N,       q_N : A_N -> Q/2Z
```

be the discriminant form.  The following is an exact necessary-and-sufficient
criterion, including all ranks outside Nikulin's stable range.

> `N` embeds primitively in `Lambda_K3` if and only if there is an even
> lattice `T` of signature `(2,20-r)` whose discriminant form is isometric to
> `-q_N`.

Indeed, a primitive sublattice of a unimodular lattice has the negative
discriminant form of its orthogonal complement.  Conversely, glue `N+T`
along the graph of an anti-isometry `A_N -> A_T`.  The result is even,
unimodular, and has signature `(3,19)`, hence is the K3 lattice; both summands
are primitive.

This gives a finite decision algorithm.

1. Check that `N` is even, nondegenerate, has signature `(1,r-1)`, and
   `r<=20`.
2. Compute the Smith normal form of a Gram matrix and the finite quadratic
   form `q_N` from its inverse modulo `2Z`.
3. Set `t=22-r`, `q=-q_N`, and decompose `q` into its finitely many primary
   Jordan forms `q_p`.  If `ell(A_{N,p})>t` for some prime `p`, reject.
4. Check Milgram's signature congruence
   
   ```text
   Br(q) = 2-(20-r) (mod 8).
   ```
   
   Here it is automatic once `q=-q_N` has been computed correctly, since
   `Br(q_N)=signature(N) (mod 8)`.
5. Apply Nikulin's even-lattice existence theorem (Theorem 1.10.1) to the
   requested invariants `(2,20-r,q)`.  Concretely, for every odd `p` at which
   `t=ell(A_{N,p})`, the determinant square class of the minimal `p`-adic
   lattice inducing `q_p` must equal the determinant square class
   `(-1)^(20-r)|A_N|`.  At `p=2` with `t=ell(A_{N,2})` and
   `q_2` not of the exceptional form `q_theta(2)+q'_2`, Nikulin's condition is
   
   ```text
   |A_N| = +/- disc(K(q_2)) in Z_2^*/(Z_2^*)^2.
   ```
   
   If `q_2=q_theta(2)+q'_2`, there are two minimal 2-adic lattices and no
   additional determinant condition.  If `t>ell(A_{N,p})`, unimodular padding
   removes that local determinant obstruction.  Only primes dividing
   `2 det(N)` occur.

Step 5 should be implemented by constructing the finite set of rank-`t`
`p`-adic determinant classes inducing `q_p`, rather than by applying only a
length inequality.  This is precisely where the borderline odd-prime and
2-adic failures live.  It is finite and necessary-and-sufficient, not merely
a sieve.

The familiar stable-range shortcut is a sufficient special case:

```text
r <= 19 and 22-r >= ell(A_N)+2.
```

Then Nikulin's primitive-embedding theorem gives an embedding, in fact a
unique embedding orbit in the usual stable hypotheses.  Equality in this
display is still inside the stable range.  When it fails, one must run the
local test above; failure of the shortcut is not failure of K3 existence.

If embedding or moduli **components** rather than mere existence are needed,
also enumerate the isometry classes of complementary lattices `T` and the
finite double cosets of anti-isometries of discriminant forms, modulo the
images of `O(N)` and `O(T)`.  Different primitive-embedding orbits need not
give the same marked carrier component.

## Why overclosures are only conditionally exhaustive

Fix primitive target normals `v_i` and let

```text
L = sum_i Z v_i.
```

Because the cone is full-dimensional, `L` has full rank.  Every full-Picard
ambient lattice containing these same vectors is then a finite-index even
overlattice

```text
L subset N subset L^vee,
```

and integrality of the facet reflections says exactly that the corresponding
isotropic subgroup of `(A_L,q_L)` is stable under all facet reflections.
Thus an exhaustive stable-overclosure census is exhaustive **for this fixed
normalization**.

But a wall hyperplane does not choose a primitive integral normal.  Replacing
the normal directions by positive rational multiples can produce a different
crystallographic root datum before any overclosure is taken.  These choices
must satisfy the integral Cartan conditions and preserve the specified sides.
A scaled dual `M^vee(c)` is one particularly useful such realization, provided
one proves:

- evenness and integrality of the target lattice;
- an orientation-preserving positive similarity of the quadratic spaces;
- positive rescaling of every wall normal; and
- integrality of every transported reflection.

Without an exhaustion theorem for all admissible normal scalings, the phrase
“wall lattice or an overclosure/scaled dual” is a sufficient search list, not
an `if and only if`.

## Exact promotion boundaries on a K3 surface

The weak theorem establishes only level 1 below.

1. **Algebraic divisor classes.**  This is automatic from `NS(X)=N`.
2. **Effective mirrors.**  If `v^2=-2`, K3 Riemann--Roch implies that exactly
   one of `v` and `-v` is effective.  For a fixed marking and ample class `a`,
   the specified sign `v` is effective exactly when `(a,v)>0`.  Hence a set
   of oriented `(-2)` mirrors is simultaneously effective exactly when the
   ample cone meets all their strict positive halfspaces.  An effective
   `(-2)` class can be reducible; effectivity does not make its hyperplane a
   nef wall.  For primitive mirrors of square below `-2`, including the
   frequent `(-4,divisibility 2)` reflections, neither sign is forced
   effective by Riemann--Roch.
3. **Ample-intersection.**  For the imported open cone `C^o`, the exact
   condition is `C^o cap Amp(X) != empty`.  Equivalently, there is a class
   `H in C^o cap Nef(X)` with `H^2>0`: if `A` is ample, then `mH+A` is ample
   and remains imported-strict for all sufficiently large `m`.  This does not
   say that any imported wall is a nef wall.
4. **Full nef cone.**  The exact statement is that `C^o` is a connected
   component of
   
   ```text
   C_N^+ minus union_{delta in N, delta^2=-2} delta^perp,
   ```
   
   with the marking chosen so that this component is the ample chamber.
   Equivalently, `C` is the closure of a complete K3 root chamber.  For a
   finite irredundant facet description, it is enough and necessary that
   every facet hyperplane be a primitive `(-2)`-root hyperplane and that no
   other `(-2)`-root hyperplane meet `C^o`; phrasing the conclusion as the
   complete-root-chamber condition avoids hidden orientation or redundancy
   assumptions.  A primitive `(-4)` facet immediately forbids full-nef
   equality, though it remains a valid divisor reflection.

None of these conditions makes a Picard reflection a Hodge isometry, makes it
preserve the ample cone, or produces a surface automorphism.  Those require a
separate discriminant-glue, period, and Torelli argument.

## Logical limits for polytope existence

- Failure for one `N` excludes only that Picard lattice.
- Exhausting all root normalizations and stable even overclosures excludes
  weak full-Picard **complex-K3** carriers for that oriented rational
  reflection representation.
- It still does not exclude the abstract real Coxeter polytope, a rational or
  Enriques carrier, a positive-characteristic carrier, or a totally-real
  eigenspace construction.  Nonrational Coxeter forms, notably golden
  `m=5` data, are outside the ordinary full-Picard integral theorem unless an
  integral realization is separately supplied.
- Algebraic-carrier failure implies abstract-polytope nonexistence only after
  a universality theorem proves that every hypothetical polytope in the
  stated class must admit one of the exhaustively classified carriers.
- Conversely, a weak carrier says nothing about compactness or finite volume.
  The Coxeter Gram feasibility and the literal Vinberg finite-volume tests
  remain independent mandatory gates.

## Primary references

- V. V. Nikulin, [*Integral symmetric bilinear forms and some of their
  applications*](https://www.mathnet.ru/eng/im1677), Math. USSR-Izv. 14
  (1980), 103--167: discriminant forms, even-lattice existence, and primitive
  embeddings (especially Theorems 1.10.1, 1.12.2, and 1.14.4).
- Vik. S. Kulikov, [*Surjectivity of the period mapping for K3
  surfaces*](https://www.mathnet.ru/eng/rm3227), Uspekhi Mat. Nauk 32:4
  (1977), 257--258.
- I. I. Pyatetskii-Shapiro and I. R. Shafarevich, [*A Torelli theorem for
  algebraic surfaces of type K3*](https://www.mathnet.ru/eng/im2021), Math.
  USSR-Izv. 5 (1971), 547--588.
- D. Burns and M. Rapoport, [*On the Torelli problem for kaehlerian K-3
  surfaces*](https://numdam.org/item/ASENS_1975_4_8_2_235_0/), Ann. Sci.
  Ecole Norm. Sup. 8 (1975), 235--273.

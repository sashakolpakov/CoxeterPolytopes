# Equivariant minimalization: full trace descent without a complement bound

**Status (2026-08-10): project-native descent theorem and exact sharpness
family.**  No Nikulin--Vinberg, ABSW, Prasad, or arithmetic-volume theorem is
used.  The result removes the minimality hypothesis from the actual-
automorphism lane whenever the carrier has nonnegative Kodaira dimension.  It
also proves that trying to bound the rank of an arbitrary nonminimal
orthogonal complement is the wrong target: even the reflection-maximal
`Gamma16` chamber admits intrinsic correspondence carriers with complement
`<-1>^r` for every `r`.

The replacement invariant is **finite-orbit equivariant minimalizability**.
If the exceptional lattice of one geometric minimal contraction has finite
numerical orbit under the Coxeter operators, finite covolume forces the
entire trace block to be orthogonal to every lattice in that orbit.  The
trace block then descends isometrically through the original contraction,
irrespective of how many exceptional curves occur.  The other orbit members
need not themselves be effective or contractible.  When the original
exceptional lattice is actually preserved, the original integral operators
themselves descend.

## 1. Setup

Let `P` be a finite-volume arithmetic Coxeter chamber in `H^n`, `n>=2`, over
a totally real field `k` of degree `d`.  Let

\[
 T_{\mathbf Q}=\operatorname {Res}_{k/\mathbf Q}V,
 \qquad Q_T=\operatorname {Tr}_{k/\mathbf Q}q,
 \qquad \operatorname {sig}(T_{\mathbf Q})=(1,d(n+1)-1).     \tag{1.1}
\]

Fix a full reflection-stable integral lattice
\(T_{\mathbf Z}\subset T_{\mathbf Q}\).  To lighten notation, write \(T\)
for \(T_{\mathbf Z}\) below, so \(T\otimes\mathbf Q=T_{\mathbf Q}\).
Existence of this integral structure is part of the intrinsic-carrier
hypothesis; the rational restriction-of-scalars space alone is not an
integral lattice.

The facet roots span `V`.  Let `Y` be a smooth projective complex surface and
suppose there is an integral isometric injection

\[
                 \iota:T\hookrightarrow \operatorname {Num}(Y)             \tag{1.2}
\]

and integral pairing isometries `u_i` such that

\[
 u_i(K_Y)=K_Y,
 \qquad u_i\iota=\iota\operatorname {Res}(s_{r_i}).          \tag{1.3}
\]

Write `Gamma=<u_i>`.  As in the minimal-carrier theorem, the minus-one
eigenspaces of the reflections and the spanning property give

\[
                              K_Y\perp\iota(T).                \tag{1.4}
\]

Let `pi:Y->X` be a birational morphism to a smooth minimal surface and put

\[
 E_\pi=\ker\bigl(\pi_*:\operatorname {Num}(Y)\to
                              \operatorname {Num}(X)\bigr).   \tag{1.5}
\]

Factoring \(\pi\) into point blowups and using total-transform exceptional
classes gives the **integral** orthogonal decomposition

\[
 \operatorname {Num}(Y)
   =\pi^*\operatorname {Num}(X)\ \perp\ E_\pi.                \tag{1.6}
\]

In particular, \(E_\pi\) is primitive and unimodular negative definite, and

\[
 \operatorname {Num}(Y)_\mathbf Q
   =\pi^*\operatorname {Num}(X)_\mathbf Q\ \perp\ E_{\pi,\mathbf Q}.
                                                                    \tag{1.7}
\]

Call `pi` **Gamma-equivariant numerically** if every `u_i` preserves
`E_pi`.  This is weaker than requiring the operators to be automorphisms, but
stronger than merely extending their matrices to `Num(Y)`.

## 2. The descent theorem

> **Theorem 2.1 (equivariant exceptional descent).**  In the setup above,
> assume `E_pi` is `Gamma`-stable.  Then
>
> \[
>                    \iota(T)\perp E_\pi.                       \tag{2.1}
> \]
>
> Consequently there is a unique isometric injection
>
> \[
>        \iota_X:T\hookrightarrow\operatorname {Num}(X),
>        \qquad \iota=\pi^*\iota_X,                             \tag{2.2}
> \]
>
> every `u_i` descends to an integral pairing isometry `u_{i,X}` fixing
> `K_X`, and
>
> \[
>       u_{i,X}\iota_X=\iota_X\operatorname {Res}(s_{r_i}),
>       \qquad K_X\perp\iota_X(T).                              \tag{2.3}
> \]
>
> Thus `X` is an intrinsic minimal canonical carrier of the same full trace
> representation and the same labelled Coxeter Gram.

**Proof.**  Stability of `E_pi` and preservation of the intersection form
also make `E_pi^\perp` stable.  Hence orthogonal projection gives a
`Gamma`-equivariant map

\[
                         p_E:T_\mathbf Q\longrightarrow E_{\pi,\mathbf Q}.
                                                                    \tag{2.4}
\]

The integral orthogonal group of a negative-definite lattice is finite.
After replacing `Gamma` by a finite-index subgroup `Gamma_0`, its action on
`E_pi` is therefore trivial.

Restriction to \(T\) maps \(\Gamma\) onto the Coxeter group \(W(P)\).
Therefore the image of \(\Gamma_0\) has finite index in \(W(P)\), and is
itself a finite-covolume hyperbolic lattice at the distinguished real place.
Its standard representation has no nonzero fixed vector. Indeed, a fixed
positive, isotropic, or negative vector would respectively force a fixed
point of hyperbolic space, a fixed ideal point, or an invariant proper
totally geodesic subspace. Each contradicts the full ideal-sphere limit set
of a finite-covolume lattice. Finally, a nonzero element of
\(\operatorname {Res}_{k/\mathbf Q}V\) cannot vanish at the distinguished
embedding: write it in a \(k\)-basis of \(V\) and use injectivity of that
field embedding. Hence

\[
                         T_\mathbf Q^{\Gamma_0}=0.               \tag{2.5}
\]

The trace form identifies invariant linear functionals on `T_Q` with fixed
vectors, so `Hom_{Gamma_0}(T_Q,Q)=0`.  Composing `p_E` with every linear
functional on `E_pi,Q`, on which `Gamma_0` acts trivially, now gives
`p_E=0`.  This proves (2.1).

The **integral** equality
\(E_\pi^\perp=\pi^*\operatorname {Num}(X)\) in (1.6), not merely the
rational splitting (1.7), gives the unique integral map (2.2). Because the
two integral orthogonal summands are stable, each \(u_i\) restricts to
\(\pi^*\operatorname {Num}(X)\) and descends integrally. Decompose

\[
                    K_Y=\pi^*K_X+K_{Y/X},\qquad K_{Y/X}\in E_\pi. \tag{2.6}
\]

The equality `u_i(K_Y)=K_Y` splits between the two stable summands, so the
descended operator fixes `K_X`.  Finally (1.4), (2.1), and (2.6) give
`K_X\perp T`.  This proves (2.3).  \(\square\)

The proof uses no irreducibility assertion for `T`.  Absence of trivial
quotients for a finite-index finite-covolume subgroup is enough.  It also
uses no bound on `rank(E_pi)` or on its discriminant.

### Theorem 2.2 (finite-orbit minimalization)

The stability premise can be weakened, and the weakened condition is exact.
For the displayed geometric contraction `pi`, the following are equivalent:

1. the numerical `Gamma`-orbit

\[
                  \mathcal E=\{u(E_\pi):u\in\Gamma\}             \tag{2.7}
\]

   is finite;
2. `iota(T)` is orthogonal to `E_pi`; and
3. the trace block descends isometrically through `pi`, that is,
   `iota(T)` is contained in `pi^*Num(X)`.

When these conditions hold, in fact

\[
                         \iota(T)\perp E
                         \quad(E\in\mathcal E).                  \tag{2.8}
\]

The descended image is orthogonal to the minimal canonical class.  Even if
the original `u_i` permute rather than stabilize the exceptional lattice,
the rational isometries

\[
     \operatorname {Res}(s_{r_i})\oplus
             1_{\iota_X(T)^\perp}                               \tag{2.9}
\]

make the target a minimal canonical carrier in the rank-theorem sense.
Consequently its minimal model has Kodaira dimension zero and

\[
                 \boxed{d(n+1)\le20}.                           \tag{2.10}
\]

**Proof.**  Form the negative-definite rational space

\[
                         F=\bigoplus_{E\in\mathcal E}E_\mathbf Q
\]

and the simultaneous orthogonal-projection map

\[
          p_{\mathcal E}:T_\mathbf Q\longrightarrow F,
          \qquad x\longmapsto(p_E(x))_{E\in\mathcal E}.         \tag{2.11}
\]

The group permutes the finitely many summands by integral isometries.  Its
image in `O(F)` is finite: after passing to the kernel of the finite
permutation action, its restriction to each negative-definite integral
summand lies in the finite group `O(E)`.  Equation (2.11) is equivariant for
the action sending the `E`-summand to the `u(E)`-summand.  The same
finite-index/no-invariant-functional argument used for (2.4) forces
`p_{mathcal E}=0`, proving (2.8).  For a chosen contraction, (2.6) then gives
canonical orthogonality on the target.  The rational orthogonal sum (2.9)
supplies the operators required by the minimal canonical-carrier theorem;
that theorem gives (2.10).

This proves `(1)=>(2)`, while `(2)<=>(3)` is exactly the orthogonal
decomposition (1.6).  Conversely, if `(2)` holds, every translate of `E_pi`
lies in the integral lattice

\[
 C=\iota(T)_\mathbf Q^\perp\cap\operatorname {Num}(Y).          \tag{2.12}
\]

The Hodge index theorem gives signature `(1,rho(Y)-1)` on `Num(Y)_Q`, while
`iota(T_Q)` already contains its sole positive direction.  Thus `C` is
negative definite.  The group preserves `C`, and `O(C)` is finite.  Hence
the orbit of `E_pi` is finite, proving `(2)=>(1)`.  \(\square\)

The geometric input is the existence of the initial contraction `pi`; its
other numerical translates need not be effective.  This distinction is
important because an arbitrary integral isometry need not preserve
effectivity, the Mori cone, or contractibility.

### Exact contraction-frame gate

The equivalence is a finite exact test for any *supplied geometric
contraction*.  Choose integral bases, let `G_Y` be the intersection matrix of
`Num(Y)`, let `J` be the inclusion matrix of the full trace lattice, and let
`C_pi` have as columns a full total-transform exceptional frame for the
actual morphism `pi`.  Define

\[
                         D_\pi=J^tG_YC_\pi.                 \tag{2.12a}
\]

Then

\[
 D_\pi=0
 \quad\Longleftrightarrow\quad
 \iota(T)\perp E_\pi
 \quad\Longleftrightarrow\quad
 \text{isometric descent through }\pi
 \quad\Longleftrightarrow\quad
 |\Gamma E_\pi|<\infty.                                   \tag{2.12b}
\]

This is a finite matrix gate **per certified contraction frame**.  It does not prove that all geometric contraction frames form a finite enumerable
set, nor may `C_pi` be replaced by an arbitrary negative-definite frame or
an incomplete list of prime exceptional curves.  Those existence and search
problems remain part of the compiler gap.

The frozen carriers give three exact behaviours.

| carrier model and chosen contraction | `rank(C_pi)` | `rank(D_pi)` | meaning |
|---|---:|---:|---|
| Wehler K3, Bugaenko rank-18 K3, and the four canonical Coble K3 covers; `pi=id` | 0 | 0 | vacuous identity-contraction zero controls |
| `Gamma16` K3 blown up at arbitrary `r` points, then contracted back to the K3 | `r` | 0 | nonminimal positive control; the trace block is the pullback block |
| toric `H6/9 = Bl_6(P2)` with `T_Q=Num_Q` | 6 | 6 | nonzero full-NS defect |
| rational Coble II/I/VI models over `P2` | `10,11,14` | `10,11,14` | nonzero full-NS defects |
| rational Coble MII model `Bl_16(P1xP1)` | 16 | 16 | nonzero full-NS defect |

For the rational rows, `J` is invertible because the trace block is the
whole rational numerical lattice.  Since `G_Y` is nondegenerate and `C_pi`
has independent columns, `D_pi` has full column rank.  These are deliberately
**defect-only negative controls**: the frozen rational models do not certify
integral canonical-preserving extensions of every simple reflection, so the
finite-orbit conclusion of Theorem 2.2 is not invoked for them.  In
particular, no statement about reflection maximality is inferred.  Their
canonical K3 covers are separate minimal zero controls.

The Type-VI and Type-MII rows, including their K3 covers, are in
characteristic three.  Their entries here test only the integral defect
identity; they do not enter the complex `20`-bound (the corresponding
positive-characteristic ceiling is `22`).  The machine replay constructs
the standard `I_(1,r)` or `U+<-1>^r` matrices and computes each rational
`D_pi` literally.  The Type-VI rank-fourteen plane frames and Type-MII
rank-sixteen `P1xP1` frame are pinned independently in the
[Coble blowdown compiler](coble_polytope_first_blowdown_incidence_bridge.md),
rather than inferred from Picard rank alone.

### Corollary 2.3 (nonminimal actual-automorphism carriers)

Suppose in addition that every `u_i` is induced by an actual automorphism of
`Y` and `kappa(Y)>=0`.  Then the minimal model of `Y` is unique.  Every
automorphism descends to it: if \(g\in\operatorname {Aut}(Y)\), then
\(h=\pi g\pi^{-1}\) is a birational self-map of the minimal surface \(X\),
hence is biregular by uniqueness of minimal models in nonnegative Kodaira
dimension, and \(\pi g=h\pi\).  Thus \(g\) preserves the contracted curves
and `E_pi` is automatically `Gamma`-stable.
Theorem 2.1 and the minimal-surface canonical-isometry theorem therefore give

\[
                 \boxed{d(n+1)\le20},\qquad
                 n\le\left\lfloor{20\over d}\right\rfloor-1. \tag{2.13}
\]

Moreover the minimal model has Kodaira dimension zero.  Thus the word
“minimal” is unnecessary in the actual-automorphism theorem on the entire
nonnegative-Kodaira-dimension lane.  Arbitrarily many equivariant preliminary
blowups do not weaken the conclusion.

This does not cover rational or ruled carriers: their minimal models are not
functorial under arbitrary Cremona actions.  Nor does it make effective nef-
wall reflections into automorphisms; the effective-wall obstruction in the
minimal theorem remains in force.

### Corollary 2.4 (rational-carrier mixing diagnostic)

Let `Y` be a rational surface carrying `T` as in (1.2)--(1.3), with `n>=2`.
No contraction to `P2` or a Hirzebruch minimal surface can have a finite
numerical `Gamma`-orbit of its exceptional lattice.  Equivalently, its
exceptional lattice can never be orthogonal to the full trace block.
Otherwise
Theorem 2.2 would inject the rank-`d(n+1)>=3` trace lattice into a minimal
rational numerical lattice of rank at most two.  Thus any such rational
carrier has an **infinite numerical operator orbit of the chosen exceptional
lattice**, unless one of the intrinsic/operator hypotheses already fails.
The translates are not asserted effective or geometrically contractible.

Therefore every proposed upgrade of the universal weak rational
eigencarrier to a canonical-reflection carrier must visibly fail at least one
of the following gates:

1. integral extension of every field reflection to the whole numerical
   lattice;
2. preservation of the canonical class;
3. a finite numerical orbit of one geometric exceptional lattice; or
4. intrinsic injection of the entire trace block.

This is an exact obstruction, not a request to bound the number of blowups.
It explains why the generic `4d(n+1)-3` rational construction cannot simply
be blown down while retaining its full Coxeter action.

## 3. Sharpness: reflection maximality cannot bound the complement

Take a generic exact-Picard K3 surface `X` with

\[
               \operatorname {NS}(X)=S=U+E_8(-1)+E_7(-1).       \tag{3.1}
\]

Its complete square-minus-two chamber is Tumarkin's finite-volume
noncompact `Gamma16` in `H16`.  The project has already certified that this
group is maximal among **all discrete reflection overgroups**, not merely
same-lattice reflection groups.

For any `r>=0`, blow up `r` distinct points:

\[
                          \pi_r:Y_r\longrightarrow X.            \tag{3.2}
\]

In the total-transform basis,

\[
 \operatorname {NS}(Y_r)=\pi_r^*S\perp\langle-1\rangle^r,
 \qquad K_{Y_r}=E_1+\cdots+E_r.                                 \tag{3.3}
\]

Extend every one of the nineteen simple-root reflections by the identity on
the exceptional block:

\[
                          \widetilde s_i=s_i\oplus1_r.            \tag{3.4}
\]

These are integral pairing isometries, preserve the exact original wall
Gram on `pi_r^*S`, and fix `K_{Y_r}`.  They are induced on the numerical
divisor lattice by integral codimension-two algebraic correspondences.
If `Z_i` realizes `s_i` on `X`, then

\[
 \widetilde Z_i={}^{t}\Gamma_{\pi_r}\circ Z_i\circ\Gamma_{\pi_r}
                  -\sum_{j=1}^r E_j\times E_j                 \tag{3.5}
\]

acts as `s_i` on `pi_r^*S` and as the identity on each `E_j`; the minus sign
uses `E_j^2=-1`.  For a K3 root `a_i`, one may take
`Z_i=\Delta_X+a_i\times a_i`.

Here correspondences act covariantly on divisor classes: \(\Gamma_{\pi_r}\)
acts by \(\pi_{r*}\) and its transpose by \(\pi_r^*\).  Thus the first term
in (3.5) sends \(\pi_r^*x\) to \(\pi_r^*s_i(x)\) and kills every \(E_j\),
while

\[
 -(E_j\times E_j)_*(E_j)=-(E_j^2)E_j=E_j.                    \tag{3.5a}
\]

It vanishes on the pullback block and on \(E_\ell\) for \(\ell\ne j\).
This proves the claimed block action and fixes the composition convention;
the certifier replays this exact correspondence matrix rather than checking
only an abstract direct sum.

It follows that

\[
      \operatorname {rank}_{\operatorname {Num}(Y_r)}
             \bigl(\iota(T)^\perp\bigr)=r,
      \qquad K_{Y_r}^2=-r                                      \tag{3.6}
\]

while the raw arithmetic chamber, its nineteen-wall Gram, and its reflection
maximality never change.  Hence there is no function of raw reflection-
maximal Coxeter data that bounds the orthogonal-complement rank **across all
intrinsic correspondence carriers**.

This is not a claim that the full nef cone of `Y_r` has only nineteen walls.
The exact original chamber survives on the orthogonal section
`pi_r^*NS(X)`: a pullback divisor is nef on `Y_r` exactly when its source is
nef on `X`.  The added exceptional directions are deliberately redundant.
That is precisely why a canonical compiler must either contract them
equivariantly or record a nonredundancy condition; reflection maximality of
`P` alone cannot detect them.

## 4. Consequence for the current Gram--Belyi and log-pair routes

The current Gram--Belyi surface stores the Gram in cross-ratios and a marked
finite map; `T` is not an intrinsic Neron--Severi summand.  It therefore
fails before Theorem 2.1 starts.  The universal rational eigencarrier embeds
the trace form, but its frozen theorem does not extend every simple
reflection as an integral canonical-preserving isometry of the whole Picard
lattice, still less prove preservation of a blowdown exceptional lattice.

The correct next compiler gate is consequently finite and geometric:

```text
raw P
  -> intrinsic full trace block and integral simple-reflection operators
  -> construct one certified geometric contraction frame C_pi
  -> compute D_pi=J^t G_Y C_pi and require D_pi=0
     (equivalently: the chosen E_pi has finite numerical operator orbit)
  -> descend automatically by Theorem 2.2
  -> apply d(n+1)<=20 on the minimal carrier.
```

For an actual automorphism carrier of nonnegative Kodaira dimension, the
middle test is automatic and the orbit is a singleton.  For arbitrary
correspondences it is substantive:
canonical-class preservation alone does not preserve effectivity or a Mori
contraction.  In the log-pair language the missing property is preservation
of the extremal exceptional face selected by the MMP, not a bound on the
number of its rays.

This advances the global theorem in one precise direction and closes another:

- **proved:** a finite numerical orbit of one geometric exceptional lattice
  forces exact full-trace descent
  and yields the existing Picard ceiling;
- **proved impossible:** reflection maximality alone bounds every
  nonminimal complement;
- **still open:** raw maximal Coxeter data force an intrinsic carrier with a
  chosen geometric exceptional lattice having finite numerical operator
  orbit, or directly force a minimal `K`-trivial carrier.

## 5. Machine replay

The symbolic descent proof is the argument in Section 2.  The exact
`Gamma16` sharpness family is frozen in
[`equivariant_minimalization_trace_descent.json`](../data/enumerations/equivariant_minimalization_trace_descent.json).
Its certifier imports the archived nineteen-wall Gram, reconstructs the
rank-seventeen determinant-two wall basis, and for sample values
`r=0,1,2,5,11` checks the block reflection identities, canonical vector,
orthogonality, signatures, and unchanged trace Gram:

```text
python3 -m scripts.certify_equivariant_minimalization_trace_descent --check-archive
python3 -m pytest -q tests/test_equivariant_minimalization_trace_descent.py
```

The proof is uniform for every `r>=0`; the samples are regression checks, not
the reason unboundedness follows.

Companion project results:

- [minimal canonical-carrier theorem](project_native_minimal_surface_automorphism_dimension_bound.md);
- [canonical-orthogonal K3-cover calibration](canonical_orthogonal_k3_cover_dimension_theorem.md);
- [`Gamma16` reflection-maximality theorem](wang_h16_maximal_reflection_no_correction.md);
- [universal rational eigencarrier](rational_surface_eigencarrier_universality.md); and
- [Gram--Belyi obstruction theorem](gram_belyi_carrier_dimension_field_bound_obstruction.md).

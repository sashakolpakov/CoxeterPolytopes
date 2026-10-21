# Independent hostile audit: exact AGCat return and global maximal arithmetic Coxeter finiteness

**Date:** 2026-08-12  
**Verdict:** **PASS after four scope/interface repairs.**

I audited the terminal composite theorem in
[`agcat_global_maximal_arithmetic_coxeter_finiteness.md`](agcat_global_maximal_arithmetic_coxeter_finiteness.md),
not merely its executable arithmetic samples.  The repaired theorem proves,
for hyperbolic dimensions (n\geq2), both:

1. an exact faithful functor from marked finite-volume arithmetic Coxeter
   chambers to a real-framed algebraic quadric--hyperplane arrangement
   category, with exact return on its image; and
2. finiteness of unmarked normalized reflection-maximal arithmetic Coxeter
   chambers and their side-reflection groups.

The proof is genuinely **AG-mediated but not AG-only**.  Integral regular
functions on the exact arrangement carrier supply the essential global
finite-field theorem.  The terminal passage also uses declared classical
hyperbolic geometry and arithmetic of orders.  It does not invoke the
Agol--Belolipetsky--Storm--Whyte finiteness conclusion, Nikulin's
fixed-field/fixed-rank finiteness conclusion, or a pre-existing global
ground-field theorem.

The owner initially had four defects visible only at the composition level.
All were repaired before the frozen bytes audited here:

1. the two global theorems now explicitly say (n\geq2), matching the
   rank-(m\geq3) native-host theorem;
2. the source is the groupoid of all marked arithmetic chambers, while the
   terminal theorem restricts to its reflection-maximal full subgroupoid;
   the target now specifies its base fields and morphisms, so “faithful
   functor” is literal rather than objectwise shorthand;
3. the fixed-((k,m)) theorem now says *admissible arithmetic Vinberg
   hosts*, exposing the nonselected definiteness used in its Minkowski
   bounds; and
4. Vinberg's structural lemma reducing arithmetic reflection groups to the
   quadratic/lattice-form model is now explicit in the import ledger.  It is
   a structural result, not a finiteness theorem.

No mathematical defect remains in the repaired composition.

## 1. Frozen theorem and dependency pins

| artifact | SHA-256 |
|---|---|
| composite owner theorem | `441deb5a379fc59c6b4c1e1e86212f5f84db4e603a4b37259200db43a4508dfd` |
| dimension-two orbicurve owner | `68ab8d85002106fde4493540ee6c6409247211cdda5364a3181000ae0135b659` |
| dimension-two hostile audit | `7be0ea8da1eb6dcc97ced1b68f40a47c9eb3ae382b47ab054e77735df4afcfc9` |
| dimension-(\geq6) field owner | `05acee48bd7a6c03e92febc32434cc38bcd3581c0e9044dafec3f07f9dbab678` |
| dimension-(\geq6) hostile audit | `e5e2754bf886d922ba50173e8afef67ff7eceffc62b706bb3c4dea4a309cac69` |
| dimension-(4,5) descent owner | `8157c0779705f5f265e75126d10361441bf2e2b504e522d917cbda848b75f019` |
| dimension-(4,5) independent audit | `74c32b41b593419951ecb6db0c43c3364134cc033e77038175156bffb385b5da` |
| compact-(\mathbf H^3) fork owner | `4a971c64460a7d234114ecb12b793e613b080945850098da37b903d2c50876d5` |
| compact-(\mathbf H^3) independent audit | `7bded0c2dc46c1652e06dd15ab691685736db4b327c6fb8a218ca84066cef996` |
| native-host owner | `7ba392dd833358379ea92d1385909a52b2dc65417a9cb4ccb856170118df02a3` |
| native-host independent audit | `0b05127f62a89b16060d9cbaacde93009a56ac056e08bb7e16e0b8f79322c3b1` |
| independent composition certifier | `8acee124164a974157ea19a32b5a8d765ad2252f15078fd9882909089bf189c7` |
| independent composition tests | `b534548ed8816e271f9a623666d756f7505e3ecb17d59730ae0100fc02564785` |

The certifier hard-pins the first eleven mathematical artifacts.  The
terminal synthesis is intentionally not hard-pinned there: it receives the
audit backlink only after this report is frozen.  The audited mathematical
owner is already immutable under its displayed pin.

## 2. The functor and return are genuinely categorical

The repaired source

\[
                  \mathsf{Cox}^{\mathrm{mark}}_{\mathrm{arith}}
\]

contains finite-sided finite-volume arithmetic Coxeter chambers with their
displayed minimal Vinberg field, distinguished real place, ordered inward
facet-root rays, and selected positive component.  Its morphisms are marked
field-compatible projective similitudes.

The target objects are, over a displayed totally real field (k),

\[
          (\mathbf P(V),Q,H_1,\ldots,H_f,\sigma_0,C),       \tag{2.1}
\]

with ordered polar hyperplanes and finite real framing.  A target morphism
is a field-isomorphism-compatible projective isomorphism preserving the
quadric, ordered hyperplanes, distinguished place, and oriented component.
Restriction to the recovered half-space intersection defines the return on
morphisms.

For (n\geq2), the objectwise inverse has no ambiguity:

* a smooth quadric (Q\subset\mathbf P(V)) determines its nondegenerate
  quadratic form up to a common scalar;
* quadric polarity sends (H_i=\mathbf P(r_i^\perp)) back to the line
  (kr_i);
* the real framing chooses the inward ray and the correct hyperbolic
  component; and
* normalized pairings recover finite angles, parallelism, ultraparallel
  distance, and therefore the identical half-space intersection.

Consequently

\[
                  R_{\rm arr}F_{\rm arr}(P)\cong P,       \tag{2.2}
\]

and the same restriction verifies functoriality on morphisms.  Two source
morphisms inducing the same marked projective map are already the same
projective similitude, so the functor is faithful.  Forgetting the ordering
is only the objectwise finite quotient by facet permutations.

The finite real decoration is essential and honestly disclosed.  This is a
genuine algebraic arrangement category with finite real/semialgebraic
framing, not a claim that a bare complex quadric remembers oriented
half-spaces.

## 3. Arithmeticity lives on the exact arrangement carrier

Let

\[
 a_{ij}=\frac{2q(r_i,r_j)}{q(r_i,r_i)},\qquad
 b_C=a_{i_1i_2}\cdots a_{i_si_1}.                         \tag{3.1}
\]

The (b_C) are regular invariant functions on the diagonal root-rescaling
quotient.  If (T_i=1-s_i), right-to-left rank-one composition gives

\[
          b_C=\operatorname {tr}(T_{i_1}\cdots T_{i_s})
              \in\mathcal O_k.                           \tag{3.2}
\]

The projective-module trace is local and glues, so (3.2) does not assume a
free root module or integral individual directed entries.  This is the
precise bridge from the AG carrier to arithmetic: Vinberg integrality and
the selected/nonselected real-place inequalities cut out integral points in
the arrangement quotient.

The compact-(\mathbf H^3) coordinate illustrates why the carrier is not
decorative.  In invariant cycle coordinates,

\[
       \beta=t_{12}+t_{13}+t_{23}-z-4\in\mathcal O_k.      \tag{3.3}
\]

For the auxiliary symmetric Gram

\[
 G=\begin{pmatrix}-2&u&a\\u&-2&b\\a&b&-2\end{pmatrix},
\]

one has (\beta=\det(G)/2=u^2+a^2+b^2+uab-4).  The sign in
(3.3) is forced because (z=-uab) in this convention.  The independent
certificate rechecks this identity over exact rational samples.  The
selected signature and nonselected definiteness give

\[
             0<\beta_0<256,\qquad \beta^\tau<0,
             \qquad k=\mathbf Q(\beta),                 \tag{3.4}
\]

which is the primitive full-field coordinate required by relative Fekete.

## 4. The global finite-field composition has no missing dimension

The branch partition is exhaustive for every (n\geq2):

| input | field mechanism | terminal output |
|---|---|---|
| noncompact, every (n\geq2) | a (k)-rational cusp vector is isotropic at every embedding, impossible at a nonselected definite place | (k=\mathbf Q) |
| compact (n\geq6) | narrow facet, face-average/Hodge forcing, finite transition arrangement charts, Fekete and polygon branches | finite set (\mathscr K_{\geq6}) |
| compact (n=4,5) | exact codimension-one/two descent or finite transition packet | reduces to (n=3,2) or a finite set |
| compact (n=3) | all-perpendicular orbicurve branch, tail, small packet, or arbitrary-label general fork | finite set (\mathscr K_3) |
| (n=2) | exact real-DM genus-zero orbifold log curve and Long--Maclachlan--Reid | finite set (\mathscr K_2) |

The iteration in dimensions four and five is in the correct order.  A
dimension-five descendant lands in dimension four or three; a dimension-four
descendant lands in dimension three or two.  Every orthogonal descent is
literal over its (k)-rational reflection eigenspaces and inherits the full
cycle field.  A failed literal descent produces one of the audited
four-node transition packets.  Nothing assumes a facet or dimension bound.

The compact-three-dimensional fork audit independently checks:

\[
 -\frac{s_m^\rho s_n^\rho}{4}\leq\beta^\tau<0,
 \qquad E=\mathbf Q(s_m,s_n),                              \tag{4.1}
\]

the correction for the unique exceptional selected place in the norm, and
the resulting uniform (m,n<10^6) cutoff.  Fixed-label Fekete finiteness
then applies.  The global incidence split includes label three directly and
does not send it through an inapplicable tail chart.

The dimension-two input is exactly scoped.  The real orbifold curve is an
algebraic carrier, while Long--Maclachlan--Reid supplies arithmetic--spectral
finiteness for congruence genus-zero Fuchsian groups.  Passing a polygon
group to a maximal arithmetic reflection overgroup preserves its field;
therefore the branch supplies finite **fields** even when the original
polygon was not reflection-maximal.

Thus

\[
  \boxed{\{k(P):P\text{ finite-volume arithmetic Coxeter},\ n\geq2\}
         \text{ is finite}.}                            \tag{4.2}
\]

Reflection maximality is not used in (4.2).

## 5. The fixed-field/fixed-rank native-host reconstruction passes

For fixed totally real (k) and (m\geq3), in the admissible Vinberg
model, the downstream proof has the following independently verified chain.

1. A normal core and finite translate sum give a full projective
   (W(P))-stable lattice (L).  Hence
   
   \[
                  W(P)\subseteq R(L)\subseteq O^+(q,L)     \tag{5.1}
   \]
   
   with finite index, and (R(L)) is a finite-covolume arithmetic
   reflection group.
2. The declared Nikulin geometric lemma selects an (m)-wall basis with
   connected nonorthogonality graph and normalized two-cycles below
   (62^2).
3. Closed projective traces are algebraic integers.  At the selected place
   their absolute values are below (62^m); at every other place they are
   at most (2^m).  Minkowski discreteness therefore leaves finitely many
   complete simple-cycle packets.
4. A spanning-tree gauge shows that a complete packet, including every
   two-cycle, determines the Cartan representation up to diagonal
   (k)-conjugacy.  Connected spanning normals prove absolute
   irreducibility, so Burnside gives
   (k[s_1,\ldots,s_m]=\operatorname {End}_k(V)), and the invariant form is
   unique up to scale.
5. For one stable projective lattice (L_0), the order
   
   \[
     \Lambda=\mathcal O_k[s_1,\ldots,s_m]
       \subset B=\operatorname {End}_{\mathcal O_k}(L_0)
   \]
   
   is full.  A conductor ideal (\mathfrak c\ne0) gives, for every stable
   full lattice,
   
   \[
                         \mathfrak cM\subseteq L\subseteq M,
                         \qquad M=BL.                     \tag{5.2}
   \]
   
   Matrix units force (M=\mathfrak aL_0); only the finite ideal class
   group remains up to homothety, and (M/\mathfrak cM) has finitely many
   submodules.

This proves finite native hosts for fixed ((k,m)) without individual
Cartan-entry integrality, freeness, principal ideals, or a hidden unit
bound.  It reconstructs and simplifies the arithmetic core behind
Nikulin's fixed-data method, as the owner explicitly says; it does not cite
that theorem's finiteness conclusion.

## 6. Dimension and maximality close in the right quantifier order

The Vinberg compact and Khovanskii--Prokhorov finite-volume theorems give
the independent universal geometric consequence

\[
                              2\leq n\leq995.             \tag{6.1}
\]

The noncompact source is often stated as nonexistence for (n>995); this
weaker consequence is all the proof uses.  It is independent of
arithmeticity, a finite field list, and fixed-data reflective-lattice
finiteness.

Combining (4.2) and (6.1) leaves only

\[
              (k,m)\in\mathscr K_{\rm Cox}\times
                         \{3,\ldots,996\},               \tag{6.2}
\]

exactly (996-3+1=994) possible ranks.  The fixed-((k,m)) theorem is
applied only after this finite product is formed.  There is no fixed-
dimension-to-uniform-dimension quantifier leak and no assumed facet bound.

For every original chamber, (5.1) holds.  Absolute reflection maximality
gives (W(P)=R(L)).  Maximality only among arithmetic reflection groups is
also enough because (R(L)) is a finite-index arithmetic reflection
overgroup.  More generally, any discrete overgroup of a finite-covolume
arithmetic lattice has finite index and is commensurable with it, so it is
arithmetic; the two maximality scopes agree here.

Finally, a generic compact path through a locally finite complete mirror
arrangement crosses finitely many mirrors.  Successive wall reflections give
a finite gallery, so a full host has one unmarked chamber orbit.  Finite
fields, finite ranks, finite hosts, and maximality therefore give the
claimed finite set of normalized chambers.  Exact return then pushes the
conclusion back from the AG carrier to the original Coxeter objects.

## 7. Dependency firewall and exact AG attribution

The following imports are real and correctly disclosed:

* Vinberg's structural reflection lemma, reducing an arithmetic reflection
  input to an admissible quadratic/lattice form;
* the genus-zero arithmetic Fuchsian theorem of
  Long--Maclachlan--Reid;
* relative Fekete theory and elementary cyclotomic norms;
* finite spherical/Lanner and four-node support calculations;
* the Vinberg--Khovanskii--Prokhorov universal dimension ceiling; and
* Nikulin's geometric connected-frame lemma with constant (62).

The following are **not** used:

* the Agol--Belolipetsky--Storm--Whyte global finiteness theorem;
* Nikulin's fixed-field/fixed-rank reflective-lattice finiteness conclusion;
* Nikulin's pre-existing global ground-field theorem or its published
  transition degree bound;
* Mostow--Prasad rigidity as a chamber-counting substitute; or
* BAB, KSBA, K-moduli, or `AG-DEG` as a claimed uniform facet-degree theorem.

The attribution “AG-mediated, not AG-only” is exact.  Algebraic geometry
enters through the real-framed quadric arrangement, its invariant quotient
and integral regular functions, the face/Hodge input in the high-dimensional
field funnel, and the real-DM orbifold-curve branch.  Arithmetic Fekete
theory proves integral-point finiteness on the contracting arrangement
charts.  Hyperbolic Coxeter geometry supplies dimension and the bounded
frame; order theory supplies the terminal native-host count.

The stronger direct moduli implication

\[
  \text{reflection maximality}\Longrightarrow
  \text{uniformly bounded Hilbert/log facet degree}       \tag{7.1}
\]

remains open and is not needed.  The theorem therefore answers the actual
functor-and-pushback question without pretending that standard surface
boundedness alone proves the result.

## 8. Executable and document checks

The independent composition certificate is
[`certify_agcat_global_maximal_arithmetic_coxeter_finiteness_independent_audit.py`](../scripts/certify_agcat_global_maximal_arithmetic_coxeter_finiteness_independent_audit.py),
with tests in
[`test_agcat_global_maximal_arithmetic_coxeter_finiteness_independent_audit.py`](../tests/test_agcat_global_maximal_arithmetic_coxeter_finiteness_independent_audit.py).
It checks the frozen dependency pins, category/scope interface, the complete
field-branch partition, (n\mapsto m=n+1), exact fork determinant/cycle
identities, local links, and TeX delimiter balance.  It returns
`all_checks_pass: true`.

The owner and independent arithmetic regression files complete with

```text
31 passed
```

and both the composite owner and this audit parse successfully with Pandoc.
These machine checks support the exact interfaces and finite calculations;
the mathematical proof is the reconstruction above and in the frozen owner
and hostile reports.

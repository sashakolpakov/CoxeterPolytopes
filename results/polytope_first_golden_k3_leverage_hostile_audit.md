# Hostile audit: polytope-first golden K3 leverage

**Verdict (2026-08-08): pass with strict scope.**  No fatal circularity was
found in the fixed-realization theorem.  The trace-overclosure, shell,
primitive-embedding, period, Weyl-frame, and compact hidden-wall arguments do
give a genuine

\[
P\longrightarrow \text{K3 divisor lattice}\longrightarrow P
\]

theorem for the stated full-trace golden species.  They do **not** give a
nonexistence theorem for compact arithmetic Coxeter polytopes in general.
The only unconditional upper bound presently obtained is on full-trace K3
**carriers**, not on polytopes.

The audited report is
[`polytope_first_golden_k3_leverage.md`](polytope_first_golden_k3_leverage.md).

## Claim-by-claim verdict

| Claim | Verdict | Exact scope or required qualification |
|---|---|---|
| finite trace overclosures | pass | fixed normalized wall realization only |
| finite compact hidden-root search | pass | all five trace-root shells must be retained for a merely trace-even overclosure |
| golden visible-shell collapse | pass | collapse to norm two and the paired-root statement require the field lattice to be even |
| Nikulin stable inequality | pass | sufficient, not necessary; applies through `H8` under the displayed length bound |
| K3 period/Weyl return | pass | returns the same complete restricted chamber; supplies no surface automorphism |
| paired roots are irreducible curves | pass | at a generic facet of a field-even, complete chamber |
| `H9` binary-complement gate is finite | pass | for each fixed rank-20 trace lattice `N` |
| no carrier from `H10` upward | pass | full-trace quadratic K3 carriers only |
| low-dimensional Lanner class is non-singleton | pass | existential full-root subchambers; the displayed seed simplex need not itself be root-complete |
| universal polytope bound | not proved | requires Missing Lemma U, including exhaustion of all normal scalings |

## 1. Finite overclosures and embedding orbits

Let `R` be the full-rank `O=Z[phi]` span of the facet roots.  Its trace form is
nondegenerate and integral.  Any same-realization integral full-trace carrier
containing those roots satisfies

\[
N_R\subseteq N'\subseteq N_R^\vee.
\]

Thus the possible additive overgroups are subgroups of the finite
discriminant group.  `O`-stability, trace evenness, field integrality, and
stability under the displayed reflections are finite exact tests.  This
proves fixed-Gram finiteness.

It does not enumerate alternative positive normal scalings of an abstract
Coxeter diagram, unequal-root realizations, non-full-Picard carriers, or
carriers outside this restriction-of-scalars construction.

The report was corrected during this audit on one important point.  Nikulin's
local discriminant criterion decides existence and complement genera; it does
not by itself enumerate literal primitive embeddings.  Embedding orbits
require complement isometry classes and the finite discriminant
anti-isometry double cosets, modulo the appropriate images of `O(N)` and
`O(T)`.  With that quotient explicitly stated, the finiteness claim is sound.

## 2. Hidden-root finiteness and the golden shell

For a visible K3 root `delta^2=-2`, writing
`h(delta,delta)=a+b phi` gives

\[
\operatorname{Tr}\!\left(-\frac{\phi}{\sqrt5}(a+b\phi)\right)
=-(a+b)=-2.
\]

Total positivity at the two real places leaves exactly

\[
1+\phi,\quad 2,\quad 3-\phi,\quad 4-2\phi,\quad 5-3\phi.
\]

If `h` is field-even, only `2` and `4-2phi` occur.  Since
`u=phi-1` and `u^2=2-phi`, multiplication by `u` pairs the two shells and
does not change their restricted hyperplane.  Comparing norms proves that a
restricted mirror has exactly the four signed representatives

\[
\pm r,\qquad \pm ur.
\]

For a merely trace-even overclosure the five-shell statement remains valid,
but the two-shell collapse must not be asserted.

If a root hyperplane meets compact `P`, its identity component lies on a
fixed-radius sphere in the positive-definite orthogonal complement of a point
of `P`.  Compactness bounds these spheres uniformly.  At the conjugate place
the form is positive definite and each of the five shell norms is fixed.
The Minkowski lattice is discrete, so only finitely many roots can meet
`P`.  This proves termination of the hidden-wall test.  The proof is
effective in principle; a general-purpose implementation of this compiler is
still listed as future work.

## 3. Stable K3 embedding and exact return

The trace lattice has

\[
\operatorname{rk}N=2(n+1),\qquad
\operatorname{sig}N=(1,2n+1),
\]

so a primitive K3 complement has signature `(2,18-2n)` and rank `20-2n`.
Nikulin's stable theorem applies when

\[
20-2n\ge \ell(A_N)+2.
\]

Consequently:

- for `n<=4`, `ell(A_N)<=rk(N)=2n+2` makes the inequality automatic;
- for `5<=n<=8`, the sufficient condition is
  `ell(A_N)<=18-2n`;
- Bugaenko at `n=8` has `ell(A_N)=2`, exactly the stable boundary.

The complement has both positive and negative directions throughout this
stable range, as required.  This is the standard sufficient range of
[Nikulin's primitive-embedding theorem](https://doi.org/10.1070/IM1980v014n01ABEH001060).

Given a primitive embedding, period surjectivity supplies a projective K3
with exact Picard lattice `N`: choose a generic period in the complement when
its period domain has positive dimension; in rank 20 the positive binary
complement is already the whole period plane.  If the complete restricted
root arrangement has `P` as its selected component, equality of all root
signs gives

\[
D\cap E_\phi=P
\]

for the ambient root chamber `D`.  The K3 Weyl group acts transitively on
root chambers, so a Weyl change of marking sends `D` to the ample chamber.
This proves the same-polytope equality with the nef eigensection.  It does
not claim that multiplication by `phi` is a Hodge isometry or a surface
automorphism.  The period and Weyl facts are also summarized in
[Huybrechts' K3 lectures](https://ncatlab.org/nlab/files/HuybrechtsLecturesOnK3.pdf).

## 4. Paired-root irreducibility

This step is valid, with the field-even qualifier now explicit.  At a generic
point of a restricted facet, shell uniqueness leaves only
`+/-r,+/-ur`.  Exact trace arithmetic gives

\[
r^2=(ur)^2=-2,\qquad (r,ur)=0.
\]

Hence the local ambient root system is `A1+A1`.  The ambient chamber chooses
one sign of each root as a simple wall.  After Weyl framing, simple walls of
the K3 ample chamber are classes of irreducible smooth rational curves.
Their intersection is zero, so the two distinct curves are disjoint.  This
argument would not justify irreducibility for an arbitrary effective
`(-2)`-class; it works because the two classes are simple chamber walls.

## 5. The `H9` and `H10` boundaries

At `n=9`, `N` has signature `(1,19)` and a primitive embedding in the
unimodular K3 lattice exists exactly when there is an even positive-definite
binary lattice `T` with a discriminant anti-isometry

\[
q_T\simeq-q_N.
\]

For fixed `N`, `det(T)=|det(N)|`; reduced positive binary forms of fixed
determinant and all anti-isometries of their finite discriminant groups form
a finite list.  Gluing produces the unique even unimodular lattice of
signature `(3,19)`.  An orientation of `T` gives the corresponding singular
K3 surface, in accord with the
[Shioda--Inose classification](https://doi.org/10.1017/CBO9780511569197.010).
The gluing data may distinguish marked carriers even when the underlying
unmarked surface is the same.

For `n>=10`, the trace rank is `2(n+1)>20`, so it cannot be the full
Neron--Severi lattice of a complex K3 surface.  More generally a degree-`e`
full-trace K3 carrier obeys

\[
e(n+1)\le20.
\]

This is a carrier obstruction.  It says nothing about a polytope that has no
carrier of this species.

## 6. The non-singleton Lanner class

The eight displayed equal-root seed Grams were replayed exactly.  Each has
identity inertia `(n,1)`, positive-definite conjugate, and positive-definite
vertex deletions; hence it defines a compact arithmetic Coxeter simplex.
The trace lattices are even, have signature `(1,2n+1)`, and satisfy the
stable embedding inequality in dimensions two, three, and four.

The passage from a seed simplex `S` to a complete chamber is sound:

1. choose a regular point in `S` and let `C` be its chamber in the full
   norm-two arrangement;
2. every seed wall is itself in that arrangement, so the signs at the chosen
   point force `C` to lie inside `S`;
3. therefore the closure of `C` is compact;
4. only finitely many mirrors meet compact `S`, so `C` has finitely many
   facets; and
5. by construction `C` is root-complete.

This gives genuine surface-free inputs and K3 outputs in three different
dimensions, so the class is nonempty and non-singleton.  It is an
existential result: the audit does not prove that each displayed simplex is
already its own full-root chamber, does not list the final walls of `C`, and
does not claim worldwide-new abstract polytopes.

## 7. Exact quantifier map

| Population | What is proved |
|---|---|
| fixed compact field-even golden realization | all trace overclosures and all K3 roots meeting `P` form finite searches |
| root-complete member, `n<=4` | a K3 carrier returning exactly `P` exists automatically |
| root-complete member, `5<=n<=8`, `ell(A_N)<=18-2n` | the same exact return exists |
| fixed rank-20 trace lattice in `H9` | K3 association is a finite binary-complement decision |
| full-trace degree-2 K3 carrier | necessarily `n<=9` |
| compact arithmetic Coxeter polytopes in general | no new dimension bound follows |

The class is surface-free because its base predicates—fixed oriented Gram,
integral `O`-lattice, admissible signatures, compactness, and
root-completeness—do not mention a surface.  In the stable ranges, K3
existence follows from a numerical lattice theorem.  Outside those ranges,
primitive embeddability is itself the exact arithmetic carrier gate; treating
it as a definition would give no independent polytope leverage.

Missing Lemma U is therefore genuinely missing.  In particular, its clause
4 quantifies over **all** admissible normal scalings, whereas the proven
compiler fixes one normalized realization.  The finite saturation search is
not an all-scaling classification.  Only after a surface-free universality
theorem supplies a carrier for every polytope in an independently defined
class can the Picard ceiling be transferred into a polytope ceiling.

## Replay

The target artifacts replayed successfully after their quantifier fixes:

```sh
PYTHONPATH=.:src python3 \
  scripts/certify_polytope_first_golden_k3_leverage.py --check-archive
PYTHONPATH=.:src pytest -q \
  tests/test_polytope_first_golden_k3_leverage.py
```

Result: certificate
`cb299551f61bb11c86c2dbb12a3b2a02edff3195566b0a9168615e2e98dd3e1e`;
`6 passed`.

The certificate checks the eight seed Grams, trace Smith data, stable
inequalities, paired-root arithmetic, and the frozen Bugaenko/silver boundary
inputs.  It does not implement the report's general overclosure and hidden-
root compiler; those parts rest on the finite lattice arguments audited
above.

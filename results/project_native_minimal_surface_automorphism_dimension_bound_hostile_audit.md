# Hostile audit: minimal-surface canonical-isometry dimension bound

**Verdict (2026-08-10): PASS after strengthening and scope repairs.** The
conditional inequality

\[
             [k:\mathbf Q](n+1)\le 20
\]

is correct over \(\mathbf C\). The proof needs less than the original draft
claimed: each restriction-of-scalars simple reflection need only extend to a
pairing isometry of \(\operatorname {Num}(X)_{\mathbf Q}\) fixing \(K_X\).
Actual automorphisms are required only for the effective- and ample-cone
statements. No universal polytope-to-surface theorem follows.

Audited theorem:
[project_native_minimal_surface_automorphism_dimension_bound.md](project_native_minimal_surface_automorphism_dimension_bound.md).

## Checks

1. **Restriction of scalars: pass.** With one real signature \((1,n)\) and
   \(d-1\) negative-definite conjugates, the trace form has rank
   \(d(n+1)\) and signature \((1,d(n+1)-1)\). The \((-1)\)-eigenspace of
   \(\operatorname {Res}(s_r)\) is the underlying rational space of \(kr\),
   hence has dimension \(d\), not one. Since the root lines span \(V\),
   canonical preservation forces \(K_X\perp T\).

2. **Hodge-index and classification step: pass.** The trace block consumes
   the unique positive numerical direction, so its orthogonal complement is
   negative definite. On a minimal surface with \(\kappa\ge0\), \(K_X\) is
   nef and \(K_X^2\ge0\); canonical orthogonality therefore gives
   \(K_X\equiv0\). The complex Kodaira-zero list and Picard ceilings are
   K3/20, Enriques/10, abelian/4, and bielliptic/2. A minimal ruled surface
   has numerical rank at most two, excluding \(n\ge2\).

3. **Positive characteristic: pass at 22, not 20.** The same numerical
   argument and the Bombieri--Mumford Kodaira-zero classification apply.
   K3 surfaces satisfy \(\rho\le b_2=22\), and rank 22 occurs; abelian
   surfaces have \(\rho\le b_2=6\), Enriques surfaces have numerical rank
   10, and (quasi-)hyperelliptic cases are smaller. A bound of 20 requires a
   compatible lift and is not asserted. Primary classification source:
   Bombieri--Mumford,
   [*Enriques' classification of surfaces in char. p, II*](https://www.dam.brown.edu/people/mumford/alg_geom/papers/1976d--EnrClass-II-NC.pdf)
   and
   [III](https://doi.org/10.1007/BF01390138).

4. **Pseudo-automorphisms: pass after proof cleanup.** Factorization of a
   birational map between smooth projective surfaces into point blowups and
   blowdowns shows that a non-biregular map or its inverse contracts a
   curve. Thus an isomorphism in codimension one is an automorphism. The
   revised proof does not rely on the overcompressed assertion that every
   common exceptional curve can individually be removed from a graph.

5. **Effective-wall obstruction: pass with a needed integrality/effectivity
   clarification.** If an automorphism sends a negative rational divisor
   class \(r\) to \(-r\), one sign being rationally effective would make both
   signs effective, contradicting positivity against an ample class. On a
   K3, Riemann--Roch makes one sign of every integral \((-2)\)-class
   effective. This conclusion does not apply to abstract canonical-preserving
   isometries or monodromy.

6. **Projected ample cone: pass only in the revised form.** Orthogonal
   projection to the distinguished Lorentzian factor sends every ample class
   into one positive component because all discarded summands are negative
   definite. For an actual-automorphism realization of a finite-covolume
   Coxeter group, the projected cone is nonempty, open, convex, and invariant.
   A lattice has full limit set; in the Klein model a proper convex invariant
   subset would be separated from part of the convex hull of that limit set.
   Hence the projection is the whole distinguished positive cone. The same
   assertion in the full rational trace hyperbolic space is false for
   \(d>1\), and the report does not make it.

7. **Wehler calibration: exact pass.** For a very general \((2,2,2)\) K3 in
   \((\mathbf P^1)^3\), Hashimoto--Lee give the displayed Picard Gram, the
   three deck involutions, and
   \(\operatorname {Aut}(X)=C_2*C_2*C_2\cong\Pi(2)\), where \(\Pi(2)\) has
   index six in \(PGL_2(\mathbf Z)\). The double-cover divisor relation gives
   \(H_i\mapsto-H_i+2H_j+2H_k\), the reflection in
   \(r_i=H_j+H_k-H_i\). Exact replay gives \(r_i^2=-4\),
   \((r_i,r_j)=4\), and ideal vertices \(\mathbf R_{>0}H_i\). Also every
   lattice square is divisible by four, so there is no \((-2)\)-vector and
   the ample cone is the positive cone. Primary source:
   Hashimoto--Lee,
   [*Free automorphism groups of K3 surfaces with Picard number 3*](https://arxiv.org/abs/2306.16147),
   Sections 1 and 5. The audited arXiv source payload had SHA-256
   `0896dc3f2e9e7f29ac30a38ecb400703ef562281faf295881658b16c85eb6085`.

8. **Auxiliary controls: pass after deleting an invalid full-lattice
   inference.** Martin's
   primary source states that a bielliptic Enriques involution acts as
   \({\rm id}_U\oplus(-{\rm id}_{E_8})\). That proves only that its action on
   the full numerical lattice is not a rational hyperplane reflection. It
   does **not** exclude a smaller invariant \(T\) on which the restriction is
   a trace reflection, because the action on \(T^\perp\) is unconstrained.
   The report now states this as an unproved calibration lane, not a no-go.
   See
   [*Automorphisms of unnodal Enriques surfaces*](https://arxiv.org/abs/1908.00049),
   Example 2.2 and Theorem A. Totaro's Hesse action is holomorphic and hence
   orientation-preserving on its full \(\mathbb H^3\) factor, so that full
   factor is not a hyperplane-reflection representation. This likewise does
   not exclude an unidentified smaller invariant trace block.

## Frozen replay

- Exact certifier:
  [certify_project_native_minimal_surface_bound.py](../scripts/certify_project_native_minimal_surface_bound.py)
- Tests:
  [test_project_native_minimal_surface_bound.py](../tests/test_project_native_minimal_surface_bound.py)
- Replay result: `5 passed`.

The certifier checks the trace signature/rank formula, the Wehler lattice
signature, all three integral reflection matrices, their level-two
congruence, the root Gram, the three ideal vertices, and the divisibility
obstruction to a \((-2)\)-vector.

## Exact remaining gap

The theorem is conditional on a full faithful trace-module injection into
one **minimal** carrier and on canonical-preserving extensions of all simple
reflections. Neither raw Coxeter maximality nor a local dessin packet proves
that existence statement. The result is therefore a valid project-native
dimension bound for the canonical-carrier lane, not yet an ABSW-type theorem
for all arithmetic Coxeter polytopes.

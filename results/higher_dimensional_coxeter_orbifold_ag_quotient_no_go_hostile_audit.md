# Hostile audit: higher-dimensional quotient-orbifold AG no-go

**Audited owner SHA-256:**
`24d34f3857417da09a531fee13675f863f2fc339e1e473b5a9a6137042648ff7`.

**Verdict: PASS for the stated direct-quotient no-go.**  The compact
\(\mathbf H^6\) control, the smooth-proper-DM obstruction, and the direct
nonzero-degree-cover obstruction are sound.  Nothing here proves a no-go for
singular AG carriers, real algebraic models, or external decorated carriers.

## Checks

1. **The source can be maximalized without reflection finiteness.**  The
   archived \(P_{6,10}\) is the index-two double of Bugaenko's compact
   arithmetic \(\mathbf H^6\) chamber, hence its reflection group \(W_0\) is
   a cocompact arithmetic lattice.  If a discrete reflection group \(R\)
   contains \(W_0\), a fundamental region for \(R\) sits inside one for
   \(W_0\); discreteness gives positive covolume, so \(R\) is a lattice and
   \([R:W_0]=\operatorname{vol}(W_0\backslash\mathbf H^6)/
   \operatorname{vol}(R\backslash\mathbf H^6)\).  The fixed-dimensional
   Kazhdan--Margulis lower bound makes these indices bounded integers.
   Choosing their maximum really excludes every strict discrete reflection
   overgroup.  Finite extension preserves cocompactness and arithmeticity.

2. **The DM-cover passage is valid.**  For
   \(\Gamma_{\max}=W_{\max}^+\), Selberg's lemma and passage to the normal
   core give a torsion-free normal finite-index subgroup \(\Gamma_0\).  The
   corresponding representable finite etale analytic cover of a hypothetical
   proper algebraic DM stack algebraizes by Riemann existence/GAGA for finite
   covers.  Trivial analytic inertia makes the cover a smooth proper
   algebraic space \(X_0\).  Chow's lemma plus resolution gives a smooth
   projective modification; topological \(\pi_1\) is unchanged by a proper
   birational modification of smooth complete complex spaces.  Consequently
   \(\Gamma_0\) would be a projective, hence Kahler, group.

3. **Carlson--Toledo applies exactly.**  Carlson--Toledo, Theorem 8.1, states
   that a cocompact discrete subgroup of \(SO(1,n)\) is not the fundamental
   group of a compact Kahler manifold for \(n>2\).  Thus it applies directly
   to the closed real-hyperbolic six-manifold
   \(\mathbf H^6/\Gamma_0\).  The owner's phrase “for \(n>3\)” is merely a
   weaker threshold, not an error in the \(n=6\) use.  Their Corollary 3.7 and
   Theorem 7.1 also give the advertised rank-at-most-two/factorization result:
   every map from a compact Kahler manifold to a closed real-hyperbolic
   \(n\)-manifold kills homology above degree two.  Hence an equal-dimensional
   projective branched cover would have degree zero, contradicting finiteness
   of positive degree.

4. **The \(\mathbf H^{10}\) stress test has the claimed scope.**  A fresh
   replay of `build_certificate(include_lineage=False)` confirms the
   finite-volume redoubling theorem, \(f(P_N)=N+12\), index \(N\), and volume
   multiplier \(N\).  The seed reflections are integral isometries of the
   signature-\((1,10)\) lattice; since their group has finite covolume, it has
   finite index in the relevant integral orthogonal lattice and is
   arithmetic.  All descendants are finite-index arithmetic reflection
   subgroups.  For \(N>1\) they are visibly nonmaximal because the seed group
   is a proper reflection overgroup.  Thus this family refutes a bound from
   coarse-sphere topology or low cohomology alone, but not a bound for the
   absolutely maximal class.

5. **The revised scope statements are exact.**

   - Section 4.4 now uses only that a holomorphic branch locus is a proper
     complex-analytic subset and hence cannot be a generic real hypersurface;
     it explicitly allows higher-codimension branching over a singular normal
     target.
   - Theorem 3.1 now excludes only a **smooth proper projective-DM / smooth
     proper KSBA-orbifold avatar** and explicitly preserves singular KSBA
     pairs and external algebraic carriers.  This is exactly the reach of the
     proof.

Finally, exact return belongs to the decorated real quotient
\(( [\mathbf H^n/W^+],\iota,\text{selected half},\text{ridge marking})\), not
to the bare unmarked quotient alone.  Section 1 describes these decorations;
the abbreviated notation in the terminal table should be read accordingly.

No owner-file edits were made in this audit.

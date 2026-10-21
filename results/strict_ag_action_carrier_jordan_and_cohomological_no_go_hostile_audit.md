# Hostile audit: strict AG-action carrier report

**Audited owner SHA-256:**
`911b8fd2d64ed3eab83599ca3a3c44ff95c05aff1aa52109a76d04af3d1b9b22`.

**Verdict: PASS.**  I find no substantive mathematical error and no use of
arithmetic-reflection finiteness hidden in the argument.  The report proves
only a conditional transfer theorem and accurately leaves uniform essential
surjectivity/bounded decorations open.

## Checks

1. **Local elementary-​2 rank.**  At a finite vertex the rank-​(n)
   spherical Coxeter group has a forest diagram; a larger bipartition class
   gives ((\mathbf Z/2)^{\lceil n/2\rceil}).  At an ideal vertex the compact
   Euclidean link has a rank-​(n-1) finite point Coxeter subgroup inside its
   affine reflection group, giving
   ((\mathbf Z/2)^{\lceil(n-1)/2\rceil}).  Finite volume excludes hyperideal
   ends without truncating facets, so a finite or ideal vertex exists.  Thus
   (r_2(W)\ge\lfloor n/2\rfloor), sharpened to
   (r_2(W)\ge\lceil n/2\rceil) in the compact case, and consequently
   (n\le2a+1) (respectively (n\le2a)).

2. **Polarized/PGL bound.**  Very ampleness makes
   (\operatorname{Aut}(X,L)\to\operatorname{PGL}_N) faithful.  The adjoint
   representation of (\operatorname{PGL}_N) on (\mathfrak{sl}_N) is
   faithful in characteristic zero.  A faithful representation of
   ((\mathbf Z/2)^r) needs at least (r) characters, hence
   (r\le N^2-1).  The stated bounds (n\le2N^2-1) and, for compact
   chambers, (n\le2N^2-2), follow.  They are deliberately crude but valid.
   The warning that a Jordan-index theorem alone does not control (r_2) is
   correct.

3. **Cohomological bound and finite image.**  Restriction of scalars gives
   rank (d(n+1)); a nondegenerate injection into a rank-​(B) lattice
   therefore gives (d(n+1)\le B).  This alone does not bound discriminants
   or the list of degree-​(d) fields.  The later finiteness statement uses
   the explicitly assumed finite constructible stratification with constant
   exact return, so it is logically valid and not smuggled in from rank
   boundedness.

4. **Triangle stress test.**  For (\Delta(2,3,m)), (m\ge7), the faithful
   Lorentz/Tits action projectivizes faithfully: its image preserves a chosen
   hyperbolic sheet, whereas the only possible nontrivial scalar kernel,
   (-I), exchanges sheets.  Linear automorphism graphs in
   (\mathbf P^2\times\mathbf P^2) have one fixed Hilbert polynomial.  The
   normalized Gram cyclic products give exactly
   (\mathbf Q(\cos^2(\pi/m))), of degree ((m-1)/2) for odd primes (m).
   The revised carrier records the real component and oriented half-spaces,
   removing the possible ambiguity of three unoriented mirror lines.  The
   report expressly does **not** claim these triangles are arithmetic or
   reflection-maximal; accordingly they refute only finiteness from bounded
   strict-action geometry by itself.

5. **(X_3/H21) scope.**  A 2-spherical twist obeys
   (T_S^r(S)\simeq S[-r]), so the displayed curve twists cannot be
   involutive Coxeter generators even though their numerical actions are
   reflections.  A regular automorphism cannot induce the same root
   reflection because it would send an effective ((-2))-class to its
   negative; birational K3 selfmaps do not evade this.  The report correctly
   limits the conclusion to the displayed lifts, explicitly allows other FM
   lifts or different AG hosts, and calls the class-wide strict species
   unproved rather than refuted.  Its rank-​20 and rank-​22 capacity statements
   agree with the companion inertia bounds.

No owner-file edits were made in this audit.

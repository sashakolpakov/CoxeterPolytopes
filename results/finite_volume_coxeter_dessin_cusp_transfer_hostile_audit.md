# Finite-volume Coxeter-dessin cusp transfer: hostile audit

**Verdict (2026-08-09): PASS after a necessary \(\beta=0\) repair and
explicit lattice/cusp hypotheses.**  The final theorem is correct in its
stated admissible-quadratic scope.  The independent replay imports no project
Python module, uses exact rational arithmetic, and ran at 0.025 GiB peak RSS
under an outer 0.5 GiB watchdog.

## 1. Defect found and repaired

The first version asserted for every \(\beta\) that

\[
 (x,r)=0\quad\Longleftrightarrow\quad
 \beta^2h^2+\|\beta z-w\|^2=\rho.
\]

That was false for \(\beta=0\): multiplying by \(\beta\) deletes the
vertical mirror equation.  An exact counterexample in signature \((1,2)\)
is

\[
 u=(1,1,0),\quad v=(1/2,-1/2,0),\quad
 h=1,\quad z=r=w=(0,0,1).
\]

Here \(x^2=1\), \((x,r)=-1\), but the boxed norm equality is nevertheless
\(1=1\).  The primary theorem and certificate were repaired before this
audit was frozen.

The exact universal identity is

\[
 \boxed{\beta^2h^2+\|\beta z-w\|^2-\rho
       =2h\beta(x,r).}
\]

It yields the completed-square equivalence only when \(\beta\ne0\).  When
\(\beta=0\), the two independent statements are

\[
 \|w\|^2=\rho,
 \qquad
 (x,r)=\frac{\alpha-\langle z,w\rangle}{h},
\]

so the vertical mirror is \(\langle z,w\rangle=\alpha\).  The hostile
certifier checks the universal identity on 300 exact rational trials and
replays the information-loss counterexample literally.

## 2. Exact finiteness theorem and its hypotheses

The argument is valid for a fixed rational lattice \(L\), with candidate
vectors in a fixed finite union of fixed-denominator cosets (in the Heegner
case, fixed \(\mu+L\subset L^\vee\)) and finitely many negative norms.  At
each cusp one needs:

1. a primitive rational/integral isotropic \(u\) and one fixed rational
   isotropic \(v\) satisfying \((u,v)=1\);
2. a finite-sided chamber and a complete maximal-parabolic neighborhood of
   rank \(n-1\), equivalently a compact horospherical cell; and
3. a fixed determinant/overlattice extension, or a proved finite list of
   finite-index extensions, so restriction of Heegner vectors retains a
   uniform denominator.

These are substantive, not cosmetic.  Without a fixed lattice/coset
denominator, fixed norm alone does not make \(\beta\), \(w\), or \(\alpha\)
discrete.  Without full parabolic rank, the vertical offset argument has no
compact \(z\)-cell.

Under the hypotheses, projection to \(K=(u,v)^\perp\) sends the fixed cosets
to finitely many cosets of a positive-definite lattice, and the
\(\alpha=(r,v)\) values form a discrete subgroup.  Therefore:

- \(\beta\ne0\) has a uniform lower bound, hence
  \(h\le\sqrt\rho/|\beta|\); removing all cusp tails above these bounds
  leaves a compact thick part, where lattice mirrors are locally finite;
- for \(\beta=0\), the fixed shell \(\|w\|^2=\rho\) contains finitely many
  projected lattice vectors, and for each \(w\), compactness of the
  horospherical cell bounds the discrete offsets \(\alpha\).

There are finitely many ideal vertices, so the union is finite and exactly
enumerable.  Neatness is not required for this stack/orbifold statement; a
neat finite-index level is needed only when the target must be an honest
variety rather than a Shimura stack.

## 3. Field split

The field claim passes in the explicitly stated scope of admissible
quadratic arithmetic lattices of simplest type.  If \(k\ne\mathbb Q\), a
nonzero \(k\)-isotropic vector remains nonzero and isotropic at every real
embedding, contradicting definiteness at a conjugate place.  Thus the
arithmetic group has \(k\)-rank zero and compact quotient.  A
finite-covolume reflection subgroup has finite index in that ambient
arithmetic lattice, so its chamber has the same compact/noncompact type.

Consequently noncompact admissible quadratic arithmetic input forces
\(k=\mathbb Q\).  Over \(\mathbb Q\), Meyer isotropy for indefinite rank
\(n+1\ge5\) shows that a compact quadratic finite-covolume reflection
chamber requires \(n\le3\).

This audit does not extend that assertion to an unspecified non-quadratic
notion of arithmetic hyperbolic lattice.

## 4. Independent toric \(H^6/9\) replay

Starting only from the nine literal vectors in \(I_{1,6}\), the hostile
certifier:

- enumerates all \(\binom96=84\) possible rank-six active sets;
- finds exactly 21 extreme rays: 18 finite, 3 ideal, and 0 hyperideal;
- proves that each ideal neighborhood is a product of affine components of
  ranks \(2+3=5=n-1\);
- checks negative semidefiniteness, corank one, and strictly positive kernel
  relations \((1,2,1)\) and \((1,1,1,1)\), which prove compactness of the
  rank-five horospherical alcove;
- independently decomposes all nine walls at each cusp using
  \(u=H-E_i\), \(v=(H+E_i)/2\); and
- obtains literally seven vertical active walls and two nonvertical walls
  at every cusp, with both nonvertical height caps satisfying \(h^2\le2\).

Thus the calibration is genuinely a complete \(18+3+0\), maximal-parabolic
finite-volume audit.  It is not an inference from finite sidedness.

## 5. Scope of the PASS

What now passes is a finite enumeration gate for fixed Heegner support at
the cusps of a known finite-volume arithmetic quadratic chamber.  It does
not prove that a desired Heegner combination is principal, that a Borcherds
product exists, that a surface section is faithful, or that an input wall
list itself has finite volume.  Those remain separate compiler and Vinberg
gates.

Artifacts:

- hostile certifier:
  `scripts/independent_finite_volume_coxeter_dessin_cusp_hostile_audit.py`;
- hostile archive:
  `data/enumerations/finite_volume_coxeter_dessin_cusp_transfer_hostile_audit.json`;
- hostile tests:
  `tests/test_finite_volume_coxeter_dessin_cusp_transfer_hostile_audit.py`.

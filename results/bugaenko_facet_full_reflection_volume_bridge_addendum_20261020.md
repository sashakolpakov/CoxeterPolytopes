# Full lattice reflection chambers: a valid but ineffective volume bridge

20 October 2026. Paper-only addendum to
`results/bugaenko_facet_binary_extension_and_collar_audit_20261020.md`
(SHA256 `774e87f7abb33680a06879a17a5be9f06e2218bf12a95023bcff5c7046ba1227`).
No proof process, root search, numerical program, Lean or TeX job ran.

**Conclusion.** Full integral reflection chambers over the fixed golden
field do have a uniform volume upper bound in dimension nine. Maximality
of the reflection group is unnecessary for this statement. Together with
the preceding collar and edge estimates it gives a finite normalized
binary reduction. The available constants make that reduction unusable
for a small classification or a new dimension bound. It does not exclude
an additional host, let alone every polytope in the facet class.

## 1. The full integral group supplies the missing hypothesis

Let `O=Z[phi]`, let R be a free admissible O-lattice of rank ten, and write

    Gamma=O^+(R),     W=<all reflections preserving R>.

Assume W has finite covolume. Gamma is the full integral group, hence is
a congruence arithmetic lattice. The set of its reflections is invariant
under conjugacy, so W is normal in Gamma. W has finite index: covolumes
give `[Gamma:W]=Vol(H9/W)/Vol(H9/Gamma)`. For a W chamber C,

    Gamma=W semidirect-product A,     A=Stab_Gamma(C).

Indeed, W acts transitively and freely on chamber interiors; multiplying
any Gamma element by a suitable W element makes it preserve C. Here C
is compact: the nonphysical positive-definite embedding makes the form
anisotropic over Q(phi), and a finite-index subgroup of the resulting
cocompact lattice remains cocompact. Consequently A is finite.

There is an elementary uniform bound on its order. Reduction modulo3
injects every finite subgroup of GL10(O) into GL10(F9). To see the
torsion-freeness of the kernel, take a nonidentity prime-order element
`g=I+3^j B`, with `j>=1` maximal and some B entry not divisible by3.
For prime order different from3 the binomial expansion has a unique
lowest-valuation term; for order3 its first term has valuation j+1 and
all later terms have larger valuation. Both contradict `g^ell=I`.
Thus

    m=|A| <= |GL10(F9)| < 9^100 = 3^200.             (1)

This argument controls the quotient in the *full integral group*. It
does not assert that W is a maximal reflection group in Isom(H9), or
that an arbitrary finite-index reflection subgroup is congruence.

## 2. The published analytic input and its precise application

Use the following established ingredients from Agol--Belolipetsky--
Storm--Whyte: the finite spherical quotient conformal-volume bound
(Theorem3.1), the orbifold spectral inequality (Theorem4.2), and the
congruence spectral gap (Theorem5.9). Their Section6 applies these to a
normal reflection subgroup with finite chamber stabilizer. The same
argument applies here because normality and congruence were established
directly above; no maximality step is needed. Primary source:
[author-hosted paper, Sections3--6](https://homepages.math.uic.edu/~kwhyte/ahfinite.pdf).

For clarity, the quotient argument is explicit. Conjugate the finite A
so its common fixed point is the center of the conformal ball. Then
`H9/Gamma=C/A` conformally embeds in `S9/A`. The proof of Theorem3.1
cuts an equivariantly embedded sphere by at most
`K=binom(10m,2)<50m^2` hyperplanes and uses at most
`2 sum_(j=0)^9 binom(K-1,j) <= 20*50^9*m^18` pieces. Therefore

    V_PC(H9/Gamma) <= 20*50^9*m^18*Vol(S9).

With `lambda_1(H9/Gamma)>=15/4`, the spectral inequality gives

    Vol(C) <= (12/5)^(9/2) * 20*50^9 * m^19 * Vol(S9)
           < U,
    U=(12/5)^(9/2) * 20*50^9 * 3^3800 * pi^5/12.   (2)

This is a genuine uniform bound for the *full* reflection chamber, not
the generally false uniform bound for arbitrary reflection subgroups.
For comparison, Belolipetsky's congruence-reflection bound applies
directly only when the reflection group itself is congruence; that
extra assumption has not been made here. See
[Definition2.1 and Section3](https://arxiv.org/pdf/1008.0934).

## 3. What follows for the marked binary extension

In the setting of the preceding audit, `R=E8+B` is the even orbit
lattice of the hypothetical compact reflection group. All its original
wall reflections preserve R, so W already has finite covolume and (2)
applies. Rebuilding the lattice from a maximal reflection overgroup is
valid but unnecessary. Extra R mirrors may still subdivide the original
metric facet; this proof does not identify the original polytope with C.

Let

    kappa=phi^4*pi^4/(24*|W(E8)|).

The independently proved collar gives `Vol(C)>=kappa*ell` for the actual
binary edge adjacent to the canonical pure-B root. The preceding
discriminant identities give, respectively,

    odd glue:   |N(det B)| <= 16*sinh(2ell)^2,
    direct sum: |N(det B)| <= 256*sinh(2ell)^2.

Hence both cases satisfy the finite explicit bound

    |N(det B)| <= 256*sinh(2U/kappa)^2.              (3)

A bounded norm does give finitely many reduced marked B, not infinitely
many unit associates. In the odd normalization
`B=[[-2phi,1],[1,2c]]`, put `D=1+4phi*c`. For every integer j, replace
the second basis vector t by

    phi^(3j)*t + (phi^(3j)-1)/(2phi)*u0.

Its pairing with u0 remains1; its determinant changes by the unit
`phi^(6j)`. Choose j so
`sqrt(|N D|)<=D<phi^6*sqrt(|N D|)` at the physical embedding.
Both embeddings of D, and therefore both coefficients of c, are then
bounded by (3). In the direct normalization use `t -> phi^j*t` and
the corresponding factor `phi^(2j)`. Thus a finite reduced list follows
without assuming that the earlier list of closed hosts was exhaustive.

## 4. Why this does not advance the dimension bound

The conformal bound (2) is enormous. More decisively, even if one somehow
knew `A=1`, the sharper direct embedding would give only

    Vol(C) <= (12/5)^(9/2)*Vol(S9),

of order10^3. The collar coefficient is about `4*10^-8`; this still
allows an edge length of order10^10 and a determinant cutoff exponential
in that number. Restoring (1) makes it far worse. No enumeration is
proposed.

The precise result is therefore a valid fixed-field finite reduction,
not a quantitatively useful exclusion. A useful next theorem must
strengthen the geometry substantially, for example by bounding the
actual marked binary edge independently of this tiny E8-sector volume,
or by forcing a second chamber symmetry uniformly. Neither such bound
nor such symmetry has been established here.

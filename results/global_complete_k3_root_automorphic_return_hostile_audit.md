# Hostile audit: global complete-root K3/automorphic return

**Verdict (2026-08-09): pass as an exact-return prototype, with strict marked
and complete-root scope; not a pass of the flagship terminal-AG-leverage
contract.**

The theorem in
[`global_complete_k3_root_automorphic_return.md`](global_complete_k3_root_automorphic_return.md)
really does remove K3 embeddability from the input.  The relevant primary
statement classifies **all** even elliptically 2-reflective hyperbolic
lattices, not merely lattices already known to occur as K3 Picard lattices.
The same source then invokes global Torelli to identify that exhaustive list
with the finite-automorphism K3 Picard list.  The automorphic `51/18` split
has the exact universe claimed: all hyperbolic 2-reflective lattices of rank
at least six, and no assertion in ranks three through five.

The result remains a synthesis of published classifications.  It neither
produces a new abstract polytope nor moves the general mixed-root or
totally-real compact boundary.

There is a decisive order-of-proof limitation.  Essential surjectivity and
the bound `n<=18` are obtained by first applying the Nikulin--Vinberg
classification to the Coxeter/lattice input and then interpreting its list
as K3 Picard lattices.  Sterk's theorem and the Cox-ring equivalence add
native algebraic-geometric structure, but do not independently imply a new
polytope bound after the return.  Thus this report closes the complete-root
denominator and supplies a benchmark for the desired inverse theorem; it
does not fulfill the user's stronger `P -> surface -> new theorem on P`
mandate.

## 1. Exact primary-source quantifiers

Gritsenko--Nikulin use signature `(n,1)` and positive norm-two roots; the
report uses the global negative form and signature `(1,n)`.  This changes no
reflection hyperplane or index statement.

In arXiv:1602.08359v2, Section 3.2, they define an even hyperbolic lattice to
be *elliptically 2-reflective* exactly when

\[
              [O(S):W^{(2)}(S)]<\infty,
\]

and state equivalently that its complete norm-two fundamental chamber has
finite volume.  Immediately afterward the source says, in full-universe
language, that Nikulin (ranks other than four) and Vinberg (rank four)
classified **all** such lattices and introduces “The list of all
elliptically 2-reflective even hyperbolic lattices of rank at least 3.”  The
list ends with one lattice in each rank 15 through 19 and the sentence that
there are none in rank at least 20.

This closes the possible loophole that the list might have assumed a K3
embedding.  Only after giving the exhaustive arithmetic list does Remark 3.1
invoke global Torelli and say that the lattices on that list give all Picard
lattices `S_X=S(-1)` of complex K3 surfaces with finite automorphism group,
and that the simple norm-two roots give all smooth rational curves.  K3
carrier existence is therefore a theorem output.

Wang, arXiv:1906.10459v2, Section 7, independently defines hyperbolic
2-reflectivity by the same finite-index condition and reproduces the exact
rank table

```text
rank  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 >=20
count 26 14  9 10  9 12 10  9  4  4  3  3  1  1  1  1  1    0.
```

His Theorem 7.2 begins with signature `(1,n)`, `n>=5`, hence rank at least
six.  It names 18 lattices with no automorphic correction and says every one
of the *other 51 lattices* has at least one correction.  Since the published
rank-six-through-rank-nineteen counts sum to 69, the theorem's universe is
exactly `51+18`, neither a sample nor a bounded-weight computation.

## 2. Same-polytope return

For every fixed marked lattice `S`, primitive K3 embedding, orientation, and
moduli-component package types form a finite list described by complement
isometry classes and discriminant anti-isometry double cosets.  The
classification proves the list nonempty.  Each package type has a nonempty
exact-Picard period locus.  Those loci generally contain continuously many
period points and K3 surfaces: finiteness refers only to the discrete package
and component types, not to the carrier objects themselves.

The K3 ample-cone theorem uses every square-minus-two mirror.  Since the
input chamber is a component of that same complete arrangement, a Weyl
change of marking sends the nef cone to the specified chamber.  This returns
the saturated lattice and the full normalized wall Gram, not merely the
nerve.  No effectiveness or hidden-root guess is left.

For an arbitrary even `S` and complete chamber `C`, the resulting equivalence
is therefore exact:

```text
C finite volume
  <=> S is one of the 118 lattices
  <=> exact-Picard K3 carrier with finite Aut
  <=> rational-polyhedral nef/effective cone
  <=> finitely generated K3 Cox ring.
```

The last two equivalences are Sterk plus Artebani--Hausen--Laface; they are
not inferred merely from the Coxeter diagram.

## 3. Rank-26 Conway hostile calibration

The ambient Conway lattice `II_(1,25)` is not a counterexample to the
rank-at-most-19 classification.  Its complete norm-two chamber has one wall
for every Leech-lattice vector, and

\[
 O^+(II_{1,25})/W_2\cong \Lambda\rtimes O(\Lambda)
\]

contains the infinite translation lattice `Lambda`.  Hence `W_2` has
infinite index; the Conway chamber is not elliptically 2-reflective and does
not have finite hyperbolic volume.  Finite-sided spherical faces of that
chamber can still produce the separate Conway-face examples in the ledger,
but the ambient chamber does not enter `Root2FV`.

This calibration is important because the words “Conway reflection
chamber” alone do not imply finite covolume.

## 4. Counting and novelty boundaries

- `118` is the number of oriented marked integral lattice/chamber types.
  Complete root chambers in one chosen positive component are
  `W_2`-conjugate.  No claim is made that forgetting the integral lattice
  leaves 118 distinct real metric-polytope isometry types.
- The compact split is independently exact: six rank-three and two rank-four
  chambers are compact.  Meyer excludes compactness from rank five onward,
  leaving 110 finite-volume noncompact types.
- The `7/12` high-rank automorphic split is a literal set intersection with
  Wang's 18 names.  It does not prove that a no-correction chamber is not
  reflective; it is already reflective by hypothesis.  It proves only that
  no `U(m)+S` has the required 2-reflective modular form.
- The 51 positive corrections and 18 obstructions are published.  Calling
  them new examples would be false.  The new project artifact is the exact
  inverse synthesis and its high-rank cross-ledger.
- The `118` count and `n<=18` bound cannot be advertised as conclusions newly
  extracted from K3 boundedness: they were needed first to prove that every
  input has a K3 carrier.  The exact equivalence is useful, but the
  flagship's independent last AG arrow remains open.

## 5. Frozen source witnesses and replay

Pinned downloads used for this audit:

```text
arXiv:1602.08359v2 source archive
80f4b3885779070a80f201f04f38a1aaaf56e65aa485bc1943596a4d78c41857
main TeX
c63ea76045f584e82443164e31dac4d470b6d188707b491f4799d3453b0bc96c

arXiv:1906.10459v2 source archive
ec3a229fc18e6a9fc85f16bd27d4f6121ebb5507a7896b161db8436d3c31fc13
main TeX
54d000cd5a1546ba92e5896f63cbebe21f305adbe0c1fccf589955c632f331b6

arXiv:2003.08985 source archive
3e445ed01446b47ded6c007b395308de34e566b736f0ba982470bc497adeb1e9
main TeX
e6eeaaa4a3f375f209762718e8ef9a361b0a95356cc8be2e9f99a6713abf5df3
```

The arithmetic replay is:

```bash
python3 scripts/certify_global_k3_root_automorphic_return.py
python3 -m pytest -q tests/test_global_k3_root_automorphic_return.py
```

It verifies the rank total, compact/noncompact split, exact `51/18` universe,
and `7/12` high-rank intersection.  It does not recompute the published
Nikulin--Vinberg or Wang classifications.

# Independent hostile audit: unconditional maximal arithmetic reflection-chamber finiteness

**Verdict (2026-08-11): PASS as an unconditional classical theorem; PASS for
the noncircular finite-field factorization; NO positive AG-only theorem.**

The owner theorem is mathematically correct in its declared scope.  For
finite-volume arithmetic reflection groups in dimensions at least two it
gives

\[
 n\le 29,\qquad [k:\mathbf Q]\le25,
\]

and, for reflection-maximal inputs, finitely many minimal Vinberg ground
fields, native pairs, and unmarked normalized chambers.  The noncompact
refinement is also correct: \(k=\mathbf Q\), and Esselmann's integral
Lorentz-lattice theorem leaves \(n\le19\) or \(n=21\).

The terminal finiteness statement already follows directly from ABSW (or
Nikulin 2007).  The owner's separate finite-native-pair factorization does
not use the ABSW finiteness conclusion: it combines the imported Vinberg
dimension bound and Belolipetsky--Linowitz degree bound with Nikulin's
bounded full-rank root frame, a finite algebraic-integer coefficient box,
Nikulin's fixed-\((k,m)\) reflective-lattice theorem, and reflection
maximality.

This result is unconditional because it has no carrier-existence
hypothesis.  It is not an algebraic-geometric derivation.  The accompanying
[AG universality--leverage no-go](unconditional_ag_absw_universality_leverage_no_go.md)
shows why the distinction is substantive: the presently universal
rational-surface carrier is a weak numerical encoding with no independent
finiteness leverage, while the strong complex-K3 routes are not universal.

## 1. Frozen owner and scope-evidence bytes

The audited terminal owner package is:

| artifact | SHA-256 |
|---|---|
| `results/unconditional_maximal_arithmetic_reflection_chamber_finiteness.md` | `6eee91512840287ce47a195a317fa4bbfd19678a95b52869ec0a45775d61e67b` |
| `scripts/certify_unconditional_maximal_arithmetic_reflection_chamber_finiteness.py` | `c32c7969f900093db388718f3ac769f669d42e43b6c2a95c76297eb34899f541` |
| `data/enumerations/unconditional_maximal_arithmetic_reflection_chamber_finiteness.json` | `2bce48c0a7156caea3f4b4890411d770450912b4d4b2d3827ce93cc7e9730acc` |
| `tests/test_unconditional_maximal_arithmetic_reflection_chamber_finiteness.py` | `4105f4f3d08dfdcbc1d98f272aae6168a99221f6e426f1bcf1fb5b50303d6702` |

The independently frozen scope memo is:

| artifact | SHA-256 |
|---|---|
| `results/unconditional_ag_absw_universality_leverage_no_go.md` | `fcc3acb8bee4751bef083cd00a179f9b542602561ba1bf81d89f374dd5523150` |

That pin matters: the owner theorem and the AG no-go answer different
questions.  The former proves the global theorem with classical arithmetic
inputs; the latter prevents that theorem from being relabelled as a new
AG-only proof.

## 2. Direct terminal route

ABSW define a maximal reflection group exactly as a discrete group generated
by the facets of a finite-volume hyperbolic polyhedron which has no proper
reflection overgroup.  Their Theorem 6.1 states that only finitely many
conjugacy classes of arithmetic maximal reflection groups occur.  Thus the
owner's terminal theorem needs no carrier and no new finite-field argument.

The passage from group conjugacy to unmarked normalized chambers is sound.
The mirrors of a discrete reflection group form a locally finite
arrangement.  A generic compact path between regular points crosses only
finitely many mirrors; reflecting successively across them gives a gallery
between the two complete chambers.  Hence all complete chambers form one
group orbit.

The two maximality conventions also agree here.  If a finite-covolume
reflection group \(W\) lies in a discrete reflection group \(W'\), then
\(W'\) has finite covolume and the covolume/index formula gives
\([W':W]<\infty\).  Arithmeticity is invariant under commensurability, so an
absolute reflection overgroup of an arithmetic \(W\) is arithmetic.

## 3. Independent audit of the finite-native-pair factorization

### 3.1 The reflective host is legitimate

Arithmeticity first supplies an admissible native \(k\)-form and an
arithmetic lattice \(L_0\).  The intersection \(W(P)\cap O(L_0)\) has finite
index in \(W(P)\).  Passing to its finite-index normal core and summing the
finitely many \(W(P)\)-translates of \(L_0\) produces a finite torsion-free
\(\mathcal O_k\)-module \(L\) preserved by every facet reflection.  This
argument does not assume that \(\mathcal O_k\) is principal.  A common
rescaling clears the pairing ideal and changes no reflection.

Consequently

\[
              W(P)\le R(L)\le O^+(L).
\]

Both endpoint groups have finite covolume, so \(W(P)\) has finite index in
\(O^+(L)\), and \(R(L)\) has finite index in \(O(L)\).  Thus \(L\) satisfies
the literal reflectivity hypothesis of Nikulin 5.2.1.  Maximality is not used
at this stage.

### 3.2 Only cyclic Cartan products are required to be integral

Apply Nikulin's constant-62 selection lemma to a complete chamber \(M(L)\)
of \(R(L)\), not to an arbitrary physical star in the original chamber.  It
selects \(m=\operatorname{rk}_{\mathcal O_k}L\) independent roots with a
connected nonorthogonality graph and finite-index span.  Put

\[
 a_{ij}={2(\delta_i,\delta_j)\over(\delta_i,\delta_i)},
 \qquad c_{ij}=a_{ij}a_{ji},
 \qquad b_C=\prod_{(i,j)\in C}a_{ij}.
\]

For projective root modules an individual \(a_{ij}\) need not lie in
\(\mathcal O_k\).  The valid Vinberg assertion is that every closed cyclic
product \(b_C\) is an algebraic integer; in particular each two-cycle
\(c_{ij}\) is integral.  The owner now states exactly this weaker and
sufficient claim.

For a simple directed cycle, symmetrizability gives

\[
                       b_C^2=\prod_{(i,j)\in C}c_{ij}.
\]

At the distinguished real place, Nikulin's frame has \(c_{ij}<62^2\).  At
every other real place, definiteness and Cauchy--Schwarz give
\(0\le\sigma(c_{ij})\le4\).  Hence, for \(|C|\le m\le30\),

\[
                         |\sigma(b_C)|\le62^m\le62^{30}.       \tag{3.1}
\]

The independent replay includes a symmetrizable rank-three Cartan packet
whose individual entries \(1/2,1/3,1/6\) are nonintegral but whose two-cycle
and oriented three-cycle products are all integral.  This is a regression
against reinstating the earlier, overly strong entrywise-integrality claim.

### 3.3 The selected cycles generate the original ground field

The selected frame has full Lorentzian rank.  Nikulin 1982, Theorem
2.1.3(2), therefore gives the same minimal ground field for its hyperbolic
subdiagram as for \(M(L)\).  On the other side,
\(W(P)\le R(L)\) has finite index.  Vinberg's least adjoint
field/ring-of-definition theorem makes the minimal Vinberg field invariant
under this finite-index passage.  Thus

\[
                   k(P)=k(M(L))
                       =\mathbf Q(b_C:C\text{ a simple frame cycle}). \tag{3.2}
\]

There is no hidden infinite generator list in (3.2).  Any repeated closed
walk decomposes multiplicatively into simple cycles and two-cycles.  A
rank-\(m\) graph has at most

\[
       S(m)=\sum_{r=2}^{m}{m!\over(m-r)!r}
\]

oriented simple-cycle slots modulo cyclic rotation.  The independent replay
checks the decomposition on every loop-free closed walk of length at most
seven on four labelled vertices and separately recomputes \(S(30)\).

Both field equalities in (3.2) are essential.  A bounded frame in the wrong
finite-index reflection group would not control \(k(P)\); the repaired owner
proof supplies precisely the missing inheritance and commensurability
bridges.

### 3.4 The algebraic-integer box is genuinely finite

Let \(B=62^{30}\).  Every generator in (3.2) has degree at most 25 and all
of its conjugates have absolute value at most \(B\).  If its minimal
polynomial has degree \(e\le25\), Vieta gives

\[
              |q_j|\le {e\choose j}B^j\le{25\choose j}B^j.
\]

There are finitely many integral coefficient vectors in these explicit
boxes, hence finitely many possible algebraic integers.  Since the number of
simple-cycle slots is bounded by \(S(30)\), only finitely many tuples and
therefore only finitely many generated fields occur.  The independent
certificate recomputes the complete 25-coordinate Vieta vector, its digest,
and a finite upper bound for the union of the degree-\(1\) through
degree-25 coefficient boxes.  It also compares those values with the frozen
owner archive rather than importing the owner's Python functions.

This use of Nikulin's bounded-frame selection lemma is not the conclusion of
Nikulin 5.2.1.  The logical order is:

1. the selection lemma plus imported \(m\le30\) and \([k:\mathbf Q]\le25\)
   proves that the set of fields and native pairs is finite;
2. only then does Nikulin 5.2.1 give finitely many reflective lattice
   similarity classes for each fixed pair.

There is therefore no ABSW-finiteness or fixed-data circularity in the
finite-pair step.

### 3.5 Fixed data and maximality close the count

Nikulin 5.2.1 applies to finite torsion-free \(\mathcal O_k\)-modules and
already includes projective ideal-class phenomena.  For each of the finitely
many \((k,m)\), it leaves finitely many reflective lattices up to similarity.
Similarity is the right quotient because a common scale cancels from every
reflection formula and normalized Coxeter Gram entry.

For a reflection-maximal input, \(R(L)\) is a discrete reflection overgroup
of \(W(P)\), so \(W(P)=R(L)\).  Under arithmetic-only maximality the same
conclusion holds because \(R(L)\) is finite-index in \(O(L)\) and hence
arithmetic.  The complete-chamber orbit argument from Section 2 then gives
at most one unmarked normalized chamber class per native similarity host.

## 4. Numerical refinements and their exact status

- Vinberg's imported crystallographic nonexistence results give \(n<30\),
  hence \(m\le30\), for every finite-volume arithmetic reflection input.
- Belolipetsky--Linowitz, Corollary 5.1, gives
  \([k:\mathbf Q]\le25\).  This is a published arithmetic input, not a
  consequence of carrier geometry and not the ABSW finiteness conclusion.
- In the noncompact case, a cusp makes the admissible native form isotropic
  over \(k\).  A nonidentity real embedding would make the same nonzero
  isotropic vector lie in a definite form, which is impossible.  Thus
  \(k=\mathbf Q\).  Esselmann's theorem for reflective integral lattices then
  gives \(n\le19\) or \(n=21\).  Dolgachev's survey states this exact result
  as Theorem 4.7; the rank-22 Borcherds lattice realizes \(n=21\).
- The weak rational-surface carrier bound \(4[k:\mathbf Q](n+1)-3\le2997\)
  is a valid corollary of the imported arithmetic bounds.  It is not the
  source of those bounds: its divisors may be noneffective numerical cuts and
  its rational Chow correspondences need not be genuine automorphisms.

## 5. Source pins and replay boundary

The following audit-time source bytes were independently rehashed:

| source | SHA-256 | locus |
|---|---|---|
| ABSW 2008 | `974ac7c5141ca1c45b515f145215ad115aba867e42551a9b103233ea52930ee4` | definition in the Introduction; Theorem 6.1 |
| Vinberg 1971 | `16a4f497688e6f968b4c63d81a6b049b3f096f1cd4f43e458b2882a5c44ded38` | least definition field/ring and cyclic products, Theorems 1, 3, 5 |
| Vinberg expanded dimension paper | `fd25f8f892eed27e4303a8f3252fd384d8954587a6c1c67f1efed2c38aabe44f` | compact and arithmetic noncompact exclusions for \(n\ge30\) |
| Nikulin 2007 | `9a820fef0b88158b283057f7d97acc89130610a9a4aacc6c53a8215e23ce874a` | Theorem 2.2 |
| Nikulin 1982 | `696652b6a67c6ed420de6c1cb61d9422e3619e28e081a78131a5323ce2edee9d` | Theorem 2.1.3(2); Theorem 5.2.1 and constant 62 |
| Nikulin 2011 | `716e6e7333efbfaa22060b60dc97e8303363a4914409ccdd56ef7a90c923c7af` | degree-transition bounds |
| Belolipetsky--Linowitz v2 | `026e0f635f71e214807e0b505bf7297673eb9e0d47a6a91b7bd6dec134eab97f` | Corollary 5.1 |
| Dolgachev survey v3 | `9f29e1e87eee55f1d66433ccbabb5be0833d3c691f9e67a1643e63f7bf44fc26` | Theorem 4.7 quoting Esselmann |

The Esselmann primary journal PDF was not locally available.  Its DOI and
journal record were checked, and the exact gap \(n\ne20\) was independently
confirmed in Dolgachev's Theorem 4.7.  The replay script treats every
temporary PDF as an optional audit-time pin: if present its hash must match,
but a clean checkout need not contain it.

The machine certificate verifies frozen repository bytes, exact owner
archive replay, the cyclic-product algebra, loop decomposition, cycle-slot
count, Vieta coefficient box, and declared logical scope.  It does not prove
the imported literature.  That boundary is explicit rather than hidden
behind a green test.

## 6. Exact answer to the AG question

An honest unconditional theorem **can** now be claimed, but only in the
following sense:

> Classical arithmetic results unconditionally give the ABSW finiteness
> theorem, \(n\le29\), \([k:\mathbf Q]\le25\), the stated noncompact
> refinement, a finite set of native pairs, and finitely many maximal
> normalized chambers.

No honest positive **AG-only** theorem with the same global quantifier is
currently proved.  Bare Neron--Severi lattices do not determine nef or
effective chambers; arbitrary divisor cuts and arbitrary Chow
correspondences give a universal but tautological encoding; exact complex
K3 root geometry is not universal; and an effective negative prime divisor
cannot itself be the reflecting class of a regular automorphism because the
reflection sends its class to its negative while automorphisms preserve the
effective cone.  A bounded mixed curve/action fundamental-domain carrier in
arbitrary characteristic remains a precise research target, not a theorem.

## 7. Independent replay

- [independent certificate generator](../scripts/audit_unconditional_maximal_arithmetic_reflection_chamber_finiteness.py)
- [independent certificate](../data/enumerations/unconditional_maximal_arithmetic_reflection_chamber_finiteness_independent_hostile_audit.json)
- [focused regressions](../tests/test_unconditional_maximal_arithmetic_reflection_chamber_finiteness_independent_hostile_audit.py)

Run

```text
python3 -m scripts.audit_unconditional_maximal_arithmetic_reflection_chamber_finiteness --check-output
python3 -m pytest -q tests/test_unconditional_maximal_arithmetic_reflection_chamber_finiteness_independent_hostile_audit.py
```

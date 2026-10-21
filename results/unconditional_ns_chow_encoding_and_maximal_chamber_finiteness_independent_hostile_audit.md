# Independent hostile audit: unconditional NS/Chow encoding and maximal chamber finiteness

## Verdict

**PASS after three scope repairs, with the two conclusions kept logically
separate.**

The terminal synthesis proves the following conjunction for finite-sided,
finite-volume arithmetic Coxeter chambers in the standard admissible Vinberg
lattice-form model, in dimensions at least two.

1. Every input has a weak rational-surface divisor-Tate carrier with an exact
   distinguished eigenslice and faithful rational Chow action.  After importing
   the classical bounds (n\le 29) and ([k:\mathbf Q]\le 25), one may use at
   most (2997) blow-ups.  Reflection maximality is not used here.
2. Reflection-maximal inputs form finitely many unmarked normalized chamber
   classes.  This is classical: it follows directly from ABSW (or Nikulin
   2007), and it also has the separately audited bounded-frame/finite-field/
   fixed-data factorization.  The weak carrier is not a premise of this
   finiteness proof.

The hostile pass initially found three presentation defects.  The opening and
terminal summary dropped the owner's lattice-form qualifier; the raw ambient
trace pairings of the transported divisors were liable to be mistaken for the
distinguished facet Gram matrix; and the machine certificate retained the
unqualified all-input string after the prose was repaired.  Before the terminal
pin, the owner restricted every weak-carrier summary to the standard admissible
Vinberg lattice-form model, identified the preserved Gram matrix as that of
(\operatorname {pr}_{E_{\sigma _0}}D_i), explicitly disclaimed equality with
the raw pairings ((D_i,D_j)_X), and repaired the archive input string and its
test.  No terminal owner byte was edited by this independent audit.

## 1. Frozen terminal packages

### Synthesis owner

| artifact | SHA-256 |
|---|---|
| [report](unconditional_ns_chow_encoding_and_maximal_chamber_finiteness.md) | `ab17efccf6961be32a45d11ae6683f7d6637c42162626653b145ce0529f2336b` |
| [certifier](../scripts/certify_unconditional_ns_chow_encoding_and_maximal_chamber_finiteness.py) | `6711f6c587c6f20c0df596ed5eaf5d385caff4cc360c0087392a920835141438` |
| [archive](../data/enumerations/unconditional_ns_chow_encoding_and_maximal_chamber_finiteness.json) | `c61c6b16d9268455ec82134c7c012b9690f40bfe4761bd9d1a9d987f67353565` |
| [tests](../tests/test_unconditional_ns_chow_encoding_and_maximal_chamber_finiteness.py) | `3f71fcb06b34a5ca52cf4c13297aaecee19ed6c31012bfc8b77aa3b97fa90814` |

### Classical owner and independent audit

| package | report | script | archive | tests |
|---|---|---|---|---|
| owner | `6eee91512840287ce47a195a317fa4bbfd19678a95b52869ec0a45775d61e67b` | `c32c7969f900093db388718f3ac769f669d42e43b6c2a95c76297eb34899f541` | `2bce48c0a7156caea3f4b4890411d770450912b4d4b2d3827ce93cc7e9730acc` | `4105f4f3d08dfdcbc1d98f272aae6168a99221f6e426f1bcf1fb5b50303d6702` |
| [independent audit](unconditional_maximal_arithmetic_reflection_chamber_finiteness_independent_hostile_audit.md) | `1f1592ef207580b765d56bdca9461429456ef9fa46eb180b0ae65d69dc68ef3c` | `e6fc5f7543093b6167e6d46612b64c7247007cc16b06f081e14260a283193400` | `87ec814ae05d29b04ac068296c4f995144264fed353fa02f607eb22841cff734` | `2ab582e72677922fa772ff57008cb17e307556f5bd595ca767f6130fe7b251ea` |

### Weak-Chow owner and independent audit

| package | report | script | archive | tests |
|---|---|---|---|---|
| owner | `df38b26e0cffd8a391834be922869153af5badf41a632e74a48b472b332e101a` | `d3a19cfe7069b811721a77a624952d92ab5bc7c5574d9e332d2ff8c7f0205eba` | `d75897bbf8d80a662a7b9d3aa3ae729052494e20633d835566db70870bb0fc9e` | `6da7b9e4590d7e531a7f3502ba2c971e806cbd3f259cb4722169553b7532c802` |
| [independent audit](universal_weak_chow_action_carrier_and_leverage_no_go_independent_hostile_audit.md) | `91345590b506ef7c7d4287be5a4d94c7563ba7613167ffbe34272e6b7fc422c1` | `98765eaae927d9265cf69591f02c31d1d0adc71653936d2199f53eb75b16a347` | `a9c19a9f2dce89c03182c9dc7ad8248a0ab75e32e5e58cd1bca6fac7683764f4` | `46be65edc68679ed6f43da75f59393c69e83630810f8c2f7404c3c40381cdc42` |

The independent replay pins every byte in all five four-artifact packages.  It
does not import any owner or branch-audit producer.

## 2. Literal weak-carrier quantifier and the lattice-form boundary

The weak theorem starts with more than an abstract real Coxeter matrix.  Its
input is an admissible totally real quadratic module over the minimal Vinberg
field (k), together with a full (W(P))-stable integral
(\mathcal O_k)-lattice and oriented facet-root lines.  At the distinguished
embedding the form has signature ((1,n)), while all other real embeddings
are negative definite.  This is exactly the standard Vinberg
quadratic/lattice-form model now stated in every owner interface.

For arithmetic reflection groups this model is the classical one used by
Vinberg and Nikulin.  The terminal synthesis nevertheless keeps the qualifier
literal rather than silently claiming a carrier for an unspecified broader
notion of arithmetic lattice.  The weak carrier itself needs neither finite
volume nor maximality once the marked admissible datum is supplied; finite
volume is used only to place the synthesis in the global arithmetic class and
to import the uniform published bounds.

The stable-lattice step is valid.  Intersect a starting arithmetic lattice
stabilizer with (W(P)), take a finite-index normal core, and sum its finitely
many (W(P))-translates.  The sum is a full finite torsion-free
(\mathcal O_k)-module preserved by every side reflection.  A totally positive
trace scalar clears the pairing ideal, so

\[
 B(x,y)=\operatorname {Tr}_{k/\mathbf Q}(\beta q(x,y))
\]

is integral and has signature ((1,R-1)), where
(R=[k:\mathbf Q](n+1)).

## 3. The (4R-3) construction and the exact slice

Rational congruence diagonalization writes the trace form as

\[
 \langle a_0,-a_1,\ldots,-a_{R-1}\rangle ,\qquad a_i>0.
\]

One negative coordinate represents (a_0) as a difference of two squares,
and four fresh negative coordinates represent every (a_i), (i>0), as a
sum of four squares.  After one common denominator is cleared,

\[
 j^t I_{1,b}j=D^2B,\qquad b\le1+4(R-1)=4R-3.
\]

The target is (\operatorname {NS}(\operatorname {Bl}_b\mathbf P^2)).
For (x\in E_{\sigma _0}), restriction of scalars gives

\[
 (x,j(r_i))_X=D^2\sigma _0(\beta)q_{\sigma _0}(x,r_i).
\]

The multiplier is positive and common to every wall.  Thus the projected
normal (\operatorname {pr}_{E_{\sigma _0}}j(r_i)), not the raw trace vector
viewed in all real summands, returns the original oriented hyperplane and
normalized Gram data.

The independent executable uses (k=\mathbf Q(\sqrt3)) and the admissible
rank-three form

\[
                       q=\langle\sqrt3,-1,-1\rangle .
\]

Its rank-six trace form has signature ((1,5)).  The replay diagonalizes the
multiplication-by-(\sqrt3) action into one Lorentzian and one negative
definite real block, checks two oriented walls on the distinguished block,
and verifies that their raw trace Gram matrix is not their projected Gram
matrix.  This is a regression against the ambiguity repaired in the owner
report.

## 4. Literal Chow composition is a group action

In an intersection-unimodular divisor basis with matrix (J), the cycle

\[
                 Z_A=\sum_i e_i^*\times A(e_i)
\]

has product-basis coefficient matrix (C_A=J^{-1}A^t).  Middle-factor
contraction gives

\[
                 C_AJC_B=C_{BA}.
\]

This is the correct order: the left displayed cycle is applied first.  Under
the standard convention that (\Gamma_g\circ\Gamma_h) means first (h),
then (g), one obtains

\[
                 \Gamma_g\circ\Gamma_h=\Gamma_{gh}.
\]

For the orthogonal projector (\Pi) and zero-complement lift (L_g), one has
(L_gL_h=L_{gh}), (L_1=\Pi), and (L_s^2=\Pi) for a reflection.  Source
self-adjointness gives a symmetric coefficient cycle, hence literal
self-transpose on the Chow summand.

The new replay independently uses a noncommuting pair in the
((2,3,\infty)) triangle, embedded in (I_{1,5}).  It checks both composition
orders: the standard order equals (C_{gh}), while the reversed order does
not.  It also checks the (m=2) and (m=3) relations, the nontrivial
parabolic word, the projector unit, and self-transpose reflection cycles.

## 5. Arithmetic bounds and the finite-box factorization

The numerical chain is exact:

\[
 n<30\Longrightarrow n+1\le30,
 \qquad [k:\mathbf Q]\le25,
\]

so

\[
 R\le25\cdot30=750,
 \qquad b\le4\cdot750-3=2997.
\]

Both inequalities on (n) and (k) are imported classical arithmetic
theorems.  The weak surface construction proves only (b\le4R-3).

The separately audited ABSW-free factorization is also sound.  For the
complete chamber of the native reflective lattice, Nikulin selects a
connected full-rank frame with

\[
 c_{ij}={4(\delta_i,\delta_j)^2\over
              (\delta_i,\delta_i)(\delta_j,\delta_j)}<62^2.
\]

Only closed cyclic Cartan products (b_C), not individual Cartan entries,
are required to be algebraic integers.  Symmetrizability gives

\[
 b_C^2=\prod_{(i,j)\in C}c_{ij},
 \qquad |\sigma(b_C)|\le62^{30}.
\]

Nikulin's full-rank field inheritance and Vinberg's finite-index invariance
identify the field generated by the simple frame cycles with the original
minimal Vinberg field.  Each generator has degree at most (25), and its
monic minimal-polynomial coefficients satisfy the finite Vieta bounds

\[
             |q_j|\le {25\choose j}(62^{30})^j.
\]

There are only

\[
 \sum_{r=2}^{30}{30!\over(30-r)!r}
 =24893912605687593731774060679329
\]

ordered simple-cycle slots modulo cyclic rotation.  Hence only finitely many
cycle tuples, fields, and native pairs occur.  The independent replay
recomputes the complete Vieta vector and cycle-slot count and cross-checks
both the classical owner and its hostile archive without importing either
producer.

## 6. From groups to unmarked chambers

There are two terminal routes.

- ABSW Theorem 6.1 directly gives finitely many conjugacy classes of
  arithmetic maximal reflection groups; Nikulin 2007 Theorem 2.2 is a second
  published terminal theorem.
- In the factorized route, finite native pairs and Nikulin 1982 Theorem
  5.2.1 give finitely many reflective-lattice similarity hosts.  Reflection
  maximality forces (W(P)=R(L)).

The group-to-chamber bridge loses no classes.  A discrete reflection
overgroup of a finite-covolume group has finite index, so it is commensurable
and arithmetic; absolute and arithmetic-only reflection maximality therefore
agree here.  The complete mirror arrangement is locally finite.  A generic
compact path between regular points crosses finitely many mirrors, and the
successive reflections form a gallery.  Thus all complete chambers form one
group orbit.  Conjugacy of groups is consequently the same finite counting
problem as unmarked normalized complete chambers.

The executable includes a finite (A_2) gallery calibration: the two simple
reflections generate six matrices, and the orbit of one regular chamber point
meets all six sign chambers exactly once.  This calibrates the orbit
conversion; the general locally finite gallery statement remains the
mathematical argument above.

## 7. Provenance separation

The proof graph has two incoming branches and no forbidden arrow:

\[
\begin{array}{c}
 \text{admissible trace datum}\to\text{weak Chow carrier}
 \to b\le4R-3,\\[2mm]
 \text{classical }(n,d)\text{ bounds}\to R\le750
 \to b\le2997,
\end{array}
\]

while

\[
 \text{classical direct theorem}
 \quad\text{or}\quad
 \bigl(\text{bounded cycles}\to\text{finite native pairs}
 \to\text{Nikulin fixed data}\to\text{maximal chamber orbit}\bigr)
 \to\text{finiteness}.
\]

There is no path from weak-carrier existence to field, dimension, native-pair,
host, or chamber finiteness.  The uniform number (2997) is a joint
corollary, not reverse leverage.  The machine audit encodes this dependency
graph and rejects a weak-carrier-to-finiteness edge.

The geometric limitation is correspondingly explicit.  The (D_i) need not
be effective, irreducible, Mori-extremal, or complete nef/movable walls, and
the Chow correspondences need not be graphs of automorphisms or birational
maps.  Blow-up Picard rank and rational matrix correspondences are deliberately
unbounded.  No positive AG-only ABSW proof follows.

## 8. Primary-source boundary

The audit independently rehashed the locally cached primary bytes and read
the following exact loci:

| source | SHA-256 | locus |
|---|---|---|
| ABSW 2008 | `974ac7c5141ca1c45b515f145215ad115aba867e42551a9b103233ea52930ee4` | Introduction and Theorem 6.1 |
| Vinberg 1971 | `16a4f497688e6f968b4c63d81a6b049b3f096f1cd4f43e458b2882a5c44ded38` | least definition field/ring and cyclic products, Theorems 1, 3, 5 |
| Vinberg expanded dimension paper | `fd25f8f892eed27e4303a8f3252fd384d8954587a6c1c67f1efed2c38aabe44f` | compact and arithmetic noncompact exclusions in dimension at least 30 |
| Nikulin 2007 | `9a820fef0b88158b283057f7d97acc89130610a9a4aacc6c53a8215e23ce874a` | Theorem 2.2 |
| Nikulin 1982 | `696652b6a67c6ed420de6c1cb61d9422e3619e28e081a78131a5323ce2edee9d` | Theorem 2.1.3(2) and Theorem 5.2.1, including (A=62) |
| Nikulin 2011 | `716e6e7333efbfaa22060b60dc97e8303363a4914409ccdd56ef7a90c923c7af` | Theorem 3.2 and the dimension-transition bounds |
| Belolipetsky--Linowitz v2 | `026e0f635f71e214807e0b505bf7297673eb9e0d47a6a91b7bd6dec134eab97f` | Corollary 5.1, degree at most 25 |

The replay treats temporary PDFs as optional audit-time pins, so a clean
checkout does not depend on `/tmp` or `/private/tmp`.  It validates their
hashes whenever they are present.  The executable certifies the composition
and arithmetic interfaces; it does not machine-prove the imported
literature.

## 9. Reproduction

- [independent script](../scripts/audit_unconditional_ns_chow_encoding_and_maximal_chamber_finiteness.py)
- [frozen archive](../data/enumerations/unconditional_ns_chow_encoding_and_maximal_chamber_finiteness_independent_hostile_audit.json)
- [focused tests](../tests/test_unconditional_ns_chow_encoding_and_maximal_chamber_finiteness_independent_hostile_audit.py)

Run

```text
python3 -m scripts.audit_unconditional_ns_chow_encoding_and_maximal_chamber_finiteness --check-output
python3 -m pytest -q tests/test_unconditional_ns_chow_encoding_and_maximal_chamber_finiteness_independent_hostile_audit.py
```

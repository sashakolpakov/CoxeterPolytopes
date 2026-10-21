# Hostile audit: canonical-orthogonal and K3-cover dimension theorem

**Verdict (2026-08-10): PASS after scope, map-domain, provenance, and TeX
repairs.** The rank theorem is valid, the complex strong-Coble corollary is
genuinely class-wide, and the thirteen frozen carriers satisfy the stated
rank hypothesis. The result is not a universal arithmetic-Coxeter dimension
theorem: existence of the canonical-orthogonal carrier remains an explicit
hypothesis. No Nikulin--Vinberg dimension theorem, ABSW, Prasad, or volume
formula enters this audit.

The audited object is the
[canonical-orthogonal/K3-cover theorem](canonical_orthogonal_k3_cover_dimension_theorem.md),
with its
[ledger](../data/enumerations/canonical_orthogonal_k3_cover_calibration.json),
[certifier](../scripts/certify_canonical_orthogonal_k3_cover_bound.py), and
[tests](../tests/test_canonical_orthogonal_k3_cover_bound.py).

## 1. Rank theorem

Let \(T=\operatorname {Res}_{k/\mathbf Q}V\), where \(V\) has dimension
\(n+1\) over a totally real field of degree \(d\), one Lorentzian real
place, and negative-definite conjugates. Then

\[
 \dim_{\mathbf Q}T=d(n+1),\qquad
 \operatorname {sig}(T)=(1,d(n+1)-1).
\]

An injective positive similarity
\(T\hookrightarrow\operatorname {Num}(X)_{\mathbf Q}\) consumes the unique
positive direction supplied by Hodge index. Its orthogonal complement is
negative definite. If its image is orthogonal to \(K_X\), minimality gives
the exhaustive alternative

    kappa=-infinity:  rho(X)<=2, excluded by d(n+1)>=3;
    kappa>=0:         K_X is nef, K_X^2>=0, hence K_X=0 in Num(X).

Surface abundance/classification then puts \(X\) in Kodaira dimension zero.
The numerical-rank maxima are \(20,10,4,2\) for complex K3, Enriques,
abelian, and bielliptic surfaces, respectively. In positive characteristic
the largest unconditional value is the K3 value \(22\). Thus

\[
 d(n+1)\leq20\quad(\mathbf C),\qquad
 d(n+1)\leq22\quad(\operatorname {char}>0).
\]

The restriction \(n\geq2\) is essential to the short exclusion of minimal
ruled surfaces. The positive-characteristic statement is over an
algebraically closed field. The machine certificate checks the rank
arithmetic; Hodge index, abundance, and surface classification remain
declared theoretical inputs rather than simulated computations.

## 2. Strong Coble cover

The words **strong Coble** must include that \(S\) is rational. For

\[
 |-K_S|=\varnothing,\qquad
 |-2K_S|=\{B_1+\cdots+B_b\},\qquad B_i\cap B_j=\varnothing,
\]

adjunction gives \(B_i^2=-4\), hence \(K_S^2=-b\). Rationality and Noether's
formula, not the linear-system conditions alone, give

\[
                         \rho(S)=10+b.                 \tag{2.1}
\]

In characteristic different from two, the cover defined by \(L=-K_S\) and
branched over \(B\in|2L|\) is smooth and satisfies

\[
 K_X=\pi^*(K_S+L)=0,\qquad
 \pi_*\mathcal O_X=\mathcal O_S\oplus\mathcal O_S(K_S).
\]

Rationality gives \(\chi(\mathcal O_X)=2\); nonempty reduced branch gives
connectedness; Serre duality then gives \(q(X)=0\). Thus \(X\) is a K3
surface. This is valid for every strong Coble surface in the stated
characteristics, not just a finite-automorphism class.

The original report incorrectly typed the extended pullback as a map with
domain \(N=\operatorname {Pic}(S)\) and then evaluated it on
\(B_i/2\notin N\). The repaired map is

\[
 \jmath:\widetilde N=\langle N,B_1/2,\ldots,B_b/2\rangle
       \longrightarrow\operatorname {NS}(X),\qquad
 \jmath|_N=\pi^*,\quad\jmath(B_i/2)=R_i.
\]

It is the restriction of rational pullback. The relations
\(\pi^*B_i=2R_i\) make its image integral, \(\pi_*\pi^*=2\) proves
injectivity, and projection formula gives

\[
             \widetilde N(2)\hookrightarrow\operatorname {NS}(X) \tag{2.2}
\]

of rank \(10+b\). Hence the cover argument itself proves \(b\leq10\) and
\(\dim\mathbf P(\operatorname {Nef}S)\leq19\) over \(\mathbf C\). It proves
only \(b\leq12\), hence the weaker ceiling \(21\), from K3 cohomology in
positive characteristic.

The latter is not the best class-wide Coble bound. Kondō records the general
theorem \(b\leq10\), in the characteristic-\(\ne2\) context, citing
Dolgachev--Kondō, Corollary 9.1.5. Thus the best known class-wide ceiling is
\(19\) in every characteristic different from two. The finite-automorphism
classification improves it to \(17\) only on that locus. The standalone
report and ledger now keep these statements separate:

| scope | boundary bound | nef projective dimension |
|---|---:|---:|
| canonical-cover cohomology over \(\mathbf C\) | \(b\leq10\) | \(\leq19\) |
| canonical-cover cohomology in characteristic \(p\ne2\) | \(b\leq12\) | \(\leq21\) |
| best general strong-Coble theorem, \(p\ne2\) | \(b\leq10\) | \(\leq19\) |
| finite-automorphism strong-Coble locus | \(b\leq8\) | \(\leq17\) |

The class-wide corollary is therefore real, but it is a new project-native
derivation over \(\mathbf C\), not a worldwide-new numerical theorem.

Primary source: Shigeyuki Kondō,
[*Coble surfaces with finite automorphism group*](https://arxiv.org/abs/2102.00591),
especially the definition and canonical cover in the introduction, the
general \(b\leq10\) statement in Section 3, and the nine-class table.

## 3. Facet correspondences and monodromy

All three Coble facet mechanisms survive direct calculation.

1. \(B_i/2\) pulls back to the ramification root \(R_i\), with \(R_i^2=-2\).
2. A physical \((-1)\)-curve \(E\) pulls back to a divisor root
   \(e=\pi^*E\), with \(e^2=-2\).
3. A boundary-disjoint \((-2)\)-curve \(C\simeq\mathbf P^1\) has split
   unramified inverse image \(c+c'\), where \(c,c'\) are orthogonal
   \((-2)\)-roots. On the invariant lattice, \(s_cs_{c'}=s_{c+c'}\).

For any divisor root \(\delta^2=-2\), the codimension-two correspondence

\[
                         \Delta_X+\delta\times\delta
\]

acts on divisor classes by
\(x\mapsto x+(x\mathbin{\cdot}\delta)\delta=s_\delta(x)\). This verifies an
actual algebraic correspondence, not merely an abstract lattice matrix.

The complex monodromy wording is also correct after reading “global
monodromy” in Huybrechts' sense: the group generated by monodromies of all
smooth proper K3 families based at \(X\) is \(O^+(H^2(X,\mathbf Z))\), and a
negative-root reflection belongs to it. Therefore every complex root or
paired-root lift is generatorwise a genuine monodromy element. This does
**not** produce one archived family, nor does it force independently selected
generator lifts to satisfy the Coxeter relations on the anti-invariant
complement. The ledger correctly records those rows as generatorwise only.

Primary source: Daniel Huybrechts,
[*Lectures on K3 Surfaces*](https://doi.org/10.1017/CBO9781316594193),
Chapter 7, Proposition 5.5.

In characteristic three the correspondences remain valid on the numerical
lattice, but no crystalline or étale Picard--Lefschetz family is certified.
The report correctly refuses to call them monodromy.

## 4. Thirteen-carrier calibration

The count is \(1+1+4+7=13\): Wehler, Bugaenko, four Coble rows, and seven
Enriques types. The calibrated trace ranks are

    Wehler       3
    Bugaenko    18 = 2*(8+1)
    Coble       11, 12, 15, 18 for b=1,2,5,8
    Enriques    10 for each of seven H9 chambers.

The action-level audit passes:

| row | strongest certified statement | hostile limit |
|---|---|---|
| Wehler | three actual K3 automorphisms | no higher-dimensional consequence |
| Bugaenko | faithful coherent action in complex K3 monodromy, via pairs of orthogonal roots | not surface automorphisms; moves nef chambers |
| two complex Coble rows | algebraic correspondences and generatorwise K3 monodromy | no coherent Coxeter lift on full \(H^2\) archived |
| two characteristic-three Coble rows | algebraic correspondences on supersingular K3 covers | no positive-characteristic monodromy claim |
| seven complex Enriques rows | paired-root correspondences and generatorwise K3 monodromy upstairs | no coherent lift archived |

The ledger hash-pins every direct frozen carrier input. It now separately
lists the theoretical sources and marks surface classification as not
machine-verified. Passing the certifier proves internal arithmetic and
provenance consistency; it is not evidence for an external classification
theorem beyond the cited source.

## 5. Repaired defects and replay

The hostile pass required four corrections:

1. remove three literal carriage-return bytes introduced by malformed rho
   strings and repair the pervasive broken inline TeX;
2. state rationality explicitly in the strong-Coble definition and identify
   where Noether's formula is used;
3. change the extended pullback domain from \(N\) to \(\widetilde N\);
4. distinguish the positive-characteristic K3-cohomological ceiling \(21\)
   from the best general strong-Coble ceiling \(19\), and remove the
   ambiguous claim that the numerical Coble theorem itself is new.

Replay:

    python3 -m py_compile scripts/certify_canonical_orthogonal_k3_cover_bound.py \
      tests/test_canonical_orthogonal_k3_cover_bound.py
    python3 scripts/certify_canonical_orthogonal_k3_cover_bound.py
    python3 -m pytest -q tests/test_canonical_orthogonal_k3_cover_bound.py

The result is PASS and 7 passed. The test suite now rejects C0 control bytes,
the old malformed TeX token, the wrong pullback domain, and any renewed
conflation of the two positive-characteristic ceilings. JSON parsing and a
Pandoc render also pass; Pandoc emits only its expected warnings for boxed
equations, equation tags, and the display array.

## 6. Exact leverage limit

The theorem proves a dimension bound **after** a whole restriction-of-scalars
trace block has entered one minimal canonical-orthogonal surface. It does not
prove any of the following arrows:

    raw maximal arithmetic Coxeter polytope
        -> bounded AG packet
        -> minimal K-trivial carrier
        -> full trace injection.

Maximality is not translated into any of those properties. Picard rank also
does not bound the defining field or its discriminant. The integral
Gram/interval packet remains necessary for field effectiveness. The Coble
lane is necessarily a special noncompact surface lane and cannot by itself
cover compact polytopes. Consequently this is a sound carrier theorem and a
useful target for the polygon-to-surface compiler, but not yet the requested
ABSW-type theorem.

# Independent hostile audit: Coxeter finiteness constant improvements

**Date:** 2026-08-12  
**Verdict:** **PASS after repair of the contraction case split.**

I audited the frozen owner package

| artifact | SHA-256 |
|---|---|
| owner report | `4d54eaeb624fbf5ad8fa6f99533c15a04dc473233291202738351a7603c53ffb` |
| owner certificate | `fca46cad8ab441c2fd7633e73e6981d353ebb58262701d69c40a6290b7e06593` |
| owner tests | `eae1ffcb12bac2e27e29e4e7826812c94a2bbc25141aa4b44c84a4c5ac425366` |

The first owner draft incorrectly said that both displayed norm ratios were
always contracting.  Counterexamples are

\[
 R_{3,3}=9/4,
 \qquad R_{3,4}=3/2,
 \qquad T_3=9/4,
 \qquad T_4=1,
 \qquad T_5=25/16.
\]

The frozen report repairs this by splitting into \(R<1\) and \(R\geq1\)
(and similarly for \(T\)).  In the contracting case, \(e\geq1\) permits
replacing \(R^e\) by \(R\).  In the noncontracting case, the proposed
necessary sieve inequality is automatic.  This repair preserves every
numeric conclusion and removes the only mathematical defect found.

## Independent checks

1. **Real cyclotomic compositum degree.**  For
   \(E=\mathbf Q(s_m,s_n)\), \(s_r=4\sin^2(\pi/r)\), direct enumeration of
   the residue classes \(a\in(\mathbf Z/\ell)^\times\),
   \(\ell=\operatorname {lcm}(m,n)\), confirms

   \[
   [E:\mathbf Q]=
   \begin{cases}
   \varphi(\ell)/2,&\gcd(m,n)>2,\\
   \varphi(\ell)/4,&\gcd(m,n)\leq2.
   \end{cases}
   \]

   I replayed the residue-kernel calculation independently for every
   \(3\leq m\leq n\leq120\), including all four sign-pair compatibility
   regimes.  I also verified \(d_m\mid D\) and \(d_n\mid D\) throughout the
   owner's complete \(3\leq m\leq n\leq840\) box.

2. **Individual label bridge.**  For \(841\leq M<10^6\), the owner checks
   \(3h(M)/4<1\) and excludes the necessary inequality

   \[
      4^{\varphi(M)/2}
      <4M^4\gamma(M)3^{\varphi(M)/2}.
   \]

   Its binary comparison is strict because
   \((4/3)^{10}>2^4\).  An independent stronger grouping,
   \((4/3)^{15}>2^6\), also excludes the entire interval without any exact
   fallback.  Direct factorization confirms that the labels below \(10^6\)
   which pass this one-label envelope have maximum \(840\).  The previously
   audited analytic theorem supplies the separate \(M\geq10^6\) tail.

3. **Fork-pair sieve.**  I independently implemented trial-division
   \(\varphi\) and \(\gamma\), the above degree formula, and exact integer
   comparison with

   \[
   4m^2n^2\gamma(m)^{D/d_m}\gamma(n)^{D/d_n}>4^D.
   \]

   It returns exactly **525** unordered pairs.  The largest surviving label
   is **420**, and the unique endpoint pair is
   \((m,n,D)=(3,420,48)\).  The SHA-256 of the independently ordered Python
   tuple list is
   `3b36049ea469acbf9fb34773e40dee5452ad4c204b1b730e5a0d6bac5557d406`.

4. **Tail sieve.**  Independent exact comparison of

   \[
        601d^4\gamma(d)^2>4^{\varphi(d)/2}
   \]

   for \(3\leq d<20000\) returns exactly **66** labels and maximum **126**.
   The owner's list agrees entry for entry.  The previously audited analytic
   endpoint argument supplies the separate \(d\geq20000\) tail.

5. **Scope.**  All filters are necessary conditions only.  The package does
   not promote any of the 525 pairs or 66 labels to an admissible field,
   integral determinant packet, Coxeter polyhedron, or new example.

The owner tests complete with `5 passed`.  No floating approximation is used
in the finite sieves.  The repaired constants

\[
                   \max(m,n)\leq420,
                   \qquad d\leq126
\]

are therefore certified consequences of the already proved norm envelopes.

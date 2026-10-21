# Hostile scope audit: exact maximal K3 eigen-intervals

**Verdict (2026-08-10): PASS for \(n=1\) and for maximality in the declared
\(\mathbf Z[t]\) full-trace host.  No \(n\ge2\) or overorder-maximal theorem
has been proved.**

## Frozen owner inputs

| artifact | SHA-256 |
|---|---|
| results/rank4_k3_exact_maximal_eigen_interval_infinite_fields.md | 9351875964fbc73da81b55ec3ff68db54af3dbe1ec90272555cff1938ab13c95 |
| scripts/certify_rank4_k3_exact_maximal_eigen_interval_infinite_fields.py | 1b1e6c60ef6e7987adac7768dc9fdf87345482e6986d8966618c10c1a5eaa48d |
| data/enumerations/rank4_k3_exact_maximal_eigen_interval_infinite_fields.json | a238ea60701c86d9c2c688ee6fa889b73c9953632a5ef1f344e2ef5441b900cd |
| tests/test_rank4_k3_exact_maximal_eigen_interval_infinite_fields.py | 9cf8b5c551764c9f456c2d204eee6524d2330c2640f7c28b13b6032913adad69 |

## Hostile attack table

| gate | attack | verdict |
|---|---|---|
| lattice | the change \((x_1,x_2,y_1,y_2)=(X,Y,Z,U-Z)\) must be unimodular and produce \(G=-2\operatorname {diag}(-1,1,1,2d)\) | PASS |
| eigenspace | the companion eigenvector must retain both real places, with signature \((1,1)\) at \(t\) and definite conjugate form | PASS |
| endpoint roots | \((0,0,-1,0)\) and \((1,-1,1,0)\) must have \(q_d\)-norm one and cut exactly \(s=0,1\) | PASS |
| hidden root | every sign-normalized root is strictly negative at \(p_0\), and weakly negative at \(p_1\) with equality only at the second endpoint | PASS |
| integral reflection | both binary reflection matrices square to one, preserve \(h_a\), lift integrally to \(G_a\), and commute with \(\tau_a\) | PASS |
| saturation | the dihedral translates tile all of \(\mathbf H^1\); an extra restricted K3 mirror would translate into the certified root-free open interval | PASS |
| host maximality | an interior integral reflection forces \(1<P<2t-1\), \(|P'|\le1\); the only trace candidates \(t-1,t\) fail their exact square equations | PASS |
| finite volume and field | the length is finite, \(\cosh^2\ell=t\), and the cyclic Cartan product \(4t\) recovers \(k\) | PASS |
| exact K3 return | the determinant-one primitive embedding, generic exact-Picard period, common Weyl framing, and divisor-product operators return the same marked interval | PASS |
| bounded polarization | \((3,-1,1,0)\) is strictly interior, primitive, and has square \(14\), while its \(\tau\)-height is \(32a-10\) | PASS |
| ambient contrast | the roots \((dk^2,0,dk^2-1,k)\) give infinitely many ambient facets; the eigenslice theorem is not ambient 2-reflectivity | PASS |

## The precise maximality statement

The positive theorem concerns the reflection host

\[
              (A^2,h_a),\qquad A=\mathbf Z[t].
\]

Every competing reflection must be \(A\)-linear and preserve this lattice.
This is a literal host-relative maximality statement, not shorthand for any
of the following:

- maximality after replacing \(A\) by its integral closure;
- maximality after passing to a commensurable trace overlattice;
- maximality inside the full rational orthogonal group \(O(L_a)\); or
- a dimension-\(\ge2\) arithmetic reflection theorem.

The full rational lattice actually has explicit extra long reflections, so
promoting the stated result to full-\(O(L_a)\) maximality would be false.

## Consequence and surviving obstruction

If a proposed field-finiteness theorem admits \(n=1\) and maximality relative
to a declared integral host, this family refutes it.  Even adding a fixed
degree-\(14\) ample polarization does not help.  The escaping integral data
are

\[
 |\det L_a|=32(a-1),\qquad
 H_a.\tau_a(H_a)=32a-10,
\]

and the polynomial/field discriminants.  Maximality controls mirror
subdivision but does not normalize the trace order, the Picard discriminant,
or the correspondence height.

The \(n\ge2\) problem must therefore either couple additional Lorentzian
directions while preserving finite volume and host maximality, or derive a
new higher-dimensional rigidity unavailable on the hyperbolic line.

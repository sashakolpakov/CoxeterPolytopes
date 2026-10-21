# Hostile scope audit: orthogonal K3 stabilization obstruction

**Verdict (2026-08-11): PASS, at the stated scope.**  The rank-six owner
artifact proves an explicit compact distinguished eigenslice on a primitive
K3 carrier, but its required mixed closing mirror has a nonintegral cyclic
Cartan product.  The same pairwise obstruction is valid in every higher
orthogonal trace stabilization.  The report correctly does **not** promote
this to an obstruction for nonorthogonal couplings, nor does it claim a
uniform higher-rank K3 embedding or compact chamber.

## Frozen owner bytes

- report: `a7d71ac262bdb709db14dafb5488ef00ce343366084783d427fef971cd6b0e01`;
- certifier: `fd0ff5530fd50117a9e29151d4487766caf05c56d59cd9dfb16f6b915567a71c`;
- archive: `b979263830c876508b69359f161c738875b597ae14418e0cd02cdb6ffe3f117b`;
- owner tests: `bca7ea7dc513d4d2d6098286b4b57eefc27ca53151d9538d568ff4413043a8ac`.

## Independent hostile checks

The standalone hostile test imports no owner arithmetic.  It reconstructs:

1. the unimodular diagonalization of
   \(2B_0\oplus2B_a\oplus2B_a\), self-adjointness of \(C_a^{\oplus3}\),
   and its quadratic field relation;
2. the determinant-one minor in the primitive K3 embedding;
3. all five norm-one roots and the strict Klein-disk location of the visible
   pentagon vertices;
4. the exact sign of \(B-tA\) on a large normalized root shell, including an
   exhaustive height-one replay;
5. the identity
   \[
      \chi_d={4(2d-1+t)\over4d-1},\qquad
      N(\chi_d)={32d\over4d-1}\notin\mathbf Z;
   \]
6. the identical norm-one mixed root in every coordinate pair of every
   orthogonal stabilization tested; and
7. exact archive replay and all nonclaim language.

The passage from nonintegral \(\chi_d\) to failure on any commensurable
restriction-of-scalars lattice is sound.  The product has eigenvalues
\(1,\lambda,\lambda^{-1}\); preservation of a \(\mathbf Z\)-lattice would
make both \(\lambda\) and \(\lambda^{-1}\) algebraic integers, forcing
\(\chi_d=\lambda+\lambda^{-1}+2\) to be integral.  This avoids the invalid
shortcut of using only the rational trace.

## Surviving problem

The no-go exhausts products and orthogonal duplications of the interval
packet.  A genuinely nonorthogonal field form could change the mixed normal
norm and repair the cyclic product.  It remains viable in principle, but it
must simultaneously preserve a primitive K3 trace embedding, a finite-volume
complete eigenslice, integral saturation, and reflection maximality.  No
classification theorem is imported and no result for all such couplings is
claimed.

The independent replay is
[here](../tests/test_rank6_k3_orthogonal_stabilization_nonarithmetic_no_go_hostile_scope.py).

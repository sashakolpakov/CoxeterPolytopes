# Outside review: quartic order exchange

20 October 2026. Full paper/source audit; no scientific execution,
worker import, Lean or TeX run. Verdict: **PASS within the literal-lattice
scope and the explicitly named arithmetic inputs.**

Reviewed source:

* `compact_biquadratic_unit_order_exchange_obstruction_20261020.md`,
  SHA256 `043f1f1f6605f6e0b80c1051fd6d288f119e6fbb65487af45d092e1669320299`.
* Complete first edge `compact_h9_biquadratic_unit_diagonal_first_edge_20261020.md`,
  SHA256 `05b6473eeadb4befa89c82af2f7002c7441ada3d1878c93d4630355a42685efe`.
* Complete second edge `compact_h9_biquadratic_unit_second_edge_20261020.md`,
  SHA256 `b3a7c8f4e2afffe13fb49dd80ee2ed03095608003fabea85e2e6b74c9df2a52b`.

The reviewer previously supplied an independent crossing-product strip
and independently checked the second finite upper. The order-exchange
construction is the author's/root's; this is an outside full review,
not a claim of three independent discoveries.

## Arithmetic and complete actual edges

I rechecked discriminant1600, PID via the ideal bound15/4 and absence
of ideals of norms2/3, the four independent unit signatures, and the
primitive pool1/2/d. The trace forms give the full B_n zero germ.
The first source's trace8/16 families and final four-embedding filters
are complete. The second source retains all six finite roles, the
complete trace strips and the dyadic congruence. Its q2 isolated row
is imprimitive. The accepted c, r+c, alpha*r+s*c are all retained,
and the latter two are positive combinations, so c alone is simple.
The finite-upper identity d(a+1)^2-a*phi^2=k^2/phi^2 is exact.

## Full faces and oriented chamber germs

Both P and Q are actual faces, not abstract positive diagrams. For Q,
projecting the old vertex along the retained positive-definite Gram
uses nonnegative coefficients; acuteness keeps every actual wall
inequality. All previous nonincident walls stay strict, and the only
released old wall beta2 becomes strict. This gives the actual whole
corner Z=(1;1,1,0...), with norm2-a<0.

The displayed Q basis has determinant -1/(alpha*phi), and its Gram
determinant is a(a-2)=(alpha*phi)^-2. Thus projection of any integral
ambient vector to Q is integral. A nonorthogonal crossing mirror
would generate a finite irreducible rank-three reflection group
containing a pi/4 pair. Only B3 qualifies. Its appropriate long-root
orbit gives a primitive norm2 root with products(-1,0) or(0,-1),
forcing Q-projection -beta-c/s or -c-beta/s, respectively. Neither
is integral. This checks arbitrary crossing mirrors, not only simple
native extensions. The analogous pi/8 classification seals P.
Intrinsic reflectors extend across the unimodular summands, proving
both FULL lower-arrangement identities.

The five A/B columns have identical block Gram diag(-a,G_P,G_Q),
unit determinants, and future time columns. I checked each positive
simple image, including Q simples -A4 and C*A3+U*A4. Consequently
D=B*A^-1 maps the full oriented germs of two actual vertices of the
same chamber. This proves chamber stabilization; mere Gram matching
would not suffice. The remaining B_(n-4) germ is fixed and cannot mix
with either sealed plane.

## Exact trace, dimension range and conclusion

Independently computing the Gram-adjoint trace gives the three terms
in equation(8). Using I+J=a(H-1) reduces the time term to
alpha*K^3-alpha*K^2+K, while the Q term is alpha*K-K^2.
Together with the P term these give

    alpha*K^3-alpha^2*K^2+(3+s)*K
      =196phi+121+s*(138phi+86).

I expanded K^2 and K^3 independently and obtained the printed values.
The full trace exceeds n+1. Thus the integral future isometry D has
infinite order. Its chamber-stabilizing powers give distinct cosets
modulo the full reflection group. Borel--Harish-Chandra plus the
nonphysical definite embedding supplies the stated uniform arithmetic
lattice; infinite index gives infinite reflection covolume.

The propagation to n>=8 uses the two complete inventories only in
dimensions n and n-2>=6. Retained B4 at the threshold has no exceptional
one-root finite extension, and the final zero germ is dimension
independent. The low-dimensional N2 compact calibration is not used
as a contradictory extrapolation. No claim about arbitrary invariant
lattices, all quartic forms, novelty, execution, or a new dimension
ceiling follows from this audit.

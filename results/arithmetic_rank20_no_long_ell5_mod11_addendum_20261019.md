# Elementary mod11 replacement for the e-label5 norm step

19 October 2026. Root-supplied alternative, checked independently by
coefficient multiplication and residue arithmetic. No scientific
program, build or Lean process ran. The source-bound report remains
unchanged:

`arithmetic_rank20_isolated_side_no_long_contact_20261019.md`,
SHA256 `af2ff3c9c4667b037ea93ed1afa8f1dcaf364c0ac41e894e3541b2bcde1df879`.

Its e-label5,N=8 argument gives an algebraic integer delta in
K15 or K20 with |Norm(delta)|=11. The following replaces the
cyclotomic residue-degree input in that ONE final step. All other
hypotheses and prior norm bounds are inherited unchanged.

The principal quotient O_K/(delta) has11 elements, so it is the
field F11. Every algebraic integer theta in K therefore has an
image in F11 satisfying its integral polynomial.

For K15 take theta=zeta15+zeta15^(-1). Dividing Phi15(zeta)
by zeta^4 and using

    zeta^3+zeta^(-3)=theta^3-3theta,
    zeta^4+zeta^(-4)=theta^4-4theta^2+2

gives

    f15(theta)=0,
    f15(X)=X^4-X^3-4X^2+4X+1.

In F11 this factors as

    f15(X)=(X^2+7X+2)*(X^2+3X+6).

The two discriminants are8 and7 modulo11. The complete square
residues are{0,1,3,4,5,9}, so neither quadratic has a root.
Thus f15 has no root in F11, a contradiction.

For K20 take theta=zeta20+zeta20^(-1). Dividing Phi20(zeta)
by zeta^4 gives

    f20(theta)=0,
    f20(X)=X^4-5X^2+5.

Modulo11,

    f20(X)=(X^2-6)*(X^2-10).

Neither6 nor10 is a square modulo11, so this polynomial also
has no root in F11. This is the second contradiction.

Consequently neither K15 nor K20 has an integral element of
absolute norm11. The e-label5 proof can use this elementary
quotient argument instead of the cyclotomic Frobenius/decomposition
law. It requires neither monogeneity nor a claimed basis for O_K:
the displayed cyclotomic theta merely belongs to the ring of
integers. No additional local or global weight conclusion is
claimed by this addendum.

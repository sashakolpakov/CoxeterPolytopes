# The 28 cusp fibrations over the full rational carrier family

## Positive theorem

Let `U_car` be the nonempty Zariski-open carrier locus from the `j`-moduli
theorem: `4a^3+27` is nonzero, the curves `C_a,Q_a,S_a` are geometrically
integral, and the complete nef chambers `P12,...,P16` persist on the iterated
blow-ups `X_n(a)`.

> **Theorem.** For every `a` in `U_car`, all 28 ideal rays occurring across
> `P12,...,P16` are semiample on `X_n(a)`. Twenty-seven give base-point-free
> pencils in their primitive class. The last class
> `F=4H-E_1-...-E_16` has exact semiample index four:
> `h0(F),h0(2F),h0(3F),h0(4F)=1,1,1,2`, and `|4F|` is base-point-free with
> connected genus-nine general fibre.

There is **no additional cusp-specific exceptional parameter**. In the
language requested by the search plan, the extra exceptional polynomial is

```text
E_cusp(a)=1.
```

This strengthens the two-fibre result at `a=0,1` to the entire
non-isotrivial carrier-open family. It is a new positive geometric
realization theorem in this project, not a claim that the abstract Vinberg
polytopes `P12,...,P16` are new. Worldwide priority for this exact
carrier-and-fibration family remains open.

## The symbolic cluster

In the affine chart `Y=1`, fifteen exact strict transforms give the selected
`S_a` branch, modulo `x^16`, as

```text
phi_S,a(x) = x^3 + a*x^7 + x^9 + (2*a^2-1)*x^11
             + (5*a+1)*x^13
             + (a-1)*(5*a^2+5*a-1)*x^15.
```

Thus every complete-ideal jet matrix is defined coefficientwise over
`ZZ[a]`. For `D=dH-sum(m_i E_i)`, its rows are exactly the negative-power
and low-total-degree coefficients after the successive blow-up substitution

```text
x^-(m_1+...+m_(k-1)) G(x,phi_<k(x)+x^(k-1)u).
```

The new builder evaluates entry-for-entry to the independent v3 two-fibre
jet builder at `a=0` and `a=1`.

## Exact rank charts

For each moving system, deterministic modular pivots choose a square maximal
minor. Its determinant is recomputed exactly in `ZZ[a]`. A polynomial second
section is reconstructed from exact integer fibres and accepted only after
coefficientwise substitution proves that it lies in the full symbolic
kernel.

| system | plane degree and multiplicities | first rank minor | second rank minor |
|---|---|---|---|
| genus 5 | `11;(3^13,2)` | `-(a+3)` | `-3` |
| genus 10 | `21;(7^5,6,4^10)` | `1` | unnecessary |
| genus 11 | `23;(7^8,5,4^7)` | `a-1` | `-7` |
| `4F` | `16;(4^16)` | `-12(a+1)` | `2(a^2-10a-13)` |

The stored second sections have respectively 33, 85, 115, and 65 nonzero
`(x,z)` monomials in the first charts, with coefficient degrees at most 4,
7, 8, and 6 in `a`. The alternate sections have 33, 115, and 65 terms and
parameter degrees at most 5, 9, and 7.

For every section, the certificate also proves:

1. Its successive strict transforms have exactly the assigned multiplicity
   sequence over `QQ(a)`. The gcds of the leading multiplicity coefficients
   are recorded centre by centre.
2. Polynomial division by each relevant wall factor among `L,C_a,Q_a,S_a`
   has nonzero remainder. The gcd of the remainder coefficients is the exact
   specialization gate preventing a common plane component.

Combining rank, independence, multiplicity, and coprimality gives:

| system | first chart gate | second chart gate | gcd |
|---|---|---|---|
| genus 5 | `a+3` | `3a+1` | `1` |
| genus 10 | `1` | -- | `1` |
| genus 11 | `a-1` | `7a-1` | `1` |
| `4F` | `a+1` | `(2a+3)(a^2-10a-13)` | `1` |

The cover is literal, not probabilistic:

```text
3(a+3)-(3a+1)=8,
7(a-1)-(7a-1)=-6,
(2a+3)(a^2-10a-13)-(a+1)(2a^2-19a-37)=-2.
```

Hence at every complex parameter at least one chart supplies a full-rank
minor and a second effective member with exact multiplicities, independent
of and coprime to the wall member. The strict transforms have no common
curve and intersection zero, so they are disjoint and span a base-point-free
pencil.

The genus-five system on `X_14(a)` pulls back to the corresponding cusp
classes on `X_15(a)` and `X_16(a)`. Together with the genus-ten and
genus-eleven systems, this accounts for the five formerly open primitive
cusps. The original 22 pencils use parameter-independent affine wall cycles,
line rulings, and the `P13` quartic construction, so they already apply
throughout `U_car`.

## The exact index-four statement

For `F,2F,3F`, the symbolic jet matrices contain constant nonzero maximal
minors `1,-2,2`. Since `Q_a,Q_a^2,Q_a^3` are universal kernel members, these
minors prove `h0(F)=h0(2F)=h0(3F)=1` for every parameter. The covered `4F`
charts prove `h0(4F)=2` and supply a coprime member with actual
multiplicities `(4^16)`, so `|4F|` is base-point-free everywhere on `U_car`.
The same Stein-factor argument as in the two-fibre theorem makes the general
fibre connected, and adjunction gives arithmetic genus nine.

## Scope boundary

The assertion is conditional only on membership in the already proved
carrier-open locus `U_car`. The factor `4a^3+27` is the explicit smoothness
factor for `C_a`; geometric integrality of `Q_a,S_a` and persistence of the
complete cone are the open conditions from the carrier theorem. This report
does not claim those carrier conditions at every parameter of the affine
line. What it proves is sharper and separate: once a parameter is a carrier,
the cusp-fibration calculation excludes no further parameter.

## Replay

The exact artifacts are:

- [`certify_rational_i1n_symbolic_cusp_family.py`](../scripts/certify_rational_i1n_symbolic_cusp_family.py)
- [`rational_i1n_symbolic_cusp_family.json`](../data/enumerations/rational_i1n_symbolic_cusp_family.json)
- [`test_rational_i1n_symbolic_cusp_family.py`](../tests/test_rational_i1n_symbolic_cusp_family.py)

The full rebuild contains two `298x298` determinants over `ZZ[a]`, so it
should remain under `caffeinate`:

```bash
PYTHONPATH=.:src caffeinate -dimsu python3 -m scripts.certify_rational_i1n_symbolic_cusp_family
PYTHONPATH=.:src python3 -m scripts.certify_rational_i1n_symbolic_cusp_family --check-archive
PYTHONPATH=.:src pytest -q tests/test_rational_i1n_symbolic_cusp_family.py
```

Frozen payload SHA-256:

```text
13975d6140b27dc0f10092e4d3e7769a7cfc7b9a19fc9de5745d2aced92dc8a6
```

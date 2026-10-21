# Hostile audit: universal Gram--Belyi Coxeter dessin

Status: **PASS after repair, with exact scope restrictions** (2026-08-10).

Audited files:

- [`gram_belyi_coxeter_dessin_theorem.md`](gram_belyi_coxeter_dessin_theorem.md),
  SHA-256 `085428e660a94a8926276db5e31dbb26d45216a169ffa6213bf1779bb17ce0c3`;
- [`certify_gram_belyi_h10_lower_bound.py`](../scripts/certify_gram_belyi_h10_lower_bound.py),
  SHA-256 `5d2d5ac188d131cdf5148d2b28d3552c588b5a497c5f85c556db8d5d23586fdc`;
- [`gram_belyi_h10_lower_bound.json`](../data/enumerations/gram_belyi_h10_lower_bound.json),
  SHA-256 `9121aee093b48e93ced3777a807bac323b5fcd990ee7dc1ffc73ece57015f8ca`;
- [`test_gram_belyi_h10_lower_bound.py`](../tests/test_gram_belyi_h10_lower_bound.py),
  SHA-256 `8b7dffd3b9910705706ee6b054b85c22bd6d683c3ea5ac3f66336936c0b1fa6d`;
- [`gram_belyi_standard_packet_calibration.md`](gram_belyi_standard_packet_calibration.md),
  SHA-256 `9afc4f0f52dfed807e7a4b96d17c0896acbbc885b76701f640c812714e82f4c4`;
- [`certify_gram_belyi_standard_packets.py`](../scripts/certify_gram_belyi_standard_packets.py),
  SHA-256 `692024c9ad2385936d75587ebc81638f0b5fee1eae6b42b3744d3d48829cdc4f`;
- [`gram_belyi_standard_packets.json`](../data/enumerations/gram_belyi_standard_packets.json),
  SHA-256 `dfabf3e170c56423603ecb92cb97dac2c2edcea14c7b64f1a088ea53031281e8`;
- [`test_gram_belyi_standard_packets.py`](../tests/test_gram_belyi_standard_packets.py),
  SHA-256 `a21336febc71a61cfb30e868fda8b59f31c238bd150756231ed13b177c0d2a24`.

The first audit found three material defects: the carrier did not explicitly
retain the distinguished real place/time sheet needed for literal return, the
closed field-degree count had an unsupported marking factor and a bad equation
reference, and the report blurred its prose proof with the much narrower
machine certificate. All three are repaired in the audited bytes above.

## 1. Bounded packet finiteness: PASS

Rodriguez's normalization is exactly the one used in the report:
`B(0),B(1)` lie in `{0,1}` and every finite critical value lies in `{0,1}`.
The source frame removes the affine `PGL2` freedom. Riemann existence then
gives only finitely many normalized Belyi polynomials of degree at most `d`.
Each polynomial has finitely many points above zero and one.

Therefore `beta_GB(P)<=B` first gives `f(f-1)<=B`, hence bounded `f` and
dimension, and then a finite list of labelled `lambda` packets and squared
Gram matrices. No arithmetic-volume finiteness theorem is hidden in this
argument. The claim is a genuine conditional finiteness theorem: it does not
provide the still-missing intrinsic bound on `beta_GB`.

Primary references used in the audit:

- [Rodriguez, *Bounding the Degree of Belyi Polynomials*](https://arxiv.org/abs/1104.2027),
  especially Definition 1 and the finiteness statement following it;
- [Khadjavi, *An Effective Version of Belyi's Theorem*](https://doi.org/10.1006/jnth.2001.2748).

### Standard-packet addendum: PASS with one exact minimum open

The added low-label table and equation (1.7) were independently checked
against the frozen exact packet certificate.  For a finite label,
`u=cos^2(pi/m)` and `lambda=u/(1+u)`.  Rodriguez's valuation theorem and the
elementary degree-`b` polynomial

```text
b^b/(a^a (b-a)^(b-a)) * x^a (1-x)^(b-a)
```

give the proved minima

```text
m=3: 5;  m=4: 3;  m=6: 7;  parallel: 2;  dotted u=2: 3.
```

The `m=2` row is correctly labelled as the separately assigned orthogonality
bit-packet, not as an intrinsic application of `lambda=u/(1+u)` at `u=0`.

For `m=5`, exact substitution gives
`lambda=(7+2sqrt(5))/29` and minimal polynomial
`29*x^2-14*x+1`.  Its 29-adic coefficient points are
`(0,0),(1,0),(2,1)`, so the Newton slopes are `0,1` and a root valuation is
`-1`.  Applying Rodriguez after the corresponding coefficient-field
embedding gives `h_Bel(lambda)>=29`.  The separate calibration report's
factored composition `F(Q(C(x)))` was checked at all marked and critical
values and really is a normalized Belyi polynomial.  It supplies only a
finite upper witness; neither the main theorem nor the calibration claims
that the exact `m=5` height is 29.

Finally, every pair contributes the baseline two.  Relative to that baseline,
the audited rows add `3,1,at least 27,5` for labels `3,4,5,6`.  Therefore

```text
beta_GB(P) >= 2*binom(f,2) + 3*N_3 + N_4 + 27*N_5 + 5*N_6
```

is correct.  Other finite labels and dotted coordinates may add further
excess; they are not silently asserted to contribute zero.  Equation (1.7)
fits the existing numbering, between (1.6) and (2.1), with no collision.

## 2. Connected nodal carrier and surface intersection: PASS

For `M=binom(f,2)`, the backbone map `x -> x^M` has `M` distinct points over
one in characteristic zero. Gluing each leaf at `lambda_ij` to a different
one of those points gives an ordinary node and joins every leaf to the same
backbone, so `C_P` is connected and nodal. The component maps agree at each
node. Since the glued curve is proper and the resulting map is quasi-finite,
`beta_P:C_P -> P1` is finite.

Consequently `S_P=C_P x P1` is a projective simple-normal-crossing surface,
and `(beta_P,id)` is finite. Pullback of `O(1,1)` is ample Cartier. On the
fundamental cycle of this reducible surface the projection formula gives

```text
A_P^2 = deg(beta_P,id) * O(1,1)^2
      = 2(M+beta_GB(P)).
```

This is not a claim that the carrier is a stable surface, nor that its field
of definition equals the Vinberg field. The repaired report explicitly makes
the latter distinction.

## 3. Exact Gram return: PASS at the stated geometric scope

The ordered four-point marking recovers each `lambda_ij`; the orthogonality
bit treats the deliberately duplicated value `lambda=1/2`; and
`u=lambda/(1-lambda)` recovers every squared off-diagonal Gram entry.

Squared entries alone do not select an algebraic conjugate, a sign convention,
or a hyperbolic time sheet. The repaired packet now explicitly retains:

- the distinguished real place;
- the acute outward-normal sign convention;
- the intended time sheet.

With those data, `-sqrt(u_ij)` (in the report's conventional Gram sign) and
the known diagonal recover the full unit-normal Gram and labelled halfspaces.
A finite-volume chamber has spanning facet normals, so the Gram determines
the hyperplane configuration up to ambient hyperbolic isometry.

The return is deliberately only geometric. It does not recover primitive
integral root lengths, divisibilities, or the integral reflection lattice;
the repaired report now states this. Such data would require another marking.

## 4. Explicit field-degree bound: PASS after correction

A transitive degree-`d` Belyi cover is represented by an ordered monodromy
pair in `S_d x S_d`, giving the crude upper bound `(d!)^2` before quotienting.
For a polynomial, infinity is its unique totally ramified pole. The three
remaining ordered source marks `0,1,lambda` lie in
`B^{-1}({0,1})`, whose underlying set has at most `2d` points. Hence the
number of source markings is at most `(2d)^3`; replacing this by `(2d)^4` is
a harmless enlargement. The corrected packet bound

```text
product_(u_ij>0) (2 d_ij)^4 (d_ij!)^2
```

therefore bounds the Galois orbit of the labelled `lambda` tuple and the
degree of `Q(u_ij)`. Adjoining all distinguished square roots costs at most
`2^M`, and the Vinberg cyclic-product field is a subfield. Since
`M<=floor(B/2)` and every `d_ij<=B`, the repaired formula

```text
[k_Vin(P):Q]
 <= 2^M product_(u_ij>0) (2 d_ij)^4 (d_ij!)^2
 <= 2^floor(B/2) ((2B)^4 (B!)^2)^floor(B/2)
```

is safe.

The intermediate `(3d)^4` version was **not false**; it was merely looser.
Indeed all four source marks lie in the union of the three special fibres,
which has at most `3d` points, so `(3d)^4` is an immediate universal upper
bound. The polynomial-specific argument above improves `3` to `2` and even
uses only three nontrivial choices. By contrast, the earlier bare `d^4`
factor had no supplied injection/counting proof and was not certificate-grade.

## 5. Machine-evidence boundary: PASS after correction

Sections 1--3 of the main report are a mathematical construction and proof;
they are not verified by either Python producer. The repaired replay section
says so explicitly. The first executable scope is the `H10` calibration:

- exact seed pairings;
- the symbolic formulas for `q_2m`, its square, and its pairing with `c`;
- the rational `p`-adic valuation computation through the archived primes;
- invocation, not reproof, of Rodriguez's general lower-bound theorem.

The second executable scope is the standard-packet calibration:

- exact rational `u -> lambda` transformations;
- exact beta-polynomial values matching Rodriguez lower bounds for the five
  proved rational minima;
- the `m=5` quadratic identity and 29-adic Newton polygon;
- every marked/critical value and the exact degree of the factored `m=5`
  upper witness;
- the explicit warning that the exact `m=5` minimum remains open.

Fresh replay:

```text
python3 scripts/certify_gram_belyi_h10_lower_bound.py --check-archive
python3 scripts/certify_gram_belyi_standard_packets.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_gram_belyi_h10_lower_bound.py \
  tests/test_gram_belyi_standard_packets.py

7 passed in 0.27s
```

## Final verdict

The repaired theorem is internally coherent and no longer overclaims what
the machine artifact proves. Its substantive output is conditional but real:

```text
arithmetic Coxeter chamber
  -> finite marked SNC Belyi-surface packet
  -> same geometric unit Gram and Vinberg field;

bounded absolute Gram--Belyi degree
  -> effective finite chamber and field list.
```

It is not yet the desired unconditional maximal-chamber finiteness theorem.
The open bridge remains a natural intrinsic geometric hypothesis that bounds
the absolute marked degree. No ABSW/Prasad input is used in the conditional
finiteness implication.

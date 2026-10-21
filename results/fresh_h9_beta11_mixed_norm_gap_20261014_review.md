# Independent review of the universal mixed-norm gap

**PASS.** The [primary argument](fresh_h9_beta11_mixed_norm_gap_20261014.md)
correctly proves that reflecting hyperplanes in the two normalized
primitive norm classes, 2 and \(2\beta\), are orthogonal or strictly
disjoint. It covers all heights and excludes both nonright intersections
and parallelism. It applies to reflections preserving this particular
integral lattice, not to arbitrary geometric caps.

The [small independent checker](../scripts/fresh_h9_beta11_mixed_norm_gap_20261014_review.py)
and [certificate](fresh_h9_beta11_mixed_norm_gap_20261014_review.json) use
standalone Fraction arithmetic. They bind the independently reviewed
[primitive norm alphabet and signatures](fresh_h9_h8_cap_height_20261014_norm_pool_review.md),
which are inherited rather than reconstructed again.

The integrality step is valid because the longer root is primitive. The
prior Bezout argument gives \(Gs\in\beta\mathcal O^{10}\) for
\(s^2=2\beta\); hence \((r,s)=\beta z\) with \(z\in\mathcal O\).
When the pairing is nonzero,

\[
 a=\frac{(r,s)^2}{\beta}=\beta z^2,
 \qquad N(a)=11N(z)^2.
\]

The two root lines cannot coincide: that would make \(\beta\) a square
in the field, whereas its norm is the nonsquare rational number eleven.
They remain linearly independent under the second embedding, since any
proportionality ratio can be read from a nonzero coordinate and then
transported back. Strict Cauchy in the conjugate positive-definite form
therefore gives \(0<a'<4\).

A nonorthogonal intersection or parallelism at the hyperbolic embedding
would give \(0<a\le4\). Thus \(0<N(a)<16\), so the nonzero integer
\(N(z)\) is \(\pm1\) and \(N(a)=11\). The positive integer trace is
less than eight and has discriminant \(T^2-44\ge0\); hence \(T=7\).
The independent checker solves this explicitly in integral coefficients:
the only candidates are \(3+\phi\) and \(4-\phi\). Each has an embedding
strictly above four. This contradicts the two simultaneous bounds,
including the proposed parallel equality.

The universal deduction is a proof review; the finite norm, trace and
coefficient exclusions are independently replayed. No Coxeter-angle
catalogue, label cutoff, bounded root search or reflectivity assumption is
used. No mathematical issue was found.

Generation and saved replay, source manifests, syntax and local links
pass. The primary and earlier frozen files were not modified.

```sh
python scripts/fresh_h9_beta11_mixed_norm_gap_20261014_review.py --verify
```

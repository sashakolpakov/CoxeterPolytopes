# A global envelope for vertex-centered H4 caps

This elementary lemma makes a finite global screen possible for a particular
infinite cap family. It does not exclude arbitrary caps or compact H9 polytopes.

Let v be a future unit timelike vertex of the compact H5 seed. Let the four
norm-two H4 simple roots be s_i, and let omega_i be their dual vectors,
`(omega_i,s_j)=delta_ij`. Write `B=Gram(H4)^(-1)`, as exactly verified in the
[standalone transverse gate](fresh_h9_h4_transverse_gate_20261013.md).
For a label `m>=3`, put

```text
t = 2 cos(pi/m),
gamma(t) = sqrt(B_ii - 2/t²),
r_i,m = t (gamma(t) v - omega_i).
```

The old H5 space is orthogonal to H4. Therefore r has square two, pairs
`−2cos(pi/m)` with s_i and zero with the other three transverse walls, and
is orthogonal to the five old walls through v. Along the released transverse
edge, its intersection is

```text
x = v - gamma(t)/B_ii * omega_i,
x² = -2/(t² B_ii) < 0.
```

This is a finite potential endpoint. The two old seed walls not through v
still require a Coxeter-angle check, and coexistence with other caps requires
all their pair checks. This formula by itself does not produce a compatible
completion.

For every `m>=3`, `t>=1` and gamma(t) is increasing. Suppose E is a proposed
escape satisfying `(v,E_old)<0`. Then

```text
(r_i,m/t,E)
 = gamma(t) (v,E_old) - (omega_i,E_H4)
 <= gamma(1) (v,E_old) - (omega_i,E_H4)
 = (r_i,3,E).
```

Thus the label-three cap is the strongest member of this ray family on that
region. A common spacelike escape for all forty label-three caps, with the
ten required old-vertex pairings negative, would exclude every label on all
forty vertex/dual-ray families. One may test this relaxation before any
compatibility filtering. The inequality even covers the entire continuous
range `1<=t<=2`, including the parallel-label limit. This file states the
reduction; it does not assert that the required common escape exists.

The compact H5 seed ensures that any nonzero vector in its future projective
cone pairs negatively with every v above. For a certificate it is simpler
to check those ten pairings directly at the saved candidate escape.

# Outside-author review of the D8 three-direction estimate

20 October 2026. **PASS** relative to the stated classical classification
and pinned actual-edge inputs. Full source reviewed:
`results/compact_d8_three_direction_component_energy_20261020.md`, SHA256
`7555245ab9ee941efbd3a6f7657a7a103d4c9c22f6ef0a63d1330c66eb9ab860`.

I did not coauthor this D8 proof. I independently rederived the contact
inventory, inverse entries and analytic argument below. Earlier work on
related component/actual-edge inputs is not presented as a second proof
of those inputs here. No scientific job, worker import, optimization,
Lean or TeX execution was performed; C28 frozen sources were untouched.

## Common opposites and complete old contacts

Deleting 4 leaves D4 on 5,6,7,8. Its positive-definite one-node
extensions permit no contact or just a 3 at one leaf. A middle 3 gives
the degree-four affine star, a leaf 4 already has zero Schur residual,
and two contacts give a nonpositive cycle. Deleting 5 leaves A4 on
1,2,3,4, so there is at most one contact to that left chain.

I checked each remaining branch independently. A right 3 at 5 forces
contact at 4 by the exact affine two-fork diagram; the A4 replacement
permits only 3 or 4 there. The 3 makes the all-3 affine triangle, so
the contact is 4. Any spectator contact gives, in the finite replacement
missing 5, the nonpositive five-node path 3,4,3,3 (the reverse of the
affine F4 path 3,3,4,3). Thus the stated profile has no omitted
spectator or old-component contact.

A right 3 at 7 or 8 cannot have a left contact: a 3 makes the full
exact all-3 cycle; an internal 4 has the retained 4,3,3 star, and an
end 4 has a retained affine F4 path. Larger bonds worsen the positive
nonpositive witness. With no left contact, the full principal tree has
arms 1,2,5, the exact affine E8 diagram. With no right contact, internal
left 3 contacts produce two-fork affine D diagrams; the two possible
end-4 extensions are affine as well. The remaining 3 at node 1 is
excluded by the connected minimal nonpositive subdiagram: it must
contain c and 5, while their unique old path has six vertices.
Spectators attach only through c and cannot shorten it. The Lanner
order bound therefore applies exactly as claimed. If there is no D8
contact, the direction-5 derivative is zero and the wall cannot close
that actual edge.

The argument uses positive definiteness of specified actual
replacements, or an exact affine principal diagram in the compact
polytope. It does not forbid an arbitrary indefinite full-star diagram.
In particular, the surviving 4,3,3 triangle is correctly retained.

For its unique profile I recovered the unit inverse entries
`C44=8, C55=10, C45=8`. The retained white residual at endpoint 5 is
`1-(1/2)(8/5)=1/5`, giving
`cosh L5=(4 sqrt2+5)/sqrt2=4+5/sqrt2>7` and hence
`tau5^2>3/4`. The pinned common-5/6 theorem applies with center 6
and neighbor 5 in D8. Thus energy at most 4/7 forces both pairs to
have distinct actual opposites, without any assumption on spectator
types or their number.

## Three-direction analytic contradiction

The inverse entries give correlations `q=sqrt(4/5)`, `k=sqrt(5/6)`
and `u4.u6=q*k`. Consequently the displayed orthonormal frame and
`u4=q*u5+t0/sqrt5` are exact. A minimizing covering vector exists;
its projections h_i are positive and below one. Monotonicity of F
justifies passing from the full tangents of the actual lengths to
`F(h4)F(h5)>=q` and `F(h5)F(h6)>=k`.

I independently derived the second denominator as
`(1+p)^2+4(1-lambda)b^2`, with `p=h5*h6>0`. Thus the concave psi
inequality and its tangent give
`psi'(p(k))=4 sqrt(1-k)/k` and
`p>=p(k)+(k sqrt(1-k)/2)b^2`.
Together with the norm bound this gives precisely
`T=(1+alpha)b^2+lambda*z^2<=lambda*(4/7)-p(k)`.
No sign of b or z is assumed; absolute values in the subsequent upper
bound preserve the inequality. If the right side were negative,
feasibility would already fail.

Expanding h5 squared while retaining the negative squared terms gives
exactly `A^2=8/105`, `B^2=4/35` and losses
`q*k*b^2+q*lambda*z^2`. Weighted Cauchy and the two loss ratios give
the stated bound in T. This step retains the extra difference penalty;
replacing the psi condition by only a product floor would lose it.

I checked the terminating-rational seals by exact hand algebra. In
particular, the bound for C reduces to
`59248/315609<47/250`, with positive cleared difference 21623.
The tight delta comparison has the displayed positive cleared
difference 26300653. The upper bounds for q*lambda*(4/7), the loss
ratio, p(q), and the increasing function on [0,1/400] have the required
strict signs. Finally `(271/12500)^2>47/100000`, and the resulting
upper coefficient is exactly
`489/1000+271/12500-71/40000=101781/200000<509/1000<p(q)`.
This contradicts the independently derived product lower bound.

## Scope and input boundary

The universal conclusion `E_D8>4/7=8/14` with arbitrary finite
spectators is proved within the named geometric/classification inputs.
Its proposed restricted-forest aggregation is compatible with the
existing D7 and smaller-component rates and the unchanged exceptional
defect accounting. It leaves the stated large A/B and D>=9 families;
there is no unrestricted compact dimension ceiling in this result.
I found no unresolved mathematical gap and request no source change.

Inputs read and rehashed:

| Source | SHA256 |
| --- | --- |
| `results/compact_dn_center_neighbor_common_opposite_energy_20261020.md` | `e371611c4903ad362654ee0d4656457140cdc5e4374b4e19a2697718c2d679b1` |
| `results/compact_universal_e6_e7_e8_component_energy_20261020.md` | `acf261d97f6eb1a11a64502e17e079ef7f418b9746fbf85dbdf6915750fe9a18` |
| `results/compact_actual_wall_reciprocal_shadow_20261020.md` | `c913ec3bd8e0f26a8c679a1ede2b98b251ad3f36efeae60486213ddfa3e8badc` |
| `results/compact_farthest_vertex_energy_obstruction_20261020.md` | `50c64e5a521d3e4e0f277a8ece02d3a29aa4ae8b71485cf3088211e587f45bff` |

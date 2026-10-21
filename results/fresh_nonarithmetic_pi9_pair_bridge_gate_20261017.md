# A complete two-wall bridge obstruction at the fixed π9 edge

For the frozen 60-wall chamber, no coherent pair of added Coxeter walls `A,B` can realize all six finite corners

- `S6+[58,59,A]` and `S6+[49,58,A]`;
- `S6+[30,59,B]` and `S6+[30,49,B]`;
- `S6+[49,A,B]` and `S6+[59,A,B]`,

where `S6=[0,1,2,3,5,6]`. These are the four old-edge terminations and two mutual corners of the proposed bridge across the exterior end near old49 and old59. Every finite label and every coefficient field is allowed. This is a restriction on that precise bridge combinatorics, not an obstruction to completing the chamber with more walls or to replacing the chosen π9 wall.

The [checker](../scripts/fresh_nonarithmetic_pi9_pair_bridge_gate_20261017.py) and [certificate](fresh_nonarithmetic_pi9_pair_bridge_gate_20261017.json) use the frozen [60-wall construction](fresh_nonarithmetic_pi8_quad_completion_20261017.md) and [complete π9 endpoint domain](fresh_nonarithmetic_pi9_endpoint_domain_20261017.md). The physical interval engine is shared with those inputs. This is an exact primary certificate with a fresh saved replay; no separate arithmetic review is claimed.

## Complete endpoint recipes

The first endpoint forces `A` to be right to `H4+I2(9)`. Its contacts with the remaining `A2` are zero, or one bond of label3,4,5 to either node. In its second endpoint,49 and58 are isolated orthogonal factors. Writing `θ_m=2cos(π/m)`, positive definiteness requires

`cost_A2 + θ_49²/2 < 2`.

For a nonzero A2 contact of label3,4,5, the possible labels to49 are respectively `{2,3,4,5}`, `{2,3}`, and `{2}`. The zero A2 profile allows every finite label to49. Thus there are fourteen nonzero-profile recipes, one zero-profile/right49 recipe, and one unbounded zero-profile tail covering all labels at least3.

Both required old-edge endpoints for `B` have type `H4+A2+I2(5)`. The inverse diagonal of the dihedral block is

`d5 = 2/(4−φ²)`.

The strict inequalities `2/3+d5>2` and `2d5>2` exclude mixing a nonzero A2 contact with a dihedral contact, and exclude labels at least4 on either dihedral node. Two simultaneous nonzero contacts to the dihedral factor also exceed the norm budget. Intersecting the two endpoint domains leaves exactly eleven recipes: zero profile with49 label2 or3; one A2 bond3/4/5 with49 right; a bond3 to30 with49 right; or a bond3 to59 with49 label2 or3.

For either endpoint frame, the old finite endpoint and old55 give the complete Lorentz2 kernel. The surviving branch has

`r(t)=p+sqrt((2−p²)/2)*(t*A0+B0/t)`, `t>1`.

The opposite orientation is obtuse to55. The endpoint on the surviving branch is proportional to `t*A0−B0/t`; its product with55 forces `t≥1`. At1, the added wall would pass through the old finite vertex. In an acute coherent extension this is impossible: that finite vertex persists strictly behind every additional wall, by its positive-definite incident Gram and a common interior point. Hence `t>1`.

Every listed49 contact gives a quadratic with one positive and one negative root; the positive root lies above1. The checker reconstructs that root, its norm, every prescribed contact, and all remaining old-wall gates. Of the fifteen finite A recipes, five survive. Of the eleven B recipes, five survive. Rejected recipes have either a strictly positive old-wall product or a strict bracket between consecutive cosine thresholds; no finite-angle cutoff is used to reject them.

## The unbounded tail and the remaining pairs

For the zero A2 profile, write `A(t)=t*A0+B0/t`. Its product with49 has coefficients `a49<0<b49`, so its negative magnitude strictly increases with `t`. Every label at least3 therefore has `t≥t3`, where `t3` is the exact label3 event.

For each of the five surviving B walls, the checker proves

`f(t)=A(t)·B=a*t+b/t`, `a<0`, `f(t3)<−2`, `f'(t3)<0`.

If `b≥0`, the derivative is always negative. If `b<0`, the derivative decreases toward `a` as `t` increases. Thus the entire tail has `A(t)·B<−2`; it cannot form either required finite mutual corner. This argument includes all real angle magnitudes at least1, hence all finite labels at least3.

Among the remaining25 pairs,23 also have mutual product strictly below−2. The only two pairs which could intersect have negative products whose magnitudes lie strictly in the gaps `(θ10,θ11)` and `(θ12,θ13)`. Neither is a Coxeter angle. The exact rational enclosures for every comparison are saved in the certificate.

The proof does not require a completeness claim about a numerical sample. The exploratory label≤1000 scan suggested the reduction, but it is not an input to this certificate.

## Scope and validation

Generation and a fresh saved replay pass. The certificate binds the frozen source and data dependencies by SHA256. There is no new chamber, new cap, or finite-volume claim. In particular, a three-wall bridge, another topology, or a changed value of the constructively chosen58–59 angle remains a legitimate next construction direction.

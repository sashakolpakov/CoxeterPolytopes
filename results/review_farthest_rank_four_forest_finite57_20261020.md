# Outside-author review of the rank-four farthest forest exclusion

20 October 2026. **PASS**, full written source/mathematical review.
No scientific program, enumeration, optimizer, Lean or TeX job was run.
The reviewer did not author or derive this assembly.

Reviewed in full:
`results/compact_farthest_rank_four_forest_exclusion_20261020.md`,
SHA256 `68ad62a619603eb7e6c2acc5689c901dde3f9d0897a7c1580b3deb948f5832b6`.
I also read its pinned small-A, B/D/F and H inputs, especially the
actual-wall alternatives, A1/A2 exception classifications and propagation
arguments. Their hashes agree with the source table: `801581c4...`,
`2eb27428...`, `3fc54634...`; the orthogonal-opposite input is
`4216d06c...`. My separate full H3/H4 review remains the detailed audit
of that component proof.

The new dihedral classification checks independently. Small energy
forces one actual opposite. For two contacts, the actual white residual
gives the displayed cosh formula, and the minimum exceeds11/8; an
individual half-tangent square then exceeds3/19>1/7. With one contact,
the bound is exactly `q(1-7k^2/16)<=1`. Its minimum
`(9d+7)/64` excludes every old label m>=8. The four remaining old labels
have the stated native-contact cutoffs. The m4/m5 cases have height
below1/4 and isolate; m6 is positive-null; m7 permits only a single
native3 contact and has Q<3/16, hence no spectator. Dotted contacts
are included by the larger-contact exclusion. No finite label tail or
multiple-contact case is omitted.

The old-I2(7)/new3 floor is exact: the mate has
`cosh^2 L>1623/1131>(109/91)^2`, with cross-product difference2652.
Thus its energy exceeds9/100. Its height square lies between0 and1/3.
The distinct old-A2/new7 defect also exceeds9/100 through the mate's
`tanh^2 L>164/541>3600/11881`, with difference884. These are different
defects and are never silently identified.

For precision, the disjoint-support wall product in Section4 is
`-sqrt(h_1 h_2)`, since h denotes squared height. The stated bounds
`0<h_i<1/2` put that product strictly in(-1/2,0). This proves all the
claimed incompatibilities, including the J short-A1 wall, whose two
supported components are themselves A1 factors and hence disjoint from
the exceptional old-I2(7) block. This is only a notation clarification;
the source's gap conclusion is unchanged.

The elementary finite-type exhaustion is complete. Cycle and star
tests have the correct Gram signs, and the rank-four path determinant
leaves exactly A4,B4,F4,H4. Under the propagated half-tangent floor1/5,
the I2 covering vector gives energy>(1/5), so it clears the required
rank-two rate. A4 uses its actual positive covering vector; the other
rank-three/four components use their intrinsic estimates. H4 is indeed
the weakest per-rank rate11/144. No common covering rate is falsely
assigned to every finite type.

All final cases exhaust the forest. The I2 defect excludes both the
A2 defect and every short A1, leaving the strict lower bound713/700.
Without it, the A2 defect or a small isolated-A1 height is excluded
by propagation. The remaining isolated-A1 debit is only2/35 and
`15/14-2/35=71/70`; the funded J case needs no debit. These inequalities
compare additive component energies at the SAME actual farthest vertex
with `||x||^2<1`.

No flaw was found. The conclusion is the stated pointwise restriction
for n>=15, with arbitrary other vertices. It supplies neither a uniform
dimension ceiling nor restrictions on larger connected finite factors.

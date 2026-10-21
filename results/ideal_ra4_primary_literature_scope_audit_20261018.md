# Primary-literature and scope audit: ideal24-cell and compact120-cell

18 October 2026. This is a bounded literature/source audit, not a new proof
or a complete bibliography. No mathematical program, optimizer, compilation
or Lean process was run. Existing frozen files were not changed.

**Entirely ideal volume integrality.** Kolpakov, *On the optimality of the
ideal right-angled24-cell*, AGT12(2012), Lemma2, printed pp1943–1944, proves

    Vol(P)=(f0-f3+4)pi^2/3.

The paper explicitly derives this from Zehrt, *The covolume of discrete
subgroups of Iso(H^{2m})*, Discrete Math.309(2009),2284–2291, Corollary4.2.
The verified statement gives pi^2/3-integrality. Equivalently,
chi_orb(P)=(f0-f3+4)/4. It does not itself prove f0-f3 is divisible by4.
This audit did not verify a stronger universal4pi^2/3-integrality theorem
of Zehrt/Kellerhals or another author, nor an entirely ideal counterexample
to that stronger divisibility. It would be incorrect to promote the weaker
formula into the stronger assertion, or to call the stronger assertion false
without additional evidence. Kolpakov's Theorem3, printed p1944, does prove
the unique minimum4pi^2/3, attained by the regular ideal24-cell.
[Published paper](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-s.pdf),
[Zehrt DOI](https://doi.org/10.1016/j.disc.2008.05.001).

**Universal decomposition into reflected24-cells.** No primary theorem or
concrete entirely ideal right-angled H4 counterexample establishing or
refuting that statement was verified in this search. Its truth, falsehood,
and precise current open-problem status remain unresolved by this audit.
Doubling/gluing constructions establish a family, not an exhaustive
classification. Mixed finite/ideal Kerckhoff–Storm examples and non-right-
angled rectified-simplex examples do not answer the entirely ideal question.
The September2026 seed-gluing theorem in Allcock–Devlin–Felikson–Kontorovich–
Whitehead, *Arithmetic Polyhedra*, concerns arithmetic ideal right-angled
polyhedra in H3; it is not an H4 decomposition theorem.
[Primary preprint, Theorems1.4–1.5](https://arxiv.org/html/2609.05349v1).

**The compact right-angled120-cell.** Ma–Zheng, *Orientable hyperbolic
4-manifolds over the120-cell*, arXiv1801.08814v4(8 June2022), Remark2.7,
printed p9, states Vol(C)=34pi^2/3 and the face vector(600,1200,720,120).
For its reflection orbifold the ordinary right-angled face-stabilizer sum is

    chi_orb=1-120/2+720/4-1200/8+600/16=17/2.

Thus the volume also equals(4pi^2/3)(17/2). Rational orbifold Euler
characteristic is compatible with integer manifold Euler characteristics;
this compact example does not settle the entirely ideal divisibility question.
[Primary paper, v4 p9](https://arxiv.org/pdf/1801.08814v4#page=9).

The right-angled realization has dihedral angle pi/2. The Davis manifold
uses the distinct compact regular120-cell realization with angle2pi/5,
as Ma–Zheng§7,p32 and Garrison–Scott(2001),p1 explicitly distinguish.
The latter polytope is not itself a Coxeter reflection chamber:2pi/5 is not
pi/m. A regular compact120-cell with all angles pi/3 is not another case:
its four incident norm-one wall normals would have diagonal1 and all
off-diagonal entries-1/2, with a negative eigenvalue, contradicting a finite
spherical corner.
[Garrison–Scott primary paper, p1](https://arxiv.org/pdf/math/0107116#page=1).

Kolpakov's thesis§4.4, printed pp73–75, poses volume and facet minimality
of the right-angled120-cell as questions and supplies partial results.
That historical question does not establish present-day openness, minimum
volume, uniqueness, or universal120-cell decomposition. No later resolution
of those general assertions was verified here; none should be reported as a
theorem on the strength of this audit.
[Author's thesis, §4.4,p73](https://homeweb.unifr.ch/kellerha/pub/KolpakovThesis_Final.pdf).

**Our equality theorem has a narrower, internal status.** The
[gap proof](ideal_ra4_second_volume_20261018.md) uses only pi^2/3-integrality
and explicitly excludes the intermediate coefficients5,6,7. The
[equality proof](ideal_ra4_second_volume_equality_20261018.md) reconstructs
an octahedral cutting section specifically at volume8pi^2/3, from total
facet defect24. It does not assume universal24-cell decomposition or stronger
volume divisibility. Rereading identified no new flaw or such hidden
assumption. It remains an internally reviewed written structural argument
with exact checks of local finite facts, not a published theorem or an
end-to-end Lean proof. Its existing Lean scope concerns earlier conditional
integer-gap logic, not the universal geometric structure. External human
mathematical review and a separate priority assessment remain necessary.

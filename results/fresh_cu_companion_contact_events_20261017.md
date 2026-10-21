# Moving companions on the deletion-free C(u) tail

Varying the exact-admissible C(u) family finds17 numerical contact events
with10 intrinsic open ends and a changed four-edge B boundary. None improves
the exact63-wall/199-finite/10-open-end control. No numerical candidate is
promoted to an exact chamber or finite-volume example.

The [C(u) theorem](fresh_dotted_c10_domain_20261017.md) retains all63 control
walls for u>=3. In the residual Lorentz4 space, the old7/30 plane is positive
definite, and C(u) lies in its orthogonal Lorentz2 plane. Choose e0^2=-1,
e1^2=1, e0.e1=0 and gamma>0 so that

    C(u)=sqrt(2)*(gamma*u*e0+sqrt(1+gamma^2*u^2)*e1).
    t=gamma*u+sqrt(1+gamma^2*u^2),  t>=t0=t(3).

For fixed native contacts to7,30,C and a fixed A2 profile, let A0 be the
projection onto span(7,30), q the residual square, and r=A.C. The complete
norm equation is

    A=A0+(r/2)*C+epsilon*sqrt(lambda)*T,
    lambda=A0^2+r^2/2-q,
    T=((t+1/t)*e0+(t-1/t)*e1)/2,

with epsilon=+/-1 when lambda>0. Therefore A=A0+t*P+Q/t. Every old native
contact, including its fixed A2 offset, is a*t+b/t+c. Contact events reduce
to quadratic equations in t. This parameterization does not impose an upper
bound on u; its floating-point roots do not constitute interval coverage.

The inherited complete finite-corner classification supplies51 symbolic
recipes at frame7/30/C. The five recipes containing an isolated high bond
are sampled at6,8,12,24, giving66 literal recipes. Eight have no real norm
branch numerically; the other58 give116 moving branches. Event targets are
all other actual old sides, at labels2..5. High event labels are not covered.
The code also samples cells separated by native products0,-2 and the center
zero; one sampled cell is not a certified interval of Coxeter compatibility.

There are1,279 tested parameters,41 coherent finite-contact events and56
coherent cell/start samples. All accepted configurations pass every old
native contact, the unchanged strict center, a complete numerical target
cone and a separate positive-definite incidence enumeration with identical
finite supports. Their counts are:

|Finite rays|Exterior rays|Open ends|Accepted samples|
|---:|---:|---:|---:|
|29|7|11|60|
|29|7|13|20|
|30|6|10|17|

The17 ten-end events give A a finite contact to old47. Their A coordinates
fall into eight numerical clusters at tolerance1e-6: changing C sometimes
decorates the same A geometry.
Each closes the whole C face as a finite triangle, while the B boundary
becomes30–52–47–A–30. The four-ray exterior component on66 remains. This
changed incidence motivates the separate [four-cycle closure search](fresh_cu_four_cycle_closure_20261017.md).

There are147 explicitly unresolved evaluations:146 large-parameter norm or
frame residuals and one cone with a possible affine support. None is accepted
or used as an exclusion. In particular, tiny floating-point coefficients and
large quadratic roots cannot certify the asymptotic tail. The accepted event
coordinates themselves pass the stated residual guards.

No old wall is deleted, so the deletion debit is zero. In any exact coherent
acute realization with this common center, Stieltjes projection would preserve
all199 old finite vertices and every old compact face. This scout supplies
neither that exact realization nor a full ambient cone certificate.

[Source](../scripts/fresh_cu_companion_contact_events_20261017.py) and
[saved data](fresh_cu_companion_contact_events_20261017.json) record the
moving coefficients, all accepted cones, blockers and unresolved cases.
Root's complete narrow source/mathematical read passed the moving-frame
algebra, native offsets, center and cone/incidence comparison within this
numerical scope; it did not repeat the event/root/cone computations.
**FINAL/FROZEN:** primary generation and a separate fresh saved replay both terminated
exit0/PASS. A separate saved-data summary confirms the17 ten-end events,
the eight numerical A clusters and all complete C-face triangles; it does
not independently reconstruct numerical roots or cones.

    python3 scripts/fresh_cu_companion_contact_events_20261017.py --verify

"""Exact silver-unit Vinberg packets for Niemeier ``N(D10 E7^2)``.

The discriminant group is ``(Z/2)^2 + Z/2 + Z/2``.  Label the D10 classes
``0,v,s,c``.  The Niemeier glue is
``{(0,0,0),(s,1,0),(c,0,1),(v,1,1)}``.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import groupby
from typing import Sequence, TypeAlias

from coxeter_polytopes.exact.silver import SILVER_ONE, SILVER_ZERO, SQRT2, Silver
from coxeter_polytopes.search.bugaenko_vinberg import _fraction_inverse, _prefix_schur_complements
from coxeter_polytopes.search.niemeier_host_triage import affine_vertex_energies, cartan_gram, highest_root_marks
from coxeter_polytopes.search.silver_a17_e7 import _determinant
from coxeter_polytopes.search.silver_a17_e7 import E7_RESIDUE_INDICES
from coxeter_polytopes.search.silver_e8_cubes import _bounded_for_caps, _mat_vec, _quadratic

Vector: TypeAlias = tuple[Silver, ...]
RationalMatrix: TypeAlias = tuple[tuple[Fraction, ...], ...]

EPSILON = SILVER_ONE + SQRT2
UNIT_NORM_CLASSES = (SILVER_ONE, Silver(2))
D10_COMPONENT = ("D", 10)
E7A_COMPONENT = ("E", 7)
E7B_COMPONENT = E7A_COMPONENT

COMPLETE_UNIT_NORM_PROOF = (
    "The coefficient epsilon=1+sqrt(2) is a unit and N(D10 E7^2) is "
    "even unimodular. The primitive pairing-ideal divisibility q(r)|2I(r), "
    "together with the dyadic parity argument, leaves valuations zero and "
    "two. Totally-positive units are squares. Hence the complete canonical "
    "primitive mirror norms are 1 and 2."
)


@lru_cache(maxsize=1)
def scaled_niemeier_basis() -> tuple[tuple[int,...],...]:
    """Twice an HNF basis, in the simple-root basis of D10+E7+E7."""
    rows=[]
    for index in range(16):
        row=[0]*24;row[index]=2;rows.append(row)
    row=[0]*24
    for index in (0,2,4,6,8,13,15,16):row[index]=1
    rows.append(row)
    for index in range(17,23):
        row=[0]*24;row[index]=2;rows.append(row)
    row=[0]*24
    for index in (0,2,4,6,9,20,22,23):row[index]=1
    rows.append(row)
    return tuple(tuple(row) for row in rows)


@lru_cache(maxsize=1)
def niemeier_gram() -> RationalMatrix:
    root_gram=[]
    for component in (D10_COMPONENT,E7A_COMPONENT,E7B_COMPONENT):
        block=cartan_gram(component)
        old=len(root_gram);rank=len(block)
        for row in root_gram:row.extend([Fraction(0)]*rank)
        root_gram.extend([[Fraction(0)]*old+list(row) for row in block])
    basis=scaled_niemeier_basis()
    return tuple(tuple(sum(Fraction(basis[i][r]*root_gram[r][s]*basis[j][s],4) for r in range(24) for s in range(24) if root_gram[r][s]) for j in range(24)) for i in range(24))


def basis_audit() -> bool:
    gram = niemeier_gram()
    return _determinant(gram)==1 and all(x.denominator==1 for row in gram for x in row) and all(gram[i][i].numerator%2==0 for i in range(24))


def _coefficient_residues(weights: Sequence[Silver], functional: Sequence[int], modulus: int) -> tuple[int,int]:
    if not all(value.is_integral for value in weights): raise ValueError("nonintegral Dynkin weight")
    return tuple(sum(c*int(getattr(v,name)) for c,v in zip(functional,weights))%modulus for name in ("a","b"))


def d10_residue(weights: Sequence[Silver]) -> tuple[int,int]:
    if len(weights)!=10: raise ValueError("D10 needs ten weights")
    # bits v,s with s=omega_9 and c=omega_10=(v+s).
    values=[]
    for name in ("a","b"):
        z=[int(getattr(w,name)) for w in weights]
        values.append(((z[0]+z[2]+z[4]+z[6]+z[9])%2,(z[8]+z[9])%2))
    return tuple(values)


def e7a_residue(weights: Sequence[Silver]) -> tuple[int,int]:
    if len(weights)!=7: raise ValueError("E7A needs seven weights")
    return tuple(sum(int(getattr(weights[i],name)) for i in E7_RESIDUE_INDICES)%2 for name in ("a","b"))


def e7b_residue(weights: Sequence[Silver]) -> tuple[int,int]:
    return e7a_residue(weights)


def glue_membership(a:Sequence[Silver],d:Sequence[Silver],e:Sequence[Silver]) -> bool:
    try: ar,dr,er=d10_residue(a),e7a_residue(d),e7b_residue(e)
    except ValueError:return False
    return all(dr[p]==(ar[p][0]^ar[p][1]) and er[p]==ar[p][0] for p in (0,1))


@dataclass(frozen=True,slots=True)
class Weight:
    square:Silver
    weights:Vector
    coordinates:Vector
    residue:tuple


@dataclass(frozen=True,slots=True)
class Enumeration:
    component:tuple[str,int]
    target:Silver
    root_norm:Silver
    cap:Silver|None
    visited:tuple[int,...]
    cap_pruned:tuple[int,...]
    schur_pruned:tuple[int,...]
    vectors:tuple[Weight,...]


def _alphabet(target:Silver,norm:Silver)->tuple[Silver,...]:
    return tuple(v for v in _bounded_for_caps((2*target,2*target)) if v.sign(0)>=0 and all((2*target-v*v).sign(p)>=0 for p in (0,1)) and (2*v/norm).is_integral)


@lru_cache(maxsize=None)
def enumerate_weights(component:tuple[str,int],target:Silver,norm:Silver,cap:Silver|None)->Enumeration:
    gram=cartan_gram(component); inverse=_fraction_inverse(gram); schur=_prefix_schur_complements(inverse);marks=highest_root_marks(component);alphabet=_alphabet(target,norm);rank=component[1]
    visited=[0]*(rank+1);cap_pruned=[0]*(rank+1);schur_pruned=[0]*(rank+1);out=[]
    residue_fn=d10_residue if component==D10_COMPONENT else e7a_residue
    def extend(prefix:Vector):
        depth=len(prefix);visited[depth]+=1
        if depth:
            if cap is not None and (2*cap-sum((marks[i]*prefix[i] for i in range(depth)),start=SILVER_ZERO)).sign(0)<0:
                cap_pruned[depth]+=1;return
            minimum=_quadratic(prefix,schur[depth-1])
            if any((target-minimum).sign(p)<0 for p in (0,1)):
                schur_pruned[depth]+=1;return
        if depth==rank:
            square=_quadratic(prefix,inverse);out.append(Weight(square,prefix,_mat_vec(inverse,prefix),residue_fn(prefix)));return
        for value in alphabet:extend((*prefix,value))
    extend(())
    return Enumeration(component,target,norm,cap,tuple(visited),tuple(cap_pruned),tuple(schur_pruned),tuple(sorted(out,key=lambda x:(x.square,x.weights))))


@dataclass(frozen=True,slots=True)
class Control:
    norm:Silver;height:Silver;control:Silver;target:Silver


@lru_cache(maxsize=None)
def controls(height_bound:int=15)->tuple[Control,...]:
    levels=set()
    for norm in UNIT_NORM_CLASSES:
        for x in _bounded_for_caps((norm*height_bound/EPSILON,norm/(-EPSILON))):
            if x.sign(0)<=0 or not (2*EPSILON*x/norm).is_integral:continue
            target=norm+EPSILON*x*x;height=EPSILON*x*x/norm
            if target.sign(1)>=0 and (Silver(height_bound)-height).sign(0)>=0:levels.add(Control(norm,height,x,target))
    return tuple(sorted(levels,key=lambda x:(x.height,x.norm,x.control)))


@dataclass(frozen=True,slots=True)
class Candidate:
    norm:Silver;height:Silver;control:Silver;target:Silver;a:Weight;d:Weight;e:Weight
    @property
    def glue_nontrivial(self):return self.a.residue!=((0,0),(0,0)) or self.d.residue!=(0,0) or self.e.residue!=(0,0)
    @property
    def three_factor_glue(self):return self.a.residue!=((0,0),(0,0)) and self.d.residue!=(0,0) and self.e.residue!=(0,0)


def primitive(root:Candidate)->bool:
    if not (root.control/SQRT2).is_integral:return True
    blocks=tuple(tuple(v/SQRT2 for v in w.weights) for w in (root.a,root.d,root.e))
    if not all(v.is_integral for block in blocks for v in block):return True
    return not glue_membership(*blocks)


@lru_cache(maxsize=None)
def enumerate_level(level:Control,caps_active:bool)->tuple[Candidate,...]:
    cap=level.control if caps_active else None
    blocks={comp:enumerate_weights(comp,level.target,level.norm,cap).vectors for comp in (D10_COMPONENT,E7A_COMPONENT,E7B_COMPONENT)}
    de=defaultdict(list)
    for d in blocks[E7A_COMPONENT]:
        for e in blocks[E7B_COMPONENT]:de[d.square+e.square].append((d,e))
    out=[]
    for a in blocks[D10_COMPONENT]:
        for d,e in de.get(level.target-a.square,()):
            if glue_membership(a.weights,d.weights,e.weights):
                root=Candidate(level.norm,level.height,level.control,level.target,a,d,e)
                if primitive(root):out.append(root)
    return tuple(sorted(out,key=lambda r:(r.a.weights,r.d.weights,r.e.weights)))


def pair(x:Candidate,y:Candidate)->Silver:
    return -EPSILON*x.control*y.control+sum((x.a.weights[i]*y.a.coordinates[i] for i in range(10)),start=SILVER_ZERO)+sum((x.d.weights[i]*y.d.coordinates[i] for i in range(7)),start=SILVER_ZERO)+sum((x.e.weights[i]*y.e.coordinates[i] for i in range(7)),start=SILVER_ZERO)


@dataclass(frozen=True,slots=True)
class Packet:
    height:Silver;controls:tuple[Control,...];candidate_count:int;rejected_count:int;roots:tuple[Candidate,...]
@dataclass(frozen=True,slots=True)
class Prefix:
    height_bound:int;controls:tuple[Control,...];packets:tuple[Packet,...];accepted_roots:tuple[Candidate,...]


@lru_cache(maxsize=1)
def height_zero()->tuple[Candidate,...]:
    return tuple(root for norm in UNIT_NORM_CLASSES for root in enumerate_level(Control(norm,SILVER_ZERO,SILVER_ZERO,norm),False))


@lru_cache(maxsize=None)
def run_prefix(height_bound:int=15,nonempty_limit:int|None=2)->Prefix:
    if not basis_audit():raise AssertionError("basis failed")
    h0=height_zero()
    if len(h0)!=3 or any(r.norm!=2 for r in h0):raise AssertionError("height zero wrong")
    cs=controls(height_bound);accepted=[];packets=[];nonempty=0
    for height,it in groupby(cs,key=lambda x:x.height):
        pcs=tuple(it);candidates=tuple(r for c in pcs for r in enumerate_level(c,nonempty>0));survivors=tuple(r for r in candidates if all(pair(r,w).sign(0)<=0 for w in accepted))
        if any(pair(x,y).sign(0)>0 for i,x in enumerate(survivors) for y in survivors[i+1:]):raise AssertionError("packet not acute")
        accepted.extend(survivors);packets.append(Packet(height,pcs,len(candidates),len(candidates)-len(survivors),survivors))
        if survivors:nonempty+=1
        if nonempty_limit is not None and nonempty==nonempty_limit:break
    return Prefix(height_bound,cs,tuple(packets),tuple(accepted))


@dataclass(frozen=True,slots=True)
class Geometry:
    rays:int;finite:int;cusps:int;hyperideal:int
def first_geometry()->Geometry:
    distribution={Fraction(0):1}
    for comp in (D10_COMPONENT,E7A_COMPONENT,E7B_COMPONENT):
        updated={}
        for energy,m in distribution.items():
            for v in affine_vertex_energies(comp):updated[energy+v]=updated.get(energy+v,0)+m
        distribution=updated
    finite=sum(m for e,m in distribution.items() if (Silver(e)-EPSILON).sign(0)<0);cusps=sum(m for e,m in distribution.items() if (Silver(e)-EPSILON).sign(0)==0);total=sum(distribution.values())
    return Geometry(total,finite,cusps,total-finite-cusps)


@dataclass(frozen=True,slots=True)
class Witness:
    weight:Fraction;norm:Silver;pairings:tuple[Silver,...];strict:bool;hyperideal:bool
@lru_cache(maxsize=1)
def strict_witness()->Witness:
    prefix=run_prefix(15);weight=Fraction(5,3);inverse=_fraction_inverse(cartan_gram(D10_COMPONENT));norm=Silver(weight*weight*inverse[0][0])-EPSILON
    values=tuple(-EPSILON*r.control+Silver(weight)*r.a.coordinates[0] for r in prefix.accepted_roots)
    return Witness(weight,norm,values,all(v.sign(0)<0 for v in values),norm.sign(0)>0)


__all__=[name for name in globals() if not name.startswith("_")]

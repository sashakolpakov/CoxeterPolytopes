#!/usr/bin/env python3
"""Explore rank-ten golden path hosts with exact bounded root pools.

This is a discovery search, not a nonreflectivity certificate.  Every
returned wall and cone is checked exactly.  The optional complete norm
alphabet does not prove that the initial A_n chamber is a chamber of every
height-zero mirror, so a negative result remains a bounded search result.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import groupby
import json
from math import isqrt
from pathlib import Path
import time

from coxeter_polytopes.exact.golden import Golden as G, ONE, ZERO, PHI
from coxeter_polytopes.exact.matrix import inertia, determinant
from coxeter_polytopes.search.golden_path_host import (
    bounded_golden_integers, spatial_inverse_entry,
    spatial_inverse_energy, _canonical_totally_positive_associate,
)
from coxeter_polytopes.verify.chamber import diagnose_chamber_cone, verify_hyperideal_extreme_ray
from sympy import divisors


def complete_norm_pool(det):
    """Canonical TP divisors of 2 det, containing every primitive root norm.

    For primitive r, adj(B) B r=det(B)r bounds the pairing ideal by det(B).
    Reflection integrality gives q | 2(B(r,L)), hence q | 2 det(B).
    Z[phi] has class number one and every TP unit is a square.
    The canonical associate satisfies sqrt(N)/phi < q <= sqrt(N)*phi;
    hence |b|=|q-q'|/sqrt(5) < sqrt(N), for q=a+b phi of norm N.
    """
    twice=2*det
    norms=set()
    for norm in divisors(abs(int(twice.norm()))):
        for b in range(-isqrt(norm)-1,isqrt(norm)+2):
            disc=5*b*b+4*norm
            s=isqrt(disc)
            if s*s != disc:
                continue
            for signed in (-s,s):
                if (signed-b)%2:
                    continue
                q=G((signed-b)//2,b)
                if q.sign()<=0 or q.sign(1)<=0 or not (twice/q).is_integral:
                    continue
                q=_canonical_totally_positive_associate(q)
                norms.add(q)
    return tuple(sorted(norms))


class Host:
    def __init__(self, dimension, coefficient):
        self.n = dimension
        self.c = coefficient
        self.qrho = self.c * (2*(self.n+1)-self.n*self.c)/(self.n+1)
        self.gram = tuple(tuple(
            (2*self.c if i == 0 else G(2)) if i == j else
            (-self.c if min(i,j) == 0 else -ONE) if abs(i-j) == 1 else ZERO
            for j in range(self.n+1)) for i in range(self.n+1))
        assert (inertia(self.gram, embedding=0).positive, inertia(self.gram, embedding=0).negative) == (self.n,1)
        assert (inertia(self.gram, embedding=1).positive, inertia(self.gram, embedding=1).negative) == (self.n+1,0)

    def pair(self, u, v):
        return 2*self.c*u[0]*v[0] + sum((2*u[i]*v[i] for i in range(1,self.n+1)),ZERO) - self.c*(u[0]*v[1]+u[1]*v[0]) - sum((u[i]*v[i+1]+u[i+1]*v[i] for i in range(1,self.n)), ZERO)

    def basis(self, i):
        return tuple(ONE if j == i else ZERO for j in range(self.n+1))

    def pair_basis(self, r):
        return (2*self.c*r[0]-self.c*r[1], self.c*(-r[0])+2*r[1]-r[2], *(2*r[i]-r[i-1]-r[i+1] for i in range(2,self.n)), 2*r[-1]-r[-2])

    @lru_cache(maxsize=None)
    def energies(self, target, q):
        alphabet = tuple(x for x in bounded_golden_integers((2*target,2*target)) if x.sign() >= 0 and (q*x/2).is_integral)
        solutions = []
        nodes = 0
        square={x:x*x for x in alphabet}
        def rec(prefix,uncentered,moment):
            nonlocal nodes
            nodes += 1
            depth=len(prefix)
            # (A_d^-1)_ij=min(i,j)-ij/(d+1), with one-based indices.
            # Maintain the two sums, avoiding a quadratic matrix product
            # at every node while retaining the exact Schur lower bound.
            slack=(depth+1)*(target-uncentered)+moment*moment
            if any(slack.sign(k)<0 for k in (0,1)):
                return
            if depth == self.n:
                if not slack:
                    solutions.append(prefix)
                return
            for x in alphabet:
                rec((*prefix,x),uncentered+2*x*moment+(depth+1)*square[x],moment+(depth+1)*x)
        rec((),ZERO,ZERO)
        return tuple(solutions), nodes

    def candidates(self, norms, height):
        controls=[]
        candidates=[]
        for q in norms:
            for x in bounded_golden_integers((-height*q/self.qrho,q/self.qrho)):
                if x.sign() <= 0:
                    continue
                h=-self.qrho*x*x/q
                target=q-self.qrho*x*x
                if (h-height).sign()>0 or target.sign(1)<0:
                    continue
                solutions,nodes=self.energies(4*target/(q*q),q)
                accepted=0
                for z in solutions:
                    rhs=(self.c*x-q*z[0]/2,*(-q*zz/2 for zz in z[1:]))
                    r=(x,*(sum((G(spatial_inverse_entry(i,j,self.n))*rhs[j] for j in range(self.n)),ZERO) for i in range(self.n)))
                    if not all(v.is_integral for v in r):
                        continue
                    assert self.pair(r,r) == q
                    if not all((2*v/q).is_integral for v in self.pair_basis(r)):
                        continue
                    candidates.append((h,q,r))
                    accepted+=1
                controls.append((q,x,h,nodes,len(solutions),accepted))
        return sorted(candidates),controls

    def search(self,norms,height):
        roots,controls=self.candidates(norms,height)
        walls=[self.basis(i) for i in range(1,self.n+1)]
        accepted=[]
        directions=set()
        for h,group in groupby(roots,key=lambda x:x[0]):
            packet=[]
            for _,q,r in group:
                direction=tuple(v/r[0] for v in r)
                if direction in directions or any(self.pair(r,w).sign()>0 for w in walls):
                    continue
                directions.add(direction)
                packet.append((h,q,r))
            assert all(self.pair(u[2],v[2]).sign()<=0 for i,u in enumerate(packet) for v in packet[i+1:])
            accepted.extend(packet)
            walls.extend(r for _,_,r in packet)
        rho=(ONE,*(self.c*G(self.n+1-i)/(self.n+1) for i in range(1,self.n+1)))
        assert self.pair(rho,rho) == self.qrho
        epsilon=G(Fraction(1,2))
        perturbation=(ZERO,*(G(Fraction(-i*(self.n+1-i),2)) for i in range(1,self.n+1)))
        for _ in range(100):
            interior=tuple(rho[i]+epsilon*perturbation[i] for i in range(self.n+1))
            if self.pair(interior,interior).sign()<0 and all(self.pair(interior,w).sign()<0 for w in walls):
                break
            epsilon/=2
        else:
            raise AssertionError('no strict timelike point found')
        geometry=None
        witness=None
        if len(walls)>=self.n+1:
            g=diagnose_chamber_cone(tuple(walls),self.gram,interior,self.n)
            geometry={"finite_vertices":len(g.finite_vertices),"cusps":len(g.cusps),"hyperideal_rays":len(g.hyperideal_rays),"opposite_sheet_rays":len(g.opposite_sheet_rays)}
            simple_ray=next((ray for ray in g.hyperideal_rays if len(ray)==self.n),None)
            if simple_ray is not None:
                w=verify_hyperideal_extreme_ray(tuple(walls),self.gram,simple_ray,self.n)
                witness={"active_facets":list(w.active_facets),"vector":[x.to_json() for x in w.ray],"norm":w.norm.to_json(),"wall_pairings":[x.to_json() for x in w.wall_pairings]}
        return {
            "dimension":self.n,"coefficient":self.c.to_json(),"height_bound":height.to_json(),
            "norm_pool":[q.to_json() for q in norms],"determinant":determinant(self.gram).to_json(),
            "control_count":len(controls),"recursion_nodes":sum(c[3] for c in controls),
            "candidate_count":len(roots),"wall_count":len(walls),"geometry":geometry,
            "hyperideal_witness":witness,
            "positive_walls":[{"height":h.to_json(),"norm":q.to_json(),"root":[v.to_json() for v in r]} for h,q,r in accepted],
            "strict_interior":[v.to_json() for v in interior],
            "all_walls":[[v.to_json() for v in r] for r in walls],
            "status":"SEARCH_PREFIX_ONLY" if geometry is None or geometry['hyperideal_rays'] or geometry['opposite_sheet_rays'] else "FINITE_VOLUME_CANDIDATE_REQUIRES_INDEPENDENT_AUDIT",
            "full_mirror_arrangement_audited":False,
        }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dimension',type=int,default=9)
    p.add_argument('--coefficient',default='2,2')
    p.add_argument('--height',type=int,default=10)
    p.add_argument('--complete-norm-alphabet',action='store_true')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    c=G(*(int(x) for x in a.coefficient.split(',')))
    norms={_canonical_totally_positive_associate(G(q)) for q in (1,2,4,8)}
    norms|={_canonical_totally_positive_associate(q*c) for q in (1,2,4)}
    host=Host(a.dimension,c)
    if a.complete_norm_alphabet:
        norms=complete_norm_pool(determinant(host.gram))
    start=time.monotonic()
    result=host.search(tuple(sorted(norms)),G(a.height))
    result['seconds']=round(time.monotonic()-start,3)
    result['primitive_mirror_norm_alphabet_complete']=a.complete_norm_alphabet
    result['initial_height_zero_chamber_completeness_audited']=False
    result['norm_pool_scope']='All canonical TP divisors of 2det(B).' if a.complete_norm_alphabet else 'Explicit exploratory norm pool.'
    result['negative_result_scope']='The bounded accepted cone is not finite volume; no host nonreflectivity consequence.'
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in ('dimension','coefficient','height_bound','control_count','recursion_nodes','candidate_count','wall_count','geometry','seconds','status')},sort_keys=True))


if __name__=='__main__':
    main()

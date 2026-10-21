#!/usr/bin/env python3
"""Literal inherited obstruction supports, used as proposal-search objectives."""
from itertools import combinations
from fresh_h9_rank5_gale_20261014_core import LocalGates, canonical


class SupportGates(LocalGates):
    def __init__(self):
        super().__init__()
        self.all_cache={}

    def all(self,mm):
        """All recognized gates on unions of two order3..5 missing faces.

        This is a complete scan of this finite pattern dictionary's specified
        carriers, not a complete list of all possible Coxeter obstructions.
        Every returned row is itself an exact, literal applicable gate.
        """
        key=tuple(mm)
        if key in self.all_cache:
            return self.all_cache[key]
        supports={}
        for a,b in combinations(mm,2):
            sa,sb=a.bit_count(),b.bit_count()
            if not (3<=sa<=5 and 3<=sb<=5):
                continue
            u=a|b
            if u in supports:
                continue
            inside=[c for c in mm if c&u==c]
            if not (a&b) and set(inside)=={a,b} and sorted((sa,sb)) in ([3,4],[3,5],[4,4],[4,5],[5,5]):
                supports[u]=dict(gate='two_lanner_blocks',blocks=[a,b],support=u)
                continue
            if u.bit_count()>9 or any(m.bit_count()==2 for m in inside):
                continue
            vv=[i for i in range(15) if u>>i&1]
            red=tuple(sorted(sum(1<<k for k,j in enumerate(vv) if m>>j&1) for m in inside))
            ck=len(vv),red
            if ck not in self.cache:
                self.cache[ck]=canonical(red,len(vv))
            can,order=self.cache[ck]
            if (len(vv),can) in self.bad:
                path,ident=self.bad[len(vv),can]
                supports[u]=dict(gate='literal_induced_pattern',source=path,pattern=ident,
                                 support=u,permutation=[vv[j] for j in order])
        out=[supports[u] for u in sorted(supports)]
        self.all_cache[key]=out
        return out
